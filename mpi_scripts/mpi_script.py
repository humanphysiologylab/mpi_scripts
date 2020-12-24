from mpi4py import MPI
import os
import json
import time
from collections import OrderedDict
import pickle
from datetime import datetime
import copy

import numpy as np
import pandas as pd
# from scipy.integrate import solve_ivp

import sklearn.metrics

import ctypes
from numba import njit

import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

from src.helpers import get_value_by_key, update_array_from_kwargs,\
    calculate_RMSE_balanced, calculate_RMSE, \
    batches_from_list, value_from_bounds, argmax_list_of_dicts,\
    calculate_composite_RMSE_V_CaT, Timer, find_index_first, flatten_iterable, strip_comments
from src.algorythm.ga import do_step


def create_genes_dict_from_config(config):
    genes_dict = {ec_name:
                      {p_name: dict(bounds=p['bounds'],
                                    gamma=p['gamma'],
                                    is_multiplier=p.get("is_multiplier", False))
                       for p_name, p in ec['params'].items() if isinstance(p, dict)}
                  for ec_name, ec in config['experimental_conditions'].items()}

    return genes_dict


def create_constants_dict_from_config(config):
    constants_dict = {
        ec_name: {p_name: value for p_name, value in ec['params'].items() if isinstance(value, (int, float))}
        for ec_name, ec in config['experimental_conditions'].items()}

    return constants_dict


def generate_organism(genes_dict, genes_m_index, state):

    genes = [value_from_bounds(gene_params['bounds'], log_scale=gene_params['is_multiplier'])
             if gene_name not in state.index else state[exp_cond_name][gene_name]
             for exp_cond_name, exp_cond in genes_dict.items() for gene_name, gene_params in exp_cond.items()]

    genes = pd.Series(data=genes, index=genes_m_index)

    organism = dict(genes=genes,
                    state=state.copy())
    return organism


def init_population(config):

    population = [generate_organism(config['runtime']['genes_dict'],
                                    config['runtime']['m_index'],
                                    config['runtime']['states_initial'])
                  for _ in range(config['runtime']['n_organisms'])]
    return population


def init_population_from_backup(backup, config):

    genes_dict = config['runtime']['genes_dict']

    assert (len(backup) == config['runtime']['n_organisms'])
    assert (len(backup[0]['genes']) == sum(map(len, genes_dict.values())))
    assert (backup[0]['state'].shape == config['runtime']['states_initial'].shape)

    population = [dict(genes=pd.Series(data=organism['genes'], index=config['runtime']['m_index']),
                       state=pd.DataFrame(data=organism['state'], columns=config['runtime']['states_initial'].columns,
                                          index=config['runtime']['states_initial'].index))
                  for organism in backup]

    return population


def run_model_ctypes(S, C, config):

    stim_period = C[config['stim_period_legend_name']]
    t_sampling = config['t_sampling']

    if 'n_beats' in config and 't_run' not in config:
        n_beats = config['n_beats']
    elif 'n_beats' not in config and 't_run' in config:
        t_run = config['t_run']
        n_beats = np.ceil(t_run / stim_period).astype(int)
        if n_beats % 2 == 0:
            n_beats += 1
        #print(stim_period, n_beats, stim_period * n_beats)
    else:
        print('Invalid config: check n_beats & t_run entries.', fflush=True)
        exit()

    tol = config.get('tol', 1e-6)

    n_samples_per_stim = int(stim_period / t_sampling)
    output = np.empty((n_samples_per_stim * n_beats + 1, len(S)))

    status = run_model_ctypes.model.run(S.values.copy(), C.values.copy(),
                                        n_beats, t_sampling, tol, output)

    output = output[-n_samples_per_stim - 1:].T  # last beat

    return status, output


def update_phenotype_state(organism, config):

    organism['phenotype'] = dict()

    legend = config['runtime']['legend']
    genes_dict = config['runtime']['genes_dict']
    constants_dict = config['runtime']['constants_dict']

    for exp_cond_name in config['experimental_conditions']:

        if exp_cond_name == 'common':
            continue

        genes_current = organism['genes'][['common',
                                           exp_cond_name]]
        constants_dict_current = {**constants_dict['common'],
                                  **constants_dict[exp_cond_name]}

        C = legend['constants'].copy()
        S = organism['state'][exp_cond_name].copy()

        for i in range(len(genes_current)):
            g_name = genes_current.index.get_level_values('g_name').to_list()[i]

            if g_name in legend['constants'].index:
                for ecn in ['common', exp_cond_name]:
                    if g_name in genes_dict[ecn]:
                        if genes_dict[ecn][g_name]['is_multiplier']:
                            C[g_name] *= genes_current[ecn, g_name]
                        else:
                            C[g_name] = genes_current[ecn, g_name]

            if g_name in legend['states'].index:
                for ecn in ['common', exp_cond_name]:
                    if g_name in genes_dict[ecn]:
                        if genes_dict[ecn][g_name]['is_multiplier']:
                            S[g_name] *= genes_current[ecn, g_name]
                        else:
                            S[g_name] = genes_current[ecn, g_name]

        for c_name, c in constants_dict_current.items():
            if c_name in legend['constants'].index:
                C[c_name] = c
            if c_name in legend['states'].index:
                S[c_name] = c

        status, res = run_model_ctypes(S, C, config)
        if (status != 2) or np.any(np.isnan(res)):
            return 1

        organism['phenotype'][exp_cond_name] = pd.DataFrame(res.T, columns=legend['states'].index)
        organism['state'][exp_cond_name] = organism['phenotype'][exp_cond_name].iloc[-1]

        for i in range(len(genes_current)):
            g_name = genes_current.index.get_level_values('g_name').to_list()[i]

            if g_name in legend['states'].index:
                for ecn in ['common', exp_cond_name]:
                    if g_name in genes_dict[ecn]:
                        if genes_dict[ecn][g_name]['is_multiplier']:
                            organism['genes'][ecn, g_name] = organism['state'][exp_cond_name][g_name]\
                                                             / legend['states'][g_name]
                        else:
                            organism['genes'][ecn, g_name] = organism['state'][exp_cond_name][g_name]

    return 0


def update_fitness(organism, config):
    loss = 0
    for exp_cond_name, exp_cond in config['experimental_conditions'].items():
        if exp_cond_name == 'common':
            continue

        # phenotype_control = exp_cond['phenotype']
        # phenotype_model = organism['phenotype'][exp_cond_name]
        # V = phenotype_model['V'].values[:len(phenotype_control)]
        #
        # volumes = legend['constants'][['Vss'] + [f'Vnonjunct{i}' for i in range(1, 4 + 1)]]
        # concentrations = phenotype_model[['Cass'] + [f'Cai{i}' for i in range(1, 4 + 1)]]
        # Cai_mean = (concentrations.values * volumes.values).sum(axis=1) / sum(volumes)
        # Cai_mean = Cai_mean[:len(phenotype_control)]
        #
        # phenotype_model = np.vstack([V, Cai_mean]).T

        phenotype_control = exp_cond['phenotype']['V']
        phenotype_model = organism['phenotype'][exp_cond_name]
        phenotype_model = phenotype_model['V'][:len(phenotype_control)]

        if config['loss'] == 'RMSE':
            loss += sklearn.metrics.mean_squared_error(phenotype_control, phenotype_model, squared=False)
        elif config['loss'] == 'MSE':
            loss += sklearn.metrics.mean_squared_error(phenotype_control, phenotype_model, squared=True)
        else:
            # loss += calculate_RMSE(phenotype_control, phenotype_model)
            loss += calculate_RMSE_balanced(phenotype_control, phenotype_model)
            # loss += calculate_composite_RMSE_V_CaT(phenotype_control, phenotype_model)
            # loss += calculate_RMSE(V, phenotype_control[:, 0])

    organism['fitness'] = -loss


def generate_bounds_gammas_mask_multipliers(genes_dict):

    bounds, gammas, mask_multipliers = map(np.array, [[x[item] for y in genes_dict.values() for x in y.values()]
                                                      for item in ['bounds', 'gamma', 'is_multiplier']])
    return bounds, gammas, mask_multipliers


def save_epoch(population, kw_output):

    dump_filename = kw_output['dump_filename']
    dump_last_filename = kw_output['dump_last_filename']
    backup_filename = kw_output['backup_filename']

    timer.start('output_sort')
    population = sorted(population, key=lambda organism: organism['fitness'], reverse=True)
    timer.end('output_sort')

    timer.start('output_prepare')
    genes = np.array([organism['genes'] for organism in population])
    fitness = np.array([organism['fitness'] for organism in population])

    dump_current = np.hstack([genes, fitness[:, None]])
    timer.end('output_prepare')

    timer.start('output_dump')
    with open(dump_filename, 'ba+') as file_dump:
        dump_current.astype(np.half).tofile(file_dump)
    timer.end('output_dump')

    np.save(dump_last_filename, dump_current)

    timer.start('output_backup')
    with open(backup_filename, "wb") as file_backup:
        pickle.dump(population, file_backup)
    timer.end('output_backup')

    print(np.array2string(dump_current[0, :],
                          formatter={'float_kind': lambda x: "%.3f" % x},
                          max_line_width=1000)[1:-1], flush=True)


@njit
def transform_genes_bounds(genes, bounds, gammas, mask_multipliers):
    assert (len(genes) == len(bounds) == len(gammas) == len(mask_multipliers))

    genes_transformed = np.zeros_like(genes)
    bounds_transformed = np.zeros_like(bounds)

    scaler_dimensional = 1 / np.sqrt(len(genes))
    for i in range(len(genes)):
        lb, ub = bounds[i]
        gene = genes[i]
        if mask_multipliers[i]:  # log10 scale
            bounds_transformed[i, 1] = np.log10(ub / lb) * 1 / (gammas[i] / scaler_dimensional)
            genes_transformed[i] = np.log10(gene)
            lb_temp = np.log10(lb)
            ub_temp = np.log10(ub)
        else:  # linear scale
            genes_transformed[i] = gene
            bounds_transformed[i, 1] = (ub - lb) * 1 / (gammas[i] / scaler_dimensional)
            lb_temp = lb
            ub_temp = ub
        genes_transformed[i] = (genes_transformed[i] - lb_temp) / (ub_temp - lb_temp) * bounds_transformed[i, 1]

    return genes_transformed, bounds_transformed


@njit
def transform_genes_bounds_back(genes_transformed, bounds_transformed, bounds_back, mask_multipliers):
    assert (len(genes_transformed) == len(bounds_transformed) == len(mask_multipliers))

    genes_back = np.zeros_like(genes_transformed)

    for i in range(len(genes_transformed)):  #log10 scale
        lb_back, ub_back = bounds_back[i]
        lb_tran, ub_tran = bounds_transformed[i]
        gene = genes_transformed[i]
        if mask_multipliers[i]:
            genes_back[i] = np.log10(lb_back) + (gene - lb_tran) / (ub_tran - lb_tran) * (
                        np.log10(ub_back) - np.log10(lb_back))
            genes_back[i] = np.power(10, genes_back[i])
        else:# linear scale
            genes_back[i] = lb_back + (gene - lb_tran) / (ub_tran - lb_tran) * (ub_back - lb_back)

    return genes_back

#### ##    ## #### ########
 ##  ###   ##  ##     ##
 ##  ####  ##  ##     ##
 ##  ## ## ##  ##     ##
 ##  ##  ####  ##     ##
 ##  ##   ###  ##     ##
#### ##    ## ####    ##

comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

if len(sys.argv) < 2:
    if comm_rank == 0:
        print(f"Usage: mpiexec -n 2 python {sys.argv[0]} config.json")
    exit()

config_filename = sys.argv[1]
config_path = os.path.dirname(os.path.realpath(config_filename))

with open(config_filename) as f:
    text = f.read()
    text = strip_comments(text)
    config = json.loads(text)

config['runtime'] = dict()

config['runtime']['genes_dict'] = create_genes_dict_from_config(config)
config['runtime']['constants_dict'] = create_constants_dict_from_config(config)

m_index_tuples = [(exp_cond_name, gene_name) for exp_cond_name, gene in config['runtime']['genes_dict'].items() for gene_name in gene]
m_index = pd.MultiIndex.from_tuples(m_index_tuples)
m_index.names = ['ec_name', 'g_name']

config['runtime']['m_index'] = m_index


 ######  ######## ##    ## ########  ########  ######
##    ##    ##     ##  ##  ##     ## ##       ##    ##
##          ##      ####   ##     ## ##       ##
##          ##       ##    ########  ######    ######
##          ##       ##    ##        ##             ##
##    ##    ##       ##    ##        ##       ##    ##
 ######     ##       ##    ##        ########  ######

filename_so = os.path.join(config_path, config["filename_so"])
# print(filename_so)
filename_so_abs = os.path.abspath(filename_so)

model = ctypes.CDLL(filename_so_abs)

# model.initialize_states_default.argtypes = [
#     np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
# ]
#
# model.initialize_states_default.restype = ctypes.c_void_p
#
#
# model.initialize_constants_default.argtypes = [
#     np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
# ]
#
# model.initialize_constants_default.restype = ctypes.c_void_p


model.run.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_double,
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')
]

model.run.restype = ctypes.c_int
run_model_ctypes.model = model

###############################################################################


legend = dict()
legend['states'] = pd.read_csv(os.path.join(config_path, config["filename_legend_states"]),
                               usecols=['name', 'value'], index_col='name')['value']  # Series
legend['constants'] = pd.read_csv(os.path.join(config_path, config["filename_legend_constants"]),
                                  usecols=['name', 'value'], index_col='name')['value']  # Series
config['runtime']['legend'] = legend

if config['n_organisms'] % comm_size != 0:
    config['runtime']['n_organisms'] = int(np.ceil(config['n_organisms'] / comm_size) * comm_size)
    if comm_rank == 0:
        print(f'# `n_organisms` is changed from {config["n_organisms"]} to {config["runtime"]["n_organisms"]}', flush=True)
else:
    config['runtime']['n_organisms'] = config['n_organisms']


for exp_cond_name, exp_cond in config['experimental_conditions'].items():

    if exp_cond_name == 'common':
        continue

    filename_phenotype = os.path.normpath(os.path.join(config_path, exp_cond['filename_phenotype']))
    # exp_cond['phenotype'] = np.loadtxt(filename_phenotype)
    exp_cond['phenotype'] = pd.read_csv(filename_phenotype)
    exp_cond['filename_phenotype'] = filename_phenotype

    filename_state = os.path.normpath(os.path.join(config_path, exp_cond['filename_state']))
    exp_cond['initial_state'] = pd.Series(np.loadtxt(filename_state), index=legend['states'].index)
    exp_cond['filename_state'] = filename_state

states_initial = pd.DataFrame(data={exp_cond_name: exp_cond['initial_state'].copy()
                                    for exp_cond_name, exp_cond in config['experimental_conditions'].items()
                                    if exp_cond_name != 'common'})

config['runtime']['states_initial'] = states_initial


output_folder_name = os.path.normpath(os.path.join(config_path, config.get("output_folder_name", "./output")))
if comm_rank == 0:
    time_suffix = datetime.now().strftime("%y%m%d_%H%M%S")
    output_folder_name = os.path.join(output_folder_name, time_suffix)
    print(f"# output_folder_name was set to {output_folder_name}", flush=True)
    config['runtime']['time_suffix'] = time_suffix
else:
    output_folder_name = None
output_folder_name = comm.bcast(output_folder_name, root=0)  # due to time call
config['runtime']["output_folder_name"] = output_folder_name
comm.Barrier()

config['runtime']['output'] = dict(output_folder_name_phenotype = os.path.join(output_folder_name, "phenotype"),
                                   dump_filename = os.path.join(output_folder_name, "dump.bin"),
                                   dump_last_filename = os.path.join(output_folder_name, "dump_last.npy"),
                                   backup_filename = os.path.join(output_folder_name, "backup.pickle"),
                                   output_folder_name_mpi = os.path.join(output_folder_name, "mpi"),
                                   config_backup_filename = os.path.join(output_folder_name, "config_backup.pickle"))

if comm_rank == 0:
    for folder in config['runtime']['output']['output_folder_name_phenotype'],\
                  config['runtime']['output']['output_folder_name_mpi']:
        os.makedirs(folder, exist_ok=True)
    with open(config['runtime']['output']['dump_filename'], "wb") as file_dump:  # create or clear and close
        pass
    with open(config['runtime']['output']['config_backup_filename'], "wb") as file_config_backup:
        pickle.dump(config, file_config_backup)
    print("# SIZE =", comm_size, flush=True)

time_start = time.time()

bounds, gammas, mask_multipliers = generate_bounds_gammas_mask_multipliers(config['runtime']['genes_dict'])
config['runtime']['bounds'] = bounds
config['runtime']['gammas'] = gammas
config['runtime']['mask_multipliers'] = mask_multipliers

config['runtime']['kw_ga'] = dict(crossover_rate=config.get('crossover_rate', 1.0),
                                  mutation_rate=config.get('mutation_rate', 0.1),
                                  gamma=config.get('gamma', 1.0))

n_orgsnisms_per_process = config['runtime']['n_organisms'] // comm_size

# organism_dummy = generate_organism(genes_dict, states_initial)
# update_phenotype_state(organism_dummy, legend, config)

genes_size = sum(map(len, config['runtime']['genes_dict'].values()))

state_size = states_initial.size
state_shape = states_initial.shape

# phenotype_size = np.concatenate(organism_dummy['phenotype']).size
# phenotype_lens = [len(phenotype) for phenotype in organism_dummy['phenotype']]
# phenotype_split_indices = np.cumsum(phenotype_lens)[:-1]

recvbuf_genes = np.empty([comm_size, n_orgsnisms_per_process * genes_size])
recvbuf_state = np.empty([comm_size, n_orgsnisms_per_process * state_size])
# recvbuf_phenotype = np.empty([comm_size, n_orgsnisms_per_process * phenotype_size])
recvbuf_fitness = np.empty([comm_size, n_orgsnisms_per_process * 1])

# del organism_dummy

#  Main cycle
if comm_rank == 0:
    initial_population_filename = config.get('initial_population_filename', None)
    if initial_population_filename:
        initial_population_filename = os.path.normpath(os.path.join(config_path, initial_population_filename))
        with open(initial_population_filename, 'rb') as f:
            backup = pickle.load(f)
        population = init_population_from_backup(backup, config)
        config['runtime']['initial_population_filename'] = initial_population_filename
        print(f"population was loaded from {initial_population_filename}", flush=True)
    else:
        population = init_population(config)
    #  population[0]['genes'] = np.ones_like(population[0]['genes'])
    #  population[-1]['genes'][-1] = 200
    population = batches_from_list(population, comm_size)
else:
    population = None

batch = comm.scatter(population, root=0)

timer = Timer()

##     ##    ###    #### ##    ##
###   ###   ## ##    ##  ###   ##
#### ####  ##   ##   ##  ####  ##
## ### ## ##     ##  ##  ## ## ##
##     ## #########  ##  ##  ####
##     ## ##     ##  ##  ##   ###
##     ## ##     ## #### ##    ##

for epoch in range(config['n_generations']):

    #  calculations
    timer.start('calc')

    index_best_per_batch = 0  # index of the best organism per batch, used for memory-optimization

    for i, organism in enumerate(batch):

        status = update_phenotype_state(organism, config)

        FLAG_INVALID = False
        msg_invalid = "Invalid organism:\n"
        if status != 0:
            FLAG_INVALID = True
            msg_invalid += f"# status = {status}\n"
        elif np.any(organism['genes'] <= bounds[:, 0]) or np.any(organism['genes'] >= bounds[:, 1]):
            FLAG_INVALID = True
            s = " ".join([f"{'*' if (gene <= bounds[i, 0] or gene >= bounds[i, 1]) else ''}{gene:.3e}"
                          for i, gene in enumerate(organism['genes'])])
            msg_invalid += f"# out of bounds: {s}\n"

            # with open("invalid.pickle", "wb") as f:
            #     pickle.dump(organism, f)
            #
            # print(msg_invalid, flush=True)
            # exit()

        if FLAG_INVALID:
            print('#', msg_invalid, end='', flush=True)
            # organism['phenotype'] = [np.empty(length) for length in phenotype_lens]  # needed for gather
            organism['fitness'] = np.NINF
            del organism['phenotype']
        else:
            update_fitness(organism, config)
            if organism['fitness'] > batch[index_best_per_batch]['fitness']:
                if 'phenotype' in batch[index_best_per_batch]:  # best guy could be invalid
                    del batch[index_best_per_batch]['phenotype']
                index_best_per_batch = i
                # print(f"{comm_rank} best is {index_best_per_batch}", flush=True)
            elif i:
                del organism['phenotype']
            #organism['genes'], bounds_transformed = transform_genes_bounds(organism['genes'], bounds, gammas,
            #                                                               n_multipliers=len(config['multipliers']))
    timer.end('calc')

     ######      ###    ######## ##     ## ######## ########
    ##    ##    ## ##      ##    ##     ## ##       ##     ##
    ##         ##   ##     ##    ##     ## ##       ##     ##
    ##   #### ##     ##    ##    ######### ######   ########
    ##    ##  #########    ##    ##     ## ##       ##   ##
    ##    ##  ##     ##    ##    ##     ## ##       ##    ##
     ######   ##     ##    ##    ##     ## ######## ##     ##

    timer.start('gather_sendbuf')
    sendbuf_genes = np.concatenate([organism['genes'] for organism in batch])
    sendbuf_state = np.concatenate([organism['state'].values.flatten() for organism in batch])
    #sendbuf_phenotype = np.concatenate([np.concatenate(organism['phenotype']) for organism in batch])
    sendbuf_fitness = np.array([organism['fitness'] for organism in batch])
    assert(not np.any(np.isnan(sendbuf_fitness)))
    timer.end('gather_sendbuf')

    timer.start('gather_allgather')
    comm.Allgatherv(sendbuf_genes, recvbuf_genes)
    comm.Allgatherv(sendbuf_state, recvbuf_state)
    #comm.Allgatherv(sendbuf_phenotype, recvbuf_phenotype)
    comm.Allgatherv(sendbuf_fitness, recvbuf_fitness)
    timer.end('gather_allgather')

    timer.start('gather_recvbuf')
    recvbuf_genes = recvbuf_genes.reshape((config['runtime']['n_organisms'], genes_size))
    recvbuf_state = recvbuf_state.reshape((config['runtime']['n_organisms'], *state_shape))
    #recvbuf_phenotype = recvbuf_phenotype.reshape((config['n_organisms'], phenotype_size))
    recvbuf_fitness = recvbuf_fitness.flatten()
    timer.end('gather_recvbuf')

    timer.start('gather_population')
    #  assert (not np.any(np.isnan(recvbuf_fitness)))

    population = [dict(genes     = recvbuf_genes[i],  # pd.Series(data=recvbuf_genes[i], index=m_index),
                       state     = recvbuf_state[i],  # pd.DataFrame(recvbuf_state[i], columns=[exp_cond_name for exp_cond_name in config['experimental_conditions'] if exp_cond_name != 'common']),
                       # phenotype = np.split(recvbuf_phenotype[i].copy(), indices_or_sections=phenotype_split_indices),
                       fitness   = recvbuf_fitness[i]) for i in range(config['runtime']['n_organisms'])]
    timer.end('gather_population')

     ######     ###    ##     ## ########
    ##    ##   ## ##   ##     ## ##
    ##        ##   ##  ##     ## ##
     ######  ##     ## ##     ## ######
          ## #########  ##   ##  ##
    ##    ## ##     ##   ## ##   ##
     ######  ##     ##    ###    ########

    index_best = argmax_list_of_dicts(population, 'fitness')

    comm_rank_best = index_best // n_orgsnisms_per_process
    index_best_batch = index_best % n_orgsnisms_per_process

    timer.start('save_phenotype')
    if comm_rank == comm_rank_best:
        # print(f"rank {comm_rank_best}: index_best = {index_best}, {index_best_batch}", flush=True)
        organism_best = batch[index_best_batch]
        for exp_cond_name in config['experimental_conditions']:
            if exp_cond_name == 'common':
                continue

            df = organism_best['phenotype'][exp_cond_name]

            # Rewrite last epoch
            filename_phenotype_save = os.path.join(config['runtime']['output']['output_folder_name_phenotype'],
                                                   f"phenotype_{exp_cond_name}.csv")
            df.to_csv(filename_phenotype_save, index=False)

            # Append last epoch to previous
            filename_phenotype_save_binary = os.path.join(config['runtime']['output']['output_folder_name_phenotype'],
                                                          f"phenotype_{exp_cond_name}.bin")
            with open(filename_phenotype_save_binary, 'ba+' if epoch else 'bw') as f:
                df.values.astype(np.float32).tofile(f)

    timer.end('save_phenotype')

    timer.start('output_sort')
    timer.end('output_sort')
    timer.start('output_prepare')
    timer.end('output_prepare')
    timer.start('output_dump')
    timer.end('output_dump')
    timer.start('output_backup')
    timer.end('output_backup')

    if comm_rank == epoch % comm_size:
        save_epoch(population, config['runtime']['output'])

     ######   ######## ##    ## ######## ######## ####  ######
    ##    ##  ##       ###   ## ##          ##     ##  ##    ##
    ##        ##       ####  ## ##          ##     ##  ##
    ##   #### ######   ## ## ## ######      ##     ##  ##
    ##    ##  ##       ##  #### ##          ##     ##  ##
    ##    ##  ##       ##   ### ##          ##     ##  ##    ##
     ######   ######## ##    ## ########    ##    ####  ######

    timer.start('gene')

    if len(population) <= 3:
        if comm_rank == 0:
            print(f"# Not enough organisms for genetic operations left: {len(population)}\nexit", flush=True)
        comm.Barrier()
        exit()

    population.sort(key=lambda organism: organism['fitness'], reverse=True)
    index_first_invalid = find_index_first(population, lambda organism: organism['fitness'] == np.NINF)
    if index_first_invalid:
        if comm_rank == 0:
            n_invalids = len(population) - index_first_invalid
            percentage_invalid = (n_invalids) / len(population) * 100
            print(f"# {n_invalids} ({percentage_invalid:.2f} %) invalids were deleted", flush=True)
        population = population[:index_first_invalid]
    elites_all = population[:config['n_elites']]  # len may be less than config['n_elites'] due to invalids
    elites_batch = elites_all[comm_rank::comm_size]
    n_elites = len(elites_batch)
    #  print(f"# {comm_rank} has {n_elites} elites", flush=True)
    #  elites_batch may be empty list

    for organism in population:
        organism['genes'], bounds_transformed = transform_genes_bounds(organism['genes'],
                                                                       bounds, gammas, mask_multipliers)
    batch = do_step(population, new_size=n_orgsnisms_per_process - n_elites,
                    elite_size=0, bounds=bounds_transformed, **config['runtime']['kw_ga'])

    batch += elites_batch

    assert (len(batch) == n_orgsnisms_per_process)

    for organism in batch:
        organism['genes'] = transform_genes_bounds_back(organism['genes'],
                                                        bounds_transformed, bounds_back=bounds,
                                                        mask_multipliers=mask_multipliers)

        organism['genes'] = pd.Series(data=organism['genes'], index=config['runtime']['m_index'])
        organism['state'] = pd.DataFrame(organism['state'],
                                         columns=[exp_cond_name for exp_cond_name in config['experimental_conditions']
                                                  if exp_cond_name != 'common'])
        organism['state'].index = states_initial.index

    timer.end('gene')

    ########  ######## ########   #######  ########  ########
    ##     ## ##       ##     ## ##     ## ##     ##    ##
    ##     ## ##       ##     ## ##     ## ##     ##    ##
    ########  ######   ########  ##     ## ########     ##
    ##   ##   ##       ##        ##     ## ##   ##      ##
    ##    ##  ##       ##        ##     ## ##    ##     ##
    ##     ## ######## ##         #######  ##     ##    ##

    if comm_rank == epoch % comm_size:
        print(f"# EPOCH {epoch}:")
        print(timer.report(sort=True), flush=True)

    filename_mpi_report = os.path.join(config['runtime']['output']['output_folder_name_mpi'],
                                       f"report_{comm_rank:04d}.csv")
    mode = 'a' if os.path.isfile(filename_mpi_report) and epoch != 0 else 'w'
    with open(filename_mpi_report, mode) as f:
        if mode == 'w':
            header = ','.join(timer.times.keys())
            f.write(header + '\n')
        f.write(','.join([str(v) for v in timer.times.values()]) + '\n')


if comm_rank == 0:
    time_end = time.time()
    print("# TIME =", time_end - time_start, flush=True)
    file_dump.close()
