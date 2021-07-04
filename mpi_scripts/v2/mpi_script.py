from mpi4py import MPI
import time
import pickle

import numpy as np

from pypoptim.helpers import find_index_first

from pypoptim.model import CardiacModel
from solmodel import SolModel
from pypoptim.algorythm.ga import GA

from io_utils import prepare_config, save_epoch
from mpi_utils import allocate_buffers, allgather

import sys

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

if comm_rank == 0:
    config_filename = sys.argv[1]
    config = prepare_config(config_filename)
else:
    config = None
config = comm.bcast(config, root=0)

if config['n_organisms'] % comm_size != 0:
    config['runtime']['n_organisms'] = int(np.ceil(config['n_organisms'] / comm_size) * comm_size)
    if comm_rank == 0:
        print(f'# `n_organisms` is changed from {config["n_organisms"]} to {config["runtime"]["n_organisms"]}',
              flush=True)
else:
    config['runtime']['n_organisms'] = config['n_organisms']


model = CardiacModel(config['runtime']['filename_so'])
SolModel.model = model
SolModel.config = config

ga_optim = GA(SolModel,
              config['runtime']['bounds'],
              config['runtime']['gammas'],
              config['runtime']['mask_multipliers'],
              keys_data_transmit=['state'])

recvbuf_genes, recvbuf_state, recvbuf_fitness = allocate_buffers(config, comm)

initial_population_filename = config.get('initial_population_filename', None)
if initial_population_filename is not None:
    raise NotImplementedError

batch = ga_optim.generate_population(config['runtime']['n_orgsnisms_per_process'])

if comm_rank == 0:
    with open(config['runtime']['output']['config_backup_filename'], "wb") as file_config_backup:
        pickle.dump(config, file_config_backup)



# timer = Timer()

##     ##    ###    #### ##    ##
###   ###   ## ##    ##  ###   ##
#### ####  ##   ##   ##  ####  ##
## ### ## ##     ##  ##  ## ## ##
##     ## #########  ##  ##  ####
##     ## ##     ##  ##  ##   ###
##     ## ##     ## #### ##    ##

for epoch in range(config['n_generations']):

    index_best_per_batch = 0  # index of the best organism per batch, used for memory-optimization

    for i, sol in enumerate(batch):
        sol.update()
        if not sol.is_valid():
            sol._y = np.inf

    batch = sorted(batch)

    population = allgather(batch, config, SolModel, recvbuf_genes, recvbuf_state, recvbuf_fitness, comm)

     ######     ###    ##     ## ########
    ##    ##   ## ##   ##     ## ##
    ##        ##   ##  ##     ## ##
     ######  ##     ## ##     ## ######
          ## #########  ##   ##  ##
    ##    ## ##     ##   ## ##   ##
     ######  ##     ##    ###    ########
    #
    # index_best = argmax_list_of_dicts(population, 'fitness')
    #
    # comm_rank_best = index_best // n_orgsnisms_per_process
    # index_best_batch = index_best % n_orgsnisms_per_process
    #
    # timer.start('save_phenotype')
    # if comm_rank == comm_rank_best:
    #
    #     organism_best = batch[index_best_batch]
    #
    #     with open(config['runtime']['output']['organism_best_filename'], 'bw') as f:
    #         pickle.dump(organism_best, f)
    #
    #     organism_best['genes'].to_csv(config['runtime']['output']['genes_best_filename'])
    #
    #     for exp_cond_name in config['experimental_conditions']:
    #         if exp_cond_name == 'common':
    #             continue
    #
    #         df = organism_best['phenotype'][exp_cond_name]
    #
    #         # Rewrite last epoch
    #         filename_phenotype_save = os.path.join(config['runtime']['output']['output_folder_name_phenotype'],
    #                                                f"phenotype_{exp_cond_name}.csv")
    #         df.to_csv(filename_phenotype_save, index=False)
    #
    #         # Append last epoch to previous
    #         filename_phenotype_save_binary = os.path.join(config['runtime']['output']['output_folder_name_phenotype'],
    #                                                       f"phenotype_{exp_cond_name}.bin")
    #         with open(filename_phenotype_save_binary, 'ba+' if epoch else 'bw') as f:
    #             df.values.astype(np.float32).tofile(f)
    #
    # timer.end('save_phenotype')
    #
    # if comm_rank == epoch % comm_size:
    #     save_epoch(population, config['runtime']['output'])


    # timer.start('gene')

    if len(population) <= 3:
        if comm_rank == 0:
            with open(config['runtime']['output']['log_filename'], "a") as file_log:
                file_log.write(f"# Not enough organisms for genetic operations left: {len(population)}\nexit\n")
        exit()

    population.sort()
    index_first_invalid = find_index_first(population, lambda sol: np.isinf(sol.y))
    if index_first_invalid:
        if comm_rank == 0:
            n_invalids = len(population) - index_first_invalid
            percentage_invalid = n_invalids / len(population) * 100
            with open(config['runtime']['output']['log_filename'], "a") as file_log:
                file_log.write(f"# {n_invalids} ({percentage_invalid:.2f} %) invalids were deleted\n")
        population = population[:index_first_invalid]
    elites_all = population[:config['n_elites']]  # len may be less than config['n_elites'] due to invalids
    elites_batch = elites_all[comm_rank::comm_size]  # elites_batch may be empty
    n_elites = len(elites_batch)
    n_mutants = config['runtime']['n_orgsnisms_per_process'] - n_elites

    mutants_batch = ga_optim.get_mutants(population, n_mutants)
    batch = elites_batch + mutants_batch

    assert (len(batch) == config['runtime']['n_orgsnisms_per_process'])


    # timer.end('gene')

    #  TODO: report