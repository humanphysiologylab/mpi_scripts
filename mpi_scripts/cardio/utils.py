import os
import pickle

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from numba import njit

from scipy.integrate import solve_ivp

import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

from src.helpers import MSE, RMSE, calculate_mean_abs_noise, calculate_RMSE_weightened, \
                        calculate_RMSE_balanced, calculate_composite_RMSE_V_CaT, autoscaling, \
                        value_from_bounds

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
    else:
        print('Invalid config: check n_beats & t_run entries.',
              file=sys.stderr, fflush=True)
        exit()

    tol = config.get('tol', 1e-6)

    n_samples_per_stim = int(stim_period / t_sampling)
    output = np.empty((n_samples_per_stim * n_beats + 1, len(S)))

    status = run_model_ctypes.model.run(S.values.copy(), C.values.copy(),
                                        n_beats, t_sampling, tol, output)

    output = output[-n_samples_per_stim - 1:].T  # last beat

    return status, output


def run_model_scipy(S, C, config):

        stim_period = C[config['stim_period_legend_name']]
        t_sampling = config['t_sampling']

        if 'n_beats' in config and 't_run' not in config:
            n_beats = config['n_beats']
        elif 'n_beats' not in config and 't_run' in config:
            t_run = config['t_run']
            n_beats = np.ceil(t_run / stim_period).astype(int)
            if n_beats % 2 == 0:
                n_beats += 1
        else:
            print('Invalid config: check n_beats & t_run entries.',
                  file=sys.stderr, fflush=True)
            exit()

        #tol = config.get('tol', 1e-6)

        n_samples_per_stim = int(stim_period / t_sampling)

        t_space = np.linspace(0, stim_period * n_beats, n_samples_per_stim * n_beats + 1, endpoint=True)
        t_span = 0, t_space[-1]

        def fun(t, y, params): # TODO: rewrite
            ydot = np.zeros_like(y)
            run_model_scipy.model.fun(t, y, ydot, params)
            return ydot

        sol = solve_ivp(fun, y0=S,
                        t_span=t_span, t_eval=t_space,
                        args=(C.values.copy(),),
                        method='LSODA',# rtol=1e-9,
                        max_step=1. * t_sampling,
                )

        #output = output[-n_samples_per_stim - 1:].T  # last beat
        output = sol.y[:, -n_samples_per_stim - 1:]
        status = sol.status

        return status, output


def update_phenotype_state(organism, config):

    organism['phenotype'] = dict()

    legend = config['runtime']['legend']
    genes_dict = config['runtime']['genes_dict']
    constants_dict = config['runtime']['constants_dict']

    for exp_cond_name in config['experimental_conditions']:

        if exp_cond_name == 'common':
            continue

        if exp_cond_name in organism['genes']:
            genes_current = organism['genes'][['common',
                                               exp_cond_name]]
        else:
            genes_current = organism['genes'][['common']]

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

        if config['use_scipy']:
            status, res = run_model_scipy(S, C, config)
            if (status != 0) or np.any(np.isnan(res)):
                return 1
        else:
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

        if 'koivumaki' in config['filename_so']:

            phenotype_control = exp_cond['phenotype'].values.copy()
            #phenotype_control[:, 1] *= 740 / 864  # vanessa Kd

            phenotype_model = organism['phenotype'][exp_cond_name]
            V = phenotype_model['V'].values[:len(phenotype_control)]

            legend = config['runtime']['legend']
            volumes = legend['constants'][['Vss'] + [f'Vnonjunct{i}' for i in range(1, 4 + 1)]]

            # concentrations = phenotype_model[['fluo_ss'] + [f'fluo_{i}' for i in range(1, 4 + 1)]]
            concentrations = phenotype_model[['Cass'] + [f'Cai{i}' for i in range(1, 4 + 1)]]

            Cai_mean = (concentrations.values * volumes.values).sum(axis=1) / sum(volumes)
            Cai_mean = Cai_mean[:len(phenotype_control)]

            phenotype_model = np.vstack([V, Cai_mean]).T

        else:

            phenotype_control = exp_cond['phenotype']['V']
            phenotype_model = organism['phenotype'][exp_cond_name]
            phenotype_model = phenotype_model['V'][:len(phenotype_control)]

        if config['loss'] == 'RMSE':
            loss += RMSE(phenotype_control, phenotype_model)

        elif config['loss'] == 'MSE':
            loss += MSE(phenotype_control, phenotype_model)

        elif config['loss'] == 'RMSE_balanced':
            loss += calculate_RMSE_balanced(phenotype_control, phenotype_model)

        elif config['loss'] == 'RMSE_weightened':

            phenotype_model   = (phenotype_model   - phenotype_control.min(axis=0)) / phenotype_control.ptp(axis=0)
            phenotype_control = (phenotype_control - phenotype_control.min(axis=0)) / phenotype_control.ptp(axis=0)

            weights = 1 / calculate_mean_abs_noise(phenotype_control)
            weights /= sum(weights)
            loss += calculate_RMSE_weightened(phenotype_control, phenotype_model, weights)

        elif config['loss'] == 'composite_RMSE_V_CaT':
            loss += calculate_composite_RMSE_V_CaT(phenotype_control, phenotype_model)

        elif config['loss'] == 'composite_RMSE_V_CaT_noisy':

            phenotype_model   = (phenotype_model   - phenotype_control.min(axis=0)) / phenotype_control.ptp(axis=0)
            phenotype_control = (phenotype_control - phenotype_control.min(axis=0)) / phenotype_control.ptp(axis=0)

            weights = 1 / calculate_mean_abs_noise(phenotype_control)
            weights /= sum(weights)

            rmse_v = RMSE(phenotype_control[:, 0], phenotype_model[:, 0])

            ca_exp_scaled, rmse_ca, coeffs = autoscaling(signal_to_scale=phenotype_control[:, 1],
                                                         signal_reference=phenotype_model[:, 1])

            loss += rmse_v * weights[0] + rmse_ca * weights[1]

        elif config['loss'] == 'CaT_RMSE':
            loss += RMSE(phenotype_control[:, 1],
                         phenotype_model[:, 1])

        elif config['loss'] == 'CaT_autoscaling':
            ca_exp_scaled, rmse_ca, coeffs = autoscaling(signal_to_scale=phenotype_control[:, 1],
                                                         signal_reference=phenotype_model[:, 1])
            loss += rmse_ca

        else:
            print(f'Unknown loss {config["loss"]}',
                  file=sys.stderr, fflush=True)
            exit()

    organism['fitness'] = -loss


def generate_bounds_gammas_mask_multipliers(genes_dict):

    bounds, gammas, mask_multipliers = map(np.array, [[x[item] for y in genes_dict.values() for x in y.values()]
                                                      for item in ['bounds', 'gamma', 'is_multiplier']])
    return bounds, gammas, mask_multipliers


def save_epoch(population, kw_output):

    dump_filename = kw_output['dump_filename']
    dump_last_filename = kw_output['dump_last_filename']
    backup_filename = kw_output['backup_filename']
    log_filename = kw_output['log_filename']

    #  timer.start('output_sort')
    population = sorted(population, key=lambda organism: organism['fitness'], reverse=True)
    #  timer.end('output_sort')

    #  timer.start('output_prepare')
    genes = np.array([organism['genes'] for organism in population])
    fitness = np.array([organism['fitness'] for organism in population])

    dump_current = np.hstack([genes, fitness[:, None]])
    #  timer.end('output_prepare')

    #  timer.start('output_dump')
    with open(dump_filename, 'ba+') as file_dump:
        dump_current.astype(np.half).tofile(file_dump)
    #  timer.end('output_dump')

    np.save(dump_last_filename, dump_current)

    #  timer.start('output_backup')
    with open(backup_filename, "wb") as file_backup:
        pickle.dump(population, file_backup)
    #  timer.end('output_backup')

    with open(log_filename, "a") as file_log:
        s = np.array2string(dump_current[0, :],
                            formatter={'float_kind': lambda x: "%.3f" % x},
                            max_line_width=1000)[1:-1] + "\n"
        file_log.write(s)


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


def plot_phenotypes(organism, config, filename_save):

    exp_cond_dict = [item for item in config['experimental_conditions'].items() if item[0] != 'common']
    phenotype_model_last = organism['phenotype']

    nrows = 2
    ncols = len(phenotype_model_last)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                             figsize=plt.figaspect(nrows / ncols),
                             sharex='col', sharey='row')

    for i_item, item in enumerate(['V', 'Cai']):

        for i_cond, (exp_cond_name, exp_cond) in enumerate(exp_cond_dict):

            ax = axes[i_item, i_cond]

            ax.set_title(f'CL = {exp_cond_name} ms')
            ax.set_xlabel('time, ms')

            xlim_right = 450
            ax.set_xticks(np.arange(0, xlim_right + 1, 100))
            ax.set_xlim(-20, xlim_right)

            ax.grid(True)

            exp = exp_cond['phenotype'][item]

            if item == 'Cai':

                ax_inset = None

                legend = config['runtime']['legend']
                volumes = legend['constants'][['Vss'] + [f'Vnonjunct{i}' for i in range(1, 4 + 1)]]
                concentrations = phenotype_model_last[exp_cond_name][['Cass'] + [f'Cai{i}' for i in range(1, 4 + 1)]]
                Cai_mean = (concentrations.values * volumes.values).sum(axis=1) / sum(volumes)
                Cai_mean = Cai_mean[:len(exp)]

                model = Cai_mean

            else:

                ax.set_yticks(np.arange(-100, 75 + 1, 25))
                ax.set_ylim(-100, 75)

                ax_inset = ax.inset_axes([0.5, 0.5, 0.4, 0.4])

                ax_inset.set_xlim(-1, 6)
                ax_inset.set_ylim(ax.get_ylim())
                ax_inset.set_yticks(ax.get_yticks()[1:-1])
                ax_inset.set_xticks(np.arange(-1, 6 + 1, 1))
                ax_inset.tick_params(axis='both', labelsize='x-small')

                ax_inset.grid(True)

                model = phenotype_model_last[exp_cond_name]['V'].values[:len(exp)]

            exp, model = map(lambda x: np.roll(x, 1), [exp, model])

            if i_cond == 0:
                ax.set_ylabel(item)

            for ml, ax_ in zip(['-', '.-'], [ax, ax_inset]):

                if ax_ is None:
                    continue

                ax_.plot(exp, ml, color='w', lw=3)
                ax_.plot(exp, ml, color='0.3', lw=2.5, label='exp')

                ax_.plot(model, ml, color='w', lw=2)
                ax_.plot(model, ml, color='C3', lw=1.25, label='model')

    fig.align_labels()
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)

    fig.suptitle(filename_save)

    plt.savefig(filename_save, dpi=300, facecolor='white')
