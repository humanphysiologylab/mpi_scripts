import os
import json
import git
import pickle
from datetime import datetime
import numpy as np
import pandas as pd

from pypoptim.helpers import strip_comments
from gene_utils import create_genes_dict_from_config, \
                       create_constants_dict_from_config, \
                       generate_bounds_gammas_mask_multipliers


def prepare_config(config_filename):

    config_path = os.path.dirname(os.path.realpath(config_filename))

    with open(config_filename) as f:
        text = f.read()
        text = strip_comments(text)
        config = json.loads(text)

    config['runtime'] = dict()

    config['runtime']['config_path'] = config_path

    config['runtime']['sha'] = git.Repo(search_parent_directories=True).head.commit.hexsha

    config['runtime']['genes_dict'] = create_genes_dict_from_config(config)
    config['runtime']['constants_dict'] = create_constants_dict_from_config(config)

    m_index_tuples = [(exp_cond_name, gene_name) for exp_cond_name, gene in config['runtime']['genes_dict'].items() for
                      gene_name in gene]
    m_index = pd.MultiIndex.from_tuples(m_index_tuples)
    m_index.names = ['ec_name', 'g_name']

    config['runtime']['m_index'] = m_index

    legend = dict()
    legend['states'] = pd.read_csv(os.path.join(config_path, config["filename_legend_states"]),
                                   usecols=['name', 'value'], index_col='name')['value']  # Series
    legend['constants'] = pd.read_csv(os.path.join(config_path, config["filename_legend_constants"]),
                                      usecols=['name', 'value'], index_col='name')['value']  # Series
    config['runtime']['legend'] = legend

    for exp_cond_name, exp_cond in config['experimental_conditions'].items():

        if exp_cond_name == 'common':
            continue

        filename_phenotype = os.path.normpath(os.path.join(config_path, exp_cond['filename_phenotype']))
        exp_cond['phenotype'] = pd.read_csv(filename_phenotype)
        exp_cond['filename_phenotype'] = filename_phenotype

        filename_state = os.path.normpath(os.path.join(config_path, exp_cond['filename_state']))
        exp_cond['initial_state'] = pd.Series(np.loadtxt(filename_state), index=legend['states'].index)
        exp_cond['filename_state'] = filename_state

        column_stim_protocol = config.get('column_stim_protocol', None)
        if column_stim_protocol is not None:
            filename_stim_protocol = os.path.normpath(os.path.join(config_path, exp_cond['filename_stim_protocol']))
            exp_cond['stim_protocol'] = pd.read_csv(
                filename_stim_protocol)  # [column_stim_protocol]  # pd.Series is returned
            exp_cond['filename_stim_protocol'] = filename_stim_protocol
        else:
            exp_cond['stim_protocol'] = None

    states_initial = pd.DataFrame(data={exp_cond_name: exp_cond['initial_state'].copy()
                                        for exp_cond_name, exp_cond in config['experimental_conditions'].items()
                                        if exp_cond_name != 'common'})

    config['runtime']['states_initial'] = states_initial

    output_folder_name = os.path.normpath(os.path.join(config_path, config.get("output_folder_name", "./output")))
    time_suffix = datetime.now().strftime("%y%m%d_%H%M%S")
    output_folder_name = os.path.join(output_folder_name, time_suffix)
    # print(f"# output_folder_name was set to {output_folder_name}", flush=True)
    config['runtime']['time_suffix'] = time_suffix
    config['runtime']["output_folder_name"] = output_folder_name

    config['runtime']['output'] = dict(output_folder_name_phenotype=os.path.join(output_folder_name, "phenotype"),
                                       dump_filename=os.path.join(output_folder_name, "dump.bin"),
                                       dump_last_filename=os.path.join(output_folder_name, "dump_last.npy"),
                                       dump_elite_filename=os.path.join(output_folder_name, "dump_elite.npy"),
                                       backup_filename=os.path.join(output_folder_name, "backup.pickle"),
                                       output_folder_name_mpi=os.path.join(output_folder_name, "mpi"),
                                       config_backup_filename=os.path.join(output_folder_name, "config_backup.pickle"),
                                       log_filename=os.path.join(output_folder_name, "runtime.log"),
                                       organism_best_filename=os.path.join(output_folder_name, "organism_best.pickle"),
                                       genes_best_filename=os.path.join(output_folder_name, "genes_best.csv"),
                                       phenotypes_plot_filename=os.path.join(output_folder_name,
                                                                             "phenotypes_plot_filename.png"))


    bounds, gammas, mask_multipliers = generate_bounds_gammas_mask_multipliers(config['runtime']['genes_dict'])
    config['runtime']['bounds'] = bounds
    config['runtime']['gammas'] = gammas
    config['runtime']['mask_multipliers'] = mask_multipliers

    config['runtime']['kw_ga'] = dict(crossover_rate=config.get('crossover_rate', 1.0),
                                      mutation_rate=config.get('mutation_rate', 0.1),
                                      gamma=config.get('gamma', 1.0))

    return config


def touch_output(config):
    for folder in config['runtime']['output']['output_folder_name_phenotype'],:
        os.makedirs(folder, exist_ok=True)
    with open(config['runtime']['output']['dump_filename'], "wb") as f:  # create or clear and close
        pass
    with open(config['runtime']['output']['dump_elite_filename'], "wb") as f:  # create or clear and close
        pass
    with open(config['runtime']['output']['log_filename'], "w") as file_log:
        file_log.write(f"# SIZE = {config['runtime']['comm_size']}\n")
        file_log.write(f"# commit {config['runtime']['sha']}\n")


def save_epoch(population, kw_output):

    dump_filename = kw_output['dump_filename']
    dump_last_filename = kw_output['dump_last_filename']
    dump_elite_filename = kw_output['dump_elite_filename']
    backup_filename = kw_output['backup_filename']
    log_filename = kw_output['log_filename']

    population = sorted(population, key=lambda organism: organism['fitness'], reverse=True)

    genes = np.array([organism['genes'] for organism in population])  # TODO: take this from buffers
    fitness = np.array([organism['fitness'] for organism in population])

    dump_current = np.hstack([genes, fitness[:, None]])

    with open(dump_filename, 'ba+') as file_dump:
        dump_current.astype(np.half).tofile(file_dump)

    np.save(dump_last_filename, dump_current)

    with open(backup_filename, "wb") as file_backup:
        pickle.dump(population, file_backup)

    with open(log_filename, "a") as file_log:
        s = np.array2string(dump_current[0, :],
                            formatter={'float_kind': lambda x: "%.3f" % x},
                            max_line_width=1000)[1:-1] + "\n"
        file_log.write(s)

    with open(dump_elite_filename, 'ba+') as f:
        dump_current[0, :].tofile(f)
