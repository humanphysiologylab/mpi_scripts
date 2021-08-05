import os
import json
import git
import pickle
from datetime import datetime
import numpy as np
import pandas as pd

from gene_utils import (
    create_genes_dict_from_config,
    create_constants_dict_from_config,
    generate_bounds_gammas_mask_multipliers,
)


def prepare_config(config_filename):

    config_path = os.path.dirname(os.path.realpath(config_filename))

    with open(config_filename) as f:
        text = f.read()
        text = strip_comments(text)
        config = json.loads(text)

    config["runtime"] = dict()

    config["runtime"]["config_path"] = config_path
    config["runtime"]["filename_so_abs"] = os.path.normpath(
        os.path.join(config_path, config["filename_so"])
    )

    config["runtime"]["sha"] = git.Repo(
        search_parent_directories=True
    ).head.commit.hexsha

    config["runtime"]["genes_dict"] = create_genes_dict_from_config(config)
    config["runtime"]["constants_dict"] = create_constants_dict_from_config(config)

    m_index_tuples = [
        (exp_cond_name, gene_name)
        for exp_cond_name, gene in config["runtime"]["genes_dict"].items()
        for gene_name in gene
    ]
    m_index = pd.MultiIndex.from_tuples(m_index_tuples)
    m_index.names = ["ec_name", "g_name"]

    config["runtime"]["m_index"] = m_index

    legend = dict()
    legend["states"] = pd.read_csv(
        os.path.join(config_path, config["filename_legend_states"]),
        usecols=["name", "value"],
        index_col="name",
    )[
        "value"
    ]  # Series
    legend["constants"] = pd.read_csv(
        os.path.join(config_path, config["filename_legend_constants"]),
        usecols=["name", "value"],
        index_col="name",
    )[
        "value"
    ]  # Series
    config["runtime"]["legend"] = legend

    for exp_cond_name, exp_cond in config["experimental_conditions"].items():

        if exp_cond_name == "common":
            continue

        filename_phenotype = os.path.normpath(
            os.path.join(config_path, exp_cond["filename_phenotype"])
        )
        exp_cond["phenotype"] = pd.read_csv(filename_phenotype)
        exp_cond["filename_phenotype"] = filename_phenotype

        filename_state = os.path.normpath(
            os.path.join(config_path, exp_cond["filename_state"])
        )
        exp_cond["initial_state"] = pd.read_csv(filename_state, index_col=0).iloc[:, -1]
        exp_cond["filename_state"] = filename_state

        column_stim_protocol = config.get("column_stim_protocol", None)
        if column_stim_protocol is not None:
            filename_stim_protocol = os.path.normpath(
                os.path.join(config_path, exp_cond["filename_stim_protocol"])
            )
            exp_cond["stim_protocol"] = pd.read_csv(
                filename_stim_protocol
            )
            exp_cond["filename_stim_protocol"] = filename_stim_protocol
        else:
            exp_cond["stim_protocol"] = None

    states_initial = pd.DataFrame(
        data={
            exp_cond_name: exp_cond["initial_state"].copy()
            for exp_cond_name, exp_cond in config["experimental_conditions"].items()
            if exp_cond_name != "common"
        }
    )

    config["runtime"]["states_initial"] = states_initial

    bounds, gammas, mask_multipliers = generate_bounds_gammas_mask_multipliers(
        config["runtime"]["genes_dict"]
    )
    config["runtime"]["bounds"] = bounds
    config["runtime"]["gammas"] = gammas
    config["runtime"]["mask_multipliers"] = mask_multipliers

    config["runtime"]["kw_ga"] = dict(
        crossover_rate=config.get("crossover_rate", 1.0),
        mutation_rate=config.get("mutation_rate", 1.0),
        selection_force=config.get("selection_force", 2),
        gamma=config.get("gamma", 1.0),
    )

    seed = config.get("seed", None)
    if seed is None:
        sq = np.random.SeedSequence()
        seed = sq.entropy
    config["runtime"]["seed"] = seed

    return config


def update_output_dict(config):

    folder = os.path.normpath(
        os.path.join(
            config["runtime"]["config_path"],
            config.get("output_folder_name", "./results"),
        )
    )
    time_suffix = datetime.now().strftime("%y%m%d_%H%M%S")
    folder = os.path.join(folder, time_suffix)
    config["runtime"]["time_suffix"] = time_suffix

    config["runtime"]["output"] = dict(
        folder=folder,
        folder_dump=os.path.join(folder, "dump"),
        folder_best=os.path.join(folder, "best"),
        folder_phenotype=os.path.join(folder, "phenotype"),
        folder_state=os.path.join(folder, "state"),
    )


def backup_config(config):
    filename = os.path.join(
        config["runtime"]["output"]["folder"], "config_backup.pickle"
    )
    with open(filename, "wb") as f:
        pickle.dump(config, f)
    config["runtime"]["output"]["config_backup"] = filename


def dump_dict(dct, folder):

    if not os.path.isdir(folder):
        os.mkdir(folder)

    for key, value in dct.items():

        filename = os.path.join(folder, key)
        if not os.path.isfile(filename):
            with open(filename, "wb") as _:
                pass

        with open(filename, "ba+") as f:
            np.asarray(value).tofile(f)


def dump_epoch(recvbuf_dict, config):
    dump_dict(recvbuf_dict, config["runtime"]["output"]["folder_dump"])


def save_sol_best(sol_best, config, save_phenotype=0):

    #  save_phenotype : 0, 1, 2

    output_dict = config["runtime"]["output"]

    genes = pd.Series(sol_best.x, index=config["runtime"]["m_index"])
    filename = os.path.join(output_dict["folder"], "sol_best.csv")
    genes.to_csv(filename)

    filename = os.path.join(output_dict["folder"], "state_best.csv")
    sol_best["state"].to_csv(filename)

    d = dict(
        genes=sol_best.x,
        state=sol_best["state"].values,
        loss=sol_best.y,
        status=sol_best.status,
    )

    folder_best = output_dict["folder_best"]
    dump_dict(d, folder_best)

    exp_cond_names = [
        exp_cond_name
        for exp_cond_name in config["experimental_conditions"]
        if exp_cond_name != "common"
    ]

    state_best_0 = {
        exp_cond_name: sol_best["phenotype"][exp_cond_name].iloc[0]
        for exp_cond_name in exp_cond_names
    }
    state_best_0 = pd.concat(state_best_0, axis=1)
    filename = os.path.join(output_dict["folder"], "state_best_0.csv")
    state_best_0.to_csv(filename)

    if save_phenotype >= 1:
        folder_phenotype = config["runtime"]["output"]["folder_phenotype"]
        if not os.path.isdir(folder_phenotype):
            os.mkdir(folder_phenotype)

        for exp_cond_name in exp_cond_names:

            df = sol_best["phenotype"][exp_cond_name]

            # Rewrite last epoch
            filename = os.path.join(folder_phenotype, f"phenotype_{exp_cond_name}.csv")
            df.to_csv(filename, index=False)

            # Append last epoch to previous
            if save_phenotype >= 2:
                filename = os.path.join(folder_phenotype, f"phenotype_{exp_cond_name}")
                if not os.path.isfile(filename):
                    with open(filename, "wb") as f:
                        pass

                with open(filename, "ba+") as f:
                    df.values.astype(np.float32).tofile(f)


def collect_results(case, dirname_results, dump_keys=None, voigt=False):

    if voigt:
        group, cell, suffix = case.split("/")
        group, cell = int(group[-1]), int(cell[-1])
    else:
        group, cell, suffix = None, None, case

    config_path = os.path.join(dirname_results, case)
    with open(os.path.join(config_path, "config_backup.pickle"), "rb") as f:
        config = pickle.load(f)

    m_index = config["runtime"]["m_index"]

    # n_genes = len(m_index)
    # n_organisms = config['runtime']['n_organisms']

    if dump_keys is None:
        dump_keys = ["dump", "best"]

    dump = {}
    for folder in dump_keys:
        dump[folder] = {}
        for key in "genes", "state", "status", "loss":
            filename = os.path.join(config_path, folder, key)
            if os.path.isfile(filename):
                dump[folder][key] = np.fromfile(filename)

    filename = os.path.join(config_path, "sol_best.csv")
    if os.path.isfile(filename):
        sol_best = pd.read_csv(filename, index_col=[0, 1]).iloc[:, -1]
    else:
        sol_best = None

    phenotype_best = {}

    for exp_cond_name in config["experimental_conditions"]:

        if exp_cond_name == "common":
            continue

        filename = os.path.join(
            config_path, "phenotype", f"phenotype_{exp_cond_name}.csv"
        )
        if os.path.isfile(filename):
            try:
                phenotype_best[exp_cond_name] = pd.read_csv(filename)
            except pd.errors.EmptyDataError as e:
                print(f"{filename} is empty")
                continue

    filename = os.path.join(config_path, "state_best.csv")
    state_best = pd.read_csv(filename, index_col=0)

    results = dict(
        trio=(group, cell, suffix),
        config=config,
        dump=dump,
        sol_best=sol_best,
        phenotype_best=phenotype_best,
        state_best=state_best,
    )

    return results


def strip_comments(code, comment_char='#'):
    lines = []
    for i, line in enumerate(code.split('\n')):
        items = line.split(comment_char)
        if len(items) >= 2:
            line = items[0]
        lines.append(line)
    code = '\n'.join(lines)
    return code
