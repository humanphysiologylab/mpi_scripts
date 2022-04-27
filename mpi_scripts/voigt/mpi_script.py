#!/usr/bin/env python

import os
import logging
import argparse

from mpi4py import MPI
from tqdm.auto import tqdm

import numpy as np

from cardiac_model import CardiacModel
from solmodel import SolModel
from pypoptim.algorythm.ga import GA

from pypoptim.helpers import argmin, is_values_inside_bounds

from io_utils import (
    prepare_config,
    update_output_dict,
    backup_config,
    dump_epoch,
    save_sol_best,
)
from mpi_utils import allocate_recvbuf, allgather, population_from_recvbuf


def mpi_script(config_filename):
    logger = logging.getLogger(__name__)

    comm = MPI.COMM_WORLD
    comm_rank = comm.Get_rank()
    comm_size = comm.Get_size()

    config = None
    if comm_rank == 0:
        config = prepare_config(config_filename)
        config["runtime"]["comm_size"] = comm_size

        for key, value in config["runtime"]["VCS"].items():
            print(f"# VCS: {key} {value['branch']} {value['commit']}")

        print(f"# size: {comm_size}")
        print(f"# seed: {config['runtime']['seed']}")

        if config["n_organisms"] % comm_size != 0:
            config["runtime"]["n_organisms"] = int(
                np.ceil(config["n_organisms"] / comm_size) * comm_size
            )
            print(
                f'# n_organisms: {config["n_organisms"]} to {config["runtime"]["n_organisms"]}',
                flush=True,
            )
        else:
            config["runtime"]["n_organisms"] = config["n_organisms"]

        update_output_dict(config)
        os.makedirs(config["runtime"]["output"]["folder"])
        print(f"# folder: {config['runtime']['output']['folder']}", flush=True)

    config = comm.bcast(config, root=0)

    recvbuf_dict = allocate_recvbuf(config, comm)

    model = CardiacModel(config["runtime"]["filename_so_abs"])
    SolModel.model = model
    SolModel.config = config

    rng = np.random.Generator(np.random.PCG64(config["runtime"]["seed"] + comm_rank))
    ga_optim = GA(
        SolModel,
        bounds=config["runtime"]["bounds"],
        gammas=config["runtime"]["gammas"],
        gamma_default=config["runtime"]["kw_ga"]["gamma"],
        mutation_rate=config["runtime"]["kw_ga"]["mutation_rate"],
        crossover_rate=config["runtime"]["kw_ga"]["crossover_rate"],
        selection_force=config["runtime"]["kw_ga"]["selection_force"],
        mask_log10_scale=config["runtime"]["mask_multipliers"],
        keys_data_transmit=["state"],
        rng=rng,
    )

    initial_population_filename = config.get("initial_population_filename", None)
    if initial_population_filename is not None:
        raise NotImplementedError

    if comm_rank == 0:
        backup_config(config)

    batch = ga_optim.generate_population(config["runtime"]["n_orgsnisms_per_process"])
    for sol in batch:
        sol["state"] = config["runtime"]["states_initial"].copy()

    if comm_rank == 0:
        pbar = tqdm(total=config["n_generations"], ascii=True)

    for epoch in range(config["n_generations"]):

        if comm_rank == 0:
            pbar.set_postfix_str("CALC")
        for i, sol in enumerate(batch):
            sol.update()
            if not (sol.is_valid() and ga_optim.is_solution_inside_bounds(sol)):
                sol._y = np.inf
                del sol.data["phenotype"]

        if comm_rank == 0:
            pbar.set_postfix_str("GATHER")

        allgather(batch, recvbuf_dict, comm)
        population = population_from_recvbuf(recvbuf_dict, SolModel, config)

        n_orgsnisms_per_process = config["runtime"]["n_orgsnisms_per_process"]
        shift = comm_rank * n_orgsnisms_per_process
        assert all(
            sol_b.is_all_equal(sol_p) for sol_b, sol_p in zip(batch, population[shift:])
        )

        if comm_rank == 0:
            pbar.set_postfix_str("SAVE")

        index_best = argmin(population)
        assert population[index_best] is min(population)
        comm_rank_best = index_best // config["runtime"]["n_orgsnisms_per_process"]
        index_best_batch = index_best % config["runtime"]["n_orgsnisms_per_process"]

        if comm_rank == comm_rank_best:
            sol_best = batch[index_best_batch]

            assert sol_best is min(batch)
            assert sol_best.is_all_equal(min(population))

            msg = f"{comm_rank} has best solution:\n{sol_best}"
            logger.debug(msg)
            save_sol_best(sol_best, config)

            assert sol_best.is_updated()
            assert sol_best.is_valid()
            assert is_values_inside_bounds(sol_best.x, config["runtime"]["bounds"])
            assert ga_optim.is_solution_inside_bounds(sol_best)

        if comm_rank == (comm_rank_best + 1) % comm_size:
            dump_epoch(recvbuf_dict, config)

        if comm_rank == 0:
            pbar.set_postfix_str("GENE")

        population = ga_optim.filter_population(population)
        # n_invalids = config['runtime']['n_organisms'] - len(population)
        # percentage_invalid = n_invalids / config['runtime']['n_organisms'] * 100
        population.sort()

        if len(population) <= 3:
            if comm_rank == 0:
                msg = f"# Not enough organisms for genetic operations left: {len(population)}"
                raise RuntimeError(msg)

        elites_all = population[: config["n_elites"]]
        # len may be less than config['n_elites'] due to invalids
        
        elites_batch = elites_all[comm_rank::comm_size]  # elites_batch may be empty
        n_elites = len(elites_batch)
        n_mutants = config["runtime"]["n_orgsnisms_per_process"] - n_elites

        mutants_batch = ga_optim.get_mutants(population, n_mutants)
        batch = elites_batch + mutants_batch

        assert len(batch) == config["runtime"]["n_orgsnisms_per_process"]

        if comm_rank == 0:
            pbar.update(1)
            pbar.refresh()

    if comm_rank == 0:
        pbar.set_postfix_str("DONE")
        pbar.refresh()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="configuration file")
    parser.add_argument("-v", "--verbosity", type=str, help="logging level")

    args = parser.parse_args()

    config_filename = args.config
    logging_level = args.verbosity
    level = dict(INFO=logging.INFO, DEBUG=logging.DEBUG).get(
        logging_level, logging.WARNING
    )
    logging.basicConfig(level=level)
    logging.getLogger("numba").setLevel(
        logging.CRITICAL
    )  # https://stackoverflow.com/a/63471108/13213091

    mpi_script(config_filename=config_filename)
