import sys
import os

from mpi4py import MPI
from tqdm.auto import tqdm

import numpy as np

from pypoptim.model import CardiacModel
from solmodel import SolModel
from pypoptim.algorythm.ga import GA

from pypoptim.helpers import argmax
from pypoptim import Timer

from io_utils import prepare_config, update_output_dict, backup_config, dump_epoch, save_sol_best
from mpi_utils import allocate_recvbuf, allgather, population_from_recvbuf


comm = MPI.COMM_WORLD
comm_rank = comm.Get_rank()
comm_size = comm.Get_size()

if len(sys.argv) < 2:
    if comm_rank == 0:
        print(f"Usage: mpiexec -n 2 python {sys.argv[0]} config.json")
    exit()

config = None
if comm_rank == 0:
    config_filename = sys.argv[1]
    config = prepare_config(config_filename)
    config['runtime']['comm_size'] = comm_size

    print(f"# commit: {config['runtime']['sha']}")
    print(f'# size: {comm_size}')

    if config['n_organisms'] % comm_size != 0:
        config['runtime']['n_organisms'] = int(np.ceil(config['n_organisms'] / comm_size) * comm_size)
        print(f'# n_organisms: {config["n_organisms"]} to {config["runtime"]["n_organisms"]}',
              flush=True)
    else:
        config['runtime']['n_organisms'] = config['n_organisms']

    update_output_dict(config)
    os.makedirs(config['runtime']['output']['folder'])
    print(f"# folder: {config['runtime']['output']['folder']}", flush=True)

config = comm.bcast(config, root=0)

recvbuf_dict = allocate_recvbuf(config, comm)

model = CardiacModel(config['runtime']['filename_so_abs'])
SolModel.model = model
SolModel.config = config

ga_optim = GA(SolModel,
              config['runtime']['bounds'],
              config['runtime']['gammas'],
              config['runtime']['mask_multipliers'],
              keys_data_transmit=['state'])

initial_population_filename = config.get('initial_population_filename', None)
if initial_population_filename is not None:
    raise NotImplementedError

if comm_rank == 0:
    backup_config(config)

batch = ga_optim.generate_population(config['runtime']['n_orgsnisms_per_process'])
for sol in batch:
    sol['state'] = config['runtime']['states_initial']

timer = Timer()

if comm_rank == 0:
    pbar = tqdm(total=config['n_generations'])

for epoch in range(config['n_generations']):

    timer.start('calc')
    for i, sol in enumerate(batch):
        sol.update()
        if not sol.is_valid():
            sol._y = np.inf
    timer.end('calc')

    timer.start('gather')
    allgather(batch, recvbuf_dict, comm)
    population = population_from_recvbuf(recvbuf_dict, SolModel, config)
    timer.end('gather')

    timer.start('save')
    index_best = argmax(population)
    comm_rank_best = index_best // config['runtime']['n_orgsnisms_per_process']
    index_best_batch = index_best % config['runtime']['n_orgsnisms_per_process']

    if comm_rank == comm_rank_best:
        sol_best = batch[index_best_batch]
        save_sol_best(sol_best, config)

    if comm_rank == (comm_rank_best + 1) % comm_size:
        dump_epoch(recvbuf_dict, config)
    timer.end('save')

    timer.start('gene')

    population = ga_optim.filter_population(population)
    # n_invalids = config['runtime']['n_organisms'] - len(population)
    # percentage_invalid = n_invalids / config['runtime']['n_organisms'] * 100
    population.sort()

    if len(population) <= 3:
        if comm_rank == 0:
            print(f"# Not enough organisms for genetic operations left: {len(population)}")
            print("# exit")
        exit()

    elites_all = population[:config['n_elites']]  # len may be less than config['n_elites'] due to invalids
    elites_batch = elites_all[comm_rank::comm_size]  # elites_batch may be empty
    n_elites = len(elites_batch)
    n_mutants = config['runtime']['n_orgsnisms_per_process'] - n_elites

    mutants_batch = ga_optim.get_mutants(population, n_mutants)
    batch = elites_batch + mutants_batch

    assert (len(batch) == config['runtime']['n_orgsnisms_per_process'])

    timer.end('gene')

    if comm_rank == 0:
        with open(os.path.join(config['runtime']['output']['folder'], 'runtime.log'), 'w') as f:
            print(timer.report(), file=f)
        pbar.update(1)

    timer.clear()

