from mpi4py import MPI
import itertools
import numpy as np

import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
from src.algorythm.ga import do_step
from src.helpers import rastrigin


def init_population(size, bounds):
    lb, ub = bounds.T
    return [{'genes': np.random.random(len(ub)) * (ub - lb) + lb, 'state': np.array(i)} for i in range(size)]


def calculate_fitness(genes):
    return -rastrigin(genes)


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

n_genes = 30
n_organisms = 2000
n_elites = 100 # max(1, int(n_organisms / 10))
n_epochs = 1000

bounds = np.array([[-5.2, 5.2]] * n_genes)

if rank == 0:
    population = init_population(n_organisms, bounds)
    #population = [population[i::size] for i in range(size)]
    file_dump = open("dump", "wb")
    time_start = time.time()
    print("SIZE =", size, flush=True)
else:
    population = None

times = dict()

for epoch in range(n_epochs):

    times['scat'] = time.time()
    #  population = comm.scatter(population, root=0)
    times['scat'] = time.time() - times['scat']

    times['calc'] = time.time()
    for i in range(len(population)):
        population[i]['fitness'] = calculate_fitness(population[i]['genes'])
    times['calc'] = time.time() - times['calc']

    times['gather'] = time.time()
    #  population = comm.gather(population, root=0)
    times['gather'] = time.time() - times['gather']

    if rank == 0:
        times['master'] = time.time()
        #population = list(itertools.chain(*population))  # flatten

        population.sort(key=lambda x: x['fitness'], reverse=True)
        genes = np.array([x['genes'] for x in population])
        fitness = np.array([x['fitness'] for x in population])

        dump_current = np.hstack([genes, fitness[:, None]])
        dump_current.tofile(file_dump)

        print(np.array2string(dump_current[0, :], precision=3, suppress_small=True,
                              max_line_width=1000)[1:-1], flush=True)

        #  Next generation
        kw_ga = dict(crossover_rate=1.0,
                     mutation_rate=0.1,
                     gamma=1.0)
        population = do_step(population, new_size=len(population),
                             elite_size=n_elites, bounds=bounds,
                             **kw_ga)
        #population = [population[i::size] for i in range(size)]
        times['master'] = time.time() - times['master']

        times['total'] = 0
        times['total'] = sum(times.values())

        continue
        print(f"# EPOCH {epoch}:")
        for t in times:
            print(f"# {t}:\t{times[t]:.6f}\t{100 * times[t] / times['total']:.2f} %")
        print()

if rank == 0:
    time_end = time.time()
    print("TIME =", time_end - time_start, flush=True)
    file_dump.close()
