import numpy as np
import pandas as pd


def allocate_buffers(config, comm):

    comm_size = comm.Get_size()

    n_orgsnisms_per_process = config['runtime']['n_organisms'] // comm_size
    config['runtime']['n_orgsnisms_per_process'] = n_orgsnisms_per_process

    genes_size = sum(map(len, config['runtime']['genes_dict'].values()))
    config['runtime']['genes_size'] = genes_size

    state_size = config['runtime']['states_initial'].size
    state_shape = config['runtime']['states_initial'].shape
    config['runtime']['state_shape'] = state_shape

    recvbuf_genes = np.empty([comm_size, n_orgsnisms_per_process * genes_size])
    recvbuf_state = np.empty([comm_size, n_orgsnisms_per_process * state_size])
    recvbuf_fitness = np.empty([comm_size, n_orgsnisms_per_process * 1])

    return recvbuf_genes, recvbuf_state, recvbuf_fitness


def allgather(batch, config, SolModel,
              recvbuf_genes, recvbuf_state, recvbuf_fitness, comm):

    sendbuf_genes = np.concatenate([sol.x for sol in batch])
    sendbuf_state = np.concatenate([sol['state'].values.flatten() for sol in batch])
    sendbuf_fitness = np.array([sol.y for sol in batch])
    assert(not np.any(np.isnan(sendbuf_fitness)))

    comm.Allgatherv(sendbuf_genes, recvbuf_genes)
    comm.Allgatherv(sendbuf_state, recvbuf_state)
    comm.Allgatherv(sendbuf_fitness, recvbuf_fitness)

    recvbuf_genes = recvbuf_genes.reshape((config['runtime']['n_organisms'], config['runtime']['genes_size']))
    recvbuf_state = recvbuf_state.reshape((config['runtime']['n_organisms'], *config['runtime']['state_shape']))
    recvbuf_fitness = recvbuf_fitness.flatten()

    population = []
    columns_state = [exp_cond_name for exp_cond_name in config['experimental_conditions'] if exp_cond_name != 'common']
    for i in range(config['runtime']['n_organisms']):
        state = pd.DataFrame(recvbuf_state[i], columns_state)
        sol = SolModel(recvbuf_genes[i], state=state)
        sol._y = recvbuf_fitness[i]
        population.append(sol)

    return population
