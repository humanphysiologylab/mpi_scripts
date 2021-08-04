import numpy as np


def calculate_n_samples_per_stim(exp_cond_name, config):
    stim_period_legend_name = config["stim_period_legend_name"]
    stim_period = config["experimental_conditions"][exp_cond_name]["params"][
        stim_period_legend_name
    ]
    t_sampling = config["t_sampling"]
    n_samples_per_stim = int(np.round(stim_period / t_sampling))
    return n_samples_per_stim
