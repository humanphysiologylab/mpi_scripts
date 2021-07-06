import numpy as np

from pypoptim.helpers import calculate_mean_abs_noise, autoscaling
from pypoptim.losses import calculate_RMSE_balanced, RMSE, calculate_RMSE_weightened
from model_utils import calculate_n_samples_per_stim


def calculate_composite_RMSE_V_CaT(x, y):
    # x -- model, y -- experiment
    assert (len(x) == len(y))

    v_model, cat_model = x.T
    v_exp, cat_exp = y.T

    cat_model = (cat_model - cat_model.min(axis=0)) / cat_model.ptp(axis=0)  # to balance V and CaT

    rmse_v = calculate_RMSE_balanced(v_model, v_exp)    # v_exp --> [0, 1]
    cat_exp_scaled, coeffs = autoscaling(signal_to_scale=cat_exp,
                                         signal_reference=cat_model)
    rmse_cat = RMSE(cat_exp_scaled, cat_model)
    rmse_total = rmse_v + rmse_cat

    return rmse_total


def _calculate_RMSE_weightened(phenotype_model, phenotype_control):

    phenotype_model = phenotype_model.to_numpy()
    phenotype_control = phenotype_control.to_numpy()

    phenotype_model = (phenotype_model - phenotype_control.min(axis=0)) / phenotype_control.ptp(axis=0)
    phenotype_control = (phenotype_control - phenotype_control.min(axis=0)) / phenotype_control.ptp(axis=0)

    weights = 1 / calculate_mean_abs_noise(phenotype_control)
    weights /= sum(weights)
    return calculate_RMSE_weightened(phenotype_control, phenotype_model, weights)


def composite_RMSE_V_CaT_noisy(phenotype_model, phenotype_control):
    phenotype_model = phenotype_model.to_numpy()
    phenotype_control = phenotype_control.to_numpy()

    phenotype_model = (phenotype_model - phenotype_control.min(axis=0)) / phenotype_control.ptp(axis=0)
    phenotype_control = (phenotype_control - phenotype_control.min(axis=0)) / phenotype_control.ptp(axis=0)

    weights = 1 / calculate_mean_abs_noise(phenotype_control)
    weights /= sum(weights)

    rmse_v = RMSE(phenotype_control[:, 0], phenotype_model[:, 0])

    ca_exp_scaled, rmse_ca, coeffs = autoscaling(signal_to_scale=phenotype_control[:, 1],
                                                 signal_reference=phenotype_model[:, 1])

    return rmse_v * weights[0] + rmse_ca * weights[1]


def align_rise(signal_model, signal_control):

    shift_model = np.argwhere(signal_model > np.mean(signal_model))[0, 0]
    shift_control = np.argwhere(signal_control > np.mean(signal_control))[0, 0]

    shift = shift_model - shift_control
    signal_control = np.roll(signal_control, shift)

    return signal_control, shift


def trim_signals(*signals):
    min_len = min(map(len, signals))
    return [signal[:min_len] for signal in signals]


def calculate_loss_cat(cat_control, cat_model):
    cat_control, _ = align_rise(signal_model=cat_model.to_numpy(),
                                signal_control=cat_control.to_numpy())
    cat_control, cat_model = trim_signals(cat_control, cat_model)
    # cat_control - uM, cat_model - mM
    loss = RMSE(cat_control, cat_model * 1000)
    return loss


def calculate_loss_ap(ap_control, ap_model):
    ap_control, ap_model = trim_signals(ap_control, ap_model)
    loss = RMSE(ap_control, ap_model)
    return loss


def calculate_loss(sol, config):

    loss = 0

    for exp_cond_name, exp_cond in config['experimental_conditions'].items():

        if exp_cond_name == 'common':
            continue

        n_samples_per_stim = calculate_n_samples_per_stim(exp_cond_name, config)

        loss_cat = 0
        loss_ap = 0
        for name in exp_cond['type']:
            df = exp_cond[name]
            if name == 'CaT':
                cat_control = df['CaT']
                cat_model = sol['phenotype'][exp_cond_name]['Ca_i'][-n_samples_per_stim:]
                loss_cat = calculate_loss_cat(cat_control, cat_model)
            if name == 'AP':
                ap_control = df['V']
                ap_model = sol['phenotype'][exp_cond_name]['V'][-n_samples_per_stim:]
                loss_ap = calculate_loss_ap(ap_control, ap_model)

        loss += loss_cat * 100 + loss_ap

    return loss


