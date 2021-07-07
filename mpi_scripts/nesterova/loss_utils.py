import numpy as np

from pypoptim.losses import RMSE
from model_utils import calculate_n_samples_per_stim

import logging
logger = logging.getLogger(__name__)

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


def calculate_DCa(cat_signal):
    return np.min(cat_signal)


def calculate_CTA(cat_signal):
    return np.ptp(cat_signal)


def calcualate_loss_default(sol, exp_cond_name, config):

    n_samples_per_stim = calculate_n_samples_per_stim(exp_cond_name, config)
    exp_cond = config['experimental_conditions'][exp_cond_name]

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

    loss = loss_cat * 100 + loss_ap
    return loss


def calculate_loss_AP_CaT_restcurves(sol, exp_cond_name, config):

    n_samples_per_stim = calculate_n_samples_per_stim(exp_cond_name, config)
    exp_cond = config['experimental_conditions'][exp_cond_name]

    loss_cat = 0
    loss_ap = 0
    for name in exp_cond['type']:
        if name == 'CaT':
            cat_params = exp_cond['CaT_params']
            cat_model = sol['phenotype'][exp_cond_name]['Ca_i'][-n_samples_per_stim:]
            DCa_control = cat_params['DCa']['mean']
            CTA_control = cat_params['CTA']['mean']
            DCa_model = calculate_DCa(cat_model)
            CTA_model = calculate_DCa(cat_model)
            loss_cat = 10 * np.abs(DCa_control - DCa_model) \
                       + np.abs(CTA_control - CTA_model)
            # explanation: DCa ~ 10 times smaller than CTA
        if name == 'AP':
            ap_control = exp_cond[name]['V']
            ap_model = sol['phenotype'][exp_cond_name]['V'][-n_samples_per_stim:]
            loss_ap = calculate_loss_ap(ap_control, ap_model)

    loss = loss_cat * 10 + loss_ap  # 10 is a magic weight
    logger.info(f'loss_cat = {loss_cat * 10}; loss_ap = {loss_ap}')
    return loss


def calculate_loss(sol, config):

    loss = 0

    for exp_cond_name, exp_cond in config['experimental_conditions'].items():

        if exp_cond_name == 'common':
            continue

        if config['loss'] == 'AP_CaT_restcurves':
            loss += calculate_loss_AP_CaT_restcurves(sol, exp_cond_name, config)

        else:  # default
            loss += calcualate_loss_default(sol, exp_cond_name, config)

    return loss
