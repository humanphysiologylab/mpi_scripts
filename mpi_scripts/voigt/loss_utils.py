import numpy as np

from pypoptim.helpers import calculate_mean_abs_noise, calculate_autoscaling
from pypoptim.losses import calculate_RMSE_balanced, RMSE, calculate_RMSE_weightened
from model_utils import calculate_n_samples_per_stim

import logging

logger = logging.getLogger(__name__)


def calculate_composite_RMSE_V_CaT(x, y):
    # x -- model, y -- experiment
    assert len(x) == len(y)

    v_model, cat_model = x.T
    v_exp, cat_exp = y.T

    cat_model = (cat_model - cat_model.min(axis=0)) / cat_model.ptp(
        axis=0
    )  # to balance V and CaT

    rmse_v = calculate_RMSE_balanced(v_model, v_exp)  # v_exp --> [0, 1]
    cat_exp_scaled, coeffs = calculate_autoscaling(
        signal_to_scale=cat_exp, signal_reference=cat_model
    )
    rmse_cat = RMSE(cat_exp_scaled, cat_model)
    rmse_total = rmse_v + rmse_cat

    return rmse_total


def _calculate_RMSE_weightened(phenotype_model, phenotype_control):

    phenotype_model = phenotype_model.to_numpy()
    phenotype_control = phenotype_control.to_numpy()

    phenotype_model = (
        phenotype_model - phenotype_control.min(axis=0)
    ) / phenotype_control.ptp(axis=0)
    phenotype_control = (
        phenotype_control - phenotype_control.min(axis=0)
    ) / phenotype_control.ptp(axis=0)

    weights = 1 / calculate_mean_abs_noise(phenotype_control)
    weights /= sum(weights)
    return calculate_RMSE_weightened(phenotype_control, phenotype_model, weights)


def composite_RMSE_V_CaT_noisy(phenotype_model, phenotype_control):
    phenotype_model = phenotype_model.to_numpy()
    phenotype_control = phenotype_control.to_numpy()

    phenotype_model = (
        phenotype_model - phenotype_control.min(axis=0)
    ) / phenotype_control.ptp(axis=0)
    phenotype_control = (
        phenotype_control - phenotype_control.min(axis=0)
    ) / phenotype_control.ptp(axis=0)

    weights = 1 / calculate_mean_abs_noise(phenotype_control)
    weights /= sum(weights)

    rmse_v = RMSE(phenotype_control[:, 0], phenotype_model[:, 0])

    ca_exp_scaled, rmse_ca, coeffs = calculate_autoscaling(
        signal_to_scale=phenotype_control[:, 1], signal_reference=phenotype_model[:, 1]
    )

    return rmse_v * weights[0] + rmse_ca * weights[1]


def align_depolarization(phenotype_model, phenotype_control):

    column_v_model = "V" if "V" in phenotype_model else "v"
    column_v_control = "V" if "V" in phenotype_control else "v"

    v_model = phenotype_model[column_v_model].to_numpy()
    v_control = phenotype_control[column_v_control].to_numpy()

    v_level = np.mean(v_model)  # 0 # mV
    shift_model = np.where(v_model > v_level)[0][0]
    shift_control = np.where(v_control > v_level)[0]

    # if not len(shift_control):
    #     loss = np.inf
    #     organism['fitness'] = -loss
    #     return

    shift_control = shift_control[0]

    shift = shift_model - shift_control

    phenotype_control = np.roll(phenotype_control, shift, axis=1)

    return phenotype_control


def calculate_V_CaT_shared(sol, config):

    columns_control = config.get("columns_control", ["V"])
    columns_model = config.get("columns_model", ["V"])

    phenotype_model_list = []
    phenotype_control_list = []

    for exp_cond_name, exp_cond in config["experimental_conditions"].items():

        if exp_cond_name == "common":
            continue

        n_samples_per_stim = calculate_n_samples_per_stim(exp_cond_name, config)

        phenotype_control = exp_cond["phenotype"][columns_control].copy()[
            -n_samples_per_stim - 1 :
        ]
        phenotype_model = sol["phenotype"][exp_cond_name][columns_model].copy()[
            -n_samples_per_stim - 1 :
        ]

        phenotype_model = phenotype_model[: len(phenotype_control)]

        phenotype_model_list.append(phenotype_model.values)
        phenotype_control_list.append(phenotype_control.values)

    cat_model_concat = np.concatenate([x[:, 1] for x in phenotype_model_list])
    cat_control_concat = np.concatenate([x[:, 1] for x in phenotype_control_list])

    cat_control_concat_scaled, (alpha, beta) = calculate_autoscaling(
        signal_to_scale=cat_control_concat, signal_reference=cat_model_concat
    )

    if alpha <= 0:
        return np.inf

    cumlen = 0
    loss = 0

    for i, x in enumerate(phenotype_control_list):
        phenotype_control_list[i][:, 1] = cat_control_concat_scaled[
            cumlen : cumlen + len(x)
        ]
        cumlen += len(x)

        phenotype_control = phenotype_control_list[i]
        phenotype_model = phenotype_model_list[i]

        phenotype_control = (
            phenotype_control - phenotype_model.min(axis=0)
        ) / phenotype_model.ptp(axis=0)
        phenotype_model = (
            phenotype_model - phenotype_model.min(axis=0)
        ) / phenotype_model.ptp(axis=0)

        weights = 1 / calculate_mean_abs_noise(phenotype_control)
        weights /= sum(weights)

        # rmse = calculate_RMSE_weightened(phenotype_control, phenotype_model, weights)
        error = (phenotype_control - phenotype_model) ** 2
        error[:10] *= 10
        error = np.sqrt(np.mean(error, axis=0))
        error *= weights
        rmse = np.sum(error)

        loss += rmse

    return loss


def calculate_loss(sol, config):

    columns_control = config.get("columns_control", ["V"])
    columns_model = config.get("columns_model", ["V"])

    loss = 0
    loss_L2 = 0

    if config.get("ridge_regression", False):
        scalers = sol.x[config["runtime"]["mask_multipliers"]]
        loss_L2 = np.sum((np.log10(scalers)) ** 2) / 10

    if config["loss"] == "V_CaT_shared":
        loss = calculate_V_CaT_shared(sol, config)

    else:

        for exp_cond_name, exp_cond in config["experimental_conditions"].items():

            if exp_cond_name == "common":
                continue

            n_samples_per_stim = calculate_n_samples_per_stim(exp_cond_name, config)

            phenotype_control = exp_cond["phenotype"][columns_control][
                -n_samples_per_stim - 1 :
            ]
            phenotype_model = sol["phenotype"][exp_cond_name][columns_model][
                -n_samples_per_stim - 1 :
            ]

            phenotype_model = phenotype_model[: len(phenotype_control)]

            if config.get("align_depolarization", False):
                phenotype_control = align_depolarization(
                    phenotype_model, phenotype_control
                )

            if config["loss"] == "RMSE":
                loss += RMSE(phenotype_control, phenotype_model)

            elif config["loss"] == "RMSE_balanced":
                loss += calculate_RMSE_balanced(phenotype_control, phenotype_model)

            elif config["loss"] == "RMSE_weightened":
                loss += _calculate_RMSE_weightened(phenotype_control, phenotype_model)

            elif config["loss"] == "composite_RMSE_V_CaT":
                loss += calculate_composite_RMSE_V_CaT(
                    phenotype_control.to_numpy(), phenotype_model.to_numpy()
                )

            elif config["loss"] == "composite_RMSE_V_CaT_noisy":
                loss += composite_RMSE_V_CaT_noisy(phenotype_model, phenotype_control)

            else:
                raise ValueError(f'Unknown loss {config["loss"]}')

    logger.debug(f"loss = {loss}; loss_L2 = {loss_L2}")

    return loss
