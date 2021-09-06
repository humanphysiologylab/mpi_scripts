def create_genes_dict_from_config(config):
    gamma_default = config.get("gamma", 1)

    genes_dict = {
        ec_name: {
            p_name: dict(
                bounds=p["bounds"],
                gamma=(gamma_default * p["gamma_multiplier"])
                if "gamma_multiplier" in p
                else None,
                is_multiplier=p.get("is_multiplier", False),
            )
            for p_name, p in ec["params"].items()
            if isinstance(p, dict)
        }
        for ec_name, ec in config["experimental_conditions"].items()
    }

    return genes_dict


def create_constants_dict_from_config(config):
    constants_dict = {
        ec_name: {
            p_name: value
            for p_name, value in ec["params"].items()
            if isinstance(value, (int, float))
        }
        for ec_name, ec in config["experimental_conditions"].items()
    }

    return constants_dict


def generate_bounds_gammas_mask_multipliers(genes_dict):
    bounds, gammas, mask_multipliers = [
        [x[item] for y in genes_dict.values() for x in y.values()]
        for item in ["bounds", "gamma", "is_multiplier"]
    ]
    return bounds, gammas, mask_multipliers


def update_S_C_from_genes(S, C, genes, exp_cond_name, config):
    legend = config["runtime"]["legend"]
    genes_dict = config["runtime"]["genes_dict"]
    constants_dict = config["runtime"]["constants_dict"]
    constants_dict_current = {
        **constants_dict["common"],
        **constants_dict[exp_cond_name],
    }

    idx_frame = genes.index.to_frame(index=False)

    for (_, idx), value in zip(idx_frame.iterrows(), genes.values):
        ec_name, g_name = idx["ec_name"], idx["g_name"]

        if ec_name in ["common", exp_cond_name]:

            if g_name in legend["constants"]:
                df = C
            elif g_name in legend["states"]:
                df = S
            else:
                msg = "No such gene in legends: {g_name}".format(g_name=g_name)
                raise KeyError(msg)

            if genes_dict[ec_name][g_name]["is_multiplier"]:
                df[g_name] *= genes[ec_name, g_name]
            else:
                df[g_name] = genes[ec_name, g_name]

    for c_name, c in constants_dict_current.items():
        if c_name in legend["constants"]:
            C[c_name] = c
        elif c_name in legend["states"]:
            S[c_name] = c
        else:
            msg = "No such constant in legends: {c_name}".format(c_name=c_name)
            raise KeyError(msg)


def update_genes_from_state(genes, state, config, exp_cond_name):
    legend = config["runtime"]["legend"]
    genes_dict = config["runtime"]["genes_dict"]

    for i, g_name in enumerate(genes.index.get_level_values("g_name")):

        if g_name in legend["states"].index:
            for ecn in ["common", exp_cond_name]:
                if g_name in genes_dict[ecn]:
                    if genes_dict[ecn][g_name]["is_multiplier"]:
                        genes[ecn, g_name] = (
                            state[exp_cond_name][g_name] / legend["states"][g_name]
                        )
                    else:
                        genes[ecn, g_name] = state[exp_cond_name][g_name]
