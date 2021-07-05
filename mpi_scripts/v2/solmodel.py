import numpy as np
import pandas as pd

from pypoptim.algorythm import Solution

from gene_utils import update_S_C_from_genes, \
                       update_genes_from_state

from model_utils import calculate_n_samples_per_stim


def calculate_loss(pred, config):
    loss = 0

    columns_control = ['V']
    columns_model = ['V']

    for exp_cond_name, exp_cond in config['experimental_conditions'].items():

        if exp_cond_name == 'common':
            continue

        n_samples_per_stim = calculate_n_samples_per_stim(exp_cond_name, config)

        phenotype_control = exp_cond['phenotype'][columns_control][-n_samples_per_stim - 1:]
        phenotype_model = pred['phenotype'][exp_cond_name][columns_model][-n_samples_per_stim - 1:]

        phenotype_model = phenotype_model[:len(phenotype_control)]

        loss += float(np.sqrt(np.mean((phenotype_control.values - phenotype_model.values) ** 2)))

    return loss


class SolModel(Solution):

    def __init__(self, x, **kwargs_data):
        super().__init__(x, **kwargs_data)
        for attr in 'model', 'config':
            if not hasattr(self, attr):
                raise AttributeError(attr)

        self._status = None
        self.__status_valid = 2

    @property
    def status(self):
        return self._status

    def update(self):

        self['phenotype'] = {}

        legend = self.config['runtime']['legend']

        genes = pd.Series(self.x, index=self.config['runtime']['m_index'])

        for exp_cond_name in self.config['experimental_conditions']:

            if exp_cond_name == 'common':
                continue

            C = legend['constants'].copy()
            S = self['state'][exp_cond_name].copy()

            update_S_C_from_genes(S, C, genes, exp_cond_name, self.config)

            stim_protocol = self.config['experimental_conditions'][exp_cond_name]['stim_protocol']
            pred = self.model.run(S, C, stim_protocol=stim_protocol,
                                  **self.config)

            self._status = self.model.status
            if self._status != self.__status_valid:
                self._x = genes.values
                self._y = np.nan
                return

            update_genes_from_state(genes=genes, state=self['state'],
                                    config=self.config, exp_cond_name=exp_cond_name)

            self['phenotype'][exp_cond_name] = pred.copy()
            self['state'][exp_cond_name] = self['phenotype'][exp_cond_name].iloc[-1]

        self._x = genes.values
        self._y = calculate_loss(self, self.config)

    def is_valid(self):
        if not self.is_updated():
            return False
        else:
            flag_valid = self._status == self.__status_valid and np.isfinite(self._y)
            if 'phenotype' not in self:  # solution was gathered via MPI
                return flag_valid
            else:
                return flag_valid and all(not np.any(np.isnan(p)) for p in self['phenotype'].values())
