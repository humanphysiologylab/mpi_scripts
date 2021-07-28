import numpy as np
import pandas as pd

from pypoptim.algorythm import Solution

from gene_utils import update_S_C_from_genes, \
                       update_genes_from_state

from loss_utils import calculate_loss

import logging
logger = logging.getLogger(__name__)

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

            self['phenotype'][exp_cond_name] = pred.copy()
            self['state'][exp_cond_name] = self['phenotype'][exp_cond_name].iloc[-1]

            update_genes_from_state(genes=genes,
                                    state=self['state'],
                                    config=self.config,
                                    exp_cond_name=exp_cond_name)

        self._x = genes.values
        self._y = calculate_loss(self, self.config)

    def is_all_equal(self, other, keys_check=None):
        if not np.allclose(self.x, other.x):
            x = np.vstack([self.x, other.x, self.x - other.x]).T
            logger.debug(f"`x`s differs: {x}")
            return False

        if self.y != other.y:
            logger.debug(f"`y`s differs: {self.y} {other.y}")
            return False

        if keys_check is None:
            keys_check = ['state']
        for key in keys_check:
            if key in self:
                if key in other:
                    if not np.allclose(self[key], other[key]):
                        logger.debug(f"`{key}`s differs")
                        return False
                else:
                    logger.debug(f"`{key}` is not in `other`")
                    return False
            else:
                if key in other:
                    logger.debug(f"`{key}` is not in `self`")
                    return False

        return True

    def is_valid(self):
        if not self.is_updated():
            return False
        else:
            flag_valid = self._status == self.__status_valid and np.isfinite(self._y)
            if 'phenotype' not in self:  # solution was gathered via MPI
                return flag_valid
            else:
                return flag_valid and all(not np.any(np.isnan(p)) for p in self['phenotype'].values())
