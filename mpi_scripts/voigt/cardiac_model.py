import os
import ctypes
import pandas as pd
import numpy as np


def wrapped_ndptr(*args, **kwargs):
    #  https://stackoverflow.com/a/37664693/13213091
    base = np.ctypeslib.ndpointer(*args, **kwargs)
    def from_param(cls, obj):
        if obj is None:
            return obj
        return base.from_param(obj)
    return type(base.__name__, (base,), {'from_param': classmethod(from_param)})

DoubleArrayType_1D = wrapped_ndptr(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
DoubleArrayType_2D = wrapped_ndptr(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')


class CardiacModel:

    def __init__(self, filename_so):

        filename_so_abs = os.path.abspath(filename_so)
        ctypes_obj = ctypes.CDLL(filename_so_abs)

        ctypes_obj.run.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # double *S
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),  # double *C
            ctypes.c_int,  # int n_beats
            ctypes.c_double,  # double t_sampling
            ctypes.c_double,  # double tol
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),  # double *output

            DoubleArrayType_2D,  # double *output_A
            DoubleArrayType_1D,  # double *output_t
            DoubleArrayType_1D,  # double *stim_protocol_Ist
            DoubleArrayType_1D,  # double *stim_protocol_t
        ]

        ctypes_obj.run.restype = ctypes.c_int
        self._run = ctypes_obj.run

        self._status = None

    @property
    def status(self):
        return self._status

    def run(self, S, C, **kwargs):

        stim_period = C[kwargs['stim_period_legend_name']]
        t_sampling = kwargs['t_sampling']

        if 'n_beats' in kwargs and 't_run' not in kwargs:
            n_beats = kwargs['n_beats']
        elif 'n_beats' not in kwargs and 't_run' in kwargs:
            t_run = kwargs['t_run']
            n_beats = np.ceil(t_run / stim_period).astype(int)
            if n_beats % 2 == 0:
                n_beats += 1
        else:
            raise ValueError('Invalid config: check n_beats & t_run entries')

        tol = kwargs.get('tol', 1e-6)

        n_samples_per_stim = int(np.round(stim_period / t_sampling))

        output = pd.DataFrame(np.empty((n_samples_per_stim * n_beats + 1, len(S))),
                              columns=S.index)

        stim_protocol = kwargs['stim_protocol']
        if stim_protocol is not None:
            stim_protocol_Ist = stim_protocol[kwargs['column_stim_protocol']].values.copy()
            if kwargs.get('sparsed_protocol', False):
                stim_protocol_t = stim_protocol['t'].values.copy()
            else:
                stim_protocol_t = None
        else:
            stim_protocol_t = None
            stim_protocol_Ist = None

        self._status = self._run(S.values.copy(), C.values.copy(),
                                 n_beats, t_sampling, tol, output.values,
                                 None, None,
                                 stim_protocol_Ist, stim_protocol_t
                                 )

        n_beats_save = kwargs.get('n_beats_save', 1)
        output = output[-n_samples_per_stim * n_beats_save - 1:]

        return output


def calculate_n_samples_per_stim(exp_cond_name, config):
    stim_period_legend_name = config["stim_period_legend_name"]
    stim_period = config["experimental_conditions"][exp_cond_name]["params"][
        stim_period_legend_name
    ]
    t_sampling = config["t_sampling"]
    n_samples_per_stim = int(np.round(stim_period / t_sampling))
    return n_samples_per_stim
