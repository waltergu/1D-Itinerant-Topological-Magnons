import numpy as np
import HamiltonianPy.ED as ED
import HamiltonianPy.Beta.TrED as TrED
from HamiltonianPy import *
from .config import *

__all__=['edconstruct','tredconstruct']

def edconstruct(parameters,basis,lattice,terms,**karg):
    config=IDFConfig(priority=DEFAULT_FERMIONIC_PRIORITY,pids=lattice.pids,map=idfmap)
    ed=ED.FED(
        dlog=       'log',
        din=        'data',
        dout=       'result/ed',
        name=       '%s_%s_%s'%(name,lattice.name,repr(basis)),
        parameters= parameters,
        map=        parametermap,
        sectors=    [basis],
        lattice=    lattice,
        config=     config,
        terms=      [term(**parameters) for term in terms],
        dtype=      np.complex128
        )
    return ed

def tredconstruct(parameters,basis,lattice,terms,**karg):
    config=IDFConfig(priority=DEFAULT_FERMIONIC_PRIORITY,pids=lattice.pids,map=idfmap)
    tred=TrED.TrFED(
        dlog=       'log',
        din=        'data',
        dout=       'result/ed',
        name=       '%s_%s_%s'%(name,lattice.name,repr(basis)),
        parameters= parameters,
        map=        parametermap,
        basis=      basis,
        lattice=    lattice,
        config=     config,
        terms=      [term(**parameters) for term in terms],
        dtype=      np.complex128
        )
    return tred
