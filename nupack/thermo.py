from .rebind import Tuple, List, Dict, Callable, forward
from .core import Structure, LRU, PairList, RawComplex, PairsMatrix, Local
from .model import Model, Ensemble
from .utility import match

import numpy, decimal, typing

################################################################################

def _count(logq):
    return int(round(decimal.Decimal(logq).exp()))

################################################################################

@forward
class ComplexSampler:
    def __init__(self, sequences, logq, scale):
        pass

    def __call__(self, engine, n, env=1, _fun_=None):
        '''Return n samples from different sets of complexes as PairLists'''
        return _fun_(self, Local.from_object(env), engine, int(n))

################################################################################

@forward
class CachedModel:
    energy_model: Model

    def __init__(self, model=None, kind='pf', bits=None, _fun_=None):
        pf = dict(pf=True, mfe=False, count=True)[kind]
        if bits is None:
            bits = (64 if pf else 32) if model is None else model.bits

        if model is None:
            model = Model(bits)

        _fun_(self, model, return_type=match(
            k for k, v in self._metadata_.items() if v.cast(Tuple[int, bool]) == (bits, pf)))

        if kind == 'count':
            self.set_beta(0)

    def set_beta(self, beta) -> None:
        '''Set the inverse temperature'''

    def capacity(self) -> int:
        '''Return maximum length of sequence that can be calculated without cache refresh'''

################################################################################

@forward
class Cache(LRU):
    def __init__(self, nbytes, ensemble, complexity=3, types=(64,), _fun_=None):
        ensemble = Ensemble.get(ensemble)
        if ensemble == Ensemble.none_nupack3: # uses same recursions as nostacking
            ensemble = Ensemble.nostacking
        key = (int(complexity), int(ensemble), *types)
        _fun_(self, int(nbytes), return_type=match(
            k for k, v in self._metadata_.items() if v.cast(object) == key))

################################################################################

@forward
class XTensor:
    slices: Tuple[Tuple[numpy.ndarray, ...], ...]

################################################################################

@forward
class PairingAction:
    def __init__(self, function=None, _fun_=None):
        _fun_(self) if function is None else _fun_(self, function)

################################################################################

@forward
class Message:
    result: float
    sequences: Tuple[str, ...]

################################################################################

def obs(fun=None):
    return fun if isinstance(fun, PairingAction) else PairingAction(fun)

################################################################################

@forward
def dynamic_program(env, strands, models, cache, observe: Callable[[Message], None], pairing) -> float:
    '''Low-level dynamic program call expecting all arguments to be specified'''

@forward
def block(env, strands, models, cache, observe: Callable[[Message], None], pairing, _fun_):
    '''Low-level dynamic program call expecting all arguments to be specified'''
    out = _fun_(env, strands, models, cache, observe, pairing)
    result = out[0].cast(float)
    names = out[1].cast(List[str])
    values = [v if isinstance(v, XTensor) else v.cast(numpy.ndarray) for v in out[2:]]
    return result, dict(zip(names, values))

@forward
def subopt(env, gap, strands, models, cache, observe: Callable[[Message], None], pairing, print_segments=False) -> Dict[Structure, Tuple[float, float]]:
    '''Low-level dynamic program call expecting all arguments to be specified'''

@forward
def subopt_stream(env, gap, strands, models, cache, callback: Callable[[Tuple[PairList, float]], bool] , observe: Callable[[Message], None], pairing) -> float:
    '''Low-level dynamic program call expecting all arguments to be specified'''

@forward
def sample(env, n, workers, strands, models, cache, observe: Callable[[Message], None], pairing) -> Tuple[List[PairList], float, int]:
    '''Low-level dynamic program call expecting all arguments to be specified'''

@forward
def pair_probability(env, strands, models, cache, observe: Callable[[Message], None], pairing) -> Tuple[numpy.ndarray, float]:
    '''Low-level dynamic program call expecting all arguments to be specified'''

@forward
def permutations(env, max_size, strands, models, cache, observe: Callable[[Message], None], pairing) -> List[Tuple[RawComplex, float]]:
    '''Low-level dynamic program call expecting all arguments to be specified'''

################################################################################


def _call(fun, strands, **kws):
    try:
        return fun(strands=strands, **kws)
    except RuntimeError:
        raise RuntimeError('Function {} failed for strands {} and options {}'.format(fun, strands, kws))

################################################################################

def options(kind, mem, model, bits, ensemble=None, count=False):
    if ensemble is not None:
        model = model.copy()
        model.ensemble = ensemble
    models = [CachedModel(model=model, kind=kind, bits=abs(b)) for b in bits]
    if count:
        for m in models:
            m.set_beta(0)
    return dict(cache=Cache(mem, model.ensemble, 3, bits) if mem else False, models=models)

################################################################################

def compute_pf(tasks, output, **kws):
    beta = kws['models'][0].energy_model.beta

    def set_pf(o, logq):
        o['free_energy'] = logq / -beta
        o['pfunc'] = decimal.Context(prec=10).create_decimal(decimal.Decimal(logq).exp())

    for k, v in tasks.items():
        o = output[k]
        if v.sample:
            nicks = numpy.cumsum(list(map(len, map(str, k))))
            strucs, logq, n = _call(sample, k, n=v.sample, workers=1, **kws)
            set_pf(o, logq)
            o['sample'] = [Structure(s, nicks) for s in strucs]

        if v.pairs is not None:
            P, logq = _call(pair_probability, k, **kws)
            assert P.max() - 1 < 1e-6, (P.max(), k)
            assert P.min() > -1e-6, (P.min(), k)
            set_pf(o, logq)
            o['pairs'] = PairsMatrix(P, sparsity=v.pairs)

        elif v.pf_matrices:
            logq, o['pf_matrices'] = _call(block, k, **kws)
            set_pf(o, logq)

        elif v.pfunc:
            set_pf(o, _call(dynamic_program, k, **kws))

################################################################################

class StructureEnergy(typing.NamedTuple):
    structure: Structure
    energy: float
    stack_energy: float

    def __repr__(self):
        return 'StructureEnergy(%r, energy=%s, stack_energy=%s)' % (self.structure, self.energy, self.stack_energy)

    def __str__(self):
        return 'StructureEnergy(%r, energy=%.2f, stack_energy=%.2f)' % (str(self.structure), self.energy, self.stack_energy)

################################################################################

def compute_mfe(tasks, output, **kws):
    for k, v in tasks.items():
        o = output[k]
        if v.subopt is not None:
            if v.subopt[1]:
                o['mfe_stack'] = _call(subopt_stream, k, gap=v.subopt[0], callback=v.subopt[1], **kws)
            else:
                strucs = _call(subopt, k, gap=v.subopt[0], **kws)
                strucs = o['subopt'] = sorted((StructureEnergy(s, *f) for s, f in strucs.items()), key=lambda p: p.stack_energy)
                stack = o['mfe_stack'] = strucs[0].stack_energy if strucs else float('inf')
                o['mfe'] = [s for s in strucs if s.stack_energy < stack + 1e-4]
        elif v.mfe_matrices:
            o['mfe_stack'], o['mfe_matrices'] = _call(block, k, **kws)
        elif v.mfe:
            o['mfe_stack'] = _call(dynamic_program, k, **kws)

################################################################################

def compute_count(tasks, output, **kws):
    for k, v in tasks.items():
        if v.ensemble_size:
            output[k]['ensemble_size'] = _count(_call(dynamic_program, k, **kws))

################################################################################
