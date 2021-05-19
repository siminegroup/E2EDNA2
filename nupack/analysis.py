import functools, itertools, collections, numpy, pandas, decimal, scipy.sparse

from .rebind import Tuple, List, Dict, Union, Callable, forward
from .core import Structure, RawComplex, Sparsity, PairsMatrix
from .utility import match, str_repr, long_output, printable, NamedTuple, \
    BaseResult, from_argument, check_instance
from .model import Model, Ensemble
from .rotation import lowest_rotation
from .constants import complexes_from_max_size, as_sequences
from . import core
from .concentration import solve_complex_concentrations
from . import thermo
from .thermo import StructureEnergy

# min(2 gb, disable it if less than 2 gb)

################################################################################

def into_blocks(matrix, lengths):
    '''split a matrix into blocks with respect to the given lengths'''
    prefix = numpy.cumsum(lengths)
    return [[matrix[p1-l1:p1, p2-l2:p2] for p2, l2 in zip(prefix, lengths)]
                                        for p1, l1 in zip(prefix, lengths)]

################################################################################

def ensemble_pair_fractions(strands, complex_data):
    '''
    Calculate ensemble pair fractions given the pair probability matrices for a set of complexes
    - strands: list of RawStrand
    - complex_data: map from RawComplex to pair probability matrix (ndarray or sparse)
    Should work with either sparse or non-sparse matrices
    '''
    strands = tuple(strands)
    assert len(set(strands)) == len(strands), 'Strands should be unique'
    # out = np.full([len(strands)]*2, 0, dtype=object)
    # P = into_blocks(out, tuple(map(len, strands)))
    P = [[0 for _ in strands] for _ in strands]
    for k, (c, x) in complex_data.items():
        blocks = into_blocks(x, tuple(map(len, k.strands)))
        for i, row in zip(k.strands, blocks):
            for j, p in zip(k.strands, row):
                P[strands.index(i)][strands.index(j)] += c * p
    P = scipy.sparse.vstack(list(map(scipy.sparse.hstack, P)))
    P /= P.sum(1)
    return P

################################################################################

def pfunc(strands, model) -> Tuple[float, float]:
    '''
    Calculate partition function of a single complex
    - model(Model): free energy model to use
    - Returns a tuple of (pfunc, free_energy)
    '''
    strands = RawComplex(strands)
    res = Specification(model).pfunc(strands).compute()
    return (res[strands].pfunc, res[strands].free_energy)

################################################################################

def pairs(strands, model, *, sparsity_fraction=1, sparsity_threshold=0) -> PairsMatrix:
    '''
    Calculate equilibrium pair probabilities of a single complex
    - model(Model): free energy model to use
    - sparsity_fraction: maximum fraction of each row to keep
    - sparsity_threshold: minimum pair probability to preserve
    - Returns a PairsMatrix object
    '''
    strands = RawComplex(strands)
    spec = Specification(model).pairs(strands,
        sparsity=Sparsity(fraction=sparsity_fraction, threshold=sparsity_threshold))
    res = spec.compute()
    return res[strands].pairs

################################################################################

def mfe(strands, model) -> List[StructureEnergy]:
    '''
    Calculate MFE structures and energies for a single complex
    - model(Model): free energy model to use
    - Returns a list of StructureEnergy objects
    '''
    strands = RawComplex(strands)
    res = Specification(model).mfe(strands).compute()
    return res[strands].mfe

################################################################################

def energy(strands, structure, model) -> float:
    '''
    Calculate energy of a single structure of a complex
    - model(Model): free energy model to use
    '''
    bonus = getattr(strands, 'energy', 0.0)
    strands = RawComplex(strands)
    return model.structure_energy(strands, structure) + bonus

################################################################################

def structure_probability(strands, structure, model) -> float:
    '''
    Calculate equilibrium probability of a structure for a given complex
    - model(Model): free energy model to use
    '''
    strands = RawComplex(strands)
    res = Specification(model).pfunc(strands).compute()
    e = energy(strands, structure, res[strands].model)
    return numpy.exp(-res[strands].model.beta * (e - res[strands].free_energy))

################################################################################

def ensemble_size(strands, model):
    '''
    Calculate number of secondary structures in a given complex ensemble
    - model(Model): free energy model to use
    '''
    strands = RawComplex(strands)
    res = Specification(model).ensemble_size(strands).compute()
    return res[strands].ensemble_size

################################################################################

def subopt(strands, energy_gap, model):
    '''
    Calculate suboptimal structures falling below a specified energy gap of the MFE
    - model(Model): free energy model to use
    - energy_gap: suboptimal structure energy gap in kcal/mol
    '''
    strands = RawComplex(strands)
    res = Specification(model).subopt(strands, gap=energy_gap).compute()
    return res[strands].subopt

################################################################################

def sample(strands, num_sample, model):
    '''
    Calculate a set of num_sample Boltzmann sampled structures
    - model(Model): free energy model to use
    - num_sample(int): number of samples
    '''
    strands = RawComplex(strands)
    res = Specification(model).sample(strands, number=num_sample).compute()
    return res[strands].sample

################################################################################

class ConcentrationResult:

    def __init__(self, strands, complexes, concentrations):
        self.strands = tuple(strands)
        self.complexes = tuple(complexes)
        self.concentrations = numpy.asarray(concentrations)
        assert len(self.complexes) == len(self.concentrations)

    @property
    def complex_concentrations(self):
        return dict(zip(self.complexes, self.concentrations))

################################################################################

class ConcentrationSolver:
    def __init__(self, strands, complex_results, *, distinguishable):
        '''
          Initialize the solver from
          - strands: an ordered list of RawStrand
        - - complex_results: a dict of RawComplex to ComplexResult
          Note that each ComplexResult must contain the partition function
          If distinguishable=True, then a RawStrand entered twice
          is equivalent to saying that the first entry is labeled differently than
          the second.
        '''

        self.strands = tuple(strands)
        self.temperature = tuple(complex_results.values())[0].model.temperature

        self.complexes = {k : v.pfunc.ln()
            for k, v in complex_results.items() if v.pfunc is not None}

        self.distinguishable = bool(distinguishable)

    def compute(self, concentrations, complexes=None, as_strands=True, **kws):
        '''
        concentrations: array
            list of molarity of each strand species
        complexes: None, int, or List[List[str]]
            None: calculate concentrations in an ensemble containing all prior calculated complexes
            int: calculate concentration in an ensemble containing complexes of only up to this size
            List[List[str]]: calculate concentrations in an ensemble of only these complexes
        as_strands: bool
            Whether initial concentrations are given for each strand or for each complex
        **kws:
            Custom options for the concentration solving algorithm (see solve_complex_concentrations)
        '''
        complexes = tuple(self.complexes.keys()) if complexes is None else complexes

        if as_strands:
            if len(concentrations) != len(self.strands):
                raise ValueError('strand number %d != concentration number %d' % (len(self.strands), len(concentrations)))
        else:
            if len(concentrations) != len(complexes):
                raise ValueError('complex number %d != concentration number %d' % (len(complexes), len(concentrations)))

        logq = [self.complexes[k] for k in complexes]
        indices = [[self.strands.index(s) for s in k] for k in complexes]
        x = solve_complex_concentrations(indices, logq, concentrations,
            kelvin=self.temperature, as_strands=as_strands,
            rotational_correction=self.distinguishable, **kws)
        return ConcentrationResult(self.strands, complexes, x)


################################################################################

class Task(NamedTuple):
    '''Inputs for thermodynamic calculations for a given complex'''

    pfunc: bool = False
    mfe: int = None
    ensemble_size: bool = False
    pairs: Sparsity = None
    sample: int = 0
    subopt: float = None
    pf_matrices: bool = False
    mfe_matrices: int = None

################################################################################

class ComplexResult(NamedTuple):
    '''Output result for analysis of a single Complex, with a lot of optional fields'''
    model: Model
    # output of pfunc
    pfunc: decimal.Decimal = None
    free_energy: float = None
    pf_matrices: dict = None
    # output of mfe
    mfe_stack: float = None
    mfe_matrices: dict = None
    # output of ensemble_size
    ensemble_size: int = None
    # output of pairs
    pairs: PairsMatrix = None
    # output of subopt
    mfe: list = None
    subopt: list = None
    # output of sample
    sample: list = None

    def structure_mfe(self) -> float:
        return self.mfe[0].energy if self.mfe else float('inf')

    def __repr__(self):
        return str_repr(self)

    def __str__(self):
        return  '{%s}' % ', '.join('{}: {}'.format(k, v) for k, v
            in zip(self._fields, self) if v is not None)

################################################################################

class Specification:
    '''
    The main class for thermodynamic computation in the complex ensemble. Computations
    are scheduled via member functions like `pfunc`, while the
    computations are then performed with `compute()`. Most methods for scheduling
    computation also take a `max_size` parameter. If this is specified as `>0`,
    all complexes of the given strands up to that size will be computed. Otherwise,
    just the complex of the given strands in that order will be computed.
    '''
    default_bits = {'pfunc': [64, -32],
                    'count': [64, -32],
                    'mfe': [32]}

    def __init__(self, model):
        '''Initialize analysis for a given free energy model of type nupack.Model'''
        self.model = check_instance(model, Model)
        self.tasks = {}

    def add(self, strands, max_size=0):
        strands = RawComplex(strands)
        for s in complexes_from_max_size(strands, max_size) if max_size else (strands,):
            r = RawComplex(s)
            yield r, self.tasks.setdefault(r, Task())

    def pfunc(self, strands, max_size=0, matrices=False):
        '''Schedule computation of complex partition function and free energy'''
        for k, v in self.add(strands, max_size):
            self.tasks[k] = v._replace(pfunc=True, pf_matrices=v.pf_matrices or matrices)
        return self

    def pairs(self, strands, max_size=0, *, sparsity=Sparsity()):
        '''Schedule computation of complex equilibrium pair probability'''
        for k, v in self.add(strands, max_size):
            self.tasks[k] = v._replace(pairs=sparsity)
        return self

    def mfe(self, strands, max_size=0, *, structures=True, matrices=False):
        '''Schedule computation of minimum free energy and, optionally, matching structures'''
        if structures:
            self.subopt(strands, max_size, gap=0)
        for k, v in self.add(strands, max_size):
            self.tasks[k] = v._replace(mfe=True, mfe_matrices=v.mfe_matrices or matrices)
        return self

    def subopt(self, strands, max_size=0, *, gap, stream=None):
        '''Schedule computation of suboptimal structures'''
        for k, v in self.add(strands, max_size):
            self.tasks[k] = v._replace(subopt=(gap if
                v.subopt is None else max(v.subopt[0], gap), stream))
        return self

    def sample(self, strands, max_size=0, *, number):
        '''Schedule computation of structures sampled from Boltzmann distribution'''
        for k, v in self.add(strands, max_size):
            self.tasks[k] = v._replace(sample=v.sample + number)
        return self

    def ensemble_size(self, strands, max_size=0):
        for k, v in self.add(strands, max_size):
            self.tasks[k] = v._replace(ensemble_size=True)
        return self

    def compute(self, pairing=None, dtypes=None, gil=True):
        '''Compute all scheduled tasks and return a dict of RawComplex to ComplexResult'''
        from nupack import config
        mem = int(float(config.cache) * 10**9) # cache bytes
        env = core.Local()

        if dtypes is None:
            bits = self.default_bits
        else:
            bits = self.default_bits.copy()
            bits.update(dtypes)

        kw = dict(env=env, pairing=thermo.obs(pairing), observe=None, gil=gil)
        out = {k: {} for k in self.tasks}

        thermo.compute_pf(self.tasks, out, **kw,
            **thermo.options('pf', mem, self.model, bits['pfunc']))
        thermo.compute_mfe(self.tasks, out, **kw,
            **thermo.options('mfe', mem, self.model, bits['mfe']))
        thermo.compute_count(self.tasks, out, **kw,
            **thermo.options('pf', mem, self.model, bits['count'], count=True))

        return {k: ComplexResult(self.model, **v) for k, v in out.items()} # maybe return caches too if requested

################################################################################
