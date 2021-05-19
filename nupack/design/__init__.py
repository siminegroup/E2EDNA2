import json, numpy, pandas
from collections import defaultdict
from decimal import Decimal
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor

from ..utility import check_instance, check_instances, Loadable, long_output, \
    printable, NamedTuple, from_argument, BaseResult, create_name, mappable, Rules
from ..model import Model
from ..core import Sequence, JSON, Domain, Structure, TargetComplex, TargetStrand, RawComplex, RawStrand, MaterialContext
from ..named import Tube, Complex, Strand, TubeResult, SetSpec, Result as AnalysisResult
from ..constants import complexes_from_max_size
from ..analysis import ComplexResult

from . import results, components

try:
    from .gui import GraphicOptimization as Optimization
except ImportError:
    from .trials import Optimization

from .core import Specification
from .components import Timer, TimeInterval, MutationInterval, \
    WriteToFileCheckpoint, StopCondition

################################################################################

@mappable
class TargetTube:
    def __init__(self, on_targets, off_targets=SetSpec(), *, name):
        '''
        Initialize tube from:
        - on_targets: a dict of TargetComplex to concentration (float)
        - off_targets: a dict containing optional keys:
            - include: an additional iterable of TargetComplex to include as offtargets
            - exclude: an additional iterable of TargetComplex to exclude
            - max_size (default 1): if >0, all non-on-target complexes of up to this size will be included as offtargets
        - name: name of the tube (str)
        '''
        self.strands = tuple(set(s for t in on_targets.keys() for s in t.strands))
        self.concentrations = numpy.zeros(len(self.strands))
        self.on_targets = dict(on_targets)
        self.name = create_name(name)

        for t, conc in on_targets.items():
            if not t.structure:
                raise ValueError('Target complex %r was not given a structure' % t)
            for s in t.strands:
                self.concentrations[self.strands.index(s)] += float(conc)

        off = SetSpec.get(off_targets).generate(self.strands, TargetComplex)
        self.off_targets = off.difference(self.on_targets)

    def __str__(self):
        return 'TargetTube({%s}, name=%r)' % (', '.join('%s: %.2e' % (k.name, v)
            for k, v in zip(self.on_targets, self.concentrations)), self.name)

    def __repr__(self):
        return '<TargetTube %s>' % self.name

    @property
    def complexes(self):
        '''All complexes in the tube (read-only)'''
        return self.off_targets.union(self.on_targets)

    def replace(self, rules):
        self.off_targets = rules.map(self.off_targets)
        self.strands = rules.map(self.strands)
        self.on_targets = {rules(k): v for k, v in self.on_targets.items()}

    def __lt__(self, other):
        if isinstance(other, TargetTube):
            return self.name < other.name
        return NotImplemented

    def __eq__(self, other):
        if not isinstance(other, TargetTube):
            return NotImplemented
        return (self.name, self.strands, self.on_targets, self.off_targets) == \
            (other.name, other.strands, other.on_targets, other.off_targets) and \
            numpy.array_equal(self.concentrations, other.concentrations)

    def __hash__(self):
        return hash(self.name)

################################################################################

def _default_format():
    return '{}'

def _get_name(x):
    return x.name

def _style(table, names, html, **formats):
    '''Format each column by the given formats'''
    if names is not None:
        table = table[table.columns[:len(names)]]
        formats.update({k: formats[v] for k, v in zip(names, table.columns) if v in formats})
        table.columns = names

    if html:
        return table.to_html(formatters=formats, index=False, na_rep='')
    else:
        return table.to_string(formatters=formats, index=False, na_rep='')

################################################################################

@mappable
class Array:
    '''
    Class representing essemtially a pandas DataFrame with a MultiIndex key
    of [tube, complex, strand, domain]. Access that DataFrame as `.named()`.

    The actual table is not held as a MultiIndex for a simpler implementation.
    '''
    table: pandas.DataFrame


    def __init__(self, tubes_or_complexes):
        '''Initialize from a list of tubes'''
        if all(isinstance(t, TargetTube) for t in tubes_or_complexes):
            rows = sorted(set((d, s, c, t) for t in tubes_or_complexes
                for c in t.on_targets for s in c.strands for d in s.domains))
            self.axes = ('domain', 'strand', 'complex', 'tube')
        elif all(isinstance(t, TargetComplex) for t in tubes_or_complexes):
            rows = sorted(set((d, s, c) for c in tubes_or_complexes for s in c.strands for d in s.domains))
            self.axes = ('domain', 'strand', 'complex')
        else:
            raise TypeError('Expected list of TargetTube or list of TargetComplex')
        self.table = pandas.DataFrame(rows, columns=self.axes)

    def replace(self, rules):
        self.table = self.table.applymap(lambda x: (
            x.substitute(rules) if hasattr(x, 'substitute') else x))

    def named(self, upper=False):
        '''Return an equivalent DataFrame using a MultiIndex'''
        df = self.table.copy()
        for a in self.axes:
            df[a] = df[a].map(_get_name)
        if upper:
            df.columns = [s[0].upper() + s[1:] for s in df.columns]
        return df

    def __str__(self):
        with long_output():
            return _style(self.named(True), None, False, weight='{:#.3g}'.format)

    def _repr_html_(self):
        with long_output():
            return _style(self.named(True), None, True, weight='{:#.3g}'.format)

    def __getattr__(self, key):
        if key.startswith('_'):
            raise AttributeError(key)
        return getattr(self.table, key)

    def mask(self, key):
        '''Return rows which match the specified key'''
        if isinstance(key, tuple):
            key = key + (len(self.axes) - len(key)) * (slice(None),)
        else:
            key = (key,) + (len(self.axes) - 1) * (slice(None),)
        mask = numpy.full(len(self.table), True)
        for a, k in zip(self.axes, key):
            if isinstance(k, slice):
                bools = numpy.full(len(self.table), False)
                bools[k] = True
                mask &= bools
            elif isinstance(k, str):
                mask &= [f.name == k for f in self.table[a]]
            else:
                mask &= (self.table[a] == k)
        return mask

################################################################################

@printable
class Weights(Array):
    '''
    Class representing objective weights on domains, strands, complexes, and tubes
    '''

    def __init__(self, tubes_or_complexes):
        '''Initialize to 1s from a list of tubes or list of complexes'''
        super().__init__(tubes_or_complexes)
        self.table['weight'] = 1.0

    def __getitem__(self, key):
        '''Get weights matching a specified key'''
        return self.table.loc[self.mask(key), 'weight']

    def __setitem__(self, key, value):
        '''Set weights matching a key to value'''
    #     assert numpy.all(value >= 0), 'Weights must be non-negative'
        self.table.loc[self.mask(key), 'weight'] = value

################################################################################

@mappable
class HardConstraint(ABC):
    '''Abstract base class for a hard constraint'''

    @abstractmethod
    def apply_hard(self, design, spec, lookup):
        pass

    @abstractmethod
    def replace(self, rules):
        pass

################################################################################

@mappable
class SoftConstraint(ABC):
    '''Abstract base class for a soft constraint'''

    @abstractmethod
    def apply_soft(self, design, spec, lookup):
        pass

    @abstractmethod
    def replace(self, rules):
        pass

################################################################################

class _DualConstraint:
    def __init__(self, left, right):
        '''
        Initialize from matching list of domains:
        - left: first list of domains
        - right: second list of domains of equal total length
        '''
        self.regions = (check_instances(left, Domain), check_instances(right, Domain))
        lengths = [sum(map(len, r)) for r in self.regions]
        assert lengths[0] == lengths[1], 'Total region lengths %r should be equal' % lengths

    def replace(self, rules):
        self.regions = tuple(map(rules.map, self.regions))

################################################################################

class Match(_DualConstraint, HardConstraint):
    """Add constraint making left and right lists of domains identical

    Args:
        left (list[Domain]): list of domains to constrain
        right (list[Domain]): list of domains to constrain (of equal total length)
    """

    def apply_hard(self, design, spec, lookup):
        spec.add_match_constraint(*[tuple(map(lookup, r)) for r in self.regions])

################################################################################

class Complementarity(_DualConstraint, HardConstraint):
    """Add constraint making left and right lists of domains complementary

    Args:
        left (list[Domain]): list of domains to constrain
        right (list[Domain]): list of domains to constrain (of equal total length)
    """
    def __init__(self, left=None, right=None, *, wobble_mutations=True):
        super().__init__(left, right)
        self.wobble_mutations = bool(wobble_mutations)

    def apply_hard(self, design, spec, lookup):
        spec.add_complementarity_constraint(*map(lookup, self.regions), allow_wobble=self.wobble_mutations)

################################################################################

class Similarity(HardConstraint, SoftConstraint):
    """Make a design element match some fraction of a reference sequence

    Args:
        domains (list[Domain]): list of domains to constrain
        reference (str): a degenerate base sequence
        limits (tuple(float, float)): min and max fraction that must be matched
        weight (float): weight to apply, if using as soft constraint
    """
    kind = 'similarity'
    name = 'Soft constraints: similarity'

    def __init__(self, domains, reference, *, limits, weight=1.0):
        self.domains = check_instances(domains, Domain)
        self.reference = Sequence(reference)
        if len(self.reference) != sum(map(len, self.domains)):
            raise ValueError('Similarity(): target and reference sequence lengths are different: %d != %d' %
                (len(self.reference), sum(map(len, self.domains))))
        self.limits = limits
        self.weight = weight

    def replace(self, rules):
        self.domains = rules.map(self.domains)

    def apply_soft(self, design, spec, lookup):
        spec.add_similarity_objective(
            tuple(map(lookup, self.domains)),
            str(self.reference), limits=self.limits, weight=self.weight)

    def apply_hard(self, design, spec, lookup):
        assert self.weight == 1, 'Weight must be unspecified for hard constraint'
        spec.add_similarity_constraint(
            tuple(map(lookup, self.domains)),
            str(self.reference), self.limits)

################################################################################

class Library(HardConstraint):
    """Constrain the concatenation of domains to be drawn from some
            alternative concatenation of library sequences

    Args:

        domains list[Domain]: domains to be concatenated in the constraint
        catalog (list[list[Sequence]]): the list of libraries whose sequences will be applied to the constraint
    """
    def __init__(self, domains, catalog):
        self.domains = check_instances(domains, Domain)
        self.catalog = [[Sequence(s) for s in c] for c in catalog]
        for c in self.catalog:
            lengths = set(map(len, c))
            if len(lengths) != 1:
                raise ValueError('Library contains concatenated sequences of different lengths %s' % lengths)
        t = sum(map(len, self.domains))
        s = sum(len(c[0]) for c in self.catalog)
        if t != s:
            raise ValueError('Source length %d does not match target length %d' % (s, t))

    def replace(self, rules):
        self.domains = rules.map(self.domains)

    def apply_hard(self, design, spec, lookup):
        names = []
        for c in self.catalog:
            name = 'library-%d' % len(spec.libraries)
            spec.add_library(name, [str(s) for s in c])
            names.append(name)
        spec.add_library_constraint(tuple(map(lookup, self.domains)), names)

################################################################################

class Window(HardConstraint):
    """Constrain the concatenation of domains to be drawn as a window
        from one of the sources

    Args:
        domains (list[Domain]): concatenated domains to constrain
        sources (list[Sequence]): list of possible source sequences

    Raises:
        ValueError: if one of the sources listed is shorter than the
            concatenation of domains and hence has no windows.
    """
    def __init__(self, domains, sources):
        self.domains = check_instances(domains, Domain)
        assert not isinstance(sources, (str, Sequence)), 'Sources should be given as a list of sequences'
        self.sources = [Sequence(s) for s in sources]

    def replace(self, rules):
        self.domains = rules.map(self.domains)

    def apply_hard(self, design, spec, lookup):
        names = []
        for s in self.sources:
            name = 'source-%d' % len(spec.sources)
            spec.add_source(name, str(s))
            names.append(name)
        spec.add_window_constraint(tuple(map(lookup, self.domains)), names)

################################################################################

class Pattern(HardConstraint, SoftConstraint):
    """Constrain given domains to exclude specified patterns

    Args:
        patterns (list[Sequence]): list of pattern sequences
        scope (list[Domain]): optional region to apply the constraint (global if not specified)
        weight (float): weight to use if using as soft constraint
    """
    kind = 'pattern'
    name = 'Soft constraints: pattern prevention'

    def __init__(self, patterns, *, scope=None, weight=1.0):
        self.patterns = [Sequence(p) for p in patterns]
        self.scope = scope
        self.weight = float(weight)

    def replace(self, rules):
        self.scope = self.scope and rules.map(self.scope)

    def apply_soft(self, design, spec, lookup):
        names = tuple(map(lookup, self.scope)) if self.scope else None
        spec.add_pattern_objective(self.patterns, names=names, weight=self.weight)

    def apply_hard(self, design, spec, lookup):
        assert self.weight == 1, 'Weight must be unspecified for hard constraint'
        scope = None if self.scope is None else tuple(map(lookup, self.scope))
        spec.add_pattern_constraints(self.patterns, scope)

################################################################################

class Diversity(HardConstraint):
    """Constrain the minimum diversity (minimum_nucleotide_types) per
    windows of word_length within each domain/strand in names

    Args:
        word (int): word length to consider
        diversity (int): number of different nucleotide types that must be present in each word
        scope (list[Domain]): optional region to apply the constraint (global if not specified)
    """

    def __init__(self, word, types, *, scope=None):
        self.word = int(word)
        self.types = int(types)
        self.scope = scope

    def replace(self, rules):
        self.scope = self.scope and rules.map(self.scope)

    def apply_hard(self, design, spec, lookup):
        names=tuple(map(lookup, self.scope)) if self.scope else None
        spec.add_diversity_constraints(self.word, self.types, names=names)

################################################################################

class SSM(SoftConstraint):
    """Degree of sequence symmetry violation for a set of complexes

    Args:
        complexes (list[TargetComplex] or None): list of complexes to constrain (global if None)
        word (int): SSM word size
        weight (float): weight to apply to the SSM objective
    """

    kind = 'sequence_symmetry'
    name = 'Soft constraints: sequence symmetry'

    def __init__(self, *, scope=None, word=4, weight=1.0):
        if scope is None:
            self.scope = None
        else:
            self.scope = check_instances(scope, TargetComplex)
        self.word = int(word)
        self.weight = float(weight)

    def replace(self, rules):
        self.scope = self.scope and rules.map(self.scope)

    def apply_soft(self, design, spec, lookup):
        if self.scope is None:
            scope = design.on_targets
        else:
            scope = self.scope
        spec.add_SSM_objective(weight=self.weight,
            word_size=self.word, names=tuple(map(lookup, scope)))

################################################################################

class EnergyMatch(SoftConstraint):
    """Constrain difference in duplex energy between sets of domains

    Args:
        domains (list[Domain]): list of domains to constrain
        energy_ref (float): constrain duplex energies toward this reference (=mean if not given)
        weight (float): weight to apply to the objective
    """

    kind = 'energy difference'
    name = 'Soft constraints: energy difference'

    def __init__(self, domains, *, energy_ref=None, weight=1.0):
        self.domains = check_instances(domains, Domain)
        self.energy_ref = energy_ref
        self.weight = float(weight)

    def replace(self, rules):
        self.domains = rules.map(self.domains)

    def apply_soft(self, design, spec, lookup):
        spec.add_energy_equalization_objective(energy=self.energy_ref,
            weight=self.weight, names=tuple(map(lookup, self.domains)))

################################################################################

class Options(NamedTuple):
    '''Algorithm options for design'''

    seed: int = 0
    f_stop: float = 0.02
    f_passive: float = 0.01
    H_split: int = -1 # defaults to 2 for RNA, 3 for DNA
    N_split: int = 12
    f_split: float = 0.99
    f_stringent: float = 0.99
    dG_clamp: float = -20
    M_bad: int = 300
    M_reseed: int = 50
    M_reopt: int = 3
    f_redecomp: float = 0.03
    f_refocus: float = 0.03
    f_sparse: float = 1e-05
    slowdown: float = 0
    log: str = ''
    decomposition_log: str = ''
    thermo_log: str = ''
    time_analysis: bool = False

    @classmethod
    def from_argument(cls, options):
        out = from_argument(cls, options)
        kws = {}
        for k, v in cls.__annotations__.items():
            try:
                kws[k] = v(getattr(out, k))
            except TypeError as e:
                raise TypeError('Incorrect type for keyword %s' % k) from e
        return cls(**kws)

################################################################################

@printable
class Mapping(BaseResult):
    '''Class for mapping from undesigned to designed entities'''

    names = ('tubes', 'complexes', 'strands', 'domains')

    def __init__(self, *, rna=True):
        self.tubes = {}
        self.complexes = {}
        self.strands = {}
        self.domains = {}
        self.is_rna = bool(rna)

    @property
    def fields(self):
        return (self.tubes, self.complexes, self.strands, self.domains)

    def domain_table(self, complements=True):
        return pandas.DataFrame([(k.name, str(v))
            for k, v in self.domains.items() if complements or not k.name.endswith('*')],
            columns=('Domain', 'Sequence'))

    def rules(self):
        return self.domains

    def __call__(self, key):
        for x in self.fields:
            y = x.get(key)
            if y is not None:
                return y
        raise KeyError(key)

    def get(self, key, default=None):
        '''Lookup a given domain'''
        return self.domains.get(key, default)

    def strand_table(self):
        return pandas.DataFrame([(k.name, str(v)) for k, v in self.strands.items()],
            columns=('Strand', 'Sequence'))

    def _fmt(self, html, x):
        if html:
            h = x.to_html(formatters={'Sequence': '1@@@@{}2@@@@'.format}, index=False)
            return h.replace('1@@@@', '<pre>').replace('2@@@@', '</pre>')
        else:
            return x.to_string(index=False)

    def __str__(self):
        with MaterialContext(self.is_rna):
            with long_output():
                return 'Domain results:\n{}\n\nStrand results:\n{}'.format(
                    self._fmt(False, self.domain_table()),
                    self._fmt(False, self.strand_table()))

    def _repr_html_(self):
        with MaterialContext(self.is_rna):
            with long_output():
                return '<b>Domain results</b>:\n{}<b>Strand results</b>:\n{}'.format(
                    self._fmt(True, self.domain_table()),
                    self._fmt(True, self.strand_table()))

################################################################################

@printable
class ComplexDefect(NamedTuple):
    '''Information on a complex's defect independent of its concentration'''
    defect: float
    normalized: float

################################################################################

@printable
class TubeDefect(NamedTuple):
    '''Information on a tube's defect'''
    defect: float
    normalized: float
    complexes: dict

################################################################################

@printable
class TargetDefect(NamedTuple):
    defect: float
    normalized: float
    concentration: float
    structural: float

################################################################################

@printable
class Defects:
    '''Complete breakdown of a design's resultant defect by tube and complex'''

    def __init__(self, defects, weighted, complexes, tubes, success, soft_constraints):
        objs = list(zip(['ensemble_defect'] + [c.kind for c in soft_constraints],
            weighted, defects, ['Weighted ensemble defect'] + [c.name for c in soft_constraints]))
        self.objectives = pandas.DataFrame(objs, columns=['kind', 'weighted', 'defect', 'name'])
        self.success = bool(success)

        self.complexes = pandas.DataFrame([(c.name, v.defect, v.normalized, c)
            for c, v in complexes.items() if len(c.structure)], columns=['complex_name', 'defect', 'normalized', 'complex'])

        self.tubes = pandas.DataFrame([(t.name, v.defect, v.normalized, t) for t, v in tubes.items()],
            columns=['tube_name', 'defect', 'normalized', 'tube'])

        self.tube_complexes = pandas.DataFrame([
            (t.name, c.name, r.structural, r.concentration, r.defect, t, c)
            for t, v in tubes.items() for c, r in v.complexes.items() if c in t.on_targets],
            columns=['tube_name', 'complex_name', 'structural', 'concentration', 'total', 'tube', 'complex'])

    @property
    def ensemble_defect(self):
        return self.objectives.defect[0]

    @property
    def weighted_ensemble_defect(self):
        return self.objectives.weighted[0]

    def _is_complex_design(self):
        return numpy.all(self.tube_complexes.concentration == 0)

    def _formatted(self, html, lns):
        g = '{:#.3g}'.format
        f = '{:#.3g}'.format

        # s =  [lns[0], _style(pandas.DataFrame([[self.ensemble_defect, self.weighted_ensemble_defect]]),
        #    ['Ensemble defect', 'Weighted ensemble defect'], html)]
        # s = [lns[0], g(self.ensemble_defect)]

        df = self.objectives.groupby('name').sum().reset_index()
        df = df.sort_values('name', key=lambda x : x != 'Weighted ensemble defect')
        if len(df) > 1:
            df.loc[len(df)] = ['Total', df.weighted.sum(), 0]
        s = [lns[0], _style(df, ['Objective type', 'Value'], html, weighted=f)]
        s += [lns[1], f(self.ensemble_defect)]
        # if len(self.objectives) > 1:
            # df = pandas.DataFrame([['Objective function + soft constraints', self.objectives.weighted.sum()]])
            # s += [lns[2], _style(df, ['Augmented objective function', 'Weighted ensemble defect + weighted penalties'], html)]

        s += [lns[2], _style(self.complexes, ['Complex', 'Complex defect (nt)', 'Normalized complex defect'],
            html, defect=g, normalized=lambda x: '{:#.3g}'.format(float(x)))]

        if not self._is_complex_design():
            s += [lns[3], _style(self.tubes, ['Tube', 'Tube defect (M)', 'Normalized tube defect'], html, defect=g, normalized=f)]
            s += [lns[4], _style(self.tube_complexes,
                ['Tube', 'On-target complex', 'Structural defect (M)', 'Concentration defect (M)', 'Total defect (M)'],
                html, concentration=g, total=g, structural=g)]

        return ''.join(s)

    def __str__(self):
        with long_output():
            return self._formatted(False, [
                'Objective function:\n', '\n\nEnsemble defect: ', '\n\n', '\n\nOn-target complex defects:\n',
                '\n\nTube defects:\n', '\nComplex contributions to tube defects:\n'
            ])

    def _repr_html_(self):
        with long_output():
            return self._formatted(True, [
                '<b>Objective function</b>:', '<b>Ensemble defect</b>: ', '<br><br><b>On-target complex defects:</b>',
                '<b>Tube defects:</b>', '<b>Complex contributions to tube defects:</b>'
            ])

################################################################################

class TargetConcentration(NamedTuple):
    actual: float
    target: float

################################################################################

class Concentrations:
    def __init__(self, tubes):
        self.table = pandas.DataFrame([(t.name, c.name, r.actual, r.target, c.nt(), t, c)
            for t, v in tubes.items() for c, r in v.items()],
            columns=['tube_name', 'complex_name', 'concentration', 'target_concentration', 'nucleotides', 'tube', 'complex'])

    def _split(self, html):
        on = self.table[self.table.target_concentration != 0]
        df = self.table.copy()
        # df['mass'] = df.target_concentration * df.nucleotides
        df['thr'] = 1e-2 * df.groupby('tube').transform('max').concentration
        off = df[(df.concentration >= df.thr) & (df.target_concentration == 0)].copy()
        empty = '\u2014' if html else '-'
        for t in set(df.tube).difference(off.tube):
            off.loc[len(off)] = [t.name, empty, float('nan'), 0, 0, t, None, 0]
        f = '{:#.3g}'.format
        return (_style(on, ['Tube', 'Complex', 'Concentration (M)', 'Target concentration (M)'], html, concentration=f, target_concentration=f),
                _style(off, ['Tube', 'Complex', 'Concentration (M)'], html, concentration=lambda x: empty if numpy.isnan(x) else f(x)))

    def __str__(self):
        with long_output():
            return 'On-target complex concentrations:\n%s\n\nSignificant off-target complex concentrations (>= 1%% max complex concentration in tube):\n%s' % self._split(False)

    def _repr_html_(self):
        with long_output():
            return '<b>On-target complex concentrations</b>:\n%s\n\n<b>Significant off-target complex concentrations (\u2265 1%% max complex concentration in tube)</b>:\n%s' % self._split(True)

################################################################################

@printable
class Result(NamedTuple):
    to_analysis: Mapping
    defects: Defects
    concentrations: Concentrations
    analysis_result: AnalysisResult
    stats: dict
    raw: 'SingleResult'
    design: 'Design'

    @property
    def ensemble_defect(self):
        return self.defects.ensemble_defect

    @property
    def domains(self):
        return self.to_analysis.domains

    def evaluate_with(self, domains=None, tubes=None, model=None,
        soft_constraints=None, defect_weights=None, objective_weight=None):
        '''
        Evaluate the same design with a differently specified set of tubes,
        weights, soft_constraints, defect_weights, or objective_weight
        '''
        # Apply non-domain customizations
        new = self.design.apply(domains=self.domains.values(), tubes=tubes, model=model,
            soft_constraints=soft_constraints, defect_weights=defect_weights,
            objective_weight=objective_weight)
        # Apply manually chosen domains
        if domains is not None:
            new = new.apply(domains)
        return new.evaluate()

    def _present(self):
        if self.defects._is_complex_design():
            return self.to_analysis, self.defects
        else:
            return self.to_analysis, self.defects, self.concentrations

    def _repr_html_(self):
        return ''.join(x._repr_html_() for x in self._present())

    def __repr__(self):
        return '<DesignResult: defect=%f>' % self.ensemble_defect

    def __str__(self):
        return '\n\n'.join(map(str, self._present()))

Result.save = Loadable.save
Result.load = Loadable.load

################################################################################

@mappable
class Design(Loadable):
    '''
    A design for a set of tubes
    - tubes (list[TargetTube]): tubes to design
    - model (nupack.Model): thermodynamic model to use
    - hard_constraints (list): hard constraints that must be satisfied at all times during the design
    - soft_constraints (list): soft constraints that are added as a weighted penalty to the objective
    - defect_weights  (nupack.design.Weights): specification for tubes, complexes, domains, and strands
    - options (nupack.design.Options): design algorithm options
    '''

    def __init__(self, tubes, model, *, wobble_mutations=False, hard_constraints=(),
        soft_constraints=(), defect_weights=None, options=None, objective_weight=1.0):

        self.options = Options.from_argument(options)
        self.defect_weights = Weights(tubes) if defect_weights is None else defect_weights

        self.tubes = check_instances(tubes, TargetTube)

        # All domains
        self.domains = set(d for t in self.tubes for s in t.strands for d in s.domains)
        self.domains.update(~d for d in tuple(self.domains))
        # All strands
        self.strands = set(s for t in self.tubes for s in t.strands)
        # All on-target complexes
        self.on_targets = set(x for t in self.tubes for x in t.on_targets)
        # All complexes
        self.complexes = set(x for t in self.tubes for x in t.complexes)

        self.hard_constraints = check_instances(hard_constraints, HardConstraint)
        self.soft_constraints = check_instances(soft_constraints, SoftConstraint)
        self.model = check_instance(model, Model)
        self.objective_weight = float(objective_weight)
        self.wobble_mutations = bool(wobble_mutations)
        self.spec = self.build_spec()

    def replace(self, rules):
        self.defect_weights = rules(self.defect_weights)
        self.tubes = rules.map(self.tubes)
        self.domains = rules.map(self.domains)
        self.strands = rules.map(self.strands)
        self.complexes = rules.map(self.complexes)
        self.hard_constraints = rules.map(self.hard_constraints)
        self.soft_constraints = rules.map(self.soft_constraints)
        self.spec = self.build_spec()

    def apply(self, domains=None, tubes=None, model=None, hard_constraints=None,
        soft_constraints=None, defect_weights=None, objective_weight=None, wobble_mutations=None, options=None):
        '''Apply new attributes to a design, returning a new one'''
        new = Design(
            tubes=self.tubes if tubes is None else tubes,
            model=self.model if model is None else model,
            soft_constraints=self.soft_constraints if soft_constraints is None else soft_constraints,
            hard_constraints=self.hard_constraints if hard_constraints is None else hard_constraints,
            defect_weights=self.defect_weights if defect_weights is None else defect_weights,
            objective_weight=self.objective_weight if objective_weight is None else objective_weight,
            options=self.options if options is None else options,
            wobble_mutations=self.wobble_mutations if wobble_mutations is None else wobble_mutations,
        )
        if domains is not None:
            assert not isinstance(domains, dict)
            if not isinstance(domains, Rules):
                find = {d.name : d for d in self.domains}
                mapping = {find[d.name]: d for d in domains}
                mapping.update({find[(~d).name]: ~d for d in domains})
                domains = Rules(mapping)
            new.replace(domains)
        return new

    def run_one(self, restart=None, *, checkpoint_condition=None, checkpoint_handler=None):
        """create the designer from the fully specified design and run and
        return results

        Arguments:
            restart: a Result object that can be used to initialize the
                sequence state of the design.
            checkpoint_condition: A 2-argument callable that takes the current design's
                statistics and timer and outputs True if a checkpoint should be emitted
                and False otherwise, e.g. if time has elapsed a certain amount
            checkpoint_handler: A single argument callable that receives a
                Result object and decides how to process it as a checkpoint,
                e.g. printing to a file

        Returns: Result: the result of a finished design
        """
        if restart is not None:
            restart = restart.raw
        try:
            raw = self.spec(restart=restart,
                checkpoint_condition=checkpoint_condition,
                checkpoint_handler=checkpoint_handler)
        except RuntimeError as e:
            if 'Space::clone' in str(e):
                raise RuntimeError('Failed to clone the constraint space. This likely means '
                    'your design specification is impossible due to 1) base-pairing or sequence '
                    'length constraints inherent in the design or 2) user-specified constraints. '
                    'The best approach is currently to check your specification or simplify your '
                    'design until it works, allowing you to pinpoint the problem')
            raise
        return self.build_result(raw)

    def run(self, trials, restart=None, *, checkpoint_condition=None, checkpoint_handler=None):
        r, c, h = ([None] * trials if x is None else list(x) for x in
            (restart, checkpoint_condition, checkpoint_handler))
        assert len(r) == trials, 'Incorrect number of restarts'
        assert len(c) == trials, 'Incorrect number of checkpoint conditions'
        assert len(h) == trials, 'Incorrect number of checkpoint handlers'

        if trials == 1:
            return [self.run_one(restart=r[0], checkpoint_condition=c[0], checkpoint_handler=h[0])]

        assert self.options.seed == 0, 'Requested multiple trials of a deterministic design'

        with ThreadPoolExecutor(trials) as ex:
            tasks = [ex.submit(self.run_one, R, checkpoint_condition=C,
                checkpoint_handler=H) for R, C, H in zip(r, c, h)]
            return [t.result() for t in tasks]

    def launch(self, trials, *, checkpoint=None, restart=None, interval=600):
        '''
        Launch a series of trials of the design optimization in a non-blocking manner
        - trials: the number of trials
        - checkpoint: the directory to save results (defaults to nowhere)
        - restart: the directory to load starting sequences (defaults to nowhere)
        - interval: how often to request a root evaluation to save results to a file (in seconds)
        '''
        out = Optimization(self, trials=trials, checkpoint=checkpoint, interval=interval)
        out.start(restart)
        return out

    def evaluate(self, env=None):
        """
        Compute the multistate test tube ensemble defect for a set of fixed
        sequences or throw and exception if the input sequences have
        variable nucleotides. Returns a Result object.

        Arguments:
        - env (Local): a C++ compute environment
        """
        for d in self.domains:
            if d.name.endswith('*'):
                d = (~d).reverse_complement(self.wobble_mutations)
                assert d.is_canonical(), 'Domain complement {} is undetermined. Pass (complement=) when making a Domain for evaluation when wobble pairs are enabled'.format(d.name)
            else:
                assert d.is_canonical(), 'Domain {} is undetermined. You may only use evaluate() for a design where all Domains are determined'.format(d.name)

        return self.build_result(self.spec.evaluate(env))

    def build_spec(self):
        '''Rebuild current raw design specification'''
        spec = Specification(self.model, wobble_mutations=self.wobble_mutations)

        for key in self.options._fields:
            if key == 'seed':
                spec.parameters.rng_seed = self.options.seed
            elif key == 'H_split':
                val = self.options.H_split
                if val < 0:
                    val = {'RNA': 2, 'DNA': 3}.get(self.model.parameters.material)
                    if val is None:
                        print('Warning: expected material to be RNA or DNA')
                        val = 2
                spec.parameters.H_split = val
            else:
                setattr(spec.parameters, key, getattr(self.options, key))

        from nupack import config
        spec.cache_bytes_of_RAM = int(float(config.cache) * 10**9)

        self.names, self.objects = {}, {}

        def f(obj):
            '''return unique object name from object'''
            name = self.objects.setdefault(obj, obj.name)
            assert name is not None, 'Every designed object must be given a name: %r' % obj
            if self.names.setdefault(name, obj) == obj:
                return name
            raise ValueError('The name %r is used to refer to different objects %r and %r'
                % (name, self.names[name], obj))

        for d in self.domains:
            s = str((~d).reverse_complement(self.wobble_mutations) if d.name.endswith('*') else d)
            spec.add_domain(f(d), s)

        for s in self.strands:
            # print('adding strand', repr(f(s)), [f(d) for d in s.domains])
            spec.add_strand(f(s), [f(d) for d in s.domains])

        for x in self.complexes:
            # print('adding complex', repr(f(x)), [f(s) for s in x.strands])
            spec.add_complex(f(x), [f(s) for s in x.strands], structure=x.structure, bonus=x.bonus)

        for t in self.tubes:
            target_concs = {f(k): v for k, v in t.on_targets.items()}
            # print('adding tube', repr(f(t)), target_concs)
            spec.add_tube(f(t), target_concs)
            include = [[f(s) for s in x.strands] for x in t.off_targets]
            spec.add_off_targets(f(t), explicit=include)

        for c in self.hard_constraints:
            c.apply_hard(self, spec, lambda k: self._lookup(k, c))

        for d, s, c, *t, w in zip(*map(self.defect_weights.table.__getitem__,
                self.defect_weights.axes + ('weight',))):
            if w != 1:
                if t:
                    spec.add_weight(w, tube=f(t[0]), complex=f(c), strand=f(s), domain=f(d))
                else:
                    spec.add_weight(w, complex=f(c), strand=f(s), domain=f(d))

        spec.add_global_objective(weight=self.objective_weight)

        for c in self.soft_constraints:
            c.apply_soft(self, spec, lambda k: self._lookup(k, c))

        return spec

    def _lookup(self, k, constraint):
        try:
            return self.objects[k]
        except KeyError as e:
            raise KeyError('Constraint %r refers to unused entity %r' % (constraint, k)) from e

    def load_result(self, path):
        '''Load a Result object from a JSON file path'''
        return self.build_result(results.Result(json=JSON.from_path(path)))

    def build_result(self, r):
        '''Build a Result object from a results.Result raw object'''
        def g(key):
            return self.names[getattr(key, 'name', key)]

        # Fix reported ensemble defect to always be unweighted
        unweighted = [
            sum(t.normalized_defect for t in r.results[0].tubes) / len(r.results[0].tubes),
            *r.results[0].defects[1:]
        ]

        # Build Defects object from the raw outputs
        defects = Defects(
            defects=unweighted,
            weighted=r.results[0].weighted_defects,
            complexes={g(x): ComplexDefect(defect=x.defect, normalized=x.normalized_defect)
                for x in r.results[0].complexes if x.name in self.names},
            tubes={g(t): TubeDefect(float(t.defect), t.normalized_defect,
                {g(x): TargetDefect(float(x.defect), float(x.normalized_defect_contribution),
                    x.concentration_defect, x.structural_defect)
                    for x in t.complexes if x.name in self.names}) for t in r.results[0].tubes},
            soft_constraints=self.soft_constraints,
            success=r.success)

        conc = Concentrations(
            {g(t): {g(x): TargetConcentration(x.concentration, target=x.target_concentration)
                for x in t.complexes if x.name in self.names} for t in r.results[0].tubes})

        # Build map of undesigned entities to designed entities
        m = Mapping(rna=self.model.material == 'RNA')
        for k, v in r.results[0].domains.items():
            m.domains[g(k)] = Domain(Sequence(v), name=g(k).name)

        for s in self.strands:
            m.strands[s] = Strand(''.join(str(m.domains[d]) for d in s.domains), name=s.name)

        for x in self.complexes:
            m.complexes[x] = Complex([m.strands[s] for s in x.strands], name=x.name, bonus=x.bonus)

        for t in self.tubes:
            strands, concs = zip(*[(m.strands[s], c)
                for s, c in zip(t.strands, t.concentrations) if c])
            m.tubes[t] = Tube(zip(strands, concs),
                SetSpec(0, include=[m.complexes.get(x) or Complex([m.strands[s] for s in x]) for x in t.complexes]),
                name=t.name)

        # lookup for the actual complexes, unnamed -_-
        logq = {x.sequence: float(x.log_partition_function) for x in r.results[0].complexes}
        pairs = {tuple(x.sequence.strands): x.pair_probabilities.todense()
            for x in r.results[0].complexes if x.pair_probabilities.nnz}

        strands = {s.name: s for s in m.strands.values()}
        complexes = {x.name: x for x in m.complexes.values()}
        for x in r.results[0].complexes:
            if x.name not in complexes:
                complexes[x.name] = Complex([strands[s] for s in x.name.split('-')])

        # Make AnalysisResult with the quantities that have been computed
        a = AnalysisResult()
        for x in r.results[0].complexes:
            x = complexes[x.name]
            pf = Decimal(logq[RawComplex(x)]).exp() / x.symmetry()
            a.complexes[x] = ComplexResult(
                pairs=pairs.get(tuple(RawComplex(x).strands)), model=r.model,
                pfunc=pf, free_energy=-1/r.model.beta * float(pf.ln()))

        for t in r.results[0].tubes:
            a.tubes[m.tubes[g(t)]] = TubeResult({complexes[x.name]:
                x.concentration for x in t.complexes}, tube=m.tubes[g(t)])

        stats = {k: getattr(r.stats, k) for k in ['num_leaf_evaluations',
            'num_reseeds', 'design_time', 'analysis_time']}

        return Result(m, defects, conc, a, stats, r, self)

    def __getstate__(self):
        out = dict(self.__dict__)
        out.pop('spec')
        return out

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.spec = self.build_spec()

################################################################################

tube_design = Design

################################################################################

def complex_design(complexes, model, *, hard_constraints=(), soft_constraints=(),
    wobble_mutations=False, defect_weights=None, options=None):
    '''
    Design a set of complexes to adopt specified structures
    - complexes (list[TargetComplex]): complexes to design
    - model (nupack.Model): thermodynamic model to use
    - hard_constraints (list): hard constraints that must be satisfied at all times during the design
    - soft_constraints (list): soft constraints that are added as a weighted penalty to the objective
    - defect_weights  (nupack.design.Weights): specification for tubes, complexes, domains, and strands
    - options (nupack.design.Options): design algorithm options
    A tube is made for each complex, and a Result is returned
    '''
    tubes = [TargetTube({c: 1e-8}, SetSpec(0), name='Tube[%s]' % c.name)
        for c in check_instances(complexes, TargetComplex)]
    return Design(tubes, model, hard_constraints=hard_constraints, wobble_mutations=wobble_mutations,
        soft_constraints=soft_constraints, defect_weights=defect_weights, options=options)

################################################################################

def structure_design(structure, strands=None, *, model, wobble_mutations=False, options=None):
    '''
    A raw utility meant to take a structure string and return a list of domain sequence strings
    - structure(str or Structure): structure to design
    - strands(List[str]): nucleotide codes of each strand (optional)
    - model(Model): free energy model
    - options(DesignOptions): specific design options
    '''
    structure = Structure(structure)
    if strands is None:
        domains = [Domain('N' * l,  name='Domain[%d]' % i) for i, l in enumerate(structure.lengths())]
    else:
        strands = strands.split('+') if isinstance(strands, str) else strands
        domains = [Domain(s, name='Domain[%d]' % i) for i, s in enumerate(strands)]
    strands = [TargetStrand([d], name='RawStrand[%d]' % i) for i, d in enumerate(domains)]
    complexes = [TargetComplex(strands, structure=structure, name='complex')]
    design = complex_design(complexes, model=model, options=options, wobble_mutations=wobble_mutations)
    result = design.run_one()
    return [str(result.domains[d]) for d in domains]

des = structure_design

################################################################################

def structure_defect(structure, strands, *, model):
    '''
    A raw utility meant to take a structure string and evaluate the multitube defect
    - structure(str or Structure): target structure
    - strands(List[str]): list of strand sequences (fully determined)
    - model(Model): free energy model
    '''
    structure = Structure(structure)
    strands = strands.split('+') if isinstance(strands, str) else strands
    domains = [Domain(str(s),  name='Domain[%d]' % i) for i, s in enumerate(strands)]
    strands = [TargetStrand([d], name='RawStrand[%d]' % i) for i, d in enumerate(domains)]
    complexes = [TargetComplex(strands, structure=structure, name='complex')]
    design = complex_design(complexes, wobble_mutations=True, model=model)
    return design.evaluate().ensemble_defect

defect = structure_defect

################################################################################
