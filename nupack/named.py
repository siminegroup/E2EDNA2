from .core import RawStrand, RawComplex, Pickleable, Sparsity, PairsMatrix, MaterialContext
from .constants import complexes_from_max_size
from .utility import create_name, check_instances, from_argument, BaseResult, printable, Loadable, long_output
from .analysis import Specification, ConcentrationSolver, ensemble_pair_fractions
from typing import Tuple, NamedTuple

import numpy, pandas, decimal

################################################################################

class Strand(Pickleable, RawStrand):
    '''
    Class representing nucleic acid strand with a mandatory name
    '''
    name: str

    def __init__(self, sequence, *, name, complement=None, _fun_=None):
        '''
        Initialize from a sequence and name, and optional reverse complement
        '''
        sequence = RawStrand(sequence)
        name = create_name(name)
        _fun_(self, sequence, complement or '', name)

    def __repr__(self):
        return '<Strand %s>' % self.name
        # return 'Strand(%r, name=%r)' % (str(self), self.name)
################################################################################

class Complex(RawComplex):
    '''
    Class representing an ordered list of connected strands with an optional name
    '''
    bonus: float

    def __init__(self, strands, *, name=None, bonus=0, _fun_=None):
        '''Initialize from a list of strands and an optional name'''
        name = getattr(strands, 'name', '') if name is None else create_name(name)
        strands = check_instances(strands, Strand)
        _fun_(self, RawComplex(strands), [s.name for s in strands],
            [~s for s in strands], name, float(bonus))

    @property
    def name(self) -> str:
        '''Manual or automatic name of complex'''

    def symmetry(self) -> int:
        '''Rotational symmetry of the labeled complex'''

    @property
    def strands(self) -> Tuple[Strand, ...]:
        '''Read only list of strands'''

    def __repr__(self):
        return '<Complex %s>' % self.name
        # return 'Complex(%r, name=%r)' % (self.strands, self.name)

################################################################################

class SetSpec(NamedTuple):
    '''
    Specification for a set of complexes that will be generated given a set of strands
    '''
    max_size: int = 1
    include: tuple = ()
    exclude: tuple = ()

    def generate(self, strands, cls=tuple):
        '''Generate a list of complexes from a given set of strands and complex type'''
        complexes = set(map(cls, self.include))
        for c in complexes_from_max_size(strands, self.max_size):
            complexes.add(cls(c))
        complexes.difference_update(map(cls, self.exclude))
        return complexes

    @classmethod
    def get(cls, obj):
        if obj is None:
            return cls()
        if isinstance(obj, cls):
            return obj
        return cls(0, obj)

################################################################################

class ComplexSet:
    '''
    Class representing a set of complexes with shared strands
    '''

    def __init__(self, strands, *, complexes=SetSpec()):
        self.strands = check_instances(strands, Strand)
        spec = SetSpec.get(complexes)
        self.complexes = spec.generate(self.strands, Complex)

    @classmethod
    def union(cls, complex_sets):
        '''Set union over an iterable of instances of ComplexSet'''
        return cls(set(s for t in complex_sets for s in t.strands),
            complexes=[c for t in complex_sets for c in t.complexes])

    def __eq__(self, other):
        if not isinstance(other, ComplexSet):
            return NotImplemented
        return self.strands == other.strands and self.complexes == other.complexes

    def __hash__(self, other):
        return hash((self.strands, self.complexes))

    def __add__(self, other):
        '''Set union with another instance (see .union())'''
        return self.union([self, other])

    def __str__(self):
        return 'ComplexSet(%s, include=%s)' % (list(self.strands), list(self.complexes))

    def __repr__(self):
        return 'ComplexSet(%r, include=%r)' % (list(self.strands), list(self.complexes))

    def __iter__(self):
        return iter(self.complexes)

################################################################################

class Tube(ComplexSet):
    '''
    Class representing a set of complexes in a tube with optional strand concentrations
    '''
    def __init__(self, strands, complexes=SetSpec(), *, name):
        self.name = create_name(name)
        strands, conc = zip(*dict(strands).items())
        super().__init__(strands, complexes=complexes)
        self.concentrations = numpy.array(conc, dtype=numpy.float64)

    def __eq__(self, other):
        if not isinstance(other, Tube):
            return NotImplemented
        return ComplexSet.__eq__(self, other) and self.concentrations == other.concentrations

    def __hash__(self):
        return hash(self.name)

    @classmethod
    def union(cls, tubes):
        '''
        Set union over an iterable of instances of Tubee
        Concentrations of the same strands are added together
        '''
        d = {}
        for t in tubes:
            for s, c in zip(t.strands, t.concentrations):
                d[s] = d.get(s, 0) + c
        return cls(d, complexes=[c for t in tubes for c in t.complexes],
            name='+'.join(t.name for t in tubes))

    def __str__(self):
        return 'Tube({%s}, name=%r)' % (', '.join('{}: {}'.format(s.name, c)
            for s, c in zip(self.strands, self.concentrations)), self.name)

    def __repr__(self):
        return '<Tube %s>' % self.name
    #     return 'Tube(%r, complexes=%r, name=%r)' % (dict(zip(self.strands,
    #         self.concentrations)), list(self.complexes), self.name)

################################################################################

@printable
class Result(BaseResult, Loadable):
    '''Class holding an analysis result'''

    names = ('tubes', 'complexes')

    def __init__(self, *, tubes=None, complexes=None):
        self.tubes = {} if tubes is None else dict(tubes)
        self.complexes = {} if complexes is None else dict(complexes)

    @property
    def fields(self):
        return (self.tubes, self.complexes)

    _unicode_cols = ['Complex', 'Pfunc', '\u0394G (kcal/mol)', 'MFE (kcal/mol)', 'Ensemble size']
    _pretty_cols = ['Complex', 'Pfunc', 'dG (kcal/mol)', 'MFE (kcal/mol)', 'Ensemble size']
    _raw_cols = ['complex', 'pfunc', 'free_energy', 'mfe', 'ensemble_size']

    _formats = [
        lambda c: c.name, # complex
        '{:.4e}'.format, # pfunc
        '{:.3f}'.format, # free energy
        '{:.3f}'.format, # mfe
        lambda i: str(i) if i < 1e10 else '{:.3e}'.format(i), # count
    ]

    _tube_fmts = [lambda x: '' if numpy.isnan(x) else '{:.3e}'.format(x)]

    def complex_frame(self, *, include_null=True, pretty: bool, unicode=True):
        '''pandas DataFrame with complexes'''
        data = [[k, v.pfunc, v.free_energy, v.mfe[0].energy if v.mfe else None, v.ensemble_size]
            for k, v in self.complexes.items()]
        if pretty:
            cols = self._unicode_cols if unicode else self._pretty_cols
        else:
            cols = self._raw_cols
        df = pandas.DataFrame(data, columns=cols)

        if include_null:
            return df

        df = df[[c for c in df.columns if not df[c].isnull().all()]]
        # sort by #strands, then alphabetically but move ( to the back
        keys = [(len(c), c.name.replace('(', '~')) for c in df[df.columns[0]]]
        return df.reindex(sorted(range(len(keys)), key=keys.__getitem__)).reset_index(drop=True)

    def tube_frame(self, *, pretty: bool):
        '''
        Return a pandas DataFrame with tubes
        - pretty: use spaces and capitalization
        '''
        if len(set(t.name for t in self.tubes)) != len(self.tubes):
            raise ValueError('Tube names are not unique')

        if self.complexes:
            complexes = self.complexes
        else:
            complexes = set(x for t in self.tubes for x in t.complexes)

        dfs = []
        for t, v in self.tubes.items():
            if pretty:
                df = pandas.DataFrame([[*x, numpy.nan] for x in v.complex_concentrations.items()],
                    columns=['Complex', t.name + ' (M)', ' '] )
            else:
                df = pandas.DataFrame(list(v.complex_concentrations.items()),
                    columns=['complex', t.name])
            # Sort by concentrations
            dfs.append(df.sort_values(by=df.columns[1], ascending=False).reset_index(drop=True))
        return pandas.concat(dfs, axis=1)

    def _material(self):
        is_rna = any(v.model.material=='RNA' for k, v in self.complexes.items())
        return MaterialContext(rna=is_rna)

    def __str__(self):
        items = []
        with long_output():
            if self.complexes:
                df = self.complex_frame(include_null=False, pretty=True, unicode=False)
                fmts = dict(zip(df.columns, self._formats))
                string = df.to_string(formatters=fmts, na_rep='')
                items.append('Complex results:\n' + string)
            if any(len(t.complexes) > 1 for t in self.tubes):
                df = self.tube_frame(pretty=True)
                f = {k: self._tube_fmts[0] for k in df.columns}
                f[df.columns[0]] = lambda c: c.name
                items.append('Concentration results:\n' + df.to_string(formatters=f, na_rep='', index=False))
        return '\n'.join(items)

    def _repr_html_(self):
        '''Return html string for output in notebooks'''
        items = []
        with long_output():
            if self.complexes:
                df = self.complex_frame(include_null=False, pretty=True, unicode=True)
                f = dict(zip(df.columns, self._formats))
                string = df.to_html(formatters=f, na_rep='', index=False)
                items.append('<b>Complex results:</b> ' + string)
            if any(len(t.complexes) > 1 for t in self.tubes):
                df = self.tube_frame(pretty=True)
                f = {k: self._tube_fmts[0] for k in df.columns}
                f[df.columns[0]] = lambda c: c.name
                string = df.to_html(formatters=f, na_rep='', index=False)
                items.append('<b>Concentration results:</b> ' + string)
        return ' '.join(items)

################################################################################

class TubeResult:
    '''
    Contains two fields:
    - complex_concentrations: a dict from complex to its equilibrium concentration
    - ensemble_pair_fractions: optionally, the equilibrium strand pairing matrix
    '''

    def __init__(self, complex_concentrations, *, tube, result=None):
        self.complex_concentrations = complex_concentrations

        if result is None or any(result[x].pairs is None for x in tube.complexes):
            self.ensemble_pair_fractions = None
        else:
            raw = ensemble_pair_fractions(tube.strands,
                {x: (complex_concentrations[x], result[x].pairs.to_sparse()) for x in tube.complexes})
            self.ensemble_pair_fractions = PairsMatrix(raw)

    def __getitem__(self, complex):
        '''Get complex concentration by complex or complex name'''
        if not isinstance(complex, str):
            return self.complex_concentrations[complex]
        try:
            return next(v for k, v
                in self.complex_concentrations.items() if k == complex)
        except StopIteration:
            raise KeyError(complex) from None

    def __str__(self):
        return str(self.complex_concentrations)

################################################################################

class Options(NamedTuple):
    '''Analysis options which apply to every complex'''

    num_sample: int = 0
    energy_gap: float = 0
    sparsity_fraction: float = 1.0
    sparsity_threshold: float = 0.0
    threads: int = 1
    cache_bytes: float = 0.0

    from_argument = classmethod(from_argument)

################################################################################

def complex_analysis(complexes, model, *, compute, options=None):
    '''
    Analyze a set of complexes for their complex ensemble thermodynamic properties
    - complexes: a ComplexSet, Tube, or list of complexes
    - compute: a list of physical quantities to compute [pfunc, mfe, pairs, sample, subopt, ensemble_size]
    - options: a dict or instance of `Options`: the options to apply to each calculation
    - model: the free energy model to use
    '''
    complexes = check_instances(set(complexes), Complex)
    spec = Specification(model=model)
    options = Options.from_argument(options)

    for c in complexes:
        for task in compute:
            if task == 'pfunc':
                spec.pfunc(c)
            elif task == 'mfe':
                spec.mfe(c)
            elif task == 'pairs':
                spec.pairs(c, sparsity=Sparsity(
                    fraction=options.sparsity_fraction, threshold=options.sparsity_threshold))
            elif task == 'sample':
                spec.sample(c, number=options.num_sample)
            elif task == 'subopt':
                spec.subopt(c, gap=options.energy_gap)
            elif task == 'ensemble_size':
                spec.ensemble_size(c)
            else:
                raise KeyError('%r is not a valid analysis task' % task)
    raw = spec.compute()

    output = {}
    for c in complexes:
        v = raw.get(RawComplex(c))
        if v is None:
            continue
        sym = c.symmetry()
        # now have to fix up all this distinguishability ... stuff
        replace = {}

        if v.mfe is not None and c.bonus != 0:
            replace['mfe'] = [s._replace(energy=s.energy + c.bonus) for s in v.mfe]
        if v.subopt is not None and c.bonus != 0:
            replace['subopt'] = [s._replace(energy=s.energy + c.bonus) for s in v.subopt]
        if v.mfe_stack is not None:
            replace['mfe_stack'] = v.mfe_stack + c.bonus
        if v.ensemble_size is not None:
            replace['ensemble_size'] = v.ensemble_size // sym
        if v.pfunc is not None:
            factor = numpy.exp(-v.model.beta * c.bonus)
            replace['pfunc'] = v.pfunc / sym * decimal.Decimal(factor)
        if v.free_energy is not None:
            replace['free_energy'] = v.free_energy + numpy.log(float(sym)) / v.model.beta + c.bonus

        output[c] = v._replace(**replace)

    return Result(complexes=output)

################################################################################

def tube_analysis(tubes, model, *, compute=(), **kws):
    '''
    Analyze a set of tubes for their complex ensemble and concentration properties
    '''
    compute = ['pfunc'] + list(compute)
    result = complex_analysis(Tube.union(tubes), compute=compute, model=model, **kws)
    for t in tubes:
        solver = ConcentrationSolver(t.strands,
            {k: result[k] for k in t.complexes}, distinguishable=False)
        solved = solver.compute(t.concentrations)
        conc = {}
        for k, v in solved.complex_concentrations.items():
            conc[k] = v
        # assert all(c is not None for c in conc)
        result.tubes[t] = TubeResult(conc, tube=t, result=result)
    return result

################################################################################

def complex_concentrations(tube, data, concentrations=None):
    '''Calculate the equilibrium concentrations for a given tube using already calculated complex results'''
    if concentrations is None:
        conc = tube.concentrations
    else:
        conc = [concentrations[strand] for strand in tube.strands]
    name = getattr(tube, 'name', 'tube')
    tube = Tube(zip(tube.strands, conc), complexes=tube.complexes, name=name)

    solver = ConcentrationSolver(tube.strands,
        {k: data[k] for k in tube.complexes}, distinguishable=False)

    res = solver.compute(tube.concentrations)
    return Result(tubes={tube: TubeResult(res.complex_concentrations, tube=tube)})

################################################################################
