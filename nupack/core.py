from typing import Tuple, List, NamedTuple
import json, numpy as np, collections
from scipy.sparse import csc_matrix
from .utility import str_repr, DUp_to_dpp, check_instances, create_name, complement_name
from .rebind import forward
from .std import Vector

##########################################################################

@forward
class JSON:
    '''C++ JSON class'''

    def load(self, string=None):
        '''Load from JSON string'''

    def dump(self, indent=0) -> str:
        '''Return JSON string with given indent'''

    def to_object(self):
        '''Convert self to a python builtin object'''
        return json.loads(self.dump())

    def dump_binary(self) -> np.ndarray:
        '''Convert self to a numpy array'''

    def load_binary(self, state):
        '''Load self from a numpy array'''

    @classmethod
    def from_object(cls, obj):
        '''Convert from a python builtin object'''
        out = cls()
        out.load(json.dumps(obj))
        return out

    @classmethod
    def from_path(cls, path):
        '''Load from a file path'''
        out = cls()
        out.load_file(str(path))
        return out

    def __getstate__(self):
        return self.dump_binary()

    def __setstate__(self, state):
        self.move_from(JSON())
        self.load_binary(memoryview(state))

    def __str__(self):
        return self.dump()

##########################################################################

class Pickleable:
    def from_json(self, json, _fun_=None):
        '''Load another instance of this class from a nupack.JSON object'''
        return _fun_(json)

    def to_json(self) -> JSON:
        '''Convert self to a nupack.JSON object'''

    def __getstate__(self):
        return self.to_json()

    def __setstate__(self, state):
        self.move_from(self.from_json(state))

##########################################################################

class SimpleType:
    def __str__(self):
        return '{}({})'.format(type(self).__name__,
            ', '.join(k + '=' + repr(getattr(self, k)) for k in self.__annotations__))

    def annotated(self):
        return {k: getattr(self, k) for k in self.__annotations__}

##########################################################################

@forward
class PairList(Pickleable):
    '''
    A container with a list of base pairs `p` such that  `p[i] == j` if
    `i` is paired to `j` and `p[i] == i` if `i` is unpaired
    '''
    values: Vector

    def __init__(self, pairs, _fun_=None):
        if isinstance(pairs, PairList):
            self.copy_from(pairs)
        elif isinstance(pairs, str):
            if 'D' in pairs or 'U' in pairs:
                pairs = DUp_to_dpp(pairs)
            _fun_(self, pairs)
        else:
            _fun_(self, list(pairs)) #  np.asarray(pairs, dtype=np.uint32)

    def dp(self, nicks=()) -> str:
        '''Return equivalent dot-parens-plus string'''

    def pairlist(self):
        '''Return numpy array of pair list indices'''
        return self.view().copy()

    def view(self):
        return self.values.cast(np.ndarray)

    def matrix(self):
        '''Return square structure matrix of 0s and 1s'''
        S = np.zeros([len(self)] * 2, dtype=np.int32)
        for i, j in enumerate(self):
            S[i, j] = 1
        return S

    structure_matrix = matrix

    def pseudoknots(self) -> List[Tuple[int,int,int,int]]:
        '''Return list of i,j,k,l pseudoknots'''

    def __repr__(self):
        return 'PairList(%r)' % self.view().tolist()

    def __str__(self):
        return str(self.view())

    def __iter__(self):
        return iter(self.pairlist())

    def __getitem__(self, i):
        return self.view()[i]

    def __len__(self):
        '''Number of nucleotides, synonym of nt()'''
        return len(self.view())

    def nt(self):
        '''Number of nucleotides, synonym of len()'''
        return len(self.view())

    def __xor__(self, other) -> int:
        '''Nucleotide distance with other pair list'''

##########################################################################

def struc_distance(structure1, structure2):
    '''Return integer nucleotide distance between two structures'''
    try:
        assert structure1.nicks == structure2.nicks
    except AttributeError:
        pass
    return PairList(structure1) ^ PairList(structure2)

##########################################################################

@forward
class Structure(PairList):
    '''
    Secondary structure including base pairing and nick information
    '''

    def __init__(self, structure, nicks=None, _fun_=None):
        '''Create Structure from either a dot-parens-plus or dpp-rle string representation'''
        if structure is None:
            _fun_(self)
        elif isinstance(structure, Structure):
            self.copy_from(structure)
        elif isinstance(structure, str):
            if 'U' in structure or 'D' in structure:
                structure = DUp_to_dpp(structure)
            _fun_(self, structure)
        else:
            assert nicks is not None
            _fun_(self, structure, list(nicks))

    def lengths(self):
        '''Return lengths of each strand as a numpy array'''
        return np.diff((0,) + tuple(self.nicks()))

    def nicks(self) -> np.ndarray:
        '''A list N of zero-based indices of each base 3′ of a nick between strands'''

    def dp(self) -> str:
        '''full dpp representation of the structure'''

    def dp_rle(self) -> str:
        '''run-length encoded dpp representation of the structure'''

    def dotparensplus(self) -> str:
        '''full dpp representation of the structure'''

    def rle_dotparensplus(self) -> str:
        '''run-length encoded dpp representation of the structure'''

    def __str__(self):
        return self.dp()

    def __repr__(self):
        return "Structure('%s')" % str(self)

##########################################################################

@forward
class Local:
    '''
    Executor for serial or shared memory parallelism
    '''
    def __init__(self, threads=None, _fun_=None):
        '''Initialize with given number of threads'''
        if threads is None:
            from nupack import config
            _fun_(self, 0 if config.parallelism else 1)
        else:
            _fun_(self, int(threads))

##########################################################################

def set_sequence_type(rna: int) -> int:
    '''
    Set whether sequence should print as RNA; return previous value
    0: print as RNA (weak control)
    1: print as RNA (strong control)
    2: print as DNA (weak control)
    3: print as DNA (strong control)
    '''

class MaterialContext:
    def __init__(self, rna):
        self.rna = 1 if rna else 3
        self.backup = 0

    def __enter__(self):
        self.backup = set_sequence_type(self.rna)

    def __exit__(self, v, t, tb):
        set_sequence_type(self.backup)

##########################################################################

@forward
class Base:
    '''
    Class for a single base: ACGT or a wildcard of those bases
U is represented the same as T
    '''
    def __init__(self, letter):
        pass

    def __str__(self):
        return chr(int(self.letter()))

    def __repr__(self):
        return 'Base(%r)' % str(self)

##########################################################################

@forward
class Sequence:
    '''
    Class for an RNA or DNA sequence, possibly containing wildcards
    '''
    def __init__(self, sequence):
        '''Initialize from a string or another Sequence'''

    def __str__(self) -> str:
        '''String of the sequence'''

    def __repr__(self):
        return 'Sequence(%r)' % str(self)

    def reverse_complement(self):
        '''Return the Watson Crick reverse complement of self'''
        return ~self

    def __invert__(self):
        '''Return the Watson Crick reverse complement of self'''

    def __getitem__(self, index, _fun_=None) -> Base:
        '''Return a given base'''
        if abs(index) >= len(self):
            raise IndexError(index)
        return _fun_(self, index % len(self))

    def __contains__(self, base) -> bool:
        '''Test if a given base is contained'''

    def __iter__(self):
        return iter(map(Base, str(self)))

    def __len__(self) -> int:
        '''Number of nucleotides'''

    def __xor__(self, other) -> int:
        '''Hamming distance'''

    def nt(self) -> int:
        '''Number of nucleotides'''

    def is_canonical(self) -> bool:
        '''Return whether each base is one of ACTGU'''

    def has_wildcard(self) -> bool:
        '''Return whether any base has multiple possibilities'''

    def is_determined(self) -> bool:
        '''Return whether no base has multiple possibilities'''

    def __getstate__(self):
        return str(self)

    def __setstate__(self, state):
        self.move_from(Sequence(state))

##########################################################################

def to_sequences(sequences_or_string):
    '''Return sequence list from a delimited string or list of strings'''

##########################################################################

def seq_distance(complex1, complex2, _fun_=None):
    '''Hamming distance between 2 sequences or complexes'''
    return _fun_(to_sequences(complex1), to_sequences(complex2)).cast(int)

##########################################################################

@forward
class RawStrand(Sequence):
    '''
    Class for a strand of RNA or DNA, not containing wildcards
    '''
    def __init__(self, sequence, *, _fun_=None):
        try:
            _fun_(self, sequence)
        except TypeError:
            raise TypeError('sequence should be specified as a string or Sequence subclass')

    def __repr__(self):
        return 'RawStrand(%r)' % str(self)

    def __invert__(self):
        '''Return Watson Crick reverse complement of self'''
        return RawStrand(Sequence.__invert__(self))

##########################################################################

@forward
class Domain(Pickleable, Sequence):
    '''
    Class for a domain of RNA or DNA, uniquely named and usually containing wildcards
    '''
    name: str

    def __init__(self, sequence, *, name, complement=None, _fun_=None):
        '''
        Initialize a Domain from a single sequence and name
        - sequence (Sequence or str)
        - name (str or list)
        - complement (Sequence or str, optional)
        '''
        sequence = Sequence(sequence)
        # assert not name.endswith('*'), 'Domain name cannot end with *; use the ~ operator for complements.'
        _fun_(self, sequence, complement or '', create_name(name))

    def reverse_complement(self, include_wobble=False):
        '''Return the complement domain, optionally including wobble complements'''

    def __repr__(self):
        return '<Domain %s>' % self.name

    def __invert__(self):
        '''Return Watson Crick reverse complement of self'''

    def substitute(self, rules):
        return rules.get(self, self)

##########################################################################

@forward
class TargetStrand(Pickleable, Sequence):
    '''
    Class for a designable strand of RNA or DNA, uniquely named and containing Domains
    '''
    name: str
    domains: Tuple[Domain, ...]

    def __init__(self, domains, *, name, _fun_=None):
        '''Initialize from list of domains'''
        _fun_(self, check_instances(domains, Domain), create_name(name))

    def __iter__(self):
        return iter(self.domains)

    def __len__(self):
        '''Number of domains'''
        return len(self.domains)

    def ndomains(self):
        '''Number of domains'''
        return len(self.domains)

    def __getitem__(self, index):
        return self.domains[index]

    def __repr__(self):
        return '<TargetStrand %s>' % self.name

    def __invert__(self):
        '''Return Watson Crick reverse complement of self'''

    def substitute(self, rules):
        return rules.get(self) or TargetStrand([d.substitute(rules) for d in self.domains], name=self.name)

##########################################################################

@forward
class RawComplex(Pickleable):
    '''
    Class for an ordered complex of RNA or DNA, not containing wildcards
    '''

    def __init__(self, strands, *, _fun_=None):
        '''
        Initialize from a list of strands
        - strands: list of Sequence or RawStrand, or RawComplex or TargetComplex
        '''
        assert not isinstance(strands, (Sequence, RawStrand)), 'strands should be a list of strands. Try using RawComplex([strand])'
        if isinstance(strands, TargetComplex):
            _fun_(self, list(strands))
        else:
            _fun_(self, strands)

    @property
    def strands(self) -> List[RawStrand]:
        '''Get strands as a python list'''

    def __len__(self) -> int:
        '''Number of strands'''

    def __xor__(self, other) -> int:
        '''Hamming distance in sequences'''

    def nstrands(self):
        '''Number of strands'''
        return len(self)

    def __getitem__(self, index, _fun_=None) -> RawStrand:
        if abs(index) >= len(self):
            raise IndexError(index)
        return _fun_(self, index % len(self))

    def __iter__(self):
        return iter(self.strands)

    def symmetry(self) -> int:
        '''Rotational symmetry factor'''

    def lowest_rotation(self):
        '''Return equivalent RawComplex in the lowest rotational order'''

    def __contains__(self, strand) -> bool:
        '''Test whether a given strand is contained'''

    def nt(self) -> int:
        '''Number of nucleotides'''

    def __str__(self):
        return '+'.join(map(str, self.strands))

    def __repr__(self):
        return 'RawComplex([%s])' % ', '.join(map(repr, map(str, self.strands)))

    # def __getstate__(self):
    #     return str(self)

    # def __setstate__(self, state):
    #     self.copy_from(RawComplex(state))

##########################################################################

@forward
class TargetComplex(Pickleable):
    '''
    Class for an ordered complex of RNA or DNA, containing designable strands
    and possibly an associated target structure
    '''

    structure: Structure
    bonus: float
    strands: Tuple[TargetStrand, ...]

    def __init__(self, strands, structure=None, *, bonus=0, name=None, _fun_=None):
        '''
        strands: list of TargetStrand, or another TargetComplex
        '''
        if isinstance(strands, TargetComplex):
            _fun_(self, strands)
            return
        name = create_name(name) if name else ''
        if structure is None or isinstance(structure, str):
            structure = Structure(structure)
        if not isinstance(structure, Structure):
            raise TypeError('Expected Structure or str: got %r' % structure)
        _fun_(self, strands, name, structure, float(bonus))

    def substitute(self, rules):
        return rules.get(self) or TargetComplex([s.substitute(rules) for s in self.strands],
            structure=self.structure, bonus=self.bonus, name=self.name)

    @property
    def name(self) -> str:
        pass

    def nstrands(self):
        '''Number of strands'''
        return len(self)

    def __len__(self) -> int:
        pass

    def __getitem__(self, index, _fun_=None) -> TargetStrand:
        if abs(index) >= len(self):
            raise IndexError(index)
        return _fun_(self, index % len(self))

    def __iter__(self):
        return iter(self.strands)

    def __contains__(self, strand) -> bool:
        '''Test whehther a given TargetStrand is contained'''

    def nt(self) -> int:
        '''Number of nucleotides'''

    def __str__(self):
        return '+'.join(map(str, self.strands))

    def __repr__(self):
        return '<TargetComplex %s>' % self.name

##########################################################################

@forward
class Fenwick:
    def __len__(self):
        pass

##########################################################################

@forward
class MemoryLimit:
    length: int
    capacity: int

##########################################################################

@forward
class LRU:
    pass

################################################################################

class Sparsity(NamedTuple):
    '''
    Specification of a sparsity level for pair probability calculation
    - fraction: the maximum fraction of elements in each row to include
    - threshold: the minimum pair probability to be considered
    '''

    fraction: float = 1.0
    threshold: float = 0.0

################################################################################

@forward
class SparsePairs:
    '''
    Raw minimal representation of a sparse pair probability matrix
    - diag: value of each diagonal element
    - values: value of each offdiagonal entry (i, j) s.t. i < j
    - rows: row i of each offdiagonal entry
    - cols: column j of each offdiagonal entry
    '''
    diag: np.ndarray
    values: np.ndarray
    rows: np.ndarray
    cols: np.ndarray

@forward
def sparse_pair_matrix(A, row_size=0, threshold=0, _fun_=None):
    '''Return matrix diagonal and off-diagonal elements'''
    return _fun_(np.asfortranarray(A), row_size, threshold).cast(SparsePairs)

################################################################################

class PairsMatrix:
    '''Class representing a possibly sparse base pairing matrix'''
    def __init__(self, full, *, sparsity=Sparsity()):
        '''Initialize from full matrix and desired sparsity between 0 and 1'''
        f, t = sparsity.fraction, sparsity.threshold

        if f < 0 or f > 1:
            raise ValueError('Sparsity fraction should be in [0:1] (is {})'.format(f))
        if t < 0 or t > 1:
            raise ValueError('Sparsity threshold should be in [0:1] (is {})'.format(t))

        if f == 1 and t == 0: # no sparsification
            self.diagonal, self.values, self.rows, self.cols = [None] * 4
            self.array = full.copy()
        else:
            p = sparse_pair_matrix(full, row_size=len(full) * f, threshold=t)
            self.diagonal, self.values, self.rows, self.cols = p.diag, p.values, p.rows, p.cols
            self.array = None

    def defect(self, structure):
        '''Calculate structure defect with respect to a Structure or PairList'''
        P = self.to_array()
        return len(P) - (structure.matrix() * P).sum()

    def to_array(self):
        '''Convert to a numpy array'''
        if self.array is not None:
            return self.array.copy()
        out = np.diag(self.diagonal)
        out[self.rows, self.cols] = self.values
        out[self.cols, self.rows] = self.values
        return out

    def to_sparse(self, cls=csc_matrix):
        '''Return equivalent sparse matrix'''
        if self.array is not None:
            return cls(self.array)
        n = len(self.diagonal)
        rows = np.concatenate([np.arange(n), self.rows, self.cols])
        cols = np.concatenate([np.arange(n), self.cols, self.rows])
        vals = np.concatenate([self.diagonal, self.values, self.values])
        return cls((vals, (rows, cols)), shape=(n, n))

    def __str__(self):
        with np.printoptions(formatter={'float': '{:.4f}'.format}):
            return str(self.to_array())
