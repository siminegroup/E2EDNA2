'''
Program-wide constants and functions for dealing with nucleic acid sequences and structures
'''
from .utility import can_convert
from . import core
from . import rotation
from .rebind import forward

import re
from typing import Callable, List, Tuple, Dict
import numpy as np

ZeroCinK = float
DefaultTemperature = float
BoltzmannConstant = float
GitBranch = str
GitRevision = str
Version = str

##########################################################################

Bases = dict(DNA=tuple('ACGT'), RNA=tuple('ACGU'))
Pairs = dict(DNA=dict('AT CG GC TA'.split()), RNA=dict('AU CG GC UA'.split()))

def reverse_complement(sequence, material='DNA'):
    get = Pairs[material].__getitem__
    return ''.join(map(get, reversed(sequence)))

##########################################################################

def random_sequence(n, material='DNA', gc=0.5):
    p = 0.5 * np.asarray((1-gc, gc, gc, 1-gc), dtype=np.float64)
    get = Bases[material].__getitem__
    return ''.join(map(get, np.random.choice(4, n, p=p)))

##########################################################################

@forward
def dp_to_pairs(dpp_string) -> Tuple[int, ...]:
    pass

##########################################################################

def as_sequences(seqs):
    if isinstance(seqs, str):
        seqs = re.split('[+, ]+', seqs)
    return tuple(i for i in seqs if i)

##########################################################################

@forward
def rotational_symmetry(indices) -> int:
    '''Rotational symmetry of a list of sequences'''

@forward
def compute_necklaces(callback: Callable[[List[int]], None], *, size, n_elements) -> int:
    '''compute necklaces and pass them to callback function'''


def complexes_from_max_size(strands, max_size):
    complexes = set()

    def add_perm(l):
        names = tuple(rotation.lowest_rotation([strands[i] for i in l]))
        complexes.add(names)

    for i in range(1, max_size+1):
        compute_necklaces(add_perm, n_elements=len(strands), size=i)

    return complexes

##########################################################################

def string_exp(exponent):
    exp = np.exp(exponent, dtype=np.float128)
    if np.isfinite(exp):
        return exp
    exponent /= np.log(10.0)
    ie = int(np.floor(exponent))
    fe = 10.0 ** (exponent - ie)
    return '{}e{}'.format(fe, '%+d' % ie) if ie else str(fe)

##########################################################################

ANALYSIS_CPU_PREFACTOR = 2.5e-9 # multiply by N^3 to get approximate CPU time to evaluate single strand of length N

def unit_evaluation_costs(n_elements, lmax) -> List[int]:
    '''Blocking and naive evaluation costs for number of elements and given lmax'''


def unit_evaluation_cost_table(max_elements, timeout) -> List[List[List[int]]]:
    '''Table of unit_evaluation_costs using a timeout'''


def unit_subblock_cost(n) -> int:
    '''subblock_cost assuming every length is 1 and there are n strands'''


def evaluation_cost(lengths, lmax):
    '''
    Return non-blocking evaluation cost of partition function computation
    lengths: lengths of each strand
    lmax: maximum complex size
    '''
    cost = 0
    def callback(v):
        nonlocal cost
        cost += sum(map(lengths.__getitem__, v))**3
    for l in range(lmax):
        compute_necklaces(callback, n_elements=len(lengths), size=l+1)
    return cost

@forward
def subblock_cost(lengths) -> int:
    '''
    Cost of calculating the lengths[:-1] x lengths[1:] subblock
    Reduces to n**3 if len(lengths) is 1
    '''

# next goal: for every m necklace in a set of n distinct strands, what is the total cost of subblocks to calculate?
# only have to consider top necklaces, they contain all of the others
# so how many of each subblock type are there?
# 1-> clearly n
# 2 ->

# next goal: for every m necklace in a set of n distinct strands, what is the total cost of naive calculation? --> #necklaces * (m n)**3
# have to add this for each necklace size

def blocking_evaluation_cost(lengths, lmax):
    '''Evaluation cost assuming every subblock is cached perfectly'''
    subblocks = set()
    def callback(v):
        '''add each subsequence of v to subblocks'''
        for i in range(len(v)):
            for j in range(i+1, len(v)+1):
                subblocks.add(tuple(v)[i:j])
    for l in range(lmax):
        compute_necklaces(callback, n_elements=len(lengths), size=l+1)
    return sum(subblock_cost([lengths[i] for i in s]) for s in subblocks)

##########################################################################

@forward
def water_molarity(T) -> float:
    '''Water molarity at a given temperature in Kelvin'''

@forward
def total_ram() -> float:
    '''Total RAM that NUPACK uses'''

@forward
def set_total_ram(n) -> None:
    '''Set total RAM that NUPACK uses'''

@forward
def total_cpu() -> int:
    '''Total CPU cores that NUPACK uses'''

@forward
def set_total_cpu(n) -> None:
    '''Set total CPU cores that NUPACK uses'''

@forward
def default_parameters_path() -> str:
    '''Default place to look for parameters'''

@forward
def set_default_parameters_path(path):
    '''Set default place to look for parameters'''

##########################################################################

def expand_sequence(string):
    words = list()
    cur = ""
    num = ""

    def add_word(x, cur, num):
        if num: num = int(num)
        else: num = 1
        words.append(cur * num)
        return x, ""

    string = string.replace(" ", "") # remove whitespace
    for x in string:
        if str.isdigit(x): num += x
        elif str.isalpha(x): cur, num = add_word(x, cur, num)
        else: raise ValueError("not valid nucleic acid sequence specification")

    add_word("", cur, num)

    return "".join(words)

##########################################################################

