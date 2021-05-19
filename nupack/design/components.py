from typing import Tuple, List, Dict, Union, Callable
import datetime, pathlib

# from .. import model # import Conditions, ParameterFile
# from .. import core
from .. import constants, utility
from ..rebind import forward

@forward
class Timer:
    """A timer for design"""
    def start(self, _fun_=None):
        """start the timer"""
        _fun_(self)

    def elapsed(self) -> float:
        """time (seconds) since timer was started"""

    def stop(self) -> float:
        """stop timer if still running and return time between start and initial stop"""

@forward
class ModelSettings:
    """The set of parameters setting the thermodynamic conditions of the
    design
    """
    # parameters: 'model.ParameterFile'
    dangle_type: str
    # conditions: 'model.Conditions'

    def __init__(self, material='RNA', ensemble='stacking', temperature=None, sodium=1.0, magnesium=0.0, _fun_=None):
        """create model settings"""
        if temperature == None:
            temperature = constants.DefaultTemperature
        _fun_(self, material, ensemble, temperature, sodium, magnesium)

@forward
class Domain:
    """A string of nucleotides"""
    name: str
    allowed_bases: str

    def __init__(self, name, sequence):
        pass

@forward
class Strand:
    """A concatenated string of domains with no strand breaks"""
    name: str
    # domain_names: 'std.Vector' # Tuple[str, ...]

    def __init__(self, name, domain_names):
        pass

@forward
class Complex:
    """An ordered list of strands, possibly with a target structure"""
    name: str
    # strands: 'std.Vector' # Tuple[str, ...]
    # structure: 'Structure'

@forward
class Tube:
    """A specification for a test tube with complexes at a given target concentration"""
    name: str
    # targets: 'std.Vector' # Tuple[Tuple[Tuple[str, ...], float], ...]


# Constraints
@forward
class DualList:
    """Holds two lists of domains that are to be constrained"""
    # left: 'std.Vector' # Tuple[str, ...]
    # right: 'std.Vector' # Tuple[str, ...]


@forward
class Pattern:
    """Pattern constraint specification"""
    name: str
    pattern: str

@forward
class Diversity:
    """Diversity constraint specification"""
    name: str
    word_length: int
    min_nucleotide_types: int

@forward
class Word:
    """Word constraint, i.e. a library or window constraint, specification"""
    # domains: 'std.Vector' # Tuple[str, ...]
    # comparisons: 'std.Vector' # Tuple[Tuple[str, ...], ...]

@forward
class Similarity:
    """Similarity constraint specification"""
    name: str
    reference: str
    range: Tuple[float, float]

@forward
class Constraints:
    """Collection of all constraint specifications"""
    # complementarity: 'std.Vector' # Tuple['DualList', ...]
    # match: 'std.Vector' # Tuple['DualList', ...]
    # pattern: 'std.Vector' # Tuple['Pattern', ...]
    # word: 'std.Vector' # Tuple['Word', ...]
    # similarity: 'std.Vector' # Tuple['Similarity', ...]

@forward
class Parameters:
    """Parameters"""
    f_stop: float
    f_sparse: float
    f_passive: float
    H_split: int
    N_split: int
    f_split: float
    f_stringent: float
    dG_clamp: float
    M_bad: int
    M_reseed: int
    M_reopt: int
    f_redecomp: float
    f_refocus: float
    rng_seed: int
    log: str = ''
    decomposition_log: str = ''
    thermo_log: str = ''

class StopCondition:
    def __init__(self):
        self.stop = False

    def cancel(self):
        self.stop = True

    def __call__(self, stats, timer, done) -> int:
        if done:
            return +1
        elif self.stop:
            return -1
        else:
            return 0


class TimeInterval(StopCondition):
    def __init__(self, interval: float):
        super().__init__()
        self.interval = interval

    def __call__(self, stats, timer, done) -> int:
        stop = super().__call__(stats, timer, done)
        if stop: return stop
        current = timer.elapsed()
        if timer.elapsed() > self.interval:
            self.last = current
            return +1
        return 0


class MutationInterval(StopCondition):
    def __init__(self, interval: int):
        super().__init__()
        self.interval = interval
        self.last = 0

    def __call__(self, stats, timer, done) -> int:
        stop = super().__call__(stats, timer, done)
        if stop: return stop
        current = stats.num_leaf_evaluations
        if current - self.last > self.interval:
            self.last = current
            return +1
        return 0


class WriteToFileCheckpoint:
    def __init__(self, path, timespec='seconds'):
        self.path = pathlib.Path(path)
        self.timespec = timespec

    def __call__(self, result):
        time = datetime.datetime.now().isoformat(timespec=self.timespec)
        self.path.mkdir(parents=True, exist_ok=True)
        path = self.path/'{}.json'.format(time)
        path.write_text(result.to_json().dump(indent=4))

