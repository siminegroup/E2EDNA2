from typing import Tuple, List, Dict, Union, Callable
from ..core import JSON, Pickleable
from scipy.sparse import csc_matrix

class Partition:
    """Partition of complexes in design between active and passive sets"""
    # mask: 'std.Vector' # Tuple[bool, ...]
    deflate: float


class Stats:
    """Run statistics of a finished design"""
    num_leaf_evaluations: int
    num_reseeds: int
    # num_redecompositions: 'std.Vector' # Tuple[int, ...]
    # offtargets_added_per_refocus: 'std.Vector' # Tuple[int, ...]
    design_time: float
    analysis_time: float
    # final_Psi: 'Partition'


class Complex:
    """Result of a complex calculation"""
    name: str
    # sequence: 'core.Complex'
    # structure: 'components.Structure'
    log_pfunc: float
    pair_probabilities: csc_matrix
    defect: float
    normalized_defect: float


class TubeComplex:
    """Result of a complex in a tube"""
    name: str
    concentration: float
    target_concentration: float
    nucleotide_defect: float
    structural_defect: float
    concentration_defect: float
    defect: float


class Tube:
    """Result of a tube calculation"""
    name: str
    nucleotide_concentration: float
    nucleotide_defect: float
    normalized_defect: float
    defect: float
    # complexes: 'std.Vector' # Tuple['TubeComplex', ...]


class Single:
    """One of the pareto optimal results at end of design"""
    domains: Dict[str, str]
    strands: Dict[str, str]
    defects: List[float]
    weighted_defects: List[float]


class Result(Pickleable):
    """Result of a design"""
    # model: 'components.ModelSettings'
    # parameters: 'components.Parameters'
    # stats: 'Stats'
    # complexes: 'std.Vector' # Tuple['Complex', ...]
    # tubes: 'std.Vector' # 'std.Vector' # Tuple['Tube', ...]
    success: bool

    def __init__(self, json=None, _fun_=None):
        """
        create a new Result, loading from JSON if specified

        Arguments:
            json: a valid nupack.JSON object specifying a Result object.
        """
        if json:
            self.move_from(self.from_json(json))
        else:
            _fun_(self)

