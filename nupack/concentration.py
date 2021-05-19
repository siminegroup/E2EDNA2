import numpy as np
from . import constants
from .rebind import forward

################################################################################

@forward
class Output:
    '''
    The output of the equilibrium concentration solving algorithm
    '''
    solution: np.ndarray
    dual_solution: np.ndarray
    objective: float
    error: float
    iters: int
    converged: bool

################################################################################

@forward
class Options:
    '''
    The options for the equilibrium concentration solving algorithm
    The tolerance is specified in terms of the maximum relative error in
    log strand concentration. If the current strand concentrations are a vector s,
    and the initial concentrations are a vector s0, this equates to:

    ```python
    converged = abs(log(s / s0)).max() < tol
    ```
    '''
    max_iters: int
    tolerance: float
    delta_min: float
    delta_max: float
    orthogonalize: bool
    method: int

    def __init__(self, delta_max=1000.0, delta_min=1e-12, max_iters=10000, tolerance=1e-08, method=1, orthogonalize=True, _fun_=None):
        '''
        Initialize solving options from specified keywords:
        - max_iters: maximum number of iterations
        - tolerance: tolerance for convergence
        - delta_min: minimum trust region radius
        - delta_max: maximum trust region radius
        - orthogonalize: whether or not to orthogonalize the input coefficient matrix
        - method: initialization
            0. fit method (not supported)
            1. coordinate descent
            2. initialization from uniform potentials
            3. initialization from given concentrations
            4. initialization from non-negative least squares solution
            5. initialization from absolute value of least squares solution
        '''
        _fun_(self)
        self.max_iters = max_iters
        self.tolerance = tolerance
        self.delta_min = delta_min
        self.delta_max = delta_max
        self.orthogonalize = orthogonalize
        self.set_method(int(method))

################################################################################

@forward
def solve(A, logb, logq, options=None, strict=True, _fun_=None):
    '''
    Solve equilibrium concentrations in any system

    - `A` (ndarray): reaction coefficient matrix (n_strands, n_complexes)
    - `logb` (ndarray): log of initial strand concentrations (n_strands)
    - `logq` (ndarray): logarithms of partition functions (n_complexes)
    - `options`: (Options) solving options
    '''
    A, logb, logq = [np.asarray(o, dtype=np.float64) for o in (A, logb, logq)]
    options = Options() if options is None else options
    out =  _fun_(np.asfortranarray(A), logb, logq, options)
    if strict:
        assert out.converged, 'Solver did not converge'
    return out

################################################################################

@forward
def solve_complexes(indices, logq, x0, options, rotational_correction=True, as_strands=True) -> Output:
    '''
    Solve equilibrium concentrations (non-dimensionalized)
    Users should generally use `complex_concentrations`.
    '''

################################################################################

@forward
def solve_complex_concentrations(indices, logq, x0, kelvin, *,
    rotational_correction=True, as_strands=True, strict=True, **kws):
    '''
    Solve equilibrium concentrations for a set of specified complexes

    - `indices` (List[List[int]]):
        for each complex, indices of each strand
    - `logq` (ndarray):
        the log partition function of each complex
    - `x0` (ndarray):
        initial concentrations for each complex
      `rotational_correction`: True if rotational correction needs to be applied here
    - `kelvin` (float): temperature in Kelvin
    '''
    logq = np.asarray(logq, dtype=np.float64)
    x0 = np.asarray(x0, dtype=np.float64) / constants.water_molarity(kelvin)
    out = solve_complexes(indices, logq, x0, Options(**kws),
        rotational_correction=rotational_correction, as_strands=as_strands)
    if strict:
        assert out.converged, 'Solver did not converge'
    return out.solution * constants.water_molarity(kelvin)

################################################################################

def coefficient_matrix(complexes, strands, dtype=None):
    A = np.zeros((len(complexes), len(strands)), dtype=dtype)
    for i, c in enumerate(complexes):
        for strand in c:
            A[i, strands.index(strand)] += 1
    return A
