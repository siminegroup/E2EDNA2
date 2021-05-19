'''
Main NUPACK module
'''

import os, numpy, scipy.sparse, typing

################################################################################

from . import info
__version__ = info.version

################################################################################

from .model import Model, ParameterFile, Conditions, Ensemble, \
    structure_energy, loop_energy

from .core import PairList, Local, Sequence, Base, Structure, JSON, Domain, \
    RawStrand, RawComplex, TargetComplex, TargetStrand, PairsMatrix, struc_distance, \
    seq_distance

from . import thermo, analysis

from .named import Strand, Complex, Result as AnalysisResult, tube_analysis, \
    Tube, complex_analysis, complex_concentrations, ComplexSet, SetSpec

from .analysis import ConcentrationSolver, ComplexResult, \
    energy, structure_probability, ensemble_size, pfunc, mfe, pairs, subopt, sample

from .design import TargetTube, complex_design, tube_design, \
    Match, Pattern, SSM, Library, Diversity, Similarity, Complementarity, \
    Window, EnergyMatch, Weights, Design, TimeInterval, \
    WriteToFileCheckpoint, StopCondition, des, defect, \
    Options as DesignOptions, Result as DesignResult

from .constants import reverse_complement, random_sequence, as_sequences
from . import concentration

################################################################################

try:
    from .cpp import document
except ImportError:
    raise ImportError('C++ module object "nupack.cpp" cannot be imported')

################################################################################

def _input_csc(A):
    return (numpy.asfortranarray(A.indices, dtype=numpy.ulonglong),
            numpy.asfortranarray(A.indptr, dtype=numpy.ulonglong),
            numpy.asfortranarray(A.data), *A.shape)

def _output_csc(A):
    data, ind, ptr, shape = A.cast(typing.Tuple[
        numpy.ndarray, numpy.ndarray, numpy.ndarray, typing.Tuple[int, int]])
    return scipy.sparse.csc_matrix((data, ind, ptr[:-1]), shape)

def render():
    from .rebind import render_module
    doc, cfg = render_module(__name__, document, monkey=[analysis])
    cfg.set_input_conversion(numpy.ndarray, numpy.asfortranarray)
    cfg.set_input_conversion(scipy.sparse.csc_matrix, _input_csc)

    cfg.set_output_conversion(numpy.ndarray, lambda a: numpy.asarray(a.cast(memoryview)))
    cfg.set_output_conversion(scipy.sparse.csc_matrix, _output_csc)
    cfg.set_output_conversion(Ensemble, lambda x: Ensemble(x.cast(int)))

    # cfg.set_input_conversion(Domain, lambda c: c.sequence)
    # cfg.set_input_conversion(RawStrand, lambda c: c.sequence)
    # cfg.set_input_conversion(RawComplex, lambda c: [s.sequence for s in c.strands])
    return doc, cfg

################################################################################

rendered_document, config = render()
config.parallelism = False
config.cache = 2.0

constants.TypeIndex = rendered_document['TypeIndex']
constants.set_default_parameters_path(os.path.join(os.path.dirname(__file__), 'parameters'))

################################################################################
