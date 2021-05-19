

class MultitubeObjective:
    """Average of all normalized tube ensemble defects"""
    __new__ = NotImplemented



class TubeObjective:
    """Normalized ensemble defect of a named test tube"""
    __new__ = NotImplemented

    tube_name: str
    tube_id: int


class ComplexObjective:
    """Normalized ensemble defect of a named complex"""
    __new__ = NotImplemented

    complex_name: str
    complex_id: int


class SSMObjective:
    """Degree of sequence symmetry violation for a set of complexes"""
    __new__ = NotImplemented

    word_size: int


class PatternObjective:
    """Fraction of windows of strands or domains """
    __new__ = NotImplemented



class SimilarityObjective:
    """Degree of similarity violation for a set of strands"""
    __new__ = NotImplemented



class EnergyEqualizationObjective:
    """Degree of binding energy inequality for a set of regions"""
    __new__ = NotImplemented



class Objective:
    """A Variant type the holds any of the other classes in design.objectives"""
    __new__ = NotImplemented



