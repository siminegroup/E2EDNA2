import enum, numpy

from .utility import nbits, match, check_instance
from .constants import as_sequences
from .core import PairList, JSON, Pickleable, RawComplex
from .rebind import forward, Dict

################################################################################

@forward
class Conditions:
    '''
    Specification of environmental conditions for the simulated system
    '''
    temperature: float
    na_molarity: float
    mg_molarity: float

    def __init__(self, celsius=37, sodium=1.0, magnesium=0.0, kelvin=None, _fun_=None):
        '''
        `kelvin`: temperature in Kelvin
        `sodium`: molarity of sodium in moles per liter
        `magnesium`: molarity of magnesium in moles per liter
        `celsius`: temperature in Celsius (overrides Kelvin if specified)
        '''
        _fun_(self)
        self.temperature = celsius + 273.15 if kelvin is None else kelvin
        self.na_molarity = sodium
        self.mg_molarity = magnesium

    def save(self):
        return dict(temperature=self.temperature, sodium=self.na_molarity, magnesium=self.mg_molarity)

################################################################################

@forward
class ParameterFile:
    '''
    Descriptor of the chosen parameter file
    '''

    path: str

    def __init__(self, name_or_file='rna'):
        '''
        Initialize from parameter set name or file name
        - 'rna'/'RNA' is translated to 'rna06'
        - 'dna'/'DNA' is translated to 'dna04'
        '''

################################################################################

@forward
class ParameterData:
    '''
    Raw free energy or enthalpy parameters loaded from parameter files
    '''

    def __init__(self, file: ParameterFile, kind: str):
        '''Initialize from file name and kind'''

    def array(self) -> numpy.ndarray:
        pass

################################################################################

@forward
class ParameterInfo:
    '''
    Unique descriptor of a parameter set
    '''
    file: ParameterFile
    kind: str
    loop_bias: float
    temperature: float

    def __init__(self, file, kind, loop_bias=0, temperature=37+273.15):
        '''Initialize from file name and other metadata'''

################################################################################

@forward
class ParameterSet(ParameterData):
    '''
    Free energy parameters loaded from parameter files
    '''
    metadata: ParameterInfo
    default_wobble_pairing: bool
    material: str

    def __init__(self, metadata: ParameterInfo):
        '''Initialize from metadata'''

################################################################################

@forward
class Ensemble(enum.IntEnum):
    '''
    Set of recursions to use in dynamic programs
    '''
    nostacking = 0
    stacking = 1
    some_nupack3 = 2
    all_nupack3 = 3
    none_nupack3 = 4

    @classmethod
    def get(cls, key):
        '''Construct Ensemble from either a string or an integer'''
        return cls(cls.lookup[key] if isinstance(key, str) else key)

    def __str__(self):
        return next(k for k, v in type(self).__members__.items() if v == self)

Ensemble.lookup = {
    'nostacking':   Ensemble.nostacking,
    'stacking':     Ensemble.stacking,
    'none-nupack3': Ensemble.none_nupack3,
    'some-nupack3': Ensemble.some_nupack3,
    'all-nupack3':  Ensemble.all_nupack3,
}

################################################################################

@forward
class Model(Pickleable):
    ensemble: Ensemble
    beta: float
    parameters: ParameterSet
    conditions: Conditions

    def __init__(self, ensemble=Ensemble.stacking, material="rna", wobble=None,
                 kelvin=None, celsius=37, sodium=1.0, magnesium=0.0, bits=64, _fun_=None):
        '''
        Construct a Model

        - `ensemble`: str, one of (nostacking, stacking, none-nupack3, some-nupack3, all-nupack3)
        - `wobble`: enable wobble pairs or not
        - `material`: ParameterFile or equivalent object
        - `kelvin`: temperature in Kelvin
        - `celsius`: temperature in Celsius (may be specified instead of Kelvin).
        - `sodium`: molarity of sodium in moles per liter
        - `magnesium`: molarity of magnesium in moles per liter
        - `bits`: number of bits in the floating type used (probably 32 or 64)
        '''
        if isinstance(material, str):
            material = ParameterFile(material)
        conditions = Conditions(celsius, sodium, magnesium, kelvin=kelvin)

        bits = nbits(bits)
        cls = match(k for k, v in self._metadata_.items() if v.cast(int) == bits)
        gu = 2 if wobble is None else int(bool(wobble))

        _fun_(self, Ensemble.get(ensemble), material, conditions, gu, return_type=cls)

    def save(self):
        out = dict(pseudoknots=False, ensemble=int(self.ensemble))
        out.update(self.conditions.save())
        out.update(self.parameters.save())
        return out

    @property
    def bits(self):
        return self._metadata_[self.type()].cast(int)

    @property
    def material(self):
        '''String describing material: usually RNA or DNA'''
        return self.parameters.material

    @property
    def temperature(self):
        '''Temperature in Kelvin'''
        return self.conditions.temperature

    def boltz(self, energy) -> float:
        '''Return the Boltzmann factor corresponding to a given free energy (kcal/mol)'''

    def dangle5(self, b, c, d) -> float:
        '''Return dangle energy on 5' side'''

    def dangle3(self, b, c, d) -> float:
        '''Return dangle energy on 5' side'''

    def stack_energies(self, loop, structure=None, _fun_=None):
        '''
        Return a dict of stacking states to their stacking energies
        - loop: a list of str or + delimited string
        - structure: a Structure or dot parens plus string
        Structure may also be given as the int index of the sequence following the nick (or -1)
        Structure does not need to be given for a non-exterior loop
        '''
        loop, structure = loop_structure(loop, structure)
        return _fun_(self, loop, structure).cast(Dict[str, float])

    def loop_energy(self, loop, structure=None, _fun_=None):
        '''
        Return a loop free energy
        - loop: a list of str or + delimited string
        - structure: a Structure or dot parens plus string
        Structure may also be given as the int index of the sequence following the nick (or -1)
        Structure does not need to be given for a non-exterior loop
        '''
        loop, structure = loop_structure(loop, structure)
        return _fun_(self, loop, structure).cast(float)

    def structure_energy(self, strands, structure, distinguishable=False):
        '''Free energy of a given secondary structure (kcal/mol)'''
        return structure_energy(as_sequences(strands), PairList(structure),
                                self, distinguishable=distinguishable)

    def structure_probability(self, strands, structure, free_energy):
        '''Probability of a given structure forming given its strands, structure, and complex free energy'''
        return self.boltz(self.structure_energy(strands, structure) - free_energy)

    def __str__(self):
        return 'Model(%r, %r, T=%r K)' % (
            str(self.ensemble), self.parameters.info.file.path, self.temperature)

################################################################################

def loop_structure(loop, structure, _fun_=None):
    '''Return parsed loop and index of nick for a specified loop structure'''
    loop = RawComplex(loop)
    if structure is None:
        return loop, -1
    if isinstance(structure, int):
        return loop, structure
    return loop, _fun_(loop.strands, PairList(structure))

################################################################################

def structure_energy(strands, structure, model, distinguishable=False, _fun_=None):
    '''Free energy of a given secondary structure in (kcal/mol)'''
    check_instance(model, Model)
    return _fun_(strands, structure, model, distinguishable).cast(float)

################################################################################

def loop_energy(sequences, *, structure=None, model):
    check_instance(model, Model)
    return model.loop_energy(sequences, structure=structure)

################################################################################
