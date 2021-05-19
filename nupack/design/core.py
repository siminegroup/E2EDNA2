from typing import Tuple, List
import json

# nupack includes
from .. import rotation  # import rotation.lowest_rotation
from .. import constants  # import compute_necklaces
from .. import core  # import Local, PairList
from ..model import Model
from ..utility import check_instance

# nupack.design includes
from . import results as rs
from . import components as cmpts
from .objectives import *
from .weights import *

def index(container, condition):
    """ returns first index in container that matches condition, otherwise -1 """
    for i, x in enumerate(container):
        if condition(x):
            return i
    return -1

def replace_or_append(container, element, condition):
    """ check container for element matching condition and overwrite if present,
        else add to container
    """
    ind = index(container, condition)
    if ind > -1:
        container[ind] = element
    else:
        container.append(element)

class Specification:
    """A representation of a multistate test tube design problem containing
    sequences, complexes, tubes, and constraints
    """
    # domains: 'std.Vector' # Tuple['cmpts.Domain', ...]
    # strands: 'std.Vector' # Tuple['cmpts.Strand', ...]
    # complexes: 'std.Vector' # Tuple["cmpts.Complex", ...]
    # tubes: 'std.Vector' # Tuple["cmpts.Tube", ...]
    # model: 'cmpts.ModelSettings'
    # constraints: 'cmpts.Constraints'
    # parameters: 'cmpts.Parameters'
    log: str

    def __init__(self, model, wobble_mutations=False, json=None, _fun_=None):
        _fun_(self, check_instance(model, Model), wobble_mutations)
        if json:
            self.from_json(json)
            return
        self.sources = dict()
        self.libraries = dict()

    def __str__(self, _fun_=None):
        """string of underlying design and sources and libraries"""
        return "\n".join([
            _fun_(self).cast(str),
            "    sources: " + str(self.sources),
            "    libraries: " + str(self.libraries)])

    def to_json(self, _fun_=None, **kwargs):
        """add Python components to C++ JSON serialization and return JSON string"""
        base = json.loads(_fun_(self).dump())
        base["sources"] = self.sources
        base["libraries"] = self.libraries
        return json.dumps(base, **kwargs)

    def from_json(self, json_string, _fun_=None):
        """load Specification from JSON string"""
        temp = json.loads(json_string)
        libraries = temp.pop("libraries", {})
        sources = temp.pop("sources", {})

        self.move_from(_fun_(core.JSON(json.dumps(temp))))
        self.sources = sources
        self.libraries = libraries

    def add_domain(self, name, sequence):
        """add a domain defined by a sequence of degenerate base codes"""
        sequence = sequence.upper() # case-insensitive
        domain = cmpts.Domain(name, sequence)
        replace_or_append(self.domains, domain, lambda x: x.name == name)

        # add generic 'N'*len(sequence) complement if not already defined
        if name[-1] != "*":
            comp_name = name + "*"
            if not any(x.name == comp_name for x in self.domains):
                self.domains.append(cmpts.Domain(comp_name, "N" * len(sequence)))

    def add_strand(self, name, domains):
        """add a strand defined by list of domains"""
        for d in domains:
            if not any(d == x.name for x in self.domains):
                raise ValueError(f"\"{d}\" is not a previously defined domain")

        strand = cmpts.Strand(name, domains)
        replace_or_append(self.strands, strand, lambda x: x.name == name)

    def add_complex(self, name, strands, structure=None, bonus=0):
        """add a complex with a name, its defining strands and optional target structure

        Args:
            name (str): a name for the complex, not the empty string
            strands (list(str) OR str): a single strand name or list of strand names
            structure (None, optional): a dot-parens-plus(-rle) structure string
        """
        comp = cmpts.Complex()
        comp.name = name
        if isinstance(strands, str):
            strands = [strands]
        comp.strands = strands
        if (structure):
            comp.structure = core.Structure(structure)
        comp.bonus = float(bonus)
        self.complexes.append(comp)

    def add_tube(self, name, on_targets):
        """add a tube with its on-targets. If a tube with the same name already
        exists, it will be overwritten.

        Args:
            name (str): tube name
            on_targets (dict(str, float)): map from complex name to target concentration
        """
        tube = cmpts.Tube()
        tube.name = name
        tube.targets = tuple(((k,), v) for k, v in on_targets.items())
        replace_or_append(self.tubes, tube, lambda t: t.name == name)

    def add_off_targets(self, name, *, max_size=0, explicit=None, exclude=None):
        """Add off-targets to existing tube

        Args:
            name (str): The name of the tube to modify. method fails if tube
                is not yet defined.
            max_size (int, optional): For programmatic generation. Causes all
                complexes up to max_size constructed from on-targets, excepting
                the on-targets and those mentioned in exclude
            explicit (None, optional): Specifications of off-targets to be
                added, either names of complexes or lists of names of their
                defining strands for complexes not used elsewhere as on-
                targets
            exclude (None, optional): Specifications of complexes to not
                add. Used to prevent programmatic generation from adding in a
                complex which is expected to conflict with the on-targets.
        """
        if not explicit:
            explicit = []
        if not exclude:
            exclude = []

        def strand_names(names):
            """if names refers to a named complex, return its strand names"""
            if isinstance(names, str):
                names = [names]

            if len(names) == 1:
                comp_matches = [i for i, x in enumerate(self.complexes) if x.name == names[0]]
                if comp_matches:
                    i, = comp_matches  # will error if multiple complexes with same name
                    names = self.complexes[i].strands.cast(List[str])

            for name in names:
                if not any(i.name == name for i in self.strands):
                    raise ValueError("{} is not an strand in any on-target".format(name))
            return tuple(names)

        explicit = [rotation.lowest_rotation(strand_names(x)) for x in explicit]
        exclude = [rotation.lowest_rotation(strand_names(x)) for x in exclude]

        tube_matches = [i for i, x in enumerate(self.tubes) if x.name == name]
        tube_ind, = tube_matches  # will error if multiple tubes with same name
        tube = self.tubes[tube_ind]

        strands = set()
        on_targets = set()

        for comp in tube.targets:
            names = strand_names(comp[0].cast(Tuple[str, ...]))
            on_targets.add(rotation.lowest_rotation(names))
            strands.update(names)

        strands = list(strands)
        off_targets = set()

        def attempt_add_offtarget(names):
            if names not in on_targets and names not in exclude:
                off_targets.add(names)

        def add_perm(l):
            """add implied complex if not an on-target or excluded"""
            names = tuple(rotation.lowest_rotation([strands[i] for i in l]))
            attempt_add_offtarget(names)

        for i in range(1, max_size+1):
            constants.compute_necklaces(add_perm, n_elements=len(strands), size=i)

        for names in explicit:
            attempt_add_offtarget(names)

        off_targets = sorted(off_targets) # repeated calls of same script generate identical output

        low_rots = [rotation.lowest_rotation(c.strands.cast(Tuple[str, ...])) for c in self.complexes]
        for o in off_targets:
            tube.targets.append((o, 0))
            if not any(o == l for l in low_rots):
                self.add_complex("", o)

        # self.tubes = self.tubes[0:tube_ind] + (tube,) + self.tubes[tube_ind+1:]


    def add_library(self, name, sequences):
        """add a named library of uniform length domains

        Args:
            name (str): the libraries name
            sequences (list(str)): the sequences (equal-length) defining the
                library

        Raises:
            ValueError: raised if all sequences are not the same length
        """
        length = len(sequences[0])
        if not all(len(s) == length for s in sequences):
            raise ValueError("all library sequences must be same length")

        self.libraries[name] = sequences

    def domain_lengths(self, domain_names):
        """sum of lengths of a list of domains with names in domain_names"""
        return sum([len(i.allowed_bases) for d in domain_names for i in self.domains if i.name == d])

    def add_library_constraint(self, domains, libraries):
        """constrain the concatenation of domains to be drawn from some
            alternative concatenation of library sequences

        Args:

            domains (list(str)): the names of the domains to be concatenated
                in the constraint

            libraries (list(str)): the list of libraries whose sequences will be applied to the
        """
        if isinstance(domains, str):
            domains = [domains]
        if isinstance(libraries, str):
            libraries = [libraries]

        # validate
        for d in domains:
            if not any(x.name == d for x in self.domains):
                raise ValueError("domain {} not defined".format(d))

        lib_seqs = [self.libraries[x] for x in libraries]
        library_len = sum(len(d[0]) for d in lib_seqs)
        domain_len = self.domain_lengths(domains)
        if domain_len != library_len:
            raise ValueError("concatenated domains are {} nt while concatenated libraries are {} nt".
                             format(domain_len, library_len))

        constraint = cmpts.Word()
        constraint.domains = domains
        constraint.comparisons = lib_seqs
        self.constraints.word.append(constraint)

    def add_source(self, name, sequence):
        """add a source sequence for window constraints to the design"""
        self.sources[name] = sequence

    def add_window_constraint(self, domains, sources):
        """constraint the concatenation of domains to be drawn as a window
        from one of the sources

        Args:
            domains (list(str)): names of concatenated domains to constrain
            sources (list(str) OR str): names of the source sequences to draw
                the windows from

        Raises:
            ValueError: if one of the sources listed is shorter than the
                concatenation of domains and hence has no windows.
        """
        if isinstance(sources, str):
            sources = [sources]

        window_length = self.domain_lengths(domains)
        # validate
        for s in sources:
            if len(self.sources[s]) < window_length:
                raise ValueError("{} is too short of a source for ({})".format(s, ",".join(domains)))

        source_windows = list()
        for s in sources:
            source = self.sources[s]
            windows = sorted(list(set(source[i:i + window_length] for i in range(len(source) - window_length + 1))))
            source_windows += windows

        constraint = cmpts.Word()
        constraint.domains = domains
        constraint.comparisons = [source_windows]
        self.constraints.word.append(constraint)

    def _add_dual_list_constraint(self, left, right, container):
        """internal method handling the common aspects of defining either a
        match or complementarity constraint

        Args:
            left (list(str)): one set of domains to concatenate
            right (list(str)): second set of domains to concatenate
            container (C++ vector of specs): a collection of specs for one or
                the other constraint type

        Raises:
            ValueError: if concatenation of left domains is a different length
                than the concatenation of the right domains, the constraint is
                ill-defined
        """
        # validate
        len_left = self.domain_lengths(left)
        len_right = self.domain_lengths(right)
        if len_left != len_right:
            raise ValueError("left domains are {} nt and right domains are {} nt"
                             .format(len_left, len_right))

        constraint = cmpts.DualList()
        constraint.left = left
        constraint.right = right
        container.append(constraint)

    def add_match_constraint(self, left, right):
        """add constraint making left and right lists of domains identical"""
        self._add_dual_list_constraint(left, right, self.constraints.match)

    def add_complementarity_constraint(self, left, right):
        """add constraint making left and right lists of domains complementary"""
        self._add_dual_list_constraint(left, right, self.constraints.complementarity)

    def add_similarity_constraint(self, domains, reference, lims):
        """make a design element match some fraction of a reference sequence

        Args:
            name (str): domain or strand name
            reference (str): a degenerate base sequence
            lims (tuple(float, float)): min and max fraction that must be matched
        """
        # validate
        # if any(d.name == name for d in self.domains):
        #     domains = [name]
        # elif any(d.name == name for d in self.strands):
        #     s, = [d for d in self.strands if d.name == name]
        #     domains = list(s.domain_names)
        # else:
        #     raise ValueError("{} has not been defined as a domain or strand".format(name))
        for s in domains:
            assert isinstance(s, str)
            assert any(d.name == s for d in self.domains)

        if self.domain_lengths(domains) != len(reference):
            raise ValueError("design element and reference sequences must have same length", self.domain_lengths(domains), len(reference))

        constraint = cmpts.Similarity()
        constraint.domains = domains
        constraint.reference = reference
        constraint.range = (min(lims), max(lims))
        self.constraints.similarity.append(constraint)

    def add_pattern_constraints(self, patterns, names=None):
        """prevent patterns either globally or on a specific set of strands
        and/or domains

        Args:
            patterns (list(str)): the patterns to prevent
            names (list(str), optional): A list of domains and strands to
                prevent the patterns on. Patterns are applied to each strand
                in the design if names is empty or defaulted
        """
        # validate
        if isinstance(patterns, str):
            patterns = [patterns]

        if names:
            if isinstance(names, str):
                names = [names]

            for n in names:
                if not (any(d.name == n for d in self.domains) ^
                        any(s.name == n for s in self.strands)):
                    raise ValueError("{} is not a domain or strand".format(n))

        # convert global
        if not names:
            names = []

        for p in patterns:
            constraint = cmpts.Pattern()
            constraint.domains = names
            constraint.pattern = p
            self.constraints.pattern.append(constraint)

    def add_diversity_constraints(self, word_length, minimum_nucleotide_types, names=None):
        """set the minimum diversity (minimum_nucleotide_types) per
        windows of word_length within each domain/strand in names
        """
        if names:
            for n in names:
                if not (any(d.name == n for d in self.domains) ^
                        any(s.name == n for s in self.strands)):
                    raise ValueError("{} is not a domain or strand".format(n))
        else:
            names = []

        constraint = cmpts.Diversity()
        constraint.domains = names
        constraint.word_length = word_length
        constraint.min_nucleotide_types = minimum_nucleotide_types
        self.constraints.diversity.append(constraint)

    # OBJECTIVES

    def add_global_objective(self, weight=None):
        '''Adds default multitube objective to design'''
        if weight is None:
            weight = 1.0
        x = MultitubeObjective()
        self.objectives.append(x)
        self.weights.add_objective_weight(weight)

    def add_tube_objective(self, name: str, weight=None):
        '''Adds objective for particular named tube to design'''
        if weight is None:
            weight = 1.0
        x = TubeObjective(name)
        self.objectives.append(x)
        self.weights.add_objective_weight(weight)

    def add_complex_objective(self, name: str, weight=None):
        '''Adds objective for particular named complex to design'''
        if weight is None:
            weight = 1.0
        x = ComplexObjective(name)
        self.objectives.append(x)
        self.weights.add_objective_weight(weight)

    def add_pattern_objective(self, patterns=None, names=None, weight=None):
        if isinstance(names, str):
            names = [names]
        if isinstance(patterns, str):
            patterns = [patterns]
        if patterns is None or not patterns:
            raise ValueError("must have at least one pattern in objective")
        if names is None:
            names = []
        if weight is None:
            weight = 1.0
        x = PatternObjective(names, patterns)
        self.objectives.append(x)
        self.weights.add_objective_weight(weight)

    def add_similarity_objective(self, names, reference_sequences, limits, weight=None):
        '''Adds objective version of similarity constraint'''
        if isinstance(names, str):
            names = [names]
        if isinstance(reference_sequences, str):
            reference_sequences = [reference_sequences]
        lo, hi = map(float, limits)
        if weight is None:
            weight = 1.0
        x = SimilarityObjective(names, reference_sequences, [[lo, hi]])
        self.objectives.append(x)
        self.weights.add_objective_weight(weight)

    def add_energy_equalization_objective(self, names, energy=None, weight=None):
        if isinstance(names, str):
            names = [names]

        if weight is None:
            weight = 1.0
        x = EnergyEqualizationObjective(names, energy)
        self.objectives.append(x)
        self.weights.add_objective_weight(weight)

    def add_SSM_objective(self, names, word_size, weight=None):
        if isinstance(names, str):
            names = [names]

        if weight is None:
            weight = 1.0
        x = SSMObjective(names, word_size)
        self.objectives.append(x)
        self.weights.add_objective_weight(weight)

    # WEIGHTS

    def add_weight(self, scale, tube=None, complex=None, strand=None, domain=None):
        weight = Weight(tube, complex, strand, domain, scale)
        self.weights.add(weight)

    def __call__(self, env=None, restart=None, checkpoint_condition=None, checkpoint_handler=None, _fun_=None):
        """create the designer from the fully specified design and run and
        return results

        Arguments:
            env (Local): a C++ compute environment
            restart: a design.results.Result object that can be used to initialize the
                sequence state of the design.
            checkpoint_condition: A 2-argument callable that takes the current design's
                statistics and timer and outputs True if a checkpoint should be emitted
                and False otherwise, e.g. if time has elapsed a certain amount
            checkpoint_handler: A single argument callable that receives a
                design.results.Result object and decides how to process it as a checkpoint,
                e.g. printing to a file

        Returns: Result: the result of a finished design
        """
        if len(self.objectives) == 0:
            self.add_global_objective()
        if not env:
            env = core.Local()
        return _fun_(self, env, checkpoint_condition, checkpoint_handler, restart, gil=False)


    def evaluate(self, env=None, _fun_=None):
        """compute the multistate test tube ensemble defect for a set of fixed
        sequences or throw and exception if the input sequences have
        variable nucleotides.

        Args:
            env (Local): a C++ compute environment

        Returns:
            Result: the result of a finished design
        """
        if len(self.objectives) == 0:
            self.add_global_objective()
        if not env:
            env = core.Local()
        return _fun_(self, env, gil=False).cast(rs.Result)


def complex_design(structure):
    """create all the necessary multitube design interface components to
    do a single complex design.
    """
    design = Specification()
    struc = core.Structure(structure)
    breaks = [0] + [+x for x in struc.nicks]
    lengths = [breaks[i+1] - breaks[i] for i in range(len(struc.nicks))]

    def domain(i):
        return f"domain {i}"

    def strand(i):
        return f"strand {i}"

    for i, l in enumerate(lengths):
        design.add_domain(domain(i), "N" * l)
        design.add_strand(strand(i), [domain(i)])

    design.add_complex("complex", [strand(i) for i in range(len(lengths))], structure)
    design.add_tube("tube", {"complex": 1})

    design.add_global_objective()

    return design
