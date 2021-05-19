import logging

from copy import deepcopy as copy
import pandas as pd

logger = logging.getLogger(__name__)
logger.setLevel("WARN")

################################################################################

def pick_first_matching(dictionary, keys=None, default=None):
    if keys is None: keys = []

    for k in keys:
        if k in dictionary:
            return dictionary[k]
    return default

################################################################################

def parse_result(data, sources):
    '''
    Fills out the sources object and returns the sorted tube names
    '''
    pfm = pick_first_matching
    # hack for now
    if "results" in data:
        data = data["results"][0]

    on_targets = sorted([copy(comp) for comp in data["complexes"] if comp["structure"]],
            key=lambda x: pfm(x, ["normalized defect", "normalized_defect"]), reverse=True)

    complexes_dict = {
        "name": [x["name"] for x in on_targets],
        "defect": [pfm(x, ["defect (nt)", "defect"]) for x in on_targets],
        "normalized_defect": [pfm(x, ["normalized defect", "normalized_defect"]) for x in on_targets],
    }

    summary_dict = {"tube_name": [], "structural_defect": [], "concentration_defect": [], "tube_defect": [], "tube_sort": []}
    on_target_dict = {"tube_name": [], "complex_name": [], "structural_defect": [], "concentration_defect": [], "complex_defect": [], "sort": []}
    concentration_dict = {"tube_name": [], "complex_name": [], "target_concentration": [], "concentration": [], "type": []}

    tubes = copy(data["tubes"])

    comp_sort = gate_number
    tube_sort = gate_tuple

    for tube in tubes:
        name = tube["name"]
        conc = pfm(tube, ["[nt] (M)", "nucleotide_concentration"])
        ontargets = [comp for comp in tube["complexes"] if pfm(comp, ["target concentration (M)", "target_concentration"]) > 0]
        complexes = sorted(copy(tube["complexes"]), key=lambda x: pfm(x, ["target concentration (M)", "target_concentration"]), reverse=True)
        # print([pfm(x, ["concentration (M)", "concentration"]) for x in complexes])

        cur_complexes = {
            "tube_name": [name]*len(complexes),
            "complex_name": [c["name"] for c in complexes],
            "target_concentration": [pfm(c, ["target concentration (M)", "target_concentration"]) for c in complexes],
            "concentration": [pfm(c, ["concentration (M)", "concentration"]) for c in complexes]
        }
        cur_complexes["type"] = ["on-target" if x > 0 else "off-target" for x in cur_complexes["target_concentration"]]

        cur_ontargets = {
            "tube_name": [name]*len(ontargets),
            "sort": [comp_sort(c['name']) for c in ontargets],
            "complex_name": [c["name"] for c in ontargets],
            "structural_defect": [pfm(c, ["structural defect (M)", "structural_defect"])/conc for c in ontargets],
            "concentration_defect": [pfm(c, ["concentration defect (M)", "concentration_defect"])/conc for c in ontargets],
            "complex_defect": [pfm(c, ["defect contribution", "normalized_defect_contribution"]) for c in ontargets]
        }

        total_complex = {
            "tube_name": name,
            "sort": comp_sort("total"),
            "tube_sort": tube_sort(name),
            "complex_name": "total",
            "structural_defect": sum(cur_ontargets["structural_defect"]),
            "concentration_defect": sum(cur_ontargets["concentration_defect"]),
            "complex_defect": pfm(tube, ["normalized defect", "normalized_defect"]),
            "tube_defect": pfm(tube, ["normalized defect", "normalized_defect"]),
        }

        for k, v in cur_complexes.items():
            concentration_dict[k] += v

        for k, v in cur_ontargets.items():
            v.append(total_complex[k])
            on_target_dict[k] += v

        for k, v in summary_dict.items():
            v.append(total_complex[k])


    # print({k: len(v) for k, v in summary_dict.items()})
    # print({k: len(v) for k, v in complexes_dict.items()})
    # print({k: len(v) for k, v in concentration_dict.items()})
    # print({k: len(v) for k, v in on_target_dict.items()})

    complexes_df = pd.DataFrame.from_dict(complexes_dict)
    summary_df = pd.DataFrame.from_dict(summary_dict).sort_values(by=["tube_sort", "tube_name"])
    conc_df = pd.DataFrame.from_dict(concentration_dict).sort_values(by=["concentration"], ascending=False)
    on_target_df = pd.DataFrame.from_dict(on_target_dict).sort_values(by=["tube_name", "sort", "complex_name"], axis=0)

    tube_names = sorted(set(on_target_df["tube_name"]), key=tube_sort)

    sources.summary = summary_df
    sources.complexes = complexes_df
    sources._per_tube_defects = on_target_df
    sources._concentrations = conc_df

    return tube_names

################################################################################

def complex_defect_contribution(ens_def, conc, target_conc, num_nucs):
    """
    computes the structural and concentration defect contributions of an on-target
    complex to the test tube ensemble defect.

    Arguments:
    ens_def -- the complex ensemble defect
    conc -- the actual concentration in the designed test tube
    target_conc -- the desired concentration in the test tube specification
    num_nucs -- the number of nucleotides in the complex; the sum of the
        lengths of the strands in the complex

    Return value:
    the structural component of the defect, the concentration component of the
            defect, and the total defect (the sum of the latter two)

    """

    str_def = ens_def * min(conc, target_conc)
    conc_def = num_nucs * max(target_conc - conc, 0)

    return str_def, conc_def, str_def + conc_def

################################################################################

def load_single_residuals(filename):
    f = open(filename)

    in_tubes = False
    in_structures = False

    temp = ""
    tube_name = ""
    tube_names = []
    complex_names = []

    ensemble_defects = {}
    num_nucs = {}
    structural_defects = []
    concentration_defects = []
    complex_defects = []
    struct_defect = None
    conc_defect = None
    combined_defect = None

    tot_struct_defect = 0
    tot_conc_defect = 0
    tot_combined_defect = 0
    num_comp = 0

    tube_conc = None
    target_conc = None
    conc = None

    def add(str_def, conc_def, comp_def, t_name, c_name):
        structural_defects.append(str_def)
        concentration_defects.append(conc_def)
        complex_defects.append(comp_def)
        tube_names.append(t_name)
        complex_names.append(c_name)

    for line in f:
        spl_line = [s.strip() for s in line.split(':')]
        if len(spl_line) > 0:
            key = spl_line[0]
            val = ""
            if len(spl_line) > 1:
                val = spl_line[-1]

            if key == 'tubes':
                in_tubes = True
                in_structures = False
            elif key == 'mean objective':
                # if num_comp > 1:
                add(tot_struct_defect, tot_conc_defect, tot_combined_defect, tube_name, "total")
                in_tubes = False
            elif key == 'structures':
                in_structures = True
            elif in_structures:
                if key == '- name':
                    temp = val
                elif key == 'defect[nt]':
                    ensemble_defects[temp] = float(val)
                elif key == "sequence":
                    num_nucs[temp] = len("".join(val.split('+')))
            elif in_tubes:
                if key == '- name':
                    temp = val
                elif key == 'complexes':
                    if tot_combined_defect > 0: # and num_comp > 1:
                        add(tot_struct_defect, tot_conc_defect, tot_combined_defect, tube_name, "total")

                    tot_struct_defect = 0
                    tot_conc_defect = 0
                    tot_combined_defect = 0
                    num_comp = 0
                    tube_name = temp
                elif key == 'nucleotide conc[M nt]':
                    tube_conc = float(val)
                elif key == "concentration":
                    conc = float(val)
                elif key == "target concentration":
                    target_conc = float(val)
                    if target_conc > 0.0:
                        num_comp += 1
                        struct_defect, conc_defect, combined_defect \
                                = complex_defect_contribution(ensemble_defects[temp],
                                conc,
                                target_conc,
                                num_nucs[temp])

                        tot_struct_defect += struct_defect / tube_conc
                        tot_conc_defect += conc_defect / tube_conc
                        tot_combined_defect += combined_defect / tube_conc

                        add(struct_defect / tube_conc, conc_defect / tube_conc,
                                combined_defect / tube_conc, tube_name, temp)

    return {"structural_defect":structural_defects,
            "concentration_defect":concentration_defects,
            "complex_defect":complex_defects,
            "tube_name":tube_names,
            "complex_name":complex_names
            }

################################################################################

def has_gatename(name):
    try:
        toks = name.split("__")
        _ = int(toks[0][1:])
        if toks[0][0] == "g":
            return True
        return False
    except:
        return False

def gate_number(name):
    if not has_gatename(name):
        return 10e9
    toks = name.split("__")
    return int(toks[0][1:])


def successive_number(tok):
    num = 10e9
    for i in range(1, len(tok)+1):
        temp = 0
        try:
            temp = int(tok[0:i])
        except:
            break

        num = temp

    return num


def gate_tuple(name):
    if not has_gatename(name):
        return 10e9, 10e9
    toks = name.split("__")
    return int(toks[0][1:]), successive_number(toks[-1])


def format_tube_name(name):
    pieces = name.split("__")
    num = float("inf")
    rest = name
    if (len(pieces)) > 1:
        try:
            num = int(pieces[0][1:])
            rest = " ".join(pieces[1:]).replace("_", " ")
        except ValueError:
            rest = " ".join(pieces).replace("_", " ")
    else:
        rest = name.replace("_", " ")


    if (num == float("inf")):
        formatted_name = rest
    else:
        formatted_name = "System {}: {}".format(num, rest)

    return formatted_name, num, rest


def tuple_name(name):
    parts = name.split('__')
    try:
        gate = int(parts[0][1:])
        rest = "∙".join(parts[1:])
        return (gate, rest)
    except ValueError:
        logger.info("bad complex name: {}".format(name))
        return (float("inf"), "∙".join(parts))
