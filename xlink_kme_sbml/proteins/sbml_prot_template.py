#!/usr/bin/env python
# %%
import sys

import numpy as np
import pandas as pd
import argparse

from prot_tools import prot_lib
from library import sbml_constants as const, sbml_xl
from Bio import SeqIO

desc = """Kai Kammer - 2021-03. 
Create a kinetic model of a protein. The final model depends on the input files provided. Providing only a fasta
file will create a model with completely random or fixed reactivities. 
Additionally providing files containing accessible surface area (SAS) and/or pKa values for the residues will
incorporate these into the lysine reactivity. 
The same is true for crosslinks: if a file containing SASD or Euclidean distances is given the crosslink 
reactivities will be adjusted accordingly. Alternatively they can either be fixed or be drawn from a normal distribution
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# not sure if just using type=str is better
parser.add_argument('input', action="store", default=None, type=argparse.FileType('r'), nargs='+',
                    help="List of input files separated by spaces. Minimum requirement is a fasta file."
                         "Further possible files: pKa file containing lysine pKa values. "
                         "asa file containing the accessible surface area. "
                         "Distance file containing crosslink distances. "
                         "Offset file containing offsets between pdb and fasta if using pKa or asa as input "
                         "(and there are actually offsets). ")
parser.add_argument('-rm', '--reactivity_monolink', action="store_true", required=True,
                    help="Possible values: r (random), f (fixed), i (input)")
parser.add_argument('-rxl', '--reactivity_crosslink', action="store_true", required=True,
                    help="Possible values: r (random) or any numeric value")
parser.add_argument('-e', '--exp_name', action="store_true", default="my_exp",
                    help="Name of the experiment to be used in the model (i.e. the protein)")
parser.add_argument('-nc', '--num_compartments', action="store_true", default=1,
                    help="Number of compartments")
parser.add_argument('--len_xl_min', action="store_true", default=5,
                    help="Minimum Crosslink length (Å) when using external distances file")
parser.add_argument('--len_xl_max', action="store_true", default=35,
                    help="Maximum Crosslink length (Å) when using external distances file")
parser.add_argument('-o', '--outname', type=argparse.FileType('w'),
                    default='sbml_model.xml', help="Name of the output file")
args = parser.parse_args()


# %%
class AllXLReactions(sbml_xl.AllXLReactionsNoDiff):

    def add_lys(self):
        lys_list = []
        for pos in self.params[const.D_POSITION_LYSINE]:
            user_data_dict = sbml_xl.get_user_data_dict(s_type=const.S_LYS, s_pos=pos, s_prot=const.S_EXP)
            s = self.min_model.sbml_model.addSpecies(user_data_dict[const.D_ID], 1, )
            s.setSpeciesType(user_data_dict[const.D_TYPE])
            s.UserData = user_data_dict
            lys_list.append(s)
        return lys_list

    def add_lys_params(self):
        lys_to_param_dict = {}
        for lys in self.species_lys:
            user_data_dict = sbml_xl.get_user_data_dict(
                s_type=const.S_K_LYS, s_precursor_list=[lys]
            )
            p_klys = self.min_model.sbml_model.addParameter(
                user_data_dict[const.D_ID],
                self.params[const.D_REACTIVITY_DICT_MONO][user_data_dict[const.D_LOCATION_LIST][0][1]]
            )  # type: tesbml.libsbml.Parameter
            p_klys.UserData = user_data_dict
            lys_to_param_dict[user_data_dict[const.D_LOCATION_ID]] = p_klys
        return lys_to_param_dict

    def add_xl_trans_params(self):
        unique_id_kon_dict = {}
        for lys in self.species_lys:
            location_id_lys = lys.UserData[const.D_LOCATION_ID]
            for mono in self.react_mono.products:
                location_id_mono = mono.UserData[const.D_LOCATION_ID]
                if location_id_lys == location_id_mono:
                    continue
                user_data_dict = sbml_xl.get_user_data_dict(
                    s_type=const.S_K_ON_XL, s_precursor_list=[lys, mono]
                )
                location_id = user_data_dict[const.D_LOCATION_ID]
                if location_id in unique_id_kon_dict:
                    continue
                pos_lys = lys.UserData[const.D_LOCATION_LIST][0][1]
                pos_mono = mono.UserData[const.D_LOCATION_LIST][0][1]
                pos_tuple = (sorted([pos_mono, pos_lys]))
                react_dict = self.params[const.D_REACTIVITY_DICT_XL]
                if pos_tuple in react_dict:
                    p_kon_xl = self.min_model.sbml_model.addParameter(
                        user_data_dict[const.D_ID],
                        react_dict[pos_tuple],
                        units=const.S_UNIT_LITRE_PER_MOLE_PER_SECOND,
                    )
                    p_kon_xl.UserData = user_data_dict
                    unique_id_kon_dict[location_id] = p_kon_xl
        return unique_id_kon_dict


class MonoOnly(AllXLReactions):
    """
    Monolink Only Reactions
    """

    def create_xl_trans(self):
        pass

    def create_xl(self):
        pass

    def add_xl_trans_params(self):
        return {}


def get_lys_pos_list(fasta_record):
    return [n + 1 for n, c in enumerate(fasta_record.seq) if c == 'K']


def get_distance(pos1, pos2, exp, df_dist, metric="SASD"):
    s_exp = "exp_name"
    s_pos1 = 'pos1'
    s_pos2 = 'pos2'
    if pos1 != pos2:
        d = df_dist[(df_dist[s_exp] == exp) & (df_dist[s_pos1] == str(pos1)) & (df_dist[s_pos2] == str(pos2))][metric]
        if len(d) > 0:
            return d.iloc[0]
        else:
            d = df_dist[(df_dist[s_exp] == exp) & (df_dist[s_pos1] == str(pos2)) & (df_dist[s_pos2] == str(pos1))][
                metric]
            if len(d) > 0:
                return d.iloc[0]
    return -1


def get_fixed_reactivity():
    return 0.01  # np.random.exponential(scale=0.1)


def get_random_reactivity():
    np.random.exponential(scale=0.1)


def get_klys_dict(lys_pos_list, df_asa):
    # first normalize in [0, 1]
    # then scale to [min_range, max_range]
    min_range = 0.0001
    max_range = 0.1
    val_array = np.array(list(asa_dict.values()))
    m = min(val_array)
    range1 = max(val_array) - m
    if lys_pos_list in asa_dict:
        asa = asa_dict[lys_pos_list]
        asa = (asa - m) / range1
        range2 = max_range - min_range
        return (asa * range2) + min_range
    print(f"Lys at pos {lys_pos_list} has no entry in asa dict. Assigning default value")
    return 0.005


# normal distribution centered around 23 A
def get_kon_xl_dict(lys_pos_list, df_dist, mu=23, sigma=5):
    if dist > 35 or dist < 5:
        return 0
    return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((dist - mu) / sigma) ** 2)


def get_fasta_record(fasta):
    with open(fasta, "r") as handle:
        record = list(SeqIO.parse(handle, "fasta"))[0]
        return record


def is_fasta_file(file_path):
    with open(file_path, "r") as file:
        first_line = file.readline()
        if first_line[0] == '>' and first_line.count("|") == 2:
            return True
    return False


def check_for_cols(file_path, cols):
    with open(file_path, "r") as file:
        first_line = file.readline()
        if all(col in first_line for col in cols):
            return True
    return False


def get_input(inp):
    file_dict = {}
    cols_dist = [prot_lib.COL_DIST_EU, prot_lib.COL_PDB_FILE, prot_lib.COL_DIST_SASD, prot_lib.COL_UXID]
    cols_dssp = [prot_lib.COL_ASA, prot_lib.COL_PDB_CHAIN_ID, prot_lib.COL_PDB_FILE, prot_lib.COL_SEC_STRUC]
    cols_pka = [prot_lib.COL_RES, prot_lib.COL_PKA, prot_lib.COL_POS, prot_lib.COL_PDB_CHAIN_ID]
    # possible files: fasta (required), pKa values from DelPhiPKa web, SASD distances from jwalk
    inp_list = inp.split()
    for inp in inp_list:
        if is_fasta_file(inp):
            file_dict['fasta'] = inp
        elif check_for_cols(inp, cols_dist):
            file_dict['dist'] = inp
        elif check_for_cols(inp, cols_dssp):
            file_dict['dssp'] = inp
        elif check_for_cols(inp, cols_pka):
            file_dict['pka'] = inp
    assert file_dict['fasta'], "ERROR: At least a fasta file is required"
    return file_dict


def get_pos_value_dict(lys_pos_list, df, value_name):
    df_g = df.groupby(prot_lib.COL_PDB_CHAIN_ID).apply(lambda x: dict(zip(x[prot_lib.COL_POS], x[value_name]))).reset_index(name=value_name)
    val_dict = dict(zip(df_g[prot_lib.COL_PDB_CHAIN_ID], df_g[value_name]))
    return val_dict


def main():
    inp_dict = get_input(args.input)
    fasta_record = get_fasta_record(inp_dict["fasta"])
    pos_lys_list = get_lys_pos_list(fasta_record)

    df_dist = None
    df_dssp = None
    df_pka = None

    if 'dist' in inp_dict:
        df_dist = pd.read_csv(inp_dict['dist'])
        split = df_dist[prot_lib.COL_UXID].str.split(":", expand=True)
        df_dist[prot_lib.COL_POS_1] = split[1]
        df_dist[prot_lib.COL_POS_2] = split[4]

    if 'dssp' in inp_dict:
        df_dssp = pd.read_csv(inp_dict['dssp'])

    if 'pka' in inp_dict:
        df_pka = pd.read_csv(inp_dict['pka'])

    if df_dssp and df_pka:
        pass
    elif df_dssp:
        react_mono_dict = get_klys_asa
    elif df_pka:
        pass
    else:
        react_mono_dict = get_random_reactivity()

    if df_dist:
        # TODO: only add to dict if dist in valid range
        # if self.params[const.D_MIN_DIST] < dist < self.params[const.D_MAX_DIST]:
            pass
    else:
        pass

    min_model = sbml_xl.MinimalModel(kinetic_params={'kh': 1e-5, 'koff': 1e-1, 'kon': 1e-3}, c_linker=111)
    params = {
        const.D_POSITION_LYSINE: pos_lys_list,
        const.S_EXP: args.exp_name,
        const.D_MIN_DIST: args.len_xl_min,
        const.D_MAX_DIST: args.len_xl_max,
        const.D_REACTIVITY_DICT_MONO: react_mono_dict,
        const.D_REACTIVITY_DICT_XL: react_xl_dict,
    }
    all_reactions = AllXLReactions(min_model, params)
    if args.num_compartments > 1:
        min_model.add_compartments(args.num_compartments - 1)
    # %%
    with open(args.output, "w") as f:
        f.write(min_model.sbml_model.toSBML())
    print(f"Write model xml to {args.output}")


if __name__ == "__main__":
    main()
