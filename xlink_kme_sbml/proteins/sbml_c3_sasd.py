#!/usr/bin/env python
# %%
import sys

sys.path.append('/home/kai/Coding/xlink_stats/src/')
import sequence_objects
import read_dssp_columns
import sequence_tools
import numpy as np
import pandas as pd

from xlink_kme_sbml.library import sbml_constants as const, sbml_xl
from Bio import SeqIO

# %%
s_eu_dist = "eucDist"
s_uxid = "uxID"
s_sasd = "SASD"
s_energy = "energy"
s_exp = "exp_name"
s_pos1 = 'pos1'
s_pos2 = 'pos2'

dir_input = '../input/'
df_dist = pd.read_csv(
    f"{dir_input}jwalk_dist_combined.csv"
)
fasta_c3 = f"{dir_input}c3.fasta"

split = df_dist[s_uxid].str.split(":", expand=True)
df_dist[s_pos1] = split[1]
df_dist[s_pos2] = split[4]


# %%
def get_lys_pos(exp='c3'):
    with open(fasta_c3, "r") as handle:
        c3_record = list(SeqIO.parse(handle, "fasta"))[0]
        if exp == 'c3':
            return [n + 1 for n, c in enumerate(c3_record.seq) if c == 'K']
        elif exp == 'c3b':
            # remove the ANA domain: 693-728
            return [n + 1 for n, c in enumerate(c3_record.seq) if c == 'K' and (n < 693 or n > 728)]


def get_asa_dict(lys_pos_list, exp):
    with open(fasta_c3, "r") as handle:
        c3_record = list(SeqIO.parse(handle, "fasta"))[0]
        uni_obj = sequence_objects.UniProtData(c3_record.id, str(c3_record.seq), {c3_record.id: c3_record.description})
        uni_obj.has_pdb_entry = True
        if exp == 'c3':
            pdb_obj = sequence_objects.UniProtData.PDBData(uni_obj, '2a73',
                                                           read_dssp_columns.get_dssp_data(dir_input, "2a73"))

        elif exp == 'c3b':
            pdb_obj = sequence_objects.UniProtData.PDBData(uni_obj, '2i07',
                                                           read_dssp_columns.get_dssp_data(dir_input, "2i07"))
        uni_obj_list = sequence_tools.biopy_seq_align([uni_obj])
        pdb_obj.is_seq_matched = True  # setting this to true starts some on match calculations in the pdb_object
        seq_obj = pdb_obj.SequenceData(pdb_obj)
        uni_to_pdb = seq_obj.uni_to_pdb_conv_dict
        return {l: pdb_obj.get_acc()[uni_to_pdb[l - 1]] for l in lys_pos_list if l in uni_to_pdb}


def get_distance(pos1, pos2, exp='c3', metric="SASD"):
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


def get_klys_asa(lys_pos, exp):
    if exp == 'c3':
        asa_dict = asa_dict_c3
    elif exp == 'c3b':
        asa_dict = asa_dict_c3b
    # first normalize in [0, 1]
    # then scale to [min_range, max_range]
    min_range = 0.0001
    max_range = 0.1
    val_array = np.array(list(asa_dict.values()))
    m = min(val_array)
    range1 = max(val_array) - m
    if lys_pos in asa_dict:
        asa = asa_dict[lys_pos]
        asa = (asa - m) / range1
        range2 = max_range - min_range
        return (asa * range2) + min_range
    print(f"Lys at pos {lys_pos} has no entry in asa dict. Assigning default value")
    return 0.005


# normal distribution centered around 23 A
def get_kon_xl(dist, mu=23, sigma=5):
    if dist > 35 or dist < 5:
        return 0
    return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((dist - mu) / sigma) ** 2)


# %%
class AllXLReactions(sbml_xl.AllXLReactionsNoDiff):

    def create_xl_trans(self):
        pass
    #
    def create_xl(self):
        pass

    def add_lys(self):
        lys_list = []
        for pos in self.params[const.D_POSITION_LYSINE]:
            user_data_dict = sbml_xl.get_user_data_dict(s_type=const.S_LYS, s_pos=pos, s_prot='C3')
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
                get_klys_asa(user_data_dict[const.D_LOCATION_LIST][0][1], exp=self.params[const.S_EXP])
                # get_fixed_reactivity()#
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
                dist = get_distance(pos_lys, pos_mono, exp=self.params[const.S_EXP])
                if dist > 0 and dist < 35:
                    p_kon_xl = self.min_model.sbml_model.addParameter(
                        user_data_dict[const.D_ID],
                        get_kon_xl(dist),
                        units=const.S_UNIT_LITRE_PER_MOLE_PER_SECOND,
                    )
                    p_kon_xl.UserData = user_data_dict
                    unique_id_kon_dict[location_id] = p_kon_xl
        return unique_id_kon_dict

    def add_xl_trans_params(self):
        return {}


# %%
pos_lys_c3 = get_lys_pos('c3')
pos_lys_c3b = get_lys_pos('c3b')
asa_dict_c3 = get_asa_dict(pos_lys_c3, exp='c3')
asa_dict_c3b = get_asa_dict(pos_lys_c3b, exp='c3b')

# %%
min_model_c3 = sbml_xl.MinimalModel(kinetic_params={'kh': 1e-5, 'koff': 1e-1, 'kon': 1e-3}, c_linker=111)
# %%
params_c3 = {
    const.D_POSITION_LYSINE: pos_lys_c3,
    const.S_EXP: 'c3',
}
all_reactions_c3 = AllXLReactions(min_model_c3, params_c3)
min_model_c3.add_compartments(2)
# %%
print("Write model xml")
with open("../output/model_c3_asa_mono_only_simple_comp.xml", "w") as f:
    f.write(min_model_c3.sbml_model.toSBML())

# %%
min_model_c3b = sbml_xl.MinimalModel(kinetic_params={'kh': 1e-5, 'koff': 1e-1, 'kon': 1e-3}, c_linker=111)
# %%
params_c3b = {
    const.D_POSITION_LYSINE: pos_lys_c3b,
    const.S_EXP: "c3b",
}
all_reactions_c3b = AllXLReactions(min_model_c3b, params_c3b)
min_model_c3b.add_compartments(2)
# %%
print("Write model xml")
with open("../output/model_c3b_asa_mono_only_simple_comp.xml", "w") as f:
    f.write(min_model_c3b.sbml_model.toSBML())

# %%
