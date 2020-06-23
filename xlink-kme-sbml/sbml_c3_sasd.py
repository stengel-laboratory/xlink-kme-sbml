# %%
import numpy as np
import pandas as pd
import sbml_xl
from Bio import SeqIO
import sys
sys.path.append('/home/kai/Coding/xlink_stats/src/')
import sequence_objects
import read_all_fasta
import read_dssp_columns
import sequence_tools

# %%
s_eu_dist = "eucDist"
s_uxid = "uxid"
s_sasd = "SASD"
s_energy = "energy"
s_exp = "exp_name"
s_pos1 = 'pos1'
s_pos2 = 'pos2'

df_dist = pd.read_csv(
    "/home/kai/Projects/c3_v2/xlink_energies/energy_dist_combined.csv"
)
fasta_c3 = "/home/kai/Projects/c3_v2/xlink_energies/c3.fasta"

split = df_dist[s_uxid].str.split(":", expand=True)
df_dist['pos1'] = split[1]
df_dist['pos2'] = split[4]

# %%
def get_lys_pos():
    with open(fasta_c3, "r") as handle:
        c3_record = list(SeqIO.parse(handle, "fasta"))[0]
        return [n+1 for n,c in enumerate(c3_record.seq) if c == 'K']

def get_asa_dict(lys_pos_list, exp):
    with open(fasta_c3, "r") as handle:
        c3_record = list(SeqIO.parse(handle, "fasta"))[0]
        uni_obj = sequence_objects.UniProtData(c3_record.id, str(c3_record.seq), {c3_record.id: c3_record.description})
        uni_obj.has_pdb_entry = True
        if exp == 'c3':
            pdb_obj = sequence_objects.UniProtData.PDBData(uni_obj, '2a73', read_dssp_columns.get_dssp_data("/media/nas2/Users/Sina/dssp/", "2a73"))
        elif exp == 'c3b':
            pdb_obj = sequence_objects.UniProtData.PDBData(uni_obj, '2i07', read_dssp_columns.get_dssp_data("/media/nas2/Users/Sina/dssp/", "2i07"))
        uni_obj_list = sequence_tools.biopy_seq_align([uni_obj])
        pdb_obj.is_seq_matched = True  # setting this to true starts some on match calculations in the pdb_object
        seq_obj = pdb_obj.SequenceData(pdb_obj)
        uni_to_pdb = seq_obj.uni_to_pdb_conv_dict
        return {l:pdb_obj.get_acc()[uni_to_pdb[l-1]] for l in lys_pos_list if l in uni_to_pdb}

def get_distance(pos1, pos2, exp='c3', metric="SASD"):
    if pos1 != pos2:
        d = df_dist[(df_dist[s_exp] == exp) & (df_dist[s_pos1] == str(pos1)) & (df_dist[s_pos2] == str(pos2))][metric]
        if len(d) > 0:
            return d.iloc[0]
        else:
            d = df_dist[(df_dist[s_exp] == exp) & (df_dist[s_pos1] == str(pos2)) & (df_dist[s_pos2] == str(pos1))][metric]
            if len(d) > 0:
                return d.iloc[0]
    return -1


def get_fixed_reactivity():
    return 0.1  # np.random.exponential(scale=0.1)

def get_klys_asa(lys_pos, exp):
    if exp == 'c3':
        asa_dict = asa_dict_c3
    elif exp == 'c3b':
        asa_dict = asa_dict_c3b
    # first normalize in [0, 1]
    # then scale to [min_range, max_range]
    min_range = 0.01
    max_range = 1
    val_array = np.array(list(asa_dict.values()))
    m = min(val_array)
    range1 = max(val_array) - m
    if lys_pos in asa_dict:
        asa = asa_dict[lys_pos]
        asa = (asa - m) / range1
        range2 = max_range - min_range
        return (asa*range2) + min_range
    print(f"Lys at pos {lys_pos} has no entry in asa dict. Assigning default value")
    return 0.1

# normal distribution centered around 23 A
def get_kon_xl(dist, mu=23, sigma=5):
    if dist > 35 or dist < 5:
        return 0
    return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((dist - mu) / sigma) ** 2)


# %%
class AllXLReactionsIMP(sbml_xl.AllXLReactions):

    def add_lys(self):
        lys_list = []
        for pos in self.params["pos_lys"]:
            user_data_dict = sbml_xl.get_user_data_dict(s_type=sbml_xl.s_str_lys,s_pos=pos, s_prot='C3')
            s = self.min_model.sbml_model.addSpecies(user_data_dict[sbml_xl.d_id], 1)
            s.setSpeciesType(user_data_dict[sbml_xl.d_type])
            s.UserData = user_data_dict
            lys_list.append(s)
        return lys_list

    def add_lys_params(self):
        mono_trans_to_param_dict = {}
        for mono_trans in self.react_mono_trans.products:
            user_data_dict = sbml_xl.get_user_data_dict(
                s_type="klys", s_precursor_list=[mono_trans]
            )
            p_klys = self.min_model.sbml_model.addParameter(
                user_data_dict[sbml_xl.d_id], get_klys_asa(user_data_dict[sbml_xl.d_location_list][0][1], exp=self.params['exp'])
            )  # type: tesbml.libsbml.Parameter
            p_klys.UserData = user_data_dict
            mono_trans_to_param_dict[mono_trans.getId()] = p_klys
        return mono_trans_to_param_dict

    def add_xl_trans_params(self):
        unique_id_kon_dict = {}
        for lys in self.species_lys:
            location_id_lys = lys.UserData[sbml_xl.d_location_id]
            for mono in self.react_mono.products:
                location_id_mono = mono.UserData[sbml_xl.d_location_id]
                if location_id_lys == location_id_mono:
                    continue
                user_data_dict = sbml_xl.get_user_data_dict(
                    s_type="kon_xl", s_precursor_list=[lys, mono]
                )
                location_id = user_data_dict[sbml_xl.d_location_id]
                if location_id in unique_id_kon_dict:
                    continue
                pos_lys = lys.UserData[sbml_xl.d_location_list][0][1]
                pos_mono = mono.UserData[sbml_xl.d_location_list][0][1]
                dist = get_distance(pos_lys, pos_mono, exp=self.params['exp'])
                if dist > 0:
                    p_kon_xl = self.min_model.sbml_model.addParameter(
                        user_data_dict[sbml_xl.d_id],
                        get_kon_xl(dist),
                        units="litre_per_mole_per_second",
                    )
                    p_kon_xl.UserData = user_data_dict
                    unique_id_kon_dict[location_id] = p_kon_xl
        return unique_id_kon_dict
# %%
pos_lys = get_lys_pos()
asa_dict_c3 = get_asa_dict(pos_lys, exp='c3')
asa_dict_c3b = get_asa_dict(pos_lys, exp='c3b')

# %%
min_model_c3 = sbml_xl.MinimalModel(c_xl=111, kh=1e-6)
# %%
params_c3 = {
    "pos_lys": pos_lys,
    "exp": 'c3',
}
all_reactions_c3 = AllXLReactionsIMP(min_model_c3, params_c3)
# %%
with open("/home/kai/Coding/xlink_kme/xlink_kme/model_c3_asa.xml", "w") as f:
    f.write(min_model_c3.sbml_model.toSBML())



# %%
min_model_c3b = sbml_xl.MinimalModel(c_xl=111, kh=1e-6)
# %%
params_c3b = {
    "pos_lys": pos_lys,
    "exp": "c3b",
}
all_reactions_c3b = AllXLReactionsIMP(min_model_c3b, params_c3b)
# %%
with open("/home/kai/Coding/xlink_kme/xlink_kme/model_c3b_asa.xml", "w") as f:
    f.write(min_model_c3b.sbml_model.toSBML())


# %%
