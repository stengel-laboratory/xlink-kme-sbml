#!/usr/bin/env python
# %%
import sys

import numpy as np
import pandas as pd
import argparse

from prot_tools import prot_lib
from xlink_kme_sbml.library import sbml_constants as const, sbml_xl

"""
2023-09 Kai-Michael Kammer
"""

desc = """
Create a kinetic crosslinking model of a protein. The final model depends on the input files provided.
"""

parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# not sure if just using type=str is better
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input files separated by spaces. A pdb file is required!"
                         ""
                         "Possible files: pKa file containing lysine pKa values. "
                         "sasa file containing the solvent-accessible surface area. "
                         "Distance file containing crosslink distances. "
                         "Model params as obtained by a previous run can be used as input if they match the system. ")
parser.add_argument('-c', '--crosslinker_conc', default=-1, type=int,
                    help="Crosslinker concentration. If none is given (i.e. the default of -1) then "
                         "half-equimolar concentration with respect to the lysines is used.")
parser.add_argument('-om', '--outname_model', type=str,
                    default='model_sbml.xml', help="Name of the sbml model output file")
parser.add_argument('-scmin', '--sasd_cutoff_min', type=int,
                    default=5, help="Minimum SASD for crosslinks")
parser.add_argument('-scmax', '--sasd_cutoff_max', type=int,
                    default=30, help="Maximum SASD for crosslinks")
args = parser.parse_args()


# %%
class AllXLReactions(sbml_xl.AllXLReactionsImplicitDiffusionSimple):

    def add_lys(self):
        df_react = self.params[const.S_REACTIVITY_DATA_MONO]
        lys_list = []
        for prot, pdb_chain, pos, pka in zip(df_react[prot_lib.COL_UNIPROT_ID],
                                             df_react[prot_lib.COL_PDB_CHAIN_ID], df_react[prot_lib.COL_POS],
                                             df_react[prot_lib.COL_VALUE_NORMALIZED]):
            user_data_dict = sbml_xl.get_user_data_dict(s_type=const.S_LYS, s_pos=pos, s_prot=prot, s_chain=pdb_chain)
            s = self.min_model.sbml_model.addSpecies(user_data_dict[const.S_ID], self.params[const.S_CONCENTRATION], )
            s.setSpeciesType(user_data_dict[const.S_TYPE])
            s.UserData = user_data_dict
            lys_list.append(s)
        return lys_list

    def add_lys_params(self):
        lys_to_param_dict = {}
        df_react = self.params[const.S_REACTIVITY_DATA_MONO]
        for lys in self.species_lys:
            user_data_dict = sbml_xl.get_user_data_dict(
                s_type=const.S_K_LYS, s_precursor_list=[lys]
            )
            pos = user_data_dict[const.S_POSITION_LIST][0]
            prot = user_data_dict[const.S_PROTEIN_LIST][0]
            chain = user_data_dict[const.S_CHAIN]
            val = df_react.loc[(df_react[prot_lib.COL_POS] == pos) & (df_react[prot_lib.COL_PDB_CHAIN_ID] == chain) & (
                        df_react[prot_lib.COL_UNIPROT_ID].str.fullmatch(prot)), prot_lib.COL_VALUE_NORMALIZED].iloc[0]
            p_klys = self.min_model.sbml_model.addParameter(
                user_data_dict[const.S_ID],  # parameter itself
                val,
                # parameter value
            )  # type: tesbml.libsbml.Parameter
            p_klys.UserData = user_data_dict
            lys_to_param_dict[user_data_dict[const.S_LOCATION_ID]] = p_klys
        return lys_to_param_dict

    def add_xl_trans_params(self):
        unique_id_kon_dict = {}
        for lys in self.species_lys:
            location_id_lys = lys.UserData[const.S_LOCATION_ID]
            for mono in self.react_mono.products:
                location_id_mono = mono.UserData[const.S_LOCATION_ID]
                if location_id_lys == location_id_mono:
                    continue
                user_data_dict = sbml_xl.get_user_data_dict(
                    s_type=const.S_K_ON_XL, s_precursor_list=[lys, mono]
                )
                location_id = user_data_dict[const.S_LOCATION_ID]
                # skip already created crosslinks
                if location_id in unique_id_kon_dict:
                    continue
                pos_lys = lys.UserData[const.S_POSITION_LIST][0]
                prot_lys = lys.UserData[const.S_PROTEIN_LIST][0]
                chain_lys = lys.UserData[const.S_CHAIN]
                pos_mono = mono.UserData[const.S_POSITION_LIST][0]
                prot_mono = mono.UserData[const.S_PROTEIN_LIST][0]
                chain_mono = mono.UserData[const.S_CHAIN]
                react_df = self.params[const.S_REACTIVITY_DATA_XL]
                xl_react = get_xl_reactivity_from_df(df_react=react_df, pos1=pos_lys, pos2=pos_mono, chain1=chain_lys,
                                                     chain2=chain_mono, prot1=prot_lys, prot2=prot_mono,
                                                     val_col=prot_lib.COL_VALUE_NORMALIZED)
                if xl_react > 0:
                    p_kon_xl = self.min_model.sbml_model.addParameter(
                        user_data_dict[const.S_ID],
                        xl_react,
                        units=const.S_UNIT_LITRE_PER_MOLE_PER_SECOND,
                    )
                    p_kon_xl.UserData = user_data_dict
                    unique_id_kon_dict[location_id] = p_kon_xl
        return unique_id_kon_dict


    # def create_xl_trans(self):
    #     pass
    # #
    # def create_xl(self):
    #     pass

class MonoReactionsOnly(AllXLReactions):
    """
    Monolink Only Reactions
    Just do nothing when overwriting the three functions below
    """

    def create_xl_trans(self):
        pass

    def create_xl(self):
        pass

    def add_xl_trans_params(self):
        return {}


def get_xl_reactivity_from_df(df_react, pos1, pos2, chain1, chain2, prot1, prot2, val_col):
    if pos1 != pos2:

        d = df_react.loc[(df_react[prot_lib.COL_POS_1] == pos1) & (df_react[prot_lib.COL_POS_2] == pos2) &
                         (df_react[prot_lib.COL_CHAIN_1] == chain1) & (df_react[prot_lib.COL_CHAIN_2] == chain2) &
                         (df_react[prot_lib.COL_UNIPROT_ID_1].str.fullmatch(prot1)) &
                         (df_react[prot_lib.COL_UNIPROT_ID_2].str.fullmatch(prot2)), val_col]
        if len(d) > 0:
            return d.iloc[0]
        else:
            d = df_react.loc[(df_react[prot_lib.COL_POS_1] == pos2) & (df_react[prot_lib.COL_POS_2] == pos1) &
                         (df_react[prot_lib.COL_CHAIN_1] == chain2) & (df_react[prot_lib.COL_CHAIN_2] == chain1) &
                         (df_react[prot_lib.COL_UNIPROT_ID_1].str.fullmatch(prot2)) &
                         (df_react[prot_lib.COL_UNIPROT_ID_2].str.fullmatch(prot1)), val_col]
            if len(d) > 0:
                return d.iloc[0]
    return -1


def get_fixed_reactivity():
    return 0.01


def get_random_reactivity(scale=0.1):
    return np.random.exponential(scale=scale)


# the range for a nhs/lysine reaction varies across orders of magnitudes
# Acetylation of protein lysines should be around a mean of 1e-3 1/M*s (range between 1e-6 to 1e-2)
# see: https://pubs.acs.org/doi/10.1021/acs.jpcb.0c02522
# while isolated NHS/lysine reactions can go up to 10 1/M*s
# see: https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1399-3011.1987.tb03319.x
# and http://www.sciencedirect.com/science/article/pii/S0003269717301112
def get_lys_reactivity_in_range(df, val_col, min_range=1e-3, max_range=1, ph=None):
    """
    Given a dataframe with some experimental data like the SASD this function will first normalize in [0, 1]
    and then scale to [min_range, max_range]. 

    :param df: dataframe containing at least lysine positions and a column with an experimental measurement
    :param val_col: name of the column containing the experimental measurement
    :param min_range: the minimal value of the scaled range
    :param max_range: the maximal value of the scaled range
    :param ph: if given, adjust lysine reactivites by pH
    """
    value_list = df[val_col].to_numpy()
    # optionally get ph scale factors
    if ph is not None:
        value_list = get_ph_scale_factor(pka=value_list, ph=ph)
    # normalize between 0 and 1
    min_val = value_list.min()
    range_val = value_list.max() - min_val
    df[prot_lib.COL_VALUE_NORMALIZED] = (value_list - min_val) / range_val
    # scale between min and max range
    range_norm = max_range - min_range
    df[prot_lib.COL_VALUE_NORMALIZED] = df[prot_lib.COL_VALUE_NORMALIZED] * range_norm + min_range
    # lysines missing from the pdb structure are not handled at the moment
    # df = _assign_missing_lysine_reactivies(df, lys_pos_list=lys_pos_list, default_value_scale=default_value_scale)
    return df

# uses the Henderson–Hasselbalch equation to approximate a pH dependent scaling factor for lysines
def get_ph_scale_factor(pka, ph):
    return (1/(1 + 10**(pka - ph)))


# not currently used; missing lysines in the pdb file are ignored
def _assign_missing_lysine_reactivies(df, lys_pos_list, default_value_scale):
    # default value for lysines not in dataframe
    reactivity_default = df[prot_lib.COL_VALUE_NORMALIZED].min() * default_value_scale
    missing_lys_dict = {}
    for lys_pos in lys_pos_list:
        if lys_pos not in list(df[prot_lib.COL_POS]):
            missing_lys_dict[lys_pos] = reactivity_default
    if missing_lys_dict:
        df_missing = pd.DataFrame(
            {prot_lib.COL_POS: missing_lys_dict.keys(), prot_lib.COL_VALUE_NORMALIZED: missing_lys_dict.values()})
        df = pd.concat([df, df_missing])
    return df


# according to https://pubs.acs.org/doi/10.1021/acs.jpcb.0c02522 the crosslinked SASD distribution
# is a linear combination of two normal distributions
# note that sigma was not given in the paper and fitted visually
def get_normal_dist(x, mu=0.0, sigma=1.0):
    return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

# correction factor adds a value of 3 (Angstrom) to both mu1 and mu2
# this corrects for the fact, that the paper uses C beta distances while jwalk calculates C alpha distances
def get_xl_dist(x, mu1=6.4, mu2=14.5, sigma1=1.15, sigma2=3.3, correction_factor=3):
    return 0.11 * get_normal_dist(x, mu=mu1+correction_factor, sigma=sigma1) + 0.89 * get_normal_dist(x, mu=mu2+correction_factor, sigma=sigma2)


def get_xl_reactivity(df, val_col, min_range=1e6, max_range=1e9):
    # first convert SASD to experimental distribution
    df[prot_lib.COL_VALUE_NORMALIZED] = df[val_col].transform(get_xl_dist)
    # then normalize between 0 and 1
    min_val = df[prot_lib.COL_VALUE_NORMALIZED].min()
    max_val = df[prot_lib.COL_VALUE_NORMALIZED].max()
    df[prot_lib.COL_VALUE_NORMALIZED] = (df[prot_lib.COL_VALUE_NORMALIZED] - min_val) / (max_val - min_val)
    # finally scale to the desired range
    range_norm = max_range - min_range
    df[prot_lib.COL_VALUE_NORMALIZED] = df[prot_lib.COL_VALUE_NORMALIZED] * range_norm + min_range
    return df


def get_input(inp_list):
    file_dict = {}
    for inp in inp_list:
        if prot_lib.is_fasta_file(inp):
            file_dict[prot_lib.COL_FASTA_FILE] = inp
        elif prot_lib.is_pdb_file(inp):
            file_dict[prot_lib.COL_PDB_FILE] = inp
        elif prot_lib.is_sasd_file(inp):
            file_dict[prot_lib.COL_DIST_SASD] = inp
        elif prot_lib.is_pka_file(inp):
            file_dict[prot_lib.COL_PKA] = inp
        elif prot_lib.is_sasa_file(inp):
            file_dict[prot_lib.COL_SASA] = inp
        elif prot_lib.is_params_file(inp):
            file_dict[prot_lib.COL_PARAM] = inp
        elif prot_lib.is_fasta_to_pdb_file(inp):
            file_dict[prot_lib.COL_OFFSET] = inp
    assert file_dict[prot_lib.COL_PDB_FILE], "ERROR: At least a pdb file is required"
    return file_dict


def get_pos_value_dict(lys_pos_list, df, value_name):
    df_g = df.groupby(prot_lib.COL_PDB_CHAIN_ID).apply(
        lambda x: dict(zip(x[prot_lib.COL_POS], x[value_name]))).reset_index(name=value_name)
    val_dict = dict(zip(df_g[prot_lib.COL_PDB_CHAIN_ID], df_g[value_name]))
    return val_dict


def get_lys_reactivity_df(chain_to_uni_id_dict=None, df_sasa=None, df_pka=None, ph=None):
    df_lsy_react_sasa = None
    df_lsy_react_pka = None
    df_lys_react_combined = None

    if df_sasa is not None:
        df_lsy_react_sasa = get_lys_reactivity_in_range(df=df_sasa,
                                                        val_col=prot_lib.COL_SASA)

    if df_pka is not None:
        df_lsy_react_pka = get_lys_reactivity_in_range(df=df_pka,
                                                       val_col=prot_lib.COL_PKA,
                                                       ph=ph)

    if df_lsy_react_sasa is not None or df_lsy_react_pka is not None:
        df_lys_react_combined = pd.concat([df_lsy_react_sasa, df_lsy_react_pka],
                                          join='inner')  # only concat common cols
        if df_lsy_react_sasa is not None and df_lsy_react_pka is not None:
            df_lys_react_combined = pd.DataFrame(
                df_lys_react_combined.groupby([prot_lib.COL_PDB_FILE, prot_lib.COL_PDB_CHAIN_ID, prot_lib.COL_POS])[
                    prot_lib.COL_VALUE_NORMALIZED].mean()).reset_index()
    if chain_to_uni_id_dict:
        df_lys_react_combined[prot_lib.COL_UNIPROT_ID] = df_lys_react_combined[prot_lib.COL_PDB_CHAIN_ID].map(chain_to_uni_id_dict)
    return df_lys_react_combined


def get_random_xl_reactivity_df(lys_pos_dict, chain_to_uni_dict):
    # assign random reactivities to a subset (skip_percentage) of all possible crosslinks
    react_list = []
    for chain_id, lys_pos_list in lys_pos_dict:
        skip_percentage = (len(lys_pos_list) * 5) / ((len(lys_pos_list) - 1) ** 2)
        pos1_list = []
        pos2_list = []
        for pos1 in lys_pos_list:
            for pos2 in lys_pos_list:
                if pos1 == pos2:
                    continue
                if skip_percentage <= np.random.uniform():
                    continue
                pos1_list.append(pos1)
                pos2_list.append(pos2)
        react_list.append(pd.DataFrame({prot_lib.COL_POS_1: pos1_list, prot_lib.COL_POS_2: pos2_list,
                                        prot_lib.COL_VALUE_NORMALIZED: get_random_reactivity(scale=5e7),
                                        prot_lib.COL_UNIPROT_ID: chain_to_uni_dict[chain_id],
                                        prot_lib.COL_PDB_CHAIN_ID: chain_id
                                        }))
    return pd.concat(react_list)


def main():
    inp_dict = get_input(args.input)
    df_sasd = None
    df_sasa = None
    df_pka = None
    ph = 7
    # fasta_records = prot_lib.get_fasta_records(inp_dict[prot_lib.COL_FASTA_FILE])
    mol_weight_prot = prot_lib.get_pdb_mol_weight(inp_dict[prot_lib.COL_PDB_FILE])
    print(f"Scaling concentration by molecular weight: {mol_weight_prot}")
    pdb_chain_to_uni_id_dict, lys_pos_dict = prot_lib.get_prot_lys_pos_dict_pdb(inp_dict[prot_lib.COL_PDB_FILE])
    if prot_lib.COL_PKA in inp_dict:
        df_pka = pd.read_csv(inp_dict[prot_lib.COL_PKA])
    if prot_lib.COL_SASA in inp_dict:
        df_sasa = pd.read_csv(inp_dict[prot_lib.COL_SASA])
    if prot_lib.COL_DIST_SASD in inp_dict:
        df_sasd = pd.read_csv(inp_dict[prot_lib.COL_DIST_SASD])
        len_before = len(df_sasd)
        df_sasd = df_sasd[(df_sasd[prot_lib.COL_DIST_SASD] >= args.sasd_cutoff_min) &
                          (df_sasd[prot_lib.COL_DIST_SASD] <= args.sasd_cutoff_max)]
        len_after = len(df_sasd)
        print(f"Removed {len_before-len_after} of {len_before} crosslinks,"
              f" leaving {len_after} crosslinks above the cutoff distance of {args.sasd_cutoff_min} and below "
              f"{args.sasd_cutoff_max}")
    df_lys_react_combined = get_lys_reactivity_df(chain_to_uni_id_dict=pdb_chain_to_uni_id_dict, df_pka=df_pka, df_sasa=df_sasa, ph=ph)
    df_lys_react_combined = df_lys_react_combined[df_lys_react_combined[prot_lib.COL_POS] > -1].reset_index(drop=True) # some pdbs contain negative sequence numbers; drop them
    num_lys = len(df_lys_react_combined)
    if df_sasd is not None:
        df_xl_react_sasd = get_xl_reactivity(df_sasd, val_col=prot_lib.COL_DIST_SASD)
    else:
        df_xl_react_sasd = get_random_xl_reactivity_df(lys_pos_dict=lys_pos_dict, chain_to_uni_dict=pdb_chain_to_uni_id_dict)
    if args.crosslinker_conc == -1:
        xl_conc = num_lys/2
        print(f"Setting half-equimolar crosslinker concentration: {num_lys}/2")
    else:
        xl_conc = args.crosslinker_conc

    # adjust lys/crosslinker concentration, so that lys concentration is 1 g/L
    # as tellurium wants concentration in mole/liter we divide 1 g/L by the molecular weight M (g/mole)
    xl_conc /= mol_weight_prot # set linker conc to 0,5 g/L * lysine number
    lys_conc = 1 / mol_weight_prot # set prot conc to 1 g/L
    
    # uncomment for constant molecular concentration instead of constanst mass concentration
    # lys_conc = 5.5e-6
    # xl_conc = num_lys*lys_conc/2

    # optionally adjust hydrolysis by pH
    # since our value already orients itself at literature values around pH 7, we leave it fixed
    kh = 2.5e-5
    # kh *= get_ph_scale_factor(pka=6, ph=7) # uncomment this line to scale hydrolysis by pH using the Henderson-Hasselbalch equation; a reasonable pKa for NHS is ~6 (D. E. Ames and T. F. Grey, The synthesis of some Nhydroxyimides, J. Chem. Soc., 1955, 631–636.); however, some further work is necessary for direct application to hydrolysis and it may make more sense to find reasonable literature values instead
    min_model_xl = sbml_xl.MinimalModel(kinetic_params={const.S_K_HYDROLYSIS: kh, const.S_K_OFF: 1e9, const.S_K_ON: 1e7}, c_linker=xl_conc, mol_weight=mol_weight_prot, num_lys=num_lys, ph=ph)
    params_xl = {
        const.S_REACTIVITY_DATA_MONO: df_lys_react_combined,
        const.S_REACTIVITY_DATA_XL: df_xl_react_sasd,
        const.S_CONCENTRATION: lys_conc,
    }
    min_model_mono = sbml_xl.MinimalModel(kinetic_params={const.S_K_HYDROLYSIS: kh, const.S_K_OFF: 1e9, const.S_K_ON: 1e7}, c_linker=xl_conc, mol_weight=mol_weight_prot, num_lys=num_lys, ph=ph)
    params_mono = {
        const.S_REACTIVITY_DATA_MONO: df_lys_react_combined,
        const.S_CONCENTRATION: lys_conc,
    }
    all_reactions = AllXLReactions(min_model_xl, params_xl)
    mono_reactions = MonoReactionsOnly(min_model_mono, params_mono)

    if prot_lib.COL_PARAM in inp_dict:
        df_params = pd.read_csv(inp_dict[prot_lib.COL_PARAM])
        min_model_xl.from_df_params(df_params)
    # %%

    with open(args.outname_model, "w") as f:
        f.write(min_model_xl.sbml_model.toSBML())
    print(f"Write model xml to {args.outname_model}")
    outname_mono_only = args.outname_model.replace(".xml", "_mono_only.xml")
    with open(outname_mono_only, "w") as f:
        f.write(min_model_mono.sbml_model.toSBML())
    print(f"Write mono-only model xml to {outname_mono_only}")
    outname_params = args.outname_model.replace(".xml", "_params.csv")
    min_model_xl.to_df_params().to_csv(outname_params, index=False)
    print(f"Write model params to {outname_params}")


if __name__ == "__main__":
    main()
