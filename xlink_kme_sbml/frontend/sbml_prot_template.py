#!/usr/bin/env python
# %%
import sys

import numpy as np
import pandas as pd
import argparse

from prot_tools import prot_lib
from xlink_kme_sbml.library import sbml_constants as const, sbml_xl

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
parser.add_argument('input', action="store", default=None, type=str, nargs='+',
                    help="List of input files separated by spaces. Minimum requirement is a fasta file. "
                         ""
                         "Further possible files: pKa file containing lysine pKa values. "
                         "sasa file containing the solvent-accessible surface area. "
                         "Distance file containing crosslink distances. "
                         "Model params as obtained by a previous run can be used as input if they match the system. ")
parser.add_argument('-c', '--crosslinker_conc', action="store_true", default=-1,
                    help="Crosslinker concentration. If none is given (i.e. the default of -1) then "
                         "equimolar concentration with respect to the lysines is used.")
parser.add_argument('-op', '--outname_params', type=str,
                    default='model_params.csv', help="Name of model parameters file")
parser.add_argument('-om', '--outname_model', type=str,
                    default='model_sbml.xml', help="Name of the sbml model output file")
args = parser.parse_args()


# %%
class AllXLReactions(sbml_xl.AllXLReactionsNoDiff):

    def add_lys(self):
        df_react = self.params[const.D_REACTIVITY_DATA_MONO]
        lys_list = []
        for prot, pdb_chain, pos, pka in zip(df_react[prot_lib.COL_UNIPROT_ID_FULL],
                                             df_react[prot_lib.COL_PDB_CHAIN_ID], df_react[prot_lib.COL_POS],
                                             df_react[prot_lib.COL_VALUE_NORMALIZED]):
            prot = prot.split("|")[1]
            user_data_dict = sbml_xl.get_user_data_dict(s_type=const.S_LYS, s_pos=pos, s_prot=prot, s_chain=pdb_chain)
            s = self.min_model.sbml_model.addSpecies(user_data_dict[const.D_ID], 1, )
            s.setSpeciesType(user_data_dict[const.D_TYPE])
            s.UserData = user_data_dict
            lys_list.append(s)
        return lys_list

    def add_lys_params(self):
        lys_to_param_dict = {}
        df_react = self.params[const.D_REACTIVITY_DATA_MONO]
        for lys in self.species_lys:
            user_data_dict = sbml_xl.get_user_data_dict(
                s_type=const.S_K_LYS, s_precursor_list=[lys]
            )
            pos = user_data_dict[const.D_POSITION_LIST][0]
            prot = user_data_dict[const.D_PROTEIN_LIST][0]
            chain = user_data_dict[const.D_CHAIN]
            val = df_react.loc[(df_react[prot_lib.COL_POS] == pos) & (df_react[prot_lib.COL_PDB_CHAIN_ID] == chain) & (
                        df_react[prot_lib.COL_UNIPROT_ID_FULL].str.contains(prot)), prot_lib.COL_VALUE_NORMALIZED].iloc[0]
            p_klys = self.min_model.sbml_model.addParameter(
                user_data_dict[const.D_ID],  # parameter itself
                val,
                # parameter value
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
                # skip already created crosslinks
                if location_id in unique_id_kon_dict:
                    continue
                pos_lys = lys.UserData[const.D_POSITION_LIST][0]
                pos_mono = mono.UserData[const.D_POSITION_LIST][0]
                react_df = self.params[const.D_REACTIVITY_DATA_XL]
                xl_react = get_xl_reactivity_from_df(df_react=react_df, pos1=pos_lys, pos2=pos_mono,
                                                     val_col=prot_lib.COL_VALUE_NORMALIZED)
                if xl_react > 0:
                    p_kon_xl = self.min_model.sbml_model.addParameter(
                        user_data_dict[const.D_ID],
                        xl_react,
                        units=const.S_UNIT_LITRE_PER_MOLE_PER_SECOND,
                    )
                    p_kon_xl.UserData = user_data_dict
                    unique_id_kon_dict[location_id] = p_kon_xl
        return unique_id_kon_dict


class MonoOnly(AllXLReactions):
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


def get_xl_reactivity_from_df(df_react, pos1, pos2, val_col):
    if pos1 != pos2:
        d = df_react.loc[(df_react[prot_lib.COL_POS_1] == pos1) & (df_react[prot_lib.COL_POS_2] == pos2), val_col]
        if len(d) > 0:
            return d.iloc[0]
        else:
            d = df_react.loc[(df_react[prot_lib.COL_POS_1] == pos2) & (df_react[prot_lib.COL_POS_2] == pos1),
                             val_col]
            if len(d) > 0:
                return d.iloc[0]
    return -1


def get_fixed_reactivity():
    return 0.01


def get_random_reactivity(scale=0.1):
    return np.random.exponential(scale=scale)


def get_lys_reactivity_in_range(lys_pos_list, df, val_col, min_range=0.0001, max_range=0.1, default_value_scale=1.05):
    """
    Given a dataframe with some experimental data like the SASD this function will first normalize in [0, 1]
    and then scale to [min_range, max_range]. A lysine position list is used to find lysines missing from the
    experimental data and assign a default for them. The default value is computed using the minimal normalized &
    scaled value multiplied by the default_value_scale parameter.

    :param lys_pos_list: a list with all lysine positions; ideally obtained from a fasta file
    :param df: dataframe containing at least lysine positions and a column with an experimental measurement
    :param val_col: name of the column containing the experimental measurement
    :param min_range: the minimal value of the scaled range
    :param max_range: the maximal value of the scaled range
    :param default_value_scale: missing value are computed by the minimal scaled value multiplied by this parameter
    """
    # normalize between 0 and 1
    m = df[val_col].min()
    range_val = df[val_col].max() - m
    df[prot_lib.COL_VALUE_NORMALIZED] = (df[val_col] - m) / range_val
    # scale between min and max range
    range_norm = max_range - min_range
    df[prot_lib.COL_VALUE_NORMALIZED] = df[prot_lib.COL_VALUE_NORMALIZED] * range_norm + min_range
    # lysines missing from the pdb structure are not handled at the moment
    # df = _assign_missing_lysine_reactivies(df, lys_pos_list=lys_pos_list, default_value_scale=default_value_scale)
    return df


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


def get_xl_dist(x, mu1=6.4, mu2=14.5, sigma1=1.15, sigma2=3.3):
    return 0.11 * get_normal_dist(x, mu=mu1, sigma=sigma1) + 0.89 * get_normal_dist(x, mu=mu2, sigma=sigma2)


def get_xl_reactivity(df, val_col, min_range=0.0001, max_range=0.1):
    df[prot_lib.COL_VALUE_NORMALIZED] = df[val_col].transform(get_xl_dist)
    range_norm = max_range - min_range
    df[prot_lib.COL_VALUE_NORMALIZED] = df[prot_lib.COL_VALUE_NORMALIZED] * range_norm + min_range
    return df


def get_input(inp_list):
    file_dict = {}
    for inp in inp_list:
        if prot_lib.is_fasta_file(inp):
            file_dict[prot_lib.COL_FASTA_FILE] = inp
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
    assert file_dict[prot_lib.COL_FASTA_FILE], "ERROR: At least a fasta file is required"
    return file_dict


def get_pos_value_dict(lys_pos_list, df, value_name):
    df_g = df.groupby(prot_lib.COL_PDB_CHAIN_ID).apply(
        lambda x: dict(zip(x[prot_lib.COL_POS], x[value_name]))).reset_index(name=value_name)
    val_dict = dict(zip(df_g[prot_lib.COL_PDB_CHAIN_ID], df_g[value_name]))
    return val_dict


def get_lys_reactivity_df(protein_id, lys_pos_list, df_sasa=None, df_pka=None):
    df_lsy_react_sasa = None
    df_lsy_react_pka = None
    df_lys_react_combined = None

    if df_sasa is not None:
        df_lsy_react_sasa = get_lys_reactivity_in_range(lys_pos_list=lys_pos_list, df=df_sasa,
                                                        val_col=prot_lib.COL_SASA)

    if df_pka is not None:
        df_lsy_react_pka = get_lys_reactivity_in_range(lys_pos_list=lys_pos_list, df=df_pka,
                                                       val_col=prot_lib.COL_PKA)

    if df_lsy_react_sasa is not None or df_lsy_react_pka is not None:
        df_lys_react_combined = pd.concat([df_lsy_react_sasa, df_lsy_react_pka],
                                          join='inner')  # only concat common cols
        if df_lsy_react_sasa is not None and df_lsy_react_pka is not None:
            df_lys_react_combined = pd.DataFrame(
                df_lys_react_combined.groupby([prot_lib.COL_PDB_FILE, prot_lib.COL_PDB_CHAIN_ID, prot_lib.COL_POS])[
                    prot_lib.COL_VALUE_NORMALIZED].mean()).reset_index()

    if df_lys_react_combined is None:
        random_react_dict = {pos: get_random_reactivity(scale=0.01) for pos in lys_pos_list}
        df_lys_react_combined = pd.DataFrame({prot_lib.COL_POS: random_react_dict.keys(),
                                              prot_lib.COL_VALUE_NORMALIZED: random_react_dict.values()})
    df_lys_react_combined[prot_lib.COL_UNIPROT_ID_FULL] = protein_id
    return df_lys_react_combined


def get_random_xl_reactivity_df(lys_pos_dict):
    # assign random reactivities to a subset (skip_percentage) of all possible crosslinks
    react_list = []
    for prot, lys_pos_list in lys_pos_dict:
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
                                        prot_lib.COL_UNIPROT_ID_FULL: prot,
                                        }))
    return pd.concat(react_list)


def main():
    inp_dict = get_input(args.input)
    df_sasd = None
    df_sasa = None
    df_pka = None
    fasta_records = prot_lib.get_fasta_records(inp_dict[prot_lib.COL_FASTA_FILE])
    lys_pos_dict = prot_lib.get_prot_lys_pos_dict(fasta_records)
    df_fasta_to_pdb = None
    if prot_lib.COL_OFFSET in inp_dict:
        df_fasta_to_pdb = pd.read_csv(inp_dict[prot_lib.COL_OFFSET])
    assert df_fasta_to_pdb is not None, "Fasta to pdb offsets file is required"
    if prot_lib.COL_PKA in inp_dict:
        df_pka = pd.read_csv(inp_dict[prot_lib.COL_PKA])
    if prot_lib.COL_SASA in inp_dict:
        df_sasa = pd.read_csv(inp_dict[prot_lib.COL_SASA])
    if prot_lib.COL_DIST_SASD in inp_dict:
        df_sasd = pd.read_csv(inp_dict[prot_lib.COL_DIST_SASD])

    df_lys_react_list = []
    for uniprot_id in df_fasta_to_pdb[prot_lib.COL_UNIPROT_ID_FULL].unique():
        lys_pos_list = lys_pos_dict[uniprot_id]
        df_lys_react_list.append(
            get_lys_reactivity_df(protein_id=uniprot_id, lys_pos_list=lys_pos_list, df_pka=df_pka, df_sasa=df_sasa))
    df_lys_react_combined = pd.concat(df_lys_react_list)

    if df_sasd is not None:
        df_xl_react_sasd = get_xl_reactivity(df_sasd, val_col=prot_lib.COL_DIST_SASD)
    else:
        df_xl_react_sasd = get_random_xl_reactivity_df(lys_pos_dict=lys_pos_dict)

    if args.crosslinker_conc == -1:
        xl_conc = len(lys_pos_list)
        f"Setting equimolar crosslinker concentration: {xl_conc}"
    else:
        xl_conc = args.crosslinker_conc

    min_model = sbml_xl.MinimalModel(kinetic_params={'kh': 1e-5, 'koff': 1e-1, 'kon': 1e-3}, c_linker=xl_conc)
    params = {
        const.D_REACTIVITY_DATA_MONO: df_lys_react_combined,
        const.D_REACTIVITY_DATA_XL: df_xl_react_sasd,
    }
    all_reactions = AllXLReactions(min_model, params)

    if prot_lib.COL_PARAM in inp_dict:
        df_params = pd.read_csv(inp_dict[prot_lib.COL_PARAM])
        min_model.from_df_params(df_params)
    # %%

    with open(args.outname_model, "w") as f:
        f.write(min_model.sbml_model.toSBML())
    print(f"Write model xml to {args.outname_model}")
    min_model.to_df_params().to_csv(args.outname_params, index=False)
    print(f"Write model params to {args.outname_params}")


if __name__ == "__main__":
    main()
