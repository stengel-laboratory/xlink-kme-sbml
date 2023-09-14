from __future__ import print_function, absolute_import

from typing import List

import pandas as pd
import tellurium as te
import xlink_kme_sbml.library.sbml_constants as const
import numpy as np
from timeit import default_timer as timer
import itertools as it

"""
2023-09 Kai-Michael Kammer
Functions to facilitate crosslink model creation
"""


def load_model(exp, exp_type):
    print(f"Loading {exp}_{exp_type} model")
    rr = te.loadSBMLModel(f"../output/model_{exp}_{exp_type}.xml")
    return rr


def load_model_path(path):
    print(f"Loading model at {path}")
    rr = te.loadSBMLModel(path)
    return rr


def write_antimony(exp, exp_type, rr):
    print("Writing Antimony")
    with open(
            f"../output/model_{exp}_antimony_{exp_type}.txt", "w"
    ) as f:
        f.write(rr.getAntimony())

    # rr.draw(layout="dot")


def simulate(rr):
    print("Starting Simulation")
    results = rr.simulate(0, 50000000)
    print("convergence", max(abs(rr.dv())))
    return results


def write_results(exp, exp_type, df_ff):
    df_ff.to_csv(f'../output/{exp}_final_frame_{exp_type}.csv', index=False)


def get_final_frame(results):
    df_res = pd.DataFrame(results, columns=results.colnames)
    return df_res.tail(1)


def load_and_simulate(exp, exp_type):
    rr = load_model(exp, exp_type)
    results = simulate(rr)
    return rr, results


def __run_sim_runner(rr, var_name, val_name, exp_name, mapper_dict, min_simulation_time, sim_vars=None,
                     scale_factor=None):
    rr.resetToOrigin()
    if sim_vars:
        for var in sim_vars.keys():
            # if hasattr(rr, var):
            setattr(rr, var, sim_vars[var])
    # crosslinker_start = rr.Crosslinker
    lys_list = [a for a in dir(rr) if const.S_LYS in a]
    lys_initial_conc = getattr(rr, lys_list[0])
    lys_no = getattr(rr, const.S_NUM_LYS, None)
    mol_weight = getattr(rr, const.S_MOLECULAR_WEIGHT, None)
    crosslinker_initial_conc = getattr(rr, const.S_CROSSLINKER)
    start = timer()
    try:
        results = rr.simulate(0, min_simulation_time, points=2)
    except RuntimeError as e:
        print(f"Simulation error for {var_name}={val_name}.")
        print(e)
        return None
    end = timer()
    minimal_convergence = max(abs(rr.dv()))
    print(
        f"Simulation for {var_name}={val_name} took {int(end - start)} seconds with a convergence of {minimal_convergence}.")
    # artificial convergence threshold
    if minimal_convergence > 1e-10:
        print("WARNING: Simulation might not have converged!")
    # while minimal_convergence > 1e-18:
    #     print(f"Simulation has not converged; increasing time by a factor of 10 for {variable}={value}")
    #     return _explore_variable_runner(rr, variable, value, exp_name, mapper_dict, min_simulation_time * 10)
    df_res = get_final_frame(results)
    df_res = prepare_df(df_res, exp_name, mapper_dict)
    if scale_factor:
        df_res[const.S_VALUE_SCALED] = df_res[const.S_VALUE] * scale_factor
    else:
        if mol_weight:
            df_res[const.S_VALUE_SCALED] = df_res[const.S_VALUE] * mol_weight
    if lys_no:
        df_res[const.S_RATIO_CROSSLINKER_LYSINE] = crosslinker_initial_conc / (lys_initial_conc * lys_no)
    rr.resetToOrigin()
    return df_res


def _explore_variable_runner(rr, variable, value, exp_name, mapper_dict, min_simulation_time, custom_vars=None,
                             scale_factor=None):
    sim_vars = {variable: value}
    if custom_vars:
        sim_vars = {**custom_vars, variable: value}
    df_res = __run_sim_runner(rr, var_name=variable, val_name=value, exp_name=exp_name, mapper_dict=mapper_dict,
                              min_simulation_time=min_simulation_time, sim_vars=sim_vars, scale_factor=scale_factor,
                              )
    if df_res is not None:
        df_res[variable] = value
    return df_res


def _explore_variable_runner_mult(rr, var_val_dict, exp_name, mapper_dict, min_simulation_time, sim_vars=None,
                                  scale_factor=None):
    var_names = list(var_val_dict.keys())
    vals = list(var_val_dict.values())
    df_res = __run_sim_runner(rr, var_name=var_names, val_name=vals, exp_name=exp_name, mapper_dict=mapper_dict,
                              min_simulation_time=min_simulation_time, sim_vars=sim_vars, scale_factor=scale_factor,
                              )
    if df_res is not None:
        for variable, value in var_val_dict.items():
            df_res[variable] = value
    return df_res


def explore_variable(rr, variable, var_range, exp_name='exp', mapper_dict=None, custom_vars=None,
                     min_simulation_time=5000000, scale_factor=None):
    df_list = []
    for value in var_range:
        df_res = _explore_variable_runner(rr=rr, variable=variable, value=value, exp_name=exp_name,
                                          mapper_dict=mapper_dict, min_simulation_time=min_simulation_time,
                                          custom_vars=custom_vars, scale_factor=scale_factor)
        if df_res is not None:
            df_list.append(df_res)
    df_final = pd.concat(df_list)
    # df_melt = pd.melt(df_final.drop(columns=['time']), id_vars=variable)
    return df_final.reset_index(drop=True)


def explore_variable_multi(rr, variable_range_dict, exp_name='exp', mapper_dict=None, custom_vars=None,
                           min_simulation_time=5000000, scale_factor=None, var_operation_dict=None):
    """
    Allows the multidimensional simulation of an arbitrary number of variables
    """

    # variables in the var range dict can either be simply overriden (var_operation_dict=None)
    # or have some mathematical operation applied on the model's default value
    def _apply_var_operation(rr, var, var_operation_dict, new_value, var_name_long=None):
        if var_operation_dict and var in var_operation_dict:
            if var_name_long:
                if hasattr(rr, var_name_long):
                    original_value = getattr(rr, var_name_long)
                else:
                    return 0
            else:
                if hasattr(rr, var):
                    original_value = getattr(rr, var)
                else:
                    return 0
            if var_operation_dict[var] == '+':
                return new_value + original_value
            elif var_operation_dict[var] == '*':
                return new_value * original_value
            elif var_operation_dict[var] == '/':
                return new_value / original_value
            elif var_operation_dict[var] == '-':
                return new_value - original_value
        return new_value

    df_list = []
    values_all_combinations = it.product(*variable_range_dict.values())
    for val_tuple in values_all_combinations:
        var_val_dict = {}
        sim_vars = {}
        rr.resetToOrigin()  # otherwise getting the original value for value operations is not reliable
        for n, var in enumerate(variable_range_dict.keys()):
            if type(var) is tuple:
                var, vars = var
                for v in vars:
                    sim_vars[v] = _apply_var_operation(rr=rr, var=var, var_operation_dict=var_operation_dict,
                                                       new_value=val_tuple[n], var_name_long=v)
            else:
                sim_vars[var] = _apply_var_operation(rr=rr, var=var, var_operation_dict=var_operation_dict,
                                                     new_value=val_tuple[n])
            var_val_dict[var] = val_tuple[n]  # contains the 'short' var names, i.e. LYS instead of [LYS_1, LYS_2, ...]
        df_res = _explore_variable_runner_mult(rr=rr, var_val_dict=var_val_dict, exp_name=exp_name,
                                               mapper_dict=mapper_dict, min_simulation_time=min_simulation_time,
                                               sim_vars={**custom_vars, **sim_vars}, scale_factor=scale_factor,
                                               )
        if df_res is not None:
            df_list.append(df_res)
    df_final = pd.concat(df_list)
    return df_final.reset_index(drop=True)


def get_mapper_dict(df_map):
    """
    Create dict from data frame in the form {simulation_name: experiment_name}
    
    """
    return dict(zip(df_map[const.S_SIM_NAME], df_map[const.S_EXP_NAME]))


def prepare_df(df: pd.DataFrame, exp_name: str, mapper_dict: dict = None) -> pd.DataFrame:
    """
    Prepare the dataframe for plotting. Expects a dataframe with only a single row, i.e. one frame of a simulation.
    Also requires a mapper dict which maps the protein names in the simulation to a UniProtID; i.e {name_sim: uni_prot_id}
    Will specifically filter for Hydrolized Monolinks and Crosslinks.
    """
    assert len(df) == 1, print(f"ERROR: Input dataframe has =! 1 row: {len(df)}")
    # filter for XLs and Monolinks
    df = df[[col for col in df.columns if const.S_REACT_MONO_HYDRO + "_" in col or const.S_REACT_XL + "_" in col]]
    # convert to long format
    df = pd.melt(df)
    # remove superfluous chars
    df[const.S_VAR] = df[const.S_VAR].str.replace("\[|\]", "", regex=True)
    # split the the species name; something like MonoHydro_C3_418 and XL_C3_904_C3_1325
    split = df[const.S_VAR].str.split("_", expand=True)
    df[const.S_LINK_TYPE] = split[0]
    if mapper_dict:
        split[1] = split[1].replace(mapper_dict)
        if const.S_REACT_XL in df[const.S_LINK_TYPE]:
            split[3] = split[3].replace(mapper_dict)
    # construct the uxID and a reversed uxID
    # first row for a crosslink and second for a monolink
    # as np.where is a if-else-statement
    if const.S_REACT_XL in df[const.S_LINK_TYPE]:
        df[const.S_UID] = np.where(df[const.S_LINK_TYPE] == const.S_REACT_XL,
                                   split[1] + ':' + split[2] + ':x:' + split[3] + ':' + split[4],
                                   split[1] + ':' + split[2])
        df[const.S_UID_REV] = np.where(df[const.S_LINK_TYPE] == const.S_REACT_XL,
                                       split[3] + ':' + split[4] + ':x:' + split[1] + ':' + split[2],
                                       split[1] + ':' + split[2])
    else:
        df[const.S_UID] = split[1] + ':' + split[2]
        df[const.S_UID_REV] = split[1] + ':' + split[2]
    df[const.S_UID_SHORT] = df[const.S_VAR].str.split('_', n=1, expand=True)[1]
    df[const.S_EXP] = exp_name
    return df


def get_rr_const(rr, constant, id_short):
    return getattr(rr, f'{constant}_{id_short}')


def get_rr_const_dict(rr, constant):
    # ugly try except workaround because Tellurium pollutes the attribute space of distinct roadrunner objects
    att_list = [a for a in dir(rr) if constant in a]
    attr_dict = {}
    # att_dict = {a: getattr(rr, a) for a in att_list} #this would work if Tellurium wasn't buggy
    for attr in att_list:
        try:
            attr_dict[attr] = getattr(rr, attr)
        except RuntimeError as re:
            print(re)
    return attr_dict


def add_initial_concentration(df, rr):
    rr.resetToOrigin()
    m_mono = df[const.S_LINK_TYPE] == const.S_REACT_MONO_HYDRO
    klys_dict = {k: get_rr_const(rr, const.S_LYS, k) for k in df.loc[m_mono, const.S_UID_SHORT].unique()}
    df_const = pd.DataFrame.from_dict(klys_dict, orient='index', columns=[const.S_CONCENTRATION]).rename_axis(
        const.S_UID_SHORT).reset_index()
    df = pd.merge(df, df_const, on=const.S_UID_SHORT, how='outer')
    df[const.S_CONCENTRATION] = df[const.S_CONCENTRATION].fillna(0)
    return df


def add_reactions_constants(df, rr):
    m_mono = df[const.S_LINK_TYPE] == const.S_REACT_MONO_HYDRO
    m_xl = df[const.S_LINK_TYPE] == const.S_REACT_XL
    klys_dict = {k: get_rr_const(rr, const.S_K_LYS, k) for k in df.loc[m_mono, const.S_UID_SHORT].unique()}
    kon_dict = {k: get_rr_const(rr, const.S_K_ON_XL, k) for k in df.loc[m_xl, const.S_UID_SHORT].unique()}
    merged_dict = {**klys_dict, **kon_dict}
    df_const = pd.DataFrame.from_dict(merged_dict, orient='index', columns=[const.S_K_GENERIC]).rename_axis(
        const.S_UID_SHORT).reset_index()
    return pd.merge(df, df_const, on=const.S_UID_SHORT)


def add_nxl(df):
    df[['prot1', 'pos1', 'prot2', 'pos2']] = df[const.S_UID_SHORT].str.split('_', expand=True)
    df_nxl1 = (df.groupby(['pos1'])[const.S_UID_SHORT].nunique() - 1).reset_index(name=const.S_NXL1)
    df_nxl2 = df.groupby(['pos2'])[const.S_UID_SHORT].nunique().reset_index(name=const.S_NXL2)
    df = pd.merge(df, df_nxl1, on='pos1')
    df = pd.merge(df, df_nxl2, on='pos2', how='outer')
    df[const.S_NXL_SUM] = df[[const.S_NXL1, const.S_NXL2]].sum(axis=1)
    return df


def add_max_connected_rate_constant(df):
    df_max_klys_p1 = df[df[const.S_LINK_TYPE] == const.S_REACT_MONO_HYDRO].groupby(['pos1'])[
        const.S_K_GENERIC].max().reset_index(
        name=const.S_K_MAX_K_LYS_P1)
    df_max_klys_p2 = df_max_klys_p1.rename(
        columns={const.S_K_MAX_K_LYS_P1: const.S_K_MAX_K_LYS_P2, 'pos1': 'pos2'})
    df_max_kon_xl_p1 = df[df[const.S_LINK_TYPE] == 'XL'].groupby('pos1')[const.S_K_GENERIC].max().reset_index(
        name=const.S_K_MAX_K_ON_XL_P1)
    df_max_kon_xl_p2 = df[df[const.S_LINK_TYPE] == 'XL'].groupby('pos2')[const.S_K_GENERIC].max().reset_index(
        name=const.S_K_MAX_K_ON_XL_P2)
    # df_max_kon_xl = pd.concat([df_max_kon_xl_p1, df_max_kon_xl_p2]).groupby('pos').max().reset_index()
    df = pd.merge(df, df_max_klys_p1, on='pos1')
    df = pd.merge(df, df_max_klys_p2, on='pos2', how='left')
    df[const.S_K_MAX_K_LYS] = df[[const.S_K_MAX_K_LYS_P1, const.S_K_MAX_K_LYS_P2]].max(axis=1)
    # df = pd.merge(df, df_max_kon_xl_mono, on=['pos1', const.S_LINK_TYPE], how='outer')
    df = pd.merge(df, df_max_kon_xl_p1, on='pos1', how='outer')
    df = pd.merge(df, df_max_kon_xl_p2, on='pos2', how='outer')
    df[const.S_K_MAX_K_ON_XL] = df[[const.S_K_MAX_K_ON_XL_P1, const.S_K_MAX_K_ON_XL_P2]].max(axis=1)
    df = df.drop(columns=[const.S_K_MAX_K_ON_XL_P1, const.S_K_MAX_K_ON_XL_P2, const.S_K_MAX_K_ON_XL_P1,
                          const.S_K_MAX_K_ON_XL_P2])
    df[const.S_K_RATE_GT_K_LYS_MAX] = df[const.S_K_GENERIC] > df[const.S_K_MAX_K_LYS]
    df[const.S_K_RATE_GT_K_ON_XL_MAX] = df[const.S_K_GENERIC] > df[const.S_K_MAX_K_ON_XL]
    return df


def get_explored_var_name(df):
    for col in df.columns:
        if col in const.var_explore_list:
            return col
    return None


def extract_aa_pos(string, split_char='_', name_prot=False, name_chain=True):
    # klys_P01024_A_264 monolink string
    # klys_P01024_A_264_P01024_A_264 crosslink string
    return_components = []
    str_comp = string.split(split_char)
    str_comp.pop(0)  # remove species name
    current_pos = 0  # 0: prot_name; 1: prot_chain; 2: pos
    while str_comp:
        next_entry = str_comp.pop(0)
        if name_prot and current_pos == 0:
            return_components.append(next_entry)
        elif name_chain and current_pos == 1:
            return_components.append(next_entry)
        elif current_pos == 2:
            return_components.append(next_entry)
        current_pos += 1
        if current_pos > 2:
            current_pos = 0
    return return_components


def extract_nums_aa_pos(string, split_char='_'):
    nums = []
    string_components = string.split(split_char)
    for comp in string_components:
        try:
            nums.append(str(int(comp)))
        except:
            None
    return nums


def get_react_dict(rr, rate_const, species_list):
    klys_dict = get_rr_const_dict(rr, rate_const)
    klys_dict_species = {}
    for k, v in klys_dict.items():
        for species in species_list:
            klys_dict_species[k.replace(rate_const, species)] = v
    return klys_dict_species


def apply_custom_params(rr, params_dict):
    rr.resetToOrigin()
    for k, v in params_dict.items():
        setattr(rr, k, v)


def get_sim_res(rr, last_frame_only=True, sim_time=1000, sim_points=100, custom_params=None, no_prot_name=True,
                scale_factor=None, exp_name='exp'):
    # params_dict = get_params_dict_equal_lys_react(rr, custom_params=custom_params)
    if custom_params:
        apply_custom_params(rr, custom_params)
    # get values for scaling to one
    # lys_list = [a for a in dir(rr) if const.S_LYS in a]
    # lys_initial_conc = getattr(rr, lys_list[0])
    mol_weight = getattr(rr, const.S_MOLECULAR_WEIGHT)
    if last_frame_only:
        points = 2
    else:
        points = sim_points
    try:
        results = rr.simulate(0, sim_time, points=points)
    except RuntimeError as e:
        print(f"Error for simulation: {e}")
        return None
    df_res = pd.DataFrame(results, columns=results.colnames)
    df_res = df_res.drop(columns=[col for col in df_res.columns if "v" in col])
    df_res.columns = df_res.columns.str.replace("\[|\]", "", regex=True)
    # klys_dict = sbml_sim_helper.get_rr_const_dict(rr, str_klys)
    # for klys in klys_dict.keys():
    #     prot_chain_pos = klys.split("_")[1:]
    #     prot_chain_pos = "_".join(prot_chain_pos)
    #     df_res[f"{str_species_mono_sum}_{prot_chain_pos}"] = df_res[f"{str_species_mono_base}_{prot_chain_pos}"] + df_res[f"{str_species_mono_hydro}_{prot_chain_pos}"]
    # melt it
    if last_frame_only:
        df_melt = pd.melt(df_res.tail(1), id_vars=const.S_TIME)
    else:
        df_melt = pd.melt(df_res, id_vars=const.S_TIME)
    # add additional columns from the variable id
    df_split_var = df_melt[const.S_VAR].str.split("_", expand=True)
    if len(df_split_var.columns) == 7:
        df_split_var.columns = [const.S_TYPE, const.S_PROT1, const.S_CHAIN_ID1, const.S_POS1, const.S_PROT2,
                                const.S_CHAIN_ID2, const.S_POS2]
    else:
        df_split_var.columns = [const.S_TYPE, const.S_PROT1, const.S_CHAIN_ID1, const.S_POS1]
    df_melt = df_melt.join(df_split_var)
    df_melt[const.S_TYPE] = df_melt[const.S_VAR].str.split("_", expand=True)[0]
    # add mono sum: monolinks (non-hydrolized) + monolinks (hydrolized)
    df_filt_mono = df_melt[df_melt[const.S_TYPE].str.fullmatch(f"{const.S_REACT_MONO}|{const.S_REACT_MONO_HYDRO}")]
    df_mono_sum = df_filt_mono.groupby([const.S_TIME, const.S_PROT1, const.S_CHAIN_ID1, const.S_POS1])[
        const.S_VALUE].sum().reset_index()
    df_mono_sum[[const.S_TYPE, const.S_VAR]] = const.S_REACT_MONO_SUM
    df_mono_sum[const.S_VAR] += "_" + df_mono_sum[const.S_PROT1] + "_" + df_mono_sum[const.S_CHAIN_ID1] + "_" + \
                                df_mono_sum[const.S_POS1]
    df_melt = pd.concat([df_melt, df_mono_sum])
    # add reaction rates
    mono_list = [m for m in df_melt[const.S_TYPE].unique() if const.S_REACT_MONO in m]
    mono_list.append(const.S_LYS)
    crosslinker_list = [m for m in df_melt[const.S_TYPE].unique() if const.S_CROSSLINKER in m]
    mono_r_dict = get_react_dict(rr, rate_const=const.S_K_LYS, species_list=mono_list)
    xl_r_dict = get_react_dict(rr, rate_const=const.S_K_ON_XL, species_list=[const.S_REACT_XL])
    crosslinker_r_dict = get_react_dict(rr, rate_const=const.S_K_HYDROLYSIS, species_list=crosslinker_list)
    df_melt[const.S_K_GENERIC] = df_melt[const.S_VAR].map({**mono_r_dict, **xl_r_dict, **crosslinker_r_dict})
    if no_prot_name:
        df_melt[const.S_VAR] = df_melt[const.S_VAR].str.replace("_[A-Z]\d+", "",
                                                                regex=True)  # remove prot name for easier readability
    if scale_factor:
        df_melt[const.S_VALUE_SCALED] = df_melt[const.S_VALUE] * scale_factor
    else:
        df_melt[const.S_VALUE_SCALED] = df_melt[const.S_VALUE] * mol_weight
    df_melt[const.S_EXP] = exp_name
    return df_melt


def get_reaction_ratio_df(df_ref: pd.DataFrame, df_list: List[pd.DataFrame], impute_missing=False):
    df_react_list = []
    df_ref = df_ref.rename(columns={const.S_VALUE: const.S_VALUE_REF, const.S_K_GENERIC: const.S_K_GENERIC_REF,
                                    const.S_VALUE_SCALED: const.S_VALUE_SCALED_REF})
    columns_to_merge = list(set(df_ref.columns) - set(const.var_supp_ref_list))
    for df in df_list:
        df_react_list.append(pd.merge(left=df_ref, right=df, on=columns_to_merge, how='outer'))
    df_concat = pd.concat(df_react_list).reset_index(drop=True)
    df_concat[const.S_REACTION_RATIO] = df_concat[const.S_VALUE] / df_concat[const.S_VALUE_REF]
    df_concat[const.S_VALUE_CHANGE] = df_concat[const.S_VALUE_SCALED] - df_concat[const.S_VALUE_SCALED_REF]
    df_concat[const.S_VALUE_INC_ABS] = abs(df_concat[const.S_VALUE_CHANGE])
    df_concat[const.S_IMPUTED] = False
    if impute_missing:
        df_concat.loc[(df_concat[const.S_VALUE] == 0) | (df_concat[const.S_VALUE_REF] == 0), const.S_IMPUTED] = True
        # df_concat[[const.S_VALUE, const.S_VALUE_SCALED, const.S_VALUE_REF, const.S_VALUE_SCALED_REF]] = df_concat[[const.S_VALUE, const.S_VALUE_SCALED, const.S_VALUE_REF, const.S_VALUE_SCALED_REF]].replace(0, 1e-5)
        median_val_scaled = df_concat[const.S_VALUE_SCALED].median()
        df_concat[[const.S_VALUE, const.S_VALUE_SCALED, const.S_VALUE_REF, const.S_VALUE_SCALED_REF]] = df_concat[
            [const.S_VALUE, const.S_VALUE_SCALED, const.S_VALUE_REF, const.S_VALUE_SCALED_REF]].replace(0,
                                                                                                        median_val_scaled / 10)
    # else:
    # df_concat = df_concat[~((df_concat[str_value] == 0) | (df_concat[str_value_ref] == 0))]
    df_concat[const.S_LOG2RATIO] = np.log2(df_concat[const.S_VALUE_SCALED] / df_concat[const.S_VALUE_SCALED_REF])
    # df_concat = df_concat[~abs(df_concat[str_log2_ratio]) == np.inf]
    return df_concat


def suppress_pos(rr, target_attr='klys', exp_name='exp', att_list=None, impute_missing=True, name_prot=False,
                 name_chain=True, sim_time=1e8, custom_params=None):
    df_list = []
    if custom_params is None:
        cp = {}
    else:
        cp = custom_params
    df_ref = get_sim_res(rr, last_frame_only=True, sim_time=sim_time, custom_params=cp, exp_name=exp_name)
    df_ref['cond_ref'] = 'S: None'
    # df_list.append(df_ref)
    if not att_list:
        # by default try ever
        att_list = [a for a in dir(rr) if target_attr in a]
    for att in att_list:
        if type(custom_params) == dict:
            cp = {att: 0, **custom_params}
        else:
            cp = {att: 0}
        df_melt_last_tp = get_sim_res(rr, last_frame_only=True, sim_time=sim_time, custom_params=cp, exp_name=exp_name)
        if df_melt_last_tp is not None:
            df_melt_last_tp['cond'] = f'S: {"_".join(extract_aa_pos(att, name_prot=name_prot, name_chain=name_chain))}'
            df_list.append(df_melt_last_tp)
    df_rr = get_reaction_ratio_df(df_ref=df_ref, df_list=df_list, impute_missing=impute_missing)
    return df_rr


def get_xl_connectivity_count_df(df, min_val=0, val_col=const.S_VALUE):
    group_by_common = [const.S_TIME, const.S_IS_SIGNIFICANT]
    if const.S_COND in df.columns:
        group_by_common.append(const.S_COND)
    if const.S_EXP in df.columns:
        group_by_common.append(const.S_EXP)
    group_by_pos1 = [const.S_PROT1, const.S_CHAIN_ID1, const.S_POS1]
    group_by_pos2 = [const.S_PROT2, const.S_CHAIN_ID2, const.S_POS2]
    df = df.copy()
    df[const.S_IS_SIGNIFICANT] = False
    df.loc[df[val_col] >= min_val, const.S_IS_SIGNIFICANT] = True
    df_pos1 = df[(df[const.S_TYPE].str.fullmatch(const.S_REACT_XL))].groupby(group_by_common + group_by_pos1)[
        const.S_VAR].count().reset_index().rename(
        columns={const.S_PROT1: const.S_PROT, const.S_CHAIN_ID1: const.S_CHAIN_ID, const.S_POS1: const.S_POS})
    df_pos2 = df[(df[const.S_TYPE].str.fullmatch(const.S_REACT_XL))].groupby(group_by_common + group_by_pos2)[
        const.S_VAR].count().reset_index().rename(
        columns={const.S_PROT2: const.S_PROT, const.S_CHAIN_ID2: const.S_CHAIN_ID, const.S_POS2: const.S_POS})
    df_concat = pd.concat([df_pos1, df_pos2])
    group_by_common.extend([const.S_PROT, const.S_CHAIN_ID, const.S_POS])
    df_count = df_concat.groupby(group_by_common)[const.S_VAR].sum().reset_index().rename(
        columns={const.S_VAR: const.S_COUNT})
    return df_count


def get_rr_significance_columns(df, min_val_abs=0, min_inc_abs=0, min_log2ratio_abs=0):
    species = [const.S_REACT_MONO_SUM, const.S_REACT_XL]
    species = '|'.join(species)
    df = df[df[const.S_TYPE].str.fullmatch(species)]
    df = df.copy()

    df[const.S_SIGNI_MIN_VAL] = False
    df.loc[df[const.S_VALUE_SCALED].abs() >= min_val_abs, const.S_SIGNI_MIN_VAL] = True

    df[const.S_SIGNI_MIN_LOG2] = False
    df.loc[df[const.S_LOG2RATIO].abs() >= min_log2ratio_abs, const.S_SIGNI_MIN_LOG2] = True

    df[const.S_SIGNI_MIN_INC] = False
    df.loc[df[const.S_VALUE_INC_ABS] >= min_inc_abs, const.S_SIGNI_MIN_INC] = True

    df[const.S_NON_SIGNI_REASON] = const.S_SIGNI_MIN_VAL + ": " + df[const.S_SIGNI_MIN_VAL].astype(
        str) + ", " + const.S_SIGNI_MIN_INC + ": " + df[const.S_SIGNI_MIN_INC].astype(
        str) + ", " + const.S_SIGNI_MIN_LOG2 + ": " + df[const.S_SIGNI_MIN_LOG2].astype(str)

    mask_overall = (df[const.S_SIGNI_MIN_VAL]) & (df[const.S_SIGNI_MIN_INC]) & (df[const.S_SIGNI_MIN_LOG2])

    df[const.S_IS_SIGNIFICANT] = False
    df.loc[mask_overall, const.S_IS_SIGNIFICANT] = True
    return df


def get_corr_df(df, x, y, group_by):
    corr_vars = [x, y]
    df_corr = (
        df.groupby(group_by)[corr_vars].corr(
            method='pearson').reset_index().  # calculate corr matrix separated by group_by
        drop(columns=[x]).  # drop superfluous columns (i.e. x and y)
        drop_duplicates(subset=group_by).  # we have dups since it's a matrix
        rename(columns={y: 'corr'})  # rename y to corr
        # fillna(0)
    )
    # if 'level_1' in df_corr.columns:
    # df_corr = df_corr.drop(columns=['level_1'])
    df = pd.merge(df, df_corr, on=group_by)
    return df


def get_xl_fraction_df(df_list, cond_list):
    df_list_frac = []
    for dfl, cond in zip(df_list, cond_list):
        df = dfl.groupby([const.S_EXP, const.S_TYPE])[const.S_VALUE_SCALED].sum().reset_index()
        df = df.pivot(columns=[const.S_TYPE], values=const.S_VALUE_SCALED, index=const.S_EXP).reset_index()
        df[const.S_REACT_MONO_XL_SUM] = df[const.S_REACT_XL] + df[const.S_REACT_MONO_SUM]
        df[const.S_REACT_XL_FRACTION] = df[const.S_REACT_XL] / df[const.S_REACT_MONO_XL_SUM]
        # df = df.drop(columns=[c for c in df if c not in [const.S_EXP, str_species_xl_fraction]])
        df = df.melt(
            id_vars=[const.S_EXP])  # .drop(columns=[str_type]).rename(columns={const.S_VALUE: str_species_xl_fraction})
        df[const.S_COND] = cond
        df_list_frac.append(df)
    df_concat = pd.concat(df_list_frac).reset_index(drop=True)
    return df_concat


def get_sorted_cats_list(series):
    return list(series.unique().sort_values())


def add_newline_to_outlist(out_list):
    out_list_n = []
    for o in out_list:
        out_list_n.append(o + '\n')
    return out_list_n


def normalize_vals(val_list, target_lower=0.0, target_upper=1.0):
    val_list = np.array(abs(val_list))  # treat negative numbers the same as positive ones
    if len(val_list) == 1:
        return target_upper  # if there is only one value return the maximum range
    # first do min max normalization leading to a range between 0 and 1
    min_val = np.min(val_list)
    max_val = np.max(val_list)
    norm_list_min_max = (val_list - min_val) / (max_val - min_val)
    # next scale to our targets
    norm_list_scaled = norm_list_min_max * (target_upper - target_lower) + target_lower
    # and finally round to 2 decimals before returning
    return np.round(norm_list_scaled, decimals=2)


def get_continous_color_scale(df, val_col, different_color_scale_mono_xl=False):
    group_by_list = [val_col + "_sign"]
    df[val_col + "_sign"] = "=="
    df[const.S_BASE_COLOR] = 0
    if (df[val_col] < 0).any():
        df.loc[df[val_col] > 0, val_col + "_sign"] = "+"
        # df.loc[df[val_col] > 0, cf.str_base_color] = 0
        df.loc[df[val_col] < 0, val_col + "_sign"] = "-"
        df.loc[df[val_col] < 0, const.S_BASE_COLOR] = 240
    col_val_norm = val_col + '_norm'
    if different_color_scale_mono_xl:
        group_by_list.append(const.S_TYPE)
    df[col_val_norm] = df.groupby(group_by_list)[val_col].transform(normalize_vals,
                                                                    **{"target_lower": 0.2, "target_upper": 1.0})
    df[const.S_HSL_COLOR] = map_colors_hsl_saturation(df[col_val_norm], df[const.S_BASE_COLOR])
    return df


# not used atm
def get_bins(df, val_col, bins=3, different_color_scale_mono_xl=False):
    col_count = 'count'
    col_bin_range = 'range'
    col_bin_start = col_bin_range + '_start'
    col_bin_end = col_bin_range + '_end'
    col_bin_start_norm = col_bin_start + '_norm'

    df = df[df[const.S_TYPE].str.fullmatch(f'{const.S_REACT_MONO_SUM}|{const.S_REACT_XL}')]
    df[col_bin_range] = pd.cut(df[val_col], bins=bins)
    df[[col_bin_start, col_bin_end]] = df[col_bin_range].astype(str).str.replace('(\(|\])', '', regex=True).str.split(
        ',', expand=True)
    df[col_bin_start], df[col_bin_end] = df[col_bin_start].astype(float), df[col_bin_end].astype(float)
    df[col_bin_start_norm] = normalize_vals(df[col_bin_start])

    df[col_count] = df[col_bin_start_norm].map(df[col_bin_start_norm].value_counts())
    df[const.S_HSL_COLOR] = map_colors_hsl_saturation(df[col_bin_start_norm])
    return df


def map_colors_hsl_saturation(saturation_list, hue_list, l=0.5):
    color_list = []
    for sat, hue in zip(saturation_list, hue_list):
        color_list.append((hue, sat, l))
    return color_list


def filter_species_min_val(df, val_col, min_val=0.05):
    df = df[df[val_col].abs() >= min_val]
    df = df[df[const.S_TYPE].str.fullmatch(f'{const.S_REACT_MONO_SUM}|{const.S_REACT_XL}')]
    return df


def get_chimerax_xl_color_script(df, pdb_id, supp_pos=None, supp_chain=None):
    df[[const.S_POS1, const.S_POS2]] = df[[const.S_POS1, const.S_POS2]].fillna(0)
    df[[const.S_POS1, const.S_POS2]] = df[[const.S_POS1, const.S_POS2]].astype(int)
    out_list = ['from chimerax.core.commands import run',
                'run(session, "close")',
                'session.logger.status("Hello Hello, someone out there???")',
                f'run(session, "open {pdb_id}")',
                # 'run(session, "show #1:lys atom")', #show side chains for lysines
                'run(session, "color #1:lys hsl(100, 0.5, 0.5)")',  # lightgreen for lysines
                'run(session, "distance style decimalPlaces 1")']

    if supp_pos:
        out_list.append(
            f'run(session, "color #1/{supp_chain}:{supp_pos} hsl(40, 1.0, 0.5)")')  # gold for suppression position

    for row in df[df[const.S_TYPE] == const.S_REACT_MONO_SUM].itertuples(index=False):
        if supp_pos != getattr(row, const.S_POS1):
            out_list.append((
                f'run(session, "color  #1/{getattr(row, const.S_CHAIN_ID1)}:{getattr(row, const.S_POS1)} hsl{getattr(row, const.S_HSL_COLOR)}")'))
    for row in df[df[const.S_TYPE] == const.S_REACT_XL].itertuples(index=False):
        out_list.append((
            f'run(session, "distance #1/{getattr(row, const.S_CHAIN_ID1)}:{getattr(row, const.S_POS1)}@CA #1/{getattr(row, const.S_CHAIN_ID2)}:{getattr(row, const.S_POS2)}@CA color hsl{getattr(row, const.S_HSL_COLOR)}")'))
    out_list = add_newline_to_outlist(out_list)
    return out_list


def get_chimerax_supp_df(df_rr, cond, val_col, min_val_abs=0, min_inc_abs=0, min_log2ratio_abs=0):
    df_rr = df_rr[df_rr['cond'] == cond]
    df_rr = get_rr_significance_columns(df_rr, min_val_abs=min_val_abs, min_inc_abs=min_inc_abs,
                                               min_log2ratio_abs=min_log2ratio_abs)
    df_rr = df_rr[df_rr[const.S_IS_SIGNIFICANT]]
    df_rr = get_continous_color_scale(df_rr, val_col=val_col, different_color_scale_mono_xl=True)
    # df_rr = df_rr.dropna(subset=[str_pos1, str_pos2])
    df_rr[[const.S_POS1, const.S_POS2]] = df_rr[[const.S_POS1, const.S_POS2]].fillna(0)
    df_rr[[const.S_POS1, const.S_POS2]] = df_rr[[const.S_POS1, const.S_POS2]].astype(int)
    return df_rr
