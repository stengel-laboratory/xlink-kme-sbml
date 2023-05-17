#!/usr/bin/env python
# %%
# fitting imports
from __future__ import print_function, absolute_import

import pandas as pd
import tellurium as te
import xlink_kme_sbml.library.sbml_constants as const
import numpy as np
from timeit import default_timer as timer
import itertools as it

# fitting imports

import csv
from scipy.optimize import differential_evolution
import random


# from multiprocessing.pool import Pool


# %%
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


# result_list_sim = []
#
#
# def log_result_sim(result):
#     # This is called whenever foo_pool(i) returns a result.
#     # result_list is modified only by the main process, not the pool workers.
#     result_list_sim.append(result)
#
#
# result_list_load = []
#
#
# def log_result_load(result):
#     # This is called whenever foo_pool(i) returns a result.
#     # result_list is modified only by the main process, not the pool workers.
#     result_list_load.append(result)


# does not work with SWIG objects which Tellurium is using
# def load_and_simulate_async(exp, exp_type, variable, var_range, mapper_dict=None):
# with Pool(len(var_range)) as pool:
#     poolObjects = []
#     for value in var_range:
#         poolObjects.append(pool.apply_async(load_model, args=(exp, exp_type), callback=log_result_load))
#     pool.close()
#     pool.join()
#     for f in poolObjects:
#         f.get()
# print(result_list_load)
# rr = load_model(exp, exp_type)
# result_list_load = [ for rr in ]
# with Pool(len(var_range)) as pool:
#     poolObjects = []
#     for i, value in enumerate(var_range):
#         # rr = load_model(exp, exp_type)
#         poolObjects.append(pool.apply_async(_explore_variable_runner, args=(rr, variable, value, exp, mapper_dict),
#                          callback=log_result_sim))
#     pool.close()
#     pool.join()
#     for f in poolObjects:
#         f.get()
# return result_list_sim





def __run_sim_runner(rr, var_name, val_name, exp_name, mapper_dict, min_simulation_time, sim_vars=None,
                     scale_factor=None, scale_to_one=False):
    rr.resetToOrigin()
    if sim_vars:
        for var in sim_vars.keys():
            # if hasattr(rr, var):
            setattr(rr, var, sim_vars[var])
    # crosslinker_start = rr.Crosslinker
    lys_list = [a for a in dir(rr) if const.S_LYS in a]
    lys_initial_conc = getattr(rr, lys_list[0])
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
    lys_no = len([lys for lys in df_res.columns if const.S_LYS in lys])
    df_res = prepare_df(df_res, exp_name, mapper_dict)
    if scale_factor:
        df_res[const.S_VALUE] *= scale_factor
    elif scale_to_one:
        df_res[const.S_VALUE] /= lys_initial_conc
    df_res[const.S_RATIO_CROSSLINKER_LYSINE] = crosslinker_initial_conc/(lys_initial_conc*lys_no)
    rr.resetToOrigin()
    return df_res

def _explore_variable_runner(rr, variable, value, exp_name, mapper_dict, min_simulation_time, custom_vars=None,
                             scale_factor=None, scale_to_one=False):
    sim_vars = {variable: value}
    if custom_vars:
        sim_vars = {**custom_vars, variable: value}
    df_res = __run_sim_runner(rr, var_name=variable, val_name=value, exp_name=exp_name, mapper_dict=mapper_dict,
                              min_simulation_time=min_simulation_time, sim_vars=sim_vars, scale_factor=scale_factor,
                              scale_to_one=scale_to_one)
    if df_res is not None:
        df_res[variable] = value
    return df_res

def _explore_variable_runner_mult(rr, var_val_dict, exp_name, mapper_dict, min_simulation_time, sim_vars=None,
                                  scale_factor=None, scale_to_one=False):
    var_names = list(var_val_dict.keys())
    vals = list(var_val_dict.values())
    df_res = __run_sim_runner(rr, var_name=var_names, val_name=vals, exp_name=exp_name, mapper_dict=mapper_dict,
                              min_simulation_time=min_simulation_time, sim_vars=sim_vars, scale_factor=scale_factor,
                              scale_to_one=scale_to_one)
    if df_res is not None:
        for variable, value in var_val_dict.items():
            df_res[variable] = value
    return df_res


def explore_variable(rr, variable, var_range, exp_name='exp', mapper_dict=None, custom_vars=None,
                     min_simulation_time=5000000, scale_factor=None, scale_to_one=False):
    df_list = []
    for value in var_range:
        df_res = _explore_variable_runner(rr=rr, variable=variable, value=value, exp_name=exp_name,
                                          mapper_dict=mapper_dict, min_simulation_time=min_simulation_time,
                                          custom_vars=custom_vars, scale_factor=scale_factor, scale_to_one=scale_to_one)
        if df_res is not None:
            df_list.append(df_res)
    df_final = pd.concat(df_list)
    # df_melt = pd.melt(df_final.drop(columns=['time']), id_vars=variable)
    return df_final.reset_index(drop=True)


def explore_variable_multi(rr, variable_range_dict, exp_name='exp', mapper_dict=None, custom_vars=None,
                           min_simulation_time=5000000, scale_factor=None, scale_to_one=False, var_operation_dict=None):
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
        rr.resetToOrigin() # otherwise getting the original value for value operations is not reliable
        for n, var in enumerate(variable_range_dict.keys()):
            if type(var) is tuple:
                var, vars = var
                for v in vars:
                    sim_vars[v] = _apply_var_operation(rr=rr, var=var, var_operation_dict=var_operation_dict, new_value=val_tuple[n], var_name_long=v)
            else:
                sim_vars[var] = _apply_var_operation(rr=rr, var=var, var_operation_dict=var_operation_dict, new_value=val_tuple[n])
            var_val_dict[var] = val_tuple[n]  # contains the 'short' var names, i.e. LYS instead of [LYS_1, LYS_2, ...]
        df_res = _explore_variable_runner_mult(rr=rr, var_val_dict=var_val_dict, exp_name=exp_name,
                                               mapper_dict=mapper_dict, min_simulation_time=min_simulation_time,
                                               sim_vars={**custom_vars, **sim_vars}, scale_factor=scale_factor,
                                               scale_to_one=scale_to_one,
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
    df_const = pd.DataFrame.from_dict(klys_dict, orient='index', columns=[const.D_CONCENTRATION]).rename_axis(
        const.S_UID_SHORT).reset_index()
    df = pd.merge(df, df_const, on=const.S_UID_SHORT, how='outer')
    df[const.D_CONCENTRATION] = df[const.D_CONCENTRATION].fillna(0)
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


# if __name__ == '__main__':
#     load_and_simulate_async('c3', 'asa', 'kh', [1e-5, 1e-4, 1e-3, 1e-2])
# %%
# load_and_simulate('c3', 'asa')
# load_and_simulate('c3b', 'asa')
# df_final_frame = write_results(exp, exp_type, results)


class ParameterEstimation(object):
    """Parameter Estimation"""

    def __init__(self, stochastic_simulation_model, bounds, data=None):
        if (data is not None):
            self.data = data

        self.model = stochastic_simulation_model
        self.bounds = bounds

    def setDataFromFile(self, FILENAME, delimiter=",", headers=True):
        """Allows the user to set the data from a File
        This data is to be compared with the simulated data in the process of parameter estimation
        
        Args:
            FILENAME: A Complete/relative readable Filename with proper permissions
            delimiter: An Optional variable with comma (",") as default value. 
                A delimiter with which the File is delimited by.
                It can be Comma (",") , Tab ("\t") or anyother thing
            headers: Another optional variable, with Boolean True as default value
                If headers are not available in the File, it can be set to False
        Returns:
            None but sets class Variable data with the data provided
        
        .. sectionauthor:: Shaik Asifullah <s.asifullah7@gmail.com>
        
        
        """
        with open(FILENAME, 'r') as dest_f:
            data_iter = csv.reader(dest_f,
                                   delimiter=",",
                                   quotechar='"')
            self.data = [data for data in data_iter]
        if (headers):
            self.data = self.data[1:]

        self.data = np.asarray(self.data, dtype=float)

    def setCustomAttr(self, attr_dict):
        self.attr_dict = attr_dict

    def run(self, func=None, params=None):
        """Allows the user to set the data from a File
        This data is to be compared with the simulated data in the process of parameter estimation
        
        Args:
            func: An Optional Variable with default value (None) which by default run differential evolution
                which is from scipy function. Users can provide reference to their defined function as argument.
            
        Returns:
            The Value of the parameter(s) which are estimated by the function provided.
        
        .. sectionauthor:: Shaik Asifullah <s.asifullah7@gmail.com>
        
        
        """
        if not params:
            params = {}
        self._parameter_names = list(self.bounds.keys())
        self._parameter_bounds = list(self.bounds.values())
        self._model_roadrunner = te.loada(self.model.model)
        x_data = self.data[:, 0]

        x_data = self.data[:, 0]
        y_data = self.data[:, 1:]
        arguments = (x_data, y_data)

        if func is None:
            result = differential_evolution(self._SSE, self._parameter_bounds, args=arguments, **params)
            return (result.x)
        else:
            result = func(self._SSE, self._parameter_bounds, args=arguments)
            return (result.x)

    def _set_theta_values(self, theta):
        """ Sets the Theta Value in the range of bounds provided to the Function.
            Not intended to be called by user.
            
        Args:
            theta: The Theta Value that is set for the function defined/provided
            
        Returns:
            None but it sets the parameter(s) to the stochastic model provided
        
        .. sectionauthor:: Shaik Asifullah <s.asifullah7@gmail.com>
        
        
        """
        for theta_i, each_theta in enumerate(self._parameter_names):
            # print(f"Set theta: {each_theta}: {theta[theta_i]}")
            setattr(self._model_roadrunner, each_theta, theta[theta_i])

    def _SSE(self, parameters, *data):
        """ Runs a simuation of SumOfSquares that get parameters and data and compute the metric.
            Not intended to be called by user.
            
        Args:
            parameters: The tuple of theta values  whose output is compared against the data provided
            data: The data provided by the user through FileName or manually 
                  which is used to compare against the simulations
            
        Returns:
            Sum of Squared Error
        
        .. sectionauthor:: Shaik Asifullah <s.asifullah7@gmail.com>
        
            
        """
        theta = parameters

        x, y = data
        sample_x, sample_y = data
        self._set_theta_values(theta)

        random.seed()
        # it is now safe to use random.randint
        # self._model.setSeed(random.randint(1000, 99999))

        self._model_roadrunner.integrator.variable_step_size = self.model.variable_step_size
        self._model_roadrunner.reset()
        if hasattr(self, 'attr_dict'):
            for att_key, att_val in self.attr_dict.items():
                setattr(self._model_roadrunner, att_key, att_val)

        simulated_data = self._model_roadrunner.simulate(self.model.from_time, self.model.to_time,
                                                         self.model.step_points)
        simulated_data = np.array(simulated_data)
        simulated_x = simulated_data[:, 0]
        simulated_y = simulated_data[:, 1:]

        SEARCH_BEGIN_INDEX = 0
        SSE_RESULT = 0

        for simulated_i in range(len(simulated_y)):
            y_i = simulated_y[simulated_i]
            # yhat_i = sample_y[simulated_i]

            x_i = simulated_x[simulated_i]
            for search_i in range(SEARCH_BEGIN_INDEX + 1, len(sample_x)):
                if (sample_x[search_i - 1] <= x_i < sample_x[search_i]):
                    yhat_i = sample_y[search_i - 1]
                    break
                SEARCH_BEGIN_INDEX += 1

            partial_result = 0
            for sse_i in range(len(y_i)):
                partial_result += (float(y_i[sse_i]) - float(yhat_i[sse_i])) ** 2
            SSE_RESULT += partial_result

        return SSE_RESULT ** 0.5


def sim_for_fake_t(rr, time):
    results = rr.simulate(0, time, points=2)
    df_tmp = pd.DataFrame(results, columns=results.colnames)
    df_tmp = df_tmp.append(df_tmp.tail(1))
    df_tmp.iloc[1, 0] = time * 0.5
    return df_tmp.reset_index(drop=True)
