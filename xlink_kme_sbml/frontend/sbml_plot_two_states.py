#!/usr/bin/env python
import pandas as pd
from library import sbml_plot_lib as lib, sbml_constants as const, sbml_sim_helper as helper
import argparse


desc = """Kai Kammer - 2020-08. 
Script to plot quantitative crosslink 2-state models. Accepts two csv files as input.\n
These files need to contain only one frame from a simulation done with Tellurium.\n
Note that the first state also acts as the reference for log2ratio calculation.\n
"""


def main():
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s1', '--state_1', action="store", required=True,
                        help="csv file of state 1")
    parser.add_argument('-n1', '--name_1', action="store", default='state_1',
                        help="Name of the first state")
    parser.add_argument('-s2', '--state_2', action="store", required=True,
                        help="csv file of state 2")
    parser.add_argument('-n2', '--name_2', action="store", default='state_2',
                        help="Name of the second state")
    parser.add_argument('-m', '--map', action="store", required=True,
                        help=f"csv file mapping simulation name to UniProtID. Columns: {const.S_SIM_NAME}, {const.S_EXP_NAME}")
    parser.add_argument('-x', '--xtract', action="store",
                        help="Optional. xTract(-like) output to compare to the simulation")
    parser.add_argument('-d', '--dist', action="store",
                        help="Optional. File containing distances for uxIDs")
    parser.add_argument('-o', '--output', action="store", default='plots',
                        help="Folder to save the plots to")
    args = parser.parse_args()

    df_s1 = pd.read_csv(args.state_1)
    df_s2 = pd.read_csv(args.state_2)
    df_map = pd.read_csv(args.map)
    assert const.S_SIM_NAME in df_map and const.S_EXP_NAME in df_map, print(
        f"The map file requires the two columns: {const.S_EXP_NAME} and {const.S_SIM_NAME}")

    mapper_dict = helper.get_mapper_dict(df_map)

    df_xtract = None
    df_dist = None
    if args.xtract:
        if ".xls" in args.xtract:
            df_xtract = pd.read_csv(args.xtract, delimiter='\t').rename(columns={'type': 'link_type'})
            # remove imputed links
            df_xtract = df_xtract[df_xtract['sign'] == '==']
        else:
            df_xtract = pd.read_csv(args.xtract)
    if args.dist:
        df_dist = pd.read_csv(args.dist)
    plot_master = lib.PlotMasterTwoStates(df_s1, df_s2, args.name_1, args.name_2, mapper_dict, df_xtract, df_dist,
                                      args.output)
    plot_master.plot_all()


if __name__ == "__main__":
    main()
