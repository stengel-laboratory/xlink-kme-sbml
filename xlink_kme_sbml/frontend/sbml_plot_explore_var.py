#!/usr/bin/env python
import pandas as pd
from library import sbml_plot_lib as lib, sbml_sim_helper as helper
import argparse


desc = """Kai Kammer - 2020-09. 
Script to plot the exploration of a variable.
"""


def main():
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', action="store", required=True,
                        help="csv file of the explored variable")
    parser.add_argument('-m', '--model', action="store", required=True,
                        help="SBML model file")
    parser.add_argument('-o', '--output', action="store", default='plots',
                        help="Folder to save the plots to")
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    rr = helper.load_model_path(args.model)

    plot_master = lib.PlotMasterVariableExplorer(df, rr, args.output)
    plot_master.plot_all()


if __name__ == "__main__":
    main()