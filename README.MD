# Kinetic Simulation Framework for Crosslinks
## Build & run kinetic crosslink simulations
### 2023-09 Kai-Michael Kammer
Kinetic simulation of the crosslink formation. Using SBML (https://sbml.org) to build the model and Tellurium for running the simulation (https://tellurium.readthedocs.io/en/latest/).

## Installation
This package requires Tellurium and simplesbml installed. Also required is the prot_tools package from https://github.com/stengel-laboratory/prot-tools  
The easiest way to install is to first follow the instructions for the prot_tools package. 
Having done so, the missing packages can be installed in the exisiting conda environment with

    conda activate prot_tools
    conda install scipy matplotlib python-libsbml jinja2 plotly
    pip install tellurium simplesbml

Note that at least on Fedora 35 Tellurium needs ncurses-libs and ncurses-compat-libs installed:

    sudo dnf install ncurses-libs ncurses-compat-libs

This package is missing on Fedora >= 37, so you'll have to symlink the following files

    sudo ln -s /usr/lib/libtinfo.so.6 /usr/lib/libtinfo.so.5
    sudo ln -s /usr/lib/libform.so.6 /usr/lib/libform.so.5
    sudo ln -s /usr/lib/libncurses.so.6 l/usr/lib/ibncurses.so.5
    sudo ln -s /usr/lib64/libtinfo.so.6 /usr/lib64/libtinfo.so.5
    sudo ln -s /usr/lib64/libform.so.6 /usr/lib64/libform.so.5
    sudo ln -s /usr/lib64/libncurses.so.6 l/usr64/lib/ibncurses.so.5

## Usage


