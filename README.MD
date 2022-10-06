# XLINK KME SBML
Kinetic simulation of the crosslink formation. Using SBML to build the model and Tellurium for the simulation.

## Installation
This package requires Tellurium, tesbml and simplesbml installed. Also required is the prot_tools package from 
The easiest way to install is to first follow the instructions for the prot_tools package. 
Having done so, the missing packages can be installed with
> conda activate prot_tools
> pip install tellurium simplesbml tesbml

Note that at least on Fedora 35 Tellurium needs ncurses-libs and ncurses-compat-libs installed:
> sudo dnf install ncurses-libs ncurses-compat-libs




