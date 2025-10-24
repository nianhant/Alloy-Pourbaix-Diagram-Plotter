import sys
import os
from pathlib import Path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from matplotlib import pyplot as plt
import itertools
from pymatgen.core import Composition

from methods.data_loader import SpeciesDataLoader
from methods.thermodynamics import PourbaixData, PourbaixAnalyzer
from methods.grid import GridMaker
from methods.visualization import GridVisualizer, PlotAccessories
from methods.utils import parse_composition, format_comp_dict
from methods.plot_util import generate_legends
from methods.set_publication_style import set_publication_style
from methods.runner import plot_pourbaix

set_publication_style()

Ni_Cu = ('Ni','Cu')
Ni_Ti = ('Ni','Ti')
Au_Pd = ('Au','Pd')
metal_list = [Ni_Cu, Ni_Ti, Au_Pd]

mu_ligand = {'NH3': -0.276037, 'Gly': -3.263014109, 'CN': 1.786800089}
T = 298.15

activity_list = [1e-4, 1e-5, 1e-6]

exp = {'NH3':0.02, 'NO2':0, 'Gly': 0.05, 'CN':0}
with_CN =  {'NH3':0.02, 'NO2':0, 'Gly': 0.05, 'CN':1e-4}
ligand_concentration_list = [exp, with_CN,
                             {'NH3':0.02, 'NO2':0, 'Gly': 0.005, 'CN':0},
                             {'NH3':0.02, 'NO2':0, 'Gly': 0.005, 'CN':1e-4},
                             {'NH3':0.02, 'NO2':0, 'Gly': 0.1, 'CN':0},
                             {'NH3':0.02, 'NO2':0, 'Gly': 0.1, 'CN':1e-4}]

data_dir = '../data'
pH_exp_range=(11.5, 13.5)
V_exp_range=(-2, 2.3)
save_fig = True
outdir = f"../figures/pourbaix_diagrams"

for m in metal_list:
    metal_1, metal_2 = m
    outdir=f'/home/x-ntian/pourbaix_paper/Accelerated-Computational-Materials-Discovery-for-Electrochemical-Nutrient-Recovery/Figures/pourbaix_diagrams/'

    for activity in activity_list:
        for ligand_concentration in ligand_concentration_list:
            plot_pourbaix(metal_1, metal_2, mu_ligand, T, activity, ligand_concentration, 
                          data_dir, pH_exp_range, V_exp_range, save_fig, outdir)
    