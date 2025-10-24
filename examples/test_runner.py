# examples/test_runner.py

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from methods.runner import process_alloys


T_list = [298.15]
activity_list = [1e-4]

mu_ligand = {'NH3': -0.276037, 'Gly': -3.263014109, 'CN': 1.786800089}


metal_1 = 'Ni'
metal_list = ['Ni','Cu']
# metal_list = list(set(['Au', 'Cu', 'Ni', 'Co',  'Mg', 'Mn', 'Zn', 'Ag', 'Cd','Sr','Pt','Pd','Fe','Cr','Cd','Pd'])) #'Fe',
aqueous_only = {'NH3': 0, 'NO2': 0, 'Gly': 0, 'CN': 0}
no_CN = {'NH3': 0.02, 'NO2': 0, 'Gly': 0.1, 'CN': 0}
with_CN = {'NH3': 0.02, 'NO2': 0, 'Gly': 0.1, 'CN': 1e-4}
experiment_concentration = {'NH3':0., 'NO2':0, 'Gly': 0.1, 'CN':0}
ligand_concentration_list = [experiment_concentration, with_CN, no_CN]

data_dir = '../data'
outdir = '../figures'
# process_alloys(metal_list, ligand_concentration_list, T_list, activity_list, mu_ligand, data_dir, outdir, save_fig=True)
process_alloys(['Ni','Ti'], ligand_concentration_list, T_list, activity_list, mu_ligand, data_dir, outdir, save_fig=True)
process_alloys(['Au','Pd'], ligand_concentration_list, T_list, activity_list, mu_ligand, data_dir, outdir, save_fig=True)

