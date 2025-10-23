# examples/test_runner.py

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from methods.runner import process_alloys

CN_list = [0]
# CN_list = [0.0001]

T_list = [298.15]
activity_list = [1e-4]

mu_ligand = {'NH3': -0.276037, 'Gly': -3.263014109, 'CN': 1.786800089}


metal_1 = 'Ni'
metal_list = ['Ni','Co']
# metal_list = ['Cu']
# metal_list = list(set(['Au', 'Cu', 'Ni', 'Co',  'Mg', 'Mn', 'Zn', 'Ag', 'Cd','Sr','Pt','Pd','Fe','Cr','Cd','Pd'])) #'Fe',
data_dir = '/storage/home/hhive1/ntian30/data/ocp/pourbaix/data'
process_alloys(metal_list, CN_list, T_list, activity_list, mu_ligand, data_dir, save_fig=True)

# metal_1 = 'Pd'
# metal_list = ['Au']
# process_alloys(metal_1, metal_list, CN_list, T_list, activity_list, mu_ligand, save_fig=True)
