# examples/test_runner.py

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from methods.runner import process_alloys

CN_list = [0]
T_list = [298.15]
activity_list = [1e-4]

mu_ligand = {'NH3': -0.276037, 'Gly': -3.263014109, 'CN': 1.786800089}

# metal_1 = 'Pd'
# metal_list = ['Au']  # Modify as needed
metal_1 = 'Ti'
metal_list = ['Ni']
metal_list = ['Cu']

process_alloys(metal_1, metal_list, CN_list, T_list, activity_list, mu_ligand, save_fig=True)
