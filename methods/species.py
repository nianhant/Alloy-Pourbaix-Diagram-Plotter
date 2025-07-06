import numpy as np
import re
from ase.symbols import string2symbols
import itertools
import pandas as pd
import sys
import os
import json

class Species:
    def __init__(self, name, phase, energy, activity=None, ligand_concentration = None):
        self.name = name
        self.formula = self.name.split('(aq)')[0]
        self.phase = phase
        self.energy = energy
        self.mass_balance = None  # Will be computed later
        self.ligand = None
        self.activity = activity
        self.alloy = False
        if self.phase == 'complex':
            self.composition, self.ligand = self.parse_complex_composition()
        else:
            self.composition = self.parse_composition()
        self.mass_balance = self.compute_mass_balance()
        
    
    def parse_composition(self):
        """Parses the chemical formula to determine element composition."""
        
        formula  = self.formula
        composition = {'charge': 0, 'H': 0, 'O': 0}
        match = re.search(r'\[(\d*)([+-])\]', formula)
        if match:
            charge = match.group(1)
            sign = match.group(2)
            composition['charge'] = int(sign + charge)

        element_only_form = formula.split('[')[0]
        symbol_list = string2symbols(element_only_form)
        n_metal = 0
        for symbol in symbol_list:
            if symbol not in composition:
                n_metal += 1
            if symbol in composition:
                composition[symbol] += 1
            else:
                composition[symbol] = 1
        if n_metal > 1:
            self.alloy = True
        return composition
    
    def parse_complex_composition(self):
        formula  = self.formula
        composition = {'H': 0,
                       'O': 0,
                       'charge': 0}
        ############## count n of electrons and protons ##############
        match = re.search(r'\[(\d*)([+-])\]', formula)
        if match:
            charge = match.group(1)
            sign = match.group(2)
            composition['charge'] = int(sign + charge)
        ############## count n of ligands ##############
        pattern = r'([A-Za-z]+)?\(?([A-Za-z0-9]+)\)?(\d*)'
        match = re.search(pattern, formula)
        if match:
            metal = match.group(1)      # Metal symbol (e.g., Ni)
            ligand = match.group(2)     # Ligand (e.g., NH3)
            ligand_count = match.group(3)  # Ligand count (e.g., 4, could be empty)
            # If the ligand count is empty, it means it's 1
            ligand_count = int(ligand_count) if ligand_count else 1
        composition[ligand] = ligand_count
        ############## count n of metal ##############
        element_only_form = formula.split('(')[0]
        symbol_list = string2symbols(element_only_form)
        for symbol in symbol_list:
            if symbol in composition:
                composition[symbol] += 1
            else:
                composition[symbol] = 1
        return composition, ligand

    def compute_mass_balance(self):
        """Computes mass balance properties for reactions."""
        
        composition = self.composition
        mass_balance = {'n_H2O': 0, 'n_charge': 0, 'n_L': {}, 'n_H': 0}
        n_L = {}

        if self.phase == 'complex':
            n_L[self.ligand] = composition[self.ligand]  

        mass_balance['n_H2O'] = composition['O']
        mass_balance['n_H'] = composition['H'] - 2 * composition['O']
        mass_balance['n_charge'] = mass_balance['n_H'] - composition['charge']
        mass_balance.update({ 'n_L': n_L})
        
        return mass_balance
