import pandas as pd
from .pymatgen_api import get_alloy_ion_formation_energy, get_alloy_solid_formation_energy
import json
import os
class SpeciesDataLoader:
    def __init__(self, metal_1, metal_2, mu_ligand, T, data_dir="data/"):
        self.metal_1 = metal_1
        self.metal_2 = metal_2
        self.mu_ligand = mu_ligand
        self.T = T
        self.data_dir = data_dir
        self.metal_reference = {metal_1: metal_1, metal_2: metal_2}

        # Load all data
        self.metal_complex_df = self.load_metal_complex_data()
        self.metal_complex_df = self.metal_complex_df[self.metal_complex_df["ligand"] != "NO2[1-]"]
        self.metal_complex_df['del_G_eV_vary_T'] = self.metal_complex_df.apply(self.calculate_G, axis=1)

        # Since Pd and Pt's complexes have controversial stability constants
        # if 'Pd' in self.metal_reference:
        #     self.metal_complex_df.loc[self.metal_complex_df["species"] == 'Pd(CN)4[2+]', "del_G_eV"] = 6.467825128
        # if 'Pt' in self.metal_reference:
        #     self.metal_complex_df.loc[self.metal_complex_df["species"] == 'Pt(CN)4[2+]', "del_G_eV"] = 5.646346039
         # harrington's stability constant = 70
        
        self.metal_complex_df["species_label"] = self.metal_complex_df.apply(self.create_species, axis=1)
        self.metal_complex, self.species_label_dict = self.filter_metal_complex()
        self.solid_eng, self.ion_eng = self.load_solid_ion_energies()
        self.stable_species = self.load_stable_species()
        
    def calculate_G(self, row):
        n_metal = int(row['n_metal'])
        G_metal = float(row['G_metal (kJ/mol)'])
        n_ligand = int(row['n_complex'])
        G_ligand = float(row['G_ligand (kJ/mol)'])

        log_B = float(row['log_B'])

        R = 8.314
        ln_log = 2.303
        kJ_per_mol_to_eV_per_atom = 0.01036427

        reactant_G = n_metal*G_metal + n_ligand * G_ligand
        del_G_rxn = -log_B * R * ln_log * self.T/1000 
        product_G = (del_G_rxn + reactant_G)* kJ_per_mol_to_eV_per_atom
        return product_G

    def load_metal_complex_data(self):
        """Loads metal complex formation energy data."""
        return pd.read_json(os.path.join(self.data_dir, "metal_complex_del_G.json"))

    def create_species(self, row):
        """Formats the species label based on charge balance."""
        def extract_ion_number(expression):
            if "[" in expression and "]" in expression:
                ion = expression[expression.find("[") + 1 : expression.find("]")]
                return int(ion.replace("+", "").replace("-", "")) * (-1 if "-" in ion else 1)
            return 0

        def format_charge(charge):
            return f"{abs(charge)}{'+' if charge > 0 else '-'}" if charge else ""

        metal_charge = extract_ion_number(row["signed_metal_ion"]) * row["n_metal"]
        ligand_charge = extract_ion_number(row["ligand"]) * row["n_complex"]
        total_charge = metal_charge + ligand_charge
        total_charge_str = format_charge(total_charge)

        species = f"[{row['signed_metal_ion'].split('[')[0]}{row['n_metal']}({row['ligand'].split('[')[0]}){row['n_complex']}]{total_charge_str}"
        species = f"{row['signed_metal_ion'].split('[')[0]}{row['n_metal']}({row['ligand'].split('[')[0]}){row['n_complex']}[{total_charge_str}]"

        return species.replace("_1", "").replace("^1", "").replace("+1", "+").replace("-1", "-").replace("[]","")

    def filter_metal_complex(self):
        """Filters data for selected metals."""
        target_df = self.metal_complex_df[
            (self.metal_complex_df["metal"] == self.metal_1) | (self.metal_complex_df["metal"] == self.metal_2)
        ]
        return target_df.set_index("species")["del_G_eV_vary_T"].to_dict(), target_df.set_index("species")["species_label"].to_dict()


    def load_solid_ion_energies(self):
        """Loads solid and ion formation energies."""
        solid_file = os.path.join(self.data_dir, f"{self.metal_1}_{self.metal_2}_solid_formation_energy.json")
        ion_file = os.path.join(self.data_dir, f"{self.metal_1}_{self.metal_2}_ion_formation_energy.json")

        solid_eng, ion_eng = {}, {}
        if os.path.exists(solid_file) and os.path.exists(ion_file):
            with open(solid_file, "r") as f:
                solid_eng = json.load(f)
            with open(ion_file, "r") as f:
                ion_eng = json.load(f)
        else:
            ion_eng = get_alloy_ion_formation_energy(self.data_dir, self.metal_1, self.metal_2)
            solid_eng = get_alloy_solid_formation_energy(self.data_dir, self.metal_1, self.metal_2)

        
        return solid_eng, ion_eng

    def load_stable_species(self):
        """Loads stable alloy species."""
        with open(os.path.join(self.data_dir, "metal_stable_species.json"), "r") as f:
            return json.load(f)