from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen.core import Composition
from .species import Species
from dataclasses import dataclass
import itertools
import numpy as np


@dataclass
class SpeciesGridResult:
    """Compact representation of the stable species at each grid point."""
    min_indices: np.ndarray
    species_list: list

    @staticmethod
    def format_species_tuple(species_tuple):
        formatted_list = [
            f"{species.formula}_{species.phase}_{species.alloy}"
            for species in species_tuple
        ]
        return tuple(sorted(formatted_list))

    @property
    def shape(self):
        return self.min_indices.shape

    def to_object_grid(self):
        species_grid = np.empty(self.shape, dtype=object)
        flat_grid = species_grid.ravel()
        for grid_idx, species_idx in enumerate(self.min_indices.ravel()):
            flat_grid[grid_idx] = self.species_list[species_idx]
        return species_grid

class PourbaixData:

    def __init__(self, species_data, activity, ligand_concentration, reference_composition):
        self.mu_ligand = species_data.mu_ligand
        self.metal_reference = species_data.metal_reference
        self.stable_species_names = species_data.stable_species
        self.ligand_concentration = ligand_concentration
        self.reference_composition = reference_composition
        self.activity = activity
        self.all_species_list = []
        for species, energy in species_data.solid_eng.items():
            solid_species = Species(species, 'solid', energy)
            self.all_species_list.append(solid_species)

        for species, energy in species_data.ion_eng.items():
            self.all_species_list.append(Species(species, 'ion', energy, activity))
            
        self.stable_species_list = self.filter_stable_species()
        
        for species, energy in species_data.metal_complex.items():
            self.stable_species_list.append(Species(species, 'complex', energy, activity = activity, ligand_concentration = ligand_concentration))
            
            self.all_species_list.append(Species(species, 'complex', energy, activity = activity, ligand_concentration = ligand_concentration))

        self.reaction_coefficients = {}
        
        self.compute_reaction_coefficients()


    def filter_stable_species(self):
        """Filters species that are in the stable species or are alloy."""
        stable_species_list = []
        for species in self.all_species_list:              
            if species.name in self.stable_species_names or all(metal in species.composition for metal in self.metal_reference):
                stable_species_list.append(species)
        return stable_species_list

    def generate_species_combinations(self):
        """Generates valid species combinations."""
        all_combinations = [list(itertools.combinations(self.stable_species_list, r)) for r in range(1, len(self.metal_reference) + 1)]
        return list(itertools.chain.from_iterable(all_combinations))
    
    def contains_required_metals(self, combo):
        required_metals = set(self.metal_reference.keys())  
        metals_in_combo = set()
        for species in combo:
            for metal in required_metals:
                if metal in species.composition:
                    metals_in_combo.add(metal)
        return metals_in_combo >= required_metals  

    # The following methods are based on pymatgen.analysis.pourbaix_diagram
    def compute_reaction_coefficients(self):
        """Calculates reaction coefficients for each species combination."""
        species_combinations = self.generate_species_combinations()
        for combo in species_combinations:
            if self.contains_required_metals(combo):
                coefficients = self.compute_reaction_coefficients_for_combo(combo)
                if coefficients:
                    self.reaction_coefficients[combo] = coefficients

                
    def compute_reaction_coefficients_for_combo(self, combo_tuple, coeff_threshold=1e-4):
        """Calculates reaction coefficients using pymatgen's reaction calculator."""
        dummy_oh = [Composition("H"), Composition("O")]
        comp_dict = {metal: 1 / len(self.metal_reference) for metal in self.metal_reference}
        prod_comp = self.reference_composition #Composition(comp_dict)
        
        composition_list = []
        for species in combo_tuple:
            metal_composition_dict = {key: species.composition[key] for key in species.composition if key in self.metal_reference}
            composition_list.append(Composition(metal_composition_dict))
        try:
            rxn = Reaction(composition_list + dummy_oh, [prod_comp])
            
            react_coeffs = [-coeff for coeff in rxn.coeffs[:len(combo_tuple)]]
            all_coeffs = [*react_coeffs, rxn.get_coeff(prod_comp)]
            if all(coeff > coeff_threshold for coeff in all_coeffs):
                return react_coeffs
            return None
        except:
            return None


class PourbaixCalculator:
    def __init__(self, species_data, grid_maker, T=298.15):
        """
        Initializes the Pourbaix Calculator.
        :param species_data: Instance of SpeciesDataLoader containing species and energy data.
        :param grid_maker: Instance of GridMaker containing pH, V, and ligand grids.
        """
        self.species_data = species_data
        self.grid_maker = grid_maker
        self.kB = 8.6173e-5  # eV/K
        self.T = T #298.15 # K
        self.mu_H2O = -2.458  # Reference for water
        self.combo_chemical_potential_dict = {}
    
    def apply_ion_correction(self, species):
#         if species.phase == 'complex':
#             species.activity = self.spcies_data.ligand_concentration[species.name]
        return species.energy + self.kB*self.T*np.log(species.activity)
  
    def formulate_coefficients(self, reactant):
        """Calculates coefficients for computing chemical potential."""
        phase = reactant.phase
        mass_balance = reactant.mass_balance

        if phase in ['ion', 'complex']:
            mu_react = self.apply_ion_correction(reactant)
        else:
            mu_react = reactant.energy

        
        coeff = self.kB * self.T * np.log(10)

        eU_coeff = mass_balance['n_charge']
        pH_coeff = mass_balance['n_H'] * coeff

        p_ligand_coeff_dict = {}
        total_ligand_mu = 0
        for ligand, n in mass_balance['n_L'].items():
            delta_mu_ligand = self.species_data.mu_ligand.get(ligand, 0)
            total_ligand_mu += n * delta_mu_ligand
            p_ligand_coeff_dict[ligand] = n * coeff

        constants = mu_react - mass_balance['n_H2O'] * self.mu_H2O - total_ligand_mu
        return eU_coeff, pH_coeff, p_ligand_coeff_dict, constants

    def compute_chemical_potential(self, eU_coeff, pH_coeff, p_ligand_coeff_dict, constants):
        """Computes the chemical potential grid for a given species."""
        pH_profile = pH_coeff * self.grid_maker.pH_values

        for ligand, coeff in p_ligand_coeff_dict.items():
            p_ligand_grid = self.grid_maker.ligand_grid_dict.get(ligand)
            pH_profile = pH_profile + coeff * p_ligand_grid

        return eU_coeff * self.grid_maker.V_values[:, np.newaxis] + pH_profile[np.newaxis, :] + constants

    def compute_combo_terms(self, combo_tuple, react_coeffs):
        """Reduces a species combination to separable V and pH terms."""
        V_coeff = 0.0
        pH_profile = np.zeros(self.grid_maker.grid_size, dtype=np.float64)
        constants = 0.0

        for i, species in enumerate(combo_tuple):
            react_coeff = react_coeffs[i]
            eU_coeff, pH_coeff, p_ligand_coeff_dict, species_constants = self.formulate_coefficients(species)
            V_coeff += react_coeff * eU_coeff
            pH_profile += react_coeff * pH_coeff * self.grid_maker.pH_values
            constants += react_coeff * species_constants

            for ligand, ligand_coeff in p_ligand_coeff_dict.items():
                p_ligand_grid = self.grid_maker.ligand_grid_dict.get(ligand)
                pH_profile += react_coeff * ligand_coeff * p_ligand_grid

        return V_coeff, pH_profile, constants

    def compute_combo_chemical_potential(self, combo_tuple, react_coeffs, row_slice=None):
        """Computes one combination's chemical potential grid, optionally for a row chunk."""
        V_coeff, pH_profile, constants = self.compute_combo_terms(combo_tuple, react_coeffs)
        V_values = self.grid_maker.V_values if row_slice is None else self.grid_maker.V_values[row_slice]
        return V_coeff * V_values[:, np.newaxis] + pH_profile[np.newaxis, :] + constants

    def compute_all_chemical_potentials(self):
        """Computes chemical potentials for all species combinations."""
        
        combo_react_coeffs_dict = self.species_data.reaction_coefficients
        for combo_tuple, react_coeffs in combo_react_coeffs_dict.items():
            self.combo_chemical_potential_dict[combo_tuple] = self.compute_combo_chemical_potential(
                combo_tuple, react_coeffs
            )

    def find_min_energy_species(self, chunk_rows=256):
        """Finds the most stable species at each grid point."""
        combo_react_coeffs_dict = self.species_data.reaction_coefficients
        species_list = list(combo_react_coeffs_dict.keys())
        if not species_list:
            raise ValueError("No valid species combinations were generated.")

        grid_shape = (self.grid_maker.grid_size, self.grid_maker.grid_size)
        min_energy = np.full(grid_shape, np.inf, dtype=np.float64)
        min_indices = np.zeros(grid_shape, dtype=np.int32)
        row_chunk = max(1, min(chunk_rows, self.grid_maker.grid_size))

        combo_terms = [
            self.compute_combo_terms(combo_tuple, combo_react_coeffs_dict[combo_tuple])
            for combo_tuple in species_list
        ]

        for combo_idx, (V_coeff, pH_profile, constants) in enumerate(combo_terms):
            for row_start in range(0, self.grid_maker.grid_size, row_chunk):
                row_stop = min(row_start + row_chunk, self.grid_maker.grid_size)
                row_slice = slice(row_start, row_stop)
                energy = (
                    V_coeff * self.grid_maker.V_values[row_slice, np.newaxis]
                    + pH_profile[np.newaxis, :]
                    + constants
                )
                chunk_min = min_energy[row_slice]
                chunk_indices = min_indices[row_slice]
                mask = energy < chunk_min
                chunk_min[mask] = energy[mask]
                chunk_indices[mask] = combo_idx

        all_species_tuples_set = set()
        for species_idx in np.unique(min_indices):
            species_tuple = species_list[species_idx]
            all_species_tuples_set.add(SpeciesGridResult.format_species_tuple(species_tuple))

        return SpeciesGridResult(min_indices, species_list), all_species_tuples_set


class PourbaixAnalyzer:
    def __init__(self, species_data, grid_maker, T):
        self.pourbaix_calculator = PourbaixCalculator(species_data, grid_maker, T)

    def analyze_and_plot(self):
        species_grid, all_species_tuples_set = self.pourbaix_calculator.find_min_energy_species()
        return species_grid, all_species_tuples_set 
