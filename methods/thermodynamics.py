from pymatgen.analysis.reaction_calculator import Reaction
from pymatgen.core import Composition
from .species import Species
import itertools
import numpy as np

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
        chemical_potential = eU_coeff * self.grid_maker.V_grid + pH_coeff * self.grid_maker.pH_grid + constants

        for ligand, coeff in p_ligand_coeff_dict.items():
            p_ligand_grid = self.grid_maker.ligand_grid_dict.get(ligand)
            chemical_potential += coeff * p_ligand_grid

        return chemical_potential

    def compute_all_chemical_potentials(self):
        """Computes chemical potentials for all species combinations."""
        
        combo_react_coeffs_dict = self.species_data.reaction_coefficients
        for combo_tuple, react_coeffs in combo_react_coeffs_dict.items():
            self.combo_chemical_potential_dict[combo_tuple] = np.zeros(
                (self.grid_maker.grid_size, self.grid_maker.grid_size), dtype=np.float64
            )
            for i, species in enumerate(combo_tuple):
                eU_coeff, pH_coeff, p_ligand_coeff_dict, constants = self.formulate_coefficients(species)
                chemical_potential = self.compute_chemical_potential(eU_coeff, pH_coeff, p_ligand_coeff_dict, constants)
                self.combo_chemical_potential_dict[combo_tuple] += chemical_potential * react_coeffs[i]

    def find_min_energy_species(self):
        """Finds the most stable species at each grid point."""
        species_list = list(self.combo_chemical_potential_dict.keys())
        energy_grids = np.array([self.combo_chemical_potential_dict[species] for species in species_list])

        min_indices = np.argmin(energy_grids, axis=0)
        min_species_grid = np.empty(min_indices.shape, dtype=object)
        all_species_tuples_set = set()
        
        for i in range(min_species_grid.shape[0]):
            for j in range(min_species_grid.shape[1]):
                species_tuple = species_list[min_indices[i, j]]
                min_species_grid[i, j] = species_tuple
                formatted_list = [f"{species.formula}_{species.phase}_{species.alloy}" for species in species_tuple]
                formatted_tuple = tuple(sorted(formatted_list))
                all_species_tuples_set.add(formatted_tuple)
        return min_species_grid, all_species_tuples_set


class PourbaixAnalyzer:
    def __init__(self, species_data, grid_maker, T):
        self.pourbaix_calculator = PourbaixCalculator(species_data, grid_maker, T)

    def analyze_and_plot(self):
        self.pourbaix_calculator.compute_all_chemical_potentials()

        species_grid, all_species_tuples_set = self.pourbaix_calculator.find_min_energy_species()
        return species_grid, all_species_tuples_set 
