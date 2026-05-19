import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba
import colorsys
import numpy as np
import re
import os
try:
    from ase.formula import Formula
except ModuleNotFoundError:
    class Formula:
        def __init__(self, formula):
            self.formula = formula

        def reduce(self):
            return (self, 1)

        def __format__(self, spec):
            return self.formula


LIGAND_ELEMENT_COUNTS = {
    'NH3': {'N': 1, 'H': 3},
    'CN': {'C': 1, 'N': 1},
    'Gly': {'C': 2, 'H': 4, 'N': 1, 'O': 2},
    'NO2': {'N': 1, 'O': 2},
}


class GridVisualizer:
    def __init__(self, grid_maker, species_data):
        """Takes a GridMaker instance and uses it for visualization."""
        self.grid_maker = grid_maker
        self.species_data = species_data
        self.metal_list = [key for key in  species_data.metal_reference.keys()]
    
    def add_H2_O2_lines(self, ax):
        PREFAC = 0.0591
        xlim = self.grid_maker.pH_range
        ylim = self.grid_maker.V_range
        h_line = np.transpose([[xlim[0], -xlim[0] * PREFAC], [xlim[1], -xlim[1] * PREFAC]])
        o_line = np.transpose([[xlim[0], -xlim[0] * PREFAC + 1.23], [xlim[1], -xlim[1] * PREFAC + 1.23]])
        
        lw = 1
        # Plot the hydrogen and oxygen lines on the axis 
        h_line_plot, = ax.plot(h_line[0], h_line[1], "b--", linewidth=lw, label= r'H$_2$O Reduction')
        o_line_plot, = ax.plot(o_line[0], o_line[1], "r--", linewidth=lw, label= r'H$_2$O Oxidation')


    def add_plot_accessories(self, ax, pH_exp_range=(11.5,13.5), V_exp_range=(-2, 2.3), ):
        PREFAC = 0.0591
        box_left = pH_exp_range[0]
        box_right = pH_exp_range[1]
        
        
        V_left_bottom = V_exp_range[0] - PREFAC * box_left
        V_right_bottom = V_exp_range[0] - PREFAC * box_right
        V_left_top = V_exp_range[1] - PREFAC * box_left
        V_right_top = V_exp_range[1] - PREFAC * box_right
        
        lw = 1
        style = 'k-'
        # Bottom side (V as a function of pH)
        bottom, = ax.plot([box_left, box_right], [V_left_bottom, V_right_bottom], style, lw=lw, 
        label = f'Exp condition\nV vs RHE={V_exp_range[0]}-{V_exp_range[1]}\npH={pH_exp_range[0]}-{pH_exp_range[1]}')  
        # Top side (V as a function of pH)
        top, = ax.plot([box_left, box_right], [V_left_top, V_right_top], style, lw=lw, label = f'{V_exp_range[1]}V vs RHE')     
        # Left side (fixed pH = box_left)
        left, = ax.plot([box_left, box_left], [V_left_bottom, V_left_top], style, lw=lw, label = f'pH={pH_exp_range[0]}')      
        # Right side (fixed pH = box_right)
        right, = ax.plot([box_right, box_right], [V_right_bottom, V_right_top], style, lw=lw,label = f'pH={pH_exp_range[1]}')     
        
        
    def format_file_name(self):
        ligand_concentration = self.species_data.ligand_concentration
        NH3 = ligand_concentration['NH3']
        Gly = ligand_concentration['Gly']
        CN = ligand_concentration['CN']
        activity = self.species_data.activity
        metal_list = "_".join(self.metal_list)
        
        reference_composition = self.species_data.reference_composition
        
        # Ensure directory exists
        # os.makedirs(outdir, exist_ok=True)
        
        # file_name = os.path.join(
        #     outdir,
        #     f"{metal_list}_alloy_{reference_composition}_NH3={NH3}M_Gly={Gly}M_CN={CN}M_activity={activity:.0e}M.png"
        # )
        file_name = f"{reference_composition}_NH3={NH3}M_Gly={Gly}M_CN={CN}M_activity={activity:.0e}M.png"
        
        return file_name

    def plot_species_distribution(self, species_grid, species_colors, ax=None, save_fig=True):
        """Plots the Pourbaix diagram using species distribution data."""
        converted_species_grid = [
            tuple(sorted(f"{species.formula}_{species.phase}_{species.alloy}" for species in species_tuple))
            for species_tuple in species_grid.flatten()
        ]

        color_values = np.array([species_colors[species_tuple] for species_tuple in converted_species_grid])
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))  # Create new figure if no axis is provided
        else:
            fig = ax.figure
        ax.set_xlim(self.grid_maker.pH_range)
        ax.set_ylim(self.grid_maker.V_range)
        ax.scatter(self.grid_maker.pH_grid.flatten(), self.grid_maker.V_grid.flatten(), 
                   c=color_values, s=1, alpha=0.75)
        
        ax.set_xlabel('pH')
        ax.set_ylabel(r'$E_{SHE}(V)$')
        self.add_H2_O2_lines(ax)
        self.add_plot_accessories(ax)
            
        return fig, ax



class PlotAccessories:
    def __init__(self,species_data): #all_species_tuples
        """Takes a GridMaker instance and uses it for visualization."""
        self.species_data = species_data

    @staticmethod
    def parse_formula_counts(formula):
        if formula in LIGAND_ELEMENT_COUNTS:
            return dict(LIGAND_ELEMENT_COUNTS[formula])

        counts = {}
        for element, number in re.findall(r'([A-Z][a-z]?)(\d*)', formula):
            counts[element] = counts.get(element, 0) + int(number or 1)
        return counts

    @staticmethod
    def get_formula_element_counts(formula):
        clean_formula = re.sub(r'\[[^\]]+\]', '', formula.replace('(aq)', ''))
        counts = {}

        for group, multiplier in re.findall(r'\(([A-Za-z0-9]+)\)(\d*)', clean_formula):
            multiplier = int(multiplier or 1)
            for element, count in PlotAccessories.parse_formula_counts(group).items():
                counts[element] = counts.get(element, 0) + count * multiplier

        clean_formula = re.sub(r'\([A-Za-z0-9]+\)\d*', '', clean_formula)
        for element, count in PlotAccessories.parse_formula_counts(clean_formula).items():
            counts[element] = counts.get(element, 0) + count

        return counts

    @staticmethod
    def split_species_key(species):
        parts = species.rsplit('_', 2)
        return parts[0], parts[1], parts[2] == 'True'

    def classify_combo(self, combo_tuple):
        phases = []
        has_oxide = False
        has_oxyhydroxide = False
        has_hydride = False
        has_ligand_complex = False
        is_alloy_combo = all(self.split_species_key(species)[2] for species in combo_tuple)

        for species in combo_tuple:
            formula, phase, _ = self.split_species_key(species)
            phases.append(phase)

            if 'complex' in phase:
                has_ligand_complex = True
                continue

            element_counts = self.get_formula_element_counts(formula)
            oxygen_count = element_counts.get('O', 0)
            hydrogen_count = element_counts.get('H', 0)

            if oxygen_count > 0 and hydrogen_count > 0:
                has_oxyhydroxide = True
            elif oxygen_count > 0:
                has_oxide = True
            elif hydrogen_count > 0:
                has_hydride = True

        if has_ligand_complex:
            return 'metal_ligand_complex'
        if all(phase != 'solid' for phase in phases):
            return 'aqueous_metal_ion'
        if has_oxyhydroxide:
            return 'metal_oxyhydroxide'
        if has_oxide:
            return 'metal_oxide'
        if has_hydride:
            return 'metal_hydride'
        if is_alloy_combo or all(phase == 'solid' for phase in phases):
            return 'metal'
        return 'aqueous_metal_ion'

    def color_intensity_for_combo(self, combo_tuple, category, used_intensities=None):
        total_metals = 0
        total_oxygen = 0
        total_other = 0

        for species in combo_tuple:
            formula, _, _ = self.split_species_key(species)
            element_counts = self.get_formula_element_counts(formula)
            for element, count in element_counts.items():
                if element in ['O', 'H', 'N', 'C']:
                    total_other += count
                else:
                    total_metals += count
            total_oxygen += element_counts.get('O', 0)

        if category == 'metal':
            intensity = 0.36
        elif category == 'metal_hydride':
            intensity = 0.42
        elif category == 'metal_oxide' and total_oxygen > 0:
            oxygen_per_metal = total_oxygen / max(total_metals, 1)
            intensity = 0.82 - 0.16 * oxygen_per_metal
        elif total_other > 0:
            metal_fraction = total_metals / (total_metals + total_other)
            intensity = 0.22 + 0.72 * metal_fraction
        else:
            intensity = 0.50

        if category == 'metal_oxyhydroxide':
            intensity += 0.08
        elif category == 'metal_ligand_complex':
            intensity += 0.10

        intensity = min(0.88, max(0.22, intensity))
        if used_intensities is None:
            return intensity

        if category == 'metal_oxide':
            used_intensities.append(intensity)
            return intensity

        while any(abs(intensity - used) < 0.055 for used in used_intensities):
            intensity += 0.065
            if intensity > 0.90:
                intensity = 0.25
        used_intensities.append(intensity)
        return intensity

    def get_color_for_category(self, combo_tuple, category, used_intensities=None):
        color_maps = {
            'metal': 'Greys',
            'metal_hydride': 'Greys',
            'metal_oxide': 'YlOrBr',
            'metal_oxyhydroxide': 'Oranges',
            'aqueous_metal_ion': 'PuBuGn',
            'metal_ligand_complex': 'RdPu',
        }
        color_map = plt.cm.get_cmap(color_maps.get(category, 'Greys'))
        intensity = self.color_intensity_for_combo(combo_tuple, category, used_intensities)
        return color_map(intensity)
        
    def count_total_phases(self, all_species_tuples):
        total_solid = 0
        
        for combo_tuple in all_species_tuples:
            if all('solid' in species for species in combo_tuple):
                total_solid += 1
        total_aq = len(all_species_tuples) - total_solid
        return total_solid, total_aq
    
    def get_color_for_label(self, all_species_tuples): 
        species_colors = {}
        used_intensities = {}
        
        for combo_tuple in all_species_tuples:  
            category = self.classify_combo(combo_tuple)
            species_colors[combo_tuple] = self.get_color_for_category(
                combo_tuple,
                category,
                used_intensities.setdefault(category, []),
            )
        return species_colors
    
    def format_formula(self, formula):
        if '[' in formula or 'aq' in formula or 'Gly' in formula:
            formatted_formula = re.sub(r"([A-Za-z\)\]])(\d+)", r"\1$_{\2}$", formula)
            formatted_formula = re.sub(r"\[([\d\+\-]+)\]", r"$^{\1}$", formatted_formula)

            formatted_formula = re.sub(r"\$_\{1\}\$", "", formatted_formula)
            formatted_formula = re.sub(r"\^\{1([+-])\}", r"^{\1}", formatted_formula)

        else:
            formula_obj = Formula(formula)
            reduced_formula = formula_obj.reduce()[0]
            formatted_formula = f'{reduced_formula:latex}'    
        
        return formatted_formula

    def format_species_combo(self, combo_tuple):
        species_list = list(combo_tuple)
        complex_label_dict = self.species_data.species_label_dict

        for i in range(len(species_list)):
            species = species_list[i]
            if 'complex' in species:
                true_charge_formula = complex_label_dict[species.split('_')[0]]
                formatted_species = self.format_formula(true_charge_formula)
            else:
                formatted_species = self.format_formula(species.split('_')[0])
            species_list[i] = formatted_species + ('(s)' if 'solid' in species else '(aq)')
            
        return '+'.join(species_list)
