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

AREA_LABEL_FONT_SIZE = 20
REGION_IMAGE_MAX_PIXELS = 2400
REGION_COLOR_STRENGTH = 0.75
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

    def plot_species_distribution(
        self, species_grid, species_colors, ax=None, save_fig=True,
        label_stable_regions=True, region_style='image'
    ):
        """Plots the Pourbaix diagram using species distribution data."""
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))  # Create new figure if no axis is provided
        else:
            fig = ax.figure
        fig.patch.set_facecolor('white')
        fig.patch.set_alpha(1.0)
        ax.patch.set_facecolor('white')
        ax.patch.set_alpha(1.0)

        ax.set_xlim(self.grid_maker.pH_range)
        ax.set_ylim(self.grid_maker.V_range)
        if region_style == 'vector':
            self.plot_vector_regions(ax, species_grid, species_colors)
        else:
            color_image = self.make_color_image(species_grid, species_colors)
            color_image = self.upscale_region_image(color_image)
            ax.imshow(
                color_image,
                extent=[
                    self.grid_maker.pH_range[0],
                    self.grid_maker.pH_range[1],
                    self.grid_maker.V_range[0],
                    self.grid_maker.V_range[1],
                ],
                origin='lower',
                aspect='auto',
                interpolation='nearest',
                alpha=1.0,
                rasterized=True,
            )
        if label_stable_regions:
            self.label_stable_regions(ax, species_grid)
        
        ax.set_xlabel('pH')
        ax.set_ylabel(r'$E_{SHE}(V)$')
        self.add_H2_O2_lines(ax)
        self.add_plot_accessories(ax)
            
        return fig, ax

    def upscale_region_image(self, color_image):
        scale = max(1, REGION_IMAGE_MAX_PIXELS // max(color_image.shape[:2]))
        if scale == 1:
            return color_image
        return np.repeat(np.repeat(color_image, scale, axis=0), scale, axis=1)

    def make_color_image(self, species_grid, species_colors):
        """Converts a compact or legacy species grid to an opaque RGB image."""
        if hasattr(species_grid, 'min_indices') and hasattr(species_grid, 'species_list'):
            palette = np.zeros((len(species_grid.species_list), 3), dtype=np.uint8)
            for species_idx, species_tuple in enumerate(species_grid.species_list):
                label = self.format_species_tuple(species_tuple)
                if label not in species_colors:
                    continue
                rgb = np.array(to_rgba(species_colors[label])[:3], dtype=np.float64)
                rgb = REGION_COLOR_STRENGTH * rgb + (1 - REGION_COLOR_STRENGTH)
                palette[species_idx] = np.clip(np.round(rgb * 255), 0, 255).astype(np.uint8)
            return palette[species_grid.min_indices]

        labels = [
            self.format_species_tuple(species_tuple)
            for species_tuple in species_grid.ravel()
        ]
        rgb_values = np.array([to_rgba(species_colors[label])[:3] for label in labels])
        rgb_values = REGION_COLOR_STRENGTH * rgb_values + (1 - REGION_COLOR_STRENGTH)
        return np.clip(np.round(rgb_values.reshape(species_grid.shape + (3,)) * 255), 0, 255).astype(np.uint8)

    def make_index_grid_and_colors(self, species_grid, species_colors):
        """Returns integer region ids and opaque colors for vector/PDF output."""
        if hasattr(species_grid, 'min_indices') and hasattr(species_grid, 'species_list'):
            used_species_indices = np.unique(species_grid.min_indices)
            remapped_grid = np.empty(species_grid.min_indices.shape, dtype=np.int32)
            colors = []
            for new_idx, species_idx in enumerate(used_species_indices):
                remapped_grid[species_grid.min_indices == species_idx] = new_idx
                species_tuple = species_grid.species_list[species_idx]
                label = self.format_species_tuple(species_tuple)
                colors.append(to_rgba(species_colors[label], alpha=1.0))
            return remapped_grid, colors

        labels = np.empty(species_grid.size, dtype=object)
        for idx, species_tuple in enumerate(species_grid.ravel()):
            labels[idx] = self.format_species_tuple(species_tuple)
        labels = labels.reshape(species_grid.shape)
        label_list = sorted(set(labels.ravel()))
        label_to_idx = {label: idx for idx, label in enumerate(label_list)}
        index_grid = np.empty(labels.shape, dtype=np.int32)
        for label, idx in label_to_idx.items():
            index_grid[labels == label] = idx
        colors = [to_rgba(species_colors[label], alpha=1.0) for label in label_list]
        return index_grid, colors

    def plot_vector_regions(self, ax, species_grid, species_colors):
        index_grid, colors = self.make_index_grid_and_colors(species_grid, species_colors)
        levels = np.arange(len(colors) + 1) - 0.5
        ax.contourf(
            self.grid_maker.pH_values,
            self.grid_maker.V_values,
            index_grid,
            levels=levels,
            colors=colors,
            antialiased=False,
        )

    @staticmethod
    def format_species_tuple(species_tuple):
        return tuple(sorted(
            f"{species.formula}_{species.phase}_{species.alloy}"
            for species in species_tuple
        ))

    def format_formula(self, formula):
        if '[' in formula or 'aq' in formula or 'Gly' in formula:
            formatted_formula = re.sub(r"([A-Za-z\)\]])(\d+)", r"\1$_{\2}$", formula)
            formatted_formula = re.sub(r"\[([\d\+\-]+)\]", r"$^{\1}$", formatted_formula)
            formatted_formula = re.sub(r"\$_\{1\}\$", "", formatted_formula)
            formatted_formula = re.sub(r"\^\{1([+-])\}", r"^{\1}", formatted_formula)
            return formatted_formula

        formula_obj = Formula(formula)
        reduced_formula = formula_obj.reduce()[0]
        return f'{reduced_formula:latex}'

    @staticmethod
    def is_solid_pair_combo(combo_tuple):
        return len(combo_tuple) == 2 and all(species.rsplit('_', 2)[1] == 'solid' for species in combo_tuple)

    @staticmethod
    def is_grey_alloy_combo(combo_tuple):
        return all(species.rsplit('_', 2)[2] == 'True' for species in combo_tuple)

    def iter_grid_labels(self, species_grid):
        if hasattr(species_grid, 'min_indices') and hasattr(species_grid, 'species_list'):
            for species_idx, species_tuple in enumerate(species_grid.species_list):
                combo_label = self.format_species_tuple(species_tuple)
                yield combo_label, species_grid.min_indices == species_idx
            return

        labels = np.empty(species_grid.size, dtype=object)
        for idx, species_tuple in enumerate(species_grid.ravel()):
            labels[idx] = self.format_species_tuple(species_tuple)
        labels = labels.reshape(species_grid.shape)
        for combo_label in set(labels.ravel()):
            yield combo_label, labels == combo_label

    def format_combo_label(self, combo_tuple):
        labels = []
        complex_label_dict = getattr(self.species_data, 'species_label_dict', {})
        for species in combo_tuple:
            formula, phase, _ = species.rsplit('_', 2)
            if formula in complex_label_dict:
                formula = complex_label_dict[formula]
            labels.append(self.format_formula(formula) + ('(s)' if phase == 'solid' else '(aq)'))
        return '+'.join(labels)

    def should_label_combo(self, combo_tuple):
        return all(species.rsplit('_', 2)[1] == 'solid' for species in combo_tuple)

    def water_line_label_angle(self, ax):
        PREFAC = 0.0591
        pH0, pH1 = self.grid_maker.pH_range
        y0 = -pH0 * PREFAC + 1.23
        y1 = -pH1 * PREFAC + 1.23
        point0 = ax.transData.transform((pH0, y0))
        point1 = ax.transData.transform((pH1, y1))
        return np.degrees(np.arctan2(point1[1] - point0[1], point1[0] - point0[0]))

    def label_stable_regions(self, ax, species_grid, min_region_fraction=0.002):
        grid_area = self.grid_maker.grid_size ** 2
        min_pixels = max(20, int(grid_area * min_region_fraction))
        label_rotation = self.water_line_label_angle(ax)

        def place_label(combo_tuple, mask):
            rows, cols = np.nonzero(mask)
            if len(rows) < min_pixels:
                return

            center_row = rows.mean()
            center_col = cols.mean()
            label_idx = np.argmin((rows - center_row) ** 2 + (cols - center_col) ** 2)
            row = rows[label_idx]
            col = cols[label_idx]
            pH = self.grid_maker.pH_values[col]
            V = self.grid_maker.V_values[row]

            ax.text(
                pH, V, self.format_combo_label(combo_tuple),
                ha='center', va='center', color='black', fontsize=AREA_LABEL_FONT_SIZE,
                rotation=label_rotation, rotation_mode='anchor',
                zorder=5,
            )

        if hasattr(species_grid, 'min_indices') and hasattr(species_grid, 'species_list'):
            for species_idx, species_tuple in enumerate(species_grid.species_list):
                combo_tuple = self.format_species_tuple(species_tuple)
                if self.should_label_combo(combo_tuple):
                    place_label(combo_tuple, species_grid.min_indices == species_idx)
            return

        for combo_tuple, mask in self.iter_grid_labels(species_grid):
            if self.should_label_combo(combo_tuple):
                place_label(combo_tuple, mask)



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
        formula, phase, alloy = species.rsplit('_', 2)
        return formula, phase, alloy == 'True'

    @staticmethod
    def is_likely_metal(element):
        return element not in {'O', 'H', 'N', 'C'}

    def is_ligand_complex_formula(self, formula):
        complex_label_dict = getattr(self.species_data, 'species_label_dict', {})
        if formula in complex_label_dict:
            return True
        return any(ligand in formula for ligand in LIGAND_ELEMENT_COUNTS)

    def classify_combo(self, combo_tuple):
        phases = []
        has_oxide = False
        has_oxyhydroxide = False
        has_hydride = False
        has_ligand_complex = False
        has_aqueous_species = False
        is_alloy_combo = all(self.split_species_key(species)[2] for species in combo_tuple)

        for species in combo_tuple:
            formula, phase, _ = self.split_species_key(species)
            phases.append(phase)
            if phase != 'solid':
                has_aqueous_species = True

            if 'complex' in phase or self.is_ligand_complex_formula(formula):
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
        if has_aqueous_species:
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

    def metal_nonmetal_ratio_for_combo(self, combo_tuple):
        total_metals = 0
        total_other = 0

        for species in combo_tuple:
            formula, _, _ = self.split_species_key(species)
            element_counts = self.get_formula_element_counts(formula)
            for element, count in element_counts.items():
                if self.is_likely_metal(element):
                    total_metals += count
                else:
                    total_other += count

        if total_other == 0:
            # print(combo_tuple, total_met  . / / .als, total_other)
            return None 
        # print(combo_tuple, total_metals, total_other)
        return total_metals / total_other
    def color_intensity_for_combo(self, combo_tuple, category, used_intensities=None):
        """Calculate intensity based on metal ratio for better visualization."""
        metal_to_nonmetal = self.metal_nonmetal_ratio_for_combo(combo_tuple)

        # Define intensity ranges for each category (min, max)
        intensity_ranges = {
            'metal_ligand_complex': (0.55, 0.75),
            'metal': (0.30, 0.80),                    # Pure metals: high variation
            'metal_hydride': (0.38, 0.65),
            'metal_oxide': (0.24, 0.95),              # Handled separately
            'metal_oxyhydroxide': (0.24, 0.95),       # Handled separately
            'aqueous_metal_ion': (0.35, 0.85),
        }
        
        min_intensity, max_intensity = intensity_ranges.get(category, (0.28, 0.80))

        # For oxides/oxyhydroxides, use full metal ratio range
        if category == 'metal' and metal_to_nonmetal is None:
            intensity = 0.75

        elif category == 'metal_hydride':
            if metal_to_nonmetal is not None:
                # PdH has ratio ~5, should be medium-light (around 0.55)
                # Higher ratio = darker
                intensity = 0.40 + 0.25 * min(metal_to_nonmetal / 10.0, 1.0)
            else:
                intensity = 0.55

        elif category in ['metal_oxide', 'metal_oxyhydroxide']:
            if metal_to_nonmetal is None:
                intensity = 0.90  # Pure oxide (shouldn't happen, but fallback)
            elif metal_to_nonmetal is not None:
                # Higher metal ratio = darker orange
                # Normalize ratio: typical range 0.3-4 → map to 0-1
                normalized_ratio = min(metal_to_nonmetal / 5.0, 1.0)
                intensity = 0.30 + 0.68 * normalized_ratio
            else:
                intensity = 0.50
        elif category == 'metal':
            if metal_to_nonmetal is not None:
                # Has some non-metal, scale accordingly
                normalized_ratio = min(metal_to_nonmetal / 10.0, 1.0)
                intensity = 0.60 + 0.30 * normalized_ratio
            else:
                intensity = 0.70
        elif category == 'aqueous_metal_ion':
            if metal_to_nonmetal is not None:
                normalized_ratio = min(metal_to_nonmetal / 5.0, 1.0)
                intensity = 0.35 + 0.40 * normalized_ratio
            else:
                intensity = 0.55
        
        # METAL LIGAND COMPLEXES - medium
        elif category == 'metal_ligand_complex':
            intensity = 0.62
        
        else:
            intensity = 0.50
        intensity = min(0.98, max(0.28, intensity))
    
        if used_intensities is None:
            return intensity

        # For oxides, don't modify (already determined by ratio)
        if category in ['metal_oxide', 'metal_oxyhydroxide']:
            used_intensities.append(intensity)
            return intensity

        # For others, find suitable intensity with spacing
        spacing = 0.055
        step = 0.065
        max_attempts = 50
        attempt = 0

        while attempt < max_attempts:
            # Check if current intensity is far enough from all used intensities
            if all(abs(intensity - used) >= spacing for used in used_intensities):
                used_intensities.append(intensity)
                return intensity
            
            # Try adjusting intensity
            if metal_to_nonmetal is not None and metal_to_nonmetal > 2:
                # High metal ratio - prefer darker
                intensity += step
                if intensity > 0.95:
                    intensity = max(0.28, intensity - 2 * step)
            else:
                # Low metal ratio - prefer lighter
                intensity -= step
                if intensity < 0.28:
                    intensity = min(0.95, intensity + 2 * step)
            
            attempt += 1

        used_intensities.append(intensity)
        return intensity

        # if category in ['metal_oxide', 'metal_oxyhydroxide']:
        #     if metal_to_nonmetal is not None:
        #         # Map metal ratio directly to intensity
        #         # Higher metal ratio → higher intensity
        #         intensity = 0.24 + 0.76 * min(metal_to_nonmetal, 1.0)
        #     else:
        #         intensity = 0.50
            
        #     intensity = min(max_intensity, max(min_intensity, intensity))
        #     if used_intensities is None:
        #         return intensity
        #     used_intensities.append(intensity)
        #     return intensity

        # # For other categories, scale intensity with metal ratio
        # if metal_to_nonmetal is not None:
        #     metal_fraction = metal_to_nonmetal / (metal_to_nonmetal + 1)
        #     # Map metal fraction to intensity range for this category
        #     intensity = min_intensity + (max_intensity - min_intensity) * metal_fraction
        # else:
        #     # No metal content, use middle of range
        #     intensity = (min_intensity + max_intensity) / 2

        # intensity = min(max_intensity, max(min_intensity, intensity))
        
        # if used_intensities is None:
        #     return intensity

        # # Find suitable intensity with spacing
        # spacing = 0.055
        # step = 0.065
        # max_attempts = 50
        # attempt = 0

        # while attempt < max_attempts:
        #     # Check if current intensity is far enough from all used intensities
        #     if all(abs(intensity - used) >= spacing for used in used_intensities):
        #         used_intensities.append(intensity)
        #         return intensity
            
        #     # Try next intensity (prefer going up for higher metal content)
        #     if metal_to_nonmetal is not None and metal_to_nonmetal > 1:
        #         # High metal ratio - try going up first
        #         intensity += step
        #         if intensity > max_intensity:
        #             intensity = max(min_intensity, intensity - 2 * step)
        #     else:
        #         # Low metal ratio - try going down first
        #         intensity -= step
        #         if intensity < min_intensity:
        #             intensity = min(max_intensity, intensity + 2 * step)
            
        #     attempt += 1

        # # Fallback
        # used_intensities.append(intensity)
        # return intensity

    # def color_intensity_for_combo(self, combo_tuple, category, used_intensities=None):
    #     metal_to_nonmetal = self.metal_nonmetal_ratio_for_combo(combo_tuple)

    #     if category == 'metal_ligand_complex':
    #         intensity = 0.62
    #     elif category == 'metal':
    #         intensity = 0.36
    #     elif category == 'metal_hydride':
    #         intensity = 0.42
    #     elif category in ['metal_oxide', 'metal_oxyhydroxide'] and metal_to_nonmetal is not None:
    #         intensity = 0.24 + 0.76 * min(metal_to_nonmetal, 1.0)
    #     elif metal_to_nonmetal is not None:
    #         metal_fraction = metal_to_nonmetal / (metal_to_nonmetal + 1)
    #         intensity = 0.22 + 0.8 * metal_fraction
    #     else:
    #         intensity = 0.95

    #     intensity = min(0.95, max(0.28, intensity))
    #     if used_intensities is None:
    #         return intensity

    #     if category in ['metal_oxide', 'metal_oxyhydroxide']:
    #         used_intensities.append(intensity)
    #         return intensity

    #     while any(abs(intensity - used) < 0.055 for used in used_intensities):
    #         intensity += 0.1#0.065
    #         if intensity > 0.90:
    #             intensity = 0.25
    #     used_intensities.append(intensity)
    #     return intensity
        
    def get_color_for_label(self, all_species_tuples):
        """Assign colors with intensity scaling based on metal ratio."""
        species_colors = {}
        used_intensities = {}
        
        # Classify all combos
        categories = {
            combo_tuple: self.classify_combo(combo_tuple)
            for combo_tuple in sorted(all_species_tuples)
        }
        
        # Handle passivation layers separately (oxides/oxyhydroxides)
        passivation_combos = [
            combo_tuple for combo_tuple, category in categories.items()
            if category in ['metal_oxide', 'metal_oxyhydroxide']
        ]
        passivation_intensities = {}
        
        if passivation_combos:
            sorted_passivation = sorted(
                passivation_combos,
                key=lambda combo_tuple: (
                    self.metal_nonmetal_ratio_for_combo(combo_tuple) or 0,
                    combo_tuple,
                ),
            )
            # Assign intensities based on metal ratio (higher ratio = higher intensity)
            if len(sorted_passivation) == 1:
                combo_tuple = sorted_passivation[0]
                passivation_intensities[combo_tuple] = self.color_intensity_for_combo(
                    combo_tuple, categories[combo_tuple]
                )
            else:
                # For multiple oxides, space them based on metal content
                metal_ratios = [
                    (self.metal_nonmetal_ratio_for_combo(combo) or 0, combo)
                    for combo in sorted_passivation
                ]
                # Sort by metal ratio (ascending)
                metal_ratios.sort(key=lambda x: x[0])
                
                # Assign linearly spaced intensities (higher ratio = higher intensity)
                for intensity, (_, combo_tuple) in zip(
                    np.linspace(0.30, 0.98, len(sorted_passivation)),
                    metal_ratios,
                ):
                    passivation_intensities[combo_tuple] = intensity
        
        # Group non-passivation combos by category
        category_combos = {}
        for combo_tuple, category in categories.items():
            if combo_tuple not in passivation_intensities:
                if category not in category_combos:
                    category_combos[category] = []
                category_combos[category].append(combo_tuple)
        
        # Sort each category by metal ratio
        for category in category_combos:
            category_combos[category].sort(
                key=lambda combo: self.metal_nonmetal_ratio_for_combo(combo) or 0
            )
        
        # Assign colors to all species
        for combo_tuple, category in categories.items():
            species_colors[combo_tuple] = self.get_color_for_category(
                combo_tuple,
                category,
                used_intensities=used_intensities.setdefault(category, []),
                fixed_intensity=passivation_intensities.get(combo_tuple),
            )
        
        return species_colors

    # def get_color_for_label(self, all_species_tuples): 
    #     species_colors = {}
    #     used_intensities = {}
    #     categories = {
    #         combo_tuple: self.classify_combo(combo_tuple)
    #         for combo_tuple in sorted(all_species_tuples)
    #     }
    #     passivation_combos = [
    #         combo_tuple for combo_tuple, category in categories.items()
    #         if category in ['metal_oxide', 'metal_oxyhydroxide']
    #     ]
    #     passivation_intensities = {}
    #     if passivation_combos:
    #         sorted_passivation = sorted(
    #             passivation_combos,
    #             key=lambda combo_tuple: (
    #                 self.metal_nonmetal_ratio_for_combo(combo_tuple) or 0,
    #                 combo_tuple,
    #             ),
    #         )
    #         if len(sorted_passivation) == 1:
    #             combo_tuple = sorted_passivation[0]
    #             passivation_intensities[combo_tuple] = self.color_intensity_for_combo(
    #                 combo_tuple, categories[combo_tuple]
    #             )
    #         else:
    #             for intensity, combo_tuple in zip(
    #                 np.linspace(0.30, 0.98, len(sorted_passivation)),
    #                 sorted_passivation,
    #             ):
    #                 passivation_intensities[combo_tuple] = intensity
        
    #     for combo_tuple, category in categories.items():
    #         species_colors[combo_tuple] = self.get_color_for_category(
    #             combo_tuple,
    #             category,
    #             used_intensities.setdefault(category, []),
    #             fixed_intensity=passivation_intensities.get(combo_tuple),
    #         )
    #     return species_colors


    def get_color_for_category(self, combo_tuple, category, used_intensities=None, fixed_intensity=None):
        color_maps = {
            'metal': 'Greys',
            'metal_hydride': 'Greys',
            'metal_oxide': 'YlOrBr',
            'metal_oxyhydroxide': 'Oranges',
            'aqueous_metal_ion': 'Blues',
            'metal_ligand_complex': 'RdPu',
        }
        color_map = plt.cm.get_cmap(color_maps.get(category, 'Greys'))
        if fixed_intensity is None:
            intensity = self.color_intensity_for_combo(combo_tuple, category, used_intensities)
        else:
            intensity = fixed_intensity
        return color_map(intensity)
        
    def count_total_phases(self, all_species_tuples):
        total_solid = 0
        
        for combo_tuple in all_species_tuples:
            if all('solid' in species for species in combo_tuple):
                total_solid += 1
        total_aq = len(all_species_tuples) - total_solid
        return total_solid, total_aq
    


    
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
