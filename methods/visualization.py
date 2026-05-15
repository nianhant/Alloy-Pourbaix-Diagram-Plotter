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
        style = 'g-'
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

        ax.set_xlim(self.grid_maker.pH_range)
        ax.set_ylim(self.grid_maker.V_range)
        if region_style == 'vector':
            self.plot_vector_regions(ax, species_grid, species_colors)
        else:
            color_image = self.make_color_image(species_grid, species_colors)
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
                rasterized=True,
            )
        if label_stable_regions:
            self.label_stable_regions(ax, species_grid)
        
        ax.set_xlabel('pH')
        ax.set_ylabel(r'$E_{SHE}(V)$')
        self.add_H2_O2_lines(ax)
        self.add_plot_accessories(ax)
            
        return fig, ax

    def make_color_image(self, species_grid, species_colors):
        """Converts a compact or legacy species grid to an RGBA image."""
        if hasattr(species_grid, 'min_indices') and hasattr(species_grid, 'species_list'):
            palette = np.zeros((len(species_grid.species_list), 4), dtype=np.uint8)
            for species_idx, species_tuple in enumerate(species_grid.species_list):
                label = self.format_species_tuple(species_tuple)
                if label not in species_colors:
                    continue
                rgba = np.array(to_rgba(species_colors[label]), dtype=np.float64)
                rgba[3] *= 0.75
                palette[species_idx] = np.clip(np.round(rgba * 255), 0, 255).astype(np.uint8)
            return palette[species_grid.min_indices]

        labels = [
            self.format_species_tuple(species_tuple)
            for species_tuple in species_grid.ravel()
        ]
        rgba_values = np.array([to_rgba(species_colors[label], alpha=0.75) for label in labels])
        return np.clip(np.round(rgba_values.reshape(species_grid.shape + (4,)) * 255), 0, 255).astype(np.uint8)

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
        for species in combo_tuple:
            formula, phase, _ = species.rsplit('_', 2)
            labels.append(self.format_formula(formula) + ('(s)' if phase == 'solid' else '(aq)'))
        return '+'.join(labels)

    def should_label_combo(self, combo_tuple):
        return self.is_solid_pair_combo(combo_tuple) or self.is_grey_alloy_combo(combo_tuple)

    def label_stable_regions(self, ax, species_grid, min_region_fraction=0.002):
        grid_area = self.grid_maker.grid_size ** 2
        min_pixels = max(20, int(grid_area * min_region_fraction))

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

            text = ax.text(
                pH, V, self.format_combo_label(combo_tuple),
                ha='center', va='center', color='black', fontsize=20,
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
        
    def count_total_phases(self, all_species_tuples):
        total_solid = 0
        
        for combo_tuple in all_species_tuples:
            if all('solid' in species for species in combo_tuple):
                total_solid += 1
        total_aq = len(all_species_tuples) - total_solid
        return total_solid, total_aq
    
    def get_color_for_label(self, all_species_tuples): 
        species_colors = {}
        total_solid, total_aq = self.count_total_phases(all_species_tuples)
        
        warmer_color_map = plt.cm.get_cmap('summer')  # Warmer colors
        cooler_color_map = plt.cm.get_cmap('cool')  # Cooler colors
        grey_color = to_rgba('grey')
        solid_index = 0
        aq_index = 0
        
        for combo_tuple in sorted(all_species_tuples):
            new_combo_tuple_key = tuple(['_'.join(species.split('_')[:-1]) for species in combo_tuple ])
            
            if all(species.split('_')[-1] == 'True' for species in combo_tuple):
                species_colors[combo_tuple] = grey_color
            elif all('solid' in species for species in combo_tuple):
                solid_index += 1
                normalized_index = solid_index/max(total_solid, 1)
                species_colors[combo_tuple] = warmer_color_map(normalized_index)
            else:
                aq_index += 1
                normalized_index = aq_index/max(total_aq, 1)
                species_colors[combo_tuple] = cooler_color_map(normalized_index)
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
