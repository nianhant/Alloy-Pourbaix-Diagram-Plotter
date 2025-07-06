from matplotlib import pyplot as plt
import itertools
from pymatgen.core import Composition

from .data_loader import SpeciesDataLoader
from .thermodynamics import PourbaixData, PourbaixAnalyzer
from .grid import GridMaker
from .visualization import GridVisualizer, PlotAccessories
from .utils import parse_composition, format_comp_dict
from .plot_util import generate_legends
from .set_publication_style import set_publication_style

def process_alloys(metal_1, metal_list, CN_list, T_list, activity_list, mu_ligand, pH_exp_range=(11.5, 13.5), V_exp_range=(-2, 2.3), save_fig=True):
    set_publication_style()
    for metal_2 in metal_list: 
        for T in T_list:
            species_data = SpeciesDataLoader(metal_1, metal_2, mu_ligand, T)
            prod_comp_dict = {
                format_comp_dict(parse_composition(alloy)): Composition({
                    metal: count / sum(parse_composition(alloy).values())
                    for metal, count in parse_composition(alloy).items()
                }) 
                for alloy in species_data.solid_eng if len(parse_composition(alloy)) >= 2
            }

            for reference_alloy, reference_composition in prod_comp_dict.items():
                species_grid_list, all_species_tuples_global = [], set()

                for CN, activity in itertools.product(CN_list, activity_list):
                    ligand_concentration = {'NH3': 0.02, 'Gly': 0.005, 'CN': CN}
                    grid_maker = GridMaker((-2, 16), (-2, 3), ligand_concentration, 400)
                    pourbaix_data = PourbaixData(species_data, activity, ligand_concentration, reference_composition)
                    analyzer = PourbaixAnalyzer(pourbaix_data, grid_maker, T)

                    species_grid, all_species_tuples_set = analyzer.analyze_and_plot()
                    species_grid_list.append((species_grid, CN, T, activity, reference_alloy))
                    all_species_tuples_global.update(all_species_tuples_set)

                fig, ax = plt.subplots(figsize=(15, 8))
                species_colors = PlotAccessories(species_data).get_color_for_label(all_species_tuples_global)

                for species_grid, CN, T, activity, reference_alloy in species_grid_list:
                    fig, ax = GridVisualizer(grid_maker, pourbaix_data).plot_species_distribution(
                        species_grid, species_colors, ax=ax, save_fig=False)

                generate_legends(ax, all_species_tuples_global, species_colors, PlotAccessories(species_data), pH_exp_range, V_exp_range)

                plt.tight_layout()
                if save_fig:
                    file_name = GridVisualizer(grid_maker, pourbaix_data).format_file_name()
                    print(file_name)
                    plt.savefig(file_name, bbox_inches='tight')
                plt.show()
