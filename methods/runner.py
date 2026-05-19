
from matplotlib import pyplot as plt
import itertools
from pymatgen.core import Composition
from pathlib import Path

from .data_loader import SpeciesDataLoader
from .thermodynamics import PourbaixData, PourbaixAnalyzer
from .grid import GridMaker
from .visualization import GridVisualizer, PlotAccessories
from .utils import parse_composition, format_comp_dict
from .set_publication_style import set_publication_style


set_publication_style()


def plot_pourbaix(metal_1, metal_2, mu_ligand, T, activity, ligand_concentration, data_dir, pH_exp_range, V_exp_range, save_fig, outdir, save_pdf=True):
    species_data = SpeciesDataLoader(metal_1, metal_2, mu_ligand, T, data_dir = data_dir)
    prod_comp_dict = {
    format_comp_dict(parse_composition(alloy)): Composition({
        metal: count / sum(parse_composition(alloy).values())
        for metal, count in parse_composition(alloy).items()
    }) 
    for alloy in species_data.solid_eng if len(parse_composition(alloy)) >= 2 and 
    all(count / sum(parse_composition(alloy).values()) == 0.5 for count in parse_composition(alloy).values())
    }
    for reference_alloy, reference_composition in prod_comp_dict.items():
        species_grid_list, all_species_tuples_global = [], set()

        grid_maker = GridMaker((-2, 16), (-2, 3), ligand_concentration, 800)
        pourbaix_data = PourbaixData(species_data, activity, ligand_concentration, reference_composition)
        analyzer = PourbaixAnalyzer(pourbaix_data, grid_maker, T)

        species_grid, all_species_tuples_set = analyzer.analyze_and_plot()
        all_species_tuples_global.update(all_species_tuples_set)

        fig, ax = plt.subplots(figsize=(8, 8))
        fig.subplots_adjust(left=0.16, right=0.96, bottom=0.14, top=0.96)
        species_colors = PlotAccessories(species_data).get_color_for_label(all_species_tuples_global)

        fig, ax = GridVisualizer(grid_maker, pourbaix_data).plot_species_distribution(
                            species_grid, species_colors, ax=ax, save_fig=False)

        if save_fig:
            file_name = GridVisualizer(grid_maker, pourbaix_data).format_file_name()
            output_path = Path(outdir) / f"{metal_1}_{metal_2}"
            output_path.mkdir(parents=True, exist_ok=True)

            file_stem = Path(file_name).stem
            output_file = output_path / f"{file_stem}.png"

            # print(file_name)
            plt.savefig(output_file, bbox_inches='tight', transparent=False, facecolor='white', dpi=600)
            if save_pdf:
                plt.close(fig)
                fig, ax = plt.subplots(figsize=(8, 8))
                fig.subplots_adjust(left=0.16, right=0.96, bottom=0.14, top=0.96)
                GridVisualizer(grid_maker, pourbaix_data).plot_species_distribution(
                    species_grid, species_colors, ax=ax, save_fig=False,
                    region_style='image',
                )
                plt.savefig(
                    output_path / f"{file_stem}.pdf",
                    bbox_inches='tight',
                    transparent=False,
                    facecolor='white',
                    dpi=600,
                )
                plt.savefig(output_path / f"{file_stem}.svg", format="svg", bbox_inches="tight")
        plt.close()
