
from matplotlib import pyplot as plt
import itertools
from pymatgen.core import Composition
from pathlib import Path

from .data_loader import SpeciesDataLoader
from .thermodynamics import PourbaixData, PourbaixAnalyzer
from .grid import GridMaker
from .visualization import GridVisualizer, PlotAccessories
from .plot_util import generate_legends
from .utils import parse_composition, format_comp_dict
from .set_publication_style import set_publication_style


set_publication_style()

BASE_FIGURE_SIZE_IN = (8.0, 8.0)
PLOT_LEFT_IN = 1.28
PLOT_BOTTOM_IN = 1.12
PLOT_WIDTH_IN = 6.40
PLOT_HEIGHT_IN = 6.56
CANVAS_PAD_IN = 0.16


def set_fixed_axes_position(fig, ax, bottom_in=PLOT_BOTTOM_IN):
    width_in, height_in = fig.get_size_inches()
    ax.set_position([
        PLOT_LEFT_IN / width_in,
        bottom_in / height_in,
        PLOT_WIDTH_IN / width_in,
        PLOT_HEIGHT_IN / height_in,
    ])


def make_fixed_plot_figure():
    fig = plt.figure(figsize=BASE_FIGURE_SIZE_IN)
    ax = fig.add_axes([0, 0, 1, 1])
    set_fixed_axes_position(fig, ax)
    return fig, ax


def expand_canvas_for_artists(fig, ax, artists):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width_in, height_in = fig.get_size_inches()

    bboxes = [
        artist.get_window_extent(renderer).transformed(fig.dpi_scale_trans.inverted())
        for artist in artists
        if artist.get_visible()
    ]
    if not bboxes:
        return

    min_y = min(bbox.y0 for bbox in bboxes)
    max_y = max(bbox.y1 for bbox in bboxes)
    max_x = max(bbox.x1 for bbox in bboxes)

    bottom_shift_in = max(0.0, CANVAS_PAD_IN - min_y)
    new_width_in = max(width_in, max_x + CANVAS_PAD_IN)
    new_height_in = max(height_in + bottom_shift_in, max_y + bottom_shift_in + CANVAS_PAD_IN)

    if new_width_in != width_in or new_height_in != height_in or bottom_shift_in:
        fig.set_size_inches(new_width_in, new_height_in)
        set_fixed_axes_position(fig, ax, PLOT_BOTTOM_IN + bottom_shift_in)
        fig.canvas.draw()


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

        fig, ax = make_fixed_plot_figure()
        plot_accessories = PlotAccessories(species_data)
        species_colors = plot_accessories.get_color_for_label(all_species_tuples_global)

        fig, ax = GridVisualizer(grid_maker, pourbaix_data).plot_species_distribution(
            species_grid, species_colors, ax=ax, save_fig=False,
            label_stable_regions=False,
            pH_exp_range=pH_exp_range, V_exp_range=V_exp_range,
        )
        legend_artists = generate_legends(
            ax, all_species_tuples_global, species_colors, plot_accessories,
            pH_exp_range, V_exp_range,
        )
        expand_canvas_for_artists(fig, ax, legend_artists)

        if save_fig:
            file_name = GridVisualizer(grid_maker, pourbaix_data).format_file_name()
            output_path = Path(outdir) / f"{metal_1}_{metal_2}"
            output_path.mkdir(parents=True, exist_ok=True)

            file_stem = Path(file_name).stem
            output_file = output_path / f"{file_stem}.png"

            # print(file_name)
            fig.savefig(
                output_file,
                transparent=False,
                facecolor='white',
                dpi=600,
            )
            if save_pdf:
                plt.close(fig)
                fig, ax = make_fixed_plot_figure()
                GridVisualizer(grid_maker, pourbaix_data).plot_species_distribution(
                    species_grid, species_colors, ax=ax, save_fig=False,
                    label_stable_regions=False, region_style='image',
                    pH_exp_range=pH_exp_range, V_exp_range=V_exp_range,
                )
                legend_artists = generate_legends(
                    ax, all_species_tuples_global, species_colors, plot_accessories,
                    pH_exp_range, V_exp_range,
                )
                expand_canvas_for_artists(fig, ax, legend_artists)
                fig.savefig(
                    output_path / f"{file_stem}.pdf",
                    transparent=False,
                    facecolor='white',
                    dpi=600,
                )
                fig.savefig(
                    output_path / f"{file_stem}.svg",
                    format="svg",
                )
        plt.close()
