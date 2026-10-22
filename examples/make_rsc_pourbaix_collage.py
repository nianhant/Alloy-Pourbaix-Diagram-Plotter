"""Generate an RSC-sized TiNi/TiCu/PdAu Pourbaix figure collage."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.ticker import MultipleLocator
from matplotlib.transforms import Bbox
from PIL import Image
from pymatgen.core import Composition

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from methods.data_loader import SpeciesDataLoader
from methods.grid import GridMaker
from methods.plot_util import categorized_legend_entries
from methods.thermodynamics import PourbaixAnalyzer, PourbaixData
from methods.utils import format_comp_dict, parse_composition
from methods.visualization import (
    REGION_COLOR_STRENGTH,
    GridVisualizer,
    PlotAccessories,
)


SYSTEMS = (
    ("Ni", "Ti", "TiNi"),
    ("Cu", "Ti", "TiCu"),
    ("Pd", "Au", "PdAu"),
)
MU_LIGAND = {"NH3": -0.276037, "Gly": -3.263014109, "CN": 1.786800089}
TEMPERATURE_K = 298.15
ION_ACTIVITY = 1e-4
LIGAND_CONCENTRATION = {"NH3": 0.02, "NO2": 0.0, "Gly": 0.1, "CN": 0.0}

PH_RANGE = (-2, 16)
V_RANGE = (-2, 3)
PH_EXP_RANGE = (10.5, 13.5)
V_EXP_RANGE = (-2, 2.0)
GRID_SIZE = 1000

DATA_DIR = REPO_ROOT / "data"
OUTPUT_DIR = REPO_ROOT / "figures" / "alloy_pourbaix_diagrams" / "rsc_collage"
OUTPUT_STEM = "TiNi_TiCu_PdAu_NH3=0.02M_Gly=0.1M_activity=1e-4M"

# Overwide side-key geometry matching the single-metal Pourbaix panel script.
CM_TO_INCH = 1 / 2.54
FIGURE_WIDTH_CM = 29.15
FIGURE_HEIGHT_CM = 6.5
FIGURE_SIZE_IN = (FIGURE_WIDTH_CM * CM_TO_INCH, FIGURE_HEIGHT_CM * CM_TO_INCH)
COLUMN_WIDTH_CM = 9.65
COLUMN_STRIDE_CM = 9.75
ROW_BOTTOM_CM = 0.6
PLOT_LEFT_OFFSET_CM = 1.0
PLOT_SIZE_CM = 4.7
KEY_LEFT_OFFSET_CM = 5.85
KEY_RIGHT_MARGIN_CM = 0.1
KEY_HEIGHT_CM = 5.2

# Displayed swatches sampled from the single-metal six-panel side-color SVG
# (Au_Ni_Cu_Pd_Pt_Ti_*_2x3_overwide_side_colors_all_axes.svg).
# Keys use ligand and stoichiometry to avoid differences in charge typography.
COMPLEX_DISPLAY_COLORS = {
    ("Ni", "Gly", 1): "#fa9db6",
    ("Ni", "Gly", 2): "#f875a8",
    ("Ni", "Gly", 3): "#e2459e",
    ("Cu", "Gly", 2): "#fba9ba",
    # Not stable in the reference six-panel figure: use a distinct RdPu shade.
    ("Cu", "Gly", 3): "#f875a8",
    ("Au", "NH3", 2): "#fa99b5",
}


def match_complex_colors(species_grid, species_colors):
    """Match reference display colors, compensating for the region white blend."""
    for index in set(species_grid.min_indices.ravel()):
        species_tuple = species_grid.species_list[index]
        complexes = [species for species in species_tuple if species.phase == "complex"]
        if not complexes:
            continue
        display_colors = []
        for species in complexes:
            metal = next(element for element in ("Ni", "Cu", "Au", "Ti", "Pd")
                         if species.composition.get(element, 0))
            key = (metal, species.ligand, int(species.composition[species.ligand]))
            if key not in COMPLEX_DISPLAY_COLORS:
                raise ValueError(f"Define a reference complex color for {key}")
            display_colors.append(to_rgba(COMPLEX_DISPLAY_COLORS[key])[:3])
        target = tuple(sum(rgb[i] for rgb in display_colors) / len(display_colors)
                       for i in range(3))
        combo = species_grid.format_species_tuple(species_tuple)
        species_colors[combo] = tuple(
            (channel - (1 - REGION_COLOR_STRENGTH)) / REGION_COLOR_STRENGTH
            for channel in target
        ) + (1.0,)


def equiatomic_reference_composition(species_data):
    """Return the unique equiatomic binary composition in a system dataset."""
    compositions = {}
    for formula in species_data.solid_eng:
        parsed = parse_composition(formula)
        total = sum(parsed.values())
        if len(parsed) != 2 or not all(count / total == 0.5 for count in parsed.values()):
            continue
        key = format_comp_dict(parsed)
        compositions[key] = Composition(
            {element: count / total for element, count in parsed.items()}
        )
    if len(compositions) != 1:
        raise ValueError(
            "Expected one equiatomic binary reference composition, found "
            f"{sorted(compositions)}."
        )
    return next(iter(compositions.values()))


def calculate_panel(metal_1, metal_2):
    """Calculate one alloy's stable-region grid and plotting metadata."""
    species_data = SpeciesDataLoader(
        metal_1, metal_2, MU_LIGAND, TEMPERATURE_K, data_dir=DATA_DIR
    )
    grid_maker = GridMaker(PH_RANGE, V_RANGE, LIGAND_CONCENTRATION, GRID_SIZE)
    pourbaix_data = PourbaixData(
        species_data,
        ION_ACTIVITY,
        LIGAND_CONCENTRATION,
        equiatomic_reference_composition(species_data),
    )
    species_grid, stable_combinations = PourbaixAnalyzer(
        pourbaix_data, grid_maker, TEMPERATURE_K
    ).analyze_and_plot()
    accessories = PlotAccessories(species_data)
    colors = accessories.get_color_for_label(stable_combinations)
    match_complex_colors(species_grid, colors)
    return grid_maker, pourbaix_data, species_grid, stable_combinations, accessories, colors


def style_panel(ax, panel_index, title):
    """Apply compact, print-legible styling to one collage panel."""
    ax.set_title(f"({chr(97 + panel_index)}) {title}", loc="left", pad=2)
    ax.set_xlabel("pH", labelpad=2)
    ax.set_ylabel(r"$E_{\mathrm{SHE}}$ (V)", labelpad=1)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.tick_params(axis="both", which="both", direction="in", length=2.5, width=0.7, pad=1.5)
    ax.set_box_aspect(1)
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)


def add_side_key(fig, column_left_cm, species_grid, grid_maker, colors, accessories, ax,
                 include_aqueous=False):
    """Draw a key ordered by category and descending region centroid height."""
    key_ax = fig.add_axes((
        (column_left_cm + KEY_LEFT_OFFSET_CM) / FIGURE_WIDTH_CM,
        ROW_BOTTOM_CM / FIGURE_HEIGHT_CM,
        (COLUMN_WIDTH_CM - KEY_LEFT_OFFSET_CM - KEY_RIGHT_MARGIN_CM)
        / FIGURE_WIDTH_CM,
        KEY_HEIGHT_CM / FIGURE_HEIGHT_CM,
    ))
    key_ax.set_axis_off()
    entries = categorized_legend_entries(
        species_grid, grid_maker.V_values, colors, accessories,
        ph_values=grid_maker.pH_values, include_aqueous=include_aqueous,
    )
    vertical_span = 0.94 if set(accessories.species_data.metal_reference) == {"Ni", "Ti"} else 0.82
    step = min(0.075, vertical_span / max(len(entries), 1))
    if include_aqueous:
        # Maintain readable line spacing; enlarge only the legend/canvas, not
        # the square scientific axes. Existing reduced legends are unchanged.
        row_points = 11.0
        height_cm = max(KEY_HEIGHT_CM, len(entries) * row_points / 72 * 2.54)
        position = key_ax.get_position()
        height = height_cm / FIGURE_HEIGHT_CM
        key_ax.set_position((position.x0, position.y0 + (position.height-height)/2,
                             position.width, height))
        step = row_points / (height_cm / 2.54 * 72)
    first_y = 0.5 + (len(entries) - 1) * step / 2
    for entry_index, (label, color, centroid) in enumerate(entries):
        y = first_y - entry_index * step
        rgba = to_rgba(color)
        blended_color = tuple(
            REGION_COLOR_STRENGTH * channel + (1 - REGION_COLOR_STRENGTH)
            for channel in rgba[:3]
        )
        key_ax.plot(
            0.04, y, marker="o", markersize=5,
            markerfacecolor=blended_color, markeredgecolor="none",
            linestyle="none", transform=key_ax.transAxes,
        )
        key_ax.text(
            0.09, y, label, transform=key_ax.transAxes,
            ha="left", va="center", fontsize=7,
        )
    return key_ax


def compact_panel_spacing(fig, panel_pairs):
    """Pack measured plot/legend groups without changing their physical size."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_inches = fig.dpi_scale_trans.inverted()
    bounds = []
    positions = []
    for ax, key_ax in panel_pairs:
        artists = [*key_ax.texts, *key_ax.lines]
        boxes = [ax.get_tightbbox(renderer)] + [
            artist.get_window_extent(renderer) for artist in artists
        ]
        bounds.append(Bbox.union(boxes).transformed(to_inches))
        positions.append([
            axis.get_position().transformed(fig.transFigure).transformed(to_inches)
            for axis in (ax, key_ax)
        ])

    margin = 0.12 / 2.54
    gap = 0.25 / 2.54
    bottom = min(box.y0 for box in bounds)
    top = max(box.y1 for box in bounds)
    height = top - bottom + 2 * margin
    width = sum(box.width for box in bounds) + gap * (len(bounds) - 1) + 2 * margin
    fig.set_size_inches(width, height)
    cursor = margin
    for pair, old_positions, box in zip(panel_pairs, positions, bounds):
        for axis, pos in zip(pair, old_positions):
            axis.set_position((
                (pos.x0 + cursor - box.x0) / width,
                (pos.y0 + margin - bottom) / height,
                pos.width / width, pos.height / height,
            ))
        cursor += box.width + gap
    print(f"Packed figure: {width * 2.54:.2f} x {height * 2.54:.2f} cm", flush=True)


def main():
    """Create TIFF, PDF, and SVG versions of the three-panel collage."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    style = {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "legend.title_fontsize": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }

    with plt.rc_context(style):
        fig = plt.figure(figsize=FIGURE_SIZE_IN, facecolor="white")
        panel_pairs = []
        for panel_index, (metal_1, metal_2, title) in enumerate(SYSTEMS):
            column_left_cm = panel_index * COLUMN_STRIDE_CM
            plot_bottom_cm = ROW_BOTTOM_CM + (KEY_HEIGHT_CM - PLOT_SIZE_CM) / 2
            ax = fig.add_axes((
                (column_left_cm + PLOT_LEFT_OFFSET_CM) / FIGURE_WIDTH_CM,
                plot_bottom_cm / FIGURE_HEIGHT_CM,
                PLOT_SIZE_CM / FIGURE_WIDTH_CM,
                PLOT_SIZE_CM / FIGURE_HEIGHT_CM,
            ))
            (
                grid_maker,
                pourbaix_data,
                species_grid,
                stable_combinations,
                accessories,
                colors,
            ) = calculate_panel(metal_1, metal_2)
            GridVisualizer(grid_maker, pourbaix_data).plot_species_distribution(
                species_grid,
                colors,
                ax=ax,
                save_fig=False,
                label_stable_regions=False,
                region_style="image",
                pH_exp_range=PH_EXP_RANGE,
                V_exp_range=V_EXP_RANGE,
            )
            style_panel(ax, panel_index, title)
            key_ax = add_side_key(
                fig, column_left_cm, species_grid, grid_maker, colors, accessories, ax
            )
            panel_pairs.append((ax, key_ax))

        compact_panel_spacing(fig, panel_pairs)
        output_base = OUTPUT_DIR / OUTPUT_STEM
        tiff_path = Path(f"{output_base}.tiff")
        fig.savefig(
            tiff_path,
            dpi=600,
            facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"},
        )
        with Image.open(tiff_path) as tiff_image:
            rgb_tiff = tiff_image.convert("RGB")
        rgb_tiff.save(tiff_path, dpi=(600, 600), compression="tiff_lzw")
        fig.savefig(Path(f"{output_base}.pdf"), facecolor="white", dpi=600)
        fig.savefig(Path(f"{output_base}.svg"), facecolor="white", dpi=600)
        plt.close(fig)


if __name__ == "__main__":
    main()
