from matplotlib.patches import Patch
import numpy as np


def categorized_legend_entries(
    species_grid, potentials, species_colors, plot_accessories, ph_values=None,
    include_aqueous=False,
):
    """Order selected regions by family, then descending mean grid potential.

    Disconnected parts of the same region contribute to one area-weighted
    centroid. Simple aqueous-ion regions are omitted unless include_aqueous.
    If ph_values is supplied, append the (pH, potential) centroid to each entry.
    """
    ranks = {
        "metal_oxide": 0, "metal_oxyhydroxide": 0,
        "metal_ligand_complex": 1, "metal": 2, "metal_hydride": 2,
    }
    if include_aqueous:
        ranks['aqueous_metal_ion'] = 3
    records = []
    for index in np.unique(species_grid.min_indices):
        combo = species_grid.format_species_tuple(species_grid.species_list[index])
        category = plot_accessories.classify_combo(combo)
        if category not in ranks:
            continue
        rows, columns = np.nonzero(species_grid.min_indices == index)
        centroid = float(np.mean(potentials[rows]))
        ph_centroid = float(np.mean(ph_values[columns])) if ph_values is not None else None
        records.append((ranks[category], -centroid, combo, ph_centroid))
    records.sort()
    return [
        (plot_accessories.format_species_combo(combo), species_colors[combo])
        + (((ph_centroid, -negative_height),) if ph_values is not None else ())
        for _, negative_height, combo, ph_centroid in records
    ]


def solid_legend_entries(all_species_tuples, species_colors, plot_accessories):
    """Return deterministic ``(label, color)`` entries for all-solid regions."""
    solid_species = [
        species
        for species in all_species_tuples
        if all(item.rsplit("_", 2)[1] == "solid" for item in species)
    ]
    solid_species.sort(key=plot_accessories.format_species_combo)
    return [
        (plot_accessories.format_species_combo(species), species_colors[species])
        for species in solid_species
    ]


def solid_legend_handles(all_species_tuples, species_colors, plot_accessories):
    """Return deterministic fill swatches for stable all-solid regions."""
    return [
        Patch(
            facecolor=color,
            edgecolor="0.25",
            linewidth=0.4,
            label=label,
        )
        for label, color in solid_legend_entries(
            all_species_tuples, species_colors, plot_accessories
        )
    ]


def generate_legends(
    ax,
    all_species_tuples_global,
    species_colors,
    plot_accessories,
    pH_exp_range,
    V_exp_range,
):
    """Add a compact legend for solid stability regions only."""
    solid_handles = solid_legend_handles(
        all_species_tuples_global, species_colors, plot_accessories
    )

    combined_legend = ax.legend(
        handles=solid_handles,
        title="Solid species",
        loc="center left", frameon=True, bbox_to_anchor=(1.01, 0.5),
        handlelength=1.0, handletextpad=0.5,
    )
    return (combined_legend,)
