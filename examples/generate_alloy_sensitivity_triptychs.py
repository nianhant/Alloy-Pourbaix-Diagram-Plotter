"""Reproduce the single-metal OAT concentration sweep for TiNi, TiCu, PdAu."""

import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba, to_hex
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from examples import make_rsc_pourbaix_collage as collage
from methods.data_loader import SpeciesDataLoader
from methods.grid import GridMaker
from methods.thermodynamics import PourbaixData, PourbaixAnalyzer
from methods.visualization import GridVisualizer, PlotAccessories, REGION_COLOR_STRENGTH

REFERENCE = Path("/global/homes/n/nianhant/data/stability_paper/single-metal-pourbaix-plotter/examples/generate_sensitivity_triptychs.py")
BASELINE = {"activity": 1e-4, "NH3": 0.02, "Gly": 0.1, "CN": 0.0}
LEVELS = {
    "activity": (1e-6, 1e-4, 1e-2),
    "NH3": (2e-3, 2e-2, 1e-1),
    "Gly": (1e-2, 1e-1, 1.0),
    "CN": (0.0, 1e-4, 1e-3),
}
LABELS = {
    "activity": r"$[\mathrm{M}]_{\mathrm{aq}}$",
    "NH3": r"$[\mathrm{NH_3}]_{\mathrm{tot}}$",
    "Gly": r"$[\mathrm{Gly}]_{\mathrm{tot}}$",
    "CN": r"$[\mathrm{CN}]_{\mathrm{tot}}$",
}
CATEGORIES = ("metal_oxide", "metal_oxyhydroxide", "metal_ligand_complex",
              "metal", "metal_hydride", "aqueous_metal_ion")
FIGURE_ROOT = ROOT / "figures/alloy_pourbaix_diagrams/sensitivity_triptychs_with_cn_sweep"
RESULT_ROOT = ROOT / "research/2026-10-02_alloy_concentration_sensitivity"


def verify_reference() -> dict:
    """Check exact constants without importing/executing the other project."""
    tree = ast.parse(REFERENCE.read_text())
    found = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in {"BASELINE", "LEVELS"}:
                    found[target.id] = ast.literal_eval(node.value)
    assert found == {"BASELINE": BASELINE, "LEVELS": LEVELS}, found
    return {"path": str(REFERENCE), "sha256": hashlib.sha256(REFERENCE.read_bytes()).hexdigest()}


def condition_key(condition: dict) -> tuple:
    """Provide an immutable cache key in the documented parameter order."""
    return tuple(condition[name] for name in BASELINE)


def calculate_case(data: SpeciesDataLoader, condition: dict, size: int) -> dict:
    """Minimize the existing model for one fully specified condition."""
    ligands = {name: condition[name] for name in ("NH3", "Gly", "CN")}
    ligands["NO2"] = 0.0
    grid = GridMaker(collage.PH_RANGE, collage.V_RANGE, ligands, size)
    thermo = PourbaixData(data, condition["activity"], ligands,
                         collage.equiatomic_reference_composition(data))
    regions, stable = PourbaixAnalyzer(thermo, grid, collage.TEMPERATURE_K).analyze_and_plot()
    return {"condition": dict(condition), "grid": grid, "thermo": thermo,
            "regions": regions, "stable": stable}


def consistent_colors(cases: list, accessories: PlotAccessories) -> dict:
    """Assign a single color per combination over the entire alloy sweep."""
    union = set().union(*(case["stable"] for case in cases))
    colors = accessories.get_color_for_label(union)
    complex_species = {}
    for case in cases:
        result = case["regions"]
        for index in np.unique(result.min_indices):
            for species in result.species_list[index]:
                if species.phase == "complex":
                    complex_species[species.formula] = species
    complex_colors = {}
    unknown = []
    for formula, species in sorted(complex_species.items()):
        metals = accessories.species_data.metal_reference
        metal = next(element for element in metals if species.composition.get(element, 0))
        key = (metal, species.ligand, int(species.composition[species.ligand]))
        if key in collage.COMPLEX_DISPLAY_COLORS:
            complex_colors[formula] = to_rgba(collage.COMPLEX_DISPLAY_COLORS[key])[:3]
        else:
            unknown.append(formula)
    for formula, intensity in zip(unknown, np.linspace(0.38, 0.68, max(len(unknown), 1))):
        raw = plt.get_cmap("RdPu")(intensity)[:3]
        complex_colors[formula] = tuple(REGION_COLOR_STRENGTH * c + 1 - REGION_COLOR_STRENGTH for c in raw)
    for case in cases:
        result = case["regions"]
        for index in np.unique(result.min_indices):
            species_tuple = result.species_list[index]
            complexes = [s for s in species_tuple if s.phase == "complex"]
            if not complexes:
                continue
            target = np.mean([complex_colors[s.formula] for s in complexes], axis=0)
            raw = (target - (1 - REGION_COLOR_STRENGTH)) / REGION_COLOR_STRENGTH
            assert np.all(raw >= 0) and np.all(raw <= 1)
            colors[result.format_species_tuple(species_tuple)] = (*raw, 1.0)
    return colors


def summarize(case: dict, baseline: dict, accessories: PlotAccessories) -> dict:
    """Compute categorical area and stable-combination changes from baseline."""
    result, grid = case["regions"], case["grid"]
    labels = [result.format_species_tuple(s) for s in result.species_list]
    base = baseline["regions"]
    assert labels == [base.format_species_tuple(s) for s in base.species_list]
    changed = result.min_indices != base.min_indices
    rhe = grid.V_values[:, None] + 0.0591 * grid.pH_values[None, :]
    window = ((grid.pH_values[None, :] >= collage.PH_EXP_RANGE[0])
              & (grid.pH_values[None, :] <= collage.PH_EXP_RANGE[1])
              & (rhe >= collage.V_EXP_RANGE[0]) & (rhe <= collage.V_EXP_RANGE[1]))
    scopes = {"full": np.ones(result.shape, dtype=bool), "experimental": window}
    stats = {}
    ranks = {name: index for index, name in enumerate(CATEGORIES)}
    category_ids = np.array([ranks[accessories.classify_combo(label)] for label in labels])
    categorized = category_ids[result.min_indices]
    for scope, mask in scopes.items():
        assert np.any(mask)
        fractions = np.bincount(categorized[mask], minlength=len(CATEGORIES)) / np.count_nonzero(mask)
        assert np.isclose(fractions.sum(), 1)
        stats[f"{scope}_changed_fraction"] = float(np.mean(changed[mask]))
        for name, fraction in zip(CATEGORIES, fractions):
            stats[f"{scope}_{name}"] = float(fraction)
        stats[f"{scope}_aqueous_containing"] = float(
            fractions[ranks["metal_ligand_complex"]] + fractions[ranks["aqueous_metal_ion"]])
    if case["condition"] == BASELINE:
        assert stats["full_changed_fraction"] == 0
    return stats


def render(alloy: str, parameter: str, cases: list, accessories: PlotAccessories,
           colors: dict, output: Path, titles: list | None = None,
           output_stem: str | None = None, include_aqueous: bool = False) -> dict:
    """Render conditions with compact side legends and optional SI titles."""
    if titles is not None:
        assert len(titles) == len(cases)
    style = {"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
             "font.size": 7, "axes.labelsize": 8, "axes.titlesize": 8,
             "xtick.labelsize": 7, "ytick.labelsize": 7,
             "pdf.fonttype": 42, "ps.fonttype": 42}
    with plt.rc_context(style):
        fig = plt.figure(figsize=collage.FIGURE_SIZE_IN, facecolor="white")
        pairs = []
        for index, case in enumerate(cases):
            left = index * collage.COLUMN_STRIDE_CM
            bottom = collage.ROW_BOTTOM_CM + (collage.KEY_HEIGHT_CM - collage.PLOT_SIZE_CM) / 2
            ax = fig.add_axes(((left + collage.PLOT_LEFT_OFFSET_CM) / collage.FIGURE_WIDTH_CM,
                               bottom / collage.FIGURE_HEIGHT_CM,
                               collage.PLOT_SIZE_CM / collage.FIGURE_WIDTH_CM,
                               collage.PLOT_SIZE_CM / collage.FIGURE_HEIGHT_CM))
            GridVisualizer(case["grid"], case["thermo"]).plot_species_distribution(
                case["regions"], colors, ax=ax, save_fig=False, label_stable_regions=False,
                pH_exp_range=collage.PH_EXP_RANGE, V_exp_range=collage.V_EXP_RANGE)
            title = (titles[index] if titles is not None else
                     f"{LABELS[parameter]} = {case['condition'][parameter]:g} M")
            collage.style_panel(ax, index, title)
            key_ax = collage.add_side_key(fig, left, case["regions"], case["grid"], colors, accessories, ax,
                                          include_aqueous=include_aqueous)
            if include_aqueous:
                assert len(key_ax.texts) == len(np.unique(case['regions'].min_indices))
            pairs.append((ax, key_ax))
        collage.compact_panel_spacing(fig, pairs)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for ax, key_ax in pairs:
            extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            assert np.isclose(extent.width * 2.54, 4.7) and np.isclose(extent.height * 2.54, 4.7)
            assert tuple(ax.get_xlim()) == collage.PH_RANGE and tuple(ax.get_ylim()) == collage.V_RANGE
            assert ax.xaxis.get_tick_params()["direction"] == ax.yaxis.get_tick_params()["direction"] == "in"
            for artist in key_ax.texts:
                box = artist.get_window_extent(renderer)
                assert box.x0 >= 0 and box.x1 <= fig.bbox.width
                assert box.y0 >= 0 and box.y1 <= fig.bbox.height
            if include_aqueous:
                boxes = [artist.get_window_extent(renderer) for artist in key_ax.texts]
                assert all(upper.y0 > lower.y1 for upper, lower in zip(boxes, boxes[1:])), 'Legend rows overlap'
        output.mkdir(parents=True, exist_ok=True)
        stem = output_stem or f"{alloy}_{parameter}_three_levels_1x3_side_colors"
        for extension in ("pdf", "svg", "tiff", "png"):
            path = output / f"{stem}.{extension}"
            kwargs = {"pil_kwargs": {"compression": "tiff_lzw"}} if extension == "tiff" else {}
            fig.savefig(path, dpi=600, facecolor="white", **kwargs)
            if extension == "tiff":
                with Image.open(path) as source:
                    image = source.convert("RGB")
                image.save(path, dpi=(600, 600), compression="tiff_lzw")
        metadata = {"alloy": alloy, "parameter": parameter,
                    "size_cm": (fig.get_size_inches() * 2.54).tolist(), "stem": stem}
        if include_aqueous:
            metadata['legend_counts'] = [len(key.texts) for _, key in pairs]
            metadata['legend_labels'] = [[text.get_text() for text in key.texts] for _, key in pairs]
        plt.close(fig)
    print(f"Saved {output / stem}", flush=True)
    return metadata


def write_csv(path: Path, records: list) -> None:
    """Write rectangular result tables with explicit headers."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    """Calculate the specified sweep, export figures, and record provenance."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-size", type=int, default=1000)
    parser.add_argument("--alloys", nargs="+", choices=[s[2] for s in collage.SYSTEMS],
                        default=[s[2] for s in collage.SYSTEMS])
    parser.add_argument("--output-dir", type=Path, default=FIGURE_ROOT)
    parser.add_argument("--results-dir", type=Path, default=RESULT_ROOT)
    args = parser.parse_args()
    reference = verify_reference()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    records, figures, palettes = [], [], {}
    for metal1, metal2, alloy in collage.SYSTEMS:
        if alloy not in args.alloys:
            continue
        data = SpeciesDataLoader(metal1, metal2, collage.MU_LIGAND, collage.TEMPERATURE_K,
                                 data_dir=collage.DATA_DIR)
        accessories = PlotAccessories(data)
        cache = {}
        conditions = [BASELINE]
        for parameter, values in LEVELS.items():
            conditions.extend([{**BASELINE, parameter: value} for value in values])
        for condition in conditions:
            key = condition_key(condition)
            if key in cache:
                continue
            print(f"Calculating {alloy} {condition}", flush=True)
            cache[key] = calculate_case(data, condition, args.grid_size)
            case = cache[key]
            audit_dir = args.results_dir / alloy
            audit_dir.mkdir(exist_ok=True)
            case_id = f"case_{len(cache):02d}"
            case["case_id"] = case_id
            np.savez_compressed(audit_dir / f"{case_id}.npz",
                                stable_indices=case["regions"].min_indices,
                                pH=case["grid"].pH_values, E_SHE=case["grid"].V_values)
            (audit_dir / f"{case_id}.json").write_text(json.dumps({
                "condition": condition,
                "species_list": [case["regions"].format_species_tuple(s)
                                 for s in case["regions"].species_list]}, indent=2))
        colors = consistent_colors(list(cache.values()), accessories)
        palettes[alloy] = [{"combination": combo, "display_color": to_hex(
            tuple(REGION_COLOR_STRENGTH * c + 1 - REGION_COLOR_STRENGTH for c in to_rgba(color)[:3]))}
            for combo, color in sorted(colors.items())]
        baseline = cache[condition_key(BASELINE)]
        for parameter, values in LEVELS.items():
            cases = [cache[condition_key({**BASELINE, parameter: value})] for value in values]
            for index, case in enumerate(cases):
                records.append({"alloy": alloy, "parameter": parameter, "panel": "abc"[index],
                                "case_id": case["case_id"], **case["condition"],
                                **summarize(case, baseline, accessories)})
            figures.append(render(alloy, parameter, cases, accessories, colors, args.output_dir / alloy))
        write_csv(args.results_dir / "sensitivity_summary.csv", records)
    hashes = {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in sorted(collage.DATA_DIR.glob("*.json"))}
    manifest = {"reference": reference, "baseline": BASELINE, "levels": LEVELS,
                "temperature_K": collage.TEMPERATURE_K, "grid_size": args.grid_size,
                "pH_range": collage.PH_RANGE, "E_SHE_range": collage.V_RANGE,
                "experimental_pH": collage.PH_EXP_RANGE, "experimental_E_RHE": collage.V_EXP_RANGE,
                "mu_ligand_eV": collage.MU_LIGAND, "data_sha256": hashes,
                "figures": figures, "palettes": palettes}
    (args.results_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Finished: {len(figures)} triptychs; {len(records)} comparison records", flush=True)


if __name__ == "__main__":
    main()
