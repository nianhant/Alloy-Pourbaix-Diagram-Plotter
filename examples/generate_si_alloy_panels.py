"""Regenerate only active SI alloy Pourbaix figures, entirely within this repo.

Parse the SI include file, preserve its figure/panel ordering, and replace
LaTeX subfigure assembly with the shared Matplotlib panel renderer. No writes
are made to the source submission. Conflicting H2O filenames/captions require
an explicit --ambiguous-baseline selection for a full run.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from examples import generate_alloy_sensitivity_triptychs as sweep
from examples import make_rsc_pourbaix_collage as collage
from methods.data_loader import SpeciesDataLoader
from methods.visualization import PlotAccessories

SUBMISSION = Path('/global/homes/n/nianhant/data/stability_paper/updated_stability_manuscript/submission')
SOURCE = SUBMISSION / 'pourbaix_figures_SI.tex'
OUTPUT = ROOT / 'figures/alloy_pourbaix_diagrams/si_regenerated'
RESULTS = ROOT / 'research/2026-10-02_si_alloy_regeneration'
ALLOY_NAMES = {('Ni', 'Ti'): 'TiNi', ('Cu', 'Ti'): 'TiCu', ('Au', 'Pd'): 'AuPd'}


def sha256(path: Path) -> str:
    """Return the SHA-256 of an existing provenance input."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory() -> list:
    """Read active alloy figure groups, ignoring commented-out LaTeX."""
    active = '\n'.join(re.split(r'(?<!\\)%', line, maxsplit=1)[0]
                       for line in SOURCE.read_text().splitlines())
    groups = []
    pattern = r'pbx_([A-Z][a-z]?)_([A-Z][a-z]?)_(H2O|NH3_Gly|NH3_Gly_CN)\.png'
    for block in re.findall(r'\\begin\{figure\}.*?\\end\{figure\}', active, re.S):
        files = re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}', block)
        if not files or not all(re.fullmatch(pattern, name) for name in files):
            continue
        matches = [re.fullmatch(pattern, name) for name in files]
        metals = matches[0].group(1, 2)
        assert all(m.group(1, 2) == metals for m in matches)
        alloy = ALLOY_NAMES.get(metals, ''.join(metals))
        labels = re.findall(r'\\label\{([^}]+)\}', block)
        groups.append({'alloy': alloy, 'metals': metals, 'sources': files,
                       'source_modes': [m.group(3) for m in matches],
                       'figure_label': labels[-1], 'panel_labels': labels[:-1],
                       'original_latex': block,
                       'ambiguous': any(m.group(3) == 'H2O' for m in matches)})
    assert len(groups) == 11 and sum(len(g['sources']) for g in groups) == 19
    return groups


def condition(mode: str) -> dict:
    """Use the SI ligand totals and the shared 1e-4 ion activity."""
    return {'activity': 1e-4, 'NH3': 0.0 if mode == 'H2O' else .02,
            'Gly': 0.0 if mode == 'H2O' else .1,
            'CN': 1e-4 if mode == 'NH3_Gly_CN' else 0.0}


def review_sheet(paths: list, output: Path, width: int = 1800) -> None:
    """Make a labeled contact sheet without changing any figure artwork."""
    previews = []
    for path in paths:
        with Image.open(path) as source:
            preview = source.convert('RGB')
        preview.thumbnail((width, 1000))
        previews.append((path.stem, preview))
    sheet = Image.new('RGB', (width + 20, sum(p.height + 30 for _, p in previews)), 'white')
    draw = ImageDraw.Draw(sheet)
    y = 0
    for label, preview in previews:
        draw.text((10, y + 4), label, fill='black')
        sheet.paste(preview, (10, y + 25))
        y += preview.height + 30
    output.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(output)


def latex_group(group: dict, stem: str, modes: list) -> str:
    """Draft a local replacement figure preserving existing cross-references."""
    descriptions = {'H2O': 'No added ligands', 'NH3_Gly': r'\ce{NH3} and Gly',
                    'NH3_Gly_CN': r'\ce{NH3}, Gly, and \ce{CN-}'}
    parts = [f'({chr(97+i)}) {descriptions[mode]}' for i, mode in enumerate(modes)]
    caption = (f"Pourbaix diagrams of {group['alloy']}. " + '; '.join(parts) + '. '
               r'The metal-ion activity is $10^{-4}$. Ligand-containing panels use '
               r'$[\ce{NH3}]_{\mathrm{tot}}=0.02\,\mathrm{M}$ and '
               r'$[\mathrm{Gly}]_{\mathrm{tot}}=0.1\,\mathrm{M}$; cyanide-containing '
               r'panels additionally use $[\ce{CN-}]_{\mathrm{tot}}=10^{-4}\,\mathrm{M}$. '
               r'The black outline marks pH 10.5--13.5 and $-2.0$ to $2.0$ V vs RHE '
               r'(converted to SHE and clipped to the plotted range).')
    width = r'\textwidth' if len(modes) > 1 else r'0.6\textwidth'
    lines = [r'\begin{figure}[!htbp]', r'    \centering',
             rf'    \includegraphics[width={width}]{{{stem}.pdf}}',
             rf'    \caption{{{caption}}}', rf"    \label{{{group['figure_label']}}}"]
    lines.extend(rf'    \phantomsubcaption\label{{{label}}}' for label in group['panel_labels'])
    return '\n'.join(lines + [r'\end{figure}'])


def main() -> None:
    """Inventory source artwork, calculate selected cases, and export locally."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alloys', nargs='+')
    parser.add_argument('--ambiguous-baseline', choices=('artwork', 'captions'))
    parser.add_argument('--inspect-only', action='store_true')
    parser.add_argument('--grid-size', type=int, default=1000)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    parser.add_argument('--results-dir', type=Path, default=RESULTS)
    args = parser.parse_args()
    groups = inventory()
    selected = [g for g in groups if not args.alloys or g['alloy'] in args.alloys]
    source_hashes = {str(SOURCE): sha256(SOURCE), str(SUBMISSION / 'SI.tex'): sha256(SUBMISSION / 'SI.tex')}
    for group in groups:
        for name in group['sources']:
            source_hashes[str(SUBMISSION / name)] = sha256(SUBMISSION / name)
    args.results_dir.mkdir(parents=True, exist_ok=True)
    (args.results_dir / 'source_inventory.json').write_text(json.dumps(groups, indent=2))
    if args.inspect_only:
        for group in groups:
            print(group['alloy'], group['source_modes'], group['figure_label'])
            review_sheet([SUBMISSION / name for name in group['sources']],
                         ROOT / '.agents/test/si_originals' / f"{group['alloy']}.png")
        return
    if any(g['ambiguous'] for g in selected) and args.ambiguous_baseline is None:
        parser.error('Select --ambiguous-baseline artwork or captions for TiCu/AuPd.')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata, snippets = [], []
    for group in selected:
        alloy = group['alloy']
        metal1, metal2 = group['metals']
        # Refuse a network/data refresh: regeneration must use local datasets.
        assert any(all((collage.DATA_DIR / f'{a}_{b}_{kind}_formation_energy.json').exists()
                       for kind in ('solid', 'ion')) for a, b in ((metal1, metal2), (metal2, metal1)))
        data = SpeciesDataLoader(metal1, metal2, collage.MU_LIGAND, collage.TEMPERATURE_K,
                                 data_dir=collage.DATA_DIR)
        accessories = PlotAccessories(data)
        modes = ['NH3_Gly' if mode == 'H2O' and args.ambiguous_baseline == 'captions' else mode
                 for mode in group['source_modes']]
        cases = []
        case_dir = args.results_dir / alloy
        case_dir.mkdir(exist_ok=True)
        for index, mode in enumerate(modes):
            settings = condition(mode)
            print(f'Calculating {alloy} ({chr(97+index)}): {settings}', flush=True)
            case = sweep.calculate_case(data, settings, args.grid_size)
            cases.append(case)
            stem = f'panel_{chr(97+index)}'
            np.savez_compressed(case_dir / f'{stem}.npz',
                                stable_indices=case['regions'].min_indices,
                                pH=case['grid'].pH_values, E_SHE=case['grid'].V_values)
            (case_dir / f'{stem}.json').write_text(json.dumps({
                'condition': settings, 'mode': mode,
                'source_image': group['sources'][index],
                'species_list': [case['regions'].format_species_tuple(s)
                                 for s in case['regions'].species_list]}, indent=2))
        colors = sweep.consistent_colors(cases, accessories)
        stem = f'SI_{alloy}_{len(cases)}panel_side_colors'
        # Panel letters alone: conditions belong in the SI caption, not titles.
        result = sweep.render(alloy, 'CN', cases, accessories, colors, args.output_dir,
                              titles=[''] * len(cases), output_stem=stem)
        result.update({'conditions': [case['condition'] for case in cases],
                       'source_images': group['sources'], 'source_modes': group['source_modes'],
                       'resolved_modes': modes, 'figure_label': group['figure_label']})
        metadata.append(result)
        snippets.append(latex_group(group, stem, modes))
    # A local, alloy-only draft; never modify submission TeX or non-alloy figures.
    (args.output_dir / 'alloy_figures_SI_replacement.tex').write_text(
        '% Local review draft only; source submission has not been edited.\n\n' + '\n\n'.join(snippets) + '\n')
    for path, digest in source_hashes.items():
        assert sha256(Path(path)) == digest, f'Submission input changed: {path}'
    manifest = {'source_hashes': source_hashes, 'grid_size': args.grid_size,
                'temperature_K': collage.TEMPERATURE_K, 'mu_ligand_eV': collage.MU_LIGAND,
                'pH_range': collage.PH_RANGE, 'E_SHE_range': collage.V_RANGE,
                'experimental_pH': collage.PH_EXP_RANGE, 'experimental_E_RHE': collage.V_EXP_RANGE,
                'ambiguous_baseline_choice': args.ambiguous_baseline, 'figures': metadata,
                'data_sha256': {str(p.relative_to(ROOT)): sha256(p) for p in collage.DATA_DIR.glob('*.json')}}
    (args.results_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(f'Finished {len(metadata)} figures; submission unchanged.', flush=True)


if __name__ == '__main__':
    main()
