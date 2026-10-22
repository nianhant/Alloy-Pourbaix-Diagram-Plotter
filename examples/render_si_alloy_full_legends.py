"""Render existing SI grids with complete legends in a separate output set."""
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from examples import generate_si_alloy_panels as si
from examples import generate_alloy_sensitivity_triptychs as sweep
from examples import make_rsc_pourbaix_collage as collage
from methods.data_loader import SpeciesDataLoader
from methods.grid import GridMaker
from methods.thermodynamics import PourbaixData, SpeciesGridResult
from methods.visualization import PlotAccessories

OUTPUT = ROOT / 'figures/alloy_pourbaix_diagrams/si_full_legends'
RESULTS = ROOT / 'research/2026-10-02_si_alloy_full_legends'


def render_set(source_results: Path, output: Path, results: Path) -> None:
    """Reuse exact stored grids, preserving concentrations and all region IDs."""
    manifest = json.loads((source_results / 'manifest.json').read_text())
    groups = {g['alloy']: g for g in json.loads((source_results / 'source_inventory.json').read_text())}
    for path, digest in manifest['data_sha256'].items():
        assert si.sha256(ROOT / path) == digest, f'Dataset changed: {path}'
    metadata, captions = [], []
    for original in manifest['figures']:
        alloy = original['alloy']
        group = groups[alloy]
        data = SpeciesDataLoader(*group['metals'], collage.MU_LIGAND, collage.TEMPERATURE_K,
                                 data_dir=collage.DATA_DIR)
        accessories = PlotAccessories(data)
        cases = []
        for index, settings in enumerate(original['conditions']):
            stem = source_results / alloy / f'panel_{chr(97+index)}'
            mapping = json.loads(stem.with_suffix('.json').read_text())
            assert mapping['condition'] == settings
            ligands = {name: settings[name] for name in ('NH3', 'Gly', 'CN')}
            ligands['NO2'] = 0.0
            grid = GridMaker(collage.PH_RANGE, collage.V_RANGE, ligands, manifest['grid_size'])
            thermo = PourbaixData(data, settings['activity'], ligands,
                                  collage.equiatomic_reference_composition(data))
            species_list = list(thermo.reaction_coefficients)
            assert [list(SpeciesGridResult.format_species_tuple(s)) for s in species_list] == mapping['species_list']
            with np.load(stem.with_suffix('.npz')) as saved:
                indices = saved['stable_indices'].copy()
                assert np.array_equal(grid.pH_values, saved['pH'])
                assert np.array_equal(grid.V_values, saved['E_SHE'])
            regions = SpeciesGridResult(indices, species_list)
            stable = {regions.format_species_tuple(species_list[i]) for i in np.unique(indices)}
            cases.append({'condition': settings, 'grid': grid, 'thermo': thermo,
                          'regions': regions, 'stable': stable})
        colors = sweep.consistent_colors(cases, accessories)
        stem = original['stem'] + '_all_species'
        print(f'Rendering {alloy}: complete legends from saved grids', flush=True)
        result = sweep.render(alloy, 'CN', cases, accessories, colors, output,
                              titles=[''] * len(cases), output_stem=stem, include_aqueous=True)
        result.update({key: original[key] for key in ('conditions', 'resolved_modes', 'source_images', 'figure_label')})
        metadata.append(result)
        captions.append(si.latex_group(group, stem, original['resolved_modes']))
    output.mkdir(parents=True, exist_ok=True)
    (output / 'alloy_figures_SI_replacement.tex').write_text(
        '% Separate complete-legend versions; previous artwork is untouched.\n\n' + '\n\n'.join(captions) + '\n')
    results.mkdir(parents=True, exist_ok=True)
    manifest.update({'source_grid_manifest': str(source_results / 'manifest.json'),
                     'figures': metadata, 'legend_order': ['solid oxides/oxyhydroxides',
                     'metal-ligand complexes', 'metals/hydrides', 'aqueous ions'],
                     'legend_row_spacing_pt': 11, 'reuses_saved_grids': True})
    (results / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    paths = [output / (f['stem'] + '.png') for f in metadata]
    for offset in range(0, len(paths), 2):
        prefix = 'alternatives' if output.name == 'caption_alternatives' else 'primary'
        si.review_sheet(paths[offset:offset+2], ROOT / '.agents/test/si_full_legends' / f'{prefix}_{offset//2+1}.png')


def main() -> None:
    """Export both existing SI interpretations and verify originals unchanged."""
    protected = {str(p): si.sha256(p) for p in si.OUTPUT.rglob('*') if p.is_file()}
    protected.update({str(p): si.sha256(p) for p in si.RESULTS.rglob('*') if p.is_file()})
    render_set(si.RESULTS, OUTPUT, RESULTS)
    render_set(si.RESULTS / 'caption_alternatives', OUTPUT / 'caption_alternatives',
               RESULTS / 'caption_alternatives')
    assert all(si.sha256(Path(path)) == digest for path, digest in protected.items())
    (RESULTS / 'preserved_originals.json').write_text(json.dumps(protected, indent=2))
    print(f'PASS: {len(protected)} original files unchanged; all legend regions included.', flush=True)


if __name__ == '__main__':
    main()
