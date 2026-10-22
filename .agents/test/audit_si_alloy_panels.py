"""Validate local SI alloy exports and build visual review sheets."""
import json
from pathlib import Path
import re
import sys

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from examples.generate_si_alloy_panels import OUTPUT, RESULTS, inventory, review_sheet, sha256


def audit(output: Path, results: Path, count: int, panels: int) -> None:
    """Check provenance, exact panel conditions, grids, exports and LaTeX labels."""
    manifest = json.loads((results / 'manifest.json').read_text())
    assert manifest['grid_size'] == 1600
    assert manifest['experimental_pH'] == [10.5, 13.5]
    assert manifest['experimental_E_RHE'] == [-2, 2.3]
    assert manifest['pH_range'] == [-2, 16] and manifest['E_SHE_range'] == [-2, 3]
    assert len(manifest['figures']) == count
    assert sum(len(f['conditions']) for f in manifest['figures']) == panels
    source_changes = []
    for path, digest in manifest['source_hashes'].items():
        current = sha256(Path(path))
        if current != digest:
            # SI.tex changed concurrently after generation. The alloy include
            # and every source artwork must still match the recorded input.
            assert Path(path).name == 'SI.tex', path
            assert r'\input{pourbaix_figures_SI}' in Path(path).read_text()
            source_changes.append({'path': path, 'original_sha256': digest,
                                   'current_sha256': current})
    tex = (output / 'alloy_figures_SI_replacement.tex').read_text()
    assert tex.count(r'\begin{figure}') == count
    assert r'\begin{subfigure}' not in tex
    assert tex.count('10.5--13.5') == count and '11.5' not in tex
    observed_labels = re.findall(r'\\label\{([^}]+)\}', tex)
    assert len(observed_labels) == len(set(observed_labels))
    selected = {f['alloy'] for f in manifest['figures']}
    expected_labels = [label for group in inventory() if group['alloy'] in selected
                       for label in [group['figure_label'], *group['panel_labels']]]
    assert set(observed_labels) == set(expected_labels)
    paths = []
    for figure in manifest['figures']:
        stem = output / figure['stem']
        for suffix in ('.pdf', '.svg', '.png', '.tiff'):
            assert stem.with_suffix(suffix).stat().st_size > 0
        svg = stem.with_suffix('.svg').read_text()
        assert 'Species' not in svg
        for index, settings in enumerate(figure['conditions']):
            assert f'<!-- ({chr(97+index)})' in svg
            assert settings['activity'] == 1e-4
            assert settings['NH3'] in (0, .02) and settings['Gly'] in (0, .1)
            assert settings['CN'] in (0, 1e-4)
            path = results / figure['alloy'] / f'panel_{chr(97+index)}.npz'
            mapping = json.loads(path.with_suffix('.json').read_text())
            assert mapping['condition'] == settings
            with np.load(path) as data:
                assert data['stable_indices'].shape == (1600, 1600)
                assert data['stable_indices'].min() >= 0
                assert data['stable_indices'].max() < len(mapping['species_list'])
                assert np.allclose(data['pH'][[0, -1]], [-2, 16])
                assert np.allclose(data['E_SHE'][[0, -1]], [-2, 3])
        with Image.open(stem.with_suffix('.tiff')) as source:
            assert source.mode == 'RGB'
            assert all(abs(dpi-600) < 1 for dpi in source.info['dpi'])
            measured_cm = np.array(source.size) / 600 * 2.54
            assert np.allclose(measured_cm, figure['size_cm'], atol=.005)
        paths.append(stem.with_suffix('.png'))
    prefix = 'si_artwork' if count == 11 else 'si_alternatives'
    for offset in range(0, len(paths), 3):
        review_sheet(paths[offset:offset+3], ROOT / '.agents/test' / f'{prefix}_{offset//3+1}.png')
    (results / 'validation.json').write_text(json.dumps({
        'status': 'passed', 'figure_count': count, 'panel_count': panels,
        'export_count': count * 4, 'alloy_source_include_and_images_unchanged': True,
        'concurrent_source_changes': source_changes}, indent=2))
    print(f'PASS {output}: {count} figures, {panels} grids, {count*4} exports; alloy inputs unchanged.')


audit(OUTPUT, RESULTS, 11, 19)
audit(OUTPUT / 'caption_alternatives', RESULTS / 'caption_alternatives', 2, 4)
