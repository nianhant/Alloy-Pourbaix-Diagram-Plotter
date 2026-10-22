"""Verify unchanged saved phase maps and apply the reviewed SI omission."""
import csv
import json
from pathlib import Path
import re
import shutil
import numpy as np

root = Path(__file__).resolve().parents[2]
data = root / 'research/2026-10-02_alloy_concentration_sensitivity'
stage = root / '.agents/test/sensitivity_pruning'
submission = Path('/global/homes/n/nianhant/data/stability_paper/updated_stability_manuscript/submission')
target = submission / 'alloy_pourbaix_sensitivity_SI.tex'
rows = list(csv.DictReader((data / 'sensitivity_summary.csv').open()))
omitted = []
for metal in ('TiNi', 'TiCu', 'PdAu'):
    for parameter in ('activity', 'NH3', 'Gly'):
        selected = [r for r in rows if r['alloy'] == metal and r['parameter'] == parameter]
        assert len(selected) == 3
        base_path = data / metal / 'case_01.npz'
        base_labels = json.loads(base_path.with_suffix('.json').read_text())['species_list']
        with np.load(base_path) as archive:
            base = archive['stable_indices'].copy()
        counts = []
        for row in selected:
            path = data / metal / (row['case_id'] + '.npz')
            assert json.loads(path.with_suffix('.json').read_text())['species_list'] == base_labels
            with np.load(path) as archive:
                counts.append(int(np.count_nonzero(archive['stable_indices'] != base)))
        assert all(abs(n / base.size - float(r['full_changed_fraction'])) < 1e-12 for n, r in zip(counts, selected))
        print(metal, parameter, 'changed grid cells:', counts)
        if not any(counts):
            omitted.append((metal, parameter))
assert omitted == [('TiNi', 'NH3'), ('PdAu', 'Gly')]
before = (stage / 'binary.before.tex').read_text()
after = (stage / 'binary.after.tex').read_text()
assert target.read_text().rstrip() == before.rstrip(), 'Concurrent edit; abort'
assert after.count(r'\begin{figure}') == after.count(r'\end{figure}') == 7
for name in re.findall(r'\\includegraphics\[[^\]]*\]\{([^}]+)\}', after):
    assert (submission / name).is_file()
for metal, parameter in omitted:
    label = f'fig:{metal}_{parameter}_sensitivity_SI'
    for path in submission.glob('*.tex'):
        text = re.sub(r'(?<!\\)%[^\n]*', '', path.read_text())
        assert not any(label in group for group in re.findall(r'\\(?:[Cc]ref|ref|autoref)\{([^}]+)\}', text)), f'Referenced omitted figure: {path}'
backup = stage / 'binary.original.tex'
assert not backup.exists()
shutil.copy2(target, backup)
assert target.read_text().rstrip() == before.rstrip()
shutil.copy2(stage / 'binary.after.tex', target)
assert target.read_bytes() == (stage / 'binary.after.tex').read_bytes()
print('PASS: two identical triptychs omitted; seven retained; original source and all artwork preserved.')
