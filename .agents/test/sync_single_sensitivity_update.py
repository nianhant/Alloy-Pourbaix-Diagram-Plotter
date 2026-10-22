"""Apply the user-approved removal of CN sweeps and concise range captions."""
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
STAGE = ROOT / 'research/2026-10-02_si_alloy_sensitivity_sync'
plan = json.loads((STAGE / 'single_update_plan.json').read_text())
target = Path(plan['source'])
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
assert digest(target) == plan['source_sha256'], 'Concurrent modification; rebase before copying.'
backup = target.with_name('pourbaix_sensitivity_SI.pre_range_only.tex.bak')
assert not backup.exists()
protected = {p: digest(p) for p in target.parent.iterdir() if p.is_file() and p != target}
shutil.copy2(target, backup)
assert digest(target) == plan['source_sha256']
shutil.copy2(STAGE / 'single_sensitivity.after.tex', target)
assert digest(target) == digest(STAGE / 'single_sensitivity.after.tex')
assert digest(backup) == plan['source_sha256']
assert all(digest(p) == checksum for p, checksum in protected.items())
all_labels = []
for name, count in [('pourbaix_sensitivity_SI.tex', 18), ('alloy_pourbaix_sensitivity_SI.tex', 9)]:
    text = (target.parent / name).read_text()
    assert text.count(r'\begin{figure}') == count
    assert 'NO2' not in text and '_CN_' not in text and 'sweep' not in text
    figures = re.findall(r'\\includegraphics\[[^\]]*\]\{([^}]+)\}', text)
    assert len(figures) == count and all((target.parent / f).is_file() for f in figures)
    captions = re.findall(r'\\caption\{(.+)', text)
    assert len(captions) == count and all('varied from' in c for c in captions)
    all_labels.extend(re.findall(r'\\label\{([^}]+)\}', text))
assert len(all_labels) == len(set(all_labels)) == 27
main = (target.parent / 'SI.tex').read_text()
assert main.count(r'\input{pourbaix_sensitivity_SI}') == 1
assert main.count(r'\input{alloy_pourbaix_sensitivity_SI}') == 1
print(json.dumps({'single_metal_figures': 18, 'alloy_figures': 9,
                  'CN_sweep_blocks_removed': 6, 'range_only_captions': 27,
                  'NO2_mentions_in_sensitivity_sections': 0,
                  'all_graphics_exist': True, 'labels_unique': True,
                  'unchanged_other_files': len(protected), 'backup': str(backup),
                  'CN_artwork_retained_on_disk': True}, indent=2))
