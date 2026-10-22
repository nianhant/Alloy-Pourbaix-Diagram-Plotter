"""Apply the reviewed, checksum-guarded SI alloy artwork update."""
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
STAGING = ROOT / 'research/2026-10-02_si_alloy_submission_sync'
plan = json.loads((STAGING / 'sync_plan.json').read_text())
target = Path(plan['source'])
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
assert digest(target) == plan['source_sha256'], 'LaTeX changed since preparation; rebase before applying.'
assert digest(STAGING / 'pourbaix_figures_SI.before.tex') == plan['source_sha256']
for item in plan['files']:
    assert digest(Path(item['source'])) == item['sha256']
    assert not Path(item['destination']).exists(), item['destination']
# Snapshot every pre-existing submission file. Only the exact LaTeX target may change.
protected = {p: digest(p) for p in target.parent.iterdir() if p.is_file() and p != target}
backup = target.with_name('pourbaix_figures_SI.pre_alloy_full_legends.tex.bak')
assert not backup.exists(), 'Do not overwrite an earlier backup.'
shutil.copy2(target, backup)
for item in plan['files']:
    shutil.copy2(item['source'], item['destination'])
assert digest(target) == plan['source_sha256'], 'Concurrent LaTeX edit detected; artwork copied but TeX not replaced.'
shutil.copy2(STAGING / 'pourbaix_figures_SI.after.tex', target)
assert digest(target) == digest(STAGING / 'pourbaix_figures_SI.after.tex')
assert digest(backup) == plan['source_sha256']
assert all(digest(path) == checksum for path, checksum in protected.items())
for item in plan['files']:
    assert digest(Path(item['destination'])) == item['sha256']
    assert digest(Path(item['source'])) == item['sha256']
active = '\n'.join(re.split(r'(?<!\\)%', line)[0] for line in target.read_text().splitlines())
names = re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+_all_species\.pdf)\}', active)
assert len(names) == 11 and all((target.parent / name).is_file() for name in names)
assert not re.search(r'pbx_(?:Ni_Ti|Cu_Ti|Ti_Fe|Ti_Mg|Ti_Mn|Ti_Zn|Pd_Ag|Pd_Cu|Pd_Ni|Pd_Pt|Au_Pd)_', active)
print(json.dumps({'copied_artwork_files': 22, 'updated_alloy_figures': 11,
                  'backup': str(backup), 'updated_tex': str(target),
                  'labels_preserved': plan['labels_preserved'],
                  'non_alloy_content_preserved': plan['non_alloy_content_preserved'],
                  'preexisting_submission_files_unchanged': len(protected),
                  'local_artwork_preserved': True, 'latex_compile_available': False}, indent=2))
