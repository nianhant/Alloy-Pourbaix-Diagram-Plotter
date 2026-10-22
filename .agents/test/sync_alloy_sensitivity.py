"""Copy the nine non-CN sensitivity triptychs and add their SI section."""
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
STAGE = ROOT / 'research/2026-10-02_si_alloy_sensitivity_sync'
plan = json.loads((STAGE / 'sync_plan.json').read_text())
target = Path(plan['source'])
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
assert digest(target) == plan['source_sha256'], 'Concurrent SI.tex edit: rebase first.'
assert digest(STAGE / 'SI.before.tex') == plan['source_sha256']
section = Path(plan['section_destination'])
backup = target.with_name('SI.pre_alloy_sensitivity.tex.bak')
assert not section.exists() and not backup.exists()
for item in plan['files']:
    assert not Path(item['destination']).exists()
    assert digest(Path(item['source'])) == item['sha256']
protected = {p: digest(p) for p in target.parent.iterdir() if p.is_file() and p != target}
shutil.copy2(target, backup)
for item in plan['files']:
    shutil.copy2(item['source'], item['destination'])
shutil.copy2(STAGE / 'alloy_pourbaix_sensitivity_SI.tex', section)
assert digest(target) == plan['source_sha256'], 'Concurrent SI edit; new files copied, main SI not replaced.'
shutil.copy2(STAGE / 'SI.after.tex', target)
assert digest(target) == digest(STAGE / 'SI.after.tex')
assert all(digest(path) == checksum for path, checksum in protected.items())
assert digest(backup) == plan['source_sha256']
for item in plan['files']:
    assert digest(Path(item['destination'])) == item['sha256']
    assert digest(Path(item['source'])) == item['sha256']
text = section.read_text()
names = re.findall(r'\\includegraphics\[[^\]]*\]\{([^}]+)\}', text)
assert len(names) == 9 and all((target.parent / name).is_file() for name in names)
assert 'NO2' not in text and '_CN_' not in text and 'sweep' not in text
labels = re.findall(r'\\label\{([^}]+)\}', text)
other = '\n'.join(p.read_text() for p in target.parent.glob('*.tex') if p != section)
assert all(f'\\label{{{label}}}' not in other for label in labels)
assert target.read_text().count(r'\input{alloy_pourbaix_sensitivity_SI}') == 1
print(json.dumps({'added_figures': 9, 'copied_PDF_TIFF_files': 18,
                  'updated_SI': str(target), 'section': str(section),
                  'backup': str(backup), 'other_submission_files_preserved': len(protected),
                  'local_figures_preserved': True, 'references_verified': True}, indent=2))
