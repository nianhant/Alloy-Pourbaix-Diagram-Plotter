"""Apply the requested SI terminology change with backup and label checks."""
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
STAGE = ROOT / 'research/2026-10-02_si_terminology'
plan = json.loads((STAGE / 'plan.json').read_text())
for item in plan:
    target = Path(item['path'])
    before = (STAGE / (item['name'] + '.before')).read_bytes()
    after = (STAGE / (item['name'] + '.after')).read_bytes()
    assert target.read_bytes() == before, 'Concurrent changes: rebase first.'
    backup = target.with_name(target.name + '.pre_binary_metal_phases.bak')
    assert not backup.exists()
    identifiers = rb'\\(?:label|input|includegraphics|[cC]ref|ref)(?:\[[^\]]*\])?\{[^}]*\}'
    assert re.findall(identifiers, before) == re.findall(identifiers, after)
for item in plan:
    target = Path(item['path'])
    assert target.read_bytes() == (STAGE / (item['name'] + '.before')).read_bytes()
    shutil.copy2(target, target.with_name(target.name + '.pre_binary_metal_phases.bak'))
    shutil.copy2(STAGE / (item['name'] + '.after'), target)
    assert target.read_bytes() == (STAGE / (item['name'] + '.after')).read_bytes()
todo = [Path(plan[0]['path']).parent / 'SI.tex']
seen = set()
while todo:
    path = todo.pop()
    if path in seen:
        continue
    seen.add(path)
    text = path.read_text()
    assert re.search(r'\balloys?\b', text, re.I) is None, path
    for name in re.findall(r'\\input\{([^}]+)\}', text):
        child = path.parent / name
        todo.append(child if child.suffix else child.with_suffix('.tex'))
print(f'Updated {len(plan)} files; verified {len(seen)} SI source files. Filenames and labels preserved; backups created.')
