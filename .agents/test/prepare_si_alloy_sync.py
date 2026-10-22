"""Read-only preparation of a narrowly scoped SI alloy figure update."""
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
SUBMISSION = Path('/global/homes/n/nianhant/data/stability_paper/updated_stability_manuscript/submission')
ARTWORK = ROOT / 'figures/alloy_pourbaix_diagrams/si_full_legends'
source = SUBMISSION / 'pourbaix_figures_SI.tex'
before = source.read_text()
draft = (ARTWORK / 'alloy_figures_SI_replacement.tex').read_text()
manifest = json.loads((ROOT / 'research/2026-10-02_si_alloy_full_legends/manifest.json').read_text())
labels = {f['figure_label'] for f in manifest['figures']}
pattern = re.compile(r'^\\begin\{figure\}.*?^\\end\{figure\}', re.M | re.S)
replacements = {}
for match in pattern.finditer(draft):
    block = match.group()
    label = next(label for label in labels if f'\\label{{{label}}}' in block)
    figure = next(f for f in manifest['figures'] if f['figure_label'] == label)
    if len(figure['conditions']) == 1:
        block = block.replace('Pourbaix diagrams of', 'Pourbaix diagram of')
    if all(c['CN'] == 0 for c in figure['conditions']):
        block = block.replace(
            r'; cyanide-containing panels additionally use $[\ce{CN-}]_{\mathrm{tot}}=10^{-4}\,\mathrm{M}$.',
            r', with no added cyanide.')
    block = block.replace(
        '(converted to SHE and clipped to the plotted range).}',
        '(converted to SHE and clipped to the plotted range). '
        'The side legends list all stable species combinations, including solids, '
        'aqueous ions, and metal--ligand complexes. The dashed blue and red lines '
        'mark the water reduction and oxidation boundaries, respectively.}')
    replacements[label] = block
assert len(replacements) == 11
seen = []
def replace(match):
    block = match.group()
    found = [label for label in labels if f'\\label{{{label}}}' in block]
    if not found:
        return block
    assert len(found) == 1
    seen.append(found[0])
    return replacements[found[0]]
after = pattern.sub(replace, before)
assert set(seen) == labels and len(seen) == 11
def active_labels(text):
    active = '\n'.join(re.split(r'(?<!\\)%', line)[0] for line in text.splitlines())
    return sorted(re.findall(r'\\label\{([^}]+)\}', active))
assert active_labels(before) == active_labels(after)
# Everything outside the selected alloy figure blocks remains byte-identical.
def mask(text):
    return pattern.sub(lambda m: '<ALLOY>' if any(f'\\label{{{l}}}' in m.group() for l in labels) else m.group(), text)
assert mask(before) == mask(after)
files = []
for figure in manifest['figures']:
    for suffix in ('.pdf', '.tiff'):
        path = ARTWORK / (figure['stem'] + suffix)
        destination = SUBMISSION / path.name
        assert path.is_file()
        files.append({'source': str(path), 'destination': str(destination),
                      'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                      'destination_exists': destination.exists()})
assert not any(f['destination_exists'] for f in files), 'Inspect existing targets before overwriting.'
print(json.dumps({'before': before, 'after': after, 'source': str(source),
                  'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                  'files': files, 'labels_preserved': True,
                  'non_alloy_content_preserved': True,
                  'latex_tools': {name: shutil.which(name) for name in ('pdflatex', 'latexmk', 'tectonic')}}))
