"""Prepare range-only single-metal sensitivity captions without CN sweeps."""
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
SUBMISSION = Path('/global/homes/n/nianhant/data/stability_paper/updated_stability_manuscript/submission')
source = SUBMISSION / 'pourbaix_sensitivity_SI.tex'
before = source.read_text()
intro = (ROOT / 'research/2026-10-02_si_alloy_sensitivity_sync/alloy_pourbaix_sensitivity_SI.tex').read_text().split('\n\n')[0]
intro = intro.replace('equiatomic TiNi, TiCu, and PdAu', 'Au, Ni, Cu, Pd, Pt, and Ti')
assert before.count('% Au sensitivity panels') == 1
body = before[before.index('% Au sensitivity panels'):]
specs = {'activity': ('nominal aqueous metal-ion concentration', r'$10^{-6}$ to $10^{-2}\,\mathrm{M}$'),
         'NH3': ('total ammonia concentration', r'$0.002$ to $0.1\,\mathrm{M}$'),
         'Gly': ('total glycine concentration', r'$0.01$ to $1.0\,\mathrm{M}$')}
removed, kept = [], []
def update(match):
    block = match.group()
    name = re.search(r'\\includegraphics\[[^\]]*\]\{([^}]+)\}', block).group(1)
    metal, parameter = re.fullmatch(r'(Au|Ni|Cu|Pd|Pt|Ti)_(activity|NH3|Gly|CN)_three_levels_1x3_side_colors.pdf', name).groups()
    labels = re.findall(r'\\label\{([^}]+)\}', block)
    assert len(labels) == 1
    if parameter == 'CN':
        removed.extend(labels)
        return ''
    kept.extend(labels)
    description, bounds = specs[parameter]
    caption = f'{metal} Pourbaix diagrams with the {description} varied from {bounds}. Other concentrations were held at their reference values.'
    return '\n'.join([r'\begin{figure}[!htbp]', r'    \centering',
        rf'    \includegraphics[width=\textwidth]{{{name}}}', rf'    \caption{{{caption}}}',
        rf'    \label{{{labels[0]}}}', r'\end{figure}'])
after = intro + '\n\n' + re.sub(r'\\begin\{figure\}.*?\\end\{figure\}', update, body, flags=re.S)
assert len(removed) == 6 and len(kept) == 18
assert 'NO2' not in after and '_CN_' not in after and 'sweep' not in after
for path in SUBMISSION.glob('*.tex'):
    text = path.read_text()
    references = re.findall(r'\\(?:[cC]ref|ref|pageref|autoref)\{([^}]+)\}', text)
    assert not any(label in [part.strip() for part in refs.split(',')] for label in removed for refs in references), path
print(json.dumps({'before': before, 'after': after, 'source': str(source),
                  'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                  'removed_labels': removed, 'retained_labels': kept}))
