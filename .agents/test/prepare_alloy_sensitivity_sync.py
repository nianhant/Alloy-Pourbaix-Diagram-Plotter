"""Prepare the alloy sensitivity SI addition without modifying submission."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUBMISSION = Path('/global/homes/n/nianhant/data/stability_paper/updated_stability_manuscript/submission')
ART = ROOT / 'figures/alloy_pourbaix_diagrams/sensitivity_triptychs_with_cn_sweep'
source = SUBMISSION / 'SI.tex'
before = source.read_text()
anchor = '\\clearpage\n\\section{Cyclic Voltammetry}'
assert before.count(anchor) == 1
addition = (r'\clearpage' + '\n' + r'\section{Alloy Pourbaix Sensitivity Analysis}' + '\n'
            + r'\label{sec:si_alloy_pourbaix_sensitivity}' + '\n'
            + r'\input{alloy_pourbaix_sensitivity_SI}' + '\n\n')
assert 'alloy_pourbaix_sensitivity_SI' not in before
after = before.replace(anchor, addition + anchor)
intro = (r'The concentration sensitivity of equiatomic TiNi, TiCu, and PdAu was examined '
         r'at 298.15~K by varying one concentration at a time. The nominal aqueous metal-ion '
         r'concentration was varied from $10^{-6}$ to $10^{-2}\,\mathrm{M}$, the total ammonia '
         r'concentration from $0.002$ to $0.1\,\mathrm{M}$, and the total glycine concentration '
         r'from $0.01$ to $1.0\,\mathrm{M}$. Other concentrations were held at their reference '
         r'values: nominal aqueous metal-ion concentration $10^{-4}\,\mathrm{M}$, '
         r'$[\ce{NH3}]_{\mathrm{tot}}=0.02\,\mathrm{M}$, and '
         r'$[\mathrm{Gly}]_{\mathrm{tot}}=0.1\,\mathrm{M}$, without added cyanide. '
         r'The nominal metal-ion concentration corresponds to the activity input under '
         r'the ideal-solution convention. Each black outline marks pH 10.5--13.5 and '
         r'$-2.0$ to $2.3$~V vs RHE, converted to SHE and clipped to the plotted range.')
specs = [('activity', 'nominal aqueous metal-ion concentration', r'$10^{-6}$ to $10^{-2}\,\mathrm{M}$'),
         ('NH3', 'total ammonia concentration', r'$0.002$ to $0.1\,\mathrm{M}$'),
         ('Gly', 'total glycine concentration', r'$0.01$ to $1.0\,\mathrm{M}$')]
blocks, files = [], []
for alloy in ('TiNi', 'TiCu', 'PdAu'):
    for parameter, description, bounds in specs:
        stem = f'{alloy}_{parameter}_three_levels_1x3_side_colors'
        caption = f'{alloy} Pourbaix diagrams with the {description} varied from {bounds}. Other concentrations were held at their reference values.'
        blocks.append('\n'.join([r'\begin{figure}[!htbp]', r'    \centering',
            rf'    \includegraphics[width=\textwidth]{{{stem}.pdf}}',
            rf'    \caption{{{caption}}}', rf'    \label{{fig:{alloy}_{parameter}_sensitivity_SI}}', r'\end{figure}']))
        for suffix in ('.pdf', '.tiff'):
            path = ART / alloy / (stem + suffix)
            destination = SUBMISSION / path.name
            assert path.is_file() and not destination.exists(), destination
            files.append({'source': str(path), 'destination': str(destination),
                          'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
text = intro + '\n\n' + '\n\n'.join(blocks) + '\n'
assert '_CN_' not in text and 'NO2' not in text
destination = SUBMISSION / 'alloy_pourbaix_sensitivity_SI.tex'
assert not destination.exists()
print(json.dumps({'before': before, 'after': after, 'section': text,
                  'source': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                  'section_destination': str(destination), 'files': files}))
