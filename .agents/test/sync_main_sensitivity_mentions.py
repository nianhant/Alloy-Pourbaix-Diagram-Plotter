"""Apply two reviewed sensitivity paragraphs with direct SI figure references."""
from pathlib import Path
import re
import shutil

stage = Path(__file__).resolve().parent / 'main_sensitivity_mentions'
submission = Path('/global/homes/n/nianhant/data/stability_paper/updated_stability_manuscript/submission')
target = submission / 'main.tex'
before = (stage / 'main.before.tex').read_text()
after = (stage / 'main.after.tex').read_text()
assert target.read_text().rstrip() == before.rstrip(), 'Concurrent edits; abort'
single = (submission / 'pourbaix_sensitivity_SI.tex').read_text()
binary = (submission / 'alloy_pourbaix_sensitivity_SI.tex').read_text()
labels = re.findall(r'\\label\{([^}]+)\}', single + binary)
needed = ['fig:Au_activity_sensitivity_SI', 'fig:Pd_NH3_sensitivity_SI',
          'fig:TiNi_activity_sensitivity_SI', 'fig:PdAu_NH3_sensitivity_SI',
          'fig:TiNi_Gly_sensitivity_SI', 'fig:TiCu_Gly_sensitivity_SI']
assert all(labels.count(label) == 1 for label in needed)
assert re.findall(r'\\label\{([^}]+)\}', single)[-1] == 'fig:Pd_NH3_sensitivity_SI'
assert re.findall(r'\\label\{([^}]+)\}', binary)[-1] == 'fig:PdAu_NH3_sensitivity_SI'
assert after.count('{') - after.count('}') == before.count('{') - before.count('}')
added = [line for line in after.splitlines() if line not in before.splitlines()]
assert len(added) == 2 and all(line.startswith(r'\textcolor{red}{') for line in added)
assert not any('sec:si_' in line for line in added)
backup = stage / 'main.original.tex'
assert not backup.exists()
shutil.copy2(target, backup)
assert target.read_text().rstrip() == before.rstrip()
shutil.copy2(stage / 'main.after.tex', target)
assert target.read_bytes() == (stage / 'main.after.tex').read_bytes()
print('PASS: two red paragraphs added; direct SI figure labels and range endpoints verified. Exact original backed up locally.')
