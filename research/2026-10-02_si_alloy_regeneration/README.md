# Local SI alloy Pourbaix regeneration

Status: complete. All 19 original active alloy diagrams were visually inspected.
Generated and visually reviewed 11 replacement figures (19 panels, 44 exports),
plus two explicitly separated caption-based alternatives (4 panels, 8 exports).
The shared renderer has optional SI title/output-name arguments; its default
sensitivity-figure behavior remains unchanged.

Scope: regenerate only the active binary-alloy Pourbaix artwork referenced by
`submission/SI.tex` through `pourbaix_figures_SI.tex`. The submission directory
is read-only for this task. No single-metal, experimental, or main-text figures
are replaced. Commented-out figure blocks are excluded.

## Regeneration protocol

Reuse the current alloy model, local energy datasets, 298.15 K, equiatomic
compositions, and ion activity 1e-4. Ligand-containing SI panels use NH3 total
0.02 M and Gly total 0.1 M, with CN total either zero or 1e-4 M, and NO2 zero.
The current 1600×1600 grid and panel styling are retained: pH -2 to 16,
E_SHE -2 to 3 V, square 4.7 cm axes, inward ticks, normal-weight text, compact
side legends sorted by category then centroid height, and no connectors.
Panel titles contain letters only. Conditions and alloy names are in captions.
The user explicitly requests the updated experimental pH interval 10.5–13.5;
the current potential window is -2 to 2.3 V vs RHE, converted to SHE and
clipped at the displayed limits. Local caption drafts use these values.

The AtomisticSkills and Pourbaix guidance informs input provenance and visual
inspection. This is an existing-model figure regeneration, not a new energy
calculation, literature study, or scientific model validation.

## Source inventory

| SI alloy | Active panels | Original modes |
| --- | ---: | --- |
| TiNi | 2 | NH3/Gly; NH3/Gly/CN |
| TiCu | 2 | H2O artwork (caption conflicts); NH3/Gly/CN |
| TiFe | 2 | NH3/Gly; NH3/Gly/CN |
| TiMg | 1 | NH3/Gly |
| TiMn | 1 | NH3/Gly |
| TiZn | 1 | NH3/Gly/CN |
| PdAg | 2 | NH3/Gly; NH3/Gly/CN |
| PdCu | 2 | NH3/Gly; NH3/Gly/CN |
| PdNi | 2 | NH3/Gly; NH3/Gly/CN |
| PdPt | 2 | NH3/Gly; NH3/Gly/CN |
| AuPd | 2 | H2O artwork (caption conflicts); NH3/Gly/CN |

TiCu(a) and AuPd(a) use filenames ending in H2O and visibly contain no
ligand-complex regions, but their original captions state NH3/Gly. This
discrepancy must be resolved explicitly with `--ambiguous-baseline artwork`
or `--ambiguous-baseline captions`. The source inventory preserves both
filenames and original LaTeX. The final manifest records the resolution.

For review, the primary output set uses `artwork`: TiCu(a) and AuPd(a) remain
ligand-free, and their local captions are corrected accordingly. Separate
`caption_alternatives/` folders contain NH3/Gly versions and matching caption
drafts. This is not a decision about which version to submit; that choice and
any sync are deferred to the user.

## Outputs and checking

Artwork: `../../figures/alloy_pourbaix_diagrams/si_regenerated/`.
Each SI figure group is assembled directly in Matplotlib, not by LaTeX
subfigures, and exported as PDF, SVG, PNG, and RGB 600-dpi TIFF.
`alloy_figures_SI_replacement.tex` is an alloy-only local draft preserving
existing figure/panel labels and updating the experimental-window captions.
It is not automatically synced. The overwide artwork follows the approved
panel style; inspect font size if reducing it to journal text width.

The generator records source and dataset hashes, verifies its source inputs
remain unchanged, and saves the computed grids and species mappings here.
Runtime rendering checks verify square axes, common limits, inward ticks,
and legend bounds. Generated figures must additionally be visually inspected.

All export/grid/caption checks passed (see `validation.json` in each result
folder). Existing alloy figure and subfigure labels are preserved exactly;
the local replacement uses one includegraphics per Matplotlib-assembled figure.
The alloy include file and all 19 source images match their recorded hashes.
`submission/SI.tex` changed concurrently after generation, while retaining the
same Pourbaix include. Its original/current hashes are recorded in the audit;
this task did not write to that file or to any part of the submission.
No LaTeX compilation was performed: `pdflatex` is unavailable in this session.

```bash
/pscratch/sd/n/nianhant/conda/envs/base-agent/bin/python examples/generate_si_alloy_panels.py --inspect-only
# For a complete run, choose the documented ambiguity resolution explicitly:
/pscratch/sd/n/nianhant/conda/envs/base-agent/bin/python examples/generate_si_alloy_panels.py --ambiguous-baseline artwork
/pscratch/sd/n/nianhant/conda/envs/base-agent/bin/python examples/generate_si_alloy_panels.py --alloys TiCu AuPd --ambiguous-baseline captions --output-dir figures/alloy_pourbaix_diagrams/si_regenerated/caption_alternatives --results-dir research/2026-10-02_si_alloy_regeneration/caption_alternatives
/pscratch/sd/n/nianhant/conda/envs/base-agent/bin/python .agents/test/audit_si_alloy_panels.py
```
