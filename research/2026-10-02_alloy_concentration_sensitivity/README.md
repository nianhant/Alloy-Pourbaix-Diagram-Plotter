# Alloy concentration sensitivity results

Completed for equiatomic TiNi, TiCu, and PdAu using the existing alloy model
at 298.15 K. The AtomisticSkills/Pourbaix guidance informed the input provenance,
reproducibility records, interpretation limits, and visual checks; no additional
DFT calculations or new thermodynamic validation were performed.

## Conditions

The generator verifies the `BASELINE` and `LEVELS` constants against the
single-metal `examples/generate_sensitivity_triptychs.py` source before running.
The source path and SHA-256 are recorded in `manifest.json`.

| Varied input | Panel (a) | Panel (b) | Panel (c) |
| --- | ---: | ---: | ---: |
| Ion activity input | 1e-6 | 1e-4 | 1e-2 |
| NH3 total / M | 0.002 | 0.02 | 0.1 |
| Gly total / M | 0.01 | 0.1 | 1 |
| CN total / M | 0 | 1e-4 | 1e-3 |

Only the named input varies. All other inputs remain at activity 1e-4,
NH3 total 0.02 M, Gly total 0.1 M, CN total 0 M, and NO2 total 0 M.
Activity is dimensionless; figure labels follow the reference's ideal-solution
molarity convention. The CN sweep has its baseline in panel (a), whereas
the other sweeps have their baseline in panel (b).

## Quantitative comparison

The table gives the percentage of grid points whose stable species combination
differs from baseline within the highlighted window: pH 10.5–13.5 and
E_RHE -2 to 2.3 V, intersected with the displayed E_SHE -2 to 3 V domain.
Each cell lists panels (a), (b), and (c), in that order.

| Alloy | Activity sweep / % | NH3 sweep / % | Gly sweep / % | CN sweep / % |
| --- | --- | --- | --- | --- |
| TiNi | 3.88, 0, 13.36 | 0, 0, 0 | 21.16, 0, 5.55 | 0, 54.14, 60.74 |
| TiCu | 10.10, 0, 2.45 | 0, 0, 0.45 | 2.45, 0, 13.89 | 0, 24.17, 34.07 |
| PdAu | 26.75, 0, 6.68 | 3.72, 0, 3.41 | 0, 0, 0 | 0, 38.69, 45.28 |

- CN addition gives the largest stable-combination changes among these tested
  sweeps for each alloy. At 1e-3 M CN, the window's aqueous-containing area
  fractions are 60.74% (TiNi), 34.36% (TiCu), and 46.93% (PdAu), compared with
  baseline fractions of 45.53%, 2.45%, and 7.33%, respectively.
- Lowering Gly to 0.01 M reduces TiNi's aqueous-containing area from 45.53%
  to 24.37%. Raising Gly to 1 M increases TiCu's from 2.45% to 16.35%.
- NH3 changes are comparatively small within this window. TiNi shows no
  resolved change at these levels; PdAu shows no resolved Gly response.
  These statements concern this model, species inventory, window, and grid,
  not a general experimental absence of ligand effects.
- Changing ion activity also matters, particularly for PdAu: the 1e-6 input
  changes 26.75% of the window relative to the 1e-4 baseline.

`sensitivity_summary.csv` contains the unrounded values, full-domain comparisons,
and category area fractions. Categories use the existing plotting priority:
ligand-complex-containing combinations first, then aqueous-ion-containing
combinations, then solid categories. A combination containing both a solid and
an aqueous species counts as aqueous-containing. Changed identity can include
one complex replacing another without increasing the total aqueous-containing
area. These are geometric sensitivity metrics, not phase yields, dissolution
probabilities, rates, or an uncertainty analysis of formation energies.

## Files and verification

Figures are under
`../../figures/alloy_pourbaix_diagrams/sensitivity_triptychs_with_cn_sweep/`,
with separate TiNi, TiCu, and PdAu directories. Each contains four 1×3
triptychs in PDF, SVG, PNG, and RGB LZW-compressed TIFF (600 dpi).
These are overwide masters (approximately 26–29 cm wide), not figures already
reduced to a journal column width; check text size at the final manuscript scale.

The figures retain square 4.7 cm axes, common ranges and ticks, inward ticks,
normal-weight text, category/centroid-sorted side legends, and no connectors.
One palette per alloy is reused across all sweeps. Existing matched single-metal
complex colors are preserved; newly appearing complexes extend the RdPu palette.
Simple aqueous-ion-only region labels remain omitted, as in the parent collage.
Only the varied input is stated in each panel title.

Validation passed for 27 distinct 1600×1600 grids, 36 comparison records,
12 triptychs, and 48 exports. Baseline comparisons are exactly zero.
All generated triptychs were visually reviewed for clipping and crowding.
Rendering assertions verify square axes, common limits, inward ticks, and
legend text within the canvas. Export checks verify RGB TIFFs and 600 dpi.

Each alloy subdirectory here contains nine `.npz` grids and corresponding
`.json` conditions/species mappings. `manifest.json` records the source settings,
dataset hashes, palettes, and figure dimensions. Reproduce from the repository root:

```bash
/pscratch/sd/n/nianhant/conda/envs/base-agent/bin/python examples/generate_alloy_sensitivity_triptychs.py
/pscratch/sd/n/nianhant/conda/envs/base-agent/bin/python .agents/test/audit_alloy_sensitivity.py
```

## Suggested caption template

Concentration sensitivity of the equiatomic [ALLOY] Pourbaix diagram at 298.15 K.
Panels (a–c) vary [INPUT] over [VALUES]; all other inputs are held at ion
activity 1e-4, NH3 total 0.02 M, Gly total 0.1 M, and CN total 0 M.
Dashed blue/red lines indicate the water reduction/oxidation boundaries.
The black outline denotes pH 10.5–13.5 and E_RHE -2 to 2.3 V, clipped to the
displayed potential domain. Side legends identify selected stable combinations.
Calculations use the existing alloy thermodynamic/speciation model on a
1600×1600 pH–potential grid.
