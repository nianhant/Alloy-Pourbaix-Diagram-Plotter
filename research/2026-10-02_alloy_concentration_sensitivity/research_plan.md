# Alloy concentration sensitivity

Reproduce the explicitly requested single-metal concentration protocol for
equiatomic TiNi, TiCu, and PdAu using the current alloy thermodynamic model.
No new energy calculations or changes to model parameters are introduced.

## Reference protocol

Source: `/global/homes/n/nianhant/data/stability_paper/single-metal-pourbaix-plotter/examples/generate_sensitivity_triptychs.py`.
The requested output folder is spelled `sensitivity_triptychs_with_cn_sweep`.
Its baseline is ion activity 1e-4, NH3 total 0.02 M, Gly total 0.1 M, CN total
0 M, NO2 total 0 M, and 298.15 K. Each sweep changes one input only:

| Input | Three levels |
| --- | --- |
| Activity | 1e-6, 1e-4, 1e-2 |
| NH3 total (M) | 0.002, 0.02, 0.1 |
| Gly total (M) | 0.01, 0.1, 1 |
| CN total (M) | 0, 1e-4, 1e-3 |

## Method and outputs

1. Verify the source constants by parsing the reference script and record its hash.
2. Use the existing local species/energy datasets and alloy chemical-potential
   minimizer. Evaluate 1600 by 1600 uniform grid points over pH -2 to 16 and
   E_SHE -2 to 3 V. Cache the repeated baseline (9 unique cases per alloy).
3. Save the discrete stable-region grids and their species mapping for audit.
4. Reuse one color dictionary for each alloy across all concentration levels;
   retain the matched single-metal complex swatches, extending RdPu for newly
   stable complexes. Preserve category/centroid-sorted side legends.
5. Generate 12 figures, each with three square, inward-ticked panels, in TIFF,
   PDF, SVG and PNG. Concentration titles identify each condition.
6. Report stable-region changes and category area fractions both over the full
   grid and the visible experimental window (pH 10.5-13.5, E_RHE -2 to 2.3 V).
   This window is intersected with the displayed E_SHE range.
7. Verify input settings, shape, tick direction, grid mappings, fraction sums,
   baseline equality, and figure layout. Inspect all generated figure previews.

## Interpretation limits

Grid area fractions are geometrical sensitivity measures, not dissolution
probabilities, kinetics, or experimentally measured phase yields. Activity is
dimensionless with the same ideal-solution molarity convention as the source.
Ligand totals enter the repository's pH-dependent speciation approximation.
No claim of new thermodynamic validation or DFT verification is made.

## Skills and infrastructure

AtomisticSkills research/coding guidance and the Pourbaix skill inform provenance,
energy consistency and visual inspection. The existing alloy workflow supplies
the required calculations; no separate MLIP/DFT workflow is needed. The artifact
planning tools described by the skill are unavailable in this environment, so
this reproducible local plan records the user-authorized protocol.
