# Alloy Pourbaix Diagram Plotter

Generate a binary alloy Pourbaix diagram with N-containing ligands.

For the paper-release example, the repository keeps one generation script with
one alloy system and one ligand condition:

```bash
python examples/make_alloy_pourbaix.py
```

The example currently generates the Ni-Ti diagram at 298.15 K with:

- metal activity: `1e-4`
- ligand concentrations: `NH3=0.02 M`, `Gly=0.005 M`, `CN=0 M`

Outputs are written under `figures/alloy_pourbaix_diagrams/`.

To generate the RSC double-column TiNi/TiCu/PdAu collage with solid-region
legends, 600 dpi TIFF output, and matching PDF/SVG files, run:

```bash
python examples/make_rsc_pourbaix_collage.py
```

The collage uses a metal-ion activity of `1e-4 M`, `NH3=0.02 M`,
`Gly=0.1 M`, and an experimental pH box spanning 10.5--13.5.
