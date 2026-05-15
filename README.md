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
