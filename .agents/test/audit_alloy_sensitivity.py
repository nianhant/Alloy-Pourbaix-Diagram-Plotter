"""Audit saved sensitivity results and make small visual review sheets."""
import csv
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "research/2026-10-02_alloy_concentration_sensitivity"
FIGURES = ROOT / "figures/alloy_pourbaix_diagrams/sensitivity_triptychs_with_cn_sweep"

manifest = json.loads((RESULTS / "manifest.json").read_text())
rows = list(csv.DictReader((RESULTS / "sensitivity_summary.csv").open()))
assert len(rows) == 36 and len(manifest["figures"]) == 12
assert manifest["grid_size"] == 1600
assert len(list(RESULTS.glob("*/*.npz"))) == 27
for path in RESULTS.glob("*/*.npz"):
    with np.load(path) as data:
        assert data["stable_indices"].shape == (1600, 1600)
        assert data["pH"].shape == data["E_SHE"].shape == (1600,)
        labels = json.loads(path.with_suffix(".json").read_text())["species_list"]
        assert data["stable_indices"].min() >= 0
        assert data["stable_indices"].max() < len(labels)
for row in rows:
    condition = {key: float(row[key]) for key in manifest["baseline"]}
    varied = row["parameter"]
    assert condition[varied] in manifest["levels"][varied]
    for key, value in condition.items():
        assert key == varied or value == manifest["baseline"][key]
    for scope in ("full", "experimental"):
        assert 0 <= float(row[f"{scope}_changed_fraction"]) <= 1
        if condition == manifest["baseline"]:
            assert float(row[f"{scope}_changed_fraction"]) == 0
for alloy in ("TiNi", "TiCu", "PdAu"):
    previews = []
    for parameter in manifest["levels"]:
        stem = FIGURES / alloy / f"{alloy}_{parameter}_three_levels_1x3_side_colors"
        for suffix in (".png", ".pdf", ".svg", ".tiff"):
            assert stem.with_suffix(suffix).stat().st_size > 0
        svg = stem.with_suffix(".svg").read_text()
        for panel in "abc":
            assert f"({panel}) {alloy}:" not in svg
        with Image.open(stem.with_suffix(".tiff")) as source:
            assert source.mode == "RGB"
            assert all(abs(dpi - 600) < 1 for dpi in source.info["dpi"])
        with Image.open(stem.with_suffix(".png")) as source:
            preview = source.convert("RGB")
        preview.thumbnail((2200, 1000))
        previews.append((parameter, preview))
    sheet = Image.new("RGB", (2220, sum(p.height + 30 for _, p in previews)), "white")
    draw = ImageDraw.Draw(sheet)
    y = 0
    for parameter, preview in previews:
        draw.text((10, y + 3), f"{alloy} / {parameter}", fill="black")
        sheet.paste(preview, (10, y + 25))
        y += preview.height + 30
    sheet.save(ROOT / ".agents/test" / f"{alloy}_sensitivity_review.png")
    for parameter in manifest["levels"]:
        selected = [r for r in rows if r["alloy"] == alloy and r["parameter"] == parameter]
        print(alloy, parameter, "experimental changed %:",
              [round(100 * float(r["experimental_changed_fraction"]), 3) for r in selected],
              "aqueous-containing %:",
              [round(100 * float(r["experimental_aqueous_containing"]), 3) for r in selected])
print("PASS: 27 grids, 36 comparison records, 12 triptychs / 48 exports")
