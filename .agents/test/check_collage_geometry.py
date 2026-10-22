"""Render and measure the collage axes at export time."""
import runpy
from pathlib import Path
from unittest.mock import patch

from matplotlib.figure import Figure
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
savefig = Figure.savefig


def measured_savefig(fig, filename, *args, **kwargs):
    if Path(filename).suffix == ".tiff":
        fig.canvas.draw()
        for ax in fig.axes:
            if not ax.axison:
                continue
            bounds = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            width, height = bounds.width * 2.54, bounds.height * 2.54
            assert abs(width - 4.7) < 1e-8 and abs(height - 4.7) < 1e-8
            assert ax.yaxis.get_tick_params()["direction"] == "in"
            assert ax.xaxis.get_tick_params()["direction"] == "in"
            assert tuple(ax.get_xlim()) == (-2, 16)
            assert tuple(ax.get_ylim()) == (-2, 3)
            assert [x for x in ax.get_xticks() if -2 <= x <= 16] == [0, 5, 10, 15]
            assert [y for y in ax.get_yticks() if -2 <= y <= 3] == [-2, -1, 0, 1, 2, 3]
            assert ax.images[0].get_array().shape[:2] == (1600, 1600)
            print(f"{ax.get_title(loc='left')}: {width:.6f} x {height:.6f} cm; all ticks inward; 1600 x 1600 grid", flush=True)
    return savefig(fig, filename, *args, **kwargs)


with patch.object(Figure, "savefig", measured_savefig):
    runpy.run_path(str(ROOT / "examples/make_rsc_pourbaix_collage.py"), run_name="__main__")

output = ROOT / "figures/alloy_pourbaix_diagrams/rsc_collage"
with Image.open(output / "TiNi_TiCu_PdAu_NH3=0.02M_Gly=0.1M_activity=1e-4M.tiff") as source:
    print(source.size, source.mode, source.info["dpi"])
    source.thumbnail((3000, 3000))
    source.save(Path(__file__).with_name("collage_preview.png"))
