from pathlib import Path

import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from methods.runner import plot_pourbaix


METALS = ("Ni", "Ti")
# METALS = ("Cu", "Ti")
METALS = ("Pd", "Au")

MU_LIGAND = {"NH3": -0.276037, "Gly": -3.263014109, "CN": 1.786800089}
TEMPERATURE_K = 298.15
ACTIVITY = 1e-4
LIGAND_CONCENTRATION = {"NH3": 0.02, "NO2": 0, "Gly": 0.005, "CN": 0}

DATA_DIR = REPO_ROOT / "data"
OUTPUT_DIR = REPO_ROOT / "figures" / "alloy_pourbaix_diagrams"

PH_EXP_RANGE = (11.5, 13.5)
V_EXP_RANGE = (-2, 2.3)


def main():
    metal_1, metal_2 = METALS
    plot_pourbaix(
        metal_1,
        metal_2,
        MU_LIGAND,
        TEMPERATURE_K,
        ACTIVITY,
        LIGAND_CONCENTRATION,
        DATA_DIR,
        PH_EXP_RANGE,
        V_EXP_RANGE,
        save_fig=True,
        outdir=OUTPUT_DIR,
        save_pdf=True,
    )


if __name__ == "__main__":
    main()
