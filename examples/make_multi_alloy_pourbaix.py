from pathlib import Path

import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from methods.runner import plot_pourbaix


# METALS = ("Ni", "Ti")
# METALS = ("Cu", "Ti")
# METALS = ("Pd", "Au")
Metal_list = [ ("Ni", "Ti"),
                ("Cu", "Ti"),
                ("Au", "Pd"),
                ("Ni", "Co"),
                ("Ni", "Cu"),
                ("Pd", "Ag"),
                ("Pd", "Cu"),
                ("Pd", "Ni"),
                ("Pd", "Pt"),
                ("Ti", "Ag"),
                ("Ti", "Au"),
                ("Ti", "Mg"),
                ("Ti", "Mn"),
                ("Ti", "Pd"),
                ("Ti", "Pt"),
                ("Ti", "Zn"),
                # ("Ti", "Ag"),
              ]
Metal_list = [ ("Ti", "Fe"),]
ligand_list = [ {"NH3": 0.02, "NO2": 0, "Gly": 0.1, "CN": 0},
                {"NH3": 0.02, "NO2": 0, "Gly": 0.1, "CN": 1e-4},
                {"NH3": 0, "NO2": 0, "Gly": 0, "CN": 0}
              ]
MU_LIGAND = {"NH3": -0.276037, "Gly": -3.263014109, "CN": 1.786800089}
TEMPERATURE_K = 298.15
ACTIVITY = 1e-4
LIGAND_CONCENTRATION = {"NH3": 0.02, "NO2": 0, "Gly": 0.1, "CN": 0}

DATA_DIR = REPO_ROOT / "data"
OUTPUT_DIR = REPO_ROOT / "figures" / "alloy_pourbaix_diagrams"
OUTPUT_DIR = f'/global/homes/n/nianhant/data/stability_paper/updated_stability_manuscript/Figures/pourbaix_diagrams'


PH_EXP_RANGE = (11.5, 13.5)
V_EXP_RANGE = (-2, 2.0)

for ligand_concentration in ligand_list:
    for metals in Metal_list:
        metal_1, metal_2 = metals
        plot_pourbaix(
            metal_1,
            metal_2,
            MU_LIGAND,
            TEMPERATURE_K,
            ACTIVITY,
            ligand_concentration,
            DATA_DIR,
            PH_EXP_RANGE,
            V_EXP_RANGE,
            save_fig=True,
            outdir=OUTPUT_DIR,
            save_pdf=False,
        )
# def main():
#     metal_1, metal_2 = METALS
#     plot_pourbaix(
#         metal_1,
#         metal_2,
#         MU_LIGAND,
#         TEMPERATURE_K,
#         ACTIVITY,
#         LIGAND_CONCENTRATION,
#         DATA_DIR,
#         PH_EXP_RANGE,
#         V_EXP_RANGE,
#         save_fig=True,
#         outdir=OUTPUT_DIR,
#         save_pdf=True,
#     )


# if __name__ == "__main__":
#     main()
