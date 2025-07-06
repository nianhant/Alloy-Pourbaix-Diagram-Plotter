import numpy as np

class GridMaker:
    def __init__(self, pH_range, V_range, ligand_concentration, grid_size):
        self.pH_range = pH_range
        self.V_range = V_range
        self.ligand_concentration = ligand_concentration
        self.grid_size = grid_size

        self.pH_values = np.linspace(pH_range[0], pH_range[1], grid_size)
        self.V_values = np.linspace(V_range[0], V_range[1], grid_size)
        self.pH_grid, self.V_grid = np.meshgrid(self.pH_values, self.V_values)

        self.ligand_grid_dict = {
            'NH3': self.generate_ligand_grid('NH3', pKa=9.25),
            'Gly': self.generate_ligand_grid('Gly', pKa1=2.35, pKa2=9.78),
            'CN': self.generate_ligand_grid('CN', pKa=9.2)
        }

    def generate_ligand_grid(self, ligand, pKa=None, pKa1=None, pKa2=None):
        """Generates ligand activity grids based on given pKa values."""
        if ligand not in self.ligand_concentration:
            raise ValueError(f"Ligand {ligand} is not defined in concentration dictionary.")

        ligand_tot = self.ligand_concentration[ligand]
        if pKa:  
            return -np.log10(ligand_tot) + np.log10(1 + 10 ** (pKa - self.pH_grid))
        elif pKa1 and pKa2:  
            return (-np.log10(ligand_tot) + (pKa1 - self.pH_grid) +
                    np.log10(1 + 10 ** (self.pH_grid - pKa1) + 1 / (10 ** (self.pH_grid - pKa2))))
        else:
            raise ValueError("Either pKa or (pKa1, pKa2) must be provided for ligand grid generation.")
