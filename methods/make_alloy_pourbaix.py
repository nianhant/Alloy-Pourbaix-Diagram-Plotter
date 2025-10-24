
from matplotlib import pyplot as plt
import itertools
from pymatgen.core import Composition

from .data_loader import SpeciesDataLoader
from .thermodynamics import PourbaixData, PourbaixAnalyzer
from .grid import GridMaker
from .visualization import GridVisualizer, PlotAccessories
from .utils import parse_composition, format_comp_dict
from .plot_util import generate_legends
from .set_publication_style import set_publication_style


set_publication_style()
for metal_1, metal_2 in list(itertools.combinations(metal_list,2)):
    print(metal_1, metal_2)

metal_1, metal_2 = 'Ni','Cu'
mu_ligand = {'NH3': -0.276037, 'Gly': -3.263014109, 'CN': 1.786800089}
T_list = 298.15
activity_list = 1e-4
ligand_concentration = {'NH3':0., 'NO2':0, 'Gly': 0.1, 'CN':0}

species_data = SpeciesDataLoader(metal_1, metal_2, mu_ligand, T, data_dir = data_dir)
prod_comp_dict = {
    format_comp_dict(parse_composition(alloy)): Composition({
        metal: count / sum(parse_composition(alloy).values())
        for metal, count in parse_composition(alloy).items()
    }) 
    for alloy in species_data.solid_eng if len(parse_composition(alloy)) >= 2 and 
    all(count / sum(parse_composition(alloy).values()) == 0.5 for count in parse_composition(alloy).values())
}
print(prod_comp_dict)