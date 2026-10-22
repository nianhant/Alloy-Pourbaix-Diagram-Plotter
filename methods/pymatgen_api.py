import json
import os
import sys
import types
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter,IonEntry
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from mp_api.client import MPRester
import re


def install_pymatgen_entries_compat_alias():
    """Support MP data serialized with the older pymatgen.core.entries path."""
    try:
        import pymatgen.entries.computed_entries as computed_entries
    except ModuleNotFoundError:
        return

    entries_module = types.ModuleType("pymatgen.core.entries")
    entries_module.computed_entries = computed_entries
    sys.modules.setdefault("pymatgen.core.entries", entries_module)
    sys.modules.setdefault("pymatgen.core.entries.computed_entries", computed_entries)

    try:
        import pymatgen.entries.compatibility as compatibility
    except ModuleNotFoundError:
        return

    sys.modules.setdefault("pymatgen.analysis.compatibility", compatibility)

    try:
        from emmet.core.vasp.calculation import PotcarSpec
    except ModuleNotFoundError:
        return

    if not hasattr(PotcarSpec, "get"):
        PotcarSpec.get = lambda self, key, default=None: self.model_dump().get(key, default)


install_pymatgen_entries_compat_alias()


def get_mpr():
    mpr_key = os.getenv("MP_API_KEY") or os.getenv("MAPI_KEY")
    if not mpr_key:
        raise ValueError("Set the MP_API_KEY environment variable before fetching Materials Project data.")
    return MPRester(mpr_key, use_document_model=False)


def save_as_json(data, data_dir, filename):
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    with open(f'{data_dir}/{filename}.json', 'w') as json_file:
        json.dump(data, json_file, indent=4)

        
def get_alloy_ion_formation_energy(data_dir, metal_1, metal_2):
    ion_formation_eng = {}
    
    mpr = get_mpr()
    pourbaix_entries = mpr.get_pourbaix_entries([metal_1, metal_2])
    for entry in pourbaix_entries:
        formula = entry.name
        if entry.phase_type == 'Ion':
            trimmed_formula = formula#.split('(aq)')[0]
            trimmed_formula = re.sub(r'\[([+-])(\d+)\]', r'[\2\1]', trimmed_formula)
            if trimmed_formula not in ion_formation_eng:
                ion_formation_eng[trimmed_formula] = entry.uncorrected_energy
            elif ion_formation_eng[trimmed_formula] > entry.uncorrected_energy:
                ion_formation_eng[trimmed_formula] = entry.uncorrected_energy
    if 'FeOH[2+]' in ion_formation_eng:
        ion_formation_eng.pop('FeOH[2+]') 
    save_as_json(ion_formation_eng, data_dir, f'{metal_1}_{metal_2}_ion_formation_energy')
    
    return ion_formation_eng


def get_alloy_solid_formation_energy(data_dir, metal_1, metal_2, dimensionally_stable = False):
    # returns per metal formation energy, not total
    solid_formation_eng = {}    
    mpr = get_mpr()
    pourbaix_entries = mpr.get_pourbaix_entries([metal_1, metal_2])
    solid_entries = [entry for entry in pourbaix_entries if entry.phase_type == "Solid"]
    
    entries_HO = [ComputedEntry("H", 0), ComputedEntry("O", 2.46)]
    solid_pd = PhaseDiagram(solid_entries + entries_HO)
    if dimensionally_stable:
        solid_entries = list(set(solid_pd.stable_entries) - set(entries_HO))
    else:
        solid_entries = list(set(solid_pd.entries) - set(entries_HO))
    for entry in solid_entries:
        if entry.phase_type == 'Solid':
            formula = entry.composition.formula
            trimmed_formula = formula.replace(" ", "")
            trimmed_formula = re.sub(r'([A-Za-z])1(?=[A-Za-z]|$)', r'\1',trimmed_formula)
            if trimmed_formula not in solid_formation_eng:
                solid_formation_eng[trimmed_formula] = entry.uncorrected_energy
            elif solid_formation_eng[trimmed_formula] > entry.uncorrected_energy: 
                solid_formation_eng[trimmed_formula] = entry.uncorrected_energy
    save_as_json(solid_formation_eng, data_dir, f'{metal_1}_{metal_2}_solid_formation_energy')
    return solid_formation_eng
