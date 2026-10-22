from pathlib import Path
import xml.etree.ElementTree as ET

path = Path('/global/homes/n/nianhant/data/stability_paper/single-metal-pourbaix-plotter/examples/figures/figures/six_panel_nh3_gly_0p1_no_cn/Au_Ni_Cu_Pd_Pt_Ti_activity_1e-4_NH3_0p02_Gly_0p1_2x3_overwide_side_colors_all_axes.svg')
root = ET.parse(path).getroot()
ns = {'s': 'http://www.w3.org/2000/svg'}
for group in root.findall('.//s:g', ns):
    if not group.get('id', '').startswith('axes_'):
        continue
    if int(group.get('id').split('_')[1]) % 2 == 0:
        print(group.get('id'))
        for child in group:
            label = ''.join(child.itertext()).strip()
            uses = child.findall('.//s:use', ns)
            if child.get('id', '').startswith('line2d_'):
                print('MARKER', [u.get('style') for u in uses])
            if child.get('id', '').startswith('text_'):
                print('LABEL', label)
