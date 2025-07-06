import math
from ase.symbols import string2symbols

def parse_composition(formula):        
    element_only_form = formula.split('[')[0]
    symbol_list = string2symbols(element_only_form)
    return {symbol: symbol_list.count(symbol) for symbol in symbol_list if symbol not in ('O', 'H')}

def format_comp_dict(comp_dict):
    values = list(comp_dict.values())  
    gcd = math.gcd(*[int(round(v)) for v in values])  
    return ''.join(f"{metal}{int(comp_dict[metal]/gcd)}" for metal in comp_dict)
