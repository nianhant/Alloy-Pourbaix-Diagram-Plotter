import matplotlib.pyplot as plt

AREA_LABEL_FONT_SIZE = 20
AXIS_LABEL_FONT_SIZE = 26
TICK_LABEL_FONT_SIZE = AXIS_LABEL_FONT_SIZE

def set_publication_style():
    plt.rcParams.update({
        'font.size': AREA_LABEL_FONT_SIZE,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'axes.labelsize': AXIS_LABEL_FONT_SIZE,
        'axes.titlesize': AXIS_LABEL_FONT_SIZE,
        'xtick.labelsize': TICK_LABEL_FONT_SIZE,
        'ytick.labelsize': TICK_LABEL_FONT_SIZE,
        'legend.fontsize': AREA_LABEL_FONT_SIZE,
        'figure.titlesize': AXIS_LABEL_FONT_SIZE,
        'savefig.dpi': 300,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })
