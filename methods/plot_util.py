import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def generate_legends(ax, all_species_tuples_global, species_colors, plot_accessories, pH_exp_range, V_exp_range):
    solid_handles = [Line2D([], [], linestyle="None", label="Solid Species")]
    aqueous_handles = [Line2D([], [], linestyle="None", label="Aqueous Species")]

    for species in all_species_tuples_global:
        label = plot_accessories.format_species_combo(species)
        handle = Line2D([0], [0], marker='o', color='w', markerfacecolor=species_colors[species], markersize=8, label=label)
        (solid_handles if '(aq)' not in label else aqueous_handles).append(handle)

    max_len = max(len(solid_handles), len(aqueous_handles))
    solid_handles.extend([Line2D([], [], linestyle="None", label="")] * (max_len - len(solid_handles)))
    aqueous_handles.extend([Line2D([], [], linestyle="None", label="")] * (max_len - len(aqueous_handles)))

    accessory_handles = [
            Line2D([0], [0], color="g", linestyle="-", linewidth=1, label=f'Exp condition\nV vs RHE={V_exp_range[0]}-{V_exp_range[1]}\npH={pH_exp_range[0]}-{pH_exp_range[1]}'),
            Line2D([0], [0], color="b", linestyle="--", linewidth=1, label=r'H$_2$O Reduction'),
            Line2D([0], [0], color="r", linestyle="--", linewidth=1, label=r'H$_2$O Oxidation')
        ]
    accessory_legend = ax.legend(handles=accessory_handles, loc="upper right", frameon=True)

    combined_legend = ax.legend(
        solid_handles + aqueous_handles,
        [h.get_label() for h in solid_handles + aqueous_handles],
        loc="center left", frameon=True, bbox_to_anchor=(1.01, 0.5), ncol=2,
        columnspacing=1.0, handletextpad=1.0
    )
    ax.add_artist(accessory_legend)
