"""
Plot showing sinusoidal field distribution across bands.
Standalone Python 3 script (not Abaqus).
"""

import numpy as np
import matplotlib.pyplot as plt
import math


def set_publication_style():
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman'],
        'font.size': 11,
        'axes.labelsize': 13,
        'axes.titlesize': 14,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'figure.figsize': (8, 5),
    })


def calculate_field_values(num_bands=5, peak_value=0.15, min_value=0.0):
    positions = []
    values = []
    num_points = num_bands + 1  # cosine reaches 0 at the virtual (num_bands+1)th point
    for i in range(num_bands):
        pos = i / (num_points - 1.0)  # 0, 0.2, 0.4, 0.6, 0.8 for 5 bands
        # Raised cosine: gradient = 0 at both center and edge (smooth step)
        field_val = (peak_value - min_value) * (1 + math.cos(math.pi * pos)) / 2.0 + min_value
        positions.append(pos)
        values.append(round(field_val, 3))
    return positions, values


def main():
    set_publication_style()

    num_bands = 5
    peak = 0.15
    min_val = 0.0

    band_positions, band_values = calculate_field_values(num_bands, peak, min_val)

    # Band edges (each band spans 20% of the region)
    band_width = 1.0 / num_bands

    # Smooth cosine curve (full symmetric profile, cosine reaches 0 at 120%)
    # cos(pi * |x| / 1.2) maps physical distance to raised cosine
    extent = 1.2  # cosine zero-crossing at 120% of physical region
    # Physical region: -100% to +100% (solid line)
    x_solid = np.linspace(-1, 1, 400)
    y_solid = (peak - min_val) * (1 + np.cos(np.pi * np.abs(x_solid) / extent)) / 2.0 + min_val
    # Extension to 120% where cosine reaches 0 (dashed line)
    x_ext_pos = np.linspace(1, extent, 50)
    x_ext_neg = np.linspace(-extent, -1, 50)
    y_ext_pos = (peak - min_val) * (1 + np.cos(np.pi * x_ext_pos / extent)) / 2.0 + min_val
    y_ext_neg = (peak - min_val) * (1 + np.cos(np.pi * np.abs(x_ext_neg) / extent)) / 2.0 + min_val

    # Colours for each band (edge to center to edge)
    colours = ['#2166ac', '#67a9cf', '#d1e5f0', '#fddbc7', '#ef8a62']

    fig, ax = plt.subplots()

    # Bar chart: mirrored bands (negative side then positive side)
    for i in range(num_bands):
        # Positive side: band edge starts at i * band_width
        pos_edge = i * band_width
        # Negative side: mirror
        neg_edge = -(i + 1) * band_width

        label = 'Band {} ({:.3f})'.format(i + 1, band_values[i])

        # Positive side bar
        ax.bar(
            pos_edge * 100,
            band_values[i],
            width=band_width * 100,
            align='edge',
            color=colours[i],
            edgecolor='grey',
            linewidth=0.8,
            alpha=0.55,
            label=label,
            zorder=2
        )
        # Negative side bar (no label to avoid duplicate legend entries)
        ax.bar(
            neg_edge * 100,
            band_values[i],
            width=band_width * 100,
            align='edge',
            color=colours[i],
            edgecolor='grey',
            linewidth=0.8,
            alpha=0.55,
            zorder=2
        )

    # Sinusoid curve (solid within physical region)
    ax.plot(
        x_solid * 100, y_solid,
        color='black', linewidth=1.8, linestyle='-',
        label='Raised cosine profile',
        zorder=3
    )
    # Dashed extension to 120% where cosine reaches 0
    ax.plot(x_ext_pos * 100, y_ext_pos, color='black', linewidth=1.2, linestyle='--', zorder=3)
    ax.plot(x_ext_neg * 100, y_ext_neg, color='black', linewidth=1.2, linestyle='--', zorder=3)

    # Sample points on the curve (both sides)
    # band_positions are in cosine space (0 to 1 = 0% to 120%); scale to physical %
    for i in range(num_bands):
        pos_x = band_positions[i] * 120
        neg_x = -band_positions[i] * 120
        # Positive side
        ax.plot(
            pos_x, band_values[i],
            marker='o', markersize=9,
            markerfacecolor=colours[i], markeredgecolor='black',
            markeredgewidth=1.2,
            zorder=4
        )
        # Negative side (skip center to avoid double-plotting)
        if i > 0:
            ax.plot(
                neg_x, band_values[i],
                marker='o', markersize=9,
                markerfacecolor=colours[i], markeredgecolor='black',
                markeredgewidth=1.2,
                zorder=4
            )
        # Annotate positive side only
        ax.annotate(
            '{:.3f}'.format(band_values[i]),
            xy=(pos_x, band_values[i]),
            xytext=(0, 12),
            textcoords='offset points',
            ha='center', fontsize=9,
            fontweight='bold'
        )

    ax.set_xlabel('Distance from Centre (%)')
    ax.set_ylabel('Field Amplitude')
    ax.set_title('Sinusoidal Predefined Field Distribution')
    ax.set_xlim(-125, 125)
    ax.set_ylim(-0.01, peak * 1.25)
    ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax.set_xticks([-120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120])

    plt.tight_layout()
    plt.savefig('field_band_distribution.pdf', dpi=300)                   
    print("Saved to field_band_distribution.pdf") 
    plt.savefig(r'C:\Users\cmb247\repos\Paper-DCM\figures\field_band_distribution.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
