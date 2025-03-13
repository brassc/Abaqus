#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Extract and plot a z-slice of logarithmic strain (LE) from an Abaqus .odb file
with a publication-quality visualization.
Compatible with Python 2.7.

Usage:
    abaqus python extract_strain.py <odb_file> [<step_name>] [<frame_number>] [<z_coord>] [<strain_component>] [<colormap>]

Arguments:
    odb_file          - Path to the Abaqus .odb file
    step_name         - Name of the step in the analysis (default: Step-1)
    frame_number      - Frame number to extract data from (default: -1, last frame)
    z_coord           - Z-coordinate value to extract the slice from (default: 0.0)
    strain_component  - Strain component to extract (default: LEMAG)
                        Options: LE11, LE22, LE33, LE12, LE13, LE23, LEMAG
    colormap          - Matplotlib colormap (default: viridis)
                        Options: viridis, plasma, inferno, magma, rainbow
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from abaqusConstants import *
from odbAccess import *

def set_publication_style():
    """Set matplotlib parameters for publication-quality figures."""
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman'],
        'mathtext.fontset': 'stix',
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 16,
        'axes.titleweight': 'bold',  # This makes titles bold
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 10,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.format': 'png',
        'savefig.bbox': 'tight',
        'axes.grid': True,
        'grid.alpha': 0.3,
        'grid.linestyle': '-',
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

def extract_le_strain_z_slice(odb_file, step_name, frame_number, z_coord, strain_component='LE11'):
    """
    Extract logarithmic strain data from a specific z-coordinate slice.
    
    Parameters:
    -----------
    odb_file : str
        Path to the Abaqus .odb file
    step_name : str
        Name of the step in the analysis
    frame_number : int
        Frame number to extract data from
    z_coord : float
        Z-coordinate value to extract the slice from
    strain_component : str
        Strain component to extract (default: LE11)
        Options: LE11, LE22, LE33, LE12, LE13, LE23, LEMAG (magnitude)
    
    Returns:
    --------
    x_coords : numpy.ndarray
        X-coordinates of nodes
    y_coords : numpy.ndarray
        Y-coordinates of nodes
    strain_values : numpy.ndarray
        Strain values at the nodes
    """
    print("Opening ODB file: %s" % odb_file)
    odb = openOdb(path=odb_file, readOnly=True)
    
    # Get the requested step and frame
    step = odb.steps[step_name]
    frame = step.frames[frame_number]
    
    # Get the field output for the logarithmic strain
    le_field = frame.fieldOutputs['LE']
    
    # Get the instance (assuming the first one if not specified)
    instance_name = odb.rootAssembly.instances.keys()[0]
    instance = odb.rootAssembly.instances[instance_name]
    
    # Collect nodes and their coordinates
    x_coords = []
    y_coords = []
    strain_values = []
    
    # Tolerance for finding nodes at the specified z-coordinate
    tol = 1e-6
    
    # Extract node positions from the instance
    nodes = instance.nodes
    
    # Create a mapping of node labels to their indices
    node_map = {}
    for i, node in enumerate(nodes):
        node_map[node.label] = i
    
    # Extract field values at nodes
    for value in le_field.values:
        # Get the node corresponding to this value
        node_label = value.nodeLabel
        
        # Skip if we can't find the node
        if node_label not in node_map:
            continue
        
        node_idx = node_map[node_label]
        node = nodes[node_idx]
        
        # Check if the node is at the specified z-coordinate
        if abs(node.coordinates[2] - z_coord) <= tol:
            x_coords.append(node.coordinates[0])
            y_coords.append(node.coordinates[1])
            
            # Extract the requested strain component
            if strain_component == 'LEMAG':
                # Calculate the magnitude of the strain tensor
                strain_tensor = [value.data[i] for i in range(len(value.data))]
                strain_mag = np.sqrt(sum(s*s for s in strain_tensor))
                strain_values.append(strain_mag)
            else:
                # Get the index of the requested component
                comp_map = {'LE11': 0, 'LE22': 1, 'LE33': 2, 'LE12': 3, 'LE13': 4, 'LE23': 5}
                comp_idx = comp_map.get(strain_component, 0)
                
                # Extract the component if available
                if comp_idx < len(value.data):
                    strain_values.append(value.data[comp_idx])
                else:
                    strain_values.append(0.0)
    
    # Close the ODB file
    odb.close()
    
    return np.array(x_coords), np.array(y_coords), np.array(strain_values)

def plot_strain_slice(x_coords, y_coords, strain_values, strain_component, cmap='viridis', 
                     title=None, fig_size=(10, 8), dpi=300, output_file=None):
    """
    Create a contour plot of the strain field at a specific z-slice with publication-quality formatting.
    
    Parameters:
    -----------
    x_coords : numpy.ndarray
        X-coordinates of nodes
    y_coords : numpy.ndarray
        Y-coordinates of nodes
    strain_values : numpy.ndarray
        Strain values at the nodes
    strain_component : str
        Strain component being plotted
    cmap : str
        Colormap name (default: 'viridis')
    title : str
        Plot title (default: None)
    fig_size : tuple
        Figure size (width, height) in inches
    dpi : int
        Resolution in dots per inch
    output_file : str
        Path to save the figure (default: None)
    """
    # Apply publication style settings
    set_publication_style()
    
    # Create a figure
    fig, ax = plt.subplots(figsize=fig_size)
    
    # Create a triangulation for irregular data points
    from matplotlib.tri import Triangulation
    triang = Triangulation(x_coords, y_coords)
    
    # Create the contour plot
    contour = ax.tricontourf(triang, strain_values, cmap=cmap, levels=20)
    
    # Add a color bar
    cbar = plt.colorbar(contour, ax=ax)
    cbar.set_label('Logarithmic Strain (%s)' % strain_component)
    
    # Set labels and title
    ax.set_xlabel('X Coordinate (mm)')
    ax.set_ylabel('Y Coordinate (mm)')
    
    if title:
        ax.set_title(title)
    else:
        ax.set_title('Z-Slice Contour Plot of %s Strain' % strain_component)
    
    # Equal aspect ratio
    ax.set_aspect('equal')
    
    # Save figure if requested
    if output_file:
        plt.savefig(output_file, bbox_inches='tight')
        print("Figure saved to: %s" % output_file)
    
    # Display the plot
    plt.tight_layout()
    plt.show()

def print_usage():
    """Print usage information."""
    print(__doc__)

def main():
    """Main function to extract and plot strain data from an Abaqus ODB file."""
    # Check if we have enough command line arguments
    if len(sys.argv) < 2:
        print_usage()
        sys.exit(1)
    
    # Get command line arguments
    odb_file = sys.argv[1]
    step_name = sys.argv[2] if len(sys.argv) > 2 else "Step-1"
    
    try:
        frame_number = int(sys.argv[3]) if len(sys.argv) > 3 else -1
    except ValueError:
        print("Invalid frame number. Using last frame (-1).")
        frame_number = -1
    
    try:
        z_coord = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    except ValueError:
        print("Invalid z-coordinate. Using z=0.0")
        z_coord = 0.0
    
    strain_component = sys.argv[5] if len(sys.argv) > 5 else "LEMAG"
    cmap = sys.argv[6] if len(sys.argv) > 6 else "viridis"
    
    # Print the parameters being used
    print("\nParameters:")
    print("  ODB File: %s" % odb_file)
    print("  Step Name: %s" % step_name)
    print("  Frame Number: %d" % frame_number)
    print("  Z-Coordinate: %f" % z_coord)
    print("  Strain Component: %s" % strain_component)
    print("  Colormap: %s\n" % cmap)
    
    # Extract data from the ODB file
    x_coords, y_coords, strain_values = extract_le_strain_z_slice(
        odb_file, step_name, frame_number, z_coord, strain_component
    )
    
    # Check if we got any data
    if len(x_coords) == 0:
        print("No nodes found at z-coordinate %f. Please try a different value." % z_coord)
        return
    
    # Plot the data
    title = "Logarithmic Strain (%s) at Z = %f" % (strain_component, z_coord)
    
    # Create output filename
    filename = os.path.basename(odb_file)
    base_name = os.path.splitext(filename)[0]
    output_file = "%s_%s_z%.3f.png" % (base_name, strain_component, z_coord)
    
    plot_strain_slice(
        x_coords, y_coords, strain_values, 
        strain_component, cmap=cmap,
        title=title, output_file=output_file
    )
    
    print("Data statistics:")
    print("  Number of nodes: %d" % len(strain_values))
    print("  Min strain: %.6f" % min(strain_values))
    print("  Max strain: %.6f" % max(strain_values))
    print("  Mean strain: %.6f" % (sum(strain_values)/len(strain_values)))

if __name__ == "__main__":
    main()