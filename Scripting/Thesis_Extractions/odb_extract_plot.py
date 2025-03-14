import os
import sys
import pickle
import gzip 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D

def set_publication_style():
	"""Set matplotlib parameters for publication-quality figures."""
	plt.rcParams.update({
		'font.family': 'serif',
		'font.serif': ['Times New Roman'],
		'mathtext.fontset': 'stix',
		'font.size': 12,
		'axes.labelsize': 14,
		'axes.titlesize': 16,
		'axes.titleweight': 'bold', # This makes titles bold
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

# for plots
set_publication_style()

# Pickle basepath
print("running...")
base_path=Path("Scripting/Thesis_Extractions/Output")
output_dir =base_path / "nofalxnohemi0.5.pkl"
model_structure_path = output_dir / "model_structure.pkl.gz"
step_1_dir = output_dir / "Step-1"
step_2_dir = output_dir / "Step-2"

if step_1_dir.exists():
    # Dictionary to store frames from each step
    step_1_frames = {}
    print("Step 1 exists")

if step_2_dir.exists():
    # Dictionary to store frames from each step
    step_2_frames = {}
    print("Step 2 exists")
# list contents of output dir

# unpickle the pickles
# unpickle model_structure.pkl.gz
print('unpickling...')
with gzip.open(model_structure_path, 'rb') as f:
    try:
        model_structure = pickle.load(f)
    except:
        f.seek(0)  # Reset file pointer to beginning
        model_structure = pickle.load(f, encoding='latin1')
    print("Successfully loaded model structure")

# Explore the model structure
print("\nModel Structure Keys:", model_structure.keys())

# Load frames from Step-1 and Step-2
step_1_frames = {}
step_2_frames = {}

# Function to load frames from a directory
def load_frames_from_dir(directory):
    frames = {}
    if directory.exists():
        # List all frame files in the directory that match the pattern
        frame_files = sorted(directory.glob("frame_00??.pkl.gz"))
        
        for frame_path in frame_files:
            frame_name = frame_path.name
            frame_number = int(frame_name.split("_")[1].split(".")[0])
            
            # Load the frame
            with gzip.open(frame_path, 'rb') as f:
                try:
                    frame_data = pickle.load(f)
                except:
                    f.seek(0)
                    frame_data = pickle.load(f, encoding='latin1')
            
            # Store the frame in the dictionary
            frames[frame_number] = frame_data
            print(f"Loaded frame {frame_number:04d}")
    
    print(f"Loaded {len(frames)} frames from {directory}")
    return frames

# Load frames from both directories
if step_1_dir.exists():
    print("\nLoading frames from Step-1...")
    step_1_frames = load_frames_from_dir(step_1_dir)
    
    # Print keys of the first frame to understand structure
    if step_1_frames:
        first_frame = list(step_1_frames.values())[0]
        print("Frame Keys:", first_frame.keys())
        
        # Look at field outputs structure
        if 'field_outputs' in first_frame and 'LE' in first_frame['field_outputs']:
            print("LE Field keys:", first_frame['field_outputs']['LE'].keys())
            
            # Check the values structure
            values = first_frame['field_outputs']['LE']['values']
            print(f"Values type: {type(values)}, shape or length: {len(values) if hasattr(values, '__len__') else 'unknown'}")

if step_2_dir.exists():
    print("\nLoading frames from Step-2...")
    step_2_frames = load_frames_from_dir(step_2_dir)

### plt model structure - NODES ONLY
# Debugging function to check node structure

def debug_node_structure(model_structure, max_nodes=5):
    """
    Debug function to examine the structure of node_coordinates
    """
    nodes = model_structure['node_coordinates']
    
    print("\n=== DEBUGGING NODE STRUCTURE ===")
    print(f"Total nodes in node_coordinates: {len(nodes)}")
    
    # Look at a few sample nodes
    sample_keys = list(nodes.keys())[:max_nodes]
    for key in sample_keys:
        print(f"Node {key} structure: {type(nodes[key])}")
        if isinstance(nodes[key], dict) and 'coordinates' in nodes[key]:
            print(f"Node {key} coordinates: {nodes[key]['coordinates']}")
        else:
            print(f"Node {key} value: {nodes[key]}")

#visualize_complete_model(model_structure, element_set_name='GM_CB', max_elements=500000)
def visualize_large_model(model_structure, element_set_name='GM_CB', sample_ratio=0.05):
    """
    Optimized visualization for very large models
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import random
    
    # Collect coordinates as before, but with random sampling
    coordinates = []  # Will store [x, y, z] coordinates
    
    elements = model_structure['element_connectivity']
    node_coordinates = model_structure['node_coordinates']
    element_set = model_structure['element_sets'].get(element_set_name, [])
    
    # Group elements by part
    part_elements = {}
    for part_name, elem_id in element_set:
        if part_name not in part_elements:
            part_elements[part_name] = []
        part_elements[part_name].append(elem_id)
    
    # Sample elements randomly
    element_count = 0
    for part_name, elem_ids in part_elements.items():
        if part_name in elements and part_name in node_coordinates:
            part_connectivity = elements[part_name]
            part_nodes = node_coordinates[part_name]
            
            # Randomly sample elements based on the sample ratio
            sample_size = int(len(elem_ids) * sample_ratio)
            if sample_size < 1:
                sample_size = 1
            
            sampled_elem_ids = random.sample(elem_ids, min(sample_size, len(elem_ids)))
            
            for elem_id in sampled_elem_ids:
                element_count += 1
                
                if elem_id in part_connectivity:
                    node_ids = part_connectivity[elem_id]['connectivity']
                    
                    # Take just one node per element to reduce point count
                    if node_ids:
                        node_id = node_ids[0]  # Just take the first node
                        
                        if node_id in part_nodes:
                            if isinstance(part_nodes[node_id], dict) and 'coordinates' in part_nodes[node_id]:
                                coords = part_nodes[node_id]['coordinates']
                                if len(coords) == 3:
                                    coordinates.append([float(coords[0]), float(coords[1]), float(coords[2])])
    
    print(f"Sampled {len(coordinates)} coordinates from {element_count} elements")
    
    if coordinates:
        # Convert to numpy array for efficient processing
        coords_array = np.array(coordinates)
        
        # Setup basic plot with reduced resolution for faster rendering
        plt.figure(figsize=(10, 8), dpi=100)
        ax = plt.subplot(111, projection='3d')
        
        # Downsample points if still too many
        max_points = 20000  # Keep plot manageable
        if len(coords_array) > max_points:
            indices = np.random.choice(len(coords_array), max_points, replace=False)
            coords_array = coords_array[indices]
            print(f"Downsampled to {len(coords_array)} points for visualization")
        
        # Plot with smaller markers for better performance
        ax.scatter(
            coords_array[:, 0], 
            coords_array[:, 1], 
            coords_array[:, 2],
            c='blue',
            marker='.',
            s=5,
            alpha=0.3
        )
        
        # Calculate center for proper view
        center = np.mean(coords_array, axis=0)
        print(f"Center of mass: {center}")
        
        # Set reasonable axis limits
        ax.set_xlim(center[0] - 100, center[0] + 100)
        ax.set_ylim(center[1] - 100, center[1] + 100)
        ax.set_zlim(center[2] - 100, center[2] + 100)
        
        # Set labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        # Set title
        plt.title(f'Model Structure (Sampled) - {element_set_name}')
        
        # Save with reasonable quality
        plt.tight_layout()
        plt.savefig('model_sampled.png', dpi=150)
        print("Plot saved as model_sampled.png")
        
        return True
    else:
        print("No coordinates found for plotting")
        return False
# Run the optimized visualization
#visualize_large_model(model_structure, sample_ratio=0.01)  # Use 1% sampling ratio

def visualize_multi_angle(model_structure, element_set_name='GM_CB', sample_ratio=0.01):
    """
    Visualization with multiple viewing angles
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import random
    
    # Collect coordinates with sampling
    coordinates = []
    
    elements = model_structure['element_connectivity']
    node_coordinates = model_structure['node_coordinates']
    element_set = model_structure['element_sets'].get(element_set_name, [])
    
    # Group elements by part
    part_elements = {}
    for part_name, elem_id in element_set:
        if part_name not in part_elements:
            part_elements[part_name] = []
        part_elements[part_name].append(elem_id)
    
    # Sample elements
    element_count = 0
    for part_name, elem_ids in part_elements.items():
        if part_name in elements and part_name in node_coordinates:
            part_connectivity = elements[part_name]
            part_nodes = node_coordinates[part_name]
            
            sample_size = int(len(elem_ids) * sample_ratio)
            if sample_size < 1:
                sample_size = 1
            
            sampled_elem_ids = random.sample(elem_ids, min(sample_size, len(elem_ids)))
            
            for elem_id in sampled_elem_ids:
                element_count += 1
                
                if elem_id in part_connectivity:
                    node_ids = part_connectivity[elem_id]['connectivity']
                    
                    # Take just one node per element
                    if node_ids:
                        node_id = node_ids[0]
                        
                        if node_id in part_nodes:
                            if isinstance(part_nodes[node_id], dict) and 'coordinates' in part_nodes[node_id]:
                                coords = part_nodes[node_id]['coordinates']
                                if len(coords) == 3:
                                    coordinates.append([float(coords[0]), float(coords[1]), float(coords[2])])
    
    print(f"Sampled {len(coordinates)} coordinates from {element_count} elements")
    
    if coordinates:
        # Convert to numpy array
        coords_array = np.array(coordinates)
        
        # Downsample if needed
        max_points = 20000
        if len(coords_array) > max_points:
            indices = np.random.choice(len(coords_array), max_points, replace=False)
            coords_array = coords_array[indices]
            print(f"Downsampled to {len(coords_array)} points for visualization")
        
        # Calculate center for proper view
        center = np.mean(coords_array, axis=0)
        print(f"Center of mass: {center}")
        
        # Define viewing angles (elev, azim) pairs
        viewing_angles = [
            (30, 0),    # Front view
            (30, 90),   # Side view
            (30, 180),  # Back view
            (30, 270),  # Other side view
            (90, 0),    # Top view
            (0, 0)      # Horizontal view
        ]
        
        # Create plots for each angle
        for i, (elev, azim) in enumerate(viewing_angles):
            # Create a new figure for each angle
            plt.figure(figsize=(10, 8), dpi=100)
            ax = plt.subplot(111, projection='3d')
            
            # Plot points
            ax.scatter(
                coords_array[:, 0], 
                coords_array[:, 1], 
                coords_array[:, 2],
                c='blue',
                marker='.',
                s=5,
                alpha=0.3
            )
            
            # Set view angle
            ax.view_init(elev=elev, azim=azim)
            
            # Set axis limits
            radius = 100  # Adjust based on your model size
            ax.set_xlim(center[0] - radius, center[0] + radius)
            ax.set_ylim(center[1] - radius, center[1] + radius)
            ax.set_zlim(center[2] - radius, center[2] + radius)
            
            # Set labels
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            
            # Set title with angle information
            plt.title(f'Model - {element_set_name} (Elev: {elev}°, Azim: {azim}°)')
            
            # Save figure
            plt.tight_layout()
            plt.savefig(f'model_angle_{i+1}_elev{elev}_azim{azim}.png', dpi=150)
            print(f"Saved view {i+1}: Elevation {elev}°, Azimuth {azim}°")
            
            # Close figure to free memory
            plt.close()
        
        print("All angle views saved successfully")
        return True
    else:
        print("No coordinates found for plotting")
        return False

# visualise LE data for the last frame of step 1
def visualize_LE_data(model_structure, frame_data, element_set_name='GM_CB', sample_ratio=0.01, component=0):
    """
    Visualize the LE strain data for a given frame, mapped to the model coordinates
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import random
    
    # Extract field data from the frame
    field_outputs = frame_data.get('field_outputs', {})
    LE_data = field_outputs.get('LE', {})
    
    if not LE_data or 'values' not in LE_data:
        print("No valid LE data found in the frame")
        return False
    
    # Get the strain values and component labels
    strain_values = LE_data['values']
    component_labels = LE_data.get('component_labels', [f'Component {i}' for i in range(6)])
    
    # Debug: Look at the strain values structure
    if len(strain_values) > 0:
        print(f"First strain value type: {type(strain_values[0])}")
        print(f"First strain value: {strain_values[0]}")
    
    # Extract model structure
    elements = model_structure['element_connectivity']
    node_coordinates = model_structure['node_coordinates']
    element_set = model_structure['element_sets'].get(element_set_name, [])
    
    # Group elements by part
    part_elements = {}
    for part_name, elem_id in element_set:
        if part_name not in part_elements:
            part_elements[part_name] = []
        part_elements[part_name].append(elem_id)
    
    # Collect coordinates and strain values
    points = []  # Will store [x, y, z, strain_value]
    
    # Sample elements
    element_count = 0
    for part_name, elem_ids in part_elements.items():
        if part_name in elements and part_name in node_coordinates:
            part_connectivity = elements[part_name]
            part_nodes = node_coordinates[part_name]
            
            # Random sampling
            sample_size = int(len(elem_ids) * sample_ratio)
            if sample_size < 1:
                sample_size = 1
            
            sampled_elem_ids = random.sample(elem_ids, min(sample_size, len(elem_ids)))
            
            for elem_id in sampled_elem_ids:
                element_count += 1
                
                if elem_id in part_connectivity:
                    # Get node coordinates
                    node_ids = part_connectivity[elem_id]['connectivity']
                    
                    # Use the first node for visualization
                    if node_ids and node_ids[0] in part_nodes:
                        node_data = part_nodes[node_ids[0]]
                        
                        if isinstance(node_data, dict) and 'coordinates' in node_data:
                            coords = node_data['coordinates']
                            
                            # Try to get a strain value for this element
                            try:
                                # This is just a heuristic - the exact mapping depends on your data
                                strain_idx = (elem_id - 1) % len(strain_values)
                                strain_value = strain_values[strain_idx]
                                
                                # Handle different types of strain values
                                if isinstance(strain_value, dict):
                                    # Extract from dictionary
                                    if 'value' in strain_value:
                                        strain_value = strain_value['value']
                                    elif 'data' in strain_value:
                                        strain_value = strain_value['data']
                                
                                # Handle numpy arrays
                                if isinstance(strain_value, np.ndarray):
                                    # If it's an array, select the component 
                                    if component < strain_value.size:
                                        strain_value = strain_value[component]
                                    else:
                                        strain_value = strain_value.mean()  # Fallback to mean
                                
                                # Final check to ensure we have a scalar value
                                strain_value = float(strain_value)
                                
                                # Store the point
                                points.append([
                                    float(coords[0]), 
                                    float(coords[1]), 
                                    float(coords[2]), 
                                    strain_value
                                ])
                            except (TypeError, ValueError, IndexError) as e:
                                # Skip this point if there's an error
                                continue
    
    print(f"Collected {len(points)} points with strain data from {element_count} elements")
    
    if points:
        # Convert to numpy array
        points_array = np.array(points)
        
        # Extract coordinates and strain values
        coordinates = points_array[:, :3]
        strains = points_array[:, 3]
        
        # Downsample if needed
        max_points = 20000
        if len(coordinates) > max_points:
            indices = np.random.choice(len(coordinates), max_points, replace=False)
            coordinates = coordinates[indices]
            strains = strains[indices]
            print(f"Downsampled to {len(coordinates)} points for visualization")
        
        # Calculate bounds for strain values
        vmin = np.min(strains)
        vmax = np.max(strains)
        print(f"Strain range: {vmin} to {vmax}")
        
        # Setup 3D plot
        fig = plt.figure(figsize=(10, 8), dpi=100)
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot the points
        scatter = ax.scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            coordinates[:, 2],
            c=strains,
            cmap='viridis',
            marker='.',
            s=5,
            alpha=0.3,
            vmin=vmin,
            vmax=vmax
        )
        
        # Calculate center for proper view
        center = np.mean(coordinates, axis=0)
        
        # Set axis limits
        radius = 100  # Adjust based on your model size
        ax.set_xlim(center[0] - radius, center[0] + radius)
        ax.set_ylim(center[1] - radius, center[1] + radius)
        ax.set_zlim(center[2] - radius, center[2] + radius)
        
        # Set labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        # Get component name
        component_name = component_labels[component] if component < len(component_labels) else f"Component {component}"
        
        # Set title
        plt.title(f'Strain Visualization - {component_name}')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax, orientation='vertical')
        cbar.set_label(component_name)
        
        # Save figure
        plt.tight_layout()
        plt.savefig(f'strain_{component_name.replace(" ", "_")}.png', dpi=150)
        print(f"Strain visualization saved as strain_{component_name.replace(' ', '_')}.png")
        
        return True
    else:
        print("No valid points with strain data found")
        return False

def visualize_LE_abaqus_style(model_structure, frame_data, element_set_name='GM_CB', sample_ratio=0.05, component=0):
    """
    Create an Abaqus-like visualization with smooth element surfaces using the working strain extraction approach
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import matplotlib.colors as colors
    import random
    
    # Extract field data
    field_outputs = frame_data.get('field_outputs', {})
    LE_data = field_outputs.get('LE', {})
    
    if not LE_data or 'values' not in LE_data:
        print("No valid LE data found in the frame")
        return False
    
    # Get strain values and component labels
    strain_values = LE_data['values']
    component_labels = LE_data.get('component_labels', [f'Component {i}' for i in range(6)])
    
    # Extract model structure
    elements = model_structure['element_connectivity']
    node_coordinates = model_structure['node_coordinates']
    element_set = model_structure['element_sets'].get(element_set_name, [])
    
    # Group elements by part
    part_elements = {}
    for part_name, elem_id in element_set:
        if part_name not in part_elements:
            part_elements[part_name] = []
        part_elements[part_name].append(elem_id)
    
    # Collection for faces and colors
    faces = []
    face_colors = []
    
    # Track element types for debugging
    element_types = set()
    
    # Process each part
    for part_name, elem_ids in part_elements.items():
        if part_name in elements and part_name in node_coordinates:
            part_connectivity = elements[part_name]
            part_nodes = node_coordinates[part_name]
            
            # Sample elements for performance
            sample_size = int(len(elem_ids) * sample_ratio)
            if sample_size < 1:
                sample_size = 1
            
            sampled_elem_ids = random.sample(elem_ids, min(sample_size, len(elem_ids)))
            print(f"Processing {len(sampled_elem_ids)} elements from part {part_name}")
            
            # Process each element
            for elem_id in sampled_elem_ids:
                if elem_id in part_connectivity:
                    elem_data = part_connectivity[elem_id]
                    node_ids = elem_data['connectivity']
                    elem_type = elem_data.get('type', '')
                    
                    # Track element types
                    element_types.add(elem_type)
                    
                    # Get strain value for this element using the working approach
                    try:
                        # This is just a heuristic - the exact mapping depends on your data
                        strain_idx = (elem_id - 1) % len(strain_values)
                        strain_value = strain_values[strain_idx]
                        
                        # Handle different types of strain values
                        if isinstance(strain_value, dict):
                            # Extract from dictionary
                            if 'value' in strain_value:
                                strain_value = strain_value['value']
                            elif 'data' in strain_value:
                                strain_value = strain_value['data']
                        
                        # Handle numpy arrays
                        if isinstance(strain_value, np.ndarray):
                            # If it's an array, select the component 
                            if component < strain_value.size:
                                strain_value = strain_value[component]
                            else:
                                strain_value = strain_value.mean()  # Fallback to mean
                        
                        # Final check to ensure we have a scalar value
                        strain_value = float(strain_value)
                    except (TypeError, ValueError, IndexError) as e:
                        # Skip this element if there's an error with strain extraction
                        continue
                    
                    # Get node coordinates
                    coords = []
                    for node_id in node_ids:
                        if node_id in part_nodes:
                            node_data = part_nodes[node_id]
                            if isinstance(node_data, dict) and 'coordinates' in node_data:
                                node_coords = node_data['coordinates']
                                coords.append([float(c) for c in node_coords])
                    
                    if len(coords) >= 3:  # Need at least 3 nodes for a face
                        # Create faces based on node count
                        if len(coords) == 4:  # Tetrahedral element
                            # Create 4 triangular faces
                            faces.append([coords[0], coords[1], coords[2]])
                            face_colors.append(strain_value)
                            faces.append([coords[0], coords[1], coords[3]])
                            face_colors.append(strain_value)
                            faces.append([coords[1], coords[2], coords[3]])
                            face_colors.append(strain_value)
                            faces.append([coords[0], coords[2], coords[3]])
                            face_colors.append(strain_value)
                        elif len(coords) == 8:  # Hexahedral element
                            # Create 6 quadrilateral faces
                            faces.append([coords[0], coords[1], coords[2], coords[3]])
                            face_colors.append(strain_value)
                            faces.append([coords[4], coords[5], coords[6], coords[7]])
                            face_colors.append(strain_value)
                            faces.append([coords[0], coords[1], coords[5], coords[4]])
                            face_colors.append(strain_value)
                            faces.append([coords[1], coords[2], coords[6], coords[5]])
                            face_colors.append(strain_value)
                            faces.append([coords[2], coords[3], coords[7], coords[6]])
                            face_colors.append(strain_value)
                            faces.append([coords[3], coords[0], coords[4], coords[7]])
                            face_colors.append(strain_value)
                        elif len(coords) == 10:  # 10-node tetrahedral (C3D10M)
                            # Use the first 4 nodes to create faces
                            faces.append([coords[0], coords[1], coords[2]])
                            face_colors.append(strain_value)
                            faces.append([coords[0], coords[1], coords[3]])
                            face_colors.append(strain_value)
                            faces.append([coords[1], coords[2], coords[3]])
                            face_colors.append(strain_value)
                            faces.append([coords[0], coords[2], coords[3]])
                            face_colors.append(strain_value)
                        else:
                            # For unknown element types, create a simple triangular face
                            faces.append([coords[0], coords[1], coords[2]])
                            face_colors.append(strain_value)
    
    print(f"Element types found: {element_types}")
    print(f"Created {len(faces)} faces for visualization")
    
    if faces:
        # Setup plot
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Create the 3D collection of polygons
        poly = Poly3DCollection(faces, alpha=0.7)
        
        # Calculate color range
        face_colors = np.array(face_colors)
        vmin = np.min(face_colors)
        vmax = np.max(face_colors)
        print(f"Strain range: {vmin} to {vmax}")
        
        # Normalize colors
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        
        # Set face colors
        poly.set_facecolor(plt.cm.viridis(norm(face_colors)))
        
        # Add the collection to the plot
        ax.add_collection3d(poly)
        
        # Find the overall bounds of the model for proper display
        all_coords = np.array([coord for face in faces for coord in face])
        x_min, y_min, z_min = np.min(all_coords, axis=0)
        x_max, y_max, z_max = np.max(all_coords, axis=0)
        
        # Set axis limits
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_zlim(z_min, z_max)
        
        # Add a colorbar
        m = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
        m.set_array([])
        cbar = plt.colorbar(m, ax=ax)
        
        component_name = component_labels[component] if component < len(component_labels) else f"Component {component}"
        cbar.set_label(component_name)
        
        # Set labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        # Set title
        plt.title(f'Strain Visualization - {component_name}')
        
        # Set good viewing angle
        ax.view_init(elev=30, azim=45)
        
        # Save figure
        plt.tight_layout()
        plt.savefig(f'strain_surface_{component_name.replace(" ", "_")}.png', dpi=300)
        print(f"Surface visualization saved as strain_surface_{component_name.replace(' ', '_')}.png")
        
        # Save alternate views
        for angle in [0, 90, 180, 270]:
            ax.view_init(elev=30, azim=angle)
            plt.savefig(f'strain_surface_{component_name.replace(" ", "_")}_angle{angle}.png', dpi=300)
            print(f"Saved view from angle {angle}°")
        
        return True
    else:
        print("No valid faces created for visualization")
        return False

# Visualize the first component (LE11) for the first frame
if step_1_frames and 0 in step_1_frames:
    visualize_LE_data(model_structure, step_1_frames[0], component=0)  # 0 = LE11
    visualize_LE_abaqus_style(model_structure, step_1_frames[0], element_set_name='GM_CB', sample_ratio=0.05, component=0)

