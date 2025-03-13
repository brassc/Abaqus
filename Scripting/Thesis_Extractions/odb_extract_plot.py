import os
import sys
import pickle
import gzip 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D




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

"""
if step_1_dir.exists():
    # Dictionary to store all frames
    # List all frame files in the directory that match the pattern
    frame_files = sorted(step_1_dir.glob("frame_00??.pkl.gz"))

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
        step_1_frames[frame_number] = frame_data
        print(f"Loaded frame {frame_number:04d}")

    print(f"Loaded {len(step_1_frames)} frames from {step_1_dir}")


# Process Step-2 if it exists
if step_2_dir.exists():
    # List all frame files in the directory that match the pattern
    frame_files = sorted(step_2_dir.glob("frame_00??.pkl.gz"))
    
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
        step_2_frames[frame_number] = frame_data
        print(f"Loaded frame {frame_number:04d} from Step-2")
    
    print(f"Loaded {len(step_2_frames)} frames from {step_2_dir}")
else:
    print(f"Step 2 directory doesn't exist: {step_2_dir}")

# Extract the data from the frames
# print the keys of the model_structure
print("Model Structure Keys:")
print(model_structure.keys())
# print the keys of the first frame
print("Frame Keys:")
print(step_1_frames[0].keys())
# print the keys of the second frame
print("Frame Keys:")
print(step_1_frames[1].keys())
# print the keys of the third frame
print("Frame Keys:")
print(step_1_frames[2].keys())
"""



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

### plt model structure 
# Function to plot the model structure correctly
def plot_model_structure(model_structure, element_set_name='GM_CB', title=None):
    """
    Plot the model structure for the specified element set.
    
    Args:
        model_structure: Model structure dictionary
        element_set_name: Name of the element set to plot
        title: Plot title
    """
    # Extract node coordinates
    nodes = model_structure['node_coordinates']
    elements = model_structure['element_connectivity']
    
    # Setup the figure
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # The element_sets structure is different than expected
    # It seems to be a dictionary of (part_name, element_id) tuples
    part_elements = {}
    element_set = model_structure['element_sets'].get(element_set_name, [])
    
    # Group elements by part
    for part_name, elem_id in element_set:
        if part_name not in part_elements:
            part_elements[part_name] = []
        part_elements[part_name].append(elem_id)
    
    print(f"Found {len(part_elements)} parts in element set {element_set_name}")
    for part_name, elem_ids in part_elements.items():
        print(f"  Part {part_name}: {len(elem_ids)} elements")
    
    # Get all nodes for visualization
    all_nodes = set()
    for part_name, elem_ids in part_elements.items():
        if part_name in elements:
            part_connectivity = elements[part_name]
            for elem_id in elem_ids:
                if elem_id in part_connectivity:
                    node_ids = part_connectivity[elem_id]['connectivity']
                    all_nodes.update(node_ids)
    
    print(f"Total unique nodes: {len(all_nodes)}")
    
    # Plot nodes
    if all_nodes:
        print("Plotting nodes...")
        x, y, z = [], [], []
        for node_id in all_nodes:
            if node_id in nodes:
                node_coords = nodes[node_id]
                x.append(node_coords[0])
                y.append(node_coords[1])
                z.append(node_coords[2])
        
        ax.scatter(x, y, z, c='blue', marker='x', alpha=0.1, s=10)
        
        # Set equal aspect ratio to keep visualization proportional
        ax.set_box_aspect([1, 1, 1])
        
        # Set labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        if title:
            plt.title(title)
        else:
            plt.title(f'Model Structure - {element_set_name}')
        
        plt.tight_layout()
        plt.savefig('model_structure.png')
        print("Plot saved as model_structure.png")
        
        return fig, ax
    else:
        print("No nodes found for plotting")
        return None, None

# Plot the model structure
#plot_model_structure(model_structure)

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

# Modified function to plot model structure
def plot_model_structure_fixed(model_structure, element_set_name='GM_CB', title=None):
    """
    Plot the model structure with fixes for the nested node structure
    """
    # Extract node coordinates
    nodes = model_structure['node_coordinates']
    elements = model_structure['element_connectivity']
    
    # Setup the figure
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # The element_sets structure is as expected
    part_elements = {}
    element_set = model_structure['element_sets'].get(element_set_name, [])
    
    # Group elements by part
    for part_name, elem_id in element_set:
        if part_name not in part_elements:
            part_elements[part_name] = []
        part_elements[part_name].append(elem_id)
    
    print(f"Found {len(part_elements)} parts in element set {element_set_name}")
    
    # Get all nodes for visualization
    all_nodes = set()
    for part_name, elem_ids in part_elements.items():
        if part_name in elements:
            part_connectivity = elements[part_name]
            for elem_id in elem_ids:
                if elem_id in part_connectivity:
                    node_ids = part_connectivity[elem_id]['connectivity']
                    all_nodes.update(node_ids)
    
    print(f"Total unique nodes: {len(all_nodes)}")
    
    # Plot nodes - account for nested node structure
    if all_nodes:
        x, y, z = [], [], []
        valid_nodes = 0
        
        for node_id in all_nodes:
            if node_id in nodes:
                try:
                    # Extract coordinates from nested dictionary
                    if isinstance(nodes[node_id], dict) and 'coordinates' in nodes[node_id]:
                        node_coords = nodes[node_id]['coordinates']
                    else:
                        node_coords = nodes[node_id]
                        
                    # Ensure we have valid 3D coordinates
                    if len(node_coords) == 3 and all(isinstance(coord, (int, float, np.number)) for coord in node_coords):
                        x.append(float(node_coords[0]))
                        y.append(float(node_coords[1]))
                        z.append(float(node_coords[2]))
                        valid_nodes += 1
                except Exception as e:
                    print(f"Error with node {node_id}: {e}")
                    continue
        
        print(f"Found {valid_nodes} valid nodes with coordinates for plotting")
        
        if valid_nodes > 0:
            scatter = ax.scatter(x, y, z, c='blue', marker='x', alpha=0.5, s=10)
            
            # Set equal aspect ratio to keep visualization proportional
            ax.set_box_aspect([1, 1, 1])
            
            # Set labels
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            
            if title:
                plt.title(title)
            else:
                plt.title(f'Model Structure - {element_set_name}')
            
            plt.tight_layout()
            plt.savefig('model_structure_fixed.png')
            print("Plot saved as model_structure_fixed.png")
            
            return fig, ax
        else:
            print("No valid nodes found for plotting")
            return None, None
    else:
        print("No nodes found for plotting")
        return None, None

# First check the node structure
#debug_node_structure(model_structure)

# Then try the fixed plotting function
#plot_model_structure_fixed(model_structure)
# Updated function to properly extract node coordinates
# More direct visualization approach with explicit control

# Try the direct visualization approach
#visualize_model_direct(model_structure)

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
visualize_large_model(model_structure, sample_ratio=0.01)  # Use 1% sampling ratio

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

# Run the multi-angle visualization
visualize_multi_angle(model_structure)

sys.exit()

# Assume we've already loaded the data
# model_structure contains node coordinates and element connectivity
# step_1_frames and step_2_frames contain the field outputs (LE) for each frame

# First, let's inspect the structure of the field outputs
def explore_frame_data(frame_data):
    """Explore the structure of frame data to find where the field values are stored"""
    print("\nExploring frame data structure:")
    print(f"Top level keys: {frame_data.keys()}")
    
    if 'field_outputs' in frame_data:
        print(f"Field output keys: {frame_data['field_outputs'].keys()}")
        
        if 'LE' in frame_data['field_outputs']:
            print(f"LE structure: {frame_data['field_outputs']['LE'].keys()}")
            
            # Print a sample of the first field output structure
            print("\nSample of LE field output structure:")
            for key, value in frame_data['field_outputs']['LE'].items():
                if isinstance(value, dict):
                    print(f"  {key}: {list(value.keys())}")
                elif isinstance(value, list) or isinstance(value, np.ndarray):
                    print(f"  {key}: Array of shape {np.array(value).shape}")
                else:
                    print(f"  {key}: {value}")

# Example usage:
frame_number = 0
if frame_number in step_1_frames:
    explore_frame_data(step_1_frames[frame_number])

# First, let's deeply inspect the model structure
def detailed_model_inspection(model_structure):
    print("\n===== DETAILED MODEL INSPECTION =====")
    
    # Check top level keys
    print(f"Model structure keys: {model_structure.keys()}")
    
    # Check node coordinates
    node_count = len(model_structure['node_coordinates'])
    print(f"Number of nodes: {node_count}")
    if node_count > 0:
        sample_node_id = list(model_structure['node_coordinates'].keys())[0]
        print(f"Sample node {sample_node_id}: {model_structure['node_coordinates'][sample_node_id]}")
    
    # Check element connectivity
    element_count = len(model_structure['element_connectivity'])
    print(f"Number of elements: {element_count}")
    if element_count > 0:
        sample_elem_id = list(model_structure['element_connectivity'].keys())[0]
        print(f"Sample element {sample_elem_id}: {model_structure['element_connectivity'][sample_elem_id]}")
    
    # Check element sets
    print(f"Element sets: {list(model_structure['element_sets'].keys())}")
    for set_name, elements in model_structure['element_sets'].items():
        print(f"  Set {set_name}: {len(elements)} elements")
        if len(elements) > 0:
            print(f"    First few elements: {elements[:5]}")
            
            # Check if these element IDs actually exist in element_connectivity
            found = 0
            for elem_id in elements[:5]:
                if elem_id in model_structure['element_connectivity']:
                    found += 1
            print(f"    Elements found in connectivity: {found}/5")
    
    # Try to understand the structure of element_connectivity
    if isinstance(model_structure['element_connectivity'], dict):
        print("Element connectivity is a dictionary")
        keys = list(model_structure['element_connectivity'].keys())
        print(f"  First few keys: {keys[:5]}")
        print(f"  Key type: {type(keys[0]) if keys else 'unknown'}")
    else:
        print(f"Element connectivity is a {type(model_structure['element_connectivity'])}")
    
    # See if the GM_CB element set contains valid elements
    if 'GM_CB' in model_structure['element_sets']:
        valid_count = 0
        for elem_id in model_structure['element_sets']['GM_CB']:
            if elem_id in model_structure['element_connectivity']:
                valid_count += 1
        print(f"Valid elements in GM_CB set: {valid_count}/{len(model_structure['element_sets']['GM_CB'])}")

# Run detailed inspection
detailed_model_inspection(model_structure)

