"""
Spinal Cord Node Set Creation Script
=====================================
Creates node sets with sinusoidal predefined field distribution for
modelling compression sites in spinal cord injury simulations.

PURPOSE
-------
To apply a spatially-varying predefined field (e.g., swelling) to a spinal
cord model, where the field intensity is highest at the compression center
and tapers off sinusoidally towards the edges.

WHAT THIS SCRIPT DOES
---------------------
1. DEFINE COMPRESSION REGION
   - Takes 3 coordinate points: center, upper limit, lower limit
   - These define an axis along the spinal cord where the field will be applied
   - Center point = location of maximum compression/swelling (field = 0.15)
   - Upper/lower points = boundaries where field tapers to minimum (field = 0.01)

2. SELECT NODES FROM INSTANCE
   - Retrieves all nodes from the specified part instance
   - No need to pre-create a node set - uses all nodes in the instance

3. CLASSIFY NODES INTO BANDS
   - Projects each node onto the compression axis
   - Assigns nodes to one of 5 bands based on distance from center
   - Nodes outside the upper/lower limits are excluded

4. CALCULATE SINUSOIDAL FIELD VALUES
   - Band 1 (center): 0.150 - peak field value
   - Band 2: 0.128
   - Band 3: 0.075
   - Band 4: 0.022
   - Band 5 (edge): 0.000 - minimum field value
   - Uses raised cosine profile (zero gradient at center and edges)

5. CREATE NODE SETS IN ABAQUS
   - Creates assembly-level node sets named <prefix>_BAND_1 through _BAND_5
   - Outputs summary table with node counts and field values
   - Saves the model database

USAGE
-----
Run in Abaqus CAE kernel:
    os.chdir('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting')
    execfile('sets_script.py')

Option 1 - Single site (set variables before execfile):
os.chdir('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting')
MODEL_NAME = 'N01-015_2026-02-09-scripting'
INSTANCE_NAME = 'PART-1_1-1'
CENTER_POINT = (-501.574E-03,-169.124924,-174.925827)
UPPER_POINT = (-887.162E-03,-168.287567,-166.848785)
LOWER_POINT = (-377.941E-03,-163.577835,-180.216278)
execfile('sets_script_old.py')

Option 2 - Multiple sites from CSV file:
    MODEL_NAME = 'Model-1'
    INSTANCE_NAME = 'PART-1_1-1'
    COORDS_FILE = 'coordinates.csv'
    execfile('sets_script.py')

CSV format (coordinates.csv):
    site_name,center_x,center_y,center_z,upper_x,upper_y,upper_z,lower_x,lower_y,lower_z
    Site1,0.0,0.0,0.0,0.0,0.0,10.0,0.0,0.0,-10.0
    Site2,5.0,0.0,5.0,5.0,0.0,15.0,5.0,0.0,-5.0

Python 2.7 compatible (Abaqus Python)
"""

from abaqus import *
from abaqusConstants import *
import math


def read_coordinates_file(filepath):
    """
    Read compression site coordinates from a CSV file.

    CSV format:
        site_name,center_x,center_y,center_z,upper_x,upper_y,upper_z,lower_x,lower_y,lower_z

    Parameters:
        filepath: path to the CSV file

    Returns:
        list of dicts: [{'name': str, 'center': tuple, 'upper': tuple, 'lower': tuple}, ...]
    """
    sites = []
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Skip header line
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue

        parts = line.split(',')
        if len(parts) < 10:
            print("Warning: Skipping invalid line: {}".format(line))
            continue

        site = {
            'name': parts[0].strip(),
            'center': (float(parts[1]), float(parts[2]), float(parts[3])),
            'upper': (float(parts[4]), float(parts[5]), float(parts[6])),
            'lower': (float(parts[7]), float(parts[8]), float(parts[9]))
        }
        sites.append(site)

    return sites


def dot_product(v1, v2):
    """Calculate dot product of two 3D vectors."""
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]


def vector_magnitude(v):
    """Calculate magnitude of a 3D vector."""
    return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)


def normalize_vector(v):
    """Normalize a 3D vector to unit length."""
    mag = vector_magnitude(v)
    if mag == 0:
        return (0, 0, 0)
    return (v[0]/mag, v[1]/mag, v[2]/mag)


def calculate_axis_projection(point, axis_vector, axis_origin):
    """
    Project a point onto an axis and return signed distance from origin along axis.

    Parameters:
        point: (x, y, z) coordinates of the point
        axis_vector: normalized (dx, dy, dz) direction vector of the axis
        axis_origin: (x, y, z) point on the axis (used as origin for distance)

    Returns:
        float: signed distance along the axis from origin
    """
    # Vector from axis origin to point
    to_point = (
        point[0] - axis_origin[0],
        point[1] - axis_origin[1],
        point[2] - axis_origin[2]
    )
    # Project onto axis (dot product with unit axis vector)
    return dot_product(to_point, axis_vector)


def get_band_index(distance, d_upper, d_lower, num_bands=5):
    """
    Determine which band (0 to num_bands-1) a node belongs to based on distance from center.

    Band 0 is the center band (highest field value).
    Band num_bands-1 is the edge band (lowest field value).

    Parameters:
        distance: signed distance from center along axis
        d_upper: distance from center to upper limit (positive)
        d_lower: distance from center to lower limit (negative)
        num_bands: number of bands to create

    Returns:
        int: band index (0 to num_bands-1), or -1 if outside region
    """
    # Check if node is within the region
    if distance > d_upper or distance < d_lower:
        return -1

    # Normalize distance to [0, 1] range based on which side of center
    if distance >= 0:
        # Upper half
        if d_upper == 0:
            normalized = 0
        else:
            normalized = abs(distance) / abs(d_upper)
    else:
        # Lower half
        if d_lower == 0:
            normalized = 0
        else:
            normalized = abs(distance) / abs(d_lower)

    # Map normalized distance to band index
    # Band 0: 0-20%, Band 1: 20-40%, etc.
    band_width = 1.0 / num_bands
    band_index = int(normalized / band_width)

    # Clamp to valid range (handle edge case where normalized == 1.0)
    if band_index >= num_bands:
        band_index = num_bands - 1

    return band_index


def calculate_field_values(num_bands=5, peak_value=0.15, min_value=0.0):
    """
    Calculate sinusoidal field values for each band.

    Uses half-cosine profile scaled between peak_value (center) and min_value (edge).

    Parameters:
        num_bands: number of bands
        peak_value: field value at center (band 0)
        min_value: field value at edge (band num_bands-1)

    Returns:
        list: field values for each band
    """
    field_values = []
    amplitude = (peak_value - min_value)

    for i in range(num_bands):
        # Position as fraction from center (0) to edge (1)
        # Band 1 = 0.0, Band 2 = 0.25, ..., Band 5 = 1.0
        if num_bands == 1:
            band_position = 0.0
        else:
            band_position = i / (num_bands - 1.0)

        # Raised cosine: gradient = 0 at both center and edge (smooth step)
        field_value = amplitude * (1 + math.cos(math.pi * band_position)) / 2.0 + min_value
        field_values.append(round(field_value, 3))

    return field_values


def create_sinusoidal_node_sets(
    center_point,
    upper_point,
    lower_point,
    num_bands=5,
    peak_field_value=0.15,
    min_field_value=0.01,
    model_name='Model-1',
    instance_name='PART-1_1-1',
    set_prefix='FIELD_BAND',
    site_index=1
):
    """
    Create node sets with sinusoidal field distribution along an axis.

    Selects all nodes from the specified instance and classifies them into bands.

    Parameters:
        center_point: (x1, y1, z1) center of the region (peak field value)
        upper_point: (x2, y2, z2) upper limit of the region
        lower_point: (x3, y3, z3) lower limit of the region
        num_bands: number of bands to create (default 5)
        peak_field_value: field value at center (default 0.15)
        min_field_value: field value at edges (default 0.01)
        model_name: name of the Abaqus model (default 'Model-1')
        instance_name: name of the part instance (default 'PART-1_1-1')
        set_prefix: prefix for node set names (default 'FIELD_BAND')

    Returns:
        dict: mapping of set names to field values
    """
    # Get the model database
    try:
        modelDB = mdb.models[model_name]
    except KeyError:
        print("Error: Model '{}' not found.".format(model_name))
        print("Available models: {}".format(list(mdb.models.keys())))
        return None
    assembly = modelDB.rootAssembly
    
    
    print('assembly loaded')
    

    # Calculate axis vector (from lower to upper)
    axis_vector = (
        upper_point[0] - lower_point[0],
        upper_point[1] - lower_point[1],
        upper_point[2] - lower_point[2]
    )
    axis_vector = normalize_vector(axis_vector)

    # print('axis vector calculated:')
    # print(axis_vector)

    

    # Calculate distances from center to upper and lower limits
    d_upper = calculate_axis_projection(upper_point, axis_vector, center_point)
    d_lower = calculate_axis_projection(lower_point, axis_vector, center_point)

    print("Axis vector: ({:.3f}, {:.3f}, {:.3f})".format(*axis_vector))
    print("Distance to upper limit: {:.3f}".format(d_upper))
    print("Distance to lower limit: {:.3f}".format(d_lower))

    

    # Get all nodes from the instance
    try:
        instance = assembly.instances[instance_name]
        all_nodes = instance.nodes
        print("Found {} nodes in instance '{}'".format(len(all_nodes), instance_name))
    except KeyError:
        print("Error: Instance '{}' not found in assembly.".format(instance_name))
        print("Available instances: {}".format(list(assembly.instances.keys())))
        return None

    

    # Initialize lists for each band
    band_nodes = [[] for _ in range(num_bands)]
    nodes_outside = 0

    # Classify each node into a band
    print("Classifying nodes into {} bands...".format(num_bands))
    for node in all_nodes:
        coords = node.coordinates
        distance = calculate_axis_projection(coords, axis_vector, center_point)
        band_idx = get_band_index(distance, d_upper, d_lower, num_bands)

        if band_idx >= 0:
            band_nodes[band_idx].append(node.label)
        else:
            nodes_outside += 1

    print("Nodes outside region: {}".format(nodes_outside))

    

    # Calculate field values for each band
    field_values = calculate_field_values(num_bands, peak_field_value, min_field_value)
    # print('field values calculated')
    

    # Create node sets for each band
    set_field_mapping = {}

    for i in range(num_bands):
        set_name = "{}_{:d}".format(set_prefix, i + 1)
        node_labels = band_nodes[i]

        if len(node_labels) > 0:
            # Create the node set
            # OLD: creates set under instance subcategory
            # assembly.SetFromNodeLabels(
            #     name=set_name,
            #     nodeLabels=((instance_name, node_labels),)
            # )
            # NEW: creates set at assembly level
            node_sequence = instance.nodes.sequenceFromLabels(node_labels)
            assembly.Set(name=set_name, nodes=node_sequence)
            set_field_mapping[set_name] = field_values[i]
            print("Created '{}': {} nodes, field value = {:.3f}".format(
                set_name, len(node_labels), field_values[i]))

            # Create predefined field for this node set
            field_name = "predefinedfield-{:d}-fieldband{:d}".format(site_index, i + 1)
            region = assembly.sets[set_name]
            modelDB.Temperature(
                name=field_name,
                createStepName='Initial',
                region=region,
                magnitudes=(field_values[i],)
            )
            print("Created predefined field '{}' = {:.3f}".format(
                field_name, field_values[i]))
        else:
            print("Warning: Band {} has no nodes, skipping.".format(i + 1))

    # Print summary table
    print("\n" + "="*50)
    print("SUMMARY: Node Sets and Field Values")
    print("="*50)
    print("{:<20} {:>10} {:>15} {:<30}".format("Node Set", "Nodes", "Field Value", "Predefined Field"))
    print("-"*75)
    for i in range(num_bands):
        set_name = "{}_{:d}".format(set_prefix, i + 1)
        field_name = "predefinedfield-{:d}-fieldband{:d}".format(site_index, i + 1)
        node_count = len(band_nodes[i])
        print("{:<20} {:>10} {:>15.3f} {:<30}".format(set_name, node_count, field_values[i], field_name))
    print("="*75)
     
    # Save the model
    # print("\nSaving model database...")
    # mdb.save()
    # print("Save complete.")

    return set_field_mapping


# =============================================================================
# SCRIPT STARTS HERE
# =============================================================================
#
# Usage: Set variables in Abaqus kernel, then execfile()
#
#   MODEL_NAME = 'Model-1'
#   INSTANCE_NAME = 'PART-1_1-1'
#   CENTER_POINT = (0.0, 0.0, 0.0)
#   UPPER_POINT = (0.0, 0.0, 10.0)
#   LOWER_POINT = (0.0, 0.0, -10.0)
#   execfile('sets_script.py')
#
# =============================================================================

print("sets_script.py loaded")



# --- Check variables are set ---
if 'MODEL_NAME' not in dir():
    print("")
    print("ERROR: MODEL_NAME not set. Usage:")
    print("")
    print("  Option 1 - Single site:")
    print("    MODEL_NAME = 'Model-1'")
    print("    INSTANCE_NAME = 'PART-1_1-1'")
    print("    CENTER_POINT = (0.0, 0.0, 0.0)")
    print("    UPPER_POINT = (0.0, 0.0, 10.0)")
    print("    LOWER_POINT = (0.0, 0.0, -10.0)")
    print("    execfile('sets_script.py')")
    print("")
    print("  Option 2 - Multiple sites from CSV:")
    print("    MODEL_NAME = 'Model-1'")
    print("    INSTANCE_NAME = 'PART-1_1-1'")
    print("    COORDS_FILE = 'coordinates.csv'")
    print("    execfile('sets_script.py')")
    print("")

elif 'COORDS_FILE' in dir():
    # --- Multiple sites from CSV file ---
    print("Model: " + MODEL_NAME)
    print("Instance: " + INSTANCE_NAME)
    print("Coordinates file: " + COORDS_FILE)

    sites = read_coordinates_file(COORDS_FILE)
    print("Found {} sites in CSV file".format(len(sites)))

    for idx, site in enumerate(sites):
        site_num = idx + 1
        print("\n" + "#"*60)
        print("SITE {}: {}".format(site_num, site['name']))
        print("#"*60)
        print("Center: {}".format(site['center']))
        print("Upper: {}".format(site['upper']))
        print("Lower: {}".format(site['lower']))

        set_prefix = "{}_BAND".format(site['name'].upper().replace(' ', '_'))

        create_sinusoidal_node_sets(
            center_point=site['center'],
            upper_point=site['upper'],
            lower_point=site['lower'],
            model_name=MODEL_NAME,
            instance_name=INSTANCE_NAME,
            set_prefix=set_prefix,
            site_index=site_num
        )

    print("\nAll {} sites processed.".format(len(sites)))

else:
    # --- Single site from variables ---
    print("Model: " + MODEL_NAME)
    print("Instance: " + INSTANCE_NAME)
    print("Center: {}".format(CENTER_POINT))
    print("Upper: {}".format(UPPER_POINT))
    print("Lower: {}".format(LOWER_POINT))

    create_sinusoidal_node_sets(
        center_point=CENTER_POINT,
        upper_point=UPPER_POINT,
        lower_point=LOWER_POINT,
        model_name=MODEL_NAME,
        instance_name=INSTANCE_NAME,
        set_prefix='FIELD_BAND',
        site_index=1
    )

    print("Done")

    
    