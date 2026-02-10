"""
Spinal Cord Node Set Creation Script - Direct .inp File Modification
=====================================================================
Creates node sets and predefined temperature fields for modelling
compression sites in spinal cord injury simulations by writing
directly into the Abaqus .inp file.

This script exists because creating assembly-level sets via the
Abaqus Python API (assembly.Set) corrupts the assembly tree into
sub-assemblies, causing the .inp writer to silently drop the sets
and predefined fields.

USAGE
-----
Run in Abaqus CAE kernel:
    os.chdir('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting')
    execfile('sets_inpmod_script.py')

Option 1 - Single site:
    MODEL_NAME = 'Model-1'
    INSTANCE_NAME = 'PART-1_1-1'
    INP_FILE = 'Job-212.inp'
    CENTER_POINT = (x, y, z)
    UPPER_POINT = (x, y, z)
    LOWER_POINT = (x, y, z)
    execfile('sets_inpmod_script.py')

Option 2 - Multiple sites from CSV:
    MODEL_NAME = 'Model-1'
    INSTANCE_NAME = 'PART-1_1-1'
    INP_FILE = 'Job-212.inp'
    COORDS_FILE = 'coordinates.csv'
    execfile('sets_inpmod_script.py')

CSV format (coordinates.csv):
    site_name,center_x,center_y,center_z,upper_x,upper_y,upper_z,lower_x,lower_y,lower_z
    Site1,0.0,0.0,0.0,0.0,0.0,10.0,0.0,0.0,-10.0

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
    Project a point onto an axis and return displacement from origin along axis.

    Parameters:
        point: (x, y, z) coordinates of the point
        axis_vector: normalized (dx, dy, dz) direction vector of the axis
        axis_origin: (x, y, z) point on the axis (used as origin for distance)

    Returns:
        float: displacement along the axis from origin
    """
    to_point = (
        point[0] - axis_origin[0],
        point[1] - axis_origin[1],
        point[2] - axis_origin[2]
    )
    return dot_product(to_point, axis_vector)


def get_band_index(distance, d_upper, d_lower, num_bands=5):
    """
    Determine which band (0 to num_bands-1) a node belongs to based on distance from center.

    Band 0 is the center band (highest field value).
    Band num_bands-1 is the edge band (lowest field value).

    Parameters:
        distance: displacement from center along axis
        d_upper: distance from center to upper limit (positive)
        d_lower: distance from center to lower limit (negative)
        num_bands: number of bands to create

    Returns:
        int: band index (0 to num_bands-1), or -1 if outside region
    """
    if distance > d_upper or distance < d_lower:
        return -1

    if distance >= 0:
        if d_upper == 0:
            normalized = 0
        else:
            normalized = abs(distance) / abs(d_upper)
    else:
        if d_lower == 0:
            normalized = 0
        else:
            normalized = abs(distance) / abs(d_lower)

    band_width = 1.0 / num_bands
    band_index = int(normalized / band_width)

    if band_index >= num_bands:
        band_index = num_bands - 1

    return band_index


def calculate_field_values(num_bands=5, peak_value=0.15, min_value=0.0):
    """
    Calculate raised cosine field values for each band.

    Uses raised cosine profile (smooth step) that reaches zero at 120% of the
    region extent. Only num_bands values are returned (the virtual point at 120%
    where the field reaches zero is not included).

    Parameters:
        num_bands: number of bands to create
        peak_value: field value at center (band 1)
        min_value: field value at the virtual zero-crossing beyond the last band

    Returns:
        list: field values for each band
    """
    field_values = []
    amplitude = (peak_value - min_value)
    num_points = num_bands + 1

    for i in range(num_bands):
        band_position = i / (num_points - 1.0)
        field_value = amplitude * (1 + math.cos(math.pi * band_position)) / 2.0 + min_value
        field_values.append(round(field_value, 3))

    return field_values


def format_node_labels(node_labels, per_line=16):
    """
    Format a list of node labels into .inp format lines (16 per line).

    Parameters:
        node_labels: list of integer node labels
        per_line: number of labels per line (default 16)

    Returns:
        str: formatted lines
    """
    lines = []
    for start in range(0, len(node_labels), per_line):
        chunk = node_labels[start:start + per_line]
        lines.append(' ' + ', '.join(str(n) for n in chunk))
    return '\n'.join(lines)


def write_sets_to_inp(inp_file, instance_name, band_nodes, field_values,
                      set_prefix, site_index, amplitude_name, num_bands):
    """
    Write node sets and predefined temperature fields directly into an .inp file.

    Inserts *Nset blocks before *End Assembly, and *Temperature blocks
    after ** PREDEFINED FIELDS in the step section.

    Parameters:
        inp_file: path to the .inp file
        instance_name: name of the part instance
        band_nodes: list of lists of node labels per band
        field_values: list of field values per band
        set_prefix: prefix for node set names
        site_index: site number for predefined field naming
        amplitude_name: amplitude name for *Temperature keyword
        num_bands: number of bands
    """
    with open(inp_file, 'r') as f:
        content = f.read()

    lines = content.split('\n')

    # Build *Nset blocks
    nset_blocks = []
    for i in range(num_bands):
        set_name = "{}_{:d}".format(set_prefix, i + 1)
        node_labels = band_nodes[i]
        if len(node_labels) > 0:
            nset_blocks.append("*Nset, nset={}, instance={}".format(set_name, instance_name))
            nset_blocks.append(format_node_labels(node_labels))

    # Build *Temperature blocks
    temp_blocks = []
    for i in range(num_bands):
        set_name = "{}_{:d}".format(set_prefix, i + 1)
        node_labels = band_nodes[i]
        if len(node_labels) > 0:
            field_name = "predefinedfield-{:d}-fieldband{:d}".format(site_index, i + 1)
            temp_blocks.append("** Name: {}   Type: Temperature".format(field_name))
            temp_blocks.append("*Temperature, amplitude={}".format(amplitude_name))
            temp_blocks.append("{}, {}".format(set_name, field_values[i]))

    # Check for duplicates
    first_set = "{}_{:d}".format(set_prefix, 1)
    if first_set in content:
        print("WARNING: '{}' already exists in .inp file. Skipping to avoid duplicates.".format(first_set))
        return

    # Find insertion points and build new content
    new_lines = []
    nsets_inserted = False
    temps_inserted = False

    for line in lines:
        # Insert *Nset blocks before *End Assembly
        if not nsets_inserted and line.strip() == '*End Assembly':
            for block_line in nset_blocks:
                new_lines.append(block_line)
            nsets_inserted = True

        new_lines.append(line)

        # Insert *Temperature blocks after ** PREDEFINED FIELDS
        if not temps_inserted and line.strip() == '** PREDEFINED FIELDS':
            # Skip the blank line after ** PREDEFINED FIELDS
            # (it will be added by the next iteration)
            for block_line in temp_blocks:
                new_lines.append(block_line)
            temps_inserted = True

    if not nsets_inserted:
        print("WARNING: Could not find '*End Assembly' in .inp file")
        print("Node sets were NOT written")
    if not temps_inserted:
        print("WARNING: Could not find '** PREDEFINED FIELDS' in .inp file")
        print("Predefined fields were NOT written")

    with open(inp_file, 'w') as f:
        f.write('\n'.join(new_lines))

    print("Written to {}".format(inp_file))
    num_written = len([b for b in band_nodes if len(b) > 0])
    if nsets_inserted:
        print("  - {} node set(s) inserted before *End Assembly".format(num_written))
    if temps_inserted:
        print("  - {} predefined field(s) inserted after ** PREDEFINED FIELDS".format(num_written))


def create_sinusoidal_node_sets(
    center_point,
    upper_point,
    lower_point,
    inp_file,
    num_bands=5,
    peak_field_value=0.15,
    min_field_value=0.0,
    model_name='Model-1',
    instance_name='PART-1_1-1',
    set_prefix='FIELD_BAND',
    site_index=1,
    amplitude_name='Amp-1-preload'
):
    """
    Classify nodes into bands and write sets/predefined fields to .inp file.

    Reads node coordinates from the Abaqus model to classify nodes, then
    writes *Nset and *Temperature keywords directly into the .inp file.
    Does NOT create any CAE objects (no assembly.Set or modelDB.Temperature).

    Parameters:
        center_point: (x, y, z) center of the region (peak field value)
        upper_point: (x, y, z) upper limit of the region
        lower_point: (x, y, z) lower limit of the region
        inp_file: path to the .inp file to modify
        num_bands: number of bands to create (default 5)
        peak_field_value: field value at center (default 0.15)
        min_field_value: field value at edges (default 0.0)
        model_name: name of the Abaqus model (default 'Model-1')
        instance_name: name of the part instance (default 'PART-1_1-1')
        set_prefix: prefix for node set names (default 'FIELD_BAND')
        site_index: site number for predefined field naming (default 1)
        amplitude_name: amplitude name for *Temperature (default 'Amp-1-preload')

    Returns:
        dict: mapping of set names to field values
    """
    # Get the model database (read-only, for node coordinates)
    try:
        modelDB = mdb.models[model_name]
    except KeyError:
        print("Error: Model '{}' not found.".format(model_name))
        print("Available models: {}".format(list(mdb.models.keys())))
        return None
    assembly = modelDB.rootAssembly

    # Calculate axis vector (from lower to upper)
    axis_vector = (
        upper_point[0] - lower_point[0],
        upper_point[1] - lower_point[1],
        upper_point[2] - lower_point[2]
    )
    axis_vector = normalize_vector(axis_vector)

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

    # Write to .inp file
    print("Writing to .inp file: {}".format(inp_file))
    write_sets_to_inp(
        inp_file=inp_file,
        instance_name=instance_name,
        band_nodes=band_nodes,
        field_values=field_values,
        set_prefix=set_prefix,
        site_index=site_index,
        amplitude_name=amplitude_name,
        num_bands=num_bands
    )

    # Print summary table
    set_field_mapping = {}
    print("\n" + "=" * 75)
    print("SUMMARY: Node Sets and Field Values")
    print("=" * 75)
    print("{:<20} {:>10} {:>15} {:<30}".format("Node Set", "Nodes", "Field Value", "Predefined Field"))
    print("-" * 75)
    for i in range(num_bands):
        set_name = "{}_{:d}".format(set_prefix, i + 1)
        field_name = "predefinedfield-{:d}-fieldband{:d}".format(site_index, i + 1)
        node_count = len(band_nodes[i])
        print("{:<20} {:>10} {:>15.3f} {:<30}".format(set_name, node_count, field_values[i], field_name))
        if node_count > 0:
            set_field_mapping[set_name] = field_values[i]
    print("=" * 75)

    return set_field_mapping


# =============================================================================
# SCRIPT STARTS HERE
# =============================================================================

print("sets_inpmod_script.py loaded")

# --- Check variables are set ---
if 'MODEL_NAME' not in dir():
    print("")
    print("ERROR: MODEL_NAME not set. Usage:")
    print("")
    print("  Option 1 - Single site:")
    print("    MODEL_NAME = 'Model-1'")
    print("    INSTANCE_NAME = 'PART-1_1-1'")
    print("    INP_FILE = 'Job-212.inp'")
    print("    CENTER_POINT = (0.0, 0.0, 0.0)")
    print("    UPPER_POINT = (0.0, 0.0, 10.0)")
    print("    LOWER_POINT = (0.0, 0.0, -10.0)")
    print("    execfile('sets_inpmod_script.py')")
    print("")
    print("  Option 2 - Multiple sites from CSV:")
    print("    MODEL_NAME = 'Model-1'")
    print("    INSTANCE_NAME = 'PART-1_1-1'")
    print("    INP_FILE = 'Job-212.inp'")
    print("    COORDS_FILE = 'coordinates.csv'")
    print("    execfile('sets_inpmod_script.py')")
    print("")

elif 'INP_FILE' not in dir():
    print("")
    print("ERROR: INP_FILE not set.")
    print("    INP_FILE = 'Job-212.inp'")
    print("")

elif 'COORDS_FILE' in dir():
    # --- Multiple sites from CSV file ---
    print("Model: " + MODEL_NAME)
    print("Instance: " + INSTANCE_NAME)
    print("INP file: " + INP_FILE)
    print("Coordinates file: " + COORDS_FILE)

    sites = read_coordinates_file(COORDS_FILE)
    print("Found {} sites in CSV file".format(len(sites)))

    for idx, site in enumerate(sites):
        site_num = idx + 1
        print("\n" + "#" * 60)
        print("SITE {}: {}".format(site_num, site['name']))
        print("#" * 60)
        print("Center: {}".format(site['center']))
        print("Upper: {}".format(site['upper']))
        print("Lower: {}".format(site['lower']))

        set_prefix = "{}_BAND".format(site['name'].upper().replace(' ', '_'))

        create_sinusoidal_node_sets(
            center_point=site['center'],
            upper_point=site['upper'],
            lower_point=site['lower'],
            inp_file=INP_FILE,
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
    print("INP file: " + INP_FILE)
    print("Center: {}".format(CENTER_POINT))
    print("Upper: {}".format(UPPER_POINT))
    print("Lower: {}".format(LOWER_POINT))

    create_sinusoidal_node_sets(
        center_point=CENTER_POINT,
        upper_point=UPPER_POINT,
        lower_point=LOWER_POINT,
        inp_file=INP_FILE,
        model_name=MODEL_NAME,
        instance_name=INSTANCE_NAME,
        set_prefix='FIELD_BAND',
        site_index=1
    )

    print("Done")
