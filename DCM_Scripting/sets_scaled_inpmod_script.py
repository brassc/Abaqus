"""
Spinal Cord Node Set Creation Script - With Cord-Diameter-Based Preload Scaling
================================================================================
Creates node sets and predefined temperature fields for modelling
compression sites in spinal cord injury simulations by writing
directly into the Abaqus .inp file.

Extends sets_inpmod_script.py with automatic per-site preload scaling
based on sagittal cord diameter measurements. The peak field value
(preload) at each site is scaled relative to a known reference
calibration point using the compression fraction:

    compression_fraction     = (upper_cord_sag_dist - indent_cord_sag_dist) / upper_cord_sag_dist
    ref_compression_fraction = (REF_UPPER_CORD_SAG_DIST - REF_INDENT_CORD_SAG_DIST) / REF_UPPER_CORD_SAG_DIST
    peak_field_value         = REFERENCE_PRELOAD * (compression_fraction / ref_compression_fraction)

    upper_cord_sag_dist  - sagittal cord AP diameter at unaffected level (mm)
    indent_cord_sag_dist - sagittal cord AP diameter at the compression/indent site (mm)

Default reference calibration (override by setting variables before execfile):
    REFERENCE_PRELOAD        = 0.5      preload that worked at the reference site
    REF_UPPER_CORD_SAG_DIST  = 6.82259  upper_cord_sag_dist at the reference site (mm)
    REF_INDENT_CORD_SAG_DIST = 4.24591  indent_cord_sag_dist at the reference site (mm)

USAGE
-----
Run in Abaqus CAE kernel:
    os.chdir('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting')
    execfile('sets_scaled_inpmod_script.py')

Option 1 - Single site (manual PEAK_FIELD_VALUE, no scaling):
    MODEL_NAME = 'Model-1'
    INSTANCE_NAME = 'PART-1_1-1'
    INP_FILE = 'Job-212.inp'
    PEAK_FIELD_VALUE = 0.5
    CENTER_POINT = (x, y, z)
    UPPER_POINT = (x, y, z)
    LOWER_POINT = (x, y, z)
    execfile('sets_scaled_inpmod_script.py')

Option 2 - Multiple sites from CSV (scaling applied automatically):
    MODEL_NAME = 'Model-1'
    INSTANCE_NAME = 'PART-1_1-1'
    INP_FILE = 'Job-212.inp'
    COORDS_FILE = 'coordinates.csv'
    execfile('sets_scaled_inpmod_script.py')

    Optionally override reference calibration before execfile:
    REFERENCE_PRELOAD        = 0.5
    REF_UPPER_CORD_SAG_DIST  = 6.82259
    REF_INDENT_CORD_SAG_DIST = 4.24591

CSV format (coordinates.csv):
    site_name,center_x,center_y,center_z,upper_x,upper_y,upper_z,lower_x,lower_y,lower_z,upper_cord_sag_dist,indent_cord_sag_dist
    Site1,0.0,0.0,0.0,0.0,0.0,10.0,0.0,0.0,-10.0,6.82259,4.24591

    Columns upper_cord_sag_dist and indent_cord_sag_dist are optional.
    If absent, falls back to PEAK_FIELD_VALUE (default 0.15).

Python 2.7 compatible (Abaqus Python)
"""

from abaqus import *
from abaqusConstants import *
import math
import os


def read_coordinates_file(filepath):
    """
    Read compression site coordinates from a CSV file.

    CSV format:
        site_name,center_x,center_y,center_z,upper_x,upper_y,upper_z,lower_x,lower_y,lower_z[,upper_cord_sag_dist,indent_cord_sag_dist]

    Parameters:
        filepath: path to the CSV file

    Returns:
        list of dicts: [{'name': str, 'center': tuple, 'upper': tuple, 'lower': tuple,
                         'upper_cord_sag_dist': float or None, 'indent_cord_sag_dist': float or None}, ...]
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
            'lower': (float(parts[7]), float(parts[8]), float(parts[9])),
            'upper_cord_sag_dist': None,
            'indent_cord_sag_dist': None
        }

        if len(parts) >= 12:
            try:
                site['upper_cord_sag_dist'] = float(parts[10])
                site['indent_cord_sag_dist'] = float(parts[11])
            except ValueError:
                print("Warning: Could not parse cord distances for site '{}'. Falling back to PEAK_FIELD_VALUE.".format(site['name']))

        sites.append(site)

    return sites


def compute_scaled_peak_value(upper_cord_sag_dist, indent_cord_sag_dist,
                               reference_preload, ref_upper_cord_sag_dist, ref_indent_cord_sag_dist):
    """
    Compute scaled peak field value from sagittal cord diameter measurements.

    Scales the preload proportionally to the compression fraction at this site
    relative to the reference calibration site:

        compression_fraction     = (upper - indent) / upper
        ref_compression_fraction = (ref_upper - ref_indent) / ref_upper
        peak_value               = reference_preload * (compression_fraction / ref_compression_fraction)

    Parameters:
        upper_cord_sag_dist      - cord AP diameter at unaffected level (mm)
        indent_cord_sag_dist     - cord AP diameter at the indent site (mm)
        reference_preload        - preload used at the reference calibration site
        ref_upper_cord_sag_dist  - upper_cord_sag_dist at the reference site (mm)
        ref_indent_cord_sag_dist - indent_cord_sag_dist at the reference site (mm)

    Returns:
        float: scaled peak field value
    """
    compression_fraction = (upper_cord_sag_dist - indent_cord_sag_dist) / upper_cord_sag_dist
    ref_compression_fraction = (ref_upper_cord_sag_dist - ref_indent_cord_sag_dist) / ref_upper_cord_sag_dist
    return reference_preload * (compression_fraction / ref_compression_fraction)


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


def generate_output_filename(inp_file):
    """
    Generate output filename from input filename.

    Example: 'Job-212.inp' -> 'Job-212_modified.inp'

    Parameters:
        inp_file: path to the input .inp file

    Returns:
        str: path to the output .inp file
    """
    base, ext = os.path.splitext(inp_file)
    return "{}_modified{}".format(base, ext)


def write_sets_to_inp(inp_file, out_file, instance_name, band_nodes, field_values,
                      set_prefix, site_index, amplitude_name, num_bands):
    """
    Write node sets and predefined temperature fields into a new .inp file.

    Reads from inp_file and writes the modified content to out_file.
    Inserts *Nset blocks before *End Assembly. For *Temperature blocks:
    - If ** PREDEFINED FIELDS exists, inserts after it
    - If not, creates a new ** PREDEFINED FIELDS section in Step-1
      (before the first section header, e.g. ** BOUNDARY CONDITIONS)

    Parameters:
        inp_file: path to the input .inp file (read only)
        out_file: path to the output .inp file (written)
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

    # Determine insertion strategy for temperature blocks
    has_predef_fields = any(l.strip() == '** PREDEFINED FIELDS' for l in lines)

    # Find insertion points and build new content
    new_lines = []
    nsets_inserted = False
    temps_inserted = False
    in_step_1 = False

    for line in lines:
        stripped = line.strip()

        # Insert *Nset blocks before *End Assembly
        if not nsets_inserted and stripped == '*End Assembly':
            for block_line in nset_blocks:
                new_lines.append(block_line)
            nsets_inserted = True

        # If no ** PREDEFINED FIELDS exists, insert a new section
        # before the first section header in Step-1
        if not temps_inserted and not has_predef_fields:
            if stripped.startswith('*Step,') and 'name=Step-1' in stripped:
                in_step_1 = True
            elif in_step_1 and stripped.startswith('** ') and len(stripped) > 3:
                # Check if this is a section header (all-uppercase text after "** ")
                header_text = stripped[3:]
                if header_text == header_text.upper() and header_text[0].isalpha():
                    new_lines.append('** ')
                    new_lines.append('** PREDEFINED FIELDS')
                    for block_line in temp_blocks:
                        new_lines.append(block_line)
                    temps_inserted = True
                    print("  - Created new ** PREDEFINED FIELDS section in Step-1")
            elif in_step_1 and stripped == '*End Step':
                # Fallback: no section headers found, insert before *End Step
                new_lines.append('** ')
                new_lines.append('** PREDEFINED FIELDS')
                for block_line in temp_blocks:
                    new_lines.append(block_line)
                temps_inserted = True
                print("  - Created new ** PREDEFINED FIELDS section before *End Step")

        new_lines.append(line)

        # If ** PREDEFINED FIELDS exists, insert temperature blocks after it
        if not temps_inserted and has_predef_fields and stripped == '** PREDEFINED FIELDS':
            for block_line in temp_blocks:
                new_lines.append(block_line)
            temps_inserted = True

    if not nsets_inserted:
        print("WARNING: Could not find '*End Assembly' in .inp file")
        print("Node sets were NOT written")
    if not temps_inserted:
        print("WARNING: Could not find Step-1 in .inp file")
        print("Predefined fields were NOT written")

    with open(out_file, 'w') as f:
        f.write('\n'.join(new_lines))

    print("Written to {}".format(out_file))
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
    out_file,
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
    writes *Nset and *Temperature keywords into a new output .inp file.
    The original input file is not modified.
    Does NOT create any CAE objects (no assembly.Set or modelDB.Temperature).

    Parameters:
        center_point: (x, y, z) center of the region (peak field value)
        upper_point: (x, y, z) upper limit of the region
        lower_point: (x, y, z) lower limit of the region
        inp_file: path to the input .inp file (read only)
        out_file: path to the output .inp file (written)
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

    # Write to output .inp file
    print("Reading from: {}".format(inp_file))
    print("Writing to: {}".format(out_file))
    write_sets_to_inp(
        inp_file=inp_file,
        out_file=out_file,
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

print("sets_scaled_inpmod_script.py loaded")

# --- Reference calibration constants ---
# Override any of these before execfile() if your reference site differs.
if 'REFERENCE_PRELOAD' not in dir():
    REFERENCE_PRELOAD = 0.3
if 'REF_UPPER_CORD_SAG_DIST' not in dir():
    REF_UPPER_CORD_SAG_DIST = 6.82259
if 'REF_INDENT_CORD_SAG_DIST' not in dir():
    REF_INDENT_CORD_SAG_DIST = 4.24591

REF_COMPRESSION_FRACTION = (REF_UPPER_CORD_SAG_DIST - REF_INDENT_CORD_SAG_DIST) / REF_UPPER_CORD_SAG_DIST
print("Reference preload:            {}".format(REFERENCE_PRELOAD))
print("Ref upper cord sag dist:      {} mm".format(REF_UPPER_CORD_SAG_DIST))
print("Ref indent cord sag dist:     {} mm".format(REF_INDENT_CORD_SAG_DIST))
print("Ref compression fraction:     {:.4f}".format(REF_COMPRESSION_FRACTION))

# --- Check variables are set ---
if 'MODEL_NAME' not in dir():
    print("")
    print("ERROR: MODEL_NAME not set. Usage:")
    print("")
    print("  Option 1 - Single site (manual peak value, no scaling):")
    print("    MODEL_NAME = 'Model-1'")
    print("    INSTANCE_NAME = 'PART-1_1-1'")
    print("    INP_FILE = 'Job-212.inp'")
    print("    PEAK_FIELD_VALUE = 0.5")
    print("    CENTER_POINT = (0.0, 0.0, 0.0)")
    print("    UPPER_POINT = (0.0, 0.0, 10.0)")
    print("    LOWER_POINT = (0.0, 0.0, -10.0)")
    print("    execfile('sets_scaled_inpmod_script.py')")
    print("")
    print("  Option 2 - Multiple sites from CSV (auto-scaled):")
    print("    MODEL_NAME = 'Model-1'")
    print("    INSTANCE_NAME = 'PART-1_1-1'")
    print("    INP_FILE = 'Job-212.inp'")
    print("    COORDS_FILE = 'coordinates.csv'")
    print("    execfile('sets_scaled_inpmod_script.py')")
    print("")
    print("  CSV columns:")
    print("    site_name,center_x,center_y,center_z,upper_x,upper_y,upper_z,")
    print("    lower_x,lower_y,lower_z,upper_cord_sag_dist,indent_cord_sag_dist")
    print("")

elif 'INP_FILE' not in dir():
    print("")
    print("ERROR: INP_FILE not set.")
    print("    INP_FILE = 'Job-212.inp'")
    print("")

elif 'COORDS_FILE' in dir():
    # --- Multiple sites from CSV file ---
    if 'PEAK_FIELD_VALUE' not in dir():
        PEAK_FIELD_VALUE = REFERENCE_PRELOAD
    OUT_FILE = generate_output_filename(INP_FILE)
    print("Model: " + MODEL_NAME)
    print("Instance: " + INSTANCE_NAME)
    print("Input INP file: " + INP_FILE)
    print("Output INP file: " + OUT_FILE)
    print("Peak field value: {}".format(PEAK_FIELD_VALUE))
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

        # Compute scaled peak value if cord distances are provided
        if site['upper_cord_sag_dist'] is not None and site['indent_cord_sag_dist'] is not None:
            site_peak_value = compute_scaled_peak_value(
                upper_cord_sag_dist=site['upper_cord_sag_dist'],
                indent_cord_sag_dist=site['indent_cord_sag_dist'],
                reference_preload=PEAK_FIELD_VALUE,
                ref_upper_cord_sag_dist=REF_UPPER_CORD_SAG_DIST,
                ref_indent_cord_sag_dist=REF_INDENT_CORD_SAG_DIST
            )
            site_compression_fraction = (site['upper_cord_sag_dist'] - site['indent_cord_sag_dist']) / site['upper_cord_sag_dist']
            print("Upper cord sag dist:       {:.5f} mm".format(site['upper_cord_sag_dist']))
            print("Indent cord sag dist:      {:.5f} mm".format(site['indent_cord_sag_dist']))
            print("Site compression fraction: {:.4f}".format(site_compression_fraction))
            print("Scaled peak field value:   {:.4f}".format(site_peak_value))
        else:
            site_peak_value = PEAK_FIELD_VALUE
            print("No cord distances in CSV. Using fallback PEAK_FIELD_VALUE: {}".format(site_peak_value))

        set_prefix = "{}_BAND".format(site['name'].upper().replace(' ', '_'))

        # First site reads from the original; subsequent sites read
        # from the output file which already contains previous sites
        read_from = INP_FILE if idx == 0 else OUT_FILE

        create_sinusoidal_node_sets(
            center_point=site['center'],
            upper_point=site['upper'],
            lower_point=site['lower'],
            inp_file=read_from,
            out_file=OUT_FILE,
            peak_field_value=site_peak_value,
            model_name=MODEL_NAME,
            instance_name=INSTANCE_NAME,
            set_prefix=set_prefix,
            site_index=site_num
        )

    print("\nAll {} sites processed.".format(len(sites)))
    print("Output: {}".format(OUT_FILE))

else:
    # --- Single site from variables ---
    if 'PEAK_FIELD_VALUE' not in dir():
        PEAK_FIELD_VALUE = 0.15
    OUT_FILE = generate_output_filename(INP_FILE)
    print("Model: " + MODEL_NAME)
    print("Instance: " + INSTANCE_NAME)
    print("Input INP file: " + INP_FILE)
    print("Output INP file: " + OUT_FILE)
    print("Peak field value: {}".format(PEAK_FIELD_VALUE))
    print("Center: {}".format(CENTER_POINT))
    print("Upper: {}".format(UPPER_POINT))
    print("Lower: {}".format(LOWER_POINT))

    create_sinusoidal_node_sets(
        center_point=CENTER_POINT,
        upper_point=UPPER_POINT,
        lower_point=LOWER_POINT,
        inp_file=INP_FILE,
        out_file=OUT_FILE,
        peak_field_value=PEAK_FIELD_VALUE,
        model_name=MODEL_NAME,
        instance_name=INSTANCE_NAME,
        set_prefix='FIELD_BAND',
        site_index=1
    )

    print("Done. Output: {}".format(OUT_FILE))
