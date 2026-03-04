"""
Spinal Cord Node Set Creation Script - With Overlap Detection and Resolution
=============================================================================
Extends sets_scaled_inpmod_script.py with automatic detection and resolution
of overlapping node assignments when multiple compression sites are processed
from a CSV file.

Overlap resolution strategy (multi-site CSV path only):
    When a node falls within the band regions of more than one site, it is
    assigned exclusively to the site that gives it the MINIMUM field value
    (based on each site's raised-cosine band profile). This is physically
    conservative: transition zones between adjacent compression levels receive
    the lower preload rather than being double-loaded.

    A significant-overlap warning is printed for any site pair where contested
    nodes fall in a high-field band (band index 0 or 1) of either site. This
    indicates that the upper/lower extents in the CSV may be configured too
    broadly.

    Same-level anterior+posterior compressions should be defined as a SINGLE
    site in the CSV - no special overlap handling is needed for that case.

Also includes cord-diameter-based preload scaling:
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
    execfile('sets_scaled_inpmod_overlap.py')

Option 1 - Single site (manual PEAK_FIELD_VALUE, no scaling):
    MODEL_NAME = 'Model-1'
    INSTANCE_NAME = 'PART-1_1-1'
    INP_FILE = 'Job-212.inp'
    PEAK_FIELD_VALUE = 0.5
    CENTER_POINT = (x, y, z)
    UPPER_POINT = (x, y, z)
    LOWER_POINT = (x, y, z)
    execfile('sets_scaled_inpmod_overlap.py')

Option 2 - Multiple sites from CSV (overlap detection + scaling applied automatically):
    MODEL_NAME = 'Model-1'
    INSTANCE_NAME = 'PART-1_1-1'
    INP_FILE = 'Job-212.inp'
    COORDS_FILE = 'coordinates_scaled.csv'
    execfile('sets_scaled_inpmod_overlap.py')

    Optionally override reference calibration before execfile:
    REFERENCE_PRELOAD        = 0.5
    REF_UPPER_CORD_SAG_DIST  = 6.82259
    REF_INDENT_CORD_SAG_DIST = 4.24591

CSV format (coordinates_scaled.csv):
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
import re


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


def generate_overlap_filenames(inp_file, reference_preload, site_names):
    """
    Generate output .inp and summary .txt filenames for the multi-site overlap path.

    Transformations applied to the input base name:
      1. 'template' (case-insensitive) is removed along with adjacent separators (- or _).
         Any resulting consecutive separators are collapsed to a single underscore.
      2. Preload level is appended, formatted as e.g. '0pt50' (decimal point -> 'pt').
      3. Site names (from CSV, in order) are appended, separated by '_'.

    The summary .txt receives the same base name as the output .inp.

    Example:
      inp_file          = 'D:/path/Job-N01-015-TEMPLATE_2STEP.inp'
      reference_preload = 0.5
      site_names        = ['Site1', 'Site2']
      -> out_inp     = 'D:/path/Job-N01-015_2STEP_0pt50_Site1_Site2.inp'
      -> out_summary = 'D:/path/Job-N01-015_2STEP_0pt50_Site1_Site2_overlap_summary.txt'

    Parameters:
        inp_file:          path to the original .inp file
        reference_preload: float, preload level to embed in filename
        site_names:        list of site name strings in CSV order

    Returns:
        tuple: (out_inp_file, out_summary_file)
    """
    base, ext = os.path.splitext(inp_file)

    # Remove 'template' along with any immediately adjacent separators (-, _)
    base = re.sub(r'[-_]*template[-_]*', '_', base, flags=re.IGNORECASE)

    # Collapse any resulting consecutive separators and strip trailing ones
    base = re.sub(r'[-_]{2,}', '_', base)
    base = base.rstrip('_-')

    # Format preload: 0.5 -> '0pt50', 0.15 -> '0pt15'
    preload_str = '{:.2f}'.format(reference_preload).replace('.', 'pt')

    # Assemble new base: original_cleaned _ preload _ site1 _ site2 ...
    new_base = '{}_{}_{}'.format(base, preload_str, '_'.join(site_names))

    out_inp     = '{}{}'.format(new_base, ext)
    out_summary = '{}_overlap_summary.txt'.format(new_base)

    return out_inp, out_summary


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


def classify_nodes_for_site(
    center_point,
    upper_point,
    lower_point,
    assembly,
    instance_name,
    num_bands=5
):
    """
    Classify all nodes in an instance into bands for one compression site.

    This is the pure classification kernel, decoupled from .inp file writing.
    Called by the multi-site overlap-resolution path (Pass 1). The single-site
    path continues to use create_sinusoidal_node_sets which contains equivalent
    logic inline.

    Parameters:
        center_point:  (x, y, z) center of the compression region
        upper_point:   (x, y, z) upper limit of the region
        lower_point:   (x, y, z) lower limit of the region
        assembly:      Abaqus rootAssembly object (resolved once by caller)
        instance_name: name of the part instance in the assembly
        num_bands:     number of bands (default 5)

    Returns:
        list of num_bands lists of integer node labels, or None if instance
        is not found.
    """
    axis_vector = (
        upper_point[0] - lower_point[0],
        upper_point[1] - lower_point[1],
        upper_point[2] - lower_point[2]
    )
    axis_vector = normalize_vector(axis_vector)

    d_upper = calculate_axis_projection(upper_point, axis_vector, center_point)
    d_lower = calculate_axis_projection(lower_point, axis_vector, center_point)

    print("  Axis vector: ({:.3f}, {:.3f}, {:.3f})".format(*axis_vector))
    print("  Distance to upper limit: {:.3f}".format(d_upper))
    print("  Distance to lower limit: {:.3f}".format(d_lower))

    try:
        instance = assembly.instances[instance_name]
        all_nodes = instance.nodes
        print("  Found {} nodes in instance '{}'".format(len(all_nodes), instance_name))
    except KeyError:
        print("Error: Instance '{}' not found in assembly.".format(instance_name))
        print("Available instances: {}".format(list(assembly.instances.keys())))
        return None

    band_nodes = [[] for _ in range(num_bands)]
    nodes_outside = 0

    for node in all_nodes:
        coords = node.coordinates
        distance = calculate_axis_projection(coords, axis_vector, center_point)
        band_idx = get_band_index(distance, d_upper, d_lower, num_bands)
        if band_idx >= 0:
            band_nodes[band_idx].append(node.label)
        else:
            nodes_outside += 1

    print("  Nodes classified: {}  outside: {}".format(
        sum([len(b) for b in band_nodes]), nodes_outside))
    return band_nodes


def resolve_overlaps(all_site_band_nodes, all_site_field_values, site_names, num_bands=5):
    """
    Resolve multi-site node conflicts using minimum-field-value assignment.

    For each node appearing in more than one site's band lists, keeps it only
    in the site+band that gives it the lowest field value. Removes it from all
    other competing sites' bands. Mutates all_site_band_nodes in-place.

    Prints a warning for each site pair where contested nodes fall in a
    high-field band (band index 0 or 1) of either site -- indicates that
    the upper/lower extents in the CSV may be configured too broadly.

    Parameters:
        all_site_band_nodes   -- list[site_idx][band_idx] -> list of int node labels
                                 (mutated in-place)
        all_site_field_values -- list[site_idx][band_idx] -> float field value
        site_names            -- list of site name strings (for warning messages)
        num_bands             -- number of bands (default 5)

    Returns:
        overlap_stats -- list of dicts, one per overlapping site pair:
            {'site_i', 'site_j', 'name_i', 'name_j',
             'total_contested', 'high_field_contested',
             'kept_by_i', 'kept_by_j'}
    """
    num_sites = len(all_site_band_nodes)

    # Step 1: Build reverse index -- node_label -> list of (site_idx, band_idx)
    node_to_assignments = {}
    for site_idx in range(num_sites):
        for band_idx in range(num_bands):
            for label in all_site_band_nodes[site_idx][band_idx]:
                if label not in node_to_assignments:
                    node_to_assignments[label] = []
                node_to_assignments[label].append((site_idx, band_idx))

    # Step 2: Collect contested nodes (claimed by more than one site)
    contested = {}
    for label, assignments in node_to_assignments.items():
        if len(assignments) > 1:
            contested[label] = assignments

    print("\nOverlap resolution: {} contested node(s) found across {} site(s).".format(
        len(contested), num_sites))

    if not contested:
        print("  No overlapping nodes detected.")
        return []

    # Step 3: Per-pair statistics accumulators keyed by canonical (i, j) pair
    pair_stats = {}

    # Step 4: Resolve each contested node
    for label, assignments in contested.items():
        # Find the assignment with the minimum field value
        best = assignments[0]
        best_val = all_site_field_values[best[0]][best[1]]
        for assignment in assignments[1:]:
            val = all_site_field_values[assignment[0]][assignment[1]]
            if val < best_val:
                best_val = val
                best = assignment

        winning_site = best[0]

        # Remove from every losing site and accumulate per-pair stats
        for assignment in assignments:
            if assignment == best:
                continue
            losing_site, losing_band = assignment
            pair = (min(winning_site, losing_site), max(winning_site, losing_site))
            if pair not in pair_stats:
                pair_stats[pair] = {
                    'total_contested': 0,
                    'high_field_contested': 0,
                    'kept_by_i': 0,
                    'kept_by_j': 0
                }
            pair_stats[pair]['total_contested'] += 1
            # High-field warning: node is in band 0 or 1 of any competing site
            for a in assignments:
                if a[1] <= 1:
                    pair_stats[pair]['high_field_contested'] += 1
                    break
            if winning_site == pair[0]:
                pair_stats[pair]['kept_by_i'] += 1
            else:
                pair_stats[pair]['kept_by_j'] += 1
            # Remove node from the losing site's band (label is unique per band)
            all_site_band_nodes[losing_site][losing_band].remove(label)

    # Step 5: Print per-pair summary table
    overlap_stats = []
    print("\n" + "-" * 70)
    print("Overlap summary (minimum-field-value assignment):")
    print("{:<14} {:<14} {:>12} {:>12} {:>8} {:>8}".format(
        "Site A", "Site B", "Contested", "High-field", "->A", "->B"))
    print("-" * 70)
    for pair in sorted(pair_stats.keys()):
        i, j = pair
        s = pair_stats[pair]
        hf = s['high_field_contested']
        warn_tag = "  *** WARNING ***" if hf > 0 else ""
        print("{:<14} {:<14} {:>12} {:>12} {:>8} {:>8}{}".format(
            site_names[i], site_names[j],
            s['total_contested'], hf,
            s['kept_by_i'], s['kept_by_j'],
            warn_tag))
        if hf > 0:
            print("  WARNING: {} high-field (inner band) node(s) contested between"
                  " '{}' and '{}'.".format(hf, site_names[i], site_names[j]))
            print("  Check that upper/lower extents in the CSV are not too broad.")
        overlap_stats.append({
            'site_i': i, 'site_j': j,
            'name_i': site_names[i], 'name_j': site_names[j],
            'total_contested': s['total_contested'],
            'high_field_contested': hf,
            'kept_by_i': s['kept_by_i'], 'kept_by_j': s['kept_by_j']
        })
    print("-" * 70)

    return overlap_stats


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

print("sets_scaled_inpmod_overlap.py loaded")

# --- Reference calibration constants ---
# Override any of these before execfile() if your reference site differs.
if 'REFERENCE_PRELOAD' not in dir():
    REFERENCE_PRELOAD = 0.5
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
    print("    execfile('sets_scaled_inpmod_overlap.py')")
    print("")
    print("  Option 2 - Multiple sites from CSV (overlap detection + auto-scaled):")
    print("    MODEL_NAME = 'Model-1'")
    print("    INSTANCE_NAME = 'PART-1_1-1'")
    print("    INP_FILE = 'Job-212.inp'")
    print("    COORDS_FILE = 'coordinates_scaled.csv'")
    print("    execfile('sets_scaled_inpmod_overlap.py')")
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
    # --- Multiple sites from CSV file (three-pass: classify, resolve, write) ---
    if 'PEAK_FIELD_VALUE' not in dir():
        PEAK_FIELD_VALUE = REFERENCE_PRELOAD
    print("Model: " + MODEL_NAME)
    print("Instance: " + INSTANCE_NAME)
    print("Input INP file: " + INP_FILE)
    print("Peak field value: {}".format(PEAK_FIELD_VALUE))
    print("Coordinates file: " + COORDS_FILE)

    sites = read_coordinates_file(COORDS_FILE)
    print("Found {} sites in CSV file".format(len(sites)))

    OUT_FILE, SUMMARY_FILE = generate_overlap_filenames(
        INP_FILE, PEAK_FIELD_VALUE, [s['name'] for s in sites])
    print("Output INP file: " + OUT_FILE)
    print("Summary file:    " + SUMMARY_FILE)

    # Resolve model and assembly once for all sites
    try:
        modelDB = mdb.models[MODEL_NAME]
    except KeyError:
        print("Error: Model '{}' not found.".format(MODEL_NAME))
        print("Available models: {}".format(list(mdb.models.keys())))
        modelDB = None

    if modelDB is not None:
        assembly = modelDB.rootAssembly

        # Accumulators populated in Pass 1
        all_site_band_nodes   = []  # [site_idx][band_idx] -> list of int node labels
        all_site_field_values = []  # [site_idx][band_idx] -> float
        all_site_peak_values  = []  # [site_idx] -> float
        all_site_prefixes     = []  # [site_idx] -> set_prefix string
        site_names            = []  # [site_idx] -> site name string

        # ----------------------------------------------------------------
        # PASS 1: Classify all sites -- no writing yet
        # ----------------------------------------------------------------
        print("\n" + "=" * 60)
        print("PASS 1: Classifying nodes for all {} sites".format(len(sites)))
        print("=" * 60)

        classification_ok = True
        for idx, site in enumerate(sites):
            site_num = idx + 1
            print("\n" + "#" * 60)
            print("SITE {}: {}".format(site_num, site['name']))
            print("#" * 60)
            print("Center: {}".format(site['center']))
            print("Upper:  {}".format(site['upper']))
            print("Lower:  {}".format(site['lower']))

            # Compute scaled peak value if cord distances are provided
            if (site['upper_cord_sag_dist'] is not None and
                    site['indent_cord_sag_dist'] is not None):
                site_peak_value = compute_scaled_peak_value(
                    upper_cord_sag_dist=site['upper_cord_sag_dist'],
                    indent_cord_sag_dist=site['indent_cord_sag_dist'],
                    reference_preload=PEAK_FIELD_VALUE,
                    ref_upper_cord_sag_dist=REF_UPPER_CORD_SAG_DIST,
                    ref_indent_cord_sag_dist=REF_INDENT_CORD_SAG_DIST
                )
                site_compression_fraction = (
                    (site['upper_cord_sag_dist'] - site['indent_cord_sag_dist'])
                    / site['upper_cord_sag_dist']
                )
                print("Upper cord sag dist:       {:.5f} mm".format(site['upper_cord_sag_dist']))
                print("Indent cord sag dist:      {:.5f} mm".format(site['indent_cord_sag_dist']))
                print("Site compression fraction: {:.4f}".format(site_compression_fraction))
                print("Scaled peak field value:   {:.4f}".format(site_peak_value))
            else:
                site_peak_value = PEAK_FIELD_VALUE
                print("No cord distances in CSV. Using fallback PEAK_FIELD_VALUE: {}".format(
                    site_peak_value))

            set_prefix = "{}_BAND".format(site['name'].upper().replace(' ', '_'))

            band_nodes = classify_nodes_for_site(
                center_point=site['center'],
                upper_point=site['upper'],
                lower_point=site['lower'],
                assembly=assembly,
                instance_name=INSTANCE_NAME
            )

            if band_nodes is None:
                print("ERROR: Node classification failed for site {}. Aborting.".format(
                    site_num))
                classification_ok = False
                break

            field_values = calculate_field_values(5, site_peak_value, 0.0)

            all_site_band_nodes.append(band_nodes)
            all_site_field_values.append(field_values)
            all_site_peak_values.append(site_peak_value)
            all_site_prefixes.append(set_prefix)
            site_names.append(site['name'])

        if classification_ok:
            # ----------------------------------------------------------------
            # PASS 2: Resolve overlaps using minimum-field-value assignment
            # ----------------------------------------------------------------
            print("\n" + "=" * 60)
            print("PASS 2: Resolving overlaps across {} sites".format(len(sites)))
            print("=" * 60)

            overlap_stats = resolve_overlaps(
                all_site_band_nodes=all_site_band_nodes,
                all_site_field_values=all_site_field_values,
                site_names=site_names
            )

            # ----------------------------------------------------------------
            # PASS 3: Write all resolved sites sequentially to .inp file
            # ----------------------------------------------------------------
            print("\n" + "=" * 60)
            print("PASS 3: Writing all {} sites to .inp file".format(len(sites)))
            print("=" * 60)

            for idx in range(len(sites)):
                site_num = idx + 1
                set_prefix = all_site_prefixes[idx]
                band_nodes = all_site_band_nodes[idx]
                field_values = all_site_field_values[idx]

                print("\n" + "-" * 60)
                print("Writing site {}: {}".format(site_num, site_names[idx]))

                # First site reads from original; subsequent sites read from
                # OUT_FILE which already contains the previous sites' sets.
                read_from = INP_FILE if idx == 0 else OUT_FILE

                write_sets_to_inp(
                    inp_file=read_from,
                    out_file=OUT_FILE,
                    instance_name=INSTANCE_NAME,
                    band_nodes=band_nodes,
                    field_values=field_values,
                    set_prefix=set_prefix,
                    site_index=site_num,
                    amplitude_name='Amp-1-preload',
                    num_bands=5
                )

                # Per-site summary table
                print("\n" + "=" * 75)
                print("SUMMARY: Site {} - {}".format(site_num, site_names[idx]))
                print("=" * 75)
                print("{:<20} {:>10} {:>15} {:<30}".format(
                    "Node Set", "Nodes", "Field Value", "Predefined Field"))
                print("-" * 75)
                for band_idx in range(5):
                    set_name = "{}_{:d}".format(set_prefix, band_idx + 1)
                    field_name = "predefinedfield-{:d}-fieldband{:d}".format(
                        site_num, band_idx + 1)
                    node_count = len(band_nodes[band_idx])
                    print("{:<20} {:>10} {:>15.3f} {:<30}".format(
                        set_name, node_count, field_values[band_idx], field_name))
                print("=" * 75)

        print("\nAll {} sites processed.".format(len(sites)))
        print("Output: {}".format(OUT_FILE))

        # ----------------------------------------------------------------
        # Write summary file
        # ----------------------------------------------------------------
        with open(SUMMARY_FILE, 'w') as sf:
            sf.write("sets_scaled_inpmod_overlap.py - Run Summary\n")
            sf.write("=" * 75 + "\n\n")
            sf.write("Model:            {}\n".format(MODEL_NAME))
            sf.write("Instance:         {}\n".format(INSTANCE_NAME))
            sf.write("Input INP file:   {}\n".format(INP_FILE))
            sf.write("Output INP file:  {}\n".format(OUT_FILE))
            sf.write("Coordinates file: {}\n".format(COORDS_FILE))
            sf.write("Sites processed:  {}\n\n".format(len(sites)))

            # Overlap resolution section
            sf.write("Overlap Resolution\n")
            sf.write("-" * 75 + "\n")
            sf.write("Column guide:\n")
            sf.write("  Contested  : total nodes found in both sites' band regions before resolution\n")
            sf.write("  High-field : contested nodes that were in an inner band (band 1 or 2,\n")
            sf.write("               highest field values) of either site -- see warnings below\n")
            sf.write("  ->A / ->B  : how many contested nodes were ultimately assigned to each\n")
            sf.write("               site after applying the minimum-field-value rule\n")
            sf.write("\n")
            sf.write("Resolution rule: each contested node is assigned exclusively to whichever\n")
            sf.write("site gives it the lower field value (based on that site's raised-cosine\n")
            sf.write("band profile). It is removed from the other site's band. This ensures\n")
            sf.write("transition zones between adjacent compression levels receive the lower\n")
            sf.write("preload rather than being double-loaded.\n")
            sf.write("\n")
            if classification_ok and overlap_stats:
                sf.write("{:<14} {:<14} {:>12} {:>12} {:>8} {:>8}\n".format(
                    "Site A", "Site B", "Contested", "High-field", "->A", "->B"))
                sf.write("-" * 70 + "\n")
                for s in overlap_stats:
                    hf = s['high_field_contested']
                    warn_tag = "  *** WARNING ***" if hf > 0 else ""
                    sf.write("{:<14} {:<14} {:>12} {:>12} {:>8} {:>8}{}\n".format(
                        s['name_i'], s['name_j'],
                        s['total_contested'], hf,
                        s['kept_by_i'], s['kept_by_j'],
                        warn_tag))
                sf.write("-" * 70 + "\n")
                sf.write("\n")
                sf.write("Notes:\n")
                for s in overlap_stats:
                    hf = s['high_field_contested']
                    total = s['total_contested']
                    ki = s['kept_by_i']
                    kj = s['kept_by_j']
                    name_i = s['name_i']
                    name_j = s['name_j']
                    sf.write("  {name_i} vs {name_j}: {total} node(s) were found within the band\n"
                             "    regions of both sites. After resolution, {ki} node(s) remained\n"
                             "    assigned to {name_i} and {kj} node(s) to {name_j} (removed from\n"
                             "    {name_i}'s bands).\n".format(
                                 name_i=name_i, name_j=name_j,
                                 total=total, ki=ki, kj=kj))
                    if hf == 0:
                        sf.write("    All contested nodes were in outer bands (low field values) --\n"
                                 "    overlap is minor and within expected boundary behaviour.\n")
                    else:
                        sf.write("    *** WARNING: {hf} of the contested node(s) were in an inner\n"
                                 "    band (high field value) of at least one site. This suggests\n"
                                 "    the upper/lower extents defined in the CSV for {name_i} and\n"
                                 "    {name_j} may overlap too broadly. Review the coordinates and\n"
                                 "    consider reducing the extents to avoid significant overlap\n"
                                 "    between these sites. ***\n".format(
                                     hf=hf, name_i=name_i, name_j=name_j))
            else:
                sf.write("  No overlapping nodes detected between any pair of sites.\n")
            sf.write("\n")

            # Per-site summary tables
            sf.write("Site Summaries\n")
            sf.write("-" * 75 + "\n")
            for idx in range(len(sites)):
                site_num = idx + 1
                set_prefix = all_site_prefixes[idx]
                band_nodes_s = all_site_band_nodes[idx]
                field_values_s = all_site_field_values[idx]
                site = sites[idx]

                sf.write("\n")
                sf.write("=" * 75 + "\n")
                sf.write("SUMMARY: Site {} - {}\n".format(site_num, site_names[idx]))
                sf.write("=" * 75 + "\n")
                sf.write("Center: {}\n".format(site['center']))
                sf.write("Upper:  {}\n".format(site['upper']))
                sf.write("Lower:  {}\n".format(site['lower']))
                if (site['upper_cord_sag_dist'] is not None and
                        site['indent_cord_sag_dist'] is not None):
                    sf.write("Upper cord sag dist:  {:.5f} mm\n".format(
                        site['upper_cord_sag_dist']))
                    sf.write("Indent cord sag dist: {:.5f} mm\n".format(
                        site['indent_cord_sag_dist']))
                sf.write("Peak field value:     {:.4f}\n\n".format(
                    all_site_peak_values[idx]))
                sf.write("{:<20} {:>10} {:>15} {:<30}\n".format(
                    "Node Set", "Nodes", "Field Value", "Predefined Field"))
                sf.write("-" * 75 + "\n")
                for band_idx in range(5):
                    set_name = "{}_{:d}".format(set_prefix, band_idx + 1)
                    field_name = "predefinedfield-{:d}-fieldband{:d}".format(
                        site_num, band_idx + 1)
                    node_count = len(band_nodes_s[band_idx])
                    sf.write("{:<20} {:>10} {:>15.3f} {:<30}\n".format(
                        set_name, node_count, field_values_s[band_idx], field_name))
                sf.write("=" * 75 + "\n")

        print("Summary written to: {}".format(SUMMARY_FILE))

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
