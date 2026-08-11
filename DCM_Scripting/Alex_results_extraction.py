"""
Alex_results_extraction.py - Extract MPS (max principal strain) and element
volume (EVOL) from the 'Cord' element set, for every frame of the step.
Outputs one CSV per ODB to the same directory as the ODB.

Run in Abaqus CAE kernel:

ODB_PATH = r'D:\\path\\to\\Job-xxx.odb'
execfile('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting\\Alex_results_extraction.py')
"""

from odbAccess import *
from abaqusConstants import *
import os

# ============================================================
# USER SETTINGS (can be overridden by setting variables before execfile())
# ============================================================
if 'ODB_PATH' not in dir():
    raise RuntimeError("Set ODB_PATH before running this script.")
if 'CORD_SET_NAME' not in dir():
    CORD_SET_NAME = 'CORD'
if 'STEP_NAME' not in dir():
    STEP_NAME = None   # None = last step
# ============================================================

OUTPUT_DIR = os.path.dirname(ODB_PATH)

odb = openOdb(path=ODB_PATH, readOnly=True)

step = odb.steps[STEP_NAME] if STEP_NAME else odb.steps.values()[-1]
print("Step: '{}' ({} frames)".format(step.name, len(step.frames)))

# 'Cord' is a NODE set, not an element set - LE/EVOL are element-based field
# outputs, so getSubset(region=...) needs an element set. Derive one: every
# element whose FULL connectivity (all nodes) lies within the Cord node set.
# Computed once (not per frame) since connectivity doesn't change over time.
assembly = odb.rootAssembly
if CORD_SET_NAME in assembly.nodeSets.keys():
    # Assembly-level set: nodes is a tuple of MeshNodeArray, one per instance spanned
    node_set = assembly.nodeSets[CORD_SET_NAME]
    if len(node_set.instances) != 1:
        raise RuntimeError("Expected Cord node set to span exactly one instance; found {}.".format(
            len(node_set.instances)))
    instance = node_set.instances[0]
    cord_node_labels = set(n.label for n in node_set.nodes[0])
else:
    # Instance-level set: nodes is a flat MeshNodeArray
    if '.' in CORD_SET_NAME:
        inst_name, sname = CORD_SET_NAME.split('.', 1)
        instance = assembly.instances[inst_name]
    else:
        matches = [inst for inst in assembly.instances.values() if CORD_SET_NAME in inst.nodeSets.keys()]
        if len(matches) == 1:
            instance = matches[0]
        elif len(matches) > 1:
            raise RuntimeError("CORD_SET_NAME '{}' found in multiple instances; "
                                "specify as 'Instance.SetName'.".format(CORD_SET_NAME))
        else:
            raise RuntimeError("CORD_SET_NAME '{}' not found at assembly or "
                                "instance level (node set).".format(CORD_SET_NAME))
        sname = CORD_SET_NAME
    node_set = instance.nodeSets[sname]
    cord_node_labels = set(n.label for n in node_set.nodes)

cord_element_labels = set()
for elem in instance.elements:
    if all(nl in cord_node_labels for nl in elem.connectivity):
        cord_element_labels.add(elem.label)

print("Cord node set '{}': {} nodes -> {} fully-contained elements in instance '{}'".format(
    CORD_SET_NAME, len(cord_node_labels), len(cord_element_labels), instance.name))

rows = []

for frame_idx in range(len(step.frames)):
    frame = step.frames[frame_idx]

    le_subset   = frame.fieldOutputs['LE'].getSubset(region=instance, position=INTEGRATION_POINT)
    evol_subset = frame.fieldOutputs['EVOL'].getSubset(region=instance)

    # Take max MPS across integration points per element.
    # Follows precedent for element-wise maximum principal strain used in
    #  threshold-based injury volume metrics e.g. Cumulative Strain Damage Measure (CSDM).
    # CSDM: quantifies the fraction of tissue volume exceeding a given MPS threshold.
    # If any integration point within an element exceeds the threshold, the full
    #  element volume is counted as at risk.
    # For C3D4 elements (1 integration point) max = average = only value.
    #  Ref: Kleiven, S. (2007). Predictors for traumatic brain injuries
    #   evaluated through accident reconstructions. Ann. Adv. Automot. Med., 51,
    #   81-92.
    elem_mps = {}
    for val in le_subset.values:
        lbl = val.elementLabel
        if lbl not in cord_element_labels:
            continue
        if lbl not in elem_mps or val.maxPrincipal > elem_mps[lbl]:
            elem_mps[lbl] = val.maxPrincipal

    elem_vol = {val.elementLabel: val.data for val in evol_subset.values
                if val.elementLabel in cord_element_labels}

    for lbl, mps in elem_mps.items():
        if lbl in elem_vol:
            rows.append((frame_idx, frame.frameValue, lbl, mps, elem_vol[lbl]))

# Write CSV
odb_basename = os.path.splitext(os.path.basename(ODB_PATH))[0]
out_path = os.path.join(OUTPUT_DIR, '{}_mps.csv'.format(odb_basename))

with open(out_path, 'w') as f:
    f.write('frame_index,frame_value,element_label,mps,volume\n')
    for r in rows:
        f.write('{},{:.6e},{},{:.6e},{:.6e}\n'.format(*r))

print("Saved {} rows -> {}".format(len(rows), out_path))
print("Add this path to id_map.csv 'csv_path' column for this participant.")

# --- Element connectivity (topology), for spatial clustering analysis
# (blob vs scattered strain distribution) downstream. One-time per
# patient/condition since connectivity doesn't change over time. This is a
# separate file with a different suffix ('_topology.csv' vs '_mps.csv') -
# purely additive, does not touch the CSV written above. No id_map.csv
# column needed for it either: downstream code derives this path from
# csv_path by swapping the suffix (same directory, same basename).
topology_path = os.path.join(OUTPUT_DIR, '{}_topology.csv'.format(odb_basename))
with open(topology_path, 'w') as f:
    f.write('element_label,nodes\n')
    for elem in instance.elements:
        if elem.label in cord_element_labels:
            nodes = ';'.join(str(n) for n in elem.connectivity)
            f.write('{},{}\n'.format(elem.label, nodes))
print("Saved topology for {} elements -> {}".format(len(cord_element_labels), topology_path))

# --- Auto-populate id_map.csv's 'csv_path' for the matching row (matched by
# ODB_PATH), removing the manual copy-paste step. Uses the stdlib csv module
# rather than pandas, since Abaqus's bundled Python isn't guaranteed to have
# pandas installed. Only the matched row's csv_path is touched - every other
# column/row is preserved exactly as read.
if 'ID_MAP_PATH' not in dir():
    ID_MAP_PATH = r'C:\Users\cmb247\repos\Abaqus\DCM_Scripting\id_map.csv'

import csv

with open(ID_MAP_PATH, 'rb') as f:
    id_map_reader = csv.DictReader(f)
    id_map_fieldnames = id_map_reader.fieldnames
    id_map_rows = list(id_map_reader)

id_map_matches = 0
for id_map_row in id_map_rows:
    if id_map_row.get('odb_path', '').strip() == ODB_PATH.strip():
        id_map_row['csv_path'] = out_path
        id_map_matches += 1

if id_map_matches:
    with open(ID_MAP_PATH, 'wb') as f:
        id_map_writer = csv.DictWriter(f, fieldnames=id_map_fieldnames)
        id_map_writer.writeheader()
        id_map_writer.writerows(id_map_rows)
    print("Updated csv_path in id_map.csv for {} matching row(s).".format(id_map_matches))
else:
    print("WARNING: ODB_PATH not found in id_map.csv 'odb_path' column - csv_path not auto-filled. "
          "Add this row to id_map.csv first, or update it manually.")

odb.close()
