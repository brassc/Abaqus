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
    ODB_PATH = r'D:\Charlotte\ABAQUS\N01-011\Pre-Op\Job-020-N01-011-PreOp-BC0pt35wEVOL\Job-020-N01-011-PreOp-BC0pt35wEVOL_0pt30_Site1_Site2_Site3_Site4.odb'
if 'CORD_SET_NAME' not in dir():
    CORD_SET_NAME = 'Cord'
if 'STEP_NAME' not in dir():
    STEP_NAME = None   # None = last step
# ============================================================

OUTPUT_DIR = os.path.dirname(ODB_PATH)

odb = openOdb(path=ODB_PATH, readOnly=True)

step = odb.steps[STEP_NAME] if STEP_NAME else odb.steps.values()[-1]
print("Step: '{}' ({} frames)".format(step.name, len(step.frames)))

# Resolve element set - try assembly level, then instance level
assembly = odb.rootAssembly
if CORD_SET_NAME in assembly.elementSets.keys():
    elem_set = assembly.elementSets[CORD_SET_NAME]
else:
    inst_name, sname = CORD_SET_NAME.split('.', 1)
    elem_set = assembly.instances[inst_name].elementSets[sname]

print("Processing set: '{}'".format(CORD_SET_NAME))

rows = []

for frame_idx in range(len(step.frames)):
    frame = step.frames[frame_idx]

    le_subset   = frame.fieldOutputs['LE'].getSubset(region=elem_set, position=INTEGRATION_POINT)
    evol_subset = frame.fieldOutputs['EVOL'].getSubset(region=elem_set)

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
        if lbl not in elem_mps or val.maxPrincipal > elem_mps[lbl]:
            elem_mps[lbl] = val.maxPrincipal

    elem_vol = {val.elementLabel: val.data for val in evol_subset.values}

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

odb.close()
