"""
Alex_results_extraction_GM_WM.py - Same per-frame MPS/EVOL extraction as
Alex_results_extraction.py (max principal strain + element volume, every
frame of the step, over the 'Cord' node set), but additionally tags every row
with which tissue the element belongs to: 'GM' (grey matter) or 'WM' (white
matter).

Unlike the Cord region (a NODE set, reconstructed into elements by full-
containment), GM and WM are resolved as ELEMENT sets directly - same approach
already used by the older single-patient pipeline, Data_Extraction/
mps_extract.py (SET_NAMES = ['PART-1-1.P44;GM', 'PART-1-1.P45;WM']). Element
sets give an unambiguous per-element material tag; node-set containment would
be ambiguous at the GM/WM interface, where boundary elements can share nodes
with both tissues.

The exact set names (e.g. 'P44;GM') vary per patient's part/instance
numbering, so this script auto-detects whichever element set's name ends in
';GM' / ';WM' (case-insensitive), instance-level first then assembly-level,
rather than hardcoding a name.

Output: '<odb_basename>_mps_GM_WM.csv' (frame_index, frame_value,
element_label, mps, volume, tissue_type), written alongside the existing
'<odb_basename>_mps.csv' - a separate file, not a replacement, so it doesn't
disturb any downstream script that already depends on the plain _mps.csv's
column schema or id_map.csv's 'csv_path' pointing to it. Downstream GM/WM
scripts should derive this path the same way '_topology.csv' is already
derived elsewhere in this codebase: csv_path.replace('_mps.csv',
'_mps_GM_WM.csv') - no id_map.csv schema change needed.

Does NOT write '_topology.csv' - that's unchanged and already produced by a
plain run of Alex_results_extraction.py; run that first for a given ODB, then
this script for the GM/WM tagging.

Does NOT auto-populate id_map.csv - unlike Alex_results_extraction.py's
'csv_path' column (which points at the plain, untagged file every other
script expects), there is no separate id_map.csv column for the GM/WM
variant; downstream code derives it from 'csv_path' as described above.

Fails gracefully the same way Alex_results_extraction_IVD.py does, for the
same reason (batching across many ODBs, some of which may have been written
by a different Abaqus version than the kernel running this script): an
openOdb() failure is caught, classified (version mismatch vs. other), and
recorded so a `for ODB_PATH in ODB_PATHS: execfile(...)` loop moves on to the
next ODB_PATH instead of dying. A missing/misnamed GM or WM element set (or
Cord set) on an ODB that DID open successfully is NOT caught - that is a real
modeling/naming problem on that specific job worth stopping on, not something
to silently skip past.

Must run under both Abaqus 2022 (Python 2.7 kernel) and Abaqus 2025 (Python
3.10 kernel) for the same reason as the IVD script - kept to syntax valid on
both (.format() not f-strings, print() function, list(...values()) instead of
indexing a Repository's .values() directly).

Run in Abaqus CAE kernel:

ODB_PATH = r'D:\\path\\to\\Job-xxx.odb'
execfile('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting\\Alex_results_extraction_GM_WM.py')

Or from the command line:
abaqus python Alex_results_extraction_GM_WM.py <odb_path>
"""

from odbAccess import *
from abaqusConstants import *
import os
import re
import sys

# ============================================================
# USER SETTINGS (can be overridden by setting variables before execfile())
# ============================================================
if 'ODB_PATH' not in dir():
    if len(sys.argv) > 1:
        ODB_PATH = sys.argv[1]
    else:
        raise RuntimeError("Set ODB_PATH before running this script.")
if 'CORD_SET_NAME' not in dir():
    CORD_SET_NAME = 'CORD'
if 'STEP_NAME' not in dir():
    STEP_NAME = None   # None = last step
# ============================================================

# This script's own absolute path - used only to print a ready-to-paste
# re-run block for ODBs that need a different Abaqus version (see the batch
# summary at the bottom). Not used anywhere in the extraction itself.
SELF_PATH = r'C:\Users\cmb247\repos\Abaqus\DCM_Scripting\Alex_results_extraction_GM_WM.py'

OUTPUT_DIR = os.path.dirname(ODB_PATH)

# Batch-loop bookkeeping - persists across repeated execfile() calls in the
# same kernel session (each loop iteration re-execfile's this whole script,
# so a plain assignment here would otherwise reset the list every time).
if 'VERSION_MISMATCH_ODB_PATHS' not in dir():
    VERSION_MISMATCH_ODB_PATHS = []
if 'OTHER_FAILED_ODB_PATHS' not in dir():
    OTHER_FAILED_ODB_PATHS = []

# Open the ODB defensively: an ODB written by a different Abaqus version than
# the one running this kernel raises here (wrong release for this ODB), and
# without a try/except that exception would propagate out of execfile() and
# kill the REST of the batch loop, not just this one ODB_PATH. Catch it,
# record it, and let the script finish normally so the loop moves on to the
# next ODB_PATH.
odb = None
try:
    odb = openOdb(path=ODB_PATH, readOnly=True)
except Exception as e:
    msg = str(e)
    is_version_mismatch = any(kw in msg.lower() for kw in ('version', 'release', 'incompatib'))
    if is_version_mismatch:
        print("SKIPPING (wrong Abaqus version for this ODB): {}".format(ODB_PATH))
        VERSION_MISMATCH_ODB_PATHS.append(ODB_PATH)
    else:
        print("SKIPPING (failed to open ODB): {}".format(ODB_PATH))
        OTHER_FAILED_ODB_PATHS.append(ODB_PATH)
    print("  {}".format(msg))

if odb is not None:
    try:
        step = odb.steps[STEP_NAME] if STEP_NAME else list(odb.steps.values())[-1]
        print("Step: '{}' ({} frames)".format(step.name, len(step.frames)))

        # 'Cord' is a NODE set, not an element set - LE/EVOL are element-based
        # field outputs, so getSubset(region=...) needs an element set. Derive
        # one: every element whose FULL connectivity (all nodes) lies within
        # the Cord node set. Computed once (not per frame) since connectivity
        # doesn't change over time. Identical logic to
        # Alex_results_extraction.py - this is still how "is this element
        # part of the cord at all" gets decided; GM/WM tissue typing below is
        # a separate, additional classification on top of this.
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

        # --------------------------------------------------------------
        # GM/WM tissue typing - resolved as ELEMENT sets (not node sets),
        # matched by name suffix since the exact set name varies per
        # patient's part numbering.
        # --------------------------------------------------------------
        def _find_elementset_by_suffix(suffix):
            """Find the single element set (instance-level first, else assembly-level) whose
            name ends with ';<suffix>' (case-insensitive) - e.g. 'P44;GM'. Suffix match (not
            an exact name) because the leading part number varies per patient's mesh."""
            pattern = re.compile(r';\s*' + suffix + r'\s*$', re.IGNORECASE)
            inst_matches = [name for name in instance.elementSets.keys() if pattern.search(name)]
            if len(inst_matches) == 1:
                return instance.elementSets[inst_matches[0]], "instance '{}'".format(inst_matches[0])
            if len(inst_matches) > 1:
                raise RuntimeError("Expected exactly one instance-level element set ending in ';{}'; "
                                    "found {}: {}".format(suffix, len(inst_matches), inst_matches))

            asm_matches = [name for name in assembly.elementSets.keys() if pattern.search(name)]
            if len(asm_matches) == 1:
                return assembly.elementSets[asm_matches[0]], "assembly '{}'".format(asm_matches[0])
            raise RuntimeError(
                "Expected exactly one element set ending in ';{}' (instance or assembly level); "
                "found {} instance-level and {} assembly-level matches.".format(
                    suffix, len(inst_matches), len(asm_matches)))

        def _element_labels_of(elem_set):
            """Flatten an OdbSet's .elements into a plain set of element labels - handles both
            assembly-level (tuple of MeshElementArray, one per spanned instance) and
            instance-level (flat MeshElementArray) sets, same distinction as the Cord node
            set resolution above. An instance-level set's .instances attribute EXISTS but is
            None (not absent/empty) - hasattr() alone doesn't catch that, so check truthiness
            via getattr() instead."""
            instances = getattr(elem_set, 'instances', None)
            if instances:
                if len(instances) != 1:
                    raise RuntimeError("Expected element set to span exactly one instance; found {}.".format(
                        len(instances)))
                return set(e.label for e in elem_set.elements[0])
            return set(e.label for e in elem_set.elements)

        gm_elem_set, gm_source = _find_elementset_by_suffix('GM')
        wm_elem_set, wm_source = _find_elementset_by_suffix('WM')
        gm_element_labels = _element_labels_of(gm_elem_set)
        wm_element_labels = _element_labels_of(wm_elem_set)

        print("GM element set: {} ({} elements)".format(gm_source, len(gm_element_labels)))
        print("WM element set: {} ({} elements)".format(wm_source, len(wm_element_labels)))

        overlap = gm_element_labels & wm_element_labels
        if overlap:
            print("WARNING: {} element(s) are in BOTH the GM and WM element sets - "
                  "check the model. Tagged 'GM' (arbitrary tie-break).".format(len(overlap)))

        cord_untagged = cord_element_labels - gm_element_labels - wm_element_labels
        if cord_untagged:
            print("WARNING: {} Cord element(s) are in NEITHER the GM nor WM element set - "
                  "tagged 'UNKNOWN'. Check CORD_SET_NAME and the GM/WM set resolution above "
                  "match the same anatomy.".format(len(cord_untagged)))

        def tissue_type_of(label):
            if label in gm_element_labels:
                return 'GM'
            if label in wm_element_labels:
                return 'WM'
            return 'UNKNOWN'

        rows = []

        for frame_idx in range(len(step.frames)):
            frame = step.frames[frame_idx]

            le_subset   = frame.fieldOutputs['LE'].getSubset(region=instance, position=INTEGRATION_POINT)
            evol_subset = frame.fieldOutputs['EVOL'].getSubset(region=instance)

            # Take max MPS across integration points per element - same
            # CSDM-style convention as Alex_results_extraction.py (see that
            # script for the full rationale/reference).
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
                    rows.append((frame_idx, frame.frameValue, lbl, mps, elem_vol[lbl], tissue_type_of(lbl)))

        # Write CSV
        odb_basename = os.path.splitext(os.path.basename(ODB_PATH))[0]
        out_path = os.path.join(OUTPUT_DIR, '{}_mps_GM_WM.csv'.format(odb_basename))

        with open(out_path, 'w') as f:
            f.write('frame_index,frame_value,element_label,mps,volume,tissue_type\n')
            for r in rows:
                f.write('{},{:.6e},{},{:.6e},{:.6e},{}\n'.format(*r))

        print("Saved {} rows -> {}".format(len(rows), out_path))
    finally:
        odb.close()

# --- Batch summary - NOT auto-printed every call (that reprints the whole
# growing list on every single ODB_PATH in the loop, which floods the console
# across a large batch). Call print_batch_summary() yourself once, after the
# LAST ODB_PATH has run, to see it. The VERSION_MISMATCH block it prints is a
# complete, self-contained snippet (import through the execfile() loop) -
# paste it directly into the Abaqus 2022 kernel as-is, nothing else to add.
def print_batch_summary():
    if not (VERSION_MISMATCH_ODB_PATHS or OTHER_FAILED_ODB_PATHS):
        print("No failures.")
        return
    print("=== Batch summary ===")
    if VERSION_MISMATCH_ODB_PATHS:
        print("{} ODB(s) need a different Abaqus version - paste this whole block into the "
              "Abaqus 2022 kernel:".format(len(VERSION_MISMATCH_ODB_PATHS)))
        print("")
        print("import os")
        print("ODB_PATHS_2022 = [")
        for p in VERSION_MISMATCH_ODB_PATHS:
            print("    r'{}',".format(p))
        print("]")
        print("for ODB_PATH in ODB_PATHS_2022:")
        print("    execfile(r'{}')".format(SELF_PATH))
        print("")
    if OTHER_FAILED_ODB_PATHS:
        print("{} ODB(s) failed for other reasons (not a version mismatch):".format(
            len(OTHER_FAILED_ODB_PATHS)))
        for p in OTHER_FAILED_ODB_PATHS:
            print("    {}".format(p))
