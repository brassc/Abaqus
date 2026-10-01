"""
Alex_results_extraction_IVD.py - Extract MPS (max principal strain) and
element volume (EVOL) from the 'P27;IVD' element set (intervertebral discs -
one combined set, not split per level), for every frame of the step. Outputs
one CSV per ODB to the same directory as the ODB.

Based on Alex_results_extraction.py, retargeted from the Cord NODE set to the
IVD ELEMENT set - IVD is already an element set (*Elset, elset=P27;IVD) so,
unlike Cord, no node->element derivation step is needed.

Separate script, not a flag on Alex_results_extraction.py: different set,
different output suffix (_ivd_mps.csv / _ivd_topology.csv, so it never
collides with the existing Cord *_mps.csv/_topology.csv), and does NOT touch
id_map.csv (its 'csv_path' column is Cord-specific; there is no IVD column).

Must run under both Abaqus 2022 (Python 2.7 kernel) and Abaqus 2025 (Python
3.10 kernel), since these will be batched separately - kept to syntax valid
on both (.format() not f-strings, print() function, list(...values())
instead of indexing a Repository's .values() directly since dict_values
isn't subscriptable on Python 3).

If an ODB was written by a different Abaqus version than the one running
this kernel, openOdb() fails for just that file - the script catches it,
records ODB_PATH, and returns normally instead of raising, so a
`for ODB_PATH in ODB_PATHS: execfile(...)` batch loop moves on to the next
ODB_PATH instead of dying on the first mismatch. See VERSION_MISMATCH_ODB_PATHS
/ OTHER_FAILED_ODB_PATHS and the batch summary printed at the bottom.

Run in Abaqus CAE kernel:

ODB_PATH = r'D:\\path\\to\\Job-xxx.odb'
execfile('C:\\Users\\cmb247\\repos\\Abaqus\\DCM_Scripting\\Alex_results_extraction_IVD.py')

Or from the command line:
abaqus python Alex_results_extraction_IVD.py <odb_path>
"""

from odbAccess import *
from abaqusConstants import *
import os
import sys

# ============================================================
# USER SETTINGS (can be overridden by setting variables before execfile())
# ============================================================
if 'ODB_PATH' not in dir():
    if len(sys.argv) > 1:
        ODB_PATH = sys.argv[1]
    else:
        raise RuntimeError("Set ODB_PATH before running this script.")
if 'IVD_SET_NAME' not in dir():
    IVD_SET_NAME = 'P27;IVD'
if 'STEP_NAME' not in dir():
    STEP_NAME = None   # None = last step
# ============================================================

OUTPUT_DIR = os.path.dirname(ODB_PATH)

# Batch-loop bookkeeping - persists across repeated execfile() calls in the
# same kernel session (each loop iteration re-execfile's this whole script,
# so a plain assignment here would otherwise reset the list every time).
#
# Two DIFFERENT version-mismatch directions need two DIFFERENT fixes, so they
# are tracked separately rather than lumped into one "VERSION_MISMATCH" bucket
# (which previously always suggested "re-run under Abaqus 2022" regardless of
# which direction the mismatch actually was - wrong and misleading whenever
# the file needed something OTHER than exactly 2022):
#   NEWER_REQUIRED - the ODB is from a MORE RECENT Abaqus release than this
#     kernel ("...must be upgraded before this output database can be
#     opened"). No command-line fix - -upgrade only converts old-format files
#     to new, never the reverse. Needs opening in an actually newer Abaqus
#     install.
#   UPGRADE_NEEDED - the ODB is from a PREVIOUS (older) Abaqus release
#     ("Run abaqus -upgrade ... to upgrade it"). Abaqus's own error message
#     already names the exact fix: `abaqus -upgrade -job <newFileName> -odb
#     <oldOdbFileName>` converts it in place - no older kernel needed at all.
if 'NEWER_REQUIRED_ODB_PATHS' not in dir():
    NEWER_REQUIRED_ODB_PATHS = []
if 'UPGRADE_NEEDED_ODB_PATHS' not in dir():
    UPGRADE_NEEDED_ODB_PATHS = []
if 'OTHER_FAILED_ODB_PATHS' not in dir():
    OTHER_FAILED_ODB_PATHS = []

# Open the ODB defensively: an ODB written by a different Abaqus version than
# the one running this kernel raises here (wrong release for this ODB), and
# without a try/except that exception would propagate out of execfile() and
# kill the REST of the batch loop, not just this one ODB_PATH. Catch it,
# classify which DIRECTION the mismatch is (see above), and let the script
# finish normally so the loop moves on to the next ODB_PATH.
odb = None
try:
    odb = openOdb(path=ODB_PATH, readOnly=True)
except Exception as e:
    msg = str(e)
    msg_lower = msg.lower()
    if 'more recent' in msg_lower:
        print("SKIPPING (ODB needs a NEWER Abaqus release than this kernel): {}".format(ODB_PATH))
        NEWER_REQUIRED_ODB_PATHS.append(ODB_PATH)
    elif 'previous release' in msg_lower or 'upgrade' in msg_lower:
        print("SKIPPING (ODB is from an OLDER Abaqus release - needs 'abaqus -upgrade', "
              "not a different kernel): {}".format(ODB_PATH))
        UPGRADE_NEEDED_ODB_PATHS.append(ODB_PATH)
    else:
        print("SKIPPING (failed to open ODB): {}".format(ODB_PATH))
        OTHER_FAILED_ODB_PATHS.append(ODB_PATH)
    print("  {}".format(msg))

if odb is not None:
    try:
        step = odb.steps[STEP_NAME] if STEP_NAME else list(odb.steps.values())[-1]
        print("Step: '{}' ({} frames)".format(step.name, len(step.frames)))

        # 'P27;IVD' is already an ELEMENT set (unlike Cord, which is a node
        # set that Alex_results_extraction.py derives elements from) - look
        # it up directly at assembly or instance level, same resolution
        # order as CORD_SET_NAME there.
        assembly = odb.rootAssembly
        if IVD_SET_NAME in assembly.elementSets.keys():
            elem_set = assembly.elementSets[IVD_SET_NAME]
            if len(elem_set.instances) != 1:
                raise RuntimeError("Expected IVD element set to span exactly one instance; found {}.".format(
                    len(elem_set.instances)))
            instance = elem_set.instances[0]
            ivd_element_labels = set(e.label for e in elem_set.elements[0])
        else:
            if '.' in IVD_SET_NAME:
                inst_name, sname = IVD_SET_NAME.split('.', 1)
                instance = assembly.instances[inst_name]
            else:
                matches = [inst for inst in assembly.instances.values() if IVD_SET_NAME in inst.elementSets.keys()]
                if len(matches) == 1:
                    instance = matches[0]
                elif len(matches) > 1:
                    raise RuntimeError("IVD_SET_NAME '{}' found in multiple instances; "
                                        "specify as 'Instance.SetName'.".format(IVD_SET_NAME))
                else:
                    raise RuntimeError("IVD_SET_NAME '{}' not found at assembly or "
                                        "instance level (element set).".format(IVD_SET_NAME))
                sname = IVD_SET_NAME
            elem_set = instance.elementSets[sname]
            ivd_element_labels = set(e.label for e in elem_set.elements)

        print("IVD element set '{}': {} elements in instance '{}'".format(
            IVD_SET_NAME, len(ivd_element_labels), instance.name))

        rows = []

        for frame_idx in range(len(step.frames)):
            frame = step.frames[frame_idx]

            le_subset   = frame.fieldOutputs['LE'].getSubset(region=instance, position=INTEGRATION_POINT)
            evol_subset = frame.fieldOutputs['EVOL'].getSubset(region=instance)

            # Take max MPS across integration points per element - same
            # convention as Alex_results_extraction.py (see that script for
            # the CSDM rationale).
            elem_mps = {}
            for val in le_subset.values:
                lbl = val.elementLabel
                if lbl not in ivd_element_labels:
                    continue
                if lbl not in elem_mps or val.maxPrincipal > elem_mps[lbl]:
                    elem_mps[lbl] = val.maxPrincipal

            elem_vol = {val.elementLabel: val.data for val in evol_subset.values
                        if val.elementLabel in ivd_element_labels}

            for lbl, mps in elem_mps.items():
                if lbl in elem_vol:
                    rows.append((frame_idx, frame.frameValue, lbl, mps, elem_vol[lbl]))

        # Write CSV
        odb_basename = os.path.splitext(os.path.basename(ODB_PATH))[0]
        out_path = os.path.join(OUTPUT_DIR, '{}_ivd_mps.csv'.format(odb_basename))

        with open(out_path, 'w') as f:
            f.write('frame_index,frame_value,element_label,mps,volume\n')
            for r in rows:
                f.write('{},{:.6e},{},{:.6e},{:.6e}\n'.format(*r))

        print("Saved {} rows -> {}".format(len(rows), out_path))

        # --- Element connectivity (topology), for spatial clustering
        # analysis downstream - same purpose/format as
        # Alex_results_extraction.py's _topology.csv, just for the IVD
        # element set instead of Cord.
        topology_path = os.path.join(OUTPUT_DIR, '{}_ivd_topology.csv'.format(odb_basename))
        with open(topology_path, 'w') as f:
            f.write('element_label,nodes\n')
            for elem in instance.elements:
                if elem.label in ivd_element_labels:
                    nodes = ';'.join(str(n) for n in elem.connectivity)
                    f.write('{},{}\n'.format(elem.label, nodes))
        print("Saved topology for {} elements -> {}".format(len(ivd_element_labels), topology_path))
    finally:
        odb.close()

# --- Batch summary - NOT auto-printed every call (that reprints the whole
# growing list on every single ODB_PATH in the loop, which floods the console
# across a large batch). Call print_batch_summary() yourself once, after the
# LAST ODB_PATH has run, to see it.
def print_batch_summary():
    if not (NEWER_REQUIRED_ODB_PATHS or UPGRADE_NEEDED_ODB_PATHS or OTHER_FAILED_ODB_PATHS):
        print("No failures.")
        return
    print("=== Batch summary ===")
    if NEWER_REQUIRED_ODB_PATHS:
        print("{} ODB(s) need a NEWER Abaqus release than this kernel - there is no command-line "
              "fix (abaqus -upgrade only converts old format to new, never the reverse); open "
              "these in whatever newer Abaqus install is available:".format(
                  len(NEWER_REQUIRED_ODB_PATHS)))
        for p in NEWER_REQUIRED_ODB_PATHS:
            print("    {}".format(p))
    if UPGRADE_NEEDED_ODB_PATHS:
        print("{} ODB(s) are from an OLDER Abaqus release - fix with 'abaqus -upgrade', NOT a "
              "different kernel (run from the command line, not inside a kernel session):".format(
                  len(UPGRADE_NEEDED_ODB_PATHS)))
        for p in UPGRADE_NEEDED_ODB_PATHS:
            print("    abaqus -upgrade -job {}_upgraded -odb \"{}\"".format(
                os.path.splitext(os.path.basename(p))[0], p))
    if OTHER_FAILED_ODB_PATHS:
        print("{} ODB(s) failed for other reasons (not a version mismatch):".format(
            len(OTHER_FAILED_ODB_PATHS)))
        for p in OTHER_FAILED_ODB_PATHS:
            print("    {}".format(p))
