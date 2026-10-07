"""
oscillation.py - Replace Step-2 (the flexion/extension moment) with a Step-3
that instead oscillates the non-cord anatomy (vertebrae/ligaments/bone) along
the C7-base -> C2-top axis, to model cord motion driven by CSF pulsation.

The C2-top and C7-base reference points are auto-detected per file via the
*Coupling lines (their *Nset names, e.g. m_Set-3/m_Set-4, vary by patient, but
the coupled surface names always contain 'c2'+'top' and 'c7'+'base'). The
C7-base RP is fixed at model level in every job checked and is kinematically
coupled to the whole vertebral/ligament surface, so driving it moves
"everything but the cord" while the cord itself is unaffected.

This reproduces (and generalises) the manual experiment already done for
patient N01-011 PreOp: Job-202-N01-011-PreOp-oscnomvt.inp. The Amp-3-osc
amplitude table and 0.76mm peak magnitude are reused verbatim from that file.

Never modifies the source .inp (only ever opened read-only) - writes into a
new, renumbered sibling job folder instead, marked '_removed' (Step-2 was
removed) and '_osc' (Step-3 oscillation added), e.g.:
  .../Pre-Op/Job-020-N01-011-PreOp-BC0pt35wEVOL/Job-020-...-BC0pt35wEVOL_0pt30_Site1_Site2_Site3_Site4.inp
  -> .../Pre-Op/Job-320-N01-011-PreOp-BC0pt35wEVOL_removed_osc/Job-320-...-BC0pt35wEVOL_removed_0pt30_Site1_Site2_Site3_Site4_osc.inp
Refuses to run if that folder or file already exists.

Set INP_PATH below to a specific .inp path for single-file mode, or to None
to run the full id_map.csv batch (every row where loading_condition ==
'Flexion'). Can also be run in the Abaqus CAE kernel via `execfile()`, or
via `abaqus python oscillation.py`.
"""
from __future__ import print_function

import csv
import math
import os
import re
import sys

# ============================================================
# CONSTANTS
# ============================================================
PEAK_OSC_MM = 0.76

# Step-3's *Static max-increment cap. Must not exceed the amplitude table's own
# sampling interval (0.02, see AMP_3_OSC_LINES below), or Abaqus could take a single
# increment large enough to skip over a peak/trough in the table entirely.
STEP3_MAX_INCREMENT = 0.02

# Copied verbatim from Job-202-N01-011-PreOp-oscnomvt.inp. Magnitude and
# table are from Sam Schaefer's work, which itself uses the oscillatory
# cervical cord motion profile described by Mikulis et al. [1,2]
#
# [1] Mikulis DJ, Wood ML, Zerdoner OAM, Poncelet BP. Oscillatory motion of
#     the normal cervical spinal cord. Radiology. 1994 Jul;192(1):117-121.
#     doi: 10.1148/radiology.192.1.8208922. PMID: 8208922.
# [2] Schaefer SD, Davies BM, Newcombe VFJ, Sutcliffe MPF. Could spinal cord
#     oscillation contribute to spinal cord injury in degenerative cervical
#     myelopathy? Brain and Spine. 2023;3:101743.
#     doi: 10.1016/j.bas.2023.101743. PMID: 37383476; PMCID: PMC10293319.
AMP_3_OSC_LINES = [
    "             0.,              0.,            0.02,         0.25723,            0.04,         0.46312,            0.06,          0.6143\n",
    "           0.08,         0.71078,             0.1,         0.75546,            0.12,         0.75357,            0.14,         0.71204\n",
    "           0.16,         0.63888,            0.18,          0.5426,             0.2,         0.43169,            0.22,         0.31417\n",
    "           0.24,         0.19723,            0.26,         0.08699,            0.28,         -0.0117,             0.3,        -0.09531\n",
    "           0.32,        -0.16161,            0.34,        -0.20963,            0.36,        -0.23951,            0.38,        -0.25232\n",
    "            0.4,        -0.24988,            0.42,        -0.23457,            0.44,        -0.20907,            0.46,        -0.17624\n",
    "           0.48,        -0.13888,             0.5,        -0.09964,            0.52,        -0.06088,            0.54,        -0.02457\n",
    "           0.56,         0.00771,            0.58,         0.03486,             0.6,         0.05618,              1.,              0.\n",
]

ID_MAP_PATH = r'C:\Users\cmb247\repos\Abaqus\DCM_Scripting\id_map.csv'

# ============================================================
# USER SETTINGS
# ============================================================
# Set to a specific .inp path to process just that one file (single-file
# mode). Set to None to instead run the full id_map.csv batch (every row
# where loading_condition == 'Flexion') - see run_batch() below.
#
# Currently hardcoded to the N01-011 Flexion PreOp job for the first test
# run, since a known-good manual result (Job-202-...-oscnomvt.inp) already
# exists for this patient to sanity-check the auto-detected reference
# points and direction vector against.
#
# >>> SET INP_PATH = None BELOW TO RUN THE FULL id_map.csv BATCH INSTEAD <<<
# INP_PATH = None
INP_PATH = r'D:\Charlotte\ABAQUS\N31-038\Job-103-N31-038-PostOpv2-BC0pt35\Job-103-N31-038-PostOpv2-BC0pt35.inp'

# INP_PATH = r'D:\Charlotte\ABAQUS\N01-014\\Job-018-N01-014-Pre-Opv10-BC0pt35\Job-018-N01-014-Pre-Opv10-BC0pt35_0pt30_Site1_Site2.inp'

#r'D:\Charlotte\ABAQUS\N01-011\Pre-Op\Job-020-N01-011-PreOp-BC0pt35wEVOL\Job-020-N01-011-PreOp-BC0pt35wEVOL_0pt30_Site1_Site2_Site3_Site4.inp'


# ============================================================
# .inp PARSING HELPERS
# ============================================================
def _parse_data_labels(lines, start_idx):
    """Read comma-separated integer labels from lines[start_idx:] until a line starting with '*'."""
    labels = []
    i = start_idx
    while i < len(lines) and not lines[i].lstrip().startswith('*'):
        parts = [p.strip() for p in lines[i].strip().rstrip(',').split(',') if p.strip()]
        labels.extend(int(p) for p in parts)
        i += 1
    return labels, i


def _parse_boundary_data(lines, start_idx):
    """Read data lines of a *Boundary block until a line starting with '*'."""
    entries = []
    i = start_idx
    while i < len(lines) and not lines[i].lstrip().startswith('*'):
        parts = [p.strip() for p in lines[i].strip().rstrip(',').split(',') if p.strip()]
        if parts:
            entries.append(tuple(parts))
        i += 1
    return entries, i


def _find_enclosing_instance(lines, line_idx):
    """Scan backward from line_idx to find the nearest still-open '*Instance, name=X'."""
    depth = 0
    for i in range(line_idx - 1, -1, -1):
        stripped = lines[i].strip()
        if re.match(r'\*End Instance', stripped, re.IGNORECASE):
            depth += 1
        elif re.match(r'\*Instance\s*,', stripped, re.IGNORECASE):
            if depth == 0:
                m = re.search(r'name=([^\s,]+)', stripped, re.IGNORECASE)
                return m.group(1) if m else None
            depth -= 1
    return None


def _resolve_nset_node(lines, nset_name):
    """Find '*Nset, nset=<nset_name>' and return (instance_name_or_None, node_label)."""
    pattern = re.compile(r'\*Nset\s*,\s*nset=' + re.escape(nset_name) + r'\b', re.IGNORECASE)
    matches = [i for i, line in enumerate(lines) if pattern.match(line.strip())]
    if len(matches) != 1:
        raise RuntimeError("Expected exactly one '*Nset, nset={0}' line, found {1}".format(
            nset_name, len(matches)))
    idx = matches[0]
    header = lines[idx]
    inst_m = re.search(r'instance=([^\s,]+)', header, re.IGNORECASE)
    instance_name = inst_m.group(1) if inst_m else _find_enclosing_instance(lines, idx)
    labels, _ = _parse_data_labels(lines, idx + 1)
    if len(labels) != 1:
        raise RuntimeError("Expected '{0}' to be a single-node (RP) set, found {1} node(s)".format(
            nset_name, len(labels)))
    return instance_name, labels[0]


def _part_block_ranges(lines):
    """Return [(start, end), ...] line index ranges for every *Part ... *End Part block.
    Used to exclude a Part's own internal node numbering when resolving a bare
    assembly-level node (instance_name is None) - a Part's node 1 is a different
    namespace from an assembly-level node 1."""
    ranges = []
    starts = [i for i, line in enumerate(lines) if re.match(r'\*Part\s*,', line.strip(), re.IGNORECASE)]
    for s in starts:
        for i in range(s + 1, len(lines)):
            if re.match(r'\*End Part\b', lines[i].strip(), re.IGNORECASE):
                ranges.append((s, i))
                break
    return ranges


def _find_node_coords(lines, instance_name, node_label):
    """Find (x, y, z) for node_label.

    If instance_name is given, searches only that *Instance ... *End Instance block's
    own *Node data (e.g. an orphan-mesh instance with its own embedded nodes).

    If instance_name is None, the *Nset had no 'instance=' and wasn't inside any open
    *Instance block, meaning it refers to a node defined directly in the *Assembly
    block itself (a standalone reference point, not part of any instanced mesh). That
    search must exclude every *Node block belonging to a *Part definition - those are
    a completely separate numbering namespace and would otherwise be matched first
    simply because Part definitions come first in the file, silently returning the
    wrong node.
    """
    if instance_name is not None:
        inst_start = None
        inst_re = re.compile(r'\*Instance\s*,.*\bname=' + re.escape(instance_name) + r'\b', re.IGNORECASE)
        for i, line in enumerate(lines):
            if inst_re.match(line.strip()):
                inst_start = i
                break
        if inst_start is None:
            raise RuntimeError("Could not find '*Instance, name={0}'".format(instance_name))
        inst_end = len(lines)
        for i in range(inst_start + 1, len(lines)):
            if re.match(r'\*End Instance', lines[i].strip(), re.IGNORECASE):
                inst_end = i
                break
        candidate_starts = [i for i in range(inst_start, inst_end)
                             if re.match(r'\*Node\b', lines[i].strip(), re.IGNORECASE)]
    else:
        part_ranges = _part_block_ranges(lines)
        candidate_starts = [i for i, line in enumerate(lines)
                             if re.match(r'\*Node\b', line.strip(), re.IGNORECASE)
                             and not any(s <= i <= e for s, e in part_ranges)]

    for nb in candidate_starts:
        i = nb + 1
        while i < len(lines) and not lines[i].lstrip().startswith('*'):
            parts = [p.strip() for p in lines[i].strip().rstrip(',').split(',') if p.strip()]
            if parts and int(parts[0]) == node_label:
                coords = [float(x) for x in parts[1:4]]
                while len(coords) < 3:
                    coords.append(0.0)
                return tuple(coords)
            i += 1
    raise RuntimeError("Could not find coordinates for node {0} in instance {1}".format(
        node_label, instance_name))


def _detect_rps(lines):
    """Auto-detect the upper ('cN-top') and lower ('cN-base') coupling reference points
    via *Coupling lines. Every model spans C2 to C7, so the surface should always be
    named 'c2-top'/'c7-base' - but at least two source files mislabel the base surface
    as 'c3-base' (a naming typo upstream, not a real anatomical difference). Matching
    any vertebra number rather than hardcoding C2/C7 means the detection still works
    despite that typo, without needing to special-case it."""
    coupling_re = re.compile(r'\*Coupling\s*,.*ref node=([^\s,]+)\s*,\s*surface=([^\s,]+)', re.IGNORECASE)
    upper, lower = None, None
    for line in lines:
        m = coupling_re.match(line.strip())
        if not m:
            continue
        ref_node, surface = m.group(1), m.group(2).lower()
        if re.search(r'c\d+-top', surface):
            if upper is not None:
                raise RuntimeError("Multiple candidate '-top' coupling reference points found")
            upper = ref_node
        elif re.search(r'c\d+-base', surface):
            if lower is not None:
                raise RuntimeError("Multiple candidate '-base' coupling reference points found")
            lower = ref_node
    if upper is None or lower is None:
        raise RuntimeError("Could not auto-detect both 'cN-top' and 'cN-base' coupling reference points "
                            "(found upper={0}, lower={1})".format(upper, lower))
    return upper, lower


def _find_step_spans(lines):
    """Return [(start_idx, end_idx), ...] for each *Step...*End Step block, in file order."""
    step_starts = [i for i, line in enumerate(lines) if re.match(r'\*Step\s*,', line.strip(), re.IGNORECASE)]
    spans = []
    for s in step_starts:
        end = None
        for i in range(s + 1, len(lines)):
            if re.match(r'\*End Step', lines[i].strip(), re.IGNORECASE):
                end = i
                break
        if end is None:
            raise RuntimeError("Could not find matching '*End Step' for step starting at line {0}".format(s + 1))
        spans.append((s, end))
    return spans


def _check_boundary_conditions(lines, step1_span, rp_lower):
    """Refuse to proceed unless the only standing BC is rp_lower's pre-Step-1 fix."""
    step1_start, step1_end = step1_span
    boundary_starts = [i for i, line in enumerate(lines[:step1_start])
                        if re.match(r'\*Boundary\b', line.strip(), re.IGNORECASE)]
    entries = []
    for b in boundary_starts:
        data, _ = _parse_boundary_data(lines, b + 1)
        entries.extend(data)
    # A single RP's fix can be written as one 'ENCASTRE' line or as several
    # per-DOF lines (e.g. 1,1 / 2,2 / ... / 6,6) - both are fine, as long as
    # every entry belongs to rp_lower and no other node set has a standing BC.
    if not entries or any(e[0].lower() != rp_lower.lower() for e in entries):
        raise RuntimeError(
            "Expected all pre-Step-1 boundary conditions to be on '{0}' only; found {1}. "
            "Refusing to reset boundary conditions in Step-3 without checking this file "
            "manually first (a bare '*Boundary, op=NEW' would silently drop these).".format(
                rp_lower, entries))
    within_step1 = [i for i in range(step1_start + 1, step1_end)
                    if re.match(r'\*Boundary\b', lines[i].strip(), re.IGNORECASE)]
    if within_step1:
        raise RuntimeError(
            "Step-1 defines its own *Boundary condition(s) at line(s) {0}; the Step-3 "
            "'op=NEW' reset would need to re-assert these but this script does not - "
            "aborting rather than silently dropping them.".format([i + 1 for i in within_step1]))


def _amplitude_insertion_point(lines, step1_start):
    """Insert point right after the last existing pre-Step-1 '*Amplitude' block, so the
    new amplitude sits with the others instead of wherever happens to be right before
    Step-1's comment header."""
    amp_starts = [i for i, line in enumerate(lines[:step1_start])
                  if re.match(r'\*Amplitude\b', line.strip(), re.IGNORECASE)]
    if not amp_starts:
        idx = step1_start
        while idx > 0 and lines[idx - 1].strip().startswith('**'):
            idx -= 1
        return idx
    i = amp_starts[-1] + 1
    while i < step1_start and not lines[i].lstrip().startswith('*'):
        i += 1
    return i


def _renumber_job(name):
    """Replace the leading digit of the job number in a 'Job-NNN...' name with '3'.

    e.g. 'Job-020-N01-011-PreOp-BC0pt35wEVOL' -> 'Job-320-N01-011-PreOp-BC0pt35wEVOL'
    """
    m = re.match(r'(Job-)(\d+)(.*)', name)
    if not m:
        raise RuntimeError("Could not find a 'Job-<number>' prefix in '{0}'".format(name))
    prefix, number, rest = m.groups()
    new_number = '3' + number[1:] if len(number) > 1 else '3'
    return prefix + new_number + rest


def _insert_removed_marker(name):
    """Insert '_removed' right before the predefined-field/site suffix (e.g.
    '_0pt30_Site1_Site2...'), so it reads '...BC0pt35wEVOL_removed_0pt30_Site...'.
    Falls back to appending '_removed' at the end if no such suffix is present
    (e.g. PostOp jobs, which have no compression-site predefined fields)."""
    m = re.search(r'_\d+pt\d+_Site', name, re.IGNORECASE)
    if m:
        idx = m.start()
        return name[:idx] + '_removed' + name[idx:]
    return name + '_removed'


def _output_path(inp_path):
    """New job folder (sibling of the source job's folder) + renumbered,
    '_removed_osc'-marked file, signalling that Step-2 was removed.

    e.g. .../Pre-Op/Job-020-N01-011-PreOp-BC0pt35wEVOL/Job-020-...-BC0pt35wEVOL_0pt30_Site1_Site2_Site3_Site4.inp
      -> .../Pre-Op/Job-320-N01-011-PreOp-BC0pt35wEVOL_removed_osc/Job-320-...-BC0pt35wEVOL_removed_0pt30_Site1_Site2_Site3_Site4_osc.inp
    """
    orig_dir = os.path.dirname(inp_path)
    orig_folder_name = os.path.basename(orig_dir)
    parent_dir = os.path.dirname(orig_dir)

    file_base, file_ext = os.path.splitext(os.path.basename(inp_path))

    new_folder_name = _insert_removed_marker(_renumber_job(orig_folder_name)) + '_osc'
    new_file_name = _insert_removed_marker(_renumber_job(file_base)) + '_osc' + file_ext

    return os.path.join(parent_dir, new_folder_name, new_file_name)


# ============================================================
# CORE TRANSFORM
# ============================================================
def add_oscillation_step(inp_path):
    with open(inp_path, 'r') as f:
        lines = f.readlines()

    rp_upper, rp_lower = _detect_rps(lines)
    inst_upper, label_upper = _resolve_nset_node(lines, rp_upper)
    inst_lower, label_lower = _resolve_nset_node(lines, rp_lower)
    coord_upper = _find_node_coords(lines, inst_upper, label_upper)
    coord_lower = _find_node_coords(lines, inst_lower, label_lower)

    vec = [coord_upper[k] - coord_lower[k] for k in range(3)]
    mag = math.sqrt(sum(c * c for c in vec))
    if mag < 1e-9:
        raise RuntimeError("Degenerate direction vector between {0} and {1}".format(rp_upper, rp_lower))
    unit = [c / mag for c in vec]
    dx, dy, dz = (PEAK_OSC_MM * c for c in unit)

    print("Job: {0}".format(inp_path))
    print("  C2-top RP  : {0} (instance={1}, node={2}) @ {3}".format(rp_upper, inst_upper, label_upper, coord_upper))
    print("  C7-base RP : {0} (instance={1}, node={2}) @ {3}  <- oscillation driver".format(
        rp_lower, inst_lower, label_lower, coord_lower))
    print("  Unit vector (C7base->C2top): ({0:.5f}, {1:.5f}, {2:.5f})".format(*unit))
    print("  Oscillation displacement components (mm, peak={0}): dx={1:.5f}, dy={2:.5f}, dz={3:.5f}".format(
        PEAK_OSC_MM, dx, dy, dz))

    step_spans = _find_step_spans(lines)
    if len(step_spans) != 2:
        raise RuntimeError(
            "Expected exactly 2 steps (Step-1, Step-2); found {0}. This file may already "
            "have a Step-3, or has an unexpected structure - check manually.".format(len(step_spans)))
    step1_span, step2_span = step_spans
    step1_start, step1_end = step1_span
    _, step2_end = step2_span

    _check_boundary_conditions(lines, step1_span, rp_lower)

    # Reuse Step-1's *Static solver-control keyword line, but give Step-3 its own
    # max-increment value: the oscillation amplitude table samples every 0.02 (in
    # normalized step time), so the max increment must not exceed that, or Abaqus
    # could legitimately take a single increment large enough to skip a peak/trough
    # in the table entirely rather than resolving the oscillation shape.
    static_keyword_line = lines[step1_start + 1]
    if not static_keyword_line.strip().startswith('*Static'):
        raise RuntimeError("Expected '*Static' as the line following Step-1's '*Step' keyword")
    step1_static_data = [p.strip() for p in lines[step1_start + 2].strip().rstrip(',').split(',')]
    if len(step1_static_data) != 4:
        raise RuntimeError("Expected Step-1's *Static data line to have 4 values "
                            "(initial, total, min, max increment), found {0}".format(step1_static_data))
    initial_inc, total_time, min_inc, _ = step1_static_data
    static_data_line = "{0}, {1}, {2}, {3}\n".format(initial_inc, total_time, min_inc, STEP3_MAX_INCREMENT)

    # Reuse Step-1's output-request block, including its '** OUTPUT REQUESTS' comment
    # header (from that header, or '*Restart' itself if there's no header, through to
    # '*End Step')
    output_start = None
    for i in range(step1_start, step1_end):
        if re.match(r'\*Restart\b', lines[i].strip(), re.IGNORECASE):
            output_start = i
            break
    if output_start is None:
        raise RuntimeError("Could not find '*Restart' (start of output requests) inside Step-1")
    while output_start > step1_start and lines[output_start - 1].strip().startswith('**'):
        output_start -= 1
    output_block = lines[output_start:step1_end]

    step3_lines = [
        "** ----------------------------------------------------------------\n",
        "** \n",
        "** STEP: Step-3\n",
        "** \n",
        "*Step, name=Step-3, nlgeom=YES\n",
        static_keyword_line,
        static_data_line,
        "** \n",
        "** BOUNDARY CONDITIONS\n",
        "** \n",
        "** Name: BC-1 Type: Displacement/Rotation\n",
        "*Boundary, op=NEW\n",
        "** Name: BC-3-osc Type: Displacement/Rotation\n",
        "*Boundary, op=NEW, amplitude=Amp-3-osc\n",
        "{0}, 1, 1, {1:.6f}\n".format(rp_lower, dx),
        "{0}, 2, 2, {1:.6f}\n".format(rp_lower, dy),
        "{0}, 3, 3, {1:.6f}\n".format(rp_lower, dz),
    ]
    step3_lines.extend(output_block)  # output_block's own leading '**' comments supply the separator
    step3_lines.append("*End Step\n")

    amp_block = ["*Amplitude, name=Amp-3-osc\n"] + list(AMP_3_OSC_LINES)
    amp_idx = _amplitude_insertion_point(lines, step1_start)

    out_lines = []
    out_lines.extend(lines[:amp_idx])
    out_lines.extend(amp_block)
    out_lines.extend(lines[amp_idx:step1_end + 1])  # unchanged, including Step-1's '*End Step'
    out_lines.extend(step3_lines)
    out_lines.extend(lines[step2_end + 1:])  # anything after Step-2's '*End Step' (usually just EOF)

    out_path = _output_path(inp_path)
    if os.path.abspath(out_path) == os.path.abspath(inp_path):
        raise RuntimeError("Refusing to write: output path resolved to the source file itself ({0})".format(inp_path))

    # out_path/out_dir are always this script's own generated '..._removed_osc' output
    # (never the source .inp, guarded above), so overwriting them to regenerate with a
    # fixed script version is fine - the source .inp is still never touched.
    out_dir = os.path.dirname(out_path)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    with open(out_path, 'w') as f:
        f.writelines(out_lines)
    print("  Wrote (new folder, source untouched): {0}".format(out_path))
    return out_path


# ============================================================
# BATCH MODE (id_map.csv, loading_condition == 'Flexion')
# ============================================================
def _read_id_map():
    if sys.version_info[0] >= 3:
        read_mode, kwargs = 'r', {'newline': ''}
    else:
        read_mode, kwargs = 'rb', {}
    with open(ID_MAP_PATH, read_mode, **kwargs) as f:
        return list(csv.DictReader(f))


def run_batch():
    rows = _read_id_map()
    flexion_rows = [r for r in rows if r.get('loading_condition', '').strip().lower() == 'flexion']
    print("Found {0} Flexion row(s) in id_map.csv".format(len(flexion_rows)))

    results = []
    for row in flexion_rows:
        odb_path = row.get('odb_path', '').strip()
        label = "P{0} {1}".format(row.get('participant', '?'), row.get('State', '?'))
        if not odb_path:
            results.append((label, 'SKIPPED', 'no odb_path'))
            continue
        inp_path = os.path.splitext(odb_path)[0] + '.inp'
        if not os.path.isfile(inp_path):
            results.append((label, 'SKIPPED', 'inp not found: {0}'.format(inp_path)))
            continue
        try:
            out_path = add_oscillation_step(inp_path)
            results.append((label, 'OK', out_path))
        except Exception as e:
            results.append((label, 'FAILED', str(e)))
        print("")

    print("=== Summary ===")
    for label, status, detail in results:
        print("{0:20s} | {1:8s} | {2}".format(label, status, detail))


# ============================================================
if 'INP_PATH' in dir() and INP_PATH:
    add_oscillation_step(INP_PATH)
else:
    run_batch()
