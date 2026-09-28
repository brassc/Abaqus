"""
remove_preload.py - Copy a PreOp job's .inp and remove its compression-site preload
(the '** PREDEFINED FIELDS' / *Temperature blocks in Step-1), leaving Step-2 (the
flexion or extension moment) untouched. Matches the shape PostOp jobs already have
naturally, since PostOp Step-1 bodies never have a predefined-fields section at all
(*Static's data line is followed straight by '** INTERACTIONS').

The *Amplitude, name=Amp-1-preload definition itself is left in the model data even
though nothing references it anymore - that's deliberate, not an oversight.

Never modifies the source .inp (only ever opened read-only) - writes into a new,
renumbered sibling job folder instead, with the compression-site suffix dropped from
the filename (it no longer applies once the preload driving those sites is gone):
  .../Pre-Op/Job-020-N01-011-PreOp-BC0pt35wEVOL/Job-020-...-BC0pt35wEVOL_0pt30_Site1_Site2_Site3_Site4.inp
  -> .../Pre-Op/Job-420-N01-011-PreOp-BC0pt35wEVOL/Job-420-N01-011-PreOp-BC0pt35wEVOL.inp
Refuses to run if that folder or file already exists.

Set INP_PATH below to a specific .inp path for single-file mode, or to None to run
the full id_map.csv batch (every PreOp row, both Flexion and Extension). Can also be
run in the Abaqus CAE kernel via execfile(), or via `abaqus python remove_preload.py`.
"""
from __future__ import print_function

import csv
import os
import re
import sys

ID_MAP_PATH = r'C:\Users\cmb247\repos\Abaqus\DCM_Scripting\id_map.csv'

# ============================================================
# USER SETTINGS
# ============================================================
# Set to a specific .inp path to process just that one file (single-file mode).
# Set to None to instead run the full id_map.csv batch (every row where State ==
# 'PreOp', both Flexion and Extension).
#
# >>> SET INP_PATH = None BELOW TO RUN THE FULL id_map.csv BATCH INSTEAD <<<
INP_PATH = None# r'D:\Charlotte\ABAQUS\N01-011\Pre-Op\Job-020-N01-011-PreOp-BC0pt35wEVOL\Job-020-N01-011-PreOp-BC0pt35wEVOL_0pt30_Site1_Site2_Site3_Site4.inp'


# ============================================================
# .inp PARSING HELPERS
# ============================================================
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


def _find_predefined_fields_span(lines, step1_span):
    """Find the ['** PREDEFINED FIELDS' header ... last SITE data line] span inside
    Step-1, including any blank '**' separator lines immediately before the header.
    Returns (start_idx, end_idx) with end_idx exclusive."""
    step1_start, step1_end = step1_span
    header_idx = None
    for i in range(step1_start, step1_end):
        if re.match(r'\*\*\s*PREDEFINED FIELDS\b', lines[i].strip(), re.IGNORECASE):
            header_idx = i
            break
    if header_idx is None:
        raise RuntimeError("No '** PREDEFINED FIELDS' section found inside Step-1 - "
                            "nothing to remove (unexpected for a PreOp job).")

    # Consume triplets: '** Name: predefinedfield-...' / '*Temperature, ...' / one data
    # line, for as long as the pattern holds. Not hardcoded to a specific next-section
    # name - stops at the first line that doesn't look like part of the block.
    i = header_idx + 1
    while i + 1 < step1_end:
        comment_line = lines[i].strip()
        keyword_line = lines[i + 1].strip()
        if not re.match(r'\*\*\s*Name:\s*predefinedfield-', comment_line, re.IGNORECASE):
            break
        if not keyword_line.startswith('*Temperature'):
            break
        i += 2
        if i < step1_end and not lines[i].lstrip().startswith('*'):
            i += 1  # the one data line following *Temperature
        else:
            break
    end_idx = i

    start_idx = header_idx
    while start_idx > step1_start and lines[start_idx - 1].strip() == '**':
        start_idx -= 1

    return start_idx, end_idx


def _renumber_job(name):
    """Replace the leading digit of the job number in a 'Job-NNN...' name with '4'.

    e.g. 'Job-020-N01-011-PreOp-BC0pt35wEVOL' -> 'Job-420-N01-011-PreOp-BC0pt35wEVOL'
    """
    m = re.match(r'(Job-)(\d+)(.*)', name)
    if not m:
        raise RuntimeError("Could not find a 'Job-<number>' prefix in '{0}'".format(name))
    prefix, number, rest = m.groups()
    new_number = '4' + number[1:] if len(number) > 1 else '4'
    return prefix + new_number + rest


def _strip_site_suffix(name):
    """Remove the compression-site suffix (e.g. '_0pt30_Site1_Site2...') from a file
    basename, since it no longer applies once the preload driving those sites is gone."""
    m = re.search(r'_\d+pt\d+_Site.*$', name, re.IGNORECASE)
    return name[:m.start()] if m else name


def _output_path(inp_path):
    """New job folder (sibling of the source job's folder) + renumbered file, with the
    compression-site suffix stripped from the filename.

    e.g. .../Pre-Op/Job-020-N01-011-PreOp-BC0pt35wEVOL/Job-020-...-BC0pt35wEVOL_0pt30_Site1_Site2_Site3_Site4.inp
      -> .../Pre-Op/Job-420-N01-011-PreOp-BC0pt35wEVOL/Job-420-N01-011-PreOp-BC0pt35wEVOL.inp
    """
    orig_dir = os.path.dirname(inp_path)
    orig_folder_name = os.path.basename(orig_dir)
    parent_dir = os.path.dirname(orig_dir)

    file_base, file_ext = os.path.splitext(os.path.basename(inp_path))

    new_folder_name = _renumber_job(orig_folder_name)
    new_file_name = _strip_site_suffix(_renumber_job(file_base)) + file_ext

    return os.path.join(parent_dir, new_folder_name, new_file_name)


# ============================================================
# CORE TRANSFORM
# ============================================================
def remove_preload(inp_path):
    with open(inp_path, 'r') as f:
        lines = f.readlines()

    step_spans = _find_step_spans(lines)
    if len(step_spans) != 2:
        raise RuntimeError(
            "Expected exactly 2 steps (Step-1, Step-2); found {0}. Check this file's "
            "structure manually.".format(len(step_spans)))
    step1_span, _ = step_spans

    pf_start, pf_end = _find_predefined_fields_span(lines, step1_span)

    print("Job: {0}".format(inp_path))
    print("  Removing '** PREDEFINED FIELDS' block: lines {0}-{1} ({2} lines)".format(
        pf_start + 1, pf_end, pf_end - pf_start))

    out_lines = lines[:pf_start] + lines[pf_end:]

    out_path = _output_path(inp_path)
    if os.path.abspath(out_path) == os.path.abspath(inp_path):
        raise RuntimeError("Refusing to write: output path resolved to the source file itself ({0})".format(inp_path))

    out_dir = os.path.dirname(out_path)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    with open(out_path, 'w') as f:
        f.writelines(out_lines)
    print("  Wrote (new folder, source untouched): {0}".format(out_path))
    return out_path


# ============================================================
# BATCH MODE (id_map.csv, State == 'PreOp', both loading conditions)
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
    preop_rows = [r for r in rows if r.get('State', '').strip().lower() == 'preop']
    print("Found {0} PreOp row(s) in id_map.csv".format(len(preop_rows)))

    results = []
    for row in preop_rows:
        odb_path = row.get('odb_path', '').strip()
        label = "P{0} {1}".format(row.get('participant', '?'), row.get('loading_condition', '?'))
        if not odb_path:
            results.append((label, 'SKIPPED', 'no odb_path'))
            continue
        inp_path = os.path.splitext(odb_path)[0] + '.inp'
        if not os.path.isfile(inp_path):
            results.append((label, 'SKIPPED', 'inp not found: {0}'.format(inp_path)))
            continue
        try:
            out_path = remove_preload(inp_path)
            results.append((label, 'OK', out_path))
        except Exception as e:
            results.append((label, 'FAILED', str(e)))
        print("")

    print("=== Summary ===")
    for label, status, detail in results:
        print("{0:20s} | {1:8s} | {2}".format(label, status, detail))


# ============================================================
if 'INP_PATH' in dir() and INP_PATH:
    remove_preload(INP_PATH)
else:
    run_batch()
