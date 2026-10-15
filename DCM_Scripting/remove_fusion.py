"""
remove_fusion.py - Copy a PostOp job's .inp and remove its fusion, i.e. the C2-C7
rigid body constraint: the '** Constraint: ...' comment plus the '*Rigid Body, ref
node=..., tie nset=...' keyword line (and any data lines under it).

Everything else is left untouched - the reference node, the nsets it used and the
output requests all stay in the file, unreferenced. Diffed against a matching
no-fusion job, the model is then functionally identical (same C2-top / C7-base
couplings, encastre BC and steps).

Never modifies the source .inp (only ever opened read-only) - writes into a new,
renumbered sibling job folder instead. The job number's leading digit is replaced
with '5', any 'rigidbody...' token is dropped from the name and '-nofusion' is
appended:
  .../Job-004-N01-011-PostOp-BC0pt35-nopreload-rigidbodyC2C7/Job-004-...-rigidbodyC2C7.inp
  -> .../Job-504-N01-011-PostOp-BC0pt35-nopreload-nofusion/Job-504-...-nofusion.inp
Refuses to overwrite the source file; re-running overwrites its own previous output.

Pure text processing, no Abaqus API. Set INP_PATH below to a path (or list of paths)
for single-file mode, or to None to run the id_map.csv batch: every row where
surgery_modelled == 'fusion', State == 'PostOp' and loading_condition in
LOADING_CONDITIONS (Flexion and Extension; Oscillation rows are skipped). A runner
.bat for the generated jobs is written to BAT_PATH. Run via `abaqus python remove_fusion.py`,
execfile() in the CAE kernel, or plain Python. Paths may also be passed on the
command line.
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
# Set to a specific .inp path (or list) for single-file mode, or None to run the
# id_map.csv batch.
INP_PATH = None  # r'D:\path\to\Job-004-...-rigidbodyC2C7.inp'

NEW_JOB_DIGIT = '5'
NAME_SUFFIX = '-nofusion'

# Runner script that submits every generated job. Job folders are written relative to
# the .bat's own folder, so it should live in the folder that contains the patient
# folders. Set to None to skip writing it.
BAT_PATH = r'D:\Charlotte\ABAQUS\multijob_nofusion.bat'

# Which id_map.csv loading conditions to process in batch mode (4 patients - P1, P5,
# P7, P8 - each with a Flexion and an Extension job).
LOADING_CONDITIONS = ('flexion', 'extension')
BAT_LOG_NAME = 'multijob_nofusion.log'
ABAQUS_CMD = 'abq2025'
CPUS = 4


# ============================================================
# HELPERS
# ============================================================
def _is_keyword(line):
    return line.startswith('*') and not line.startswith('**')


def _find_rigid_body_spans(lines):
    """Return [(start, end), ...] (end exclusive) for each '*Rigid Body' block,
    including a '** Constraint:' comment line immediately before it."""
    spans = []
    for i, line in enumerate(lines):
        if not re.match(r'\*rigid body\b', line.strip(), re.IGNORECASE):
            continue
        end = i + 1
        while end < len(lines) and not lines[end].startswith('*'):
            end += 1  # data lines (none expected, but harmless)
        start = i
        if i > 0 and re.match(r'\*\*\s*Constraint:', lines[i - 1].strip(), re.IGNORECASE):
            start = i - 1
        spans.append((start, end))
    return spans


def _new_name(name):
    m = re.match(r'(Job-)(\d+)(.*)', name)
    if not m:
        raise RuntimeError("Could not find a 'Job-<number>' prefix in '{0}'".format(name))
    prefix, number, rest = m.groups()
    number = NEW_JOB_DIGIT + number[1:]
    rest = re.sub(r'[-_]?rigidbody\w*', '', rest, flags=re.IGNORECASE)
    return prefix + number + rest + NAME_SUFFIX


def _output_path(inp_path):
    orig_dir = os.path.dirname(os.path.abspath(inp_path))
    parent_dir = os.path.dirname(orig_dir)
    file_base, file_ext = os.path.splitext(os.path.basename(inp_path))
    new_base = _new_name(file_base)
    return os.path.join(parent_dir, _new_name(os.path.basename(orig_dir)), new_base + file_ext), new_base


# ============================================================
# CORE TRANSFORM
# ============================================================
def remove_fusion(inp_path):
    with open(inp_path, 'r') as f:
        lines = f.readlines()

    print("Job: {0}".format(inp_path))
    spans = _find_rigid_body_spans(lines)
    if not spans:
        raise RuntimeError("No '*Rigid Body' constraint found - nothing to remove.")

    drop = set()
    for a, b in spans:
        print("  Removing lines {0}-{1}: {2}".format(a + 1, b, lines[b - 1].strip() if b - a == 1
                                                      else ' | '.join(l.strip() for l in lines[a:b])))
        drop.update(range(a, b))
    out_lines = [l for i, l in enumerate(lines) if i not in drop]

    out_path, new_base = _output_path(inp_path)
    if os.path.abspath(out_path) == os.path.abspath(inp_path):
        raise RuntimeError("Refusing to write: output path resolved to the source file itself ({0})".format(inp_path))

    # keep the '** Job name:' header comment consistent with the new file name
    old_base = os.path.splitext(os.path.basename(inp_path))[0]
    for i, l in enumerate(out_lines[:60]):
        if l.startswith('** Job name:') and old_base in l:
            out_lines[i] = l.replace(old_base, new_base, 1)
            break

    out_dir = os.path.dirname(out_path)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    with open(out_path, 'w') as f:
        f.writelines(out_lines)
    print("  Wrote (new folder, source untouched): {0}".format(out_path))
    return out_path


# ============================================================
# BATCH MODE (id_map.csv, fusion-modelled PostOp Flexion/Extension rows)
# ============================================================
def _read_id_map():
    if sys.version_info[0] >= 3:
        read_mode, kwargs = 'r', {'newline': ''}
    else:
        read_mode, kwargs = 'rb', {}
    with open(ID_MAP_PATH, read_mode, **kwargs) as f:
        return list(csv.DictReader(f))


def run_batch():
    rows = [r for r in _read_id_map()
            if r.get('surgery_modelled', '').strip().lower() == 'fusion'
            and r.get('State', '').strip().lower() == 'postop'
            and r.get('loading_condition', '').strip().lower() in LOADING_CONDITIONS]
    print("Found {0} fusion PostOp row(s) in id_map.csv".format(len(rows)))

    results = []
    jobs = []
    for row in rows:
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
            out_path = remove_fusion(inp_path)
            results.append((label, 'OK', out_path))
            jobs.append((row.get('id', '').strip(), label, out_path))
        except Exception as e:
            results.append((label, 'FAILED', str(e)))
        print("")

    print("=== Summary ===")
    for label, status, detail in results:
        print("{0:20s} | {1:8s} | {2}".format(label, status, detail))
    return jobs


# ============================================================
# .bat RUNNER
# ============================================================
def _rel_dir(path, start):
    """Path of `path` relative to `start` (for use after %~dp0); absolute if on another drive."""
    try:
        rel = os.path.relpath(path, start)
    except ValueError:
        return None
    return rel


def write_bat(jobs, bat_path):
    """Write a .bat that runs each (patient_id, label, inp_path) job in turn and logs
    pass/fail from the job's own .sta file ('COMPLETED SUCCESSFULLY')."""
    if not jobs:
        print("No jobs - not writing .bat")
        return
    bat_dir = os.path.dirname(os.path.abspath(bat_path))

    patients = []
    for pid, _, _ in jobs:
        if pid not in patients:
            patients.append(pid)

    out = [
        '@echo off',
        'REM ============================================================',
        'REM Runs the fusion-removed (no-fusion) jobs generated by',
        'REM DCM_Scripting\\remove_fusion.py for the fusion PostOp rows in id_map.csv.',
        'REM',
        'REM Each DIR_N is relative to this script\'s own folder; each block',
        'REM returns to it afterward via "%~dp0".',
        'REM',
        'REM Pass/fail for each job is logged to {0} (next to this'.format(BAT_LOG_NAME),
        'REM script), judged by whether the job\'s own .sta file contains',
        'REM "COMPLETED SUCCESSFULLY" - the abaqus exit code isn\'t reliable.',
        'REM',
        'REM ask_delete=OFF - overwrite old job files from a previous run of the',
        'REM same job name instead of prompting (which would hang the script).',
        'REM ============================================================',
        '',
        'set CPUS={0}'.format(CPUS),
        'set LOGFILE=%~dp0{0}'.format(BAT_LOG_NAME),
        'echo No-fusion batch run started: %DATE% %TIME% > "%LOGFILE%"',
        '',
    ]
    for n, pid in enumerate(patients, 1):
        out.append('set PATIENT_{0}={1}'.format(n, pid))
    out.append('')

    for k, (pid, label, inp_path) in enumerate(jobs, 1):
        job_dir = os.path.dirname(os.path.abspath(inp_path))
        rel = _rel_dir(job_dir, bat_dir)
        dir_value = rel if rel is not None else job_dir
        job_name = os.path.splitext(os.path.basename(inp_path))[0]
        out += ['REM --- {0} {1} ---'.format(label, pid),
                'set DIR_{0}={1}'.format(k, dir_value),
                'set JOBNAME_{0}={1}'.format(k, job_name),
                '']

    out.append('REM ============================================================')
    out.append('')
    for k, (pid, _, _) in enumerate(jobs, 1):
        p = patients.index(pid) + 1
        out += ['cd /d "%~dp0%DIR_{0}%"'.format(k),
                'call {0} job=%JOBNAME_{1}% input=%JOBNAME_{1}%.inp interactive cpus=%CPUS% ask_delete=OFF'.format(ABAQUS_CMD, k),
                'call :LogResult %JOBNAME_{0}% %PATIENT_{1}%'.format(k, p),
                'cd /d "%~dp0"',
                '']
    out += ['goto :eof',
            '',
            ':LogResult',
            'findstr /C:"COMPLETED SUCCESSFULLY" "%~1.sta" >nul',
            'if %ERRORLEVEL% EQU 0 (',
            '    echo %~1 finished successfully >> "%LOGFILE%"',
            ') else (',
            '    echo %~2 job %~1 failed >> "%LOGFILE%"',
            ')',
            'goto :eof',
            '']

    bat_parent = os.path.dirname(bat_path)
    if bat_parent and not os.path.isdir(bat_parent):
        os.makedirs(bat_parent)
    with open(bat_path, 'w') as f:
        f.write('\n'.join(out))  # text mode writes CRLF on Windows
    print("Wrote runner: {0} ({1} job(s))".format(bat_path, len(jobs)))


def _patient_id_from_path(p):
    m = re.search(r'N\d\d-\d{3}', p)
    return m.group(0) if m else 'unknown'


# ============================================================
if __name__ == '__main__' and len(sys.argv) > 1:
    INP_PATH = sys.argv[1:]

if INP_PATH:
    jobs = []
    for p in ([INP_PATH] if isinstance(INP_PATH, str) else INP_PATH):
        out_path = remove_fusion(p)
        jobs.append((_patient_id_from_path(p), os.path.basename(p), out_path))
else:
    jobs = run_batch()

if BAT_PATH:
    write_bat(jobs, BAT_PATH)
