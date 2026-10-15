"""
make_extension_jobs.py - Create the Extension version of each Flexion '-wfusion' job.

For every Flexion .inp in FLEXION_INPS this writes a sibling job folder next to the
source folder, with:
  - job number + 1             Job-502-...  ->  Job-503-...
  - 'BC0pt35' -> 'BCminus0pt175' in the folder name, file name and '** Job name:' header
    (if the name has no BC token, '-BCminus0pt175' is inserted before '-wfusion')
  - the BC-2 rotation (UR1, dof 4) changed to -0.175:
        ** Name: BC-2 Type: Displacement/Rotation
        *Boundary, amplitude=Amp-2
        m_Set-4, 4, 4, 0.35        ->        m_Set-4, 4, 4, -0.175
Everything else is copied unchanged. The one exception to "source is read-only": if a
Flexion source's BC-2 isn't 0.35 (FIX_FLEXION = True), that single line is corrected in
the source, with the original kept as <name>.inp.orig.

Existing target files are not overwritten unless OVERWRITE = True.

Pure text processing - run with `abaqus python make_extension_jobs.py`, execfile()
in the CAE kernel, or plain Python.
"""
from __future__ import print_function

import os
import re

# ============================================================
# USER SETTINGS
# ============================================================
ROOT = r'D:\Charlotte\ABAQUS'

# Flexion -wfusion jobs, relative to ROOT.
FLEXION_INPS = [
    r'N01-016\Job-502-N01-016-Post-Op-BC0pt35-wfusion\Job-502-N01-016-Post-Op-BC0pt35-wfusion.inp',
    r'N01-017\Job-511-N01-017-PostOp-BC0pt35-wfusion\Job-511-N01-017-PostOp-BC0pt35-wfusion.inp',
    r'N01-012\Job-507-N01-012-PostOp-BC0pt35-wfusion\Job-507-N01-012-PostOp-BC0pt35-wfusion.inp',
    r'N31-038\Job-503-N31-038-PostOpv2-BC0pt35-wfusion\Job-503-N31-038-PostOpv2-BC0pt35-wfusion.inp',
    r'N31-048\Job-501-N31-048-PostOp-wfusion\Job-501-N31-048-PostOp-wfusion.inp',
    r'N01-021\Job-503-N01-021-PostOp-BC0pt35-wfusion\Job-503-N01-021-PostOp-BC0pt35-wfusion.inp',
    r'N31-053\Job-501-N31-053-PostOp-BC0pt35-wfusion\Job-501-N31-053-PostOp-BC0pt35-wfusion.inp',
    r'N01-018\Job-503-N01-018-PostOp-BC0pt35-wfusion\Job-503-N01-018-PostOp-BC0pt35-wfusion.inp',
]

FLEXION_VALUE = '0.35'
EXTENSION_VALUE = '-0.175'
OVERWRITE = False

# If a Flexion source's BC-2 isn't FLEXION_VALUE, fix that one line in the source
# (original saved alongside as <name>.inp.orig). Set False to only warn.
FIX_FLEXION = True

# Set False to skip creating/fixing .inp files (e.g. to only regenerate the .bat).
CREATE_JOBS = True

# Runner script for all the jobs (each patient's Flexion then Extension job), in the
# same form as multijob_nofusion.bat from remove_fusion.py. None to skip.
BAT_PATH = r'D:\Charlotte\ABAQUS\multijob_lamiwfusion.bat'
BAT_LOG_NAME = 'multijob_lamiwfusion.log'
ABAQUS_CMD = 'abq2025'
CPUS = 4


# ============================================================
# HELPERS
# ============================================================
def _ext_name(name):
    """Job-502-...-BC0pt35-wfusion -> Job-503-...-BCminus0pt175-wfusion"""
    m = re.match(r'(Job-)(\d+)(.*)', name)
    if not m:
        raise RuntimeError("No 'Job-<number>' prefix in '{0}'".format(name))
    prefix, number, rest = m.groups()
    number = str(int(number) + 1).zfill(len(number))
    if re.search(r'BC0pt35', rest):
        rest = rest.replace('BC0pt35', 'BCminus0pt175')
    elif re.search(r'-wfusion$', rest):
        rest = re.sub(r'-wfusion$', '-BCminus0pt175-wfusion', rest)
    else:
        raise RuntimeError("Don't know where to put the BC token in '{0}'".format(name))
    return prefix + number + rest


def _set_bc2(lines, value=None):
    """Change the BC-2 dof-4 value in place (to `value`, default EXTENSION_VALUE).
    Returns (line_idx, old_value)."""
    value = EXTENSION_VALUE if value is None else value
    for i, line in enumerate(lines):
        if not re.match(r'\*\*\s*Name:\s*BC-2\b', line.strip(), re.IGNORECASE):
            continue
        for j in range(i + 1, min(i + 8, len(lines))):
            m = re.match(r'^(\s*\S+,\s*4,\s*4,\s*)(\S+)(\s*)$', lines[j])
            if m:
                old = m.group(2)
                eol = '\n' if lines[j].endswith('\n') else ''
                lines[j] = m.group(1) + value + eol
                return j, old
        raise RuntimeError("Found '** Name: BC-2' but no '<set>, 4, 4, <value>' data line after it")
    raise RuntimeError("No '** Name: BC-2' block found")


def make_extension(rel_path):
    src = os.path.join(ROOT, rel_path)
    if not os.path.isfile(src):
        raise RuntimeError("source not found: {0}".format(src))

    src_dir = os.path.dirname(src)
    parent = os.path.dirname(src_dir)
    base = os.path.splitext(os.path.basename(src))[0]
    new_folder = _ext_name(os.path.basename(src_dir))
    new_base = _ext_name(base)
    out_path = os.path.join(parent, new_folder, new_base + '.inp')

    if os.path.abspath(out_path) == os.path.abspath(src):
        raise RuntimeError("target resolves to the source file")
    if os.path.exists(out_path) and not OVERWRITE:
        raise RuntimeError("target already exists (set OVERWRITE = True): {0}".format(out_path))

    with open(src, 'r') as f:
        lines = f.readlines()

    # The Flexion source must itself be at FLEXION_VALUE. If it isn't (e.g. it was copied
    # from an Extension model), correct that one line in the source, keeping a backup.
    flex_lines = list(lines)
    _, flex_old = _set_bc2(flex_lines, FLEXION_VALUE)
    if flex_old != FLEXION_VALUE:
        if FIX_FLEXION:
            backup = src + '.orig'
            if not os.path.exists(backup):
                with open(backup, 'w') as f:
                    f.writelines(lines)
            with open(src, 'w') as f:
                f.writelines(flex_lines)
            print("  Flexion source BC-2 was {0} -> set to {1} (backup: {2})".format(
                flex_old, FLEXION_VALUE, os.path.basename(backup)))
        else:
            print("  *** WARNING: Flexion source BC-2 is {0}, expected {1}".format(flex_old, FLEXION_VALUE))

    idx, old = _set_bc2(lines)
    warn = ''

    for i, l in enumerate(lines[:60]):
        if l.startswith('** Job name:') and base in l:
            lines[i] = l.replace(base, new_base, 1)
            break

    out_dir = os.path.dirname(out_path)
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    with open(out_path, 'w') as f:
        f.writelines(lines)
    print("  BC-2 line {0}: {1} -> {2}{3}".format(idx + 1, old, EXTENSION_VALUE, warn))
    return out_path


def write_bat(rel_inps, bat_path):
    """One DIR_n/JOBNAME_n pair per job (Flexion, then its Extension partner), a run
    block per job, and pass/fail logging from each job's .sta file."""
    bat_dir = os.path.dirname(os.path.abspath(bat_path))
    jobs = []  # (patient_id, label, abs_inp)
    for rel in rel_inps:
        src = os.path.join(ROOT, rel)
        d = os.path.dirname(src)
        base = os.path.splitext(os.path.basename(src))[0]
        pid = re.search(r'N\d\d-\d{3}', rel).group(0)
        ext = os.path.join(os.path.dirname(d), _ext_name(os.path.basename(d)), _ext_name(base) + '.inp')
        jobs.append((pid, 'Flexion', src))
        jobs.append((pid, 'Extension', ext))

    patients = []
    for pid, _, _ in jobs:
        if pid not in patients:
            patients.append(pid)

    out = [
        '@echo off',
        'REM ============================================================',
        'REM Runs the lami-with-fusion (counterfactual) jobs generated for the',
        'REM lami patients: each patient\'s Flexion job then its Extension job.',
        'REM Extension jobs come from DCM_Scripting\\make_extension_jobs.py.',
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
        'echo Lami-wfusion batch run started: %DATE% %TIME% > "%LOGFILE%"',
        '',
    ]
    for n, pid in enumerate(patients, 1):
        out.append('set PATIENT_{0}={1}'.format(n, pid))
    out.append('')
    for k, (pid, label, inp) in enumerate(jobs, 1):
        job_dir = os.path.dirname(os.path.abspath(inp))
        try:
            dir_value = os.path.relpath(job_dir, bat_dir)
        except ValueError:
            dir_value = job_dir
        out += ['REM --- {0} {1} ---'.format(pid, label),
                'set DIR_{0}={1}'.format(k, dir_value),
                'set JOBNAME_{0}={1}'.format(k, os.path.splitext(os.path.basename(inp))[0]),
                '']
    out += ['REM ============================================================', '']
    for k, (pid, _, _) in enumerate(jobs, 1):
        out += ['cd /d "%~dp0%DIR_{0}%"'.format(k),
                'call {0} job=%JOBNAME_{1}% input=%JOBNAME_{1}%.inp interactive cpus=%CPUS% ask_delete=OFF'.format(ABAQUS_CMD, k),
                'call :LogResult %JOBNAME_{0}% %PATIENT_{1}%'.format(k, patients.index(pid) + 1),
                'cd /d "%~dp0"',
                '']
    out += ['goto :eof', '', ':LogResult',
            'findstr /C:"COMPLETED SUCCESSFULLY" "%~1.sta" >nul',
            'if %ERRORLEVEL% EQU 0 (',
            '    echo %~1 finished successfully >> "%LOGFILE%"',
            ') else (',
            '    echo %~2 job %~1 failed >> "%LOGFILE%"',
            ')',
            'goto :eof', '']
    with open(bat_path, 'w') as f:
        f.write('\n'.join(out))  # text mode writes CRLF on Windows
    print("Wrote runner: {0} ({1} job(s))".format(bat_path, len(jobs)))


# ============================================================
results = []
for rel in (FLEXION_INPS if CREATE_JOBS else []):
    print(os.path.basename(rel))
    try:
        results.append((os.path.basename(rel), 'OK', make_extension(rel)))
    except Exception as e:
        results.append((os.path.basename(rel), 'FAILED', str(e)))
    print("")

print("=== Summary ===")
for name, status, detail in results:
    print("{0:8s} | {1}\n           {2}".format(status, name, detail))

if BAT_PATH:
    write_bat(FLEXION_INPS, BAT_PATH)
