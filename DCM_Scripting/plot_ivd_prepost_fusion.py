"""
plot_ivd_prepost_fusion.py - Compare intervertebral disc (IVD) strain between
PreOp-NoPreload and PostOp, to evaluate the effect of fusion surgery on the
IVD (adjacent-segment strain is a known fusion concern - does fusing levels
push more strain into the IVD, vs. decompression-only?).

Metric: % IVD volume above threshold, using a FIXED, COHORT-POOLED threshold
per percentile - T90/T95/T97/T99 are each computed ONCE from every patient's,
every condition's, every state's element data pooled together
(mps_common.volume_weighted_percentile applied to the whole pooled cache),
giving one shared MPS cutoff per percentile. Each patient/condition/state's
value is then % of ITS OWN volume exceeding that SAME shared cutoff
(mps_common.pct_volume_above) - same convention as the original Cord analysis
in "Alex_results_multipatient_plot - compare_threshold.py"'s first plot
(T90/T95/T99 pooled across the whole cohort, then % volume above per patient).

This replaces an earlier version of this script that computed each
percentile independently PER job (per patient/condition/state) rather than
pooling first - that gave each state its own, different MPS cutoff by
construction, which couldn't answer "did more volume cross a shared
threshold after surgery" (the cutoff itself moved between states, not just
the volume above it). Pooling first, then measuring against one fixed value,
is what makes PreOp vs PostOp comparable at all.

A patient-specific-but-state-pooled threshold (pool only that one patient's
PreOp+PostOp data, not the whole cohort) is a plausible alternative if
individual baseline IVD strain varies enough across patients to wash out
patient-specific signal under one cohort-wide cutoff - not implemented here,
flagged for later if the cohort-pooled version doesn't show anything useful.

Color = percentile, marker shape = loading condition (Flexion square,
Extension triangle), fill = state (PreOp-NoPreload solid, PostOp hollow), x
ordered by change in mJOA (postop - preop, ascending), fusion patients (P1,
P5, P7, P8) shaded orange across the whole column - same convention as
plot_prepost_sortedbymjoachange.py's combined single-row plot. The raw peak
is still saved in the summary CSV for reference, just not plotted.

Reads '_ivd_mps.csv' (from Alex_results_extraction_IVD.py), derived from
id_map.csv's existing 'csv_path' column by suffix swap (csv_path.replace(
'_mps.csv', '_ivd_mps.csv')) - same convention '_topology.csv' already uses
elsewhere in this codebase. As of this script's creation, IVD extraction is
still in progress (4 PostOp jobs - N01-011 and N01-017, both conditions -
pending an Abaqus-version fix), so this script will report those as
missing/skipped rather than failing outright; re-run once extraction for them
completes.

Caches the reduced per-element (mps, volume) data (element-wise peak over all
frames) to cache_ivd_prepost_peak.csv, so a later %-above-threshold plot can
reuse it without re-reading the raw '_ivd_mps.csv' files.

Run: python plot_ivd_prepost_fusion.py
"""

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from mps_common import volume_weighted_percentile, pct_volume_above, PLOT_STYLE

# ============================================================
# USER SETTINGS
# ============================================================
ID_MAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'id_map.csv')
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

STATE_PREOP = 'PreOp-NoPreload'
STATE_POSTOP = 'PostOp'
STATE_FILLED = {STATE_PREOP: True, STATE_POSTOP: False}   # solid vs hollow marker

CONDITION_MARKERS = {'flexion': 's', 'extension': '^'}

PERCENTILES = {'t90': 0.90, 't95': 0.95, 't97': 0.97, 't99': 0.99}
PERCENTILE_COLORS = {'t90': '#548235', 't95': '#2e75b6', 't97': '#c00000', 't99': '#7030a0'}

# Fusion patients - shaded orange across the whole column (both states), same
# convention as plot_prepost_sortedbymjoachange.py's combined single-row plot.
FUSION_PARTICIPANTS = {'P1', 'P5', 'P7', 'P8'}
# ============================================================

plt.rcParams.update(PLOT_STYLE)


def reduce_to_peak(df):
    peak_mps = df.groupby('element_label')['mps'].max()
    volume = df.groupby('element_label')['volume'].first()
    return pd.DataFrame({'mps': peak_mps, 'volume': volume}).reset_index()


# ============================================================
# Cache: reduced (participant, loading_condition, state) -> [element_label, mps, volume]
# ============================================================
own_cache_path = os.path.join(OUT_DIR, 'cache_ivd_prepost_peak.csv')

id_map = pd.read_csv(ID_MAP_PATH, skipinitialspace=True)

if os.path.isfile(own_cache_path):
    print("Loading cached reduced data: {}".format(own_cache_path))
    cache_df = pd.read_csv(own_cache_path)
else:
    print("No cache at {} yet - building it.".format(own_cache_path))
    rows = id_map[
        (id_map['State'].astype(str).str.strip().isin([STATE_PREOP, STATE_POSTOP])) &
        (id_map['loading_condition'].astype(str).str.strip().str.lower().isin(['flexion', 'extension']))
    ]

    cache_parts = []
    missing = []
    for _, row in rows.iterrows():
        participant = int(row['participant'])
        condition = str(row['loading_condition']).strip()
        state = str(row['State']).strip()
        csv_path = str(row.get('csv_path', '')).strip()
        if not csv_path or csv_path.lower() == 'nan':
            missing.append((participant, condition, state))
            continue
        ivd_csv_path = csv_path.replace('_mps.csv', '_ivd_mps.csv')
        if not os.path.isfile(ivd_csv_path):
            missing.append((participant, condition, state))
            continue
        raw = pd.read_csv(ivd_csv_path)
        # Completeness check: Alex_results_extraction_IVD.py opens the file in
        # 'w' mode and writes rows in a loop, so a file read WHILE extraction
        # is still running is a clean-but-truncated read (fewer complete rows,
        # not a parse error) - pandas wouldn't raise, it would just silently
        # understate peak/percentile strain. id_map.csv's 'last_frame_idx'
        # column (from the original Cord extraction, same ODB/step, so same
        # frame count) lets us catch this: if the IVD file's last frame_index
        # doesn't reach it, extraction for this job isn't finished yet.
        expected_last_frame = row.get('last_frame_idx', None)
        if expected_last_frame not in (None, '') and not pd.isna(expected_last_frame):
            if raw['frame_index'].max() < int(expected_last_frame):
                print("  P{} ({}, {}): '_ivd_mps.csv' exists but looks incomplete "
                      "(max frame_index {} < expected {}) - still extracting, skipping for now.".format(
                          participant, condition, state, raw['frame_index'].max(), int(expected_last_frame)))
                missing.append((participant, condition, state))
                continue
        reduced = reduce_to_peak(raw)
        reduced.insert(0, 'state', state)
        reduced.insert(0, 'loading_condition', condition)
        reduced.insert(0, 'participant', participant)
        cache_parts.append(reduced)

    if missing:
        print("Skipping {} job(s) missing '_ivd_mps.csv' (IVD extraction not done yet for these):".format(
            len(missing)))
        for p, c, s in missing:
            print("  P{} ({}, {})".format(p, c, s))

    if not cache_parts:
        raise SystemExit("No IVD data loaded - check that Alex_results_extraction_IVD.py has been run.")

    cache_df = pd.concat(cache_parts, ignore_index=True)
    cache_df.to_csv(own_cache_path, index=False)
    print("Cached reduced data: {}".format(own_cache_path))

# ============================================================
# mJOA change ordering - same convention as plot_prepost_sortedbymjoachange.py
# ============================================================
preop_mjoa_by_participant = {}
for _, row in id_map.iterrows():
    if str(row.get('State', '')).strip().lower() != 'preop':
        continue
    p_label = 'P{}'.format(int(row['participant']))
    preop_mjoa_by_participant[p_label] = row.get('mJOA', '')

postop_mjoa_by_participant = {}
for _, row in id_map.iterrows():
    if str(row.get('State', '')).strip().lower() != 'postop':
        continue
    p_label = 'P{}'.format(int(row['participant']))
    postop_mjoa_by_participant[p_label] = row.get('mJOA', '')

mjoa_delta_by_participant = {}
for p_label, preop_val in preop_mjoa_by_participant.items():
    postop_val = postop_mjoa_by_participant.get(p_label)
    if preop_val in ('', None) or postop_val in ('', None):
        continue
    try:
        mjoa_delta_by_participant[p_label] = float(postop_val) - float(preop_val)
    except (TypeError, ValueError):
        continue

# ============================================================
# Cohort-pooled thresholds - computed ONCE from every patient/condition/state's
# element data pooled together, not per job. This is the fixed MPS cutoff each
# job's volume is measured against below.
# ============================================================
COHORT_THRESHOLDS = {name: volume_weighted_percentile(cache_df, p=p) for name, p in PERCENTILES.items()}
print("Cohort-pooled thresholds: " + "  ".join(
    "{}={:.4f}".format(name.upper(), val) for name, val in COHORT_THRESHOLDS.items()))

# ============================================================
# Summary: % IVD volume above each cohort-pooled threshold, per
# participant/condition/state - plus the raw peak, kept for reference only
# (one row per participant/condition/state, not per percentile, since it
# doesn't depend on p).
# ============================================================
records = []
peak_records = []
for (participant, condition, state), grp in cache_df.groupby(['participant', 'loading_condition', 'state']):
    p_label = 'P{}'.format(participant)
    for name, threshold_val in COHORT_THRESHOLDS.items():
        records.append({
            'participant': p_label,
            'loading_condition': condition,
            'state': state,
            'percentile': name,
            'pct_above': pct_volume_above(grp, threshold_val),
        })
    peak_records.append({
        'participant': p_label, 'loading_condition': condition, 'state': state,
        'peak_ivd_mps': grp['mps'].max(),
    })
summary = pd.DataFrame(records)
peak_summary = pd.DataFrame(peak_records)

# Only patients with both a pre-op and post-op mJOA value (needed for the x-axis ordering).
summary = summary[summary['participant'].isin(mjoa_delta_by_participant)]
peak_summary = peak_summary[peak_summary['participant'].isin(mjoa_delta_by_participant)]

summary_path = os.path.join(OUT_DIR, 'multipatient_ivd_summary_percentiles_prepost_sortedbymjoachange.csv')
summary.merge(peak_summary, on=['participant', 'loading_condition', 'state']).to_csv(summary_path, index=False)
print(summary.to_string(index=False))

participants = sorted(summary['participant'].unique(),
                       key=lambda p: (mjoa_delta_by_participant[p], int(p[1:])))
x_pos = {p: i for i, p in enumerate(participants)}

# ============================================================
# Plot: one SEPARATE figure per percentile (not one combined multi-color
# plot - with 4 percentiles x 2 states x 2 conditions all overlaid, individual
# patients like P1 were hard to pick out). Each figure: marker shape = loading
# condition, fill = state (PreOp-NoPreload solid, PostOp hollow), fusion
# patients shaded orange across the whole column, x ordered by mJOA change.
# Same per-threshold-figure pattern already used elsewhere in this codebase
# (see "Alex_results_multipatient_plot - compare_threshold.py"'s per-patient
# blob-distribution plots).
# ============================================================
from matplotlib.patches import Patch

state_handles = [
    Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='black', label='PreOp (no preload)'),
    Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='none', label='PostOp'),
]
condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                      for cond, marker in CONDITION_MARKERS.items()]
fusion_handle = [Patch(facecolor='orange', alpha=0.2, label='Fusion')]
blank = Line2D([0], [0], linestyle='none', marker='None', label='')

all_handles = (
    [Line2D([0], [0], linestyle='none', marker='None', label='State')] + state_handles +
    [blank] +
    [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + condition_handles +
    [blank] + fusion_handle
)

plot_paths = []
for name, p in PERCENTILES.items():
    color = PERCENTILE_COLORS[name]
    col_df = summary[summary['percentile'] == name]

    fig, ax = plt.subplots(figsize=(9, 5.5))

    for participant, grp in col_df.groupby('participant'):
        for state, state_grp in grp.groupby('state'):
            if len(state_grp) == 2:
                xp = x_pos[participant]
                ax.vlines(xp, state_grp['pct_above'].min(), state_grp['pct_above'].max(),
                          color=color, linewidth=1.0, alpha=0.5, zorder=2)
    for (condition, state), grp in col_df.groupby(['loading_condition', 'state']):
        marker = CONDITION_MARKERS.get(condition.strip().lower(), 'o')
        filled = STATE_FILLED.get(state, True)
        xs = [x_pos[p] for p in grp['participant']]
        if filled:
            ax.scatter(xs, grp['pct_above'], color=color, marker=marker, s=60, zorder=3)
        else:
            ax.scatter(xs, grp['pct_above'], facecolors='none', edgecolors=color, marker=marker,
                       s=60, linewidths=1.4, zorder=3)

    for fusion_p in FUSION_PARTICIPANTS:
        if fusion_p in x_pos:
            xp = x_pos[fusion_p]
            ax.axvspan(xp - 0.5, xp + 0.5, color='orange', alpha=0.2, zorder=0)

    ax.set_xticks(range(len(participants)))
    ax.set_xticklabels(['{}\n({:+.0f})'.format(p, mjoa_delta_by_participant[p]) for p in participants])
    ax.annotate('Δ mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
    ax.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')
    ax.set_ylabel('% IVD volume above threshold')
    ax.set_title('IVD strain ({} = {:.4f}, cohort-pooled): PreOp (no preload) vs PostOp, by fusion status'.format(
        name.upper(), COHORT_THRESHOLDS[name]))
    ax.legend(handles=all_handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

    fig.tight_layout()

    plot_path = os.path.join(OUT_DIR, 'multipatient_ivd_plot_{}_prepost_sortedbymjoachange.pdf'.format(name))
    fig.savefig(plot_path, bbox_inches='tight')
    plt.close(fig)
    plot_paths.append(plot_path)

print()
print("Summary saved: {}".format(summary_path))
for p in plot_paths:
    print("Plot saved: {}".format(p))

# ============================================================
# Patient-wise version: instead of one cohort-pooled threshold shared by
# everyone, each patient gets their OWN threshold - pooling that one
# patient's PreOp-NoPreload + PostOp data (both conditions) together. This is
# the "patient-specific-but-state-pooled" alternative flagged in the
# docstring above: still comparable PreOp vs PostOp WITHIN a patient (same
# fixed cutoff both times), but no longer forces every patient to be judged
# against one shared cohort-wide cutoff - useful if baseline IVD strain
# varies enough across patients that a cohort threshold would sit miles
# above some patients' whole distribution and miles below others'.
#
# Separate output files (same plots, '_patientwise' suffix) - not a
# replacement for the cohort-pooled plots above, a second thing to compare
# against them. Reuses cache_df, mjoa_delta_by_participant,
# PERCENTILES/PERCENTILE_COLORS, CONDITION_MARKERS, STATE_FILLED,
# FUSION_PARTICIPANTS, all_handles, OUT_DIR, Line2D, pd, plt, os already
# loaded/defined above - does not modify anything above this point.
# ============================================================
PATIENT_THRESHOLDS = {}
for participant, pgrp in cache_df.groupby('participant'):
    p_label = 'P{}'.format(participant)
    PATIENT_THRESHOLDS[p_label] = {name: volume_weighted_percentile(pgrp, p=p) for name, p in PERCENTILES.items()}

pw_records = []
pw_peak_records = []
for (participant, condition, state), grp in cache_df.groupby(['participant', 'loading_condition', 'state']):
    p_label = 'P{}'.format(participant)
    for name, threshold_val in PATIENT_THRESHOLDS[p_label].items():
        pw_records.append({
            'participant': p_label,
            'loading_condition': condition,
            'state': state,
            'percentile': name,
            'threshold_value': threshold_val,
            'pct_above': pct_volume_above(grp, threshold_val),
        })
    pw_peak_records.append({
        'participant': p_label, 'loading_condition': condition, 'state': state,
        'peak_ivd_mps': grp['mps'].max(),
    })
pw_summary = pd.DataFrame(pw_records)
pw_peak_summary = pd.DataFrame(pw_peak_records)

pw_summary = pw_summary[pw_summary['participant'].isin(mjoa_delta_by_participant)]
pw_peak_summary = pw_peak_summary[pw_peak_summary['participant'].isin(mjoa_delta_by_participant)]

pw_summary_path = os.path.join(
    OUT_DIR, 'multipatient_ivd_summary_percentiles_prepost_sortedbymjoachange_patientwise.csv')
pw_summary.merge(pw_peak_summary, on=['participant', 'loading_condition', 'state']).to_csv(
    pw_summary_path, index=False)
print()
print(pw_summary.to_string(index=False))

pw_plot_paths = []
for name, p in PERCENTILES.items():
    color = PERCENTILE_COLORS[name]
    col_df = pw_summary[pw_summary['percentile'] == name]

    fig, ax = plt.subplots(figsize=(9, 5.5))

    for participant, grp in col_df.groupby('participant'):
        for state, state_grp in grp.groupby('state'):
            if len(state_grp) == 2:
                xp = x_pos[participant]
                ax.vlines(xp, state_grp['pct_above'].min(), state_grp['pct_above'].max(),
                          color=color, linewidth=1.0, alpha=0.5, zorder=2)
    for (condition, state), grp in col_df.groupby(['loading_condition', 'state']):
        marker = CONDITION_MARKERS.get(condition.strip().lower(), 'o')
        filled = STATE_FILLED.get(state, True)
        xs = [x_pos[p] for p in grp['participant']]
        if filled:
            ax.scatter(xs, grp['pct_above'], color=color, marker=marker, s=60, zorder=3)
        else:
            ax.scatter(xs, grp['pct_above'], facecolors='none', edgecolors=color, marker=marker,
                       s=60, linewidths=1.4, zorder=3)

    for fusion_p in FUSION_PARTICIPANTS:
        if fusion_p in x_pos:
            xp = x_pos[fusion_p]
            ax.axvspan(xp - 0.5, xp + 0.5, color='orange', alpha=0.2, zorder=0)

    ax.set_xticks(range(len(participants)))
    ax.set_xticklabels(['{}\n({:+.0f})'.format(p, mjoa_delta_by_participant[p]) for p in participants])
    ax.annotate('Δ mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
    ax.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')
    ax.set_ylabel('% IVD volume above threshold')
    ax.set_title('IVD strain ({}, patient-specific threshold): PreOp (no preload) vs PostOp, '
                  'by fusion status'.format(name.upper()))
    ax.legend(handles=all_handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

    fig.tight_layout()

    pw_plot_path = os.path.join(
        OUT_DIR, 'multipatient_ivd_plot_{}_prepost_sortedbymjoachange_patientwise.pdf'.format(name))
    fig.savefig(pw_plot_path, bbox_inches='tight')
    plt.close(fig)
    pw_plot_paths.append(pw_plot_path)

print()
print("Summary saved: {}".format(pw_summary_path))
for p in pw_plot_paths:
    print("Plot saved: {}".format(p))

# ============================================================
# Fusion-patient PreOp vs PostOp diff, T95 ONLY, BOTH threshold versions
# (cohort-pooled and patient-wise) side by side in one Obsidian-ready
# markdown table. Flexion and Extension kept FULLY SEPARATE (never averaged
# together - they are distinct loading modes; averaging them hides cases
# where one condition increases PostOp while the other decreases, as
# happened for P5 Flexion at patient-wise T90). T95 is the one combination
# (of cohort-pooled/patient-wise x T90/T95/T97/T99) where all of P5/P7/P8
# increase PostOp in BOTH conditions individually under the patient-wise
# version, with P1 the clear exception - if that changes (e.g. after further
# data fixes), update PRINT_PERCENTILE below rather than printing every
# combination again. Reuses pw_summary, summary, FUSION_PARTICIPANTS, pd
# already loaded/defined above - does not modify anything above this point.
# ============================================================
PRINT_PERCENTILE = 't95'


def _fusion_pivot(df, percentile):
    fusion_df = df[(df['participant'].isin(FUSION_PARTICIPANTS)) & (df['percentile'] == percentile)]
    piv = fusion_df.pivot_table(index=['participant', 'loading_condition'],
                                 columns='state', values='pct_above').reset_index()
    piv['delta'] = piv[STATE_POSTOP] - piv[STATE_PREOP]
    return piv.set_index(['participant', 'loading_condition'])


pw_piv = _fusion_pivot(pw_summary, PRINT_PERCENTILE)
glob_piv = _fusion_pivot(summary, PRINT_PERCENTILE)
combined = pw_piv.join(glob_piv, lsuffix='_pw', rsuffix='_global').reset_index()
combined = combined.sort_values(['participant', 'loading_condition'])

print()
print("Fusion-patient PreOp vs PostOp diff - {} (Obsidian-ready markdown):".format(PRINT_PERCENTILE.upper()))
print()
print("| Participant | Condition | PreOp PW (%) | PostOp PW (%) | Delta PW (pp) "
      "| PreOp Global (%) | PostOp Global (%) | Delta Global (pp) |")
print("|---|---|---|---|---|---|---|---|")
for _, r in combined.iterrows():
    print("| {} | {} | {:.2f} | {:.2f} | {:+.2f} | {:.2f} | {:.2f} | {:+.2f} |".format(
        r['participant'], r['loading_condition'],
        r[STATE_PREOP + '_pw'], r[STATE_POSTOP + '_pw'], r['delta_pw'],
        r[STATE_PREOP + '_global'], r[STATE_POSTOP + '_global'], r['delta_global']))
