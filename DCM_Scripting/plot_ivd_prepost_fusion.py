"""
plot_ivd_prepost_fusion.py - % IVD volume above T95, PreOp-NoPreload vs
PostOp, per patient/condition, fusion patients highlighted. Two threshold
versions, both plotted: global (cohort-pooled) and patientwise (per-patient,
2nd y-axis shows each patient's own threshold).

Reads '_ivd_mps.csv' via id_map.csv's csv_path. Caches peak-reduced data to
cache_ivd_prepost_peak.csv.

Outputs: 2 summary CSVs, 2 PDF plots, 1 markdown table (fusion patients) to stdout.

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

# T90/T97/T99 commented out, not deleted - everything downstream reads this dict.
PERCENTILES = {
    # 't90': 0.90,
    't95': 0.95,
    # 't97': 0.97,
    # 't99': 0.99,
}
PERCENTILE_COLORS = {
    # 't90': '#548235',
    't95': '#2e75b6',
    # 't97': '#c00000',
    # 't99': '#7030a0',
}

# Shaded orange across the whole column in every plot.
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
        # Catches a file still being written (truncated, not a parse error) -
        # compare against id_map.csv's 'last_frame_idx' from the Cord extraction.
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
# mJOA change ordering (x-axis)
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
# Cohort-pooled thresholds - one shared cutoff per percentile, pooled across all patients.
# ============================================================
COHORT_THRESHOLDS = {name: volume_weighted_percentile(cache_df, p=p) for name, p in PERCENTILES.items()}
print("Cohort-pooled thresholds: " + "  ".join(
    "{}={:.4f}".format(name.upper(), val) for name, val in COHORT_THRESHOLDS.items()))

# ============================================================
# Summary: % IVD volume above each cohort-pooled threshold, per job. Raw peak kept for reference.
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
# Plot: one figure per percentile. Marker shape = condition, fill = state.
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
# Patient-wise version: each patient's own threshold, pooled from their own
# PreOp-NoPreload + PostOp data. Separate output files ('_patientwise' suffix).
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

    # Starts at 0 explicitly - matplotlib's autoscale otherwise pads slightly below 0.
    ax.set_ylim(bottom=0)
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

    # --- Separate figure: same plot plus a second y-axis showing each
    # patient's own threshold value, saved to a distinct '_w_thresholds' file.
    fig2, ax2 = plt.subplots(figsize=(9, 5.5))

    for participant, grp in col_df.groupby('participant'):
        for state, state_grp in grp.groupby('state'):
            if len(state_grp) == 2:
                xp = x_pos[participant]
                ax2.vlines(xp, state_grp['pct_above'].min(), state_grp['pct_above'].max(),
                           color=color, linewidth=1.0, alpha=0.5, zorder=2)
    for (condition, state), grp in col_df.groupby(['loading_condition', 'state']):
        marker = CONDITION_MARKERS.get(condition.strip().lower(), 'o')
        filled = STATE_FILLED.get(state, True)
        xs = [x_pos[p] for p in grp['participant']]
        if filled:
            ax2.scatter(xs, grp['pct_above'], color=color, marker=marker, s=60, zorder=3)
        else:
            ax2.scatter(xs, grp['pct_above'], facecolors='none', edgecolors=color, marker=marker,
                        s=60, linewidths=1.4, zorder=3)
    for fusion_p in FUSION_PARTICIPANTS:
        if fusion_p in x_pos:
            xp = x_pos[fusion_p]
            ax2.axvspan(xp - 0.5, xp + 0.5, color='orange', alpha=0.2, zorder=0)
    ax2.set_ylim(bottom=0)

    # Right axis: each patient's own threshold value, flat per column, NaN
    # gaps so adjacent (unrelated) patients aren't connected. Pinned to 0 too
    # - both axes share the same plot box, so bottom=0 on both is what
    # actually aligns the two zero points.
    threshold_color = '#c00000'
    step_xs, step_ys = [], []
    for participant in participants:
        xp = x_pos[participant]
        val = PATIENT_THRESHOLDS[participant][name]
        step_xs.extend([xp - 0.5, xp + 0.5, float('nan')])
        step_ys.extend([val, val, float('nan')])
    ax2b = ax2.twinx()
    # Background element - behind the data markers (zorder 2-3), not competing with them.
    ax2b.plot(step_xs, step_ys, color=threshold_color, linewidth=1.2, alpha=0.4, zorder=1, solid_capstyle='butt')
    ax2b.set_ylim(bottom=0)
    ax2b.set_ylabel('Patient-specific {} threshold (MPS)'.format(name.upper()), color=threshold_color)
    ax2b.tick_params(axis='y', labelcolor=threshold_color)
    # re-show right spine for twin axis - thin/neutral, not a bold red bar;
    # the colored tick labels/ylabel do the association work.
    ax2b.spines['right'].set_visible(True)
    ax2b.spines['right'].set_linewidth(0.8)

    ax2.set_xticks(range(len(participants)))
    ax2.set_xticklabels(['{}\n({:+.0f})'.format(p, mjoa_delta_by_participant[p]) for p in participants])
    ax2.annotate('Δ mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                 xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
    ax2.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')
    ax2.set_ylabel('% IVD volume above threshold')
    ax2.set_title('IVD strain ({}, patient-specific threshold): PreOp (no preload) vs PostOp, '
                   'by fusion status'.format(name.upper()))

    threshold_line_handle = [Line2D([0], [0], color=threshold_color, linewidth=1.0,
                                     label='{} threshold value'.format(name.upper()))]
    ax2.legend(handles=all_handles + [blank] + threshold_line_handle,
               loc='center left', bbox_to_anchor=(1.12, 0.5), frameon=False)

    fig2.tight_layout()

    pw_thresh_plot_path = os.path.join(
        OUT_DIR,
        'multipatient_ivd_plot_{}_prepost_sortedbymjoachange_patientwise_w_thresholds.pdf'.format(name))
    fig2.savefig(pw_thresh_plot_path, bbox_inches='tight')
    plt.close(fig2)
    pw_plot_paths.append(pw_thresh_plot_path)

print()
print("Summary saved: {}".format(pw_summary_path))
for p in pw_plot_paths:
    print("Plot saved: {}".format(p))

# ============================================================
# Fusion-patient PreOp vs PostOp diff, both threshold versions, as one
# Obsidian-ready markdown table. Flexion/Extension kept separate - never average them.
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
combined['threshold_pw'] = combined['participant'].map(
    lambda p: PATIENT_THRESHOLDS[p][PRINT_PERCENTILE])
combined['threshold_global'] = COHORT_THRESHOLDS[PRINT_PERCENTILE]

print()
print("Fusion-patient PreOp vs PostOp diff - {} (Obsidian-ready markdown):".format(PRINT_PERCENTILE.upper()))
print()
print("| Participant | Condition | Threshold PW | PreOp PW (%) | PostOp PW (%) | Delta PW (pp) "
      "| Threshold Global | PreOp Global (%) | PostOp Global (%) | Delta Global (pp) |")
print("|---|---|---|---|---|---|---|---|---|---|")
for _, r in combined.iterrows():
    print("| {} | {} | {:.4f} | {:.2f} | {:.2f} | {:+.2f} | {:.4f} | {:.2f} | {:.2f} | {:+.2f} |".format(
        r['participant'], r['loading_condition'], r['threshold_pw'],
        r[STATE_PREOP + '_pw'], r[STATE_POSTOP + '_pw'], r['delta_pw'],
        r['threshold_global'],
        r[STATE_PREOP + '_global'], r[STATE_POSTOP + '_global'], r['delta_global']))

# ============================================================
# Standalone table: just the threshold VALUES (patient-wise per fusion
# patient, plus the single global value for reference) - separate from the
# PreOp/PostOp %-above table above.
# ============================================================
print()
print("Threshold values - {} (Obsidian-ready markdown):".format(PRINT_PERCENTILE.upper()))
print()
print("| Participant | Threshold PW (MPS) | Threshold Global (MPS) |")
print("|---|---|---|")
for fusion_p in sorted(FUSION_PARTICIPANTS, key=lambda p: int(p[1:])):
    print("| {} | {:.4f} | {:.4f} |".format(
        fusion_p, PATIENT_THRESHOLDS[fusion_p][PRINT_PERCENTILE], COHORT_THRESHOLDS[PRINT_PERCENTILE]))
