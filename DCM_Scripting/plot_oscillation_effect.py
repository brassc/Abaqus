"""
plot_oscillation_effect.py - For the CSF-pulsation oscillation (Step-3 of each
id_map.csv loading_condition='Oscillation' job), compare two different % cord
volume above threshold metrics, per patient:

  Cumulative - each element's own PEAK mps over the whole oscillation history
  (the usual 'peak' reduction used everywhere else in this codebase), i.e. the
  % of cord volume that exceeds the threshold AT SOME POINT during the cycle,
  regardless of whether different elements peak at different phases.

  At-peak    - for each single frame, the % of cord volume exceeding the
  threshold SIMULTANEOUSLY at that frame; then the worst (maximum) frame
  across the whole oscillation history.

Cumulative >= At-peak always, since it's a union over all frames vs. a single
frame's snapshot. The gap between them is the interesting number: it shows how
much of the "cumulative" exceeding volume is actually never under threshold-
exceeding strain all at once - i.e. how spatially/temporally spread the
oscillation's strain excursions are, rather than one coincident high-strain
event.

Thresholds: the Cumulative-vs-At-peak plot uses 0.10/0.15 (same as
plot_preload_effect.py's Plot A - the oscillation jobs all have the preload
applied). The delta (At-peak - Baseline) plot additionally includes 0.02/0.05:
0.10/0.15 alone showed a near-zero delta for most patients there, since the
preload-only baseline already saturates those thresholds, leaving little
headroom for oscillation's small (0.76mm) incremental displacement to push
more volume across - the lower pair checks whether that incremental effect is
visible below that saturation point.

Caches the small derived per-frame/cumulative summary (not the full raw
per-element-per-frame data) to cache_oscillation_cumulative_vs_peak.csv, so
re-running after changing thresholds/colors is fast.

Run: python plot_oscillation_effect.py
"""

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from mps_common import pct_volume_above, PLOT_STYLE

# ============================================================
# USER SETTINGS - edit these and re-run for fast iteration
# ============================================================
ID_MAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'id_map.csv')
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

STATE_OSCILLATION = 'Oscillation'

METRIC_CUMULATIVE = 'cumulative'
METRIC_AT_PEAK = 'at_peak'
METRIC_TITLES = {METRIC_CUMULATIVE: 'Cumulative (element-wise peak over cycle)',
                  METRIC_AT_PEAK: 'At-peak (worst single frame)'}
METRIC_FILLED = {METRIC_CUMULATIVE: True, METRIC_AT_PEAK: False}   # solid vs hollow marker

# All thresholds computed/cached. The first (cumulative vs at-peak) plot only
# displays t0p10/t0p15 (unchanged); the second (delta) plot displays all seven
# - 0.10/0.15 showed a near-zero delta almost everywhere (baseline already
# saturates them for most patients), so the 0.01-0.05 sweep is here to check
# whether oscillation's incremental effect is visible below that saturation
# point. Low thresholds use a sequential green ramp (light->dark as the
# threshold increases) - orange was too close to the 0.15 red to tell apart
# at a glance; 0.10/0.15 keep the blue/red used throughout the other
# multipatient plots in this codebase.
CACHE_THRESHOLDS = {'t0p01': 0.01, 't0p02': 0.02, 't0p03': 0.03, 't0p04': 0.04, 't0p05': 0.05,
                     't0p10': 0.10, 't0p15': 0.15}
CACHE_THRESHOLD_COLORS = {
    't0p01': '#c7e9c0', 't0p02': '#a1d99b', 't0p03': '#74c476', 't0p04': '#31a354', 't0p05': '#006d2c',
    't0p10': '#2e75b6', 't0p15': '#c00000',
}

CUM_PEAK_DISPLAY_THRESHOLDS = {'t0p10': 0.10, 't0p15': 0.15}
DELTA_DISPLAY_THRESHOLDS = CACHE_THRESHOLDS
# ============================================================

plt.rcParams.update(PLOT_STYLE)


def reduce_to_peak(df):
    peak_mps = df.groupby('element_label')['mps'].max()
    volume = df.groupby('element_label')['volume'].first()
    return pd.DataFrame({'mps': peak_mps, 'volume': volume}).reset_index()


# ============================================================
# Cache: small derived summary - per-frame pct_above series (for at-peak) +
# the single cumulative pct_above value, per participant per threshold. NOT
# the full raw per-element-per-frame data (much larger, unnecessary to keep).
# ============================================================
own_cache_path = os.path.join(OUT_DIR, 'cache_oscillation_cumulative_vs_peak.csv')

id_map = pd.read_csv(ID_MAP_PATH, skipinitialspace=True)

if os.path.isfile(own_cache_path):
    print("Loading cached summary: {}".format(own_cache_path))
    cache_df = pd.read_csv(own_cache_path)
else:
    print("No cache at {} yet - building it.".format(own_cache_path))
    osc_rows = id_map[id_map['loading_condition'].astype(str).str.strip().str.lower() == STATE_OSCILLATION.lower()]
    if osc_rows.empty:
        raise SystemExit("No rows with loading_condition == '{}' found in id_map.csv.".format(STATE_OSCILLATION))

    cache_rows = []
    missing = []
    for _, row in osc_rows.iterrows():
        participant = int(row['participant'])
        csv_path = str(row.get('csv_path', '')).strip()
        if not csv_path or csv_path.lower() == 'nan' or not os.path.isfile(csv_path):
            missing.append(participant)
            continue
        raw = pd.read_csv(csv_path)

        # At-peak: per-frame pct_above, for every frame of the oscillation history.
        for frame_idx, frame_df in raw.groupby('frame_index'):
            frame_value = frame_df['frame_value'].iloc[0]
            for name, val in CACHE_THRESHOLDS.items():
                cache_rows.append({
                    'participant': participant, 'metric': METRIC_AT_PEAK,
                    'frame_index': frame_idx, 'frame_value': frame_value,
                    'threshold': name, 'pct_above': pct_volume_above(frame_df, val),
                })

        # Cumulative: element-wise peak over the whole history, then one pct_above per threshold.
        peak_reduced = reduce_to_peak(raw)
        for name, val in CACHE_THRESHOLDS.items():
            cache_rows.append({
                'participant': participant, 'metric': METRIC_CUMULATIVE,
                'frame_index': None, 'frame_value': None,
                'threshold': name, 'pct_above': pct_volume_above(peak_reduced, val),
            })

    if missing:
        print("Skipping {} patient(s) with no csv_path set (or file not found):".format(len(missing)))
        for p in missing:
            print("  P{}".format(p))

    if not cache_rows:
        raise SystemExit("No data loaded - check id_map.csv.")

    cache_df = pd.DataFrame(cache_rows)
    cache_df.to_csv(own_cache_path, index=False)
    print("Cached summary: {}".format(own_cache_path))

target_participants = sorted(cache_df['participant'].unique())

# Pre-op mJOA - pulled directly from the Oscillation rows themselves (same
# clinical value regardless of which simulation variant carries it).
preop_mjoa_by_participant = {}
for _, row in id_map.iterrows():
    if str(row.get('loading_condition', '')).strip().lower() != STATE_OSCILLATION.lower():
        continue
    p_label = 'P{}'.format(int(row['participant']))
    preop_mjoa_by_participant[p_label] = row.get('mJOA', '')

participants = sorted(('P{}'.format(p) for p in target_participants),
                       key=lambda p: (preop_mjoa_by_participant.get(p, float('inf')), int(p[1:])))
x_pos = {p: i for i, p in enumerate(participants)}

# ============================================================
# Summary: cumulative vs at-peak (worst frame) % cord volume above threshold
# ============================================================
records = []
for p in target_participants:
    p_label = 'P{}'.format(p)
    for name in CUM_PEAK_DISPLAY_THRESHOLDS:
        cum_row = cache_df[(cache_df['participant'] == p) & (cache_df['metric'] == METRIC_CUMULATIVE)
                            & (cache_df['threshold'] == name)]
        peak_rows = cache_df[(cache_df['participant'] == p) & (cache_df['metric'] == METRIC_AT_PEAK)
                              & (cache_df['threshold'] == name)]
        if cum_row.empty or peak_rows.empty:
            continue
        cumulative_pct = cum_row['pct_above'].iloc[0]
        at_peak_row = peak_rows.loc[peak_rows['pct_above'].idxmax()]
        records.append({
            'participant': p_label,
            'threshold': name,
            'pct_above_cumulative': cumulative_pct,
            'pct_above_at_peak': at_peak_row['pct_above'],
            'at_peak_frame_index': at_peak_row['frame_index'],
            'at_peak_frame_value': at_peak_row['frame_value'],
            'gap_cumulative_minus_at_peak': cumulative_pct - at_peak_row['pct_above'],
        })
summary = pd.DataFrame(records)

summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary_oscillation_cumulative_vs_peak.csv')
summary.to_csv(summary_path, index=False)
print(summary.to_string(index=False))

print()
print("Mean/median gap (cumulative - at-peak) by threshold:")
print(summary.groupby('threshold')['gap_cumulative_minus_at_peak'].agg(['mean', 'median']).to_string())

# ============================================================
# Plot: Cumulative (solid) vs At-peak (hollow), connected per patient per
# threshold - no loading-condition split (single oscillation job/patient).
# ============================================================
long_df = pd.concat([
    summary[['participant', 'threshold', 'pct_above_cumulative']].rename(
        columns={'pct_above_cumulative': 'pct_above'}).assign(metric=METRIC_CUMULATIVE),
    summary[['participant', 'threshold', 'pct_above_at_peak']].rename(
        columns={'pct_above_at_peak': 'pct_above'}).assign(metric=METRIC_AT_PEAK),
], ignore_index=True)

fig, ax = plt.subplots(figsize=(8, 5.5))

for name in CUM_PEAK_DISPLAY_THRESHOLDS:
    col_df = long_df[long_df['threshold'] == name]
    color = CACHE_THRESHOLD_COLORS[name]
    for participant, grp in col_df.groupby('participant'):
        if len(grp) == 2:
            xp = x_pos[participant]
            ax.vlines(xp, grp['pct_above'].min(), grp['pct_above'].max(), color=color, linewidth=1.0, alpha=0.5, zorder=2)
    for metric, grp in col_df.groupby('metric'):
        filled = METRIC_FILLED.get(metric, True)
        xs = [x_pos[p] for p in grp['participant']]
        if filled:
            ax.scatter(xs, grp['pct_above'], color=color, marker='o', s=55, zorder=3)
        else:
            ax.scatter(xs, grp['pct_above'], facecolors='none', edgecolors=color, marker='o',
                       s=55, linewidths=1.3, zorder=3)

threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=CACHE_THRESHOLD_COLORS[name],
                             label='{:.2f}'.format(val))
                     for name, val in CUM_PEAK_DISPLAY_THRESHOLDS.items()]
metric_handles = [
    Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='black', label='Cumulative (ever exceeds)'),
    Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='none', label='At-peak (worst single frame)'),
]
blank = Line2D([0], [0], linestyle='none', marker='None', label='')

all_handles = (
    [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + threshold_handles +
    [blank] +
    [Line2D([0], [0], linestyle='none', marker='None', label='Metric')] + metric_handles
)
ax.legend(handles=all_handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

ax.set_xticks(range(len(participants)))
ax.set_xticklabels(['{}\n({})'.format(p, preop_mjoa_by_participant.get(p, '?')) for p in participants])
ax.annotate('mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
            xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
ax.set_xlabel('Participant (ordered by pre-op mJOA, ascending)')
ax.set_ylabel('% cord volume above threshold')
ax.set_title('Oscillation: cumulative (ever-exceeds) vs. at-peak (simultaneous) % cord volume')
fig.tight_layout()

plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot_oscillation_cumulative_vs_peak_sortedbypreopmJOA.pdf')
fig.savefig(plot_path, bbox_inches='tight')
plt.close(fig)

print()
print("Summary saved: {}".format(summary_path))
print("Plot saved: {}".format(plot_path))

# ============================================================
# Second plot: Delta (At-peak worst single frame MINUS Baseline) % cord
# volume above threshold - a single diffed value per patient/threshold rather
# than paired absolute points. Baseline is the pre-oscillation starting state
# (frame_index==0 of the same Step-3, i.e. preload only, before any
# oscillation displacement) - already in cache_df as the frame_index==0 rows
# of the METRIC_AT_PEAK per-frame series, so this needs no new raw-CSV reads.
# Reuses cache_df, target_participants, participants, x_pos,
# preop_mjoa_by_participant, PLOT_THRESHOLDS/COLORS, OUT_DIR, Line2D, pd, plt,
# os already loaded/defined above - does not modify anything above this
# point.
# ============================================================
delta_records = []
for p in target_participants:
    p_label = 'P{}'.format(p)
    for name in DELTA_DISPLAY_THRESHOLDS:
        peak_rows = cache_df[(cache_df['participant'] == p) & (cache_df['metric'] == METRIC_AT_PEAK)
                              & (cache_df['threshold'] == name)]
        base_row = peak_rows[peak_rows['frame_index'] == 0]
        if peak_rows.empty or base_row.empty:
            continue
        at_peak_pct = peak_rows['pct_above'].max()
        baseline_pct = base_row['pct_above'].iloc[0]
        delta_records.append({
            'participant': p_label,
            'threshold': name,
            'pct_above_baseline': baseline_pct,
            'pct_above_at_peak': at_peak_pct,
            'delta_at_peak_minus_baseline': at_peak_pct - baseline_pct,
        })
delta_summary = pd.DataFrame(delta_records)

delta_summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary_oscillation_delta_at_peak_vs_baseline.csv')
delta_summary.to_csv(delta_summary_path, index=False)
print(delta_summary.to_string(index=False))

print()
print("Mean/median delta (at-peak - baseline) by threshold:")
print(delta_summary.groupby('threshold')['delta_at_peak_minus_baseline'].agg(['mean', 'median']).to_string())

delta_min = delta_summary['delta_at_peak_minus_baseline'].min()
delta_max = delta_summary['delta_at_peak_minus_baseline'].max()
delta_pad = 0.1 * max(abs(delta_min), abs(delta_max), 1e-9)
delta_ylim = (delta_min - delta_pad, delta_max + delta_pad)

delta_fig, delta_ax = plt.subplots(figsize=(8, 5.5))

for name in DELTA_DISPLAY_THRESHOLDS:
    col_df = delta_summary[delta_summary['threshold'] == name]
    color = CACHE_THRESHOLD_COLORS[name]
    xs = [x_pos[p] for p in col_df['participant']]
    # alpha<1 so overlapping points (common here - many patients cluster near
    # delta=0) show as visibly darker/stacked instead of hiding each other.
    delta_ax.scatter(xs, col_df['delta_at_peak_minus_baseline'], color=color, marker='o', s=55,
                      alpha=0.55, edgecolors='none', zorder=3)

delta_ax.axhline(0, color='gray', linestyle='--', linewidth=0.8, zorder=1)
delta_ax.set_ylim(delta_ylim)
delta_ax.set_xticks(range(len(participants)))
delta_ax.set_xticklabels(['{}\n({})'.format(p, preop_mjoa_by_participant.get(p, '?')) for p in participants])
delta_ax.annotate('mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                   xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
delta_ax.set_xlabel('Participant (ordered by pre-op mJOA, ascending)')
delta_ax.set_ylabel('Δ % cord volume above threshold\n(At-peak - Baseline)')
delta_ax.set_title('Oscillation: Δ % cord volume above threshold, worst single frame vs. pre-oscillation baseline')

delta_threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=CACHE_THRESHOLD_COLORS[name],
                                   label='{:.2f}'.format(val))
                            for name, val in DELTA_DISPLAY_THRESHOLDS.items()]
delta_ax.legend(handles=[Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + delta_threshold_handles,
                loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

delta_fig.tight_layout()

delta_plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot_oscillation_delta_at_peak_vs_baseline_sortedbypreopmJOA.pdf')
delta_fig.savefig(delta_plot_path, bbox_inches='tight')
plt.close(delta_fig)

print()
print("Summary saved: {}".format(delta_summary_path))
print("Plot saved: {}".format(delta_plot_path))
