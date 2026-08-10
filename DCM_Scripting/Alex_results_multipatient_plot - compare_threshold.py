"""
Alex_results_multipatient_plot.py - Pool MPS data across all patients AND
both loading conditions to compute cohort-wide volume-weighted 90th/95th/99th
percentile MPS thresholds (T90, T95, T99), then plot, per patient, the % of
cord volume exceeding each threshold for flexion and extension separately
(paired points per patient per threshold, connected by a grey line - same
style as Data_Extraction/mps_group_plot.py). Marker shape encodes loading
condition, color encodes threshold. X axis is the anonymized participant
number.

Reads DCM_Scripting/id_map.csv (gitignored - real patient IDs, never
committed) to find each patient's raw per-frame extraction CSV per loading
condition (from Alex_results_extraction.py, also never committed - see
'csv_path' column). Only the small aggregated summary, keyed by anonymized
participant number (e.g. 'P6') and loading condition, and the resulting
plot are written into the repo.

Run: python Alex_results_multipatient_plot.py
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

# 'last'  - use each element's MPS at the last frame of the step
# 'peak'  - use each element's peak (max-ever) MPS across all frames
# Decide using Alex_results_plotting.py's single-patient time-history
# diagnostic before running this at full cohort scale.
FRAME_MODE = 'peak'

THRESHOLD_PERCENTILES = {'t90': 0.90, 't95': 0.95, 't99': 0.99}
THRESHOLD_COLORS = {'t90': '#548235', 't95': '#2e75b6', 't99': '#c00000'}

# Fixed MPS thresholds (not derived from the data) for a second comparison plot
MANUAL_THRESHOLDS = {'t0p05': 0.05, 't0p10': 0.10, 't0p15': 0.15, 't0p20': 0.20}
MANUAL_THRESHOLD_COLORS = {'t0p05': '#548235', 't0p10': '#2e75b6', 't0p15': '#c00000', 't0p20': '#7030a0'}

CONDITION_MARKERS = {'flexion': 's', 'extension': '^'}

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
# ============================================================

plt.rcParams.update(PLOT_STYLE)


def reduce_to_frame_mode(df, mode):
    if mode == 'last':
        last_frame_idx = df['frame_index'].max()
        return df[df['frame_index'] == last_frame_idx][['element_label', 'mps', 'volume']]
    elif mode == 'peak':
        peak_mps = df.groupby('element_label')['mps'].max()
        volume = df.groupby('element_label')['volume'].first()
        return pd.DataFrame({'mps': peak_mps, 'volume': volume}).reset_index()
    else:
        raise ValueError("FRAME_MODE must be 'last' or 'peak', got '{}'".format(mode))


def build_summary(per_patient, thresholds):
    records = []
    for (participant, condition), df in sorted(per_patient.items()):
        record = {'participant': 'P{}'.format(participant), 'loading_condition': condition}
        for name, val in thresholds.items():
            record['pct_above_{}'.format(name)] = pct_volume_above(df, val)
        records.append(record)
    return pd.DataFrame(records)


def make_comparison_plot(summary, thresholds, threshold_colors, threshold_labels, plot_path):
    # x = participant (ordered by number). For each threshold, paired
    # flexion/extension points are connected by a thin line in that
    # threshold's color; marker shape encodes loading condition, marker/line
    # color encodes threshold - following mps_group_plot.py's paired-point
    # style, extended to multiple thresholds at once.
    participants = sorted(summary['participant'].unique(), key=lambda p: int(p[1:]))
    x_pos = {p: i for i, p in enumerate(participants)}

    fig, ax = plt.subplots(figsize=(7, 5))

    for name in thresholds:
        col = 'pct_above_{}'.format(name)
        color = threshold_colors[name]
        for participant, grp in summary.groupby('participant'):
            if len(grp) == 2:
                xp = x_pos[participant]
                ax.vlines(xp, grp[col].min(), grp[col].max(), color=color, linewidth=1.0, alpha=0.5, zorder=2)
        for condition, grp in summary.groupby('loading_condition'):
            marker = CONDITION_MARKERS.get(condition.strip().lower(), 'o')
            xs = [x_pos[p] for p in grp['participant']]
            ax.scatter(xs, grp[col], color=color, marker=marker, s=55, zorder=3)

    # Two separate legends: color -> threshold, marker shape -> loading condition
    threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=threshold_colors[name],
                                 label=threshold_labels[name])
                          for name in thresholds]
    condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                          for cond, marker in CONDITION_MARKERS.items()]

    threshold_legend = ax.legend(handles=threshold_handles, loc='upper left', frameon=False, title='Threshold')
    ax.add_artist(threshold_legend)
    ax.legend(handles=condition_handles, loc='upper right', frameon=False, title='Loading condition')

    ax.set_xticks(range(len(participants)))
    ax.set_xticklabels(participants)
    ax.set_xlabel('Participant')
    ax.set_ylabel('% cord volume above threshold')
    fig.tight_layout()

    fig.savefig(plot_path, bbox_inches='tight')
    plt.close(fig)


id_map = pd.read_csv(ID_MAP_PATH, skipinitialspace=True)

per_patient = {}   # (participant, loading_condition) -> reduced DataFrame
missing = []
for _, row in id_map.iterrows():
    csv_path = str(row.get('csv_path', '')).strip()
    condition = str(row.get('loading_condition', '')).strip()
    if not csv_path or csv_path.lower() == 'nan' or not os.path.isfile(csv_path):
        missing.append(row)
        continue
    df = pd.read_csv(csv_path)
    per_patient[(int(row['participant']), condition)] = reduce_to_frame_mode(df, FRAME_MODE)

if missing:
    print("Skipping {} row(s) with no csv_path set in id_map.csv (or file not found):".format(len(missing)))
    for row in missing:
        cond = str(row.get('loading_condition', '')).strip() or '?'
        print("  P{} ({})".format(int(row['participant']), cond))

if not per_patient:
    raise SystemExit("No patient data loaded - fill in csv_path in id_map.csv first.")

# Pool ALL patients' AND both conditions' reduced (element, mps, volume) rows
# to compute cohort-wide volume-weighted T90/T95/T99. Elements weighted by
# volume, not counted equally - correct across patients with different cord
# sizes/mesh densities.
pooled = pd.concat(per_patient.values(), ignore_index=True)
thresholds = {name: volume_weighted_percentile(pooled, p=p) for name, p in THRESHOLD_PERCENTILES.items()}
print("Cohort-pooled thresholds (FRAME_MODE='{}'): ".format(FRAME_MODE) +
      "  ".join("{}={:.4f}".format(name.upper(), val) for name, val in thresholds.items()))

records = []
for (participant, condition), df in sorted(per_patient.items()):
    record = {'participant': 'P{}'.format(participant), 'loading_condition': condition}
    for name, val in thresholds.items():
        record['pct_above_{}'.format(name)] = pct_volume_above(df, val)
    records.append(record)
summary = pd.DataFrame(records)

summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary_compare_thresholds.csv')
summary.to_csv(summary_path, index=False)
print(summary.to_string(index=False))

# Plot: x = participant (ordered by number). For each threshold, paired
# flexion/extension points are connected by a thin line in that threshold's
# color; marker shape encodes loading condition, marker/line color encodes
# threshold - following mps_group_plot.py's paired-point style, extended to
# three thresholds at once.
participants = sorted(summary['participant'].unique(), key=lambda p: int(p[1:]))
x_pos = {p: i for i, p in enumerate(participants)}

fig, ax = plt.subplots(figsize=(7, 5))

for name in thresholds:
    col = 'pct_above_{}'.format(name)
    color = THRESHOLD_COLORS[name]
    for participant, grp in summary.groupby('participant'):
        if len(grp) == 2:
            xp = x_pos[participant]
            ax.vlines(xp, grp[col].min(), grp[col].max(), color=color, linewidth=1.0, alpha=0.5, zorder=2)
    for condition, grp in summary.groupby('loading_condition'):
        marker = CONDITION_MARKERS.get(condition.strip().lower(), 'o')
        xs = [x_pos[p] for p in grp['participant']]
        ax.scatter(xs, grp[col], color=color, marker=marker, s=55, zorder=3)

# Two separate legends: color -> threshold, marker shape -> loading condition
threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=THRESHOLD_COLORS[name],
                             label='{} ({:.4f})'.format(name.upper(), thresholds[name]))
                      for name in thresholds]
condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                      for cond, marker in CONDITION_MARKERS.items()]

threshold_legend = ax.legend(handles=threshold_handles, loc='upper left', frameon=False, title='Threshold')
ax.add_artist(threshold_legend)
ax.legend(handles=condition_handles, loc='upper right', frameon=False, title='Loading condition')

ax.set_xticks(range(len(participants)))
ax.set_xticklabels(participants)
ax.set_xlabel('Participant')
ax.set_ylabel('% cord volume above threshold')
fig.tight_layout()

plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot_compare_thresholds.pdf')
fig.savefig(plot_path, bbox_inches='tight')
plt.close(fig)

print("Summary saved: {}".format(summary_path))
print("Plot saved: {}".format(plot_path))

# ============================================================
# Second plot: manual fixed MPS thresholds (0.05 / 0.10 / 0.15 / 0.20)
# instead of percentile-derived ones. Reuses per_patient loaded above.
# ============================================================
manual_summary = build_summary(per_patient, MANUAL_THRESHOLDS)
manual_summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary_compare_thresholds_manual_thresholds.csv')
manual_summary.to_csv(manual_summary_path, index=False)
print(manual_summary.to_string(index=False))

manual_labels = {name: '{:.2f}'.format(val) for name, val in MANUAL_THRESHOLDS.items()}
manual_plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot_compare_thresholds_manual_thresholds.pdf')
make_comparison_plot(manual_summary, MANUAL_THRESHOLDS, MANUAL_THRESHOLD_COLORS, manual_labels, manual_plot_path)

print("Summary saved: {}".format(manual_summary_path))
print("Plot saved: {}".format(manual_plot_path))
