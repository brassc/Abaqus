"""
Alex_results_multipatient_plot.py - Pool MPS data across all patients AND
both loading conditions to compute a single cohort-wide volume-weighted 95th
percentile MPS threshold (T95), then plot, per patient, the % of cord volume
exceeding T95 for flexion and extension separately (paired points per
patient, connected by a grey line - same style as
Data_Extraction/mps_group_plot.py). X axis is the anonymized participant
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

THRESHOLD_PERCENTILE = 0.95   # T95

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
# to compute one cohort-wide volume-weighted T95. Elements weighted by
# volume, not counted equally - correct across patients with different cord
# sizes/mesh densities.
pooled = pd.concat(per_patient.values(), ignore_index=True)
t95 = volume_weighted_percentile(pooled, p=THRESHOLD_PERCENTILE)
print("Cohort-pooled threshold (FRAME_MODE='{}'): T95={:.4f}".format(FRAME_MODE, t95))

records = []
for (participant, condition), df in sorted(per_patient.items()):
    records.append({
        'participant':       'P{}'.format(participant),
        'loading_condition': condition,
        'pct_above_t95':     pct_volume_above(df, t95),
    })
summary = pd.DataFrame(records)

summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary.csv')
summary.to_csv(summary_path, index=False)
print(summary.to_string(index=False))

# Plot: x = participant (ordered by number), paired flexion/extension points
# connected by a grey vertical line, following mps_group_plot.py's style.
participants = sorted(summary['participant'].unique(), key=lambda p: int(p[1:]))
x_pos = {p: i for i, p in enumerate(participants)}

fig, ax = plt.subplots(figsize=(6, 4.5))

for participant, grp in summary.groupby('participant'):
    if len(grp) == 2:
        ax.vlines(x_pos[participant], grp['pct_above_t95'].min(), grp['pct_above_t95'].max(),
                  color='#aaaaaa', linewidth=1.2, zorder=2)

style = {
    'flexion':   {'color': '#DD8452', 'marker': 's', 'label': 'Flexion'},
    'extension': {'color': '#55A868', 'marker': '^', 'label': 'Extension'},
}
for condition, grp in summary.groupby('loading_condition'):
    s = style.get(condition.strip().lower(), {'color': 'grey', 'marker': 'o', 'label': condition})
    xs = [x_pos[p] for p in grp['participant']]
    ax.scatter(xs, grp['pct_above_t95'], label=s['label'], color=s['color'], marker=s['marker'], s=70, zorder=3)

ax.set_xticks(range(len(participants)))
ax.set_xticklabels(participants)
ax.set_xlabel('Participant')
ax.set_ylabel('% cord volume >= T95 ({:.4f})'.format(t95))
ax.legend(frameon=False)
fig.tight_layout()

plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot.pdf')
fig.savefig(plot_path, bbox_inches='tight')
plt.close(fig)

print("Summary saved: {}".format(summary_path))
print("Plot saved: {}".format(plot_path))
