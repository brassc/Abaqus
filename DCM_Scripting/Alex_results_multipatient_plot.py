"""
Alex_results_multipatient_plot.py - Pool MPS data across all patients to
compute cohort-wide volume-weighted 90th/95th/99th percentile MPS thresholds
(T90, T95, T99), then plot, per patient, the % of their own cord volume that
exceeds each threshold. X axis is the anonymized participant number.

Reads DCM_Scripting/id_map.csv (gitignored - real patient IDs, never
committed) to find each patient's raw per-frame extraction CSV (from
Alex_results_extraction.py, also never committed - see mps.csv location
in id_map.csv 'csv_path' column). Only the small aggregated summary,
keyed by anonymized participant number (e.g. 'P6'), and the resulting
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

per_patient = {}
missing = []
for _, row in id_map.iterrows():
    csv_path = str(row.get('csv_path', '')).strip()
    if not csv_path or csv_path.lower() == 'nan' or not os.path.isfile(csv_path):
        missing.append(row)
        continue
    df = pd.read_csv(csv_path)
    per_patient[int(row['participant'])] = reduce_to_frame_mode(df, FRAME_MODE)

if missing:
    print("Skipping {} participant(s) with no csv_path set in id_map.csv (or file not found):".format(len(missing)))
    for row in missing:
        print("  P{}".format(int(row['participant'])))

if not per_patient:
    raise SystemExit("No patient data loaded - fill in csv_path in id_map.csv first.")

# Pool ALL patients' reduced (element, mps, volume) rows to compute cohort-wide
# volume-weighted T90/T95/T99. Elements weighted by volume, not counted equally -
# correct across patients with different cord sizes/mesh densities.
pooled = pd.concat(per_patient.values(), ignore_index=True)
t90 = volume_weighted_percentile(pooled, p=0.90)
t95 = volume_weighted_percentile(pooled, p=0.95)
t99 = volume_weighted_percentile(pooled, p=0.99)
print("Cohort-pooled thresholds (FRAME_MODE='{}'): T90={:.4f}  T95={:.4f}  T99={:.4f}".format(
    FRAME_MODE, t90, t95, t99))

records = []
for participant, df in sorted(per_patient.items()):
    records.append({
        'participant':    'P{}'.format(participant),
        'pct_above_t90':  pct_volume_above(df, t90),
        'pct_above_t95':  pct_volume_above(df, t95),
        'pct_above_t99':  pct_volume_above(df, t99),
    })
summary = pd.DataFrame(records)

summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary.csv')
summary.to_csv(summary_path, index=False)
print(summary.to_string(index=False))

fig, ax = plt.subplots(figsize=(6, 4.5))
x = range(len(summary))
ax.scatter(x, summary['pct_above_t90'], label='% volume >= T90 ({:.4f})'.format(t90),
           color='#548235', marker='s', s=60, zorder=3)
ax.scatter(x, summary['pct_above_t95'], label='% volume >= T95 ({:.4f})'.format(t95),
           color='#2e75b6', marker='o', s=60, zorder=3)
ax.scatter(x, summary['pct_above_t99'], label='% volume >= T99 ({:.4f})'.format(t99),
           color='#c00000', marker='^', s=60, zorder=3)
ax.set_xticks(list(x))
ax.set_xticklabels(summary['participant'])
ax.set_xlabel('Participant')
ax.set_ylabel('% cord volume above threshold')
ax.legend(frameon=False)
fig.tight_layout()

plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot.pdf')
fig.savefig(plot_path, bbox_inches='tight')
plt.close(fig)

print("Summary saved: {}".format(summary_path))
print("Plot saved: {}".format(plot_path))
