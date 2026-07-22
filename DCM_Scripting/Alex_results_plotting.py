"""
Alex_results_plotting.py - For ONE representative patient's extraction
CSV (from Alex_results_extraction.py), plot how the % of cord volume above
the 95th/99th percentile MPS threshold changes over the simulation frames.

Thresholds (T95, T99) are computed once from the LAST frame's volume-weighted
distribution, then tracked across all frames. Use this to decide whether
last-frame extraction is adequate, or whether the exceedance peaks mid-
simulation and drops off by the end (in which case switch to a
peak-over-simulation-time reduction for the full-cohort analysis).

Run: python Alex_results_plotting.py
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
CSV_PATH = r'D:\Charlotte\ABAQUS\N31-038\Job-009-N31-038-PreOpv11-BC0pt35\Job-009-N31-038-PreOpv11-BC0pt35_0pt30_site1_site2_site3_mps.csv'
# ============================================================

plt.rcParams.update(PLOT_STYLE)

df = pd.read_csv(CSV_PATH)

last_frame_idx = df['frame_index'].max()
last_frame_df = df[df['frame_index'] == last_frame_idx]

t95 = volume_weighted_percentile(last_frame_df, p=0.95)
t99 = volume_weighted_percentile(last_frame_df, p=0.99)
print("Thresholds from last frame ({}): T95={:.4f}  T99={:.4f}".format(last_frame_idx, t95, t99))

records = []
for frame_idx, grp in df.groupby('frame_index'):
    records.append({
        'frame_index': frame_idx,
        'frame_value': grp['frame_value'].iloc[0],
        'pct_above_t95': pct_volume_above(grp, t95),
        'pct_above_t99': pct_volume_above(grp, t99),
    })
history = pd.DataFrame(records).sort_values('frame_index')

peak_t95_idx = history['pct_above_t95'].idxmax()
peak_t99_idx = history['pct_above_t99'].idxmax()
print("Peak %% above T95 at frame {} (last frame is {})".format(
    history.loc[peak_t95_idx, 'frame_index'], last_frame_idx))
print("Peak %% above T99 at frame {} (last frame is {})".format(
    history.loc[peak_t99_idx, 'frame_index'], last_frame_idx))
if peak_t95_idx == history.index[-1] and peak_t99_idx == history.index[-1]:
    print("-> Exceedance peaks at the final frame: last-frame extraction (FRAME_MODE='last') looks adequate.")
else:
    print("-> Exceedance peaks before the final frame: consider peak-over-simulation-time (FRAME_MODE='peak').")

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(history['frame_value'], history['pct_above_t95'], label='% volume >= T95 ({:.4f})'.format(t95), color='#2e75b6')
ax.plot(history['frame_value'], history['pct_above_t99'], label='% volume >= T99 ({:.4f})'.format(t99), color='#c00000')
ax.set_xlabel('Frame value')
ax.set_ylabel('% cord volume above threshold')
ax.legend(frameon=False)
fig.tight_layout()

out_dir = os.path.dirname(CSV_PATH)
out_path = os.path.join(out_dir, 'time_history_diagnostic.pdf')
fig.savefig(out_path, bbox_inches='tight')
plt.close(fig)

print("Plot saved: {}".format(out_path))
