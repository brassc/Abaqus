"""
mps_analysis.py - Analyse output from mps_extract.py.
Run: python mps_analysis.py
"""

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# USER SETTINGS
# ============================================================
CSV_PATH  = r'D:\Charlotte\ABAQUS\N01-011\Pre-Op\Job-021-N01-011-PreOp-BCminus0pt175wEVOL\Job-021-N01-011-PreOp-BCminus0pt175wEVOL_0pt30_Site1_Site2_Site3_Site4_mps.csv'
THRESHOLD = 0.1
SET_NAMES = ['PART-1-1.P44;GM', 'PART-1-1.P45;WM']
# ============================================================

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'],
                     'font.size': 10, 'figure.dpi': 300})

df = pd.read_csv(CSV_PATH)

# MPS95: volume-weighted 95th percentile. Elements sorted by MPS; MPS95 is the
# strain value below which 95% of tissue volume lies. Physically correct when
# combining sets with different mesh densities — each element weighted by volume,
# not counted equally. Threshold-free alternative to CSDM.
# Ref: Kleiven, S. (2007). Predictors for traumatic brain injuries evaluated
# through accident reconstructions. Ann. Adv. Automot. Med., 51, 81-92.
def volume_weighted_percentile(grp, p=0.95):
    s = grp.sort_values('mps')
    cumvol = s['volume'].cumsum()
    return s.loc[cumvol >= p * s['volume'].sum(), 'mps'].iloc[0]

records = []
for (set_name, frame_idx), grp in df.groupby(['set_name', 'frame_index']):
    total = grp['volume'].sum()
    above = grp.loc[grp['mps'] >= THRESHOLD, 'volume'].sum()
    records.append({'set_name': set_name, 'frame_index': frame_idx,
                    'frame_value': grp['frame_value'].iloc[0],
                    'total_volume': total, 'volume_above': above,
                    'pct_above': 100.0 * above / total if total > 0 else 0.0,
                    'mps95': volume_weighted_percentile(grp)})

summary = pd.DataFrame(records).sort_values(['set_name', 'frame_index'])

# Combined summary across all sets per frame. MPS95 and volume stats are computed
# from the full element distribution (all sets pooled), not averaged from per-set values.
combined_records = []
for frame_idx, grp in df.groupby('frame_index'):
    total = grp['volume'].sum()
    above = grp.loc[grp['mps'] >= THRESHOLD, 'volume'].sum()
    combined_records.append({'set_name': 'COMBINED', 'frame_index': frame_idx,
                             'frame_value': grp['frame_value'].iloc[0],
                             'total_volume': total, 'volume_above': above,
                             'pct_above': 100.0 * above / total if total > 0 else 0.0,
                             'mps95': volume_weighted_percentile(grp)})

combined = pd.DataFrame(combined_records).sort_values('frame_index')
summary_full = pd.concat([summary, combined], ignore_index=True)

# Terminal summary (last frame)
print("\nLast-frame summary (threshold = {}):".format(THRESHOLD))
last = summary_full.loc[summary_full.groupby('set_name')['frame_index'].idxmax()]
for _, row in last.iterrows():
    print("  {}  total={:.4e}  above={:.4e}  ({:.1f}%)  MPS95={:.4f}".format(
        row['set_name'], row['total_volume'], row['volume_above'], row['pct_above'],
        row['mps95']))

out_dir = os.path.dirname(CSV_PATH)

# Summary CSV
summary_full.to_csv(os.path.join(out_dir, 'mps_summary.csv'), index=False)

# # Time-history of volume above threshold (CSDM-style) - commented out, use MPS95 plot instead
# fig, ax = plt.subplots(figsize=(6, 4))
# for set_name, grp in summary.groupby('set_name'):
#     ax.plot(grp['frame_value'], grp['volume_above'], label=set_name)
# ax.set_xlabel('Frame value')
# ax.set_ylabel('Volume above MPS >= {}'.format(THRESHOLD))
# ax.legend()
# fig.tight_layout()
# fig.savefig(os.path.join(out_dir, 'mps_volume_threshold.pdf'), bbox_inches='tight')
# plt.close(fig)

# MPS95 time-history
plot_config = {
    SET_NAMES[0]: ('Grey Matter Only',  ':',  '#9dc3e6'),
    SET_NAMES[1]: ('White Matter Only', '--', '#2e75b6'),
    'COMBINED':   ('Combined MPS95',    '-',  '#1a3a5c'),
}
fig, ax = plt.subplots(figsize=(6, 4))
for set_name, (label, ls, colour) in plot_config.items():
    grp = summary_full[summary_full['set_name'] == set_name].sort_values('frame_index')
    ax.plot(grp['frame_value'], grp['mps95'], label=label, linestyle=ls, color=colour)
ax.set_xlabel('Frame value')
ax.set_ylabel('MPS95')
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'mps95_time_history.pdf'), bbox_inches='tight')
plt.close(fig)

# Histogram at last frame
sets = df['set_name'].unique()
fig, axes = plt.subplots(1, len(sets), figsize=(5 * len(sets), 4), squeeze=False)
last_frame = df['frame_index'].max()
for i, set_name in enumerate(sets):
    mps_vals = df.loc[(df['set_name'] == set_name) & (df['frame_index'] == last_frame), 'mps']
    axes[0][i].hist(mps_vals, bins=50, color='steelblue', edgecolor='none')
    axes[0][i].axvline(THRESHOLD, color='crimson', linestyle='--', label='Threshold')
    axes[0][i].set_xlabel('MPS')
    axes[0][i].set_ylabel('Element count')
    axes[0][i].set_title(set_name)
    axes[0][i].legend()
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'mps_histogram.pdf'), bbox_inches='tight')
plt.close(fig)

print("Done. Outputs saved to: {}".format(out_dir))
