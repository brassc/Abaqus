"""
mps_group_plot.py - Aggregate pre-op MPS results across patients and plot
volume fraction above threshold vs pre-op mJOA score.

Requires patient_metadata.csv in the same directory.
Run: python mps_group_plot.py
"""

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# USER SETTINGS
# ============================================================
METADATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'patient_metadata.csv')
THRESHOLD     = 0.1
# ============================================================

plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'],
                     'font.size': 10, 'figure.dpi': 300})

meta = pd.read_csv(METADATA_PATH)
meta = meta[meta['model_type'] == 'pre-op'].reset_index(drop=True)

results = []
for _, row in meta.iterrows():
    df = pd.read_csv(row['csv_path'])
    last_frame = df['frame_index'].max()
    grp = df[df['frame_index'] == last_frame]
    total = grp['volume'].sum()
    above = grp.loc[grp['mps'] >= THRESHOLD, 'volume'].sum()
    results.append({'anon_id':        row['anon_id'],
                    'mjoa_pre':       row['mjoa_pre'],
                    'vol_frac_above': above / total if total > 0 else 0.0})

results = pd.DataFrame(results)

fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(results['mjoa_pre'], results['vol_frac_above'], color='steelblue', zorder=3)
for _, r in results.iterrows():
    ax.annotate(r['anon_id'], (r['mjoa_pre'], r['vol_frac_above']),
                textcoords='offset points', xytext=(5, 3), fontsize=8)
ax.set_xlabel('Pre-operative mJOA score')
ax.set_ylabel('Volume fraction above MPS threshold ({})'.format(THRESHOLD))
ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
fig.tight_layout()

out = os.path.join(os.path.dirname(METADATA_PATH), 'mps_group_mjoa.pdf')
fig.savefig(out, bbox_inches='tight')
plt.close(fig)

print(results.to_string(index=False))
print("\nPlot saved: {}".format(out))
