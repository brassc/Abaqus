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

plt.rcParams.update({
    'font.family':       'Arial',
    'font.size':         11,
    'axes.linewidth':    0.8,
    'xtick.major.size':  4,
    'ytick.major.size':  4,
    'axes.spines.top':   False,
    'axes.spines.right': False,
})

meta = pd.read_csv(METADATA_PATH)
meta = meta[meta['model_type'] == 'pre-op'].reset_index(drop=True)

results = []
missing = []
for _, row in meta.iterrows():
    if not os.path.isfile(row['csv_path']):
        missing.append(row)
        continue
    df = pd.read_csv(row['csv_path'])
    last_frame = df['frame_index'].max()
    grp = df[df['frame_index'] == last_frame]
    total = grp['volume'].sum()
    above = grp.loc[grp['mps'] >= THRESHOLD, 'volume'].sum()
    results.append({'anon_id':            row['anon_id'],
                    'mjoa_pre':           row['mjoa_pre'],
                    'loading_condition':  row['loading_condition'],
                    'vol_frac_above':     above / total if total > 0 else 0.0})

if missing:
    print("\nWARNING: Missing CSV files — paste the following into the Abaqus kernel for each:")
    for row in missing:
        odb = row['csv_path'].replace('_mps.csv', '.odb')
        print("\n  --- {} ({}) ---".format(row['anon_id'], row['loading_condition']))
        print("ODB_PATH = r'{}'".format(odb))
        print("execfile('C:\\\\Users\\\\cmb247\\\\repos\\\\Abaqus\\\\DCM_Scripting\\\\Data_Extraction\\\\mps_extract.py')")
    if not results:
        raise SystemExit("No data loaded. Exiting.")

results = pd.DataFrame(results)

style = {
    'flexion':   {'color': '#DD8452', 'marker': 's', 'label': 'Preload + Flexion'},
    'extension': {'color': '#55A868', 'marker': '^', 'label': 'Preload + Extension'},
}

fig, ax = plt.subplots(figsize=(5.5, 4.5))

# Grey vertical line connecting flexion and extension per patient
for anon_id, grp in results.groupby('anon_id'):
    if len(grp) == 2:
        ax.vlines(grp['mjoa_pre'].iloc[0],
                  grp['vol_frac_above'].min(), grp['vol_frac_above'].max(),
                  color='#aaaaaa', linewidth=1.2, zorder=2)

# Scatter points
for condition, grp in results.groupby('loading_condition'):
    s = style.get(condition, {'color': 'grey', 'marker': 'o', 'label': condition})
    ax.scatter(grp['mjoa_pre'], grp['vol_frac_above'],
               label=s['label'], color=s['color'], marker=s['marker'], s=70, zorder=3)

# Patient label once per patient, above the higher of the two points
for anon_id, grp in results.groupby('anon_id'):
    y_max = grp['vol_frac_above'].max()
    mjoa  = grp['mjoa_pre'].iloc[0]
    ax.annotate(anon_id, xy=(mjoa, y_max), xytext=(0, 7),
                textcoords='offset points', ha='center', fontsize=9, color='#444444')

ax.set_ylim(0, grp['vol_frac_above'].max()*1.25)
ax.set_xlim(0, 18)
ax.set_xlabel('Pre-operative mJOA score')
ax.set_ylabel('Volume fraction above MPS threshold ({})'.format(THRESHOLD))
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()

out = os.path.join(os.path.dirname(METADATA_PATH), 'mps_group_mjoa.pdf')
fig.savefig(out, bbox_inches='tight')
plt.close(fig)

print(results.to_string(index=False))
print("\nPlot saved: {}".format(out))
