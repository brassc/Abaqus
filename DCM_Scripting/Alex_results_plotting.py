"""
Alex_results_plotting.py - For each patient's extraction CSV (from
Alex_results_extraction.py) with a csv_path filled in id_map.csv, plot how
the % of cord volume above the 95th/99th percentile MPS threshold changes
over the simulation frames.

Thresholds (T95, T99) are computed once from the LAST frame's volume-weighted
distribution, then tracked across all frames. Use this to decide whether
last-frame extraction is adequate, or whether the exceedance peaks mid-
simulation and drops off by the end (in which case switch to a
peak-over-simulation-time reduction for the full-cohort analysis).

The frame index where each threshold's exceedance peaks is logged back into
id_map.csv (gitignored - real patient IDs never leave that file) for every
row processed, in one write at the end.

CSV_PATH = None (default) - process every id_map.csv row with a csv_path set
CSV_PATH = '<path>'       - process only that one CSV (original single-patient behavior)

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
if 'CSV_PATH' not in dir():
    CSV_PATH = None   # None = batch mode over every id_map.csv row with a csv_path
ID_MAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'id_map.csv')
# ============================================================

plt.rcParams.update(PLOT_STYLE)

id_map = pd.read_csv(ID_MAP_PATH, skipinitialspace=True)
id_map.columns = [c.strip() for c in id_map.columns]

for col in ('peak_frame_t95', 'peak_frame_t99', 'last_frame_idx'):
    if col not in id_map.columns:
        id_map[col] = pd.NA
    # Force object dtype - column may be inferred as float64/string on read
    # (numeric-looking values + blanks), which rejects int assignment under
    # pandas's strict dtype casting.
    id_map[col] = id_map[col].astype('object')

if CSV_PATH is not None:
    targets = id_map[id_map['csv_path'].astype(str).str.strip() == CSV_PATH.strip()]
    if targets.empty:
        raise SystemExit("CSV_PATH not found in id_map.csv 'csv_path' column - "
                          "add this participant's row first.")
else:
    csv_path_str = id_map['csv_path'].astype(str).str.strip()
    targets = id_map[(csv_path_str != '') & (csv_path_str.str.lower() != 'nan')]

processed = 0
skipped = []

for idx, row in targets.iterrows():
    csv_path = str(row['csv_path']).strip()
    label = 'P{} ({}, {})'.format(int(row['participant']), row.get('loading_condition', '?'),
                                   row.get('State', '?'))
    if not csv_path or csv_path.lower() == 'nan' or not os.path.isfile(csv_path):
        skipped.append(label)
        continue

    print("\n=== {} ===".format(label))

    df = pd.read_csv(csv_path)

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
    peak_frame_t95 = int(history.loc[peak_t95_idx, 'frame_index'])
    peak_frame_t99 = int(history.loc[peak_t99_idx, 'frame_index'])
    print("Peak %% above T95 at frame {} (last frame is {})".format(peak_frame_t95, last_frame_idx))
    print("Peak %% above T99 at frame {} (last frame is {})".format(peak_frame_t99, last_frame_idx))
    if peak_t95_idx == history.index[-1] and peak_t99_idx == history.index[-1]:
        print("-> Exceedance peaks at the final frame: last-frame extraction (FRAME_MODE='last') looks adequate.")
    else:
        print("-> Exceedance peaks before the final frame: consider peak-over-simulation-time (FRAME_MODE='peak').")

    id_map.loc[idx, 'peak_frame_t95'] = peak_frame_t95
    id_map.loc[idx, 'peak_frame_t99'] = peak_frame_t99
    id_map.loc[idx, 'last_frame_idx'] = int(last_frame_idx)
    processed += 1

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(history['frame_value'], history['pct_above_t95'],
            label='% volume >= T95 ({:.4f})'.format(t95), color='#2e75b6')
    ax.plot(history['frame_value'], history['pct_above_t99'],
            label='% volume >= T99 ({:.4f})'.format(t99), color='#c00000')
    ax.set_xlabel('Frame value')
    ax.set_ylabel('% cord volume above threshold')
    ax.set_title(label)
    ax.legend(frameon=False)
    fig.tight_layout()

    out_dir = os.path.dirname(csv_path)
    out_path = os.path.join(out_dir, 'time_history_diagnostic.pdf')
    fig.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
    print("Plot saved: {}".format(out_path))

id_map.to_csv(ID_MAP_PATH, index=False)
print("\nLogged peak frames to id_map.csv for {} row(s).".format(processed))

if skipped:
    print("Skipped {} row(s) with no csv_path set (or file not found):".format(len(skipped)))
    for label in skipped:
        print("  {}".format(label))
