"""
plot_prepost_sortedbymjoachange.py - Standalone, fast replot of
multipatient_mps_plot_prepost_sortedbymjoachange.pdf (the PreOp-NoPreload vs
PostOp, mJOA-change-ordered threshold comparison) from the cached per-element
(mps, volume) data written by "Alex_results_multipatient_plot - compare_threshold.py".

Reads ONLY id_map.csv (for mJOA) and the small cache CSV - never the raw
per-frame extraction CSVs - so changing DELTA_THRESHOLDS/colors below and
re-running this script takes seconds instead of redoing the full pipeline.

Run: python plot_prepost_sortedbymjoachange.py
Requires the cache to exist first - run the full
"Alex_results_multipatient_plot - compare_threshold.py" at least once (or
after the underlying per-patient data changes) to (re)generate it.
"""

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from mps_common import pct_volume_above, PLOT_STYLE

# ============================================================
# USER SETTINGS - edit these and re-run for fast iteration
# ============================================================
ID_MAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'id_map.csv')
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Must match the FRAME_MODE the cache was generated with (selects the cache file).
FRAME_MODE = 'peak'

DELTA_THRESHOLDS = {'t0p01': 0.01, 't0p02': 0.02, 't0p03': 0.03}
DELTA_THRESHOLD_COLORS = {
    't0p01': '#2a78d6',   # blue
    't0p02': '#eb6834',   # orange
    't0p03': '#1baf7a',   # aqua
}

CONDITION_MARKERS = {'flexion': 's', 'extension': '^'}
DELTA_TOP_STATE = 'preop-nopreload'     # top row: PreOp WITHOUT simulated preload
DELTA_BOTTOM_STATE = 'postop'           # bottom row: PostOp (unchanged)
DELTA_STATE_TITLES = {DELTA_TOP_STATE: 'PreOp (no preload)', DELTA_BOTTOM_STATE: 'PostOp'}

# Fusion patients (surgical detail - only shaded on the PostOp row, since
# fusion is a post-operative property and has no PreOp meaning).
DELTA_FUSION_PARTICIPANTS = {'P1', 'P5', 'P7', 'P8'}
# ============================================================

plt.rcParams.update(PLOT_STYLE)

delta_cache_path = os.path.join(OUT_DIR, 'cache_prepost_sortedbymjoachange_{}.csv'.format(FRAME_MODE))
if not os.path.isfile(delta_cache_path):
    raise SystemExit(
        "No cache found at {}. Run "
        "\"Alex_results_multipatient_plot - compare_threshold.py\" once first "
        "(with FRAME_MODE='{}') to generate it.".format(delta_cache_path, FRAME_MODE))

delta_cache_df = pd.read_csv(delta_cache_path)
state_per_patient = {
    (int(p), c, s): grp[['element_label', 'mps', 'volume']].reset_index(drop=True)
    for (p, c, s), grp in delta_cache_df.groupby(['participant', 'loading_condition', 'state'])
}

id_map = pd.read_csv(ID_MAP_PATH, skipinitialspace=True)

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

delta_records = []
for (delta_participant, delta_condition, delta_state), delta_df in sorted(state_per_patient.items()):
    delta_p_label = 'P{}'.format(delta_participant)
    if delta_p_label not in mjoa_delta_by_participant:
        continue
    if delta_state.strip().lower() not in (DELTA_TOP_STATE, DELTA_BOTTOM_STATE):
        continue
    for delta_threshold_name, delta_threshold_val in DELTA_THRESHOLDS.items():
        delta_records.append({
            'participant': delta_p_label,
            'loading_condition': delta_condition,
            'state': delta_state,
            'threshold': delta_threshold_name,
            'pct_above': pct_volume_above(delta_df, delta_threshold_val),
        })
delta_summary = pd.DataFrame(delta_records)

delta_summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary_prepost_sortedbymjoachange.csv')
delta_summary.to_csv(delta_summary_path, index=False)
print(delta_summary.to_string(index=False))

delta_participants = sorted(mjoa_delta_by_participant.keys(),
                             key=lambda p: (mjoa_delta_by_participant[p], int(p[1:])))
delta_x_pos = {p: i for i, p in enumerate(delta_participants)}

# Common y-limits across BOTH rows, so PreOp and PostOp are directly
# comparable at a glance instead of each auto-scaling to its own data range.
delta_y_max = delta_summary['pct_above'].max()
delta_ylim = (0, delta_y_max * 1.1 if delta_y_max > 0 else 1)

delta_fig, (delta_ax_preop, delta_ax_postop) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
delta_ax_by_state = {DELTA_TOP_STATE: delta_ax_preop, DELTA_BOTTOM_STATE: delta_ax_postop}

for delta_state_name, delta_ax in delta_ax_by_state.items():
    delta_state_sub = delta_summary[delta_summary['state'].astype(str).str.strip().str.lower() == delta_state_name]
    for delta_threshold_name in DELTA_THRESHOLDS:
        delta_sub = delta_state_sub[delta_state_sub['threshold'] == delta_threshold_name]
        delta_color = DELTA_THRESHOLD_COLORS[delta_threshold_name]
        for delta_participant_label, grp in delta_sub.groupby('participant'):
            if len(grp) == 2:
                xp = delta_x_pos[delta_participant_label]
                delta_ax.vlines(xp, grp['pct_above'].min(), grp['pct_above'].max(),
                                 color=delta_color, linewidth=1.0, alpha=0.5, zorder=2)
        for delta_condition_name, grp in delta_sub.groupby('loading_condition'):
            marker = CONDITION_MARKERS.get(delta_condition_name.strip().lower(), 'o')
            xs = [delta_x_pos[p] for p in grp['participant']]
            delta_ax.scatter(xs, grp['pct_above'], color=delta_color, marker=marker, s=55, zorder=3)
    if delta_state_name == DELTA_BOTTOM_STATE:
        for delta_fusion_p in DELTA_FUSION_PARTICIPANTS:
            if delta_fusion_p in delta_x_pos:
                xp = delta_x_pos[delta_fusion_p]
                delta_ax.axvspan(xp - 0.5, xp + 0.5, color='orange', alpha=0.2, zorder=0)
    delta_ax.set_ylim(delta_ylim)
    delta_ax.set_ylabel('% cord volume\nabove threshold')
    delta_ax.set_title(DELTA_STATE_TITLES[delta_state_name])

delta_threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=DELTA_THRESHOLD_COLORS[name],
                                   label='{:.2f}'.format(val))
                            for name, val in DELTA_THRESHOLDS.items()]
delta_condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                            for cond, marker in CONDITION_MARKERS.items()]
delta_fusion_handle = [Patch(facecolor='orange', alpha=0.2, label='Fusion (PostOp)')]
delta_blank_handle = Line2D([0], [0], linestyle='none', marker='None', label='')

delta_all_handles = (
    [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + delta_threshold_handles +
    [delta_blank_handle] +
    [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + delta_condition_handles +
    [delta_blank_handle] + delta_fusion_handle
)
delta_fig.legend(handles=delta_all_handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)

delta_ax_preop.tick_params(labelbottom=False)
delta_ax_postop.set_xticks(range(len(delta_participants)))
delta_ax_postop.set_xticklabels(
    ['{}\n(delta {:+.0f})'.format(p, mjoa_delta_by_participant[p]) for p in delta_participants])
delta_ax_postop.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')
delta_fig.suptitle('% cord volume above MPS threshold, PreOp (no preload) vs PostOp')
delta_fig.tight_layout()

delta_plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot_prepost_sortedbymjoachange.pdf')
delta_fig.savefig(delta_plot_path, bbox_inches='tight')
plt.close(delta_fig)

print("Summary saved: {}".format(delta_summary_path))
print("Plot saved: {}".format(delta_plot_path))
