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

DELTA_THRESHOLDS = {'t0p01': 0.01, 't0p015': 0.015, 't0p02': 0.02, 't0p025': 0.025}
DELTA_THRESHOLD_COLORS = {
    't0p01': '#2a78d6',   # blue
    't0p015': '#1baf7a',  # aqua
    't0p02': '#eb6834',   # orange
    't0p025': '#c00000',  # red
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
                                   label='{:g}'.format(val))
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
    ['{}\n({:+.0f})'.format(p, mjoa_delta_by_participant[p]) for p in delta_participants])
delta_ax_postop.annotate('Δ mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                          xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
delta_ax_postop.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')
delta_fig.suptitle('% cord volume above MPS threshold, PreOp (no preload) vs PostOp')
delta_fig.tight_layout()

delta_plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot_prepost_sortedbymjoachange.pdf')
delta_fig.savefig(delta_plot_path, bbox_inches='tight')
plt.close(delta_fig)

print("Summary saved: {}".format(delta_summary_path))
print("Plot saved: {}".format(delta_plot_path))

# ============================================================
# Second plot: delta (PostOp - PreOp no preload) % cord volume above
# threshold, split into two rows by operation type - Fusion (top) vs
# Decompression-only (bottom). Reuses state_per_patient, id_map,
# mjoa_delta_by_participant, DELTA_THRESHOLDS/COLORS, CONDITION_MARKERS,
# DELTA_FUSION_PARTICIPANTS, OUT_DIR, Line2D, pd, plt, os,
# pct_volume_above already loaded/defined above.
# ============================================================
OP_GROUP_TITLES = {'fusion': 'Laminectomy with Fusion', 'decompression': 'Laminectomy Only'}

op_delta_records = []
op_participants_present = sorted(set(p for (p, _, _) in state_per_patient))
for op_participant in op_participants_present:
    op_p_label = 'P{}'.format(op_participant)
    if op_p_label not in mjoa_delta_by_participant:
        continue
    op_group = 'fusion' if op_p_label in DELTA_FUSION_PARTICIPANTS else 'decompression'
    for op_condition in ('Flexion', 'Extension'):
        op_pre_df = state_per_patient.get((op_participant, op_condition, 'PreOp-NoPreload'))
        op_post_df = state_per_patient.get((op_participant, op_condition, 'PostOp'))
        if op_pre_df is None or op_post_df is None:
            continue
        for op_threshold_name, op_threshold_val in DELTA_THRESHOLDS.items():
            op_pct_pre = pct_volume_above(op_pre_df, op_threshold_val)
            op_pct_post = pct_volume_above(op_post_df, op_threshold_val)
            op_delta_records.append({
                'participant': op_p_label,
                'group': op_group,
                'loading_condition': op_condition,
                'threshold': op_threshold_name,
                'pct_above_pre': op_pct_pre,
                'pct_above_post': op_pct_post,
                'delta': op_pct_post - op_pct_pre,
            })
op_delta_summary = pd.DataFrame(op_delta_records)

op_summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary_delta_by_operation_2row.csv')
op_delta_summary.to_csv(op_summary_path, index=False)
print(op_delta_summary.to_string(index=False))

# Numeric evaluation aid: mean/median delta per threshold per group, to help
# judge which threshold best separates Fusion from Decompression-only.
print()
print("Mean/median delta by threshold and group:")
op_agg = op_delta_summary.groupby(['threshold', 'group'])['delta'].agg(['mean', 'median']).reset_index()
print(op_agg.to_string(index=False))

# Each row gets its own x-ordering (ascending by mJOA change, tie-broken by
# participant number) since Fusion/Decompression are disjoint patient sets.
op_group_participants = {}
op_group_x_pos = {}
for op_group_name in ('fusion', 'decompression'):
    op_grp_participants = sorted(
        (p for p in mjoa_delta_by_participant
         if p in op_delta_summary.loc[op_delta_summary['group'] == op_group_name, 'participant'].unique()),
        key=lambda p: (mjoa_delta_by_participant[p], int(p[1:])))
    op_group_participants[op_group_name] = op_grp_participants
    op_group_x_pos[op_group_name] = {p: i for i, p in enumerate(op_grp_participants)}

# Shared y-limits across both rows (symmetric padding around the combined
# min/max delta) so Fusion and Decompression are directly comparable.
op_delta_min = op_delta_summary['delta'].min()
op_delta_max = op_delta_summary['delta'].max()
op_delta_pad = 0.1 * max(abs(op_delta_min), abs(op_delta_max), 1e-9)
op_delta_ylim = (op_delta_min - op_delta_pad, op_delta_max + op_delta_pad)

op_fig, (op_ax_fusion, op_ax_decomp) = plt.subplots(2, 1, figsize=(9, 8))
op_ax_by_group = {'fusion': op_ax_fusion, 'decompression': op_ax_decomp}

for op_group_name, op_ax in op_ax_by_group.items():
    op_x_pos = op_group_x_pos[op_group_name]
    op_group_sub = op_delta_summary[op_delta_summary['group'] == op_group_name]
    for op_threshold_name in DELTA_THRESHOLDS:
        op_threshold_sub = op_group_sub[op_group_sub['threshold'] == op_threshold_name]
        op_color = DELTA_THRESHOLD_COLORS[op_threshold_name]
        for op_participant_label, op_grp in op_threshold_sub.groupby('participant'):
            if len(op_grp) == 2:
                op_xp = op_x_pos[op_participant_label]
                op_ax.vlines(op_xp, op_grp['delta'].min(), op_grp['delta'].max(),
                             color=op_color, linewidth=1.0, alpha=0.5, zorder=2)
        for op_condition_name, op_grp in op_threshold_sub.groupby('loading_condition'):
            op_marker = CONDITION_MARKERS.get(op_condition_name.strip().lower(), 'o')
            op_xs = [op_x_pos[p] for p in op_grp['participant']]
            op_ax.scatter(op_xs, op_grp['delta'], color=op_color, marker=op_marker, s=55, zorder=3)
    op_ax.axhline(0, color='gray', linestyle='--', linewidth=0.8, zorder=1)
    op_ax.set_ylim(op_delta_ylim)
    op_ax.set_ylabel('Δ % cord volume\nabove threshold\n(PostOp - PreOp)')
    op_ax.set_title(OP_GROUP_TITLES[op_group_name])
    op_ax.set_xticks(range(len(op_group_participants[op_group_name])))
    op_ax.set_xticklabels(
        ['{}\n({:+.0f})'.format(p, mjoa_delta_by_participant[p])
         for p in op_group_participants[op_group_name]])
    op_ax.annotate('Δ mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                    xytext=(-12, -26), textcoords='offset points', ha='right', va='center')

op_ax_decomp.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')

op_threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=DELTA_THRESHOLD_COLORS[name],
                                label='{:g}'.format(val))
                         for name, val in DELTA_THRESHOLDS.items()]
op_condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                         for cond, marker in CONDITION_MARKERS.items()]
op_blank_handle = Line2D([0], [0], linestyle='none', marker='None', label='')

op_all_handles = (
    [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + op_threshold_handles +
    [op_blank_handle] +
    [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + op_condition_handles
)
op_fig.legend(handles=op_all_handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)

op_fig.suptitle('Δ % cord volume above MPS threshold (PostOp - PreOp no preload), by operation type')
op_fig.tight_layout()

op_plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot_delta_by_operation_2row.pdf')
op_fig.savefig(op_plot_path, bbox_inches='tight')
plt.close(op_fig)

print()
print("Summary saved: {}".format(op_summary_path))
print("Plot saved: {}".format(op_plot_path))

# ============================================================
# Third plot: same delta (PostOp - PreOp no preload) data as the second
# plot above, but as a single combined panel - all 12 patients in one row
# ordered by mJOA change, with the 4 Fusion patients' x-region shaded
# (same orange/axvspan convention as the first plot on this page). Built
# to compare directly against the 2-row version above and pick whichever
# reads better. Reuses op_delta_summary, mjoa_delta_by_participant,
# DELTA_THRESHOLDS/COLORS, CONDITION_MARKERS, DELTA_FUSION_PARTICIPANTS,
# OUT_DIR, Line2D, Patch, pd, plt, os already loaded/defined above.
# ============================================================
op_combined_participants = sorted(mjoa_delta_by_participant.keys(),
                                   key=lambda p: (mjoa_delta_by_participant[p], int(p[1:])))
op_combined_x_pos = {p: i for i, p in enumerate(op_combined_participants)}

op_combined_fig, op_combined_ax = plt.subplots(figsize=(9, 5.5))

for op_threshold_name in DELTA_THRESHOLDS:
    op_threshold_sub = op_delta_summary[op_delta_summary['threshold'] == op_threshold_name]
    op_color = DELTA_THRESHOLD_COLORS[op_threshold_name]
    for op_participant_label, op_grp in op_threshold_sub.groupby('participant'):
        if len(op_grp) == 2:
            op_xp = op_combined_x_pos[op_participant_label]
            op_combined_ax.vlines(op_xp, op_grp['delta'].min(), op_grp['delta'].max(),
                                   color=op_color, linewidth=1.0, alpha=0.5, zorder=2)
    for op_condition_name, op_grp in op_threshold_sub.groupby('loading_condition'):
        op_marker = CONDITION_MARKERS.get(op_condition_name.strip().lower(), 'o')
        op_xs = [op_combined_x_pos[p] for p in op_grp['participant']]
        op_combined_ax.scatter(op_xs, op_grp['delta'], color=op_color, marker=op_marker, s=55, zorder=3)

for op_fusion_p in DELTA_FUSION_PARTICIPANTS:
    if op_fusion_p in op_combined_x_pos:
        op_xp = op_combined_x_pos[op_fusion_p]
        op_combined_ax.axvspan(op_xp - 0.5, op_xp + 0.5, color='orange', alpha=0.2, zorder=0)

op_combined_ax.axhline(0, color='gray', linestyle='--', linewidth=0.8, zorder=1)
op_combined_ax.set_ylim(op_delta_ylim)
op_combined_ax.set_ylabel('Δ % cord volume above threshold\n(PostOp - PreOp)')
op_combined_ax.set_xticks(range(len(op_combined_participants)))
op_combined_ax.set_xticklabels(
    ['{}\n({:+.0f})'.format(p, mjoa_delta_by_participant[p]) for p in op_combined_participants])
op_combined_ax.annotate('Δ mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                         xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
op_combined_ax.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')

op_combined_fusion_handle = [Patch(facecolor='orange', alpha=0.2, label='Fusion')]
op_combined_all_handles = (
    [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + op_threshold_handles +
    [op_blank_handle] +
    [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + op_condition_handles +
    [op_blank_handle] + op_combined_fusion_handle
)
op_combined_fig.legend(handles=op_combined_all_handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)

op_combined_fig.suptitle('Δ % cord volume above MPS threshold (PostOp - PreOp no preload), by operation type')
op_combined_fig.tight_layout()

op_combined_plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot_delta_by_operation_combined.pdf')
op_combined_fig.savefig(op_combined_plot_path, bbox_inches='tight')
plt.close(op_combined_fig)

op_combined_summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary_delta_by_operation_combined.csv')
op_delta_summary.to_csv(op_combined_summary_path, index=False)

print()
print("Summary saved: {}".format(op_combined_summary_path))
print("Plot saved: {}".format(op_combined_plot_path))
