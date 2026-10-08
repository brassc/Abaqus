"""
DCM_main.py - Stage 1: PreOp (no preload) vs PostOp.

Part A: % cord volume above threshold. Part B: blob-size distribution.
Part C: whole-cord LMM. Part D: GM vs WM LMM. F: IVD strain (T95) by fusion status. Output
-> preopnopreloadvspostop_results/ (QQ plots and blob scalar summaries in
diagnostic_plots/).

Run: python DCM_main.py (first run is slow - builds caches from raw data).
"""

import math
import os
import sys
from collections import defaultdict

# Fixes a reproducibility quirk: matplotlib writes its PDF backend's
# /CreationDate from the current time on every run, so re-running this
# script with identical data still produces byte-different PDFs (git sees
# every figure as "modified" even when nothing actually changed). Setting
# SOURCE_DATE_EPOCH pins that timestamp instead - set before matplotlib
# import, and setdefault so an externally-set value isn't clobbered.
os.environ.setdefault('SOURCE_DATE_EPOCH', '1700000000')

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import to_rgba
from matplotlib.ticker import LogLocator, NullFormatter

# ============================================================
# Shared plotting/volume helpers
# ============================================================


def pct_volume_above(df, threshold, mps_col='mps', vol_col='volume'):
    """% of total volume in df with mps_col >= threshold."""
    total = df[vol_col].sum()
    if total <= 0:
        return 0.0
    above = df.loc[df[mps_col] >= threshold, vol_col].sum()
    return 100.0 * above / total


def volume_weighted_percentile(df, p=0.95, mps_col='mps', vol_col='volume'):
    """Volume-weighted percentile: the MPS value below which fraction p of
    total volume lies. Elements weighted by volume, not counted equally -
    correct when combining regions/patients with different mesh densities.

    Ref: Kleiven, S. (2007). Predictors for traumatic brain injuries
    evaluated through accident reconstructions. Ann. Adv. Automot. Med., 51, 81-92.
    """
    s = df.sort_values(mps_col)
    cumvol = s[vol_col].cumsum()
    return s.loc[cumvol >= p * s[vol_col].sum(), mps_col].iloc[0]


PLOT_STYLE = {
    'font.family':       'Arial',
    'font.size':         11,
    'axes.linewidth':    0.8,
    'xtick.major.size':  4,
    'ytick.major.size':  4,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'figure.dpi':        300,
}

# ============================================================
# USER SETTINGS - edit these and re-run for fast iteration
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ID_MAP_PATH = os.path.join(SCRIPT_DIR, 'id_map.csv')

RESULTS_DIR = os.path.join(SCRIPT_DIR, 'preopnopreloadvspostop_results')
DIAG_DIR = os.path.join(RESULTS_DIR, 'diagnostic_plots')
CACHE_DIR = os.path.join(RESULTS_DIR, 'cache')
os.makedirs(DIAG_DIR, exist_ok=True)
os.makedirs(CACHE_DIR, exist_ok=True)
OUT_DIR = RESULTS_DIR   # kept as OUT_DIR below so the ported blocks are unchanged

# Must match the FRAME_MODE the caches were generated with (selects the cache files).
FRAME_MODE = 'peak'

DELTA_THRESHOLDS = {'t0p01': 0.01, 't0p015': 0.015, 't0p02': 0.02, 't0p025': 0.025}
DELTA_THRESHOLD_COLORS = {
    't0p01': '#2a78d6',   # blue
    't0p015': '#1baf7a',  # aqua
    't0p02': '#eb6834',   # orange
    't0p025': '#c00000',  # red
}

PPBLOB_THRESHOLDS = {'t0p01': 0.01, 't0p015': 0.015, 't0p02': 0.02, 't0p025': 0.025, 't0p03': 0.03, 't0p035': 0.035}
PPBLOB_THRESHOLD_COLORS = {
    't0p01': '#2a78d6',    # blue
    't0p015': '#1baf7a',   # aqua
    't0p02': '#eb6834',    # orange
    't0p025': '#e34948',   # red
    't0p03': '#4a3aa7',    # violet
    't0p035': '#8c8c8c',   # grey
}

CONDITION_MARKERS = {'flexion': 's', 'extension': '^'}
BLOB_CONDITION_LINESTYLES = {'flexion': '-', 'extension': '--'}

DELTA_TOP_STATE = 'preop-nopreload'     # top row: PreOp WITHOUT simulated preload
DELTA_BOTTOM_STATE = 'postop'           # bottom row: PostOp (unchanged)
DELTA_STATE_TITLES = {DELTA_TOP_STATE: 'PreOp (no preload)', DELTA_BOTTOM_STATE: 'PostOp'}

# Fusion patients (surgical detail - only shaded on the PostOp row, since
# fusion is a post-operative property and has no PreOp meaning).
DELTA_FUSION_PARTICIPANTS = {'P1', 'P5', 'P7', 'P8'}
# ============================================================

plt.rcParams.update(PLOT_STYLE)

BLOB_FACE_SHARING_MIN_NODES = 4   # >=4 shared nodes approximates shared face (all elements are C3D8)


class _BlobUnionFind(object):
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def find_blobs(exceeding_labels, edges):
    uf = _BlobUnionFind()
    for lbl in exceeding_labels:
        uf.find(lbl)   # register as singleton
    for a, b in edges:
        if a in exceeding_labels and b in exceeding_labels:
            uf.union(a, b)

    groups = defaultdict(list)
    for lbl in exceeding_labels:
        groups[uf.find(lbl)].append(lbl)
    return list(groups.values())


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


def load_adjacency_edges(topology_path, min_shared_nodes=BLOB_FACE_SHARING_MIN_NODES):
    elem_nodes = {}
    with open(topology_path) as f:
        next(f)  # header
        for line in f:
            label_str, nodes_str = line.strip().split(',', 1)
            elem_nodes[int(label_str)] = set(int(n) for n in nodes_str.split(';'))

    node_to_elems = defaultdict(list)
    for elem, nodes in elem_nodes.items():
        for n in nodes:
            node_to_elems[n].append(elem)

    shared_count = defaultdict(int)
    for n, elems in node_to_elems.items():
        elems = sorted(elems)
        for i in range(len(elems)):
            for j in range(i + 1, len(elems)):
                shared_count[(elems[i], elems[j])] += 1

    return [pair for pair, cnt in shared_count.items() if cnt >= min_shared_nodes]


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

delta_participants = sorted(mjoa_delta_by_participant.keys(),
                             key=lambda p: (mjoa_delta_by_participant[p], int(p[1:])))
delta_x_pos = {p: i for i, p in enumerate(delta_participants)}

# ============================================================
# Rows this script needs data for: PreOp-NoPreload/PostOp, Flexion/Extension,
# with a populated csv_path (set by Alex_results_extraction.py).
# ============================================================
_target_rows = id_map[
    id_map['State'].astype(str).str.strip().str.lower().isin((DELTA_TOP_STATE, DELTA_BOTTOM_STATE)) &
    id_map['loading_condition'].astype(str).str.strip().str.lower().isin(('flexion', 'extension')) &
    id_map['csv_path'].astype(str).str.strip().astype(bool)
]

# --- mps/volume cache: build if missing, else read (shared in-memory by
# Parts A, B and C below) ---
mps_cache_path = os.path.join(CACHE_DIR, 'cache_prepost_sortedbymjoachange_{}.csv'.format(FRAME_MODE))
if os.path.isfile(mps_cache_path):
    mps_cache_df = pd.read_csv(mps_cache_path)
else:
    print("No mps/volume cache found - building it from raw per-frame CSVs (slow, one-time)...")
    _mps_rows = []
    for _, _row in _target_rows.iterrows():
        _csv_path = str(_row['csv_path']).strip()
        if not os.path.isfile(_csv_path):
            print("  SKIPPING (file not found): {}".format(_csv_path))
            continue
        _raw_df = pd.read_csv(_csv_path)
        _reduced = reduce_to_frame_mode(_raw_df, FRAME_MODE)
        _reduced = _reduced.copy()
        _reduced.insert(0, 'state', str(_row['State']).strip())
        _reduced.insert(0, 'loading_condition', str(_row['loading_condition']).strip())
        _reduced.insert(0, 'participant', int(_row['participant']))
        _mps_rows.append(_reduced)
    mps_cache_df = pd.concat(_mps_rows, ignore_index=True)
    mps_cache_df.to_csv(mps_cache_path, index=False)
    print("Cached mps/volume data: {}".format(mps_cache_path))

state_per_patient = {
    (int(p), c, s): grp[['element_label', 'mps', 'volume']].reset_index(drop=True)
    for (p, c, s), grp in mps_cache_df.groupby(['participant', 'loading_condition', 'state'])
}

# --- blob adjacency cache: build if missing, else read ---
adjacency_cache_path = os.path.join(CACHE_DIR, 'cache_prepost_blob_adjacency_{}.csv'.format(FRAME_MODE))
if os.path.isfile(adjacency_cache_path):
    adjacency_cache_df = pd.read_csv(adjacency_cache_path)
else:
    print("No blob adjacency cache found - building it from '_topology.csv' files (slow, one-time)...")
    _adjacency_rows = []
    for _, _row in _target_rows.iterrows():
        _csv_path = str(_row['csv_path']).strip()
        _topology_path = _csv_path.replace('_mps.csv', '_topology.csv')
        if not os.path.isfile(_topology_path):
            print("  SKIPPING (no '_topology.csv', re-run Alex_results_extraction.py): {}".format(_topology_path))
            continue
        for _elem_a, _elem_b in load_adjacency_edges(_topology_path):
            _adjacency_rows.append({
                'participant':       int(_row['participant']),
                'loading_condition': str(_row['loading_condition']).strip(),
                'state':             str(_row['State']).strip(),
                'elem_a':            _elem_a,
                'elem_b':            _elem_b,
            })
    adjacency_cache_df = pd.DataFrame(_adjacency_rows)
    adjacency_cache_df.to_csv(adjacency_cache_path, index=False)
    print("Cached blob adjacency edges: {}".format(adjacency_cache_path))

ppblob_adjacency_cache = {
    (int(p), c, s): list(zip(grp['elem_a'], grp['elem_b']))
    for (p, c, s), grp in adjacency_cache_df.groupby(['participant', 'loading_condition', 'state'])
}

print("")
print("=" * 70)
print("PART A: % cord volume above threshold, sorted by change in mJOA")
print("=" * 70)

# ============================================================
# First plot: PreOp-NoPreload vs PostOp, % cord volume above threshold,
# sorted by change in mJOA (non-delta - raw values per state).
# ============================================================
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
        for delta_fusion_p in sorted(DELTA_FUSION_PARTICIPANTS, key=lambda p: int(p[1:])):
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

# 2-row figure (by operation type) - COMMENTED OUT, not part of the
# current reduced output set. op_delta_ylim/op_threshold_handles/
# op_condition_handles/op_blank_handle below are kept active - the combined
# plot (third plot) still needs them.
# op_summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary_delta_by_operation_2row.csv')
# op_delta_summary.to_csv(op_summary_path, index=False)
# print(op_delta_summary.to_string(index=False))
#
# # Numeric evaluation aid: mean/median delta per threshold per group, to help
# # judge which threshold best separates Fusion from Decompression-only.
# print()
# print("Mean/median delta by threshold and group:")
# op_agg = op_delta_summary.groupby(['threshold', 'group'])['delta'].agg(['mean', 'median']).reset_index()
# print(op_agg.to_string(index=False))
#
# # Each row gets its own x-ordering (ascending by mJOA change, tie-broken by
# # participant number) since Fusion/Decompression are disjoint patient sets.
# op_group_participants = {}
# op_group_x_pos = {}
# for op_group_name in ('fusion', 'decompression'):
#     op_grp_participants = sorted(
#         (p for p in mjoa_delta_by_participant
#          if p in op_delta_summary.loc[op_delta_summary['group'] == op_group_name, 'participant'].unique()),
#         key=lambda p: (mjoa_delta_by_participant[p], int(p[1:])))
#     op_group_participants[op_group_name] = op_grp_participants
#     op_group_x_pos[op_group_name] = {p: i for i, p in enumerate(op_grp_participants)}

# Shared y-limits across both rows (symmetric padding around the combined
# min/max delta) so Fusion and Decompression are directly comparable. Still
# needed below by the combined plot.
op_delta_min = op_delta_summary['delta'].min()
op_delta_max = op_delta_summary['delta'].max()
op_delta_pad = 0.1 * max(abs(op_delta_min), abs(op_delta_max), 1e-9)
op_delta_ylim = (op_delta_min - op_delta_pad, op_delta_max + op_delta_pad)

# op_fig, (op_ax_fusion, op_ax_decomp) = plt.subplots(2, 1, figsize=(9, 8))
# op_ax_by_group = {'fusion': op_ax_fusion, 'decompression': op_ax_decomp}
#
# for op_group_name, op_ax in op_ax_by_group.items():
#     op_x_pos = op_group_x_pos[op_group_name]
#     op_group_sub = op_delta_summary[op_delta_summary['group'] == op_group_name]
#     for op_threshold_name in DELTA_THRESHOLDS:
#         op_threshold_sub = op_group_sub[op_group_sub['threshold'] == op_threshold_name]
#         op_color = DELTA_THRESHOLD_COLORS[op_threshold_name]
#         for op_participant_label, op_grp in op_threshold_sub.groupby('participant'):
#             if len(op_grp) == 2:
#                 op_xp = op_x_pos[op_participant_label]
#                 op_ax.vlines(op_xp, op_grp['delta'].min(), op_grp['delta'].max(),
#                              color=op_color, linewidth=1.0, alpha=0.5, zorder=2)
#         for op_condition_name, op_grp in op_threshold_sub.groupby('loading_condition'):
#             op_marker = CONDITION_MARKERS.get(op_condition_name.strip().lower(), 'o')
#             op_xs = [op_x_pos[p] for p in op_grp['participant']]
#             op_ax.scatter(op_xs, op_grp['delta'], color=op_color, marker=op_marker, s=55, zorder=3)
#     op_ax.axhline(0, color='gray', linestyle='--', linewidth=0.8, zorder=1)
#     op_ax.set_ylim(op_delta_ylim)
#     op_ax.set_ylabel('Δ % cord volume\nabove threshold\n(PostOp - PreOp)')
#     op_ax.set_title(OP_GROUP_TITLES[op_group_name])
#     op_ax.set_xticks(range(len(op_group_participants[op_group_name])))
#     op_ax.set_xticklabels(
#         ['{}\n({:+.0f})'.format(p, mjoa_delta_by_participant[p])
#          for p in op_group_participants[op_group_name]])
#     op_ax.annotate('Δ mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
#                     xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
#
# op_ax_decomp.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')

# Legend handles - kept active, the combined plot (third plot) reuses these.
op_threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=DELTA_THRESHOLD_COLORS[name],
                                label='{:g}'.format(val))
                         for name, val in DELTA_THRESHOLDS.items()]
op_condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                         for cond, marker in CONDITION_MARKERS.items()]
op_blank_handle = Line2D([0], [0], linestyle='none', marker='None', label='')

# op_all_handles = (
#     [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + op_threshold_handles +
#     [op_blank_handle] +
#     [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + op_condition_handles
# )
# op_fig.legend(handles=op_all_handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)
#
# op_fig.suptitle('Δ % cord volume above MPS threshold (PostOp - PreOp no preload), by operation type')
# op_fig.tight_layout()
#
# op_plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot_delta_by_operation_2row.pdf')
# op_fig.savefig(op_plot_path, bbox_inches='tight')
# plt.close(op_fig)
#
# print()
# print("Summary saved: {}".format(op_summary_path))
# print("Plot saved: {}".format(op_plot_path))

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

for op_fusion_p in sorted(DELTA_FUSION_PARTICIPANTS, key=lambda p: int(p[1:])):
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

print("")
print("=" * 70)
print("PART B: Per-patient blob-size distribution, PreOp-NoPreload vs PostOp")
print("=" * 70)

# ============================================================
# Blob-finding: per (threshold, state, condition), pool blobs from every
# patient with adjacency data. (ppblob_adjacency_cache/state_per_patient
# already built above.)
# ============================================================
ppblob_records = []

for ppblob_threshold_name, ppblob_threshold_val in PPBLOB_THRESHOLDS.items():
    for (ppblob_participant, ppblob_condition, ppblob_state), ppblob_edges in ppblob_adjacency_cache.items():
        ppblob_key = (ppblob_participant, ppblob_condition, ppblob_state)
        if ppblob_key not in state_per_patient:
            continue
        ppblob_df_patient = state_per_patient[ppblob_key]
        ppblob_vol_lookup = dict(zip(ppblob_df_patient['element_label'], ppblob_df_patient['volume']))
        ppblob_exceeding = set(
            ppblob_df_patient.loc[ppblob_df_patient['mps'] >= ppblob_threshold_val, 'element_label'])
        if not ppblob_exceeding:
            continue
        ppblob_state_norm = ppblob_state.strip().lower()
        for ppblob in find_blobs(ppblob_exceeding, ppblob_edges):
            ppblob_volume = sum(ppblob_vol_lookup[lbl] for lbl in ppblob)
            ppblob_r = (3.0 * ppblob_volume / (4.0 * math.pi)) ** (1.0 / 3.0)
            ppblob_records.append({
                'participant':       'P{}'.format(ppblob_participant),
                'loading_condition': ppblob_condition,
                'state':             ppblob_state_norm,
                'threshold':         ppblob_threshold_name,
                'n_elements':        len(ppblob),
                'volume':            ppblob_volume,
                'r':                 ppblob_r,
            })

ppblob_df = pd.DataFrame(ppblob_records)
ppblob_summary_path = os.path.join(OUT_DIR, 'multipatient_blob_distribution_prepost_summary.csv')
ppblob_df.to_csv(ppblob_summary_path, index=False)
print("PrePost blob summary saved: {}".format(ppblob_summary_path))

# ============================================================
# Per-patient faceted blob-size distribution, one panel per threshold (6,
# laid out 2x3), produced as TWO complete grids - one for PreOp-NoPreload,
# one for PostOp - rather than pooling across patients or mixing states
# into one panel.
# ============================================================
pp_patients_sorted = sorted(ppblob_df['participant'].unique(),
                             key=lambda p: (preop_mjoa_by_participant.get(p, float('inf')), int(p[1:])))
pp_patient_colors = {p: plt.cm.tab10(i % 10) for i, p in enumerate(pp_patients_sorted)}

# Common x/y range across BOTH states and ALL thresholds, so every panel
# (and the two states' grids) is directly comparable.
pp_r_min = ppblob_df['r'].min()
pp_r_max = ppblob_df['r'].max()
pp_xlim = (pp_r_min * 0.9, pp_r_max * 1.1)
pp_ylim = (0, 102)

# Each state's own total cord volume per (participant, condition) - fixed
# regardless of threshold, so a curve's final height shows the real % of
# that patient's cord exceeding the threshold. PreOp-NoPreload and PostOp
# are different meshes for the same patient, so this is kept per-state
# rather than shared, the same reasoning as the removed pooled version.
pp_total_cord_vol_by_state = {}
for (pp_p, pp_c, pp_s), pp_df in state_per_patient.items():
    pp_s_norm = pp_s.strip().lower()
    pp_total_cord_vol_by_state.setdefault(pp_s_norm, {})[(pp_p, pp_c)] = pp_df['volume'].sum()

PP_LOG_MAJOR_LOCATOR = LogLocator(base=10.0)
PP_LOG_NULL_FORMATTER = NullFormatter()


def pp_plot_threshold(ax, threshold_name, state_name):
    sub = ppblob_df[(ppblob_df['threshold'] == threshold_name) & (ppblob_df['state'] == state_name)]
    pp_total_cord_vol = pp_total_cord_vol_by_state.get(state_name, {})
    for (pp_participant, pp_condition), grp in sub.groupby(['participant', 'loading_condition']):
        grp_sorted = grp.sort_values('r')
        pp_participant_num = int(pp_participant[1:])   # 'P6' -> 6, to match state_per_patient's int key
        pp_total_vol = pp_total_cord_vol.get((pp_participant_num, pp_condition), 0.0)
        if pp_total_vol <= 0:
            continue
        pp_cum_pct = list(100.0 * grp_sorted['volume'].cumsum() / pp_total_vol)
        # Extend to the shared axis edges (0% at the left, held flat at the
        # final value out to the right) so the curve doesn't stop short and
        # appear to float in the middle of the panel.
        pp_rs = [pp_xlim[0]] + list(grp_sorted['r']) + [pp_xlim[1]]
        pp_ys = [0.0] + pp_cum_pct + [pp_cum_pct[-1]]
        pp_color = pp_patient_colors[pp_participant]
        pp_linestyle = BLOB_CONDITION_LINESTYLES.get(pp_condition.strip().lower(), ':')
        ax.plot(pp_rs, pp_ys, color=pp_color, linestyle=pp_linestyle, linewidth=1.2, drawstyle='steps-post')
    ax.set_xscale('log')
    ax.set_xlim(pp_xlim)
    ax.set_ylim(pp_ylim)
    ax.xaxis.set_major_locator(PP_LOG_MAJOR_LOCATOR)
    ax.xaxis.set_minor_formatter(PP_LOG_NULL_FORMATTER)
    ax.set_title('Threshold {} = {:g}'.format(threshold_name.upper(), PPBLOB_THRESHOLDS[threshold_name]))
    ax.set_xlabel('MPS concentration effective radius r (mm)')
    ax.set_ylabel('Cumulative % of total cord volume above MPS threshold')


pp_patient_handles = [Line2D([0], [0], color=pp_patient_colors[p], linestyle='-',
                              label='{} (preop mJOA {})'.format(p, preop_mjoa_by_participant.get(p, '?')))
                       for p in pp_patients_sorted]
pp_condition_handles = [Line2D([0], [0], color='black', linestyle=ls, label=cond.capitalize())
                         for cond, ls in BLOB_CONDITION_LINESTYLES.items()]

for pp_state_name in (DELTA_TOP_STATE, DELTA_BOTTOM_STATE):
    pp_state_suffix = pp_state_name.replace('-', '')   # 'preop-nopreload' -> 'preopnopreload'

    # Six individual per-threshold plots - COMMENTED OUT, only the combined
    # grid below is part of the current reduced output set.
    # for pp_threshold_name in PPBLOB_THRESHOLDS:
    #     pp_fig, pp_ax = plt.subplots(figsize=(7, 5.5))
    #     pp_plot_threshold(pp_ax, pp_threshold_name, pp_state_name)
    #     pp_fig.suptitle(DELTA_STATE_TITLES[pp_state_name])
    #     pp_fig.legend(handles=pp_patient_handles + pp_condition_handles,
    #                   loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)
    #     pp_fig.tight_layout()
    #     pp_plot_path = os.path.join(
    #         OUT_DIR, 'multipatient_mps_plot_blob_distribution_prepost_{}_{}.pdf'.format(
    #             pp_state_suffix, pp_threshold_name))
    #     pp_fig.savefig(pp_plot_path, bbox_inches='tight')
    #     plt.close(pp_fig)
    #     print("Plot saved: {}".format(pp_plot_path))

    # --- Combined 2x3 grid (6 panels) ---
    pp_grid_fig, pp_grid_axes = plt.subplots(2, 3, figsize=(18, 10))
    for pp_ax_grid, pp_threshold_name in zip(pp_grid_axes.flat, PPBLOB_THRESHOLDS):
        pp_plot_threshold(pp_ax_grid, pp_threshold_name, pp_state_name)

    pp_grid_fig.suptitle(DELTA_STATE_TITLES[pp_state_name])
    pp_grid_fig.legend(handles=pp_patient_handles + pp_condition_handles,
                        loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)
    pp_grid_fig.tight_layout()
    pp_grid_plot_path = os.path.join(
        OUT_DIR, 'multipatient_mps_plot_blob_distribution_prepost_{}_grid.pdf'.format(pp_state_suffix))
    pp_grid_fig.savefig(pp_grid_plot_path, bbox_inches='tight')
    plt.close(pp_grid_fig)
    print("Plot saved: {}".format(pp_grid_plot_path))

# --- Aggregated per-(participant, condition, threshold, state) stats ---
ppblob_stats = ppblob_df.groupby(['participant', 'loading_condition', 'threshold', 'state']).agg(
    n_blobs=('r', 'size'),
    median_r=('r', 'median'),
    max_r=('r', 'max'),
    total_volume=('volume', 'sum'),
    max_blob_volume=('volume', 'max'),
).reset_index()
ppblob_stats['pct_volume_in_largest_blob'] = (
    100.0 * ppblob_stats['max_blob_volume'] / ppblob_stats['total_volume'])
ppblob_stats_path = os.path.join(OUT_DIR, 'multipatient_blob_distribution_prepost_stats.csv')
ppblob_stats.to_csv(ppblob_stats_path, index=False)
print("PrePost blob stats saved: {}".format(ppblob_stats_path))

# ============================================================
# Two scalar-summary plots - COMMENTED OUT, not part of the current reduced
# output set.
# ============================================================


# def make_ppblob_scalar_plot(metric_col, ylabel, title, out_suffix):
#     sub_all = ppblob_stats[ppblob_stats['participant'].isin(delta_participants)]
#     metric_max = sub_all[metric_col].max()
#     ylim = (0, metric_max * 1.1 if metric_max > 0 else 1)
#
#     fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
#     ax_by_state = {DELTA_TOP_STATE: ax_top, DELTA_BOTTOM_STATE: ax_bottom}
#
#     for state_name, ax in ax_by_state.items():
#         state_sub = sub_all[sub_all['state'].astype(str).str.strip().str.lower() == state_name]
#         for threshold_name in PPBLOB_THRESHOLDS:
#             threshold_sub = state_sub[state_sub['threshold'] == threshold_name]
#             color = PPBLOB_THRESHOLD_COLORS[threshold_name]
#             for participant_label, grp in threshold_sub.groupby('participant'):
#                 if len(grp) == 2:
#                     xp = delta_x_pos[participant_label]
#                     ax.vlines(xp, grp[metric_col].min(), grp[metric_col].max(),
#                               color=color, linewidth=1.0, alpha=0.5, zorder=2)
#             for condition_name, grp in threshold_sub.groupby('loading_condition'):
#                 marker = CONDITION_MARKERS.get(condition_name.strip().lower(), 'o')
#                 xs = [delta_x_pos[p] for p in grp['participant']]
#                 ax.scatter(xs, grp[metric_col], color=color, marker=marker, s=55, zorder=3)
#         ax.set_ylim(ylim)
#         ax.set_ylabel(ylabel)
#         ax.set_title(DELTA_STATE_TITLES[state_name])
#
#     threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=PPBLOB_THRESHOLD_COLORS[name],
#                                  label='{:g}'.format(val))
#                           for name, val in PPBLOB_THRESHOLDS.items()]
#     condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
#                           for cond, marker in CONDITION_MARKERS.items()]
#     blank_handle = Line2D([0], [0], linestyle='none', marker='None', label='')
#     all_handles = (
#         [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + threshold_handles +
#         [blank_handle] +
#         [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + condition_handles
#     )
#     fig.legend(handles=all_handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)
#
#     ax_top.tick_params(labelbottom=False)
#     ax_bottom.set_xticks(range(len(delta_participants)))
#     ax_bottom.set_xticklabels(
#         ['{}\n({:+.0f})'.format(p, mjoa_delta_by_participant[p]) for p in delta_participants])
#     ax_bottom.annotate('Δ mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
#                         xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
#     ax_bottom.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')
#     fig.suptitle(title)
#     fig.tight_layout()
#
#     plot_path = os.path.join(DIAG_DIR, 'multipatient_blob_distribution_prepost_plot_{}.pdf'.format(out_suffix))
#     fig.savefig(plot_path, bbox_inches='tight')
#     plt.close(fig)
#     print("Plot saved: {}".format(plot_path))
#
#
# make_ppblob_scalar_plot('pct_volume_in_largest_blob', '% of exceeding volume\nin largest blob',
#                          '% of exceeding volume in the largest blob, PreOp (no preload) vs PostOp',
#                          'pctlargest')
# make_ppblob_scalar_plot('n_blobs', 'Number of blobs',
#                          'Number of distinct blobs, PreOp (no preload) vs PostOp',
#                          'nblobs')

print("")
print("=" * 70)
print("PART C: Whole-cord LMM, PreOp-NoPreload vs PostOp, threshold=0.015")
print("=" * 70)

# ============================================================
# Linear mixed-effects tests (R's lme4/lmerTest via rpy2 - Kenward-Roger-
# corrected t-tests, not Satterthwaite or the asymptotic z-tests statsmodels
# gives - N=12 patients is small enough that the bias correction applied to the
# fixed-effect covariance matrix in KR could be important). 
# Patient is a random intercept.
# ============================================================
os.environ.setdefault('R_HOME', r'C:\Program Files\R\R-4.6.1')
os.environ.setdefault('R_LIBS_USER', os.path.join(os.path.expanduser('~'), 'Documents', 'R', 'win-library', '4.6'))
os.environ['PATH'] = os.path.join(os.environ['R_HOME'], 'bin', 'x64') + os.pathsep + os.environ['PATH']
# If Git's own sh.exe is on PATH (e.g. running under Git Bash), R's "CMD
# config" routes through its Unix-style config.sh, which shells out to
# 'make' - not installed here, which crashes rpy2's import. Strip Git's
# bin dirs so R falls back to its native Windows config path instead.
os.environ['PATH'] = os.pathsep.join(
    p for p in os.environ['PATH'].split(os.pathsep) if 'Git' not in p)

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

importr('lme4')
importr('lmerTest')


def r_table_to_markdown_from_df(df):
    """Renders an already-converted pandas DataFrame (from an R coefficient
    table, with a 'Term' column) as a GitHub/Obsidian-style Markdown table.
    Escapes '|' - lme4's own column name 'Pr(>|t|)' contains two of them,
    which would otherwise corrupt the table's column structure."""
    def esc(v):
        s = '{:.4g}'.format(v) if isinstance(v, float) else str(v)
        return s.replace('|', '\\|')

    cols = [esc(c) for c in df.columns]
    lines = ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join(['---'] * len(cols)) + '|']
    for _, row in df.iterrows():
        lines.append('| ' + ' | '.join(esc(v) for v in row) + ' |')
    return '\n'.join(lines)


ro.r('''
    get_coef_df <- function(model, ddf = "Kenward-Roger") {
        df <- as.data.frame(coef(summary(model, ddf = ddf)))
        df <- cbind(Term = rownames(df), df)
        rownames(df) <- NULL
        df
    }
    get_varcorr_df <- function(model) {
        as.data.frame(VarCorr(model))
    }
''')


def format_random_effects_md(varcorr_df):
    """Random-effects table + ICC (patient variance / total variance)."""
    patient_var = resid_var = None
    lines = ['| Group | Variance | Std.Dev. |', '|---|---|---|']
    for _, row in varcorr_df.iterrows():
        grp = row['grp']
        lines.append('| {} | {:.4g} | {:.4g} |'.format(grp, row['vcov'], row['sdcor']))
        if grp == 'patient':
            patient_var = row['vcov']
        elif grp == 'Residual':
            resid_var = row['vcov']
    table = '\n'.join(lines)
    if patient_var is not None and resid_var is not None and (patient_var + resid_var) > 0:
        icc = patient_var / (patient_var + resid_var)
        flag = ' (singular/near-zero)' if icc < 0.01 else ''
        icc_line = '\n\n**ICC = {:.3f}**{}'.format(icc, flag)
    else:
        icc_line = ''
    return "**Random effects**\n\n" + table + icc_line


NAVY = '#003f5c'    # raw data points
TEAL = '#58a4b0'    # box fill

# Whole-cord data has no tissue_type column, just whole-Cord mps/volume -
# mps_cache_df built above (shared with Parts A/B) is exactly this: already
# filtered to PreOp-NoPreload/PostOp, Flexion/Extension.
WHOLE_CORD_STATE_LABELS = {'preop-nopreload': 'PreOp (no preload)', 'postop': 'PostOp'}


def _fit_one_condition_state_model(df, condition, threshold):
    """pct_above ~ state + (1 | patient) on ONE condition's whole-cord data
    alone. Tries the LMM first; falls back to a paired t-test if the random
    intercept is singular. Returns (df, pre_mean, post_mean, p_value,
    model_kind, coef_df, varcorr_df_or_None)."""
    cond_df = df[df['loading_condition'] == condition]
    rows = []
    for (participant, state_norm), grp in cond_df.groupby(['participant', 'state_norm']):
        rows.append({'patient': 'P{}'.format(int(participant)), 'state': WHOLE_CORD_STATE_LABELS[state_norm],
                     'pct_above': pct_volume_above(grp, threshold)})
    long_df = pd.DataFrame(rows).dropna(subset=['pct_above'])

    wide = long_df.pivot(index='patient', columns='state', values='pct_above')
    common = wide.dropna().index
    long_df = long_df[long_df['patient'].isin(common)]

    print("--- Whole-cord data, {} (patient x state, N={}) ---".format(condition, len(common)))
    print(long_df.to_string(index=False))

    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv['wc_data'] = ro.conversion.py2rpy(long_df)

    print("")
    print("--- Whole cord ({}): pct_above ~ state + (1 | patient) ---".format(condition))
    ro.r('''
        wc_data$patient <- factor(wc_data$patient)
        wc_data$state   <- factor(wc_data$state, levels = c("PreOp (no preload)", "PostOp"))
        model_wc <- lmerTest::lmer(pct_above ~ state + (1 | patient), data = wc_data)
        print(summary(model_wc, ddf = "Kenward-Roger"))
        singular_wc <- isSingular(model_wc)
    ''')
    singular = bool(ro.r('singular_wc')[0])

    pre_vals = wide.loc[common, 'PreOp (no preload)'].values
    post_vals = wide.loc[common, 'PostOp'].values
    pre_mean, post_mean = float(pre_vals.mean()), float(post_vals.mean())

    if not singular:
        print("  ==> Model used: LMM (random intercept not singular)")
        ro.r('''
            fe_wc <- fixef(model_wc)
            pre_mean_wc <- as.numeric(fe_wc['(Intercept)'])
            post_mean_wc <- as.numeric(fe_wc['(Intercept)'] + fe_wc['statePostOp'])
            p_wc <- summary(model_wc, ddf = "Kenward-Roger")$coefficients['statePostOp', 'Pr(>|t|)']
            coef_wc <- get_coef_df(model_wc)
            varcorr_wc <- get_varcorr_df(model_wc)
        ''')
        p_val = float(ro.r('p_wc')[0])
        with localconverter(ro.default_converter + pandas2ri.converter):
            coef_df = ro.conversion.rpy2py(ro.r('coef_wc'))
            varcorr_df = ro.conversion.rpy2py(ro.r('varcorr_wc'))
        return long_df, pre_mean, post_mean, p_val, 'lmm', coef_df, varcorr_df

    print("  ==> Model used: paired t-test (LMM random intercept was singular - falling back)")
    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv['pre_vals'] = ro.FloatVector(pre_vals)
        ro.globalenv['post_vals'] = ro.FloatVector(post_vals)
    ro.r('''
        tt <- t.test(post_vals, pre_vals, paired = TRUE)
        print(tt)
        ttest_df <- data.frame(
            Term = "PostOp - PreOp (no preload)", Estimate = as.numeric(tt$estimate),
            CI_low = tt$conf.int[1], CI_high = tt$conf.int[2],
            df = tt$parameter, t = tt$statistic, `Pr(>|t|)` = tt$p.value,
            check.names = FALSE
        )
    ''')
    p_val = float(ro.r('tt$p.value')[0])
    with localconverter(ro.default_converter + pandas2ri.converter):
        coef_df = ro.conversion.rpy2py(ro.r('ttest_df'))
    return long_df, pre_mean, post_mean, p_val, 'ttest', coef_df, None


def run_wholecord_state(df, tag, threshold, ymax=None):
    """PreOp-NoPreload vs PostOp, whole cord (not tissue-split), run
    completely separately for Flexion and Extension (2 independent LMMs).
    Tries the LMM first per condition, falling back to a paired t-test if
    singular. Saves a 2-panel boxplot (-> OUT_DIR) and per-condition QQ
    plots (-> DIAG_DIR), returns the Markdown section covering both
    conditions."""
    df = df.copy()
    df['state_norm'] = df['state'].astype(str).str.strip().str.lower()
    df = df[df['state_norm'].isin(WHOLE_CORD_STATE_LABELS)]

    results = {}
    for condition in ('Flexion', 'Extension'):
        wc_df, pre_mean, post_mean, p_val, model_kind, coef_df, varcorr_df = \
            _fit_one_condition_state_model(df, condition, threshold)
        results[condition] = {'df': wc_df, 'pre_mean': pre_mean, 'post_mean': post_mean, 'p': p_val,
                               'model_kind': model_kind, 'coef_df': coef_df, 'varcorr_df': varcorr_df}

        # QQ coordinates drawn via matplotlib (not R's own pdf() device) so
        # output is reproducible under SOURCE_DATE_EPOCH - R's base graphics
        # device doesn't honor that env var and would re-stamp its own
        # creation date every run regardless. Same approach as GMvsWM's
        # save_gmvswm_qq_grid.
        qq_path = os.path.join(DIAG_DIR, 'wholecord_state_qq_{}_{}.pdf'.format(tag, condition.lower()))
        if model_kind == 'lmm':
            ro.r('wc_qq <- qqnorm(resid(model_wc), plot.it=FALSE)')
            q1, q3 = list(ro.r('as.numeric(quantile(resid(model_wc), c(0.25, 0.75)))'))
            ro.r('sw_wc <- shapiro.test(resid(model_wc))')
            sw_label = 'residuals'
            qq_title = "Whole-cord ({}) residual Q-Q".format(condition)
        else:
            # Paired t-test's actual assumption: the per-patient differences
            # are ~normal - not "model residuals" (there's no model here).
            wide_qq = wc_df.pivot(index='patient', columns='state', values='pct_above')
            diffs = (wide_qq['PostOp'] - wide_qq['PreOp (no preload)']).values
            with localconverter(ro.default_converter + pandas2ri.converter):
                ro.globalenv['wc_diffs'] = ro.FloatVector(diffs)
            ro.r('wc_qq <- qqnorm(wc_diffs, plot.it=FALSE)')
            q1, q3 = list(ro.r('as.numeric(quantile(wc_diffs, c(0.25, 0.75)))'))
            ro.r('sw_wc <- shapiro.test(wc_diffs)')
            sw_label = 'paired differences'
            qq_title = "Whole-cord ({}) paired-difference Q-Q".format(condition)

        qq_theoretical = list(ro.r('wc_qq$x'))
        qq_sample = list(ro.r('wc_qq$y'))
        nq1, nq3 = list(ro.r('qnorm(c(0.25, 0.75))'))
        qq_slope = (q3 - q1) / (nq3 - nq1)
        qq_intercept = q1 - qq_slope * nq1

        qq_fig, qq_ax = plt.subplots(figsize=(5, 5))
        qq_ax.scatter(qq_theoretical, qq_sample, color=NAVY, s=20, alpha=0.8)
        qq_line_x = [min(qq_theoretical), max(qq_theoretical)]
        qq_line_y = [qq_slope * x + qq_intercept for x in qq_line_x]
        qq_ax.plot(qq_line_x, qq_line_y, color='red')
        qq_ax.set_xlabel('Theoretical Quantiles')
        qq_ax.set_ylabel('Sample Quantiles')
        qq_ax.set_title(qq_title)
        qq_fig.tight_layout()
        qq_fig.savefig(qq_path, bbox_inches='tight')
        plt.close(qq_fig)

        results[condition]['shapiro'] = (float(ro.r('sw_wc$statistic')[0]), float(ro.r('sw_wc$p.value')[0]), sw_label)
        print("  Plot saved: {}".format(qq_path))

    if ymax is None:
        ymax = max(max(r['df']['pct_above'].max(), r['pre_mean'], r['post_mean'])
                   for r in results.values()) * 1.2

    fig, axes = plt.subplots(1, 2, figsize=(9, 5.5), sharey=True)
    for ax, condition in zip(axes, ('Flexion', 'Extension')):
        r = results[condition]
        pre_vals = r['df'][r['df']['state'] == 'PreOp (no preload)']['pct_above']
        post_vals = r['df'][r['df']['state'] == 'PostOp']['pct_above']
        box = ax.boxplot([pre_vals.values, post_vals.values], positions=[0, 1], widths=0.35,
                          showfliers=False, patch_artist=True, zorder=2)
        for patch in box['boxes']:
            patch.set_facecolor(to_rgba(TEAL, 0.4))
            patch.set_edgecolor('black')
            patch.set_linewidth(0.5)
        for part in ('whiskers', 'caps', 'medians'):
            for line in box[part]:
                line.set_color('black')
                line.set_linewidth(0.5)
        wide = r['df'].pivot(index='patient', columns='state', values='pct_above')
        for _, row in wide.iterrows():
            ax.plot([0, 1], [row['PreOp (no preload)'], row['PostOp']],
                    color='grey', linewidth=0.6, alpha=0.5, zorder=2.5)

        marker = CONDITION_MARKERS[condition.lower()]
        ax.scatter([0] * len(pre_vals), pre_vals.values, color=NAVY, s=25, alpha=0.7, marker=marker, zorder=3)
        ax.scatter([1] * len(post_vals), post_vals.values, color=NAVY, s=25, alpha=0.7, marker=marker, zorder=3)
        ax.scatter([0, 1], [r['pre_mean'], r['post_mean']], marker='d', s=80, facecolor='red',
                   edgecolor='black', linewidth=1.5, zorder=5)
        # for patient_id, row in wide.iterrows():
        #     ax.text(1.06, row['PostOp'], patient_id, fontsize=6, va='center', ha='left', color=NAVY, zorder=4)

        p_label = 'p < 0.001' if r['p'] < 0.001 else 'p = {:.3f}'.format(r['p'])
        bracket_y = 0.35 * ymax if condition == 'Extension' else (65 / 70) * ymax
        tick = (1.5 / 70) * ymax
        ax.plot([0, 0, 1, 1], [bracket_y - tick, bracket_y, bracket_y, bracket_y - tick],
                color='black', linewidth=1.2, zorder=6)
        ax.text(0.5, bracket_y + (1 / 70) * ymax, p_label, ha='center', va='bottom', fontsize=10)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['PreOp\n(no preload)', 'PostOp'])
        model_label = 'LMM' if r['model_kind'] == 'lmm' else 'paired t-test'
        # ax.set_xlabel('{}\nn={} ({})'.format(condition, r['df']['patient'].nunique(), model_label))
        ax.set_xlabel('{}'.format(condition))
        ax.set_ylim(0-(0.05*ymax), ymax*1.1)
        ax.set_xlim(-0.5, 1.3)

    axes[0].set_ylabel('% of whole cord volume above threshold ({:g})'.format(threshold))
    legend_handles = [
        Line2D([0], [0], marker='d', linestyle='', markerfacecolor='red', markeredgecolor='black',
                       label='Model Prediction'),
        Line2D([0], [0], marker=CONDITION_MARKERS['flexion'], linestyle='', markerfacecolor=NAVY,
               markeredgecolor=NAVY, alpha=0.7, label='Flexion'),
        Line2D([0], [0], marker=CONDITION_MARKERS['extension'], linestyle='', markerfacecolor=NAVY,
               markeredgecolor=NAVY, alpha=0.7, label='Extension'),
        
    ]
    axes[1].legend(handles=legend_handles, loc='upper right', frameon=False)
    fig.suptitle('Whole Cord: PreOp (no preload) vs PostOp')
    fig.tight_layout()
    box_path = os.path.join(OUT_DIR, 'wholecord_state_boxplot_{}.pdf'.format(tag))
    fig.savefig(box_path, bbox_inches='tight')
    plt.close(fig)
    print("  Plot saved: {}".format(box_path))

    

    sections = []
    for condition in ('Flexion', 'Extension'):
        r = results[condition]
        sw_stat, sw_p, sw_label = r['shapiro']
        sw_flag = ' (deviates from normality)' if sw_p < 0.05 else ''
        sw_line = "\n\n**Shapiro-Wilk ({})**: W = {:.4g}, p = {:.4g}{}".format(sw_label, sw_stat, sw_p, sw_flag)

        if r['model_kind'] == 'lmm':
            model_note = "**Model used: linear mixed-effects model** (random intercept not singular)."
            hypothesis = (
                "$$y_i = \\beta_0 + \\beta_{\\text{statePostOp}}\\,\\mathbb{1}[\\text{state}_i=\\text{PostOp}] "
                "+ u_i + \\varepsilon_i$$\n\n$$H_0:\\ \\beta_{\\text{statePostOp}} = 0$$"
            )
            extra = "\n\n{}{}".format(format_random_effects_md(r['varcorr_df']), sw_line)
        else:
            model_note = ("**Model used: paired t-test** (LMM random-intercept variance was singular - "
                           "adds nothing over a plain paired comparison).")
            hypothesis = (
                "$$H_0:\\ \\mu_{\\Delta} = 0, \\quad \\Delta_i = \\text{pct\\_above}_{i,\\text{PostOp}} "
                "- \\text{pct\\_above}_{i,\\text{PreOp}}$$"
            )
            extra = sw_line

        sections.append("""\
#### Whole cord ({condition}): PreOp (no preload) vs PostOp

{model_note}

{hypothesis}

No difference in % of whole cord volume above MPS $={threshold:g}$ between PreOp \
(no preload) and PostOp, within {condition} (not pooled with {other}).

{coef_table}{extra}
""".format(condition=condition, other='Extension' if condition == 'Flexion' else 'Flexion',
           model_note=model_note, hypothesis=hypothesis, threshold=threshold,
           coef_table=r_table_to_markdown_from_df(r['coef_df']), extra=extra))

    return "\n".join(sections)


print("")
print("=== Whole cord: PreOp (no preload) vs PostOp, threshold=0.015 ===")
MID_THRESHOLD = 0.015
wholecord_section_t015 = run_wholecord_state(mps_cache_df, 'prepost_t0p015', threshold=MID_THRESHOLD)

WHOLECORD_PREAMBLE = """\
Tries a mixed model first (patient random intercept); falls back to a \
paired t-test if that random intercept is singular - which model actually \
ran is stated explicitly under each condition below (re-decided every run, \
not assumed from a past diagnostic). Flexion/Extension kept as separate \
tests, not pooled.
"""

summary_md_wholecord_t015 = (
    "### Whole cord - PreOp (no preload) vs PostOp (threshold = {:.3f})\n\n".format(MID_THRESHOLD)
    + WHOLECORD_PREAMBLE + "\n" + wholecord_section_t015
)
summary_md_wholecord_t015_path = os.path.join(OUT_DIR, 'summary_wholecord_prepost_t0p015.md')
with open(summary_md_wholecord_t015_path, 'w', encoding='utf-8') as f:
    f.write(summary_md_wholecord_t015)
print("")
print("Summary saved: {}".format(summary_md_wholecord_t015_path))

print("")
print("=" * 70)
print("PART D: GM/WM tissue boundary enrichment, PreOp (no preload)")
print("=" * 70)

# ============================================================
# Does high-strain volume sit in GM or WM, and/or at the GM/WM boundary?
# One number per (tissue, region) combo - % of that subset's own volume
# above threshold - for GM x Interior, GM x Boundary, WM x Interior,
# WM x Boundary. Needs '_mps_GM_WM.csv' and '_topology.csv' per job (from
# Alex_results_extraction_GM_WM.py); missing files are skipped. Shared by
# Stage 2's PreOp-with-preload run below - defined once here (out_dir/
# cache_dir passed explicitly per call) rather than duplicated per stage.
# ============================================================
PREOP_THRESHOLDS = {'t0p10': 0.10, 't0p15': 0.15}


def find_boundary_elements(edges, tissue_type_by_element):
    adjacency = defaultdict(set)
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)

    boundary = set()
    for elem, neighbors in adjacency.items():
        t = tissue_type_by_element.get(elem)
        if t is None:
            continue
        for n in neighbors:
            nt = tissue_type_by_element.get(n)
            if nt is not None and nt != t:
                boundary.add(elem)
                break
    return boundary


def pct_above_for_subset(sub_df, thresholds):
    """{threshold: pct_above} - % of sub_df's own volume above each
    threshold. NaN if sub_df has zero volume."""
    total = sub_df['volume'].sum()
    per_threshold = {}
    for name, val in thresholds.items():
        if total <= 0:
            per_threshold[name] = float('nan')
            continue
        above = sub_df.loc[sub_df['mps'] >= val, 'volume'].sum()
        per_threshold[name] = 100.0 * above / total
    return per_threshold


def compute_tissue_region_pct_above(job_df, thresholds):
    """{(tissue, region): {threshold: pct_above}} for the 4 combinations of
    tissue (GM/WM) x region (Interior/Boundary) - each normalized against
    its own subset's volume, so all 4 are directly comparable."""
    is_boundary = job_df['is_boundary'].astype(bool)
    result = {}
    for tissue in ('GM', 'WM'):
        tissue_mask = job_df['tissue_type'] == tissue
        result[(tissue, 'Interior')] = pct_above_for_subset(job_df[tissue_mask & ~is_boundary], thresholds)
        result[(tissue, 'Boundary')] = pct_above_for_subset(job_df[tissue_mask & is_boundary], thresholds)
    return result


def build_job_df(csv_path):
    """Per-element mps/volume/tissue_type/is_boundary for one job, or None if
    '_mps_GM_WM.csv' is missing. is_boundary is meaningless unless
    has_topology is True."""
    gmwm_path = csv_path.replace('_mps.csv', '_mps_GM_WM.csv')
    if not os.path.isfile(gmwm_path):
        return None

    raw = pd.read_csv(gmwm_path)
    tissue_type_by_element = raw.groupby('element_label')['tissue_type'].first().to_dict()
    job_df = reduce_to_frame_mode(raw, 'peak')
    job_df['tissue_type'] = job_df['element_label'].map(tissue_type_by_element)

    topology_path = csv_path.replace('_mps.csv', '_topology.csv')
    has_topology = os.path.isfile(topology_path)
    job_df['has_topology'] = has_topology
    if has_topology:
        edges = load_adjacency_edges(topology_path)
        boundary_elements = find_boundary_elements(edges, tissue_type_by_element)
        job_df['is_boundary'] = job_df['element_label'].isin(boundary_elements)
    else:
        job_df['is_boundary'] = False

    return job_df


def make_x_order(participants):
    return sorted(participants, key=lambda p: (preop_mjoa_by_participant.get(p, float('inf')), int(p[1:])))


TISSUE_COLOR = {'GM': '#2e75b6', 'WM': '#c00000'}
REGION_FILLED = {'Interior': True, 'Boundary': False}

# Standard grouped-bar "dodge" layout: the 4 columns use GROUP_WIDTH of each
# participant's 1.0-wide slot, evenly spaced; the rest (1 - GROUP_WIDTH) is a
# guaranteed gutter to the next participant, sized as 2x the within-group gap
# so groups read as visually separate regardless of column count.
GROUP_WIDTH = 0.6
CATEGORY_ORDER = [('GM', 'Interior'), ('GM', 'Boundary'), ('WM', 'Interior'), ('WM', 'Boundary')]
_n = len(CATEGORY_ORDER)
CATEGORY_X_OFFSET = {cat: GROUP_WIDTH * (i / (_n - 1) - 0.5) for i, cat in enumerate(CATEGORY_ORDER)}


def plot_tissue_boundary(df, threshold_val, plot_path, title, has_condition):
    """One plot per threshold, 4 columns per participant: GM (blue) left,
    WM (red) right; within each, Interior (solid) left of Boundary (hollow).
    Vline connects each column's Flexion/Extension pair."""
    plot_df = df.dropna(subset=['pct_above'])
    if plot_df.empty:
        print("  Nothing to plot (zero volume for every job).")
        return

    participants = make_x_order(plot_df['participant'].unique())
    base_x = {p: i for i, p in enumerate(participants)}

    def xpos(p, tissue, region):
        return base_x[p] + CATEGORY_X_OFFSET[(tissue, region)]

    # Figure width scales with participant count so the wider spacing above
    # renders as real physical space, not just a stretched data range.
    fig, ax = plt.subplots(figsize=(max(9.5, 0.9 * len(participants) + 2), 5.5))

    if has_condition:
        for (p, tissue, region), grp in plot_df.groupby(['participant', 'tissue', 'region']):
            if len(grp) == 2:
                ax.vlines(xpos(p, tissue, region), grp['pct_above'].min(), grp['pct_above'].max(),
                          color=TISSUE_COLOR[tissue], linewidth=1.0, alpha=0.5, zorder=2)

    for (tissue, region), grp_tr in plot_df.groupby(['tissue', 'region']):
        color = TISSUE_COLOR[tissue]
        condition_groups = grp_tr.groupby('loading_condition') if has_condition else [(None, grp_tr)]
        for condition, grp in condition_groups:
            marker = CONDITION_MARKERS.get(str(condition).strip().lower(), 'o') if has_condition else 'o'
            xs = [xpos(p, tissue, region) for p in grp['participant']]
            if REGION_FILLED[region]:
                ax.scatter(xs, grp['pct_above'], color=color, marker=marker, s=55, alpha=0.85, zorder=3)
            else:
                ax.scatter(xs, grp['pct_above'], facecolors='none', edgecolors=color, marker=marker,
                           s=55, linewidths=1.4, zorder=3)

    ax.set_xticks([base_x[p] for p in participants])
    ax.set_xticklabels(['{}\n({})'.format(p, preop_mjoa_by_participant.get(p, '?')) for p in participants])
    ax.annotate('mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
    ax.set_xlabel('Participant (ordered by pre-op mJOA, ascending)')
    ax.set_ylabel('% of subset volume above threshold ({:.2f})'.format(threshold_val))
    ax.set_title(title)

    tissue_handles = [
        Line2D([0], [0], marker='o', linestyle='', color=TISSUE_COLOR['GM'],
               markerfacecolor=TISSUE_COLOR['GM'], label='GM'),
        Line2D([0], [0], marker='o', linestyle='', color=TISSUE_COLOR['WM'],
               markerfacecolor=TISSUE_COLOR['WM'], label='WM'),
    ]
    region_handles = [
        Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='black', label='Interior'),
        Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='none', label='Boundary'),
    ]
    blank = Line2D([0], [0], linestyle='none', marker='None', label='')
    handles = (
        [Line2D([0], [0], linestyle='none', marker='None', label='Tissue')] + tissue_handles +
        [blank, Line2D([0], [0], linestyle='none', marker='None', label='Region')] + region_handles
    )
    if has_condition:
        condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                              for cond, marker in CONDITION_MARKERS.items()]
        handles += [blank, Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + condition_handles
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

    fig.tight_layout()
    fig.savefig(plot_path, bbox_inches='tight')
    plt.close(fig)
    print("  Plot saved: {}".format(plot_path))


def load_or_build_combined_df(rows, tag, id_fn, cache_dir):
    """Pooled per-element data for this dataset, cached to
    'cache_gm_wm_boundary_<tag>.csv'. Delete the cache to force a rebuild."""
    id_cols = list(id_fn(rows.iloc[0]).keys()) if not rows.empty else []
    cache_path = os.path.join(cache_dir, 'cache_gm_wm_boundary_{}.csv'.format(tag))

    if os.path.isfile(cache_path):
        print("Loading cached per-element data: {}".format(cache_path))
        return pd.read_csv(cache_path), id_cols

    print("No cache yet - building {} (slow: reads raw CSVs + topology).".format(cache_path))
    parts, missing = [], []
    for _, row in rows.iterrows():
        ids = id_fn(row)
        csv_path = str(row.get('csv_path', '')).strip()
        job_df = build_job_df(csv_path) if csv_path and csv_path.lower() != 'nan' else None
        if job_df is None:
            missing.append(ids)
            continue
        for k, v in ids.items():
            job_df[k] = v
        parts.append(job_df)

    if missing:
        print("Skipping {} job(s) missing '_mps_GM_WM.csv':".format(len(missing)))
        for ids in missing:
            print("  {}".format(', '.join(str(v) for v in ids.values())))

    if not parts:
        return pd.DataFrame(), id_cols

    combined_df = pd.concat(parts, ignore_index=True)
    combined_df.to_csv(cache_path, index=False)
    print("Cached per-element data: {}".format(cache_path))
    return combined_df, id_cols


def run_dataset(rows, thresholds, tag, title_prefix, has_condition, id_fn, out_dir, cache_dir):
    combined_df, id_cols = load_or_build_combined_df(rows, tag, id_fn, cache_dir)

    records, missing_topology = [], []
    for key, job_df in (combined_df.groupby(id_cols) if not combined_df.empty else []):
        ids = dict(zip(id_cols, key if isinstance(key, tuple) else (key,)))

        if not bool(job_df['has_topology'].iloc[0]):
            missing_topology.append(ids)
            continue
        for (tissue, region), per in compute_tissue_region_pct_above(job_df, thresholds).items():
            for name in thresholds:
                records.append(dict(ids, threshold=name, tissue=tissue, region=region, pct_above=per[name]))

    if missing_topology:
        print("Skipping {} job(s) missing '_topology.csv':".format(len(missing_topology)))
        for ids in missing_topology:
            print("  {}".format(', '.join(str(v) for v in ids.values())))

    summary = pd.DataFrame(records)
    if not summary.empty:
        path = os.path.join(out_dir, 'multipatient_gm_wm_tissue_boundary_{}.csv'.format(tag))
        summary.to_csv(path, index=False)
        print(summary.to_string(index=False))
        print("Summary saved: {}".format(path))
        for name, val in thresholds.items():
            plot_path = os.path.join(
                out_dir, 'multipatient_gm_wm_tissue_boundary_{}_{}_sortedbypreopmJOA.pdf'.format(tag, name))
            plot_tissue_boundary(summary[summary['threshold'] == name], val, plot_path,
                                  '{}: high-strain volume by tissue and region (threshold {:.2f})'.format(
                                      title_prefix, val),
                                  has_condition)
    else:
        print("No data available yet.")

    return combined_df, summary


print("=== PreOp (no preload) ===")
preop_nopreload_rows = id_map[
    (id_map['State'].astype(str).str.strip().str.lower() == 'preop-nopreload') &
    (id_map['loading_condition'].astype(str).str.strip().str.lower().isin(['flexion', 'extension']))
]
preop_nopreload_elements, preop_nopreload_summary = run_dataset(
    preop_nopreload_rows, PREOP_THRESHOLDS, 'preop_nopreload', 'PreOp (no preload)',
    has_condition=True,
    id_fn=lambda row: {'participant': 'P{}'.format(int(row['participant'])),
                        'loading_condition': str(row['loading_condition']).strip()},
    out_dir=RESULTS_DIR, cache_dir=CACHE_DIR)

print("")
print("=== PostOp ===")
postop_rows = id_map[
    (id_map['State'].astype(str).str.strip().str.lower() == 'postop') &
    (id_map['loading_condition'].astype(str).str.strip().str.lower().isin(['flexion', 'extension']))
]
postop_elements, postop_summary = run_dataset(
    postop_rows, PREOP_THRESHOLDS, 'postop', 'PostOp',
    has_condition=True,
    id_fn=lambda row: {'participant': 'P{}'.format(int(row['participant'])),
                        'loading_condition': str(row['loading_condition']).strip()},
    out_dir=RESULTS_DIR, cache_dir=CACHE_DIR)

print("")
print("=" * 70)
print("PART D (cont.): GM vs WM LMM, PreOp (no preload) and PostOp, threshold=0.015")
print("=" * 70)

# ============================================================
# GMvsWM: pct_above ~ tissue + (1 | patient), GM vs WM, fit completely
# separately for Flexion and Extension (not pooled with condition as a
# covariate - the two loading modes differ too much in magnitude). Reuses
# preop_nopreload_elements (per-element data from Part F's run_dataset
# call above), pct_above_for_subset, and the rpy2/lme4/get_coef_df/
# get_varcorr_df/r_table_to_markdown_from_df/format_random_effects_md/
# NAVY/TEAL machinery already set up in Part C.
# ============================================================
PREAMBLE = """\
Patient is a random intercept; loading condition (Flexion/Extension) is kept \
as its own main-effect covariate rather than averaged away - they differ \
hugely in magnitude, so averaging would blend two different mechanical \
regimes into one number. Kenward-Roger-corrected t-tests (R `lme4`/`lmerTest`/`pbkrtest`) \
- not Satterthwaite (df-only correction) or asymptotic z - since N=12 patients is small \
enough that Kenward-Roger's extra bias-correction to the fixed-effect covariance matrix \
matters (recommended below ~30 clusters).
"""


def _fit_one_condition_tissue_model(elements_df, condition, threshold):
    """Fits pct_above ~ tissue + (1 | patient) on ONE condition's data alone.
    Returns (df, gm_mean, wm_mean, p_value, coef_df, varcorr_df, sw_stat,
    sw_p, resid_values) - residuals/Shapiro extracted here, immediately
    after fitting, since the R model object is overwritten on the next
    condition's fit."""
    cond_df = elements_df[elements_df['loading_condition'] == condition]
    rows = []
    for (participant, tissue), grp in cond_df.groupby(['participant', 'tissue_type']):
        rows.append({'patient': participant, 'tissue': tissue,
                     'pct_above': pct_above_for_subset(grp, {'t': threshold})['t']})
    df = pd.DataFrame(rows).dropna(subset=['pct_above'])

    print("--- GMvsWM data, {} (patient x tissue, N={}) ---".format(condition, df['patient'].nunique()))
    print(df.to_string(index=False))

    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv['gmvswm_data'] = ro.conversion.py2rpy(df)

    print("")
    print("--- GMvsWM ({}): pct_above ~ tissue + (1 | patient) ---".format(condition))
    ro.r('''
        gmvswm_data$patient <- factor(gmvswm_data$patient)
        gmvswm_data$tissue  <- factor(gmvswm_data$tissue, levels = c("GM", "WM"))
        gmvswm_model <- lmerTest::lmer(pct_above ~ tissue + (1 | patient), data = gmvswm_data)
        print(summary(gmvswm_model, ddf = "Kenward-Roger"))
    ''')
    ro.r('''
        gmvswm_fe <- fixef(gmvswm_model)
        gm_mean <- as.numeric(gmvswm_fe['(Intercept)'])
        wm_mean <- as.numeric(gmvswm_fe['(Intercept)'] + gmvswm_fe['tissueWM'])
        gmvswm_p <- summary(gmvswm_model, ddf = "Kenward-Roger")$coefficients['tissueWM', 'Pr(>|t|)']
        coef_df <- get_coef_df(gmvswm_model)
        varcorr_df <- get_varcorr_df(gmvswm_model)
    ''')
    gm_mean, wm_mean = ro.r('gm_mean')[0], ro.r('wm_mean')[0]
    gmvswm_p = ro.r('gmvswm_p')[0]
    with localconverter(ro.default_converter + pandas2ri.converter):
        coef_df = ro.conversion.rpy2py(ro.r('coef_df'))
        varcorr_df = ro.conversion.rpy2py(ro.r('varcorr_df'))

    ro.r('gmvswm_sw <- shapiro.test(resid(gmvswm_model))')
    sw_stat = float(ro.r('gmvswm_sw$statistic')[0])
    sw_p = float(ro.r('gmvswm_sw$p.value')[0])

    # Q-Q coordinates + qqline's reference-line slope/intercept, computed in
    # R (no scipy dependency in this codebase) for later combined plotting.
    ro.r('gmvswm_qq <- qqnorm(resid(gmvswm_model), plot.it=FALSE)')
    qq_theoretical = list(ro.r('gmvswm_qq$x'))
    qq_sample = list(ro.r('gmvswm_qq$y'))
    q1, q3 = list(ro.r('as.numeric(quantile(resid(gmvswm_model), c(0.25, 0.75)))'))
    nq1, nq3 = list(ro.r('qnorm(c(0.25, 0.75))'))
    qq_slope = (q3 - q1) / (nq3 - nq1)
    qq_intercept = q1 - qq_slope * nq1

    return (df, gm_mean, wm_mean, gmvswm_p, coef_df, varcorr_df, sw_stat, sw_p,
            qq_theoretical, qq_sample, qq_slope, qq_intercept)


def run_gmvswm_lmm(elements_df, tag, title_prefix, threshold, out_dir, ymax=None):
    """Fits GM vs WM completely separately for Flexion and Extension. Saves
    a 2-panel boxplot (Flexion left, Extension right, shared y-axis, with
    spaghetti lines connecting each patient's GM/WM pair). Residual QQ plots
    are NOT drawn here - the caller combines them across states into one
    grid (see save_gmvswm_qq_grid). Returns (markdown_section, results) -
    results carries 'resid'/'shapiro' per condition for that combined grid."""
    results = {}
    for condition in ('Flexion', 'Extension'):
        (df, gm_mean, wm_mean, p_val, coef_df, varcorr_df, sw_stat, sw_p,
         qq_theoretical, qq_sample, qq_slope, qq_intercept) = _fit_one_condition_tissue_model(
            elements_df, condition, threshold)
        results[condition] = {'df': df, 'gm_mean': gm_mean, 'wm_mean': wm_mean, 'p': p_val,
                               'coef_df': coef_df, 'varcorr_df': varcorr_df,
                               'shapiro': (sw_stat, sw_p),
                               'qq': (qq_theoretical, qq_sample, qq_slope, qq_intercept)}

    if ymax is None:
        ymax = max(max(r['df']['pct_above'].max(), r['gm_mean'], r['wm_mean']) for r in results.values()) * 1.2

    fig, axes = plt.subplots(1, 2, figsize=(9, 5.5), sharey=True)
    for ax, condition in zip(axes, ('Flexion', 'Extension')):
        r = results[condition]
        gm_vals = r['df'][r['df']['tissue'] == 'GM']['pct_above']
        wm_vals = r['df'][r['df']['tissue'] == 'WM']['pct_above']
        box = ax.boxplot([gm_vals.values, wm_vals.values], positions=[0, 1], widths=0.35,
                          showfliers=False, patch_artist=True, zorder=2)
        for patch in box['boxes']:
            patch.set_facecolor(to_rgba(TEAL, 0.4))
            patch.set_edgecolor('black')
            patch.set_linewidth(0.5)
        for part in ('whiskers', 'caps', 'medians'):
            for line in box[part]:
                line.set_color('black')
                line.set_linewidth(0.5)

        # Spaghetti lines connecting each patient's GM/WM pair - the LMM's
        # random intercept is modeling exactly this within-patient pairing.
        wide = r['df'].pivot(index='patient', columns='tissue', values='pct_above')
        for _, row in wide.iterrows():
            if pd.notna(row.get('GM')) and pd.notna(row.get('WM')):
                ax.plot([0, 1], [row['GM'], row['WM']], color='grey', linewidth=0.6, alpha=0.5, zorder=2.5)

        marker = CONDITION_MARKERS[condition.lower()]
        ax.scatter([0] * len(gm_vals), gm_vals.values, color=NAVY, s=25, alpha=0.7, marker=marker, zorder=3)
        ax.scatter([1] * len(wm_vals), wm_vals.values, color=NAVY, s=25, alpha=0.7, marker=marker, zorder=3)
        ax.scatter([0, 1], [r['gm_mean'], r['wm_mean']], marker='d', s=80, facecolor='red',
                   edgecolor='black', linewidth=1.5, zorder=5)

        p_label = 'p < 0.001' if r['p'] < 0.001 else 'p = {:.3f}'.format(r['p'])
        bracket_y = 0.3 * ymax if condition == 'Extension' else (65 / 70) * ymax
        tick = (1.5 / 70) * ymax
        ax.plot([0, 0, 1, 1], [bracket_y - tick, bracket_y, bracket_y, bracket_y - tick],
                color='black', linewidth=1.2, zorder=6)
        ax.text(0.5, bracket_y + (1 / 70) * ymax, p_label, ha='center', va='bottom', fontsize=10)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['GM', 'WM'])
        ax.set_xlabel(condition)
        ax.set_ylim(0 - (0.05 * ymax), ymax * 1.05)

    axes[0].set_ylabel('% of tissue volume above threshold ({:g})'.format(threshold))
    legend_handles = [
        Line2D([0], [0], marker='d', linestyle='', markerfacecolor='red', markeredgecolor='black',
                       label='Model Prediction'),
        Line2D([0], [0], marker=CONDITION_MARKERS['flexion'], linestyle='', markerfacecolor=NAVY,
               markeredgecolor=NAVY, alpha=0.7, label='Flexion'),
        Line2D([0], [0], marker=CONDITION_MARKERS['extension'], linestyle='', markerfacecolor=NAVY,
               markeredgecolor=NAVY, alpha=0.7, label='Extension'),
        
    ]
    axes[1].legend(handles=legend_handles, loc='upper right', frameon=False)
    fig.suptitle('{}: Grey Matter vs. White Matter'.format(title_prefix))
    fig.tight_layout()
    box_path = os.path.join(out_dir, 'lmm_gmvswm_boxplot_{}.pdf'.format(tag))
    fig.savefig(box_path, bbox_inches='tight')
    plt.close(fig)
    print("  Plot saved: {}".format(box_path))
    sys.exit()
    sections = []
    for condition in ('Flexion', 'Extension'):
        r = results[condition]
        sections.append("""\
#### GM vs WM ({condition})

$$y_i = \\beta_0 + \\beta_{{\\text{{tissueWM}}}}\\,\\mathbb{{1}}[\\text{{tissue}}_i=\\text{{WM}}] + u_i + \\varepsilon_i$$

- $y_i$: `pct_above` for patient $i$'s tissue observation (GM or WM), {condition} only
- $\\beta_0$: intercept - expected `pct_above` for GM (the reference level)
- $\\beta_{{\\text{{tissueWM}}}}$: fixed effect of WM vs GM
- $u_i \\sim \\mathcal{{N}}(0,\\tau^2)$: random intercept per patient
- $\\varepsilon_i \\sim \\mathcal{{N}}(0,\\sigma^2)$: residual error

`lmer(pct_above ~ tissue + (1 | patient))`, {condition} data only

$$H_0:\\ \\beta_{{\\text{{tissueWM}}}} = 0$$

No difference in % of cord volume above MPS $={threshold:g}$ between grey and white \
matter, within {condition} ({pct_lbl}).

{coef_table}

{random_effects}

**Shapiro-Wilk (residuals)**: W = {sw_stat:.4g}, p = {sw_p:.4g}{sw_flag}
""".format(condition=condition, threshold=threshold,
           pct_lbl='not pooled with Extension' if condition == 'Flexion' else 'not pooled with Flexion',
           coef_table=r_table_to_markdown_from_df(r['coef_df']),
           random_effects=format_random_effects_md(r['varcorr_df']),
           sw_stat=r['shapiro'][0], sw_p=r['shapiro'][1],
           sw_flag=' (residuals deviate from normality)' if r['shapiro'][1] < 0.05 else ''))

    return "\n".join(sections), results


def save_gmvswm_qq_grid(entries, grid_path):
    """entries: list of (row_label, column_label, qq_tuple, sw_stat, sw_p),
    qq_tuple = (theoretical, sample, line_slope, line_intercept) from R's
    qqnorm()/qqline() (computed in _fit_one_condition_tissue_model - no
    scipy dependency in this codebase). Grid rows = distinct row_label (in
    order of first appearance); each row's own entries are laid out into
    columns POSITIONALLY, in the order they appear for that row -
    column_label is just that panel's own title text, not matched across
    rows, so different rows may show different column labels/counts (e.g.
    different threshold sets per state). Combines what would otherwise be
    one QQ PDF per (row, column) into one figure. Name kept 'gmvswm' for
    history - also reused by Part K's oscillation LMM grid."""
    row_labels = list(dict.fromkeys(e[0] for e in entries))
    entries_by_row = {row_label: [e for e in entries if e[0] == row_label] for row_label in row_labels}
    n_cols = max(len(v) for v in entries_by_row.values())

    fig, axes = plt.subplots(len(row_labels), n_cols, figsize=(4.5 * n_cols, 4.5 * len(row_labels)), squeeze=False)
    for row_idx, row_label in enumerate(row_labels):
        row_entries = entries_by_row[row_label]
        for col_idx in range(n_cols):
            ax = axes[row_idx][col_idx]
            if col_idx >= len(row_entries):
                ax.axis('off')
                continue
            _, column_label, (qq_theoretical, qq_sample, qq_slope, qq_intercept), sw_stat, sw_p = \
                row_entries[col_idx]
            ax.scatter(qq_theoretical, qq_sample, color=NAVY, s=20, alpha=0.8)
            line_x = [min(qq_theoretical), max(qq_theoretical)]
            line_y = [qq_slope * x + qq_intercept for x in line_x]
            ax.plot(line_x, line_y, color='red')
            ax.set_xlabel('Theoretical Quantiles')
            ax.set_ylabel('Sample Quantiles')
            ax.set_title('{} ({})\nShapiro-Wilk: W={:.3g}, p={:.3g}'.format(row_label, column_label, sw_stat, sw_p))

    fig.tight_layout()
    fig.savefig(grid_path, bbox_inches='tight')
    plt.close(fig)
    print("Plot saved: {}".format(grid_path))


gmvswm_section_nopreload_t015, gmvswm_results_nopreload = run_gmvswm_lmm(
    preop_nopreload_elements, 'preop_nopreload_t0p015', 'PreOp without Preload (threshold=0.015)',
    threshold=0.015, out_dir=RESULTS_DIR)

summary_md_nopreload_t015 = (
    "### GM/WM linear mixed-effects model - PreOp without preload (threshold = 0.015)\n\n"
    + PREAMBLE + "\n" + gmvswm_section_nopreload_t015
)
summary_md_nopreload_t015_path = os.path.join(OUT_DIR, 'lmm_summary_preop_nopreload_t0p015.md')
with open(summary_md_nopreload_t015_path, 'w', encoding='utf-8') as f:
    f.write(summary_md_nopreload_t015)
print("")
print("Summary saved: {}".format(summary_md_nopreload_t015_path))

# PostOp GMvsWM LMM - COMMENTED OUT. At threshold=0.015, PostOp/GM is
# heavily zero-inflated (boxes collapse onto the axis, Extension/GM has
# effectively no spread) - the output doesn't give a meaningful comparison.
# PreOp-no-preload above is unaffected and still runs.
# gmvswm_section_postop_t015, gmvswm_results_postop = run_gmvswm_lmm(
#     postop_elements, 'postop_t0p015', 'PostOp (threshold=0.015)', threshold=0.015, out_dir=RESULTS_DIR)
#
# summary_md_postop_t015 = (
#     "### GM/WM linear mixed-effects model - PostOp (threshold = 0.015)\n\n"
#     + PREAMBLE + "\n" + gmvswm_section_postop_t015
# )
# summary_md_postop_t015_path = os.path.join(OUT_DIR, 'lmm_summary_postop_t0p015.md')
# with open(summary_md_postop_t015_path, 'w', encoding='utf-8') as f:
#     f.write(summary_md_postop_t015)
# print("")
# print("Summary saved: {}".format(summary_md_postop_t015_path))

# QQ grid - PreOp-no-preload only now (PostOp disabled above).
gmvswm_qq_entries = []
for gmvswm_row_label, gmvswm_state_results in (('PreOp (no preload)', gmvswm_results_nopreload),):
    for gmvswm_condition in ('Flexion', 'Extension'):
        gmvswm_r = gmvswm_state_results[gmvswm_condition]
        gmvswm_qq_entries.append((gmvswm_row_label, gmvswm_condition, gmvswm_r['qq'],
                                   gmvswm_r['shapiro'][0], gmvswm_r['shapiro'][1]))

gmvswm_qq_grid_path = os.path.join(DIAG_DIR, 'lmm_gmvswm_residual_qq_combined.pdf')
save_gmvswm_qq_grid(gmvswm_qq_entries, gmvswm_qq_grid_path)

print("")
print("=" * 70)
print("PART F: IVD strain (T95), PreOp-NoPreload vs PostOp, by fusion status")
print("=" * 70)

# ============================================================
# IVD strain, % IVD volume above T95, PreOp-NoPreload vs PostOp, per
# patient/condition, fusion patients highlighted. Two threshold versions:
# global (cohort-pooled) and patientwise (per-patient). Reads '_ivd_mps.csv'
# via id_map.csv's csv_path (written by Alex_results_extraction_IVD.py).
# Reuses preop_mjoa_by_participant/postop_mjoa_by_participant/
# mjoa_delta_by_participant/CONDITION_MARKERS/DELTA_FUSION_PARTICIPANTS
# already built above.
# ============================================================
IVD_STATE_PREOP = 'PreOp-NoPreload'
IVD_STATE_POSTOP = 'PostOp'
IVD_STATE_FILLED = {IVD_STATE_PREOP: True, IVD_STATE_POSTOP: False}   # solid vs hollow marker

# T90/T97/T99 commented out, not deleted - everything downstream reads this dict.
IVD_PERCENTILES = {
    # 't90': 0.90,
    't95': 0.95,
    # 't97': 0.97,
    # 't99': 0.99,
}
IVD_PERCENTILE_COLORS = {
    # 't90': '#548235',
    't95': '#2e75b6',
    # 't97': '#c00000',
    # 't99': '#7030a0',
}

# --- Cache: reduced (participant, loading_condition, state) -> [element_label, mps, volume] ---
ivd_own_cache_path = os.path.join(CACHE_DIR, 'cache_ivd_prepost_peak.csv')

if os.path.isfile(ivd_own_cache_path):
    print("Loading cached reduced data: {}".format(ivd_own_cache_path))
    ivd_cache_df = pd.read_csv(ivd_own_cache_path)
else:
    print("No cache at {} yet - building it.".format(ivd_own_cache_path))
    ivd_rows = id_map[
        (id_map['State'].astype(str).str.strip().isin([IVD_STATE_PREOP, IVD_STATE_POSTOP])) &
        (id_map['loading_condition'].astype(str).str.strip().str.lower().isin(['flexion', 'extension']))
    ]

    ivd_cache_parts = []
    ivd_missing = []
    for _, ivd_row in ivd_rows.iterrows():
        ivd_participant = int(ivd_row['participant'])
        ivd_condition = str(ivd_row['loading_condition']).strip()
        ivd_state = str(ivd_row['State']).strip()
        ivd_csv_path = str(ivd_row.get('csv_path', '')).strip()
        if not ivd_csv_path or ivd_csv_path.lower() == 'nan':
            ivd_missing.append((ivd_participant, ivd_condition, ivd_state))
            continue
        ivd_ivd_csv_path = ivd_csv_path.replace('_mps.csv', '_ivd_mps.csv')
        if not os.path.isfile(ivd_ivd_csv_path):
            ivd_missing.append((ivd_participant, ivd_condition, ivd_state))
            continue
        ivd_raw = pd.read_csv(ivd_ivd_csv_path)
        # Catches a file still being written (truncated, not a parse error) -
        # compare against id_map.csv's 'last_frame_idx' from the Cord extraction.
        ivd_expected_last_frame = ivd_row.get('last_frame_idx', None)
        if ivd_expected_last_frame not in (None, '') and not pd.isna(ivd_expected_last_frame):
            if ivd_raw['frame_index'].max() < int(ivd_expected_last_frame):
                print("  P{} ({}, {}): '_ivd_mps.csv' exists but looks incomplete "
                      "(max frame_index {} < expected {}) - still extracting, skipping for now.".format(
                          ivd_participant, ivd_condition, ivd_state,
                          ivd_raw['frame_index'].max(), int(ivd_expected_last_frame)))
                ivd_missing.append((ivd_participant, ivd_condition, ivd_state))
                continue
        ivd_reduced = reduce_to_frame_mode(ivd_raw, 'peak')
        ivd_reduced = ivd_reduced.copy()
        ivd_reduced.insert(0, 'state', ivd_state)
        ivd_reduced.insert(0, 'loading_condition', ivd_condition)
        ivd_reduced.insert(0, 'participant', ivd_participant)
        ivd_cache_parts.append(ivd_reduced)

    if ivd_missing:
        print("Skipping {} job(s) missing '_ivd_mps.csv' (IVD extraction not done yet for these):".format(
            len(ivd_missing)))
        for ivd_p, ivd_c, ivd_s in ivd_missing:
            print("  P{} ({}, {})".format(ivd_p, ivd_c, ivd_s))

    if not ivd_cache_parts:
        raise SystemExit("No IVD data loaded - check that Alex_results_extraction_IVD.py has been run.")

    ivd_cache_df = pd.concat(ivd_cache_parts, ignore_index=True)
    ivd_cache_df.to_csv(ivd_own_cache_path, index=False)
    print("Cached reduced data: {}".format(ivd_own_cache_path))

# ============================================================
# Cohort-pooled thresholds - one shared cutoff per percentile, pooled across all patients.
# ============================================================
IVD_COHORT_THRESHOLDS = {name: volume_weighted_percentile(ivd_cache_df, p=p) for name, p in IVD_PERCENTILES.items()}
print("Cohort-pooled thresholds: " + "  ".join(
    "{}={:.4f}".format(name.upper(), val) for name, val in IVD_COHORT_THRESHOLDS.items()))

# ============================================================
# Summary: % IVD volume above each cohort-pooled threshold, per job. Raw peak kept for reference.
# ============================================================
ivd_records = []
ivd_peak_records = []
for (ivd_participant, ivd_condition, ivd_state), ivd_grp in ivd_cache_df.groupby(
        ['participant', 'loading_condition', 'state']):
    ivd_p_label = 'P{}'.format(ivd_participant)
    for ivd_name, ivd_threshold_val in IVD_COHORT_THRESHOLDS.items():
        ivd_records.append({
            'participant': ivd_p_label,
            'loading_condition': ivd_condition,
            'state': ivd_state,
            'percentile': ivd_name,
            'pct_above': pct_volume_above(ivd_grp, ivd_threshold_val),
        })
    ivd_peak_records.append({
        'participant': ivd_p_label, 'loading_condition': ivd_condition, 'state': ivd_state,
        'peak_ivd_mps': ivd_grp['mps'].max(),
    })
ivd_summary = pd.DataFrame(ivd_records)
ivd_peak_summary = pd.DataFrame(ivd_peak_records)

# Only patients with both a pre-op and post-op mJOA value (needed for the x-axis ordering).
ivd_summary = ivd_summary[ivd_summary['participant'].isin(mjoa_delta_by_participant)]
ivd_peak_summary = ivd_peak_summary[ivd_peak_summary['participant'].isin(mjoa_delta_by_participant)]

ivd_summary_path = os.path.join(OUT_DIR, 'multipatient_ivd_summary_percentiles_prepost_sortedbymjoachange.csv')
ivd_summary.merge(ivd_peak_summary, on=['participant', 'loading_condition', 'state']).to_csv(
    ivd_summary_path, index=False)
print(ivd_summary.to_string(index=False))

ivd_participants = sorted(ivd_summary['participant'].unique(),
                           key=lambda p: (mjoa_delta_by_participant[p], int(p[1:])))
ivd_x_pos = {p: i for i, p in enumerate(ivd_participants)}

# ============================================================
# Plot: one figure per percentile (just T95). Marker shape = condition, fill = state.
# ============================================================
ivd_state_handles = [
    Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='black', label='PreOp (no preload)'),
    Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='none', label='PostOp'),
]
ivd_condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                          for cond, marker in CONDITION_MARKERS.items()]
ivd_fusion_handle = [Patch(facecolor='orange', alpha=0.2, label='Fusion')]
ivd_blank = Line2D([0], [0], linestyle='none', marker='None', label='')

ivd_all_handles = (
    [Line2D([0], [0], linestyle='none', marker='None', label='State')] + ivd_state_handles +
    [ivd_blank] +
    [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + ivd_condition_handles +
    [ivd_blank] + ivd_fusion_handle
)

ivd_plot_paths = []
for ivd_name, ivd_p in IVD_PERCENTILES.items():
    ivd_color = IVD_PERCENTILE_COLORS[ivd_name]
    ivd_col_df = ivd_summary[ivd_summary['percentile'] == ivd_name]

    ivd_fig, ivd_ax = plt.subplots(figsize=(9, 5.5))

    for ivd_participant_label, ivd_grp in ivd_col_df.groupby('participant'):
        for ivd_state, ivd_state_grp in ivd_grp.groupby('state'):
            if len(ivd_state_grp) == 2:
                ivd_xp = ivd_x_pos[ivd_participant_label]
                ivd_ax.vlines(ivd_xp, ivd_state_grp['pct_above'].min(), ivd_state_grp['pct_above'].max(),
                               color=ivd_color, linewidth=1.0, alpha=0.5, zorder=2)
    for (ivd_condition, ivd_state), ivd_grp in ivd_col_df.groupby(['loading_condition', 'state']):
        ivd_marker = CONDITION_MARKERS.get(ivd_condition.strip().lower(), 'o')
        ivd_filled = IVD_STATE_FILLED.get(ivd_state, True)
        ivd_xs = [ivd_x_pos[p] for p in ivd_grp['participant']]
        if ivd_filled:
            ivd_ax.scatter(ivd_xs, ivd_grp['pct_above'], color=ivd_color, marker=ivd_marker, s=60, zorder=3)
        else:
            ivd_ax.scatter(ivd_xs, ivd_grp['pct_above'], facecolors='none', edgecolors=ivd_color,
                            marker=ivd_marker, s=60, linewidths=1.4, zorder=3)

    for ivd_fusion_p in sorted(DELTA_FUSION_PARTICIPANTS, key=lambda p: int(p[1:])):
        if ivd_fusion_p in ivd_x_pos:
            ivd_xp = ivd_x_pos[ivd_fusion_p]
            ivd_ax.axvspan(ivd_xp - 0.5, ivd_xp + 0.5, color='orange', alpha=0.2, zorder=0)

    ivd_ax.set_xticks(range(len(ivd_participants)))
    ivd_ax.set_xticklabels(['{}\n({:+.0f})'.format(p, mjoa_delta_by_participant[p]) for p in ivd_participants])
    ivd_ax.annotate('Δ mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                     xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
    ivd_ax.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')
    ivd_ax.set_ylabel('% IVD volume above threshold')
    ivd_ax.set_title('IVD strain ({} = {:.4f}, cohort-pooled): PreOp (no preload) vs PostOp, by fusion status'.format(
        ivd_name.upper(), IVD_COHORT_THRESHOLDS[ivd_name]))
    ivd_ax.legend(handles=ivd_all_handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

    ivd_fig.tight_layout()

    ivd_plot_path = os.path.join(OUT_DIR, 'multipatient_ivd_plot_{}_prepost_sortedbymjoachange.pdf'.format(ivd_name))
    ivd_fig.savefig(ivd_plot_path, bbox_inches='tight')
    plt.close(ivd_fig)
    ivd_plot_paths.append(ivd_plot_path)

print()
print("Summary saved: {}".format(ivd_summary_path))
for _p in ivd_plot_paths:
    print("Plot saved: {}".format(_p))

# ============================================================
# Patient-wise version: each patient's own threshold, pooled from their own
# PreOp-NoPreload + PostOp data. Separate output files ('_patientwise' suffix).
# ============================================================
IVD_PATIENT_THRESHOLDS = {}
for ivd_participant, ivd_pgrp in ivd_cache_df.groupby('participant'):
    ivd_p_label = 'P{}'.format(ivd_participant)
    IVD_PATIENT_THRESHOLDS[ivd_p_label] = {
        name: volume_weighted_percentile(ivd_pgrp, p=p) for name, p in IVD_PERCENTILES.items()}

ivd_pw_records = []
ivd_pw_peak_records = []
for (ivd_participant, ivd_condition, ivd_state), ivd_grp in ivd_cache_df.groupby(
        ['participant', 'loading_condition', 'state']):
    ivd_p_label = 'P{}'.format(ivd_participant)
    for ivd_name, ivd_threshold_val in IVD_PATIENT_THRESHOLDS[ivd_p_label].items():
        ivd_pw_records.append({
            'participant': ivd_p_label,
            'loading_condition': ivd_condition,
            'state': ivd_state,
            'percentile': ivd_name,
            'threshold_value': ivd_threshold_val,
            'pct_above': pct_volume_above(ivd_grp, ivd_threshold_val),
        })
    ivd_pw_peak_records.append({
        'participant': ivd_p_label, 'loading_condition': ivd_condition, 'state': ivd_state,
        'peak_ivd_mps': ivd_grp['mps'].max(),
    })
ivd_pw_summary = pd.DataFrame(ivd_pw_records)
ivd_pw_peak_summary = pd.DataFrame(ivd_pw_peak_records)

ivd_pw_summary = ivd_pw_summary[ivd_pw_summary['participant'].isin(mjoa_delta_by_participant)]
ivd_pw_peak_summary = ivd_pw_peak_summary[ivd_pw_peak_summary['participant'].isin(mjoa_delta_by_participant)]

ivd_pw_summary_path = os.path.join(
    OUT_DIR, 'multipatient_ivd_summary_percentiles_prepost_sortedbymjoachange_patientwise.csv')
ivd_pw_summary.merge(ivd_pw_peak_summary, on=['participant', 'loading_condition', 'state']).to_csv(
    ivd_pw_summary_path, index=False)
print()
print(ivd_pw_summary.to_string(index=False))

ivd_pw_plot_paths = []
for ivd_name, ivd_p in IVD_PERCENTILES.items():
    ivd_color = IVD_PERCENTILE_COLORS[ivd_name]
    ivd_col_df = ivd_pw_summary[ivd_pw_summary['percentile'] == ivd_name]

    ivd_fig, ivd_ax = plt.subplots(figsize=(9, 5.5))

    for ivd_participant_label, ivd_grp in ivd_col_df.groupby('participant'):
        for ivd_state, ivd_state_grp in ivd_grp.groupby('state'):
            if len(ivd_state_grp) == 2:
                ivd_xp = ivd_x_pos[ivd_participant_label]
                ivd_ax.vlines(ivd_xp, ivd_state_grp['pct_above'].min(), ivd_state_grp['pct_above'].max(),
                               color=ivd_color, linewidth=1.0, alpha=0.5, zorder=2)
    for (ivd_condition, ivd_state), ivd_grp in ivd_col_df.groupby(['loading_condition', 'state']):
        ivd_marker = CONDITION_MARKERS.get(ivd_condition.strip().lower(), 'o')
        ivd_filled = IVD_STATE_FILLED.get(ivd_state, True)
        ivd_xs = [ivd_x_pos[p] for p in ivd_grp['participant']]
        if ivd_filled:
            ivd_ax.scatter(ivd_xs, ivd_grp['pct_above'], color=ivd_color, marker=ivd_marker, s=60, zorder=3)
        else:
            ivd_ax.scatter(ivd_xs, ivd_grp['pct_above'], facecolors='none', edgecolors=ivd_color,
                            marker=ivd_marker, s=60, linewidths=1.4, zorder=3)

    for ivd_fusion_p in sorted(DELTA_FUSION_PARTICIPANTS, key=lambda p: int(p[1:])):
        if ivd_fusion_p in ivd_x_pos:
            ivd_xp = ivd_x_pos[ivd_fusion_p]
            ivd_ax.axvspan(ivd_xp - 0.5, ivd_xp + 0.5, color='orange', alpha=0.2, zorder=0)

    # Starts at 0 explicitly - matplotlib's autoscale otherwise pads slightly below 0.
    ivd_ax.set_ylim(bottom=0)
    ivd_ax.set_xticks(range(len(ivd_participants)))
    ivd_ax.set_xticklabels(['{}\n({:+.0f})'.format(p, mjoa_delta_by_participant[p]) for p in ivd_participants])
    ivd_ax.annotate('Δ mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                     xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
    ivd_ax.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')
    ivd_ax.set_ylabel('% IVD volume above threshold')
    ivd_ax.set_title('IVD strain ({}, patient-specific threshold): PreOp (no preload) vs PostOp, '
                      'by fusion status'.format(ivd_name.upper()))
    ivd_ax.legend(handles=ivd_all_handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

    ivd_fig.tight_layout()

    ivd_pw_plot_path = os.path.join(
        OUT_DIR, 'multipatient_ivd_plot_{}_prepost_sortedbymjoachange_patientwise.pdf'.format(ivd_name))
    ivd_fig.savefig(ivd_pw_plot_path, bbox_inches='tight')
    plt.close(ivd_fig)
    ivd_pw_plot_paths.append(ivd_pw_plot_path)

# Dual-axis variant (adds each patient's own threshold value on a second
# y-axis, '_w_thresholds' suffix) - NOT part of the current output set.

print()
print("Summary saved: {}".format(ivd_pw_summary_path))
for _p in ivd_pw_plot_paths:
    print("Plot saved: {}".format(_p))

# ============================================================
# Fusion-patient PreOp vs PostOp diff, both threshold versions, as one
# Obsidian-ready markdown table (printed, not saved to a file - matches the
# source script). Flexion/Extension kept separate - never averaged.
# ============================================================
IVD_PRINT_PERCENTILE = 't95'


def _ivd_fusion_pivot(df, percentile):
    fusion_df = df[(df['participant'].isin(DELTA_FUSION_PARTICIPANTS)) & (df['percentile'] == percentile)]
    piv = fusion_df.pivot_table(index=['participant', 'loading_condition'],
                                 columns='state', values='pct_above').reset_index()
    piv['delta'] = piv[IVD_STATE_POSTOP] - piv[IVD_STATE_PREOP]
    return piv.set_index(['participant', 'loading_condition'])


ivd_pw_piv = _ivd_fusion_pivot(ivd_pw_summary, IVD_PRINT_PERCENTILE)
ivd_glob_piv = _ivd_fusion_pivot(ivd_summary, IVD_PRINT_PERCENTILE)
ivd_combined = ivd_pw_piv.join(ivd_glob_piv, lsuffix='_pw', rsuffix='_global').reset_index()
ivd_combined = ivd_combined.sort_values(['participant', 'loading_condition'])
ivd_combined['threshold_pw'] = ivd_combined['participant'].map(
    lambda p: IVD_PATIENT_THRESHOLDS[p][IVD_PRINT_PERCENTILE])
ivd_combined['threshold_global'] = IVD_COHORT_THRESHOLDS[IVD_PRINT_PERCENTILE]

print()
print("Fusion-patient PreOp vs PostOp diff - {} (Obsidian-ready markdown):".format(IVD_PRINT_PERCENTILE.upper()))
print()
print("| Participant | Condition | Threshold PS | PreOp PS (%) | PostOp PS (%) | Delta PS (pp) "
      "| Threshold Global | PreOp Global (%) | PostOp Global (%) | Delta Global (pp) |")
print("|---|---|---|---|---|---|---|---|---|---|")
for _, _r in ivd_combined.iterrows():
    print("| {} | {} | {:.4f} | {:.2f} | {:.2f} | {:+.2f} | {:.4f} | {:.2f} | {:.2f} | {:+.2f} |".format(
        _r['participant'], _r['loading_condition'], _r['threshold_pw'],
        _r[IVD_STATE_PREOP + '_pw'], _r[IVD_STATE_POSTOP + '_pw'], _r['delta_pw'],
        _r['threshold_global'],
        _r[IVD_STATE_PREOP + '_global'], _r[IVD_STATE_POSTOP + '_global'], _r['delta_global']))

# Direction tally, not a hypothesis test - n/4 patients with PostOp > PreOp.
print()
print("Direction summary - {} (fusion patients only, Obsidian-ready markdown):".format(IVD_PRINT_PERCENTILE.upper()))
print()
print("| Condition | Threshold | Increased (PostOp > PreOp) | Median delta (pp) |")
print("|---|---|---|---|")
for ivd_direction_condition, ivd_direction_grp in ivd_combined.groupby('loading_condition'):
    for ivd_delta_col, ivd_delta_label in (('delta_pw', 'Patient-specific'), ('delta_global', 'Global')):
        ivd_n_increased = int((ivd_direction_grp[ivd_delta_col] > 0).sum())
        ivd_n_total = len(ivd_direction_grp)
        ivd_median_delta = ivd_direction_grp[ivd_delta_col].median()
        print("| {} | {} | {}/{} | {:+.2f} |".format(
            ivd_direction_condition, ivd_delta_label, ivd_n_increased, ivd_n_total, ivd_median_delta))

# ============================================================
# Standalone table: just the threshold VALUES (patient-wise per fusion
# patient, plus the single global value for reference).
# ============================================================
print()
print("Threshold values - {} (Obsidian-ready markdown):".format(IVD_PRINT_PERCENTILE.upper()))
print()
print("| Participant | Threshold PS (MPS) | Threshold Global (MPS) |")
print("|---|---|---|")
for ivd_fusion_p in sorted(DELTA_FUSION_PARTICIPANTS, key=lambda p: int(p[1:])):
    print("| {} | {:.4f} | {:.4f} |".format(
        ivd_fusion_p, IVD_PATIENT_THRESHOLDS[ivd_fusion_p][IVD_PRINT_PERCENTILE], IVD_COHORT_THRESHOLDS[IVD_PRINT_PERCENTILE]))

# ============================================================
# ============================================================
# STAGE 2: PreOp (with preload)
# ============================================================
# ============================================================
STAGE2_RESULTS_DIR = os.path.join(SCRIPT_DIR, 'preop_results')
STAGE2_DIAG_DIR = os.path.join(STAGE2_RESULTS_DIR, 'diagnostic_plots')
STAGE2_CACHE_DIR = os.path.join(STAGE2_RESULTS_DIR, 'cache')
os.makedirs(STAGE2_DIAG_DIR, exist_ok=True)
os.makedirs(STAGE2_CACHE_DIR, exist_ok=True)

print("")
print("=" * 70)
print("PART G: Per-patient blob-size distribution, PreOp (with preload)")
print("=" * 70)

# ============================================================
# Per-patient faceted blob-size distribution (spatial clustering), PreOp
# WITH preload only (State == 'preop', not 'preop-nopreload'). Grid only
# (2x2, all 4 manual thresholds) - no individual per-threshold plots.
# ============================================================
MANUAL_THRESHOLDS = {'t0p05': 0.05, 't0p10': 0.10, 't0p15': 0.15, 't0p20': 0.20}

s2blob_per_patient = {}   # (participant, loading_condition) -> reduced DataFrame, PreOp (with preload) only
s2blob_missing = []
for _, s2blob_row in id_map.iterrows():
    if str(s2blob_row.get('State', '')).strip().lower() != 'preop':
        continue
    s2blob_condition = str(s2blob_row.get('loading_condition', '')).strip()
    if s2blob_condition.strip().lower() not in ('flexion', 'extension'):
        continue
    s2blob_csv_path = str(s2blob_row.get('csv_path', '')).strip()
    if not s2blob_csv_path or s2blob_csv_path.lower() == 'nan' or not os.path.isfile(s2blob_csv_path):
        s2blob_missing.append(s2blob_row)
        continue
    s2blob_raw = pd.read_csv(s2blob_csv_path)
    s2blob_per_patient[(int(s2blob_row['participant']), s2blob_condition)] = reduce_to_frame_mode(
        s2blob_raw, FRAME_MODE)

if s2blob_missing:
    print("Skipping {} row(s) with no csv_path set in id_map.csv (or file not found):".format(len(s2blob_missing)))
    for s2blob_row in s2blob_missing:
        print("  P{} ({})".format(int(s2blob_row['participant']),
                                   str(s2blob_row.get('loading_condition', '')).strip() or '?'))

if not s2blob_per_patient:
    raise SystemExit("No PreOp (with preload) patient data loaded - fill in csv_path in id_map.csv first.")

# Blob adjacency cache, PreOp (with preload) only
s2blob_adjacency_cache = {}
s2blob_missing_topology = []
for _, s2blob_row in id_map.iterrows():
    if str(s2blob_row.get('State', '')).strip().lower() != 'preop':
        continue
    s2blob_participant = int(s2blob_row['participant'])
    s2blob_condition = str(s2blob_row.get('loading_condition', '')).strip()
    s2blob_csv_path = str(s2blob_row.get('csv_path', '')).strip()
    s2blob_key = (s2blob_participant, s2blob_condition)
    if s2blob_key not in s2blob_per_patient or not s2blob_csv_path or s2blob_csv_path.lower() == 'nan':
        continue
    s2blob_topology_path = s2blob_csv_path.replace('_mps.csv', '_topology.csv')
    if not os.path.isfile(s2blob_topology_path):
        s2blob_missing_topology.append(s2blob_key)
        continue
    s2blob_adjacency_cache[s2blob_key] = load_adjacency_edges(s2blob_topology_path)

if s2blob_missing_topology:
    print("Skipping {} patient/condition(s) missing '_topology.csv' for blob analysis:".format(
        len(s2blob_missing_topology)))
    for s2blob_participant, s2blob_condition in s2blob_missing_topology:
        print("  P{} ({})".format(s2blob_participant, s2blob_condition))

# For each (threshold, condition), pool blobs from every patient with adjacency data
s2blob_records = []
for s2blob_threshold_name, s2blob_threshold_val in MANUAL_THRESHOLDS.items():
    for (s2blob_participant, s2blob_condition), s2blob_edges in s2blob_adjacency_cache.items():
        s2blob_df_patient = s2blob_per_patient[(s2blob_participant, s2blob_condition)]
        s2blob_vol_lookup = dict(zip(s2blob_df_patient['element_label'], s2blob_df_patient['volume']))
        s2blob_exceeding = set(
            s2blob_df_patient.loc[s2blob_df_patient['mps'] >= s2blob_threshold_val, 'element_label'])
        if not s2blob_exceeding:
            continue
        for s2blob in find_blobs(s2blob_exceeding, s2blob_edges):
            s2blob_volume = sum(s2blob_vol_lookup[lbl] for lbl in s2blob)
            s2blob_r = (3.0 * s2blob_volume / (4.0 * math.pi)) ** (1.0 / 3.0)
            s2blob_records.append({
                'participant':       'P{}'.format(s2blob_participant),
                'loading_condition': s2blob_condition,
                'threshold':         s2blob_threshold_name,
                'n_elements':        len(s2blob),
                'volume':            s2blob_volume,
                'r':                 s2blob_r,
            })

s2blob_df = pd.DataFrame(s2blob_records)
s2blob_summary_path = os.path.join(STAGE2_RESULTS_DIR, 'multipatient_blob_distribution_summary.csv')
s2blob_df.to_csv(s2blob_summary_path, index=False)
print("Blob summary saved: {}".format(s2blob_summary_path))

s2blob_patients_sorted = sorted(s2blob_df['participant'].unique(),
                                 key=lambda p: (preop_mjoa_by_participant.get(p, float('inf')), int(p[1:])))
s2blob_patient_colors = {p: plt.cm.tab10(i % 10) for i, p in enumerate(s2blob_patients_sorted)}

# Common x/y range across ALL thresholds, so every panel is directly comparable.
s2blob_r_min = s2blob_df['r'].min()
s2blob_r_max = s2blob_df['r'].max()
s2blob_xlim = (s2blob_r_min * 0.9, s2blob_r_max * 1.1)
s2blob_ylim = (0, 102)

# Each patient's TOTAL CORD VOLUME (not exceeding-volume subset) - fixed
# regardless of threshold, so a curve's final height shows the real % of
# that patient's cord exceeding the threshold.
s2blob_total_cord_vol = {key: df['volume'].sum() for key, df in s2blob_per_patient.items()}


def s2blob_plot_threshold(ax, threshold_name):
    sub = s2blob_df[s2blob_df['threshold'] == threshold_name]
    for (s2blob_p, s2blob_c), grp in sub.groupby(['participant', 'loading_condition']):
        grp_sorted = grp.sort_values('r')
        s2blob_p_num = int(s2blob_p[1:])   # 'P6' -> 6, to match s2blob_per_patient's int key
        s2blob_total_vol = s2blob_total_cord_vol.get((s2blob_p_num, s2blob_c), 0.0)
        if s2blob_total_vol <= 0:
            continue
        s2blob_cum_pct = list(100.0 * grp_sorted['volume'].cumsum() / s2blob_total_vol)
        s2blob_rs = [s2blob_xlim[0]] + list(grp_sorted['r']) + [s2blob_xlim[1]]
        s2blob_ys = [0.0] + s2blob_cum_pct + [s2blob_cum_pct[-1]]
        s2blob_color = s2blob_patient_colors[s2blob_p]
        s2blob_linestyle = BLOB_CONDITION_LINESTYLES.get(s2blob_c.strip().lower(), ':')
        ax.plot(s2blob_rs, s2blob_ys, color=s2blob_color, linestyle=s2blob_linestyle, linewidth=1.2,
                 drawstyle='steps-post')
    ax.set_xscale('log')
    ax.set_xlim(s2blob_xlim)
    ax.set_ylim(s2blob_ylim)
    ax.xaxis.set_major_locator(PP_LOG_MAJOR_LOCATOR)
    ax.xaxis.set_minor_formatter(PP_LOG_NULL_FORMATTER)
    ax.set_title('Threshold {} = {:.2f}'.format(threshold_name.upper(), MANUAL_THRESHOLDS[threshold_name]))
    ax.set_xlabel('MPS concentration effective radius r (mm)')
    ax.set_ylabel('Cumulative % of total cord volume above MPS threshold')


s2blob_patient_handles = [Line2D([0], [0], color=s2blob_patient_colors[p], linestyle='-',
                                  label='{} (preop mJOA {})'.format(p, preop_mjoa_by_participant.get(p, '?')))
                          for p in s2blob_patients_sorted]
s2blob_condition_handles = [Line2D([0], [0], color='black', linestyle=ls, label=cond.capitalize())
                             for cond, ls in BLOB_CONDITION_LINESTYLES.items()]

s2blob_grid_fig, s2blob_grid_axes = plt.subplots(2, 2, figsize=(13, 10))
for s2blob_ax_grid, s2blob_threshold_name in zip(s2blob_grid_axes.flat, MANUAL_THRESHOLDS):
    s2blob_plot_threshold(s2blob_ax_grid, s2blob_threshold_name)

s2blob_grid_fig.legend(handles=s2blob_patient_handles + s2blob_condition_handles,
                        loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)
s2blob_grid_fig.tight_layout()
s2blob_grid_plot_path = os.path.join(
    STAGE2_RESULTS_DIR, 'multipatient_mps_plot_blob_distribution_perpatient_grid_sortedbypreopmJOA.pdf')
s2blob_grid_fig.savefig(s2blob_grid_plot_path, bbox_inches='tight')
plt.close(s2blob_grid_fig)
print("Plot saved: {}".format(s2blob_grid_plot_path))

print("")
print("=" * 70)
print("PART H: GM/WM tissue boundary enrichment, PreOp (with preload)")
print("=" * 70)

# Reuses PREOP_THRESHOLDS/find_boundary_elements/pct_above_for_subset/
# compute_tissue_region_pct_above/build_job_df/make_x_order/TISSUE_COLOR/
# REGION_FILLED/plot_tissue_boundary/load_or_build_combined_df/run_dataset
# defined in Stage 1's Part F above - only out_dir/cache_dir differ.
print("=== PreOp (with preload) ===")
preop_rows = id_map[
    (id_map['State'].astype(str).str.strip().str.lower() == 'preop') &
    (id_map['loading_condition'].astype(str).str.strip().str.lower().isin(['flexion', 'extension']))
]
preop_elements, preop_summary = run_dataset(
    preop_rows, PREOP_THRESHOLDS, 'preop', 'PreOp (with preload)',
    has_condition=True,
    id_fn=lambda row: {'participant': 'P{}'.format(int(row['participant'])),
                        'loading_condition': str(row['loading_condition']).strip()},
    out_dir=STAGE2_RESULTS_DIR, cache_dir=STAGE2_CACHE_DIR)

print("")
print("=" * 70)
print("PART H (cont.): GM vs WM LMM, PreOp (with preload), thresholds=0.10/0.15")
print("=" * 70)

# Reuses _fit_one_condition_tissue_model/run_gmvswm_lmm/save_gmvswm_qq_grid/
# PREAMBLE defined in Stage 1's Part F above - only out_dir/diag_dir differ.
# Run at both 0.10 and 0.15 (not 0.015 as in Stage 1) - PreOp WITH preload
# has much higher-magnitude strain, so the lower Stage-1 cutoff isn't
# appropriate; 0.10/0.15 matches the tissue-boundary-enrichment thresholds
# already used for this same cohort above.
gmvswm_results_preop_by_threshold = {}
for gmvswm_threshold_name, gmvswm_threshold_val in (('t0p10', 0.10), ('t0p15', 0.15)):
    gmvswm_section_preop, gmvswm_results_preop = run_gmvswm_lmm(
        preop_elements, 'preop_{}'.format(gmvswm_threshold_name),
        'PreOp with Preload (threshold={:g})'.format(gmvswm_threshold_val),
        threshold=gmvswm_threshold_val, out_dir=STAGE2_RESULTS_DIR)
    gmvswm_results_preop_by_threshold[gmvswm_threshold_name] = gmvswm_results_preop

    summary_md_preop = (
        "### GM/WM linear mixed-effects model - PreOp with preload (threshold = {:g})\n\n".format(
            gmvswm_threshold_val)
        + PREAMBLE + "\n" + gmvswm_section_preop
    )
    summary_md_preop_path = os.path.join(STAGE2_RESULTS_DIR, 'lmm_summary_preop_{}.md'.format(gmvswm_threshold_name))
    with open(summary_md_preop_path, 'w', encoding='utf-8') as f:
        f.write(summary_md_preop)
    print("")
    print("Summary saved: {}".format(summary_md_preop_path))

# Combined QQ grid across both thresholds (rows = threshold, cols = condition).
gmvswm_qq_entries_stage2 = []
for gmvswm_threshold_name, gmvswm_threshold_label in (('t0p10', 'Threshold 0.10'), ('t0p15', 'Threshold 0.15')):
    for gmvswm_condition in ('Flexion', 'Extension'):
        gmvswm_r = gmvswm_results_preop_by_threshold[gmvswm_threshold_name][gmvswm_condition]
        gmvswm_qq_entries_stage2.append((gmvswm_threshold_label, gmvswm_condition, gmvswm_r['qq'],
                                          gmvswm_r['shapiro'][0], gmvswm_r['shapiro'][1]))

gmvswm_qq_grid_path_stage2 = os.path.join(STAGE2_DIAG_DIR, 'lmm_gmvswm_residual_qq_combined.pdf')
save_gmvswm_qq_grid(gmvswm_qq_entries_stage2, gmvswm_qq_grid_path_stage2)

print("")
print("=" * 70)
print("PART I: Effect of simulated compression-site preload (PreOp)")
print("=" * 70)

# ============================================================
# Compare PreOp WITH the simulated compression-site preload (State='preop')
# against PreOp WITHOUT it (State='preop-nopreload'), patient by patient, for
# participants with both variants - isolates the effect of the preload
# modeling assumption itself, independent of surgery. Both states are
# pre-operative, so patients are ordered by pre-op mJOA (not mJOA change,
# which needs a post-op value that doesn't apply here). No stats - purely
# descriptive, same as the source script. Reuses s2blob_per_patient (PreOp-
# with-preload, built in Part G above) and Stage 1's state_per_patient
# (PreOp-no-preload half) - no new raw-CSV reads.
# ============================================================
PRELOAD_STATE_PRELOAD = 'PreOp'
PRELOAD_STATE_NO_PRELOAD = 'PreOp-NoPreload'
PRELOAD_STATE_FILLED = {PRELOAD_STATE_PRELOAD: True, PRELOAD_STATE_NO_PRELOAD: False}   # solid vs hollow marker

PRELOAD_PLOT_A_THRESHOLDS = {'t0p10': 0.10, 't0p15': 0.15}
PRELOAD_PLOT_A_THRESHOLD_COLORS = {'t0p10': '#2e75b6', 't0p15': '#c00000'}

preload_with_preload = s2blob_per_patient   # (participant, loading_condition) -> reduced df, from Part G
preload_no_preload = {
    (p, c): df for (p, c, s), df in state_per_patient.items() if s.strip().lower() == 'preop-nopreload'
}

preload_participants_with = {p for (p, c) in preload_with_preload}
preload_participants_without = {p for (p, c) in preload_no_preload}
preload_target_participants = sorted(preload_participants_with & preload_participants_without)
if not preload_target_participants:
    raise SystemExit("No participants have both PreOp and PreOp-NoPreload reduced data.")
print("Participants with both PreOp and PreOp-NoPreload data: {}".format(
    ', '.join('P{}'.format(p) for p in preload_target_participants)))

preload_participants = sorted(('P{}'.format(p) for p in preload_target_participants),
                              key=lambda p: (preop_mjoa_by_participant.get(p, float('inf')), int(p[1:])))
preload_x_pos = {p: i for i, p in enumerate(preload_participants)}

# ============================================================
# Plot A: paired threshold comparison, PreOp (solid) vs PreOp-NoPreload (hollow)
# ============================================================
preload_records = []
for preload_state_name, preload_state_data in ((PRELOAD_STATE_PRELOAD, preload_with_preload),
                                                (PRELOAD_STATE_NO_PRELOAD, preload_no_preload)):
    for (preload_p, preload_c), preload_df in sorted(preload_state_data.items()):
        if preload_p not in preload_target_participants:
            continue
        preload_record = {'participant': 'P{}'.format(preload_p), 'loading_condition': preload_c,
                           'state': preload_state_name}
        for preload_name, preload_val in PRELOAD_PLOT_A_THRESHOLDS.items():
            preload_record['pct_above_{}'.format(preload_name)] = pct_volume_above(preload_df, preload_val)
        preload_records.append(preload_record)
preload_summary = pd.DataFrame(preload_records)

preload_summary_path = os.path.join(
    STAGE2_RESULTS_DIR, 'multipatient_mps_summary_preload_effect_sortedbypreopmJOA.csv')
preload_summary.to_csv(preload_summary_path, index=False)
print(preload_summary.to_string(index=False))

preload_fig, preload_ax = plt.subplots(figsize=(8, 5.5))

for preload_name in PRELOAD_PLOT_A_THRESHOLDS:
    preload_col = 'pct_above_{}'.format(preload_name)
    preload_color = PRELOAD_PLOT_A_THRESHOLD_COLORS[preload_name]
    # Connect Flexion<->Extension pairs within the SAME state only.
    for (_, _), preload_grp in preload_summary.groupby(['participant', 'state']):
        if len(preload_grp) == 2:
            preload_xp = preload_x_pos[preload_grp['participant'].iloc[0]]
            preload_ax.vlines(preload_xp, preload_grp[preload_col].min(), preload_grp[preload_col].max(),
                               color=preload_color, linewidth=1.0, alpha=0.5, zorder=2)
    for (preload_cond, preload_state), preload_grp in preload_summary.groupby(['loading_condition', 'state']):
        preload_marker = CONDITION_MARKERS.get(preload_cond.strip().lower(), 'o')
        preload_filled = PRELOAD_STATE_FILLED.get(preload_state, True)
        preload_xs = [preload_x_pos[p] for p in preload_grp['participant']]
        if preload_filled:
            preload_ax.scatter(preload_xs, preload_grp[preload_col], color=preload_color, marker=preload_marker,
                                s=55, zorder=3)
        else:
            preload_ax.scatter(preload_xs, preload_grp[preload_col], facecolors='none', edgecolors=preload_color,
                                marker=preload_marker, s=55, linewidths=1.3, zorder=3)

preload_threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=PRELOAD_PLOT_A_THRESHOLD_COLORS[name],
                                     label='{:.2f}'.format(val))
                              for name, val in PRELOAD_PLOT_A_THRESHOLDS.items()]
preload_condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                              for cond, marker in CONDITION_MARKERS.items()]
preload_state_handles = [
    Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='black', label='PreOp'),
    Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='none', label='PreOp-NoPreload'),
]
preload_blank = Line2D([0], [0], linestyle='none', marker='None', label='')

preload_all_handles = (
    [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + preload_threshold_handles +
    [preload_blank] +
    [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + preload_condition_handles +
    [preload_blank] +
    [Line2D([0], [0], linestyle='none', marker='None', label='State')] + preload_state_handles
)
preload_ax.legend(handles=preload_all_handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

preload_ax.set_xticks(range(len(preload_participants)))
preload_ax.set_xticklabels(
    ['{}\n({})'.format(p, preop_mjoa_by_participant.get(p, '?')) for p in preload_participants])
preload_ax.set_xlabel('Participant (ordered by pre-op mJOA, ascending)')
preload_ax.set_ylabel('% cord volume above threshold')
preload_ax.set_title('Effect of simulated compression-site preload (PreOp)')
preload_fig.tight_layout()

preload_plot_path = os.path.join(STAGE2_RESULTS_DIR, 'multipatient_mps_plot_preload_effect_sortedbypreopmJOA.pdf')
preload_fig.savefig(preload_plot_path, bbox_inches='tight')
plt.close(preload_fig)

print("Summary saved: {}".format(preload_summary_path))
print("Plot saved: {}".format(preload_plot_path))

# ============================================================
# ============================================================
# STAGE 3: Oscillation (PreOp only)
# ============================================================
# ============================================================
STAGE3_RESULTS_DIR = os.path.join(SCRIPT_DIR, 'oscillation_results')
STAGE3_DIAG_DIR = os.path.join(STAGE3_RESULTS_DIR, 'diagnostic_plots')
STAGE3_CACHE_DIR = os.path.join(STAGE3_RESULTS_DIR, 'cache')
os.makedirs(STAGE3_DIAG_DIR, exist_ok=True)
os.makedirs(STAGE3_CACHE_DIR, exist_ok=True)

print("")
print("=" * 70)
print("PART J: Oscillation, cumulative vs at-peak vs baseline (PreOp only)")
print("=" * 70)

# ============================================================
# For the CSF-pulsation oscillation (Step-3 of each id_map.csv State=='PreOp',
# loading_condition='Oscillation' job - PreOp oscillation only, NOT the
# separate PostOp-oscillation rows id_map.csv also has), compare two
# different % cord volume above threshold metrics, per patient:
#   Cumulative - each element's own PEAK mps over the whole oscillation
#   history, i.e. % of cord volume that exceeds the threshold AT SOME POINT.
#   At-peak    - for each frame, % of cord volume exceeding the threshold
#   SIMULTANEOUSLY; then the worst (max) frame across the whole history.
# Cumulative >= At-peak always; the gap shows how spatially/temporally
# spread the oscillation's strain excursions are. Thresholds 0.10/0.15
# (oscillation jobs all have the preload applied, same magnitude as Stage 2).
# Reuses id_map/preop_mjoa_by_participant/pct_volume_above/
# reduce_to_frame_mode defined above.
# ============================================================
OSC_STATE_OSCILLATION = 'Oscillation'
OSC_STATES = ('PreOp',)   # PostOp dropped entirely (zero-inflated, Shapiro-Wilk fails badly at
# every threshold tried, 0.005/0.01/0.015 - see Part K for detail); PreOp-NoPreload has no
# Oscillation rows at all in id_map either way.

OSC_METRIC_CUMULATIVE = 'cumulative'
OSC_METRIC_AT_PEAK = 'at_peak'
OSC_METRIC_FILLED = {OSC_METRIC_CUMULATIVE: True, OSC_METRIC_AT_PEAK: False}   # solid vs hollow marker

# PostOp strain magnitude is much lower than PreOp-with-preload (same
# reason Stage 1 needed 0.015 instead of 0.15 for PostOp's GM-vs-WM LMM) -
# 0.10/0.15 read as all-zero for PostOp oscillation. Thresholds are
# therefore per-state, not shared. t0p005 is cached pre-emptively (cheap -
# same raw data, one more threshold) as a ready fallback if 0.01/0.015 also
# turn out degenerate for PostOp; swap it into OSC_THRESHOLDS_BY_STATE if so.
OSC_ALL_THRESHOLDS = {'t0p005': 0.005, 't0p01': 0.01, 't0p015': 0.015, 't0p10': 0.10, 't0p15': 0.15}
OSC_THRESHOLDS_BY_STATE = {
    'PreOp': {'t0p10': 0.10, 't0p15': 0.15},
    'PostOp': {'t0p005': 0.005, 't0p01': 0.01, 't0p015': 0.015},
}
OSC_CACHE_THRESHOLD_COLORS = {
    't0p005': '#548235',  # green (fallback threshold, distinct from the others)
    't0p01': '#2a78d6',   # blue (matches Stage 1 convention)
    't0p015': '#1baf7a',  # aqua
    't0p10': '#2e75b6',   # blue (existing oscillation convention)
    't0p15': '#c00000',   # red (existing oscillation convention)
}

# --- Cache: small derived summary (per-frame + cumulative pct_above per
# participant per threshold) - not the full raw per-element-per-frame data.
# Covers BOTH PreOp and PostOp oscillation rows (id_map has both; PreOp-
# NoPreload has none - oscillation is built on top of the preload-bearing
# Step-1, so a no-preload variant doesn't exist in this dataset). ---
osc_own_cache_path = os.path.join(STAGE3_CACHE_DIR, 'cache_oscillation_cumulative_vs_peak.csv')

# Validated, not just trusted - a cache built before a threshold/state set
# changed in code would otherwise be loaded silently and produce empty
# results for whatever it's missing (exactly what happened when PostOp's
# thresholds were added after PreOp's cache already existed).
osc_required_thresholds = set().union(*(set(d.keys()) for d in OSC_THRESHOLDS_BY_STATE.values()))
osc_cache_valid = False
if os.path.isfile(osc_own_cache_path):
    print("Loading cached summary: {}".format(osc_own_cache_path))
    osc_cache_df = pd.read_csv(osc_own_cache_path)
    osc_cache_valid = (osc_required_thresholds <= set(osc_cache_df['threshold'].unique()) and
                       set(OSC_STATES) <= set(osc_cache_df['state'].unique()))
    if not osc_cache_valid:
        print("  Cache is missing required state(s)/threshold(s) (need thresholds {}, states {}) - "
              "rebuilding.".format(sorted(osc_required_thresholds), OSC_STATES))

if not osc_cache_valid:
    print("No (valid) cache at {} yet - building it.".format(osc_own_cache_path))
    osc_rows = id_map[
        (id_map['loading_condition'].astype(str).str.strip().str.lower() == OSC_STATE_OSCILLATION.lower()) &
        (id_map['State'].astype(str).str.strip().isin(OSC_STATES))
    ]
    if osc_rows.empty:
        raise SystemExit("No rows with loading_condition == '{}' and State in {} found in "
                          "id_map.csv.".format(OSC_STATE_OSCILLATION, OSC_STATES))

    osc_cache_rows = []
    osc_missing = []
    for _, osc_row in osc_rows.iterrows():
        osc_participant = int(osc_row['participant'])
        osc_state = str(osc_row['State']).strip()
        osc_csv_path = str(osc_row.get('csv_path', '')).strip()
        if not osc_csv_path or osc_csv_path.lower() == 'nan' or not os.path.isfile(osc_csv_path):
            osc_missing.append((osc_participant, osc_state))
            continue
        osc_raw = pd.read_csv(osc_csv_path)

        # At-peak: per-frame pct_above, for every frame of the oscillation history.
        for osc_frame_idx, osc_frame_df in osc_raw.groupby('frame_index'):
            osc_frame_value = osc_frame_df['frame_value'].iloc[0]
            for osc_name, osc_val in OSC_ALL_THRESHOLDS.items():
                osc_cache_rows.append({
                    'participant': osc_participant, 'state': osc_state, 'metric': OSC_METRIC_AT_PEAK,
                    'frame_index': osc_frame_idx, 'frame_value': osc_frame_value,
                    'threshold': osc_name, 'pct_above': pct_volume_above(osc_frame_df, osc_val),
                })

        # Cumulative: element-wise peak over the whole history, then one pct_above per threshold.
        osc_peak_reduced = reduce_to_frame_mode(osc_raw, 'peak')
        for osc_name, osc_val in OSC_ALL_THRESHOLDS.items():
            osc_cache_rows.append({
                'participant': osc_participant, 'state': osc_state, 'metric': OSC_METRIC_CUMULATIVE,
                'frame_index': None, 'frame_value': None,
                'threshold': osc_name, 'pct_above': pct_volume_above(osc_peak_reduced, osc_val),
            })

    if osc_missing:
        print("Skipping {} patient/state row(s) with no csv_path set (or file not found):".format(len(osc_missing)))
        for osc_p, osc_s in osc_missing:
            print("  P{} ({})".format(osc_p, osc_s))

    if not osc_cache_rows:
        raise SystemExit("No data loaded - check id_map.csv.")

    osc_cache_df = pd.DataFrame(osc_cache_rows)
    osc_cache_df.to_csv(osc_own_cache_path, index=False)
    print("Cached summary: {}".format(osc_own_cache_path))

osc_participants = sorted(('P{}'.format(p) for p in osc_cache_df['participant'].unique()),
                          key=lambda p: (preop_mjoa_by_participant.get(p, float('inf')), int(p[1:])))
osc_x_pos = {p: i for i, p in enumerate(osc_participants)}

# Alternative ordering (change in mJOA, postop - preop) for the cumulative-
# vs-at-peak plot only - same precedent as Stage 1's "sortedbymjoachange"
# plots. Restricted to participants with both pre- and post-op mJOA, same
# as that precedent.
osc_mjoa_change_participants = sorted(
    (p for p in osc_participants if p in mjoa_delta_by_participant),
    key=lambda p: (mjoa_delta_by_participant[p], int(p[1:])))
osc_mjoa_change_x_pos = {p: i for i, p in enumerate(osc_mjoa_change_participants)}

# ============================================================
# Summary + Plot 1: cumulative vs at-peak (worst frame) % cord volume above
# threshold. Plot 2: Delta (at-peak - baseline). Both run once per state
# (PreOp, PostOp), separate output files - same convention as the GMvsWM
# boxplots (combine diagnostics, keep the headline plots per-cohort).
# ============================================================
for osc_state in OSC_STATES:
    osc_state_df = osc_cache_df[osc_cache_df['state'] == osc_state]
    osc_state_tag = osc_state.lower()
    osc_display_thresholds = OSC_THRESHOLDS_BY_STATE[osc_state]
    osc_target_participants = sorted(osc_state_df['participant'].unique())
    if not osc_target_participants:
        print("No oscillation data for state '{}' - skipping.".format(osc_state))
        continue

    osc_records = []
    for osc_p in osc_target_participants:
        osc_p_label = 'P{}'.format(osc_p)
        for osc_name in osc_display_thresholds:
            osc_cum_row = osc_state_df[(osc_state_df['participant'] == osc_p) &
                                        (osc_state_df['metric'] == OSC_METRIC_CUMULATIVE) &
                                        (osc_state_df['threshold'] == osc_name)]
            osc_peak_rows = osc_state_df[(osc_state_df['participant'] == osc_p) &
                                          (osc_state_df['metric'] == OSC_METRIC_AT_PEAK) &
                                          (osc_state_df['threshold'] == osc_name)]
            if osc_cum_row.empty or osc_peak_rows.empty:
                continue
            osc_cumulative_pct = osc_cum_row['pct_above'].iloc[0]
            osc_at_peak_row = osc_peak_rows.loc[osc_peak_rows['pct_above'].idxmax()]
            osc_records.append({
                'participant': osc_p_label,
                'threshold': osc_name,
                'pct_above_cumulative': osc_cumulative_pct,
                'pct_above_at_peak': osc_at_peak_row['pct_above'],
                'at_peak_frame_index': osc_at_peak_row['frame_index'],
                'at_peak_frame_value': osc_at_peak_row['frame_value'],
                'gap_cumulative_minus_at_peak': osc_cumulative_pct - osc_at_peak_row['pct_above'],
            })
    osc_summary = pd.DataFrame(osc_records)

    osc_summary_path = os.path.join(
        STAGE3_RESULTS_DIR, 'multipatient_mps_summary_oscillation_cumulative_vs_peak_{}.csv'.format(osc_state_tag))
    osc_summary.to_csv(osc_summary_path, index=False)
    print("")
    print("=== Oscillation cumulative vs at-peak: {} ===".format(osc_state))
    print(osc_summary.to_string(index=False))

    print()
    print("Mean/median gap (cumulative - at-peak) by threshold:")
    print(osc_summary.groupby('threshold')['gap_cumulative_minus_at_peak'].agg(['mean', 'median']).to_string())

    osc_long_df = pd.concat([
        osc_summary[['participant', 'threshold', 'pct_above_cumulative']].rename(
            columns={'pct_above_cumulative': 'pct_above'}).assign(metric=OSC_METRIC_CUMULATIVE),
        osc_summary[['participant', 'threshold', 'pct_above_at_peak']].rename(
            columns={'pct_above_at_peak': 'pct_above'}).assign(metric=OSC_METRIC_AT_PEAK),
    ], ignore_index=True)

    osc_fig, osc_ax = plt.subplots(figsize=(8, 5.5))

    # Sorted by change in mJOA (not pre-op mJOA) - same precedent as Stage
    # 1's "sortedbymjoachange" plots. Restricted to participants with both
    # pre- and post-op mJOA.
    osc_cvp_long_df = osc_long_df[osc_long_df['participant'].isin(osc_mjoa_change_participants)]

    for osc_name in osc_display_thresholds:
        osc_col_df = osc_cvp_long_df[osc_cvp_long_df['threshold'] == osc_name]
        osc_color = OSC_CACHE_THRESHOLD_COLORS[osc_name]
        for osc_participant, osc_grp in osc_col_df.groupby('participant'):
            if len(osc_grp) == 2:
                osc_xp = osc_mjoa_change_x_pos[osc_participant]
                osc_ax.vlines(osc_xp, osc_grp['pct_above'].min(), osc_grp['pct_above'].max(),
                               color=osc_color, linewidth=1.0, alpha=0.5, zorder=2)
        for osc_metric, osc_grp in osc_col_df.groupby('metric'):
            osc_filled = OSC_METRIC_FILLED.get(osc_metric, True)
            osc_xs = [osc_mjoa_change_x_pos[p] for p in osc_grp['participant']]
            if osc_filled:
                osc_ax.scatter(osc_xs, osc_grp['pct_above'], color=osc_color, marker='o', s=55, zorder=3)
            else:
                osc_ax.scatter(osc_xs, osc_grp['pct_above'], facecolors='none', edgecolors=osc_color, marker='o',
                                s=55, linewidths=1.3, zorder=3)

    osc_threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=OSC_CACHE_THRESHOLD_COLORS[name],
                                     label='{:g}'.format(val))
                              for name, val in osc_display_thresholds.items()]
    osc_metric_handles = [
        Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='black',
               label='Cumulative (ever exceeds)'),
        Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='none',
               label='At-peak (worst single frame)'),
    ]
    osc_blank = Line2D([0], [0], linestyle='none', marker='None', label='')

    osc_all_handles = (
        [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + osc_threshold_handles +
        [osc_blank] +
        [Line2D([0], [0], linestyle='none', marker='None', label='Metric')] + osc_metric_handles
    )
    osc_ax.legend(handles=osc_all_handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

    osc_ax.set_xticks(range(len(osc_mjoa_change_participants)))
    osc_ax.set_xticklabels(
        ['{}\n({:+.0f})'.format(p, mjoa_delta_by_participant[p]) for p in osc_mjoa_change_participants])
    osc_ax.annotate('Δ mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                     xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
    osc_ax.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')
    osc_ax.set_ylabel('% cord volume above threshold')
    osc_ax.set_title('Oscillation ({}): cumulative (ever-exceeds) vs. at-peak (simultaneous) % cord volume'.format(
        osc_state))
    osc_fig.tight_layout()

    osc_plot_path = os.path.join(
        STAGE3_RESULTS_DIR,
        'multipatient_mps_plot_oscillation_cumulative_vs_peak_sortedbymjoachange_{}.pdf'.format(osc_state_tag))
    osc_fig.savefig(osc_plot_path, bbox_inches='tight')
    plt.close(osc_fig)

    print()
    print("Summary saved: {}".format(osc_summary_path))
    print("Plot saved: {}".format(osc_plot_path))

    # --- Plot 2: Delta (At-peak worst single frame MINUS Baseline). Baseline
    # is the pre-oscillation starting state (frame_index==0 of the same
    # Step-3, i.e. preload only, before any oscillation displacement). ---
    osc_delta_records = []
    for osc_p in osc_target_participants:
        osc_p_label = 'P{}'.format(osc_p)
        for osc_name in osc_display_thresholds:
            osc_peak_rows = osc_state_df[(osc_state_df['participant'] == osc_p) &
                                          (osc_state_df['metric'] == OSC_METRIC_AT_PEAK) &
                                          (osc_state_df['threshold'] == osc_name)]
            osc_base_row = osc_peak_rows[osc_peak_rows['frame_index'] == 0]
            if osc_peak_rows.empty or osc_base_row.empty:
                continue
            osc_at_peak_pct = osc_peak_rows['pct_above'].max()
            osc_baseline_pct = osc_base_row['pct_above'].iloc[0]
            osc_delta_records.append({
                'participant': osc_p_label,
                'threshold': osc_name,
                'pct_above_baseline': osc_baseline_pct,
                'pct_above_at_peak': osc_at_peak_pct,
                'delta_at_peak_minus_baseline': osc_at_peak_pct - osc_baseline_pct,
            })
    osc_delta_summary = pd.DataFrame(osc_delta_records)

    osc_delta_summary_path = os.path.join(
        STAGE3_RESULTS_DIR,
        'multipatient_mps_summary_oscillation_delta_at_peak_vs_baseline_{}.csv'.format(osc_state_tag))
    osc_delta_summary.to_csv(osc_delta_summary_path, index=False)
    print(osc_delta_summary.to_string(index=False))

    print()
    print("Mean/median delta (at-peak - baseline) by threshold:")
    print(osc_delta_summary.groupby('threshold')['delta_at_peak_minus_baseline'].agg(['mean', 'median']).to_string())

    osc_delta_min = osc_delta_summary['delta_at_peak_minus_baseline'].min()
    osc_delta_max = osc_delta_summary['delta_at_peak_minus_baseline'].max()
    osc_delta_pad = 0.1 * max(abs(osc_delta_min), abs(osc_delta_max), 1e-9)
    osc_delta_ylim = (osc_delta_min - osc_delta_pad, osc_delta_max + osc_delta_pad)

    osc_delta_fig, osc_delta_ax = plt.subplots(figsize=(8, 5.5))

    for osc_name in osc_display_thresholds:
        osc_col_df = osc_delta_summary[osc_delta_summary['threshold'] == osc_name]
        osc_color = OSC_CACHE_THRESHOLD_COLORS[osc_name]
        osc_xs = [osc_x_pos[p] for p in osc_col_df['participant']]
        # alpha<1 so overlapping points (common here - many patients cluster near
        # delta=0) show as visibly darker/stacked instead of hiding each other.
        osc_delta_ax.scatter(osc_xs, osc_col_df['delta_at_peak_minus_baseline'], color=osc_color, marker='o', s=55,
                              alpha=1, edgecolors='none', zorder=3)

    osc_delta_ax.axhline(0, color='gray', linestyle='--', linewidth=0.8, zorder=1)
    osc_delta_ax.set_ylim(osc_delta_ylim)
    osc_delta_ax.set_xticks(range(len(osc_participants)))
    osc_delta_ax.set_xticklabels(
        ['{}\n({})'.format(p, preop_mjoa_by_participant.get(p, '?')) for p in osc_participants])
    osc_delta_ax.annotate('mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                           xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
    osc_delta_ax.set_xlabel('Participant (ordered by pre-op mJOA, ascending)')
    osc_delta_ax.set_ylabel('Δ % cord volume above threshold\n(At-peak - Baseline)')
    osc_delta_ax.set_title(
        'Oscillation ({}): Δ % cord volume above threshold, worst single frame vs. pre-oscillation baseline'.format(
            osc_state))

    osc_delta_threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=OSC_CACHE_THRESHOLD_COLORS[name],
                                           label='{:g}'.format(val))
                                    for name, val in osc_display_thresholds.items()]
    osc_delta_ax.legend(
        handles=[Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + osc_delta_threshold_handles,
        loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

    osc_delta_fig.tight_layout()

    osc_delta_plot_path = os.path.join(
        STAGE3_RESULTS_DIR,
        'multipatient_mps_plot_oscillation_delta_at_peak_vs_baseline_sortedbypreopmJOA_{}.pdf'.format(osc_state_tag))
    osc_delta_fig.savefig(osc_delta_plot_path, bbox_inches='tight')
    plt.close(osc_delta_fig)

    print()
    print("Summary saved: {}".format(osc_delta_summary_path))
    print("Plot saved: {}".format(osc_delta_plot_path))

print("")
print("=" * 70)
print("PART K: Oscillation vs no-oscillation LMM (PreOp only)")
print("=" * 70)

# ============================================================
# Tests the effect of oscillation on MPS directly: pct_above ~ oscillation +
# (1 | patient), where "oscillation" compares the pre-oscillation baseline
# (frame 0 of Step-3, i.e. preload only, before any cyclic displacement)
# against the CUMULATIVE metric (% of cord that exceeds the threshold at
# ANY point during the step - not at-peak/worst-frame, which the plots
# above use instead). Fit per threshold (0.10, 0.15), same LMM-first/
# paired-t-test-fallback convention as Part C. Reuses get_coef_df/
# get_varcorr_df/r_table_to_markdown_from_df/format_random_effects_md/NAVY/
# TEAL/PREAMBLE/save_gmvswm_qq_grid already set up above.
#
# PostOp dropped entirely (OSC_STATES is PreOp-only above) - at 0.01/0.015
# (and 0.005), the paired differences are heavily zero-inflated (most
# patients exactly 0.0, a few nonzero) - Shapiro-Wilk fails badly (W around
# 0.55, p<0.0001) and the paired t-test result isn't trustworthy.
# ============================================================


def _fit_osc_exposure_model(osc_state_df, threshold_name):
    """Fits pct_above ~ oscillation + (1 | patient) for ONE (state,
    threshold) combination. 'oscillation' is a factor: 'No oscillation'
    (baseline, frame_index==0 of the at-peak series) vs 'Oscillation'
    (cumulative - exceeds threshold at any point in the step). Returns
    (long_df, no_osc_mean, osc_mean, p_value, model_kind, coef_df,
    varcorr_df_or_None, sw_stat, sw_p, sw_label, qq_tuple)."""
    threshold_df = osc_state_df[osc_state_df['threshold'] == threshold_name]
    baseline = threshold_df[(threshold_df['metric'] == OSC_METRIC_AT_PEAK) &
                             (threshold_df['frame_index'] == 0)].set_index('participant')['pct_above']
    cumulative = threshold_df[threshold_df['metric'] == OSC_METRIC_CUMULATIVE].set_index(
        'participant')['pct_above']
    common = sorted(set(baseline.index) & set(cumulative.index))

    long_rows = []
    for p in common:
        long_rows.append({'patient': 'P{}'.format(int(p)), 'oscillation': 'No oscillation',
                           'pct_above': baseline.loc[p]})
        long_rows.append({'patient': 'P{}'.format(int(p)), 'oscillation': 'Oscillation',
                           'pct_above': cumulative.loc[p]})
    long_df = pd.DataFrame(long_rows)

    print("--- Oscillation data, threshold {} (patient x oscillation, N={}) ---".format(threshold_name, len(common)))
    print(long_df.to_string(index=False))

    # Guard: zero variance (every value identical, e.g. all 0.0 - nothing in
    # either group ever reaches this threshold) crashes lmer()/summary() with
    # a "not a positive definite matrix" error before isSingular() can even
    # run. Not a fallback-worthy case like an ordinary singular fit - there
    # is no model to fit at all, so report it descriptively instead.
    if long_df['pct_above'].max() - long_df['pct_above'].min() < 1e-9:
        no_osc_mean = float(long_df.loc[long_df['oscillation'] == 'No oscillation', 'pct_above'].mean())
        osc_mean = float(long_df.loc[long_df['oscillation'] == 'Oscillation', 'pct_above'].mean())
        coef_df = pd.DataFrame([{'Term': 'n/a', 'Note': 'All values identical ({:g}) - zero variance, '
                                  'no model fit.'.format(no_osc_mean)}])
        return (long_df, no_osc_mean, osc_mean, None, 'degenerate', coef_df, None,
                None, None, 'n/a', None)

    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv['osc_data'] = ro.conversion.py2rpy(long_df)

    print("")
    print("--- Oscillation (threshold {}): pct_above ~ oscillation + (1 | patient) ---".format(threshold_name))
    ro.r('''
        osc_data$patient <- factor(osc_data$patient)
        osc_data$oscillation <- factor(osc_data$oscillation, levels = c("No oscillation", "Oscillation"))
        model_osc <- lmerTest::lmer(pct_above ~ oscillation + (1 | patient), data = osc_data)
        print(summary(model_osc, ddf = "Kenward-Roger"))
        singular_osc <- isSingular(model_osc)
    ''')
    singular = bool(ro.r('singular_osc')[0])

    wide = long_df.pivot(index='patient', columns='oscillation', values='pct_above')
    no_osc_vals = wide['No oscillation'].values
    osc_vals = wide['Oscillation'].values
    no_osc_mean, osc_mean = float(no_osc_vals.mean()), float(osc_vals.mean())

    if not singular:
        print("  ==> Model used: LMM (random intercept not singular)")
        ro.r('''
            fe_osc <- fixef(model_osc)
            no_osc_mean_r <- as.numeric(fe_osc['(Intercept)'])
            osc_mean_r <- as.numeric(fe_osc['(Intercept)'] + fe_osc['oscillationOscillation'])
            p_osc <- summary(model_osc, ddf = "Kenward-Roger")$coefficients['oscillationOscillation', 'Pr(>|t|)']
            coef_osc <- get_coef_df(model_osc)
            varcorr_osc <- get_varcorr_df(model_osc)
        ''')
        p_val = float(ro.r('p_osc')[0])
        with localconverter(ro.default_converter + pandas2ri.converter):
            coef_df = ro.conversion.rpy2py(ro.r('coef_osc'))
            varcorr_df = ro.conversion.rpy2py(ro.r('varcorr_osc'))
        model_kind = 'lmm'
        ro.r('osc_qq <- qqnorm(resid(model_osc), plot.it=FALSE)')
        q1, q3 = list(ro.r('as.numeric(quantile(resid(model_osc), c(0.25, 0.75)))'))
        ro.r('sw_osc <- shapiro.test(resid(model_osc))')
        sw_label = 'residuals'
    else:
        print("  ==> Model used: paired t-test (LMM random intercept was singular - falling back)")
        with localconverter(ro.default_converter + pandas2ri.converter):
            ro.globalenv['no_osc_vals'] = ro.FloatVector(no_osc_vals)
            ro.globalenv['osc_vals'] = ro.FloatVector(osc_vals)
        ro.r('''
            tt_osc <- t.test(osc_vals, no_osc_vals, paired = TRUE)
            print(tt_osc)
            ttest_osc_df <- data.frame(
                Term = "Oscillation - No oscillation", Estimate = as.numeric(tt_osc$estimate),
                CI_low = tt_osc$conf.int[1], CI_high = tt_osc$conf.int[2],
                df = tt_osc$parameter, t = tt_osc$statistic, `Pr(>|t|)` = tt_osc$p.value,
                check.names = FALSE
            )
        ''')
        p_val = float(ro.r('tt_osc$p.value')[0])
        with localconverter(ro.default_converter + pandas2ri.converter):
            coef_df = ro.conversion.rpy2py(ro.r('ttest_osc_df'))
        varcorr_df = None
        model_kind = 'ttest'
        diffs = osc_vals - no_osc_vals
        with localconverter(ro.default_converter + pandas2ri.converter):
            ro.globalenv['osc_diffs'] = ro.FloatVector(diffs)
        ro.r('osc_qq <- qqnorm(osc_diffs, plot.it=FALSE)')
        q1, q3 = list(ro.r('as.numeric(quantile(osc_diffs, c(0.25, 0.75)))'))
        ro.r('sw_osc <- shapiro.test(osc_diffs)')
        sw_label = 'paired differences'

    qq_theoretical = list(ro.r('osc_qq$x'))
    qq_sample = list(ro.r('osc_qq$y'))
    nq1, nq3 = list(ro.r('qnorm(c(0.25, 0.75))'))
    qq_slope = (q3 - q1) / (nq3 - nq1)
    qq_intercept = q1 - qq_slope * nq1
    sw_stat = float(ro.r('sw_osc$statistic')[0])
    sw_p = float(ro.r('sw_osc$p.value')[0])

    return (long_df, no_osc_mean, osc_mean, p_val, model_kind, coef_df, varcorr_df,
            sw_stat, sw_p, sw_label, (qq_theoretical, qq_sample, qq_slope, qq_intercept))


def run_osc_lmm(osc_state_df, state_label, tag, thresholds):
    """Fits pct_above ~ oscillation + (1|patient) separately for each
    threshold in `thresholds` (per-state - PreOp uses 0.10/0.15, PostOp
    uses 0.01/0.015, since PostOp's strain magnitude is much lower) - not
    pooled, consistent with every other LMM in this codebase treating
    different thresholds/conditions as separate tests. Saves an N-panel
    boxplot (one per threshold, shared y-axis, spaghetti lines per patient).
    Returns (markdown_section, results)."""
    results = {}
    for osc_threshold_name in thresholds:
        (long_df, no_osc_mean, osc_mean, p_val, model_kind, coef_df, varcorr_df,
         sw_stat, sw_p, sw_label, qq) = _fit_osc_exposure_model(osc_state_df, osc_threshold_name)
        results[osc_threshold_name] = {
            'df': long_df, 'no_osc_mean': no_osc_mean, 'osc_mean': osc_mean, 'p': p_val,
            'model_kind': model_kind, 'coef_df': coef_df, 'varcorr_df': varcorr_df,
            'shapiro': (sw_stat, sw_p, sw_label), 'qq': qq,
        }

    ymax = max(max(r['df']['pct_above'].max(), r['no_osc_mean'], r['osc_mean']) for r in results.values()) * 1.2
    ymax = ymax if ymax > 0 else 1.0

    fig, axes = plt.subplots(1, len(thresholds), figsize=(4.5 * len(thresholds), 5.5), sharey=True, squeeze=False)
    axes = axes[0]
    for ax, osc_threshold_name in zip(axes, thresholds):
        r = results[osc_threshold_name]
        no_osc_vals = r['df'][r['df']['oscillation'] == 'No oscillation']['pct_above']
        osc_vals = r['df'][r['df']['oscillation'] == 'Oscillation']['pct_above']
        box = ax.boxplot([no_osc_vals.values, osc_vals.values], positions=[0, 1], widths=0.35,
                          showfliers=False, patch_artist=True, zorder=2)
        for patch in box['boxes']:
            patch.set_facecolor(to_rgba(TEAL, 0.4))
            patch.set_edgecolor('black')
            patch.set_linewidth(0.5)
        for part in ('whiskers', 'caps', 'medians'):
            for line in box[part]:
                line.set_color('black')
                line.set_linewidth(0.5)

        wide = r['df'].pivot(index='patient', columns='oscillation', values='pct_above')
        for _, row in wide.iterrows():
            ax.plot([0, 1], [row['No oscillation'], row['Oscillation']],
                    color='grey', linewidth=0.6, alpha=0.5, zorder=2.5)

        ax.scatter([0] * len(no_osc_vals), no_osc_vals.values, color=NAVY, s=25, alpha=0.7, marker='o', zorder=3)
        ax.scatter([1] * len(osc_vals), osc_vals.values, color=NAVY, s=25, alpha=0.7, marker='o', zorder=3)
        ax.scatter([0, 1], [r['no_osc_mean'], r['osc_mean']], marker='d', s=80, facecolor='red',
                   edgecolor='black', linewidth=1.5, zorder=5)

        if r['model_kind'] == 'degenerate':
            model_label = 'no variance'
            ax.text(0.5, 0.9 * ymax, 'All values identical', ha='center', va='top', fontsize=9, style='italic')
        else:
            p_label = 'p < 0.001' if r['p'] < 0.001 else 'p = {:.3f}'.format(r['p'])
            bracket_y, tick = 0.83 * ymax, 0.02 * ymax
            ax.plot([0, 0, 1, 1], [bracket_y - tick, bracket_y, bracket_y, bracket_y - tick],
                    color='black', linewidth=1.2, zorder=6)
            ax.text(0.5, bracket_y + 0.015 * ymax, p_label, ha='center', va='bottom', fontsize=10)
            model_label = 'LMM' if r['model_kind'] == 'lmm' else 'paired t-test'

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['No\noscillation', 'Oscillation'])
        ax.set_xlabel('Threshold {:g}\nn={} ({})'.format(
            thresholds[osc_threshold_name], r['df']['patient'].nunique(), model_label))
        ax.set_ylim(0, ymax)

    axes[0].set_ylabel('% cord volume above threshold')
    mean_handle = [Line2D([0], [0], marker='d', linestyle='', markerfacecolor='red', markeredgecolor='black',
                           label='Mean / estimate')]
    axes[0].legend(handles=mean_handle, loc='upper left', frameon=False)
    fig.suptitle('Oscillation ({}): No oscillation vs Oscillation'.format(state_label))
    fig.tight_layout()
    box_path = os.path.join(STAGE3_RESULTS_DIR, 'lmm_oscillation_boxplot_{}.pdf'.format(tag))
    fig.savefig(box_path, bbox_inches='tight')
    plt.close(fig)
    print("  Plot saved: {}".format(box_path))

    sections = []
    for osc_threshold_name in thresholds:
        osc_threshold_val = thresholds[osc_threshold_name]
        r = results[osc_threshold_name]

        if r['model_kind'] == 'degenerate':
            sections.append("""\
#### Threshold {threshold:g}

All patients show % cord volume above MPS ${threshold:g}$ = {value:g} for both No oscillation \
and Oscillation, for {state} - zero variance, no model fit (nothing to test).
""".format(threshold=osc_threshold_val, state=state_label, value=r['no_osc_mean']))
            continue

        sw_stat, sw_p, sw_label = r['shapiro']
        sw_flag = ' (deviates from normality)' if sw_p < 0.05 else ''
        sw_line = "\n\n**Shapiro-Wilk ({})**: W = {:.4g}, p = {:.4g}{}".format(sw_label, sw_stat, sw_p, sw_flag)

        if r['model_kind'] == 'lmm':
            model_note = "**Model used: linear mixed-effects model** (random intercept not singular)."
            hypothesis = (
                "$$y_i = \\beta_0 + \\beta_{\\text{oscillation}}\\,"
                "\\mathbb{1}[\\text{oscillation}_i=\\text{Oscillation}] + u_i + \\varepsilon_i$$"
                "\n\n$$H_0:\\ \\beta_{\\text{oscillation}} = 0$$"
            )
            extra = "\n\n{}{}".format(format_random_effects_md(r['varcorr_df']), sw_line)
        else:
            model_note = ("**Model used: paired t-test** (LMM random-intercept variance was singular - "
                           "adds nothing over a plain paired comparison).")
            hypothesis = (
                "$$H_0:\\ \\mu_{\\Delta} = 0, \\quad \\Delta_i = \\text{pct\\_above}_{i,\\text{Oscillation}} "
                "- \\text{pct\\_above}_{i,\\text{No oscillation}}$$"
            )
            extra = sw_line

        sections.append("""\
#### Threshold {threshold:g}

{model_note}

{hypothesis}

No difference in % of cord volume that exceeds MPS $={threshold:g}$ at any point during the \
oscillation step, compared to the pre-oscillation baseline (frame 0), for {state}.

{coef_table}{extra}
""".format(threshold=osc_threshold_val, state=state_label, model_note=model_note, hypothesis=hypothesis,
           coef_table=r_table_to_markdown_from_df(r['coef_df']), extra=extra))

    return "\n".join(sections), results


osc_lmm_results_by_state = {}
osc_lmm_sections_by_state = {}
for osc_state in OSC_STATES:
    osc_state_df_lmm = osc_cache_df[osc_cache_df['state'] == osc_state]
    osc_section, osc_results = run_osc_lmm(
        osc_state_df_lmm, osc_state, osc_state.lower(), OSC_THRESHOLDS_BY_STATE[osc_state])
    osc_lmm_results_by_state[osc_state] = osc_results
    osc_lmm_sections_by_state[osc_state] = osc_section

# Tabulated at-a-glance summary, one row per (state, threshold), before the
# detailed per-threshold sections below.
osc_table_lines = [
    "| State | Threshold | No oscillation (mean) | Oscillation (mean) | Delta | Model | p-value | "
    "Shapiro-Wilk p |",
    "|---|---|---|---|---|---|---|---|",
]
for osc_state in OSC_STATES:
    for osc_threshold_name, osc_threshold_val in OSC_THRESHOLDS_BY_STATE[osc_state].items():
        r = osc_lmm_results_by_state[osc_state][osc_threshold_name]
        if r['model_kind'] == 'degenerate':
            osc_table_lines.append("| {} | {:g} | {:g} | {:g} | 0 | n/a (zero variance) | n/a | n/a |".format(
                osc_state, osc_threshold_val, r['no_osc_mean'], r['osc_mean']))
            continue
        osc_model_label = 'LMM' if r['model_kind'] == 'lmm' else 'paired t-test'
        osc_p_str = '<0.001' if r['p'] < 0.001 else '{:.3f}'.format(r['p'])
        osc_sw_p_str = '{:.3f}'.format(r['shapiro'][1])
        osc_table_lines.append("| {} | {:g} | {:.2f} | {:.2f} | {:+.2f} | {} | {} | {} |".format(
            osc_state, osc_threshold_val, r['no_osc_mean'], r['osc_mean'],
            r['osc_mean'] - r['no_osc_mean'], osc_model_label, osc_p_str, osc_sw_p_str))
osc_table_md = "\n".join(osc_table_lines)

OSC_PREAMBLE = """\
Patient is a random intercept; thresholds are fit as separate models, not \
pooled. Oscillation has no Flexion/Extension split (it is its own single \
loading mode, not crossed with condition), unlike every other LMM in this \
codebase - so there is no loading-condition covariate here. Kenward-Roger-corrected \
t-tests (R `lme4`/`lmerTest`/`pbkrtest`) - not Satterthwaite or asymptotic z - since \
N=12 patients is small enough that Kenward-Roger's extra bias-correction to the \
fixed-effect covariance matrix matters (recommended below ~30 clusters).
"""

summary_md_osc = (
    "### Oscillation vs no-oscillation - effect on MPS\n\n" + OSC_PREAMBLE + "\n" +
    osc_table_md + "\n\n" +
    "\n".join("### {}\n\n{}".format(osc_state, osc_lmm_sections_by_state[osc_state])
              for osc_state in OSC_STATES)
)
summary_md_osc_path = os.path.join(STAGE3_RESULTS_DIR, 'lmm_summary_oscillation.md')
with open(summary_md_osc_path, 'w', encoding='utf-8') as f:
    f.write(summary_md_osc)
print("")
print("Summary saved: {}".format(summary_md_osc_path))

# Combined QQ grid across both states - each state's own threshold columns
# (different per state), positional per row. Degenerate (zero-variance)
# results have no residuals/QQ to show, so they're skipped (that panel
# renders blank rather than erroring).
osc_qq_entries = []
for osc_state in OSC_STATES:
    for osc_threshold_name, osc_threshold_val in OSC_THRESHOLDS_BY_STATE[osc_state].items():
        r = osc_lmm_results_by_state[osc_state][osc_threshold_name]
        if r['qq'] is None:
            continue
        osc_qq_entries.append((osc_state, 'Threshold {:g}'.format(osc_threshold_val), r['qq'],
                               r['shapiro'][0], r['shapiro'][1]))

osc_qq_grid_path = os.path.join(STAGE3_DIAG_DIR, 'lmm_oscillation_residual_qq_combined.pdf')
save_gmvswm_qq_grid(osc_qq_entries, osc_qq_grid_path)
