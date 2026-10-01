"""
plot_prepost_blob_distribution.py - Standalone, fast plot of blob-size
distribution (spatial clustering, not just magnitude), PreOp-NoPreload vs
PostOp, at 5 low thresholds (0.01-0.05). A curve's final height reflects
magnitude (% of cord volume exceeding the threshold); its shape/shift
reflects spatial concentration (a few large blobs vs many small scattered
ones - the "chicken pox" pattern).

Fully self-contained: needs only id_map.csv plus the raw per-patient
extraction files already on disk (the '_mps.csv'/'_topology.csv' pairs
Alex_results_extraction.py writes next to each ODB). If either of its two
caches is missing, this script builds it itself:
  - cache_prepost_sortedbymjoachange_{FRAME_MODE}.csv (element_label, mps,
    volume per participant/condition/state - read from each row's raw
    '_mps.csv' and reduced via FRAME_MODE)
  - cache_prepost_blob_adjacency_{FRAME_MODE}.csv (adjacency edges per
    participant/condition/state - parsed from each row's '_topology.csv',
    threshold-independent, the slow part of blob analysis)
Both are also the exact format plot_prepost_sortedbymjoachange.py (for the
first one) and "Alex_results_multipatient_plot - compare_threshold.py"
(for both) read/write - whichever script runs first builds them, every
other script just reads them back.

Run: python plot_prepost_blob_distribution.py
First run (no caches yet) is slow - reads every PreOp-NoPreload/PostOp
job's raw per-frame CSV and parses its topology. Every run after that,
with both caches present, is fast.
"""

import math
import os
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter

from mps_common import PLOT_STYLE

# ============================================================
# USER SETTINGS - edit these and re-run for fast iteration
# ============================================================
ID_MAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'id_map.csv')
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Must match the FRAME_MODE the caches were generated with (selects the cache files).
FRAME_MODE = 'peak'

PPBLOB_THRESHOLDS = {'t0p01': 0.01, 't0p02': 0.02, 't0p03': 0.03, 't0p04': 0.04, 't0p05': 0.05}
PPBLOB_THRESHOLD_COLORS = {
    't0p01': '#2a78d6',   # blue
    't0p02': '#eb6834',   # orange
    't0p03': '#1baf7a',   # aqua
    't0p04': '#e34948',   # red
    't0p05': '#4a3aa7',   # violet
}

CONDITION_MARKERS = {'flexion': 's', 'extension': '^'}
BLOB_CONDITION_LINESTYLES = {'flexion': '-', 'extension': '--'}

DELTA_TOP_STATE = 'preop-nopreload'     # top row / blue curves: PreOp WITHOUT simulated preload
DELTA_BOTTOM_STATE = 'postop'           # bottom row / orange curves: PostOp
DELTA_STATE_TITLES = {DELTA_TOP_STATE: 'PreOp (no preload)', DELTA_BOTTOM_STATE: 'PostOp'}
# ============================================================

plt.rcParams.update(PLOT_STYLE)


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


BLOB_FACE_SHARING_MIN_NODES = 4   # >=4 shared nodes approximates shared face (all elements are C3D8)


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

# Rows this script needs data for: PreOp-NoPreload/PostOp, Flexion/Extension,
# with a populated csv_path (set by Alex_results_extraction.py).
_target_rows = id_map[
    id_map['State'].astype(str).str.strip().str.lower().isin((DELTA_TOP_STATE, DELTA_BOTTOM_STATE)) &
    id_map['loading_condition'].astype(str).str.strip().str.lower().isin(('flexion', 'extension')) &
    id_map['csv_path'].astype(str).str.strip().astype(bool)
]

# --- mps/volume cache: build if missing, else read ---
mps_cache_path = os.path.join(OUT_DIR, 'cache_prepost_sortedbymjoachange_{}.csv'.format(FRAME_MODE))
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
adjacency_cache_path = os.path.join(OUT_DIR, 'cache_prepost_blob_adjacency_{}.csv'.format(FRAME_MODE))
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
# Blob-finding: per (threshold, state, condition), pool blobs from every
# patient with adjacency data. Same logic as the removed Seventh-plot
# section in the main script, just sourced from the two caches instead of
# re-parsing topology.
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
# Seventh plot: per-patient faceted blob-size distribution, one panel per
# threshold (5, laid out 2x3 with one blank) - EXACTLY the existing
# PreOp-only grid plot's style (one curve per (patient, condition), color =
# patient via tab10, linestyle = condition), produced as TWO complete
# grids - one for PreOp-NoPreload, one for PostOp - rather than pooling
# across patients or mixing states into one panel.
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
    ax.set_title('Threshold {} = {:.2f}'.format(threshold_name.upper(), PPBLOB_THRESHOLDS[threshold_name]))
    ax.set_xlabel('MPS concentration effective radius r (mm)')
    ax.set_ylabel('Cumulative % of total cord volume above MPS threshold')


pp_patient_handles = [Line2D([0], [0], color=pp_patient_colors[p], linestyle='-',
                              label='{} (preop mJOA {})'.format(p, preop_mjoa_by_participant.get(p, '?')))
                       for p in pp_patients_sorted]
pp_condition_handles = [Line2D([0], [0], color='black', linestyle=ls, label=cond.capitalize())
                         for cond, ls in BLOB_CONDITION_LINESTYLES.items()]

for pp_state_name in (DELTA_TOP_STATE, DELTA_BOTTOM_STATE):
    pp_state_suffix = pp_state_name.replace('-', '')   # 'preop-nopreload' -> 'preopnopreload'

    # --- Five individual per-threshold plots ---
    for pp_threshold_name in PPBLOB_THRESHOLDS:
        pp_fig, pp_ax = plt.subplots(figsize=(7, 5.5))
        pp_plot_threshold(pp_ax, pp_threshold_name, pp_state_name)
        pp_fig.suptitle(DELTA_STATE_TITLES[pp_state_name])
        pp_fig.legend(handles=pp_patient_handles + pp_condition_handles,
                      loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)
        pp_fig.tight_layout()
        pp_plot_path = os.path.join(
            OUT_DIR, 'multipatient_mps_plot_blob_distribution_prepost_{}_{}.pdf'.format(
                pp_state_suffix, pp_threshold_name))
        pp_fig.savefig(pp_plot_path, bbox_inches='tight')
        plt.close(pp_fig)
        print("Plot saved: {}".format(pp_plot_path))

    # --- Combined 2x3 grid (5 panels + 1 blank) ---
    pp_grid_fig, pp_grid_axes = plt.subplots(2, 3, figsize=(18, 10))
    for pp_ax_grid, pp_threshold_name in zip(pp_grid_axes.flat, PPBLOB_THRESHOLDS):
        pp_plot_threshold(pp_ax_grid, pp_threshold_name, pp_state_name)
    pp_grid_axes.flat[-1].axis('off')   # 6th slot unused (5 thresholds, 2x3 grid)

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
# Eighth plot: two scalar-summary plots (not just the CSV), same 2-row
# layout as plot_prepost_sortedbymjoachange.py's first plot (top =
# PreOp-NoPreload, bottom = PostOp, x = participant ordered by mJOA change
# ascending) - swaps that plot's y-metric (% cord volume above threshold)
# for pct_volume_in_largest_blob and n_blobs from the aggregated stats
# table above.
# ============================================================


def make_ppblob_scalar_plot(metric_col, ylabel, title, out_suffix):
    sub_all = ppblob_stats[ppblob_stats['participant'].isin(delta_participants)]
    metric_max = sub_all[metric_col].max()
    ylim = (0, metric_max * 1.1 if metric_max > 0 else 1)

    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    ax_by_state = {DELTA_TOP_STATE: ax_top, DELTA_BOTTOM_STATE: ax_bottom}

    for state_name, ax in ax_by_state.items():
        state_sub = sub_all[sub_all['state'].astype(str).str.strip().str.lower() == state_name]
        for threshold_name in PPBLOB_THRESHOLDS:
            threshold_sub = state_sub[state_sub['threshold'] == threshold_name]
            color = PPBLOB_THRESHOLD_COLORS[threshold_name]
            for participant_label, grp in threshold_sub.groupby('participant'):
                if len(grp) == 2:
                    xp = delta_x_pos[participant_label]
                    ax.vlines(xp, grp[metric_col].min(), grp[metric_col].max(),
                              color=color, linewidth=1.0, alpha=0.5, zorder=2)
            for condition_name, grp in threshold_sub.groupby('loading_condition'):
                marker = CONDITION_MARKERS.get(condition_name.strip().lower(), 'o')
                xs = [delta_x_pos[p] for p in grp['participant']]
                ax.scatter(xs, grp[metric_col], color=color, marker=marker, s=55, zorder=3)
        ax.set_ylim(ylim)
        ax.set_ylabel(ylabel)
        ax.set_title(DELTA_STATE_TITLES[state_name])

    threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=PPBLOB_THRESHOLD_COLORS[name],
                                 label='{:.2f}'.format(val))
                          for name, val in PPBLOB_THRESHOLDS.items()]
    condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                          for cond, marker in CONDITION_MARKERS.items()]
    blank_handle = Line2D([0], [0], linestyle='none', marker='None', label='')
    all_handles = (
        [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + threshold_handles +
        [blank_handle] +
        [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + condition_handles
    )
    fig.legend(handles=all_handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)

    ax_top.tick_params(labelbottom=False)
    ax_bottom.set_xticks(range(len(delta_participants)))
    ax_bottom.set_xticklabels(
        ['{}\n({:+.0f})'.format(p, mjoa_delta_by_participant[p]) for p in delta_participants])
    ax_bottom.annotate('Δ mJOA', xy=(0, 0), xycoords=('axes fraction', 'axes fraction'),
                        xytext=(-12, -26), textcoords='offset points', ha='right', va='center')
    ax_bottom.set_xlabel('Participant (ordered by change in mJOA, postop - preop, ascending)')
    fig.suptitle(title)
    fig.tight_layout()

    plot_path = os.path.join(OUT_DIR, 'multipatient_blob_distribution_prepost_plot_{}.pdf'.format(out_suffix))
    fig.savefig(plot_path, bbox_inches='tight')
    plt.close(fig)
    print("Plot saved: {}".format(plot_path))


make_ppblob_scalar_plot('pct_volume_in_largest_blob', '% of exceeding volume\nin largest blob',
                         '% of exceeding volume in the largest blob, PreOp (no preload) vs PostOp',
                         'pctlargest')
make_ppblob_scalar_plot('n_blobs', 'Number of blobs',
                         'Number of distinct blobs, PreOp (no preload) vs PostOp',
                         'nblobs')
