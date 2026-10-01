"""
plot_prepost_blob_distribution.py - Standalone, fast plot of blob-size
distribution (spatial clustering, not just magnitude), PreOp-NoPreload vs
PostOp, at 5 low thresholds (0.01-0.05). A curve's final height reflects
magnitude (% of cord volume exceeding the threshold); its shape/shift
reflects spatial concentration (a few large blobs vs many small scattered
ones - the "chicken pox" pattern).

Reads ONLY id_map.csv (for mJOA) and two small cache CSVs written by
"Alex_results_multipatient_plot - compare_threshold.py" - never the raw
per-frame extraction CSVs or the '_topology.csv' connectivity files - so
editing PPBLOB_THRESHOLDS/colors below and re-running takes seconds:
  - cache_prepost_sortedbymjoachange_{FRAME_MODE}.csv (element_label, mps,
    volume per participant/condition/state - also used by
    plot_prepost_sortedbymjoachange.py)
  - cache_prepost_blob_adjacency_{FRAME_MODE}.csv (adjacency edges per
    participant/condition/state - threshold-independent, the expensive part
    of blob analysis, parsed from topology ONCE in the main script)

Run: python plot_prepost_blob_distribution.py
Requires both caches to exist first - run
"Alex_results_multipatient_plot - compare_threshold.py" at least once (or
after the underlying per-patient data changes) to (re)generate them.
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
PPBLOB_STATE_COLORS = {DELTA_TOP_STATE: '#2a78d6', DELTA_BOTTOM_STATE: '#eb6834'}   # blue / orange
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


# --- Load caches ---
mps_cache_path = os.path.join(OUT_DIR, 'cache_prepost_sortedbymjoachange_{}.csv'.format(FRAME_MODE))
adjacency_cache_path = os.path.join(OUT_DIR, 'cache_prepost_blob_adjacency_{}.csv'.format(FRAME_MODE))
for _cache_path in (mps_cache_path, adjacency_cache_path):
    if not os.path.isfile(_cache_path):
        raise SystemExit(
            "No cache found at {}. Run "
            "\"Alex_results_multipatient_plot - compare_threshold.py\" once first "
            "(with FRAME_MODE='{}') to generate it.".format(_cache_path, FRAME_MODE))

mps_cache_df = pd.read_csv(mps_cache_path)
state_per_patient = {
    (int(p), c, s): grp[['element_label', 'mps', 'volume']].reset_index(drop=True)
    for (p, c, s), grp in mps_cache_df.groupby(['participant', 'loading_condition', 'state'])
}

adjacency_cache_df = pd.read_csv(adjacency_cache_path)
ppblob_adjacency_cache = {
    (int(p), c, s): list(zip(grp['elem_a'], grp['elem_b']))
    for (p, c, s), grp in adjacency_cache_df.groupby(['participant', 'loading_condition', 'state'])
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
ppblob_pooled = defaultdict(list)   # (threshold_name, state, condition) -> list of (r, volume)

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
            ppblob_pooled[(ppblob_threshold_name, ppblob_state_norm, ppblob_condition)].append(
                (ppblob_r, ppblob_volume))
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

# Each state's own total cord volume per condition (not shared across
# states - PreOp-NoPreload and PostOp are different meshes for the same
# patient, so a shared denominator would let mesh differences masquerade
# as a magnitude change).
ppblob_cord_vol_by_state_condition = defaultdict(float)
for (ppblob_p, ppblob_c, ppblob_s) in ppblob_adjacency_cache:
    if (ppblob_p, ppblob_c, ppblob_s) not in state_per_patient:
        continue
    ppblob_s_norm = ppblob_s.strip().lower()
    ppblob_cord_vol_by_state_condition[(ppblob_s_norm, ppblob_c)] += (
        state_per_patient[(ppblob_p, ppblob_c, ppblob_s)]['volume'].sum())

# ============================================================
# Seventh plot: cumulative blob-size distribution grid, one panel per
# threshold (5, laid out 2x3 with one blank). Pooled across all patients -
# NOT per-patient (would be 12 x 2 x 2 = 48 curves per panel). Color =
# state, linestyle = condition.
# ============================================================
ppblob_all_r = [r for blobs in ppblob_pooled.values() for r, _ in blobs]
if not ppblob_all_r:
    raise SystemExit("No blobs found at any threshold - check the caches have data for both states.")
ppblob_r_min = min(ppblob_all_r) * 0.9
ppblob_r_max = max(ppblob_all_r) * 1.1

PP_LOG_MAJOR_LOCATOR = LogLocator(base=10.0)
PP_LOG_NULL_FORMATTER = NullFormatter()


def ppblob_plot_threshold(ax, threshold_name):
    for ppblob_state_name in (DELTA_TOP_STATE, DELTA_BOTTOM_STATE):
        for ppblob_condition_name in ('Flexion', 'Extension'):
            blobs = ppblob_pooled.get((threshold_name, ppblob_state_name, ppblob_condition_name))
            if not blobs:
                continue
            blobs_sorted = sorted(blobs, key=lambda b: b[0])
            total_vol = ppblob_cord_vol_by_state_condition[(ppblob_state_name, ppblob_condition_name)]
            if total_vol <= 0:
                continue
            rs = [ppblob_r_min] + [b[0] for b in blobs_sorted] + [ppblob_r_max]
            cum_vol = 0.0
            cum_pct = [0.0]
            for _, v in blobs_sorted:
                cum_vol += v
                cum_pct.append(100.0 * cum_vol / total_vol)
            cum_pct.append(cum_pct[-1])
            color = PPBLOB_STATE_COLORS[ppblob_state_name]
            linestyle = BLOB_CONDITION_LINESTYLES.get(ppblob_condition_name.strip().lower(), ':')
            ax.plot(rs, cum_pct, color=color, linestyle=linestyle, linewidth=1.5, drawstyle='steps-post')
    ax.set_xscale('log')
    ax.set_xlim(ppblob_r_min, ppblob_r_max)
    ax.set_ylim(0, 100)
    ax.xaxis.set_major_locator(PP_LOG_MAJOR_LOCATOR)
    ax.xaxis.set_minor_formatter(PP_LOG_NULL_FORMATTER)
    ax.set_title('Threshold {} = {:.2f}'.format(threshold_name.upper(), PPBLOB_THRESHOLDS[threshold_name]))
    ax.set_xlabel('MPS concentration effective radius r (mm)')
    ax.set_ylabel('Cumulative % of total cord volume above MPS threshold')


ppblob_state_handles = [Line2D([0], [0], color=PPBLOB_STATE_COLORS[s], linestyle='-',
                                label=DELTA_STATE_TITLES[s])
                         for s in (DELTA_TOP_STATE, DELTA_BOTTOM_STATE)]
ppblob_condition_handles = [Line2D([0], [0], color='black', linestyle=ls, label=cond.capitalize())
                             for cond, ls in BLOB_CONDITION_LINESTYLES.items()]

# --- Five individual per-threshold plots ---
for ppblob_threshold_name in PPBLOB_THRESHOLDS:
    ppblob_fig, ppblob_ax = plt.subplots(figsize=(7, 5.5))
    ppblob_plot_threshold(ppblob_ax, ppblob_threshold_name)
    ppblob_fig.legend(handles=ppblob_state_handles + ppblob_condition_handles,
                       loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)
    ppblob_fig.tight_layout()
    ppblob_plot_path = os.path.join(
        OUT_DIR, 'multipatient_mps_plot_blob_distribution_prepost_{}.pdf'.format(ppblob_threshold_name))
    ppblob_fig.savefig(ppblob_plot_path, bbox_inches='tight')
    plt.close(ppblob_fig)
    print("Plot saved: {}".format(ppblob_plot_path))

# --- Combined 2x3 grid (5 panels + 1 blank) ---
ppblob_grid_fig, ppblob_grid_axes = plt.subplots(2, 3, figsize=(18, 10))
for ppblob_ax_grid, ppblob_threshold_name in zip(ppblob_grid_axes.flat, PPBLOB_THRESHOLDS):
    ppblob_plot_threshold(ppblob_ax_grid, ppblob_threshold_name)
ppblob_grid_axes.flat[-1].axis('off')   # 6th slot unused (5 thresholds, 2x3 grid)

ppblob_grid_fig.legend(handles=ppblob_state_handles + ppblob_condition_handles,
                        loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)
ppblob_grid_fig.tight_layout()
ppblob_grid_plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot_blob_distribution_prepost_grid.pdf')
ppblob_grid_fig.savefig(ppblob_grid_plot_path, bbox_inches='tight')
plt.close(ppblob_grid_fig)
print("Plot saved: {}".format(ppblob_grid_plot_path))

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
