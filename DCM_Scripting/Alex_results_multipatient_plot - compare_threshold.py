"""
Alex_results_multipatient_plot.py - Pool MPS data across all patients AND
both loading conditions to compute cohort-wide volume-weighted 90th/95th/99th
percentile MPS thresholds (T90, T95, T99), then plot, per patient, the % of
cord volume exceeding each threshold for flexion and extension separately
(paired points per patient per threshold, connected by a grey line - same
style as Data_Extraction/mps_group_plot.py). Marker shape encodes loading
condition, color encodes threshold. X axis is the anonymized participant
number.

Reads DCM_Scripting/id_map.csv (gitignored - real patient IDs, never
committed) to find each patient's raw per-frame extraction CSV per loading
condition (from Alex_results_extraction.py, also never committed - see
'csv_path' column). Only the small aggregated summary, keyed by anonymized
participant number (e.g. 'P6') and loading condition, and the resulting
plot are written into the repo.

Run: python Alex_results_multipatient_plot.py
"""

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from mps_common import volume_weighted_percentile, pct_volume_above, PLOT_STYLE

# ============================================================
# USER SETTINGS
# ============================================================
ID_MAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'id_map.csv')

# 'last'  - use each element's MPS at the last frame of the step
# 'peak'  - use each element's peak (max-ever) MPS across all frames
# Decide using Alex_results_plotting.py's single-patient time-history
# diagnostic before running this at full cohort scale.
FRAME_MODE = 'peak'

THRESHOLD_PERCENTILES = {'t90': 0.90, 't95': 0.95, 't99': 0.99}
THRESHOLD_COLORS = {'t90': '#548235', 't95': '#2e75b6', 't99': '#c00000'}

# Fixed MPS thresholds (not derived from the data) for a second comparison plot
MANUAL_THRESHOLDS = {'t0p05': 0.05, 't0p10': 0.10, 't0p15': 0.15, 't0p20': 0.20}
MANUAL_THRESHOLD_COLORS = {'t0p05': '#548235', 't0p10': '#2e75b6', 't0p15': '#c00000', 't0p20': '#7030a0'}

CONDITION_MARKERS = {'flexion': 's', 'extension': '^'}

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
# ============================================================

plt.rcParams.update(PLOT_STYLE)


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


def build_summary(per_patient, thresholds):
    records = []
    for (participant, condition), df in sorted(per_patient.items()):
        record = {'participant': 'P{}'.format(participant), 'loading_condition': condition}
        for name, val in thresholds.items():
            record['pct_above_{}'.format(name)] = pct_volume_above(df, val)
        records.append(record)
    return pd.DataFrame(records)


def make_comparison_plot(summary, thresholds, threshold_colors, threshold_labels, plot_path):
    # x = participant (ordered by number). For each threshold, paired
    # flexion/extension points are connected by a thin line in that
    # threshold's color; marker shape encodes loading condition, marker/line
    # color encodes threshold - following mps_group_plot.py's paired-point
    # style, extended to multiple thresholds at once.
    participants = sorted(summary['participant'].unique(),
                           key=lambda p: (preop_mjoa_by_participant.get(p, float('inf')), int(p[1:])))
    x_pos = {p: i for i, p in enumerate(participants)}

    fig, ax = plt.subplots(figsize=(7, 5))

    for name in thresholds:
        col = 'pct_above_{}'.format(name)
        color = threshold_colors[name]
        for participant, grp in summary.groupby('participant'):
            if len(grp) == 2:
                xp = x_pos[participant]
                ax.vlines(xp, grp[col].min(), grp[col].max(), color=color, linewidth=1.0, alpha=0.5, zorder=2)
        for condition, grp in summary.groupby('loading_condition'):
            marker = CONDITION_MARKERS.get(condition.strip().lower(), 'o')
            xs = [x_pos[p] for p in grp['participant']]
            ax.scatter(xs, grp[col], color=color, marker=marker, s=55, zorder=3)

    # Single combined legend, top-left: Threshold entries, then Loading
    # condition entries directly underneath (blank marker/label entries act
    # as section headers) - keeps both groups together instead of split
    # across corners.
    threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=threshold_colors[name],
                                 label=threshold_labels[name])
                          for name in thresholds]
    condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                          for cond, marker in CONDITION_MARKERS.items()]
    blank_handle = Line2D([0], [0], linestyle='none', marker='None', label='')

    all_handles = (
        [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + threshold_handles +
        [blank_handle] +
        [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + condition_handles
    )
    fig.legend(handles=all_handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)

    ax.set_xticks(range(len(participants)))
    ax.set_xticklabels(['{}\n({})'.format(p, preop_mjoa_by_participant.get(p, '?'))
                         for p in participants])
    ax.set_xlabel('Participant (ordered by pre-op mJOA, ascending)')
    ax.set_ylabel('% cord volume above threshold')
    fig.tight_layout()

    fig.savefig(plot_path, bbox_inches='tight')
    plt.close(fig)


id_map = pd.read_csv(ID_MAP_PATH, skipinitialspace=True)

# PRE-OP mJOA per participant - used below to order patients on every plot
# by ascending pre-op mJOA (participant number as tie-break) and to label
# each patient with it. Explicitly filtered to State == 'PreOp' rows only -
# post-op mJOA can differ from pre-op for the same patient (e.g. N01-017:
# preop 11, postop 10), so this must not just take whichever row comes
# first in id_map.csv.
preop_mjoa_by_participant = {}
for _, row in id_map.iterrows():
    if str(row.get('State', '')).strip().lower() != 'preop':
        continue
    p_label = 'P{}'.format(int(row['participant']))
    preop_mjoa_by_participant[p_label] = row.get('mJOA', '')

per_patient = {}   # (participant, loading_condition) -> reduced DataFrame, PreOp only
missing = []
for _, row in id_map.iterrows():
    # PreOp only - without this filter, a PostOp row for the same
    # (participant, condition) silently overwrites the PreOp entry below,
    # since the dict key doesn't include state.
    state = str(row.get('State', '')).strip().lower()
    if state and state != 'preop':
        continue
    condition = str(row.get('loading_condition', '')).strip()
    if condition.strip().lower() not in ('flexion', 'extension'):
        continue
    csv_path = str(row.get('csv_path', '')).strip()
    if not csv_path or csv_path.lower() == 'nan' or not os.path.isfile(csv_path):
        missing.append(row)
        continue
    df = pd.read_csv(csv_path)
    per_patient[(int(row['participant']), condition)] = reduce_to_frame_mode(df, FRAME_MODE)

if missing:
    print("Skipping {} row(s) with no csv_path set in id_map.csv (or file not found):".format(len(missing)))
    for row in missing:
        cond = str(row.get('loading_condition', '')).strip() or '?'
        print("  P{} ({})".format(int(row['participant']), cond))

if not per_patient:
    raise SystemExit("No patient data loaded - fill in csv_path in id_map.csv first.")

# Pool ALL patients' AND both conditions' reduced (element, mps, volume) rows
# to compute cohort-wide volume-weighted T90/T95/T99. Elements weighted by
# volume, not counted equally - correct across patients with different cord
# sizes/mesh densities.
pooled = pd.concat(per_patient.values(), ignore_index=True)
thresholds = {name: volume_weighted_percentile(pooled, p=p) for name, p in THRESHOLD_PERCENTILES.items()}
print("Cohort-pooled thresholds (FRAME_MODE='{}'): ".format(FRAME_MODE) +
      "  ".join("{}={:.4f}".format(name.upper(), val) for name, val in thresholds.items()))

records = []
for (participant, condition), df in sorted(per_patient.items()):
    record = {'participant': 'P{}'.format(participant), 'loading_condition': condition}
    for name, val in thresholds.items():
        record['pct_above_{}'.format(name)] = pct_volume_above(df, val)
    records.append(record)
summary = pd.DataFrame(records)

summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary_compare_thresholds.csv')
summary.to_csv(summary_path, index=False)
print(summary.to_string(index=False))

# Plot: x = participant (ordered by number). For each threshold, paired
# flexion/extension points are connected by a thin line in that threshold's
# color; marker shape encodes loading condition, marker/line color encodes
# threshold - following mps_group_plot.py's paired-point style, extended to
# three thresholds at once.
participants = sorted(summary['participant'].unique(),
                       key=lambda p: (preop_mjoa_by_participant.get(p, float('inf')), int(p[1:])))
x_pos = {p: i for i, p in enumerate(participants)}

fig, ax = plt.subplots(figsize=(7, 5))

for name in thresholds:
    col = 'pct_above_{}'.format(name)
    color = THRESHOLD_COLORS[name]
    for participant, grp in summary.groupby('participant'):
        if len(grp) == 2:
            xp = x_pos[participant]
            ax.vlines(xp, grp[col].min(), grp[col].max(), color=color, linewidth=1.0, alpha=0.5, zorder=2)
    for condition, grp in summary.groupby('loading_condition'):
        marker = CONDITION_MARKERS.get(condition.strip().lower(), 'o')
        xs = [x_pos[p] for p in grp['participant']]
        ax.scatter(xs, grp[col], color=color, marker=marker, s=55, zorder=3)

# Single combined legend, top-left: Threshold entries, then Loading
# condition entries directly underneath (blank marker/label entries act as
# section headers) - keeps both groups together instead of split across corners.
threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=THRESHOLD_COLORS[name],
                             label='{} ({:.4f})'.format(name.upper(), thresholds[name]))
                      for name in thresholds]
condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                      for cond, marker in CONDITION_MARKERS.items()]
blank_handle = Line2D([0], [0], linestyle='none', marker='None', label='')

all_handles = (
    [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + threshold_handles +
    [blank_handle] +
    [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + condition_handles
)
fig.legend(handles=all_handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)

ax.set_xticks(range(len(participants)))
ax.set_xticklabels(['{}\n({})'.format(p, preop_mjoa_by_participant.get(p, '?'))
                     for p in participants])
ax.set_xlabel('Participant (ordered by pre-op mJOA, ascending)')
ax.set_ylabel('% cord volume above threshold')
fig.tight_layout()

plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot_compare_thresholds_sortedbypreopmJOA.pdf')
fig.savefig(plot_path, bbox_inches='tight')
plt.close(fig)

print("Summary saved: {}".format(summary_path))
print("Plot saved: {}".format(plot_path))

# ============================================================
# Second plot: manual fixed MPS thresholds (0.05 / 0.10 / 0.15 / 0.20)
# instead of percentile-derived ones. Reuses per_patient loaded above.
# ============================================================
manual_summary = build_summary(per_patient, MANUAL_THRESHOLDS)
manual_summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary_compare_thresholds_manual_thresholds.csv')
manual_summary.to_csv(manual_summary_path, index=False)
print(manual_summary.to_string(index=False))

manual_labels = {name: '{:.2f}'.format(val) for name, val in MANUAL_THRESHOLDS.items()}
manual_plot_path = os.path.join(
    OUT_DIR, 'multipatient_mps_plot_compare_thresholds_manual_thresholds_sortedbypreopmJOA.pdf')
make_comparison_plot(manual_summary, MANUAL_THRESHOLDS, MANUAL_THRESHOLD_COLORS, manual_labels, manual_plot_path)

print("Summary saved: {}".format(manual_summary_path))
print("Plot saved: {}".format(manual_plot_path))

# ============================================================
# Third plot: cumulative blob-size distribution (spatial clustering).
# Distinguishes a few large contiguous high-strain regions ("blob") from
# many small scattered ones ("chicken pox"), for all 4 manual thresholds and
# both loading conditions at once. Blobs are connected components of
# exceeding-threshold elements (face-sharing adjacency, approximated as
# >=4 shared nodes),
# pooled across ALL patients per (threshold, condition) - same pooling
# convention already used above for the percentile/manual thresholds - so
# this produces 4 thresholds x 2 conditions = 8 curves, not one per patient.
#
# Requires the companion '<basename>_topology.csv' files (element
# connectivity) next to each '_mps.csv', derived from csv_path by suffix
# swap - re-run Alex_results_extraction.py per patient/condition to
# generate them if missing (older extractions won't have them yet).
#
# Reuses id_map, per_patient, MANUAL_THRESHOLDS, MANUAL_THRESHOLD_COLORS,
# OUT_DIR, Line2D already loaded/defined above - does not modify anything
# above this point.
# ============================================================
import math
from collections import defaultdict

BLOB_FACE_SHARING_MIN_NODES = 4   # >=4 shared nodes approximates shared face (all elements are C3D8)
BLOB_CONDITION_LINESTYLES = {'flexion': '-', 'extension': '--'}


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


# Cache adjacency edges per (participant, condition) - reused across all 4 thresholds
blob_adjacency_cache = {}
blob_missing_topology = []

for _, row in id_map.iterrows():
    # PreOp only - same reasoning as per_patient above: without this, a
    # PostOp row for the same (participant, condition) could pair PostOp
    # topology (different mesh) with PreOp field data for that key.
    blob_state = str(row.get('State', '')).strip().lower()
    if blob_state and blob_state != 'preop':
        continue
    blob_participant = int(row['participant'])
    blob_condition = str(row.get('loading_condition', '')).strip()
    blob_csv_path = str(row.get('csv_path', '')).strip()
    blob_key = (blob_participant, blob_condition)
    if blob_key not in per_patient or not blob_csv_path or blob_csv_path.lower() == 'nan':
        continue
    blob_topology_path = blob_csv_path.replace('_mps.csv', '_topology.csv')
    if not os.path.isfile(blob_topology_path):
        blob_missing_topology.append(blob_key)
        continue
    blob_adjacency_cache[blob_key] = load_adjacency_edges(blob_topology_path)

if blob_missing_topology:
    print("Skipping {} patient/condition(s) missing '_topology.csv' for blob analysis "
          "(re-run Alex_results_extraction.py to generate it):".format(len(blob_missing_topology)))
    for blob_participant, blob_condition in blob_missing_topology:
        print("  P{} ({})".format(blob_participant, blob_condition))

# For each (threshold, condition), pool blobs from every patient with adjacency data
blob_records = []
blob_pooled = defaultdict(list)   # (threshold_name, condition) -> list of (r, volume)

for blob_threshold_name, blob_threshold_val in MANUAL_THRESHOLDS.items():
    for (blob_participant, blob_condition), blob_edges in blob_adjacency_cache.items():
        blob_df_patient = per_patient[(blob_participant, blob_condition)]
        blob_vol_lookup = dict(zip(blob_df_patient['element_label'], blob_df_patient['volume']))
        blob_exceeding = set(blob_df_patient.loc[blob_df_patient['mps'] >= blob_threshold_val, 'element_label'])
        if not blob_exceeding:
            continue
        for blob in find_blobs(blob_exceeding, blob_edges):
            blob_volume = sum(blob_vol_lookup[lbl] for lbl in blob)
            blob_r = (3.0 * blob_volume / (4.0 * math.pi)) ** (1.0 / 3.0)
            blob_pooled[(blob_threshold_name, blob_condition)].append((blob_r, blob_volume))
            blob_records.append({
                'participant':       'P{}'.format(blob_participant),
                'loading_condition': blob_condition,
                'threshold':         blob_threshold_name,
                'n_elements':        len(blob),
                'volume':            blob_volume,
                'r':                 blob_r,
            })

blob_df = pd.DataFrame(blob_records)
blob_summary_path = os.path.join(OUT_DIR, 'multipatient_blob_distribution_summary.csv')
blob_df.to_csv(blob_summary_path, index=False)
print("Blob summary saved: {}".format(blob_summary_path))

# ============================================================
# Fourth plot: per-patient faceted breakdown of the same blob distribution,
# one curve per (patient, condition). Each patient's curve is normalized to
# THEIR OWN total cord volume (fixed regardless of threshold), so a curve's
# final height shows the real % of that patient's cord exceeding the
# threshold - comparable across patients and across thresholds, rather than
# every curve being forced to reach 100% regardless of magnitude. Produces
# one figure per threshold
# (4 total) plus a single 2x2 grid image combining all four. Reuses blob_df,
# MANUAL_THRESHOLDS, BLOB_CONDITION_LINESTYLES, OUT_DIR, Line2D, pd, plt, os
# already loaded/defined above - does not modify anything above this point.
# ============================================================
from matplotlib.ticker import LogLocator, NullFormatter

PP_LOG_MAJOR_LOCATOR = LogLocator(base=10.0)
PP_LOG_NULL_FORMATTER = NullFormatter()

pp_patients_sorted = sorted(blob_df['participant'].unique(),
                             key=lambda p: (preop_mjoa_by_participant.get(p, float('inf')), int(p[1:])))
pp_patient_colors = {p: plt.cm.tab10(i % 10) for i, p in enumerate(pp_patients_sorted)}

# Common x/y range across ALL thresholds, so every panel is directly
# comparable - otherwise each panel auto-scales to its own data subset,
# which makes cross-threshold shifts misleading rather than informative.
pp_r_min = blob_df['r'].min()
pp_r_max = blob_df['r'].max()
pp_xlim = (pp_r_min * 0.9, pp_r_max * 1.1)
pp_ylim = (0, 102)

# Each patient's TOTAL CORD VOLUME (not their exceeding-volume subset) -
# fixed regardless of threshold, so a curve's final height shows the real %
# of that patient's cord exceeding the threshold, comparable across
# thresholds and patients instead of every curve being forced to 100%.
pp_total_cord_vol = {key: df['volume'].sum() for key, df in per_patient.items()}


def pp_plot_threshold(ax, threshold_name):
    sub = blob_df[blob_df['threshold'] == threshold_name]
    for (pp_participant, pp_condition), grp in sub.groupby(['participant', 'loading_condition']):
        grp_sorted = grp.sort_values('r')
        pp_participant_num = int(pp_participant[1:])   # 'P6' -> 6, to match per_patient's int key
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
        ax.plot(pp_rs, pp_ys, color=pp_color, linestyle=pp_linestyle, linewidth=1.2,
                drawstyle='steps-post')
    ax.set_xscale('log')
    ax.set_xlim(pp_xlim)
    ax.set_ylim(pp_ylim)
    # Major ticks only at clean powers of ten (10^-1, 10^0, 10^1, ...) - log
    # axes otherwise default to also labeling minor ticks, which gets
    # cluttered/unreadable once several decades are spanned.
    ax.xaxis.set_major_locator(PP_LOG_MAJOR_LOCATOR)
    ax.xaxis.set_minor_formatter(PP_LOG_NULL_FORMATTER)
    ax.set_title('Threshold {} = {:.2f}'.format(threshold_name.upper(), MANUAL_THRESHOLDS[threshold_name]))
    ax.set_xlabel('MPS concentration effective radius r (mm)')
    ax.set_ylabel('Cumulative % of total cord volume above MPS threshold')


pp_patient_handles = [Line2D([0], [0], color=pp_patient_colors[p], linestyle='-',
                              label='{} (preop mJOA {})'.format(p, preop_mjoa_by_participant.get(p, '?')))
                       for p in pp_patients_sorted]
pp_condition_handles = [Line2D([0], [0], color='black', linestyle=ls, label=cond.capitalize())
                         for cond, ls in BLOB_CONDITION_LINESTYLES.items()]

# --- Four individual per-threshold plots ---
for pp_threshold_name in MANUAL_THRESHOLDS:
    pp_fig, pp_ax = plt.subplots(figsize=(7, 5.5))
    pp_plot_threshold(pp_ax, pp_threshold_name)

    pp_fig.legend(handles=pp_patient_handles + pp_condition_handles,
                  loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)

    pp_fig.tight_layout()
    pp_plot_path = os.path.join(
        OUT_DIR, 'multipatient_mps_plot_blob_distribution_perpatient_{}_sortedbypreopmJOA.pdf'.format(
            pp_threshold_name))
    pp_fig.savefig(pp_plot_path, bbox_inches='tight')
    plt.close(pp_fig)
    print("Plot saved: {}".format(pp_plot_path))

# --- Combined 2x2 grid of all four thresholds in one image ---
pp_grid_fig, pp_grid_axes = plt.subplots(2, 2, figsize=(13, 10))
for pp_ax_grid, pp_threshold_name in zip(pp_grid_axes.flat, MANUAL_THRESHOLDS):
    pp_plot_threshold(pp_ax_grid, pp_threshold_name)

pp_grid_fig.legend(handles=pp_patient_handles + pp_condition_handles,
                    loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)
pp_grid_fig.tight_layout()
pp_grid_plot_path = os.path.join(
    OUT_DIR, 'multipatient_mps_plot_blob_distribution_perpatient_grid_sortedbypreopmJOA.pdf')
pp_grid_fig.savefig(pp_grid_plot_path, bbox_inches='tight')
plt.close(pp_grid_fig)
print("Plot saved: {}".format(pp_grid_plot_path))

# ============================================================
# Fifth plot: manual-threshold comparison (0.05/0.10/0.15/0.20) INCLUDING
# both PreOp and PostOp states, per patient - a separate, new plot, not a
# modification of the PreOp-only manual-thresholds plot earlier in this
# file. Same visual encoding as that plot (color = threshold, marker shape =
# loading condition), plus marker FILL for state: solid = PreOp, hollow =
# PostOp. Flexion/Extension pairs are connected only within the same state.
# Single combined legend to the right of the plot (not overlapping it),
# grouped into Threshold / Loading condition / State sections.
#
# Builds its own (participant, condition, state) -> reduced DataFrame dict
# from id_map's 'State' column (independent of per_patient above, which is
# PreOp-only) - reuses id_map, reduce_to_frame_mode, pct_volume_above,
# MANUAL_THRESHOLDS, MANUAL_THRESHOLD_COLORS, CONDITION_MARKERS, OUT_DIR,
# Line2D, pd, plt, os already loaded/defined above - does not modify
# anything above this point.
# ============================================================
STATE_FILLED = {'preop': True, 'postop': False}   # solid marker vs hollow marker

state_per_patient = {}   # (participant, loading_condition, state) -> reduced DataFrame, ALL states
state_missing = []
for _, row in id_map.iterrows():
    state_condition = str(row.get('loading_condition', '')).strip()
    if state_condition.strip().lower() not in ('flexion', 'extension'):
        # Excludes id_map's 'Oscillation' rows - this dict feeds both the
        # fifth and sixth plots below, neither of which defines a
        # marker/style for a third loading condition, so Oscillation rows
        # would otherwise silently plot with the scatter default marker 'o'
        # mixed in among the Flexion/Extension points.
        continue
    state_csv_path = str(row.get('csv_path', '')).strip()
    state_state = str(row.get('State', '')).strip() or 'PreOp'
    if not state_csv_path or state_csv_path.lower() == 'nan' or not os.path.isfile(state_csv_path):
        state_missing.append((int(row['participant']), state_condition, state_state))
        continue
    state_df_raw = pd.read_csv(state_csv_path)
    state_per_patient[(int(row['participant']), state_condition, state_state)] = reduce_to_frame_mode(
        state_df_raw, FRAME_MODE)

if state_missing:
    print("Skipping {} row(s) with no csv_path set in id_map.csv (or file not found) "
          "for PreOp/PostOp plot:".format(len(state_missing)))
    for sp, sc, ss in state_missing:
        print("  P{} ({}, {})".format(sp, sc, ss))

state_records = []
for (state_participant, state_condition, state_state), state_df in sorted(state_per_patient.items()):
    state_record = {'participant': 'P{}'.format(state_participant), 'loading_condition': state_condition,
                     'state': state_state}
    for name, val in MANUAL_THRESHOLDS.items():
        state_record['pct_above_{}'.format(name)] = pct_volume_above(state_df, val)
    state_records.append(state_record)
state_summary = pd.DataFrame(state_records)

state_summary_path = os.path.join(
    OUT_DIR, 'multipatient_mps_summary_compare_thresholds_manual_thresholds_prepost.csv')
state_summary.to_csv(state_summary_path, index=False)
print(state_summary.to_string(index=False))

state_participants = sorted(state_summary['participant'].unique(),
                             key=lambda p: (preop_mjoa_by_participant.get(p, float('inf')), int(p[1:])))
state_x_pos = {p: i for i, p in enumerate(state_participants)}

# Plot only shows T0.10/T0.15 (summary CSV above still has all 4 thresholds
# for reference) - keeps this comparison plot readable.
STATE_PLOT_THRESHOLDS = {'t0p10': MANUAL_THRESHOLDS['t0p10'], 't0p15': MANUAL_THRESHOLDS['t0p15']}

state_fig, state_ax = plt.subplots(figsize=(9, 5.5))

for state_threshold_name in STATE_PLOT_THRESHOLDS:
    state_col = 'pct_above_{}'.format(state_threshold_name)
    state_color = MANUAL_THRESHOLD_COLORS[state_threshold_name]
    # Connect Flexion<->Extension pairs within the SAME state only
    for (_, _), grp in state_summary.groupby(['participant', 'state']):
        if len(grp) == 2:
            xp = state_x_pos[grp['participant'].iloc[0]]
            state_ax.vlines(xp, grp[state_col].min(), grp[state_col].max(),
                             color=state_color, linewidth=1.0, alpha=0.5, zorder=2)
    for (condition, state), grp in state_summary.groupby(['loading_condition', 'state']):
        marker = CONDITION_MARKERS.get(condition.strip().lower(), 'o')
        filled = STATE_FILLED.get(state.strip().lower(), True)
        xs = [state_x_pos[p] for p in grp['participant']]
        if filled:
            state_ax.scatter(xs, grp[state_col], color=state_color, marker=marker, s=55, zorder=3)
        else:
            state_ax.scatter(xs, grp[state_col], facecolors='none', edgecolors=state_color, marker=marker,
                              s=55, linewidths=1.3, zorder=3)

# Single combined legend to the right of the plot (not overlapping), grouped
# into Threshold / Loading condition / State sections via blank header entries.
state_threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=MANUAL_THRESHOLD_COLORS[name],
                                   label='{:.2f}'.format(val))
                            for name, val in STATE_PLOT_THRESHOLDS.items()]
state_condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                            for cond, marker in CONDITION_MARKERS.items()]
state_state_handles = [
    Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='black', label='PreOp'),
    Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='none', label='PostOp'),
]
state_blank = Line2D([0], [0], linestyle='none', marker='None', label='')

state_all_handles = (
    [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + state_threshold_handles +
    [state_blank] +
    [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + state_condition_handles +
    [state_blank] +
    [Line2D([0], [0], linestyle='none', marker='None', label='State')] + state_state_handles
)
state_ax.legend(handles=state_all_handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

state_ax.set_xticks(range(len(state_participants)))
state_ax.set_xticklabels(['{}\n({})'.format(p, preop_mjoa_by_participant.get(p, '?'))
                           for p in state_participants])
# Explicitly "pre-op" here since this plot shows PostOp points too, and
# post-op mJOA can differ from pre-op for the same patient - the x-axis
# ordering/label always uses pre-op mJOA regardless of a point's own state.
state_ax.set_xlabel('Participant (ordered by pre-op mJOA, ascending)')
state_ax.set_ylabel('% cord volume above threshold')
state_fig.tight_layout()

state_plot_path = os.path.join(
    OUT_DIR, 'multipatient_mps_plot_compare_thresholds_manual_thresholds_prepost_sortedbypreopmJOA.pdf')
state_fig.savefig(state_plot_path, bbox_inches='tight')
plt.close(state_fig)

print("Summary saved: {}".format(state_summary_path))
print("Plot saved: {}".format(state_plot_path))

# ============================================================
# Sixth plot: PreOp-NoPreload vs PostOp, 5 low manual thresholds
# (0.01/0.02/0.03/0.04/0.05), 2 rows x 1 column (top = PreOp WITHOUT the
# simulated preload step, bottom = PostOp, sharing the same x order) -
# patients ordered by CHANGE in mJOA (postop - preop), ascending (most
# worsened on the left, most improved on the right), rather than raw pre-op
# mJOA. mJOA itself is a clinical property of the patient, not of which
# simulation variant is plotted, so the delta/ordering still comes from the
# regular PreOp/PostOp mJOA values regardless of the top row using the
# NoPreload model. Only patients with BOTH a pre-op and post-op mJOA value
# can have a delta, so anyone missing either is left out of this plot only.
#
# Reuses state_per_patient (built for the fifth plot - already restricted to
# Flexion/Extension, all States including PreOp-NoPreload), id_map,
# preop_mjoa_by_participant, CONDITION_MARKERS, OUT_DIR, Line2D, pd, plt, os,
# pct_volume_above already loaded/defined above - does not modify anything
# above this point.
# ============================================================
from matplotlib.patches import Patch

DELTA_THRESHOLDS = {'t0p01': 0.01, 't0p02': 0.02, 't0p03': 0.03}#, 't0p04': 0.04, 't0p05': 0.05}
DELTA_THRESHOLD_COLORS = {
    't0p01': '#2a78d6',   # blue
    't0p02': '#eb6834',   # orange
    't0p03': '#1baf7a',   # aqua
    #'t0p04': '#1c5cab',
    #'t0p05': '#104281',
}

DELTA_TOP_STATE = 'preop-nopreload'     # top row: PreOp WITHOUT simulated preload
DELTA_BOTTOM_STATE = 'postop'           # bottom row: PostOp (unchanged)
DELTA_STATE_TITLES = {DELTA_TOP_STATE: 'PreOp (no preload)', DELTA_BOTTOM_STATE: 'PostOp'}

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

# Cache the per-element (mps, volume) data this plot depends on - the
# threshold-independent, expensive-to-recompute part (reading + reducing
# every patient's raw per-frame CSV) - so plot_prepost_sortedbymjoachange.py
# can replot with different DELTA_THRESHOLDS/colors without re-running this
# whole script. Keyed by FRAME_MODE since 'peak' vs 'last' changes the values.
delta_cache_path = os.path.join(
    OUT_DIR, 'cache_prepost_sortedbymjoachange_{}.csv'.format(FRAME_MODE))
delta_cache_parts = []
for (cache_participant, cache_condition, cache_state), cache_df in state_per_patient.items():
    if cache_state.strip().lower() not in (DELTA_TOP_STATE, DELTA_BOTTOM_STATE):
        continue
    cache_part = cache_df[['element_label', 'mps', 'volume']].copy()
    cache_part.insert(0, 'state', cache_state)
    cache_part.insert(0, 'loading_condition', cache_condition)
    cache_part.insert(0, 'participant', cache_participant)
    delta_cache_parts.append(cache_part)
pd.concat(delta_cache_parts, ignore_index=True).to_csv(delta_cache_path, index=False)
print("Cached reduced data for standalone replotting: {}".format(delta_cache_path))

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

# Fusion patients (surgical detail - only shaded on the PostOp row, since
# fusion is a post-operative property and has no PreOp meaning).
DELTA_FUSION_PARTICIPANTS = {'P1', 'P5', 'P7', 'P8'}   

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
# Blob adjacency CACHE (PreOp-NoPreload vs PostOp) - not a plot itself.
# Parsing each job's '_topology.csv' and finding shared-node element pairs
# (load_adjacency_edges) is the slow part of blob analysis and is entirely
# threshold-independent, so it's computed ONCE here and written to a CSV
# cache. plot_prepost_blob_distribution.py (a standalone script, same
# fast-iteration pattern as plot_prepost_sortedbymjoachange.py) reads this
# cache plus the existing cache_prepost_sortedbymjoachange_{FRAME_MODE}.csv
# to do the actual blob-finding and plotting in seconds, without re-running
# this whole pipeline or re-parsing topology every time thresholds/colors
# change.
#
# Needs a SEPARATE adjacency cache keyed by (participant, condition, state)
# - PreOp-NoPreload and PostOp are different jobs/meshes for the same
# patient, so their topology must be loaded and kept separately (the
# existing blob_adjacency_cache above is keyed (participant, condition)
# only, PreOp-specific).
#
# Reuses id_map, state_per_patient (built for the Fifth/Sixth plots -
# already covers every State including PreOp-NoPreload/PostOp, already
# Flexion/Extension only), load_adjacency_edges, DELTA_TOP_STATE/
# BOTTOM_STATE, OUT_DIR, pd, os already loaded/defined above - does not
# modify anything above this point.
# ============================================================
# state_per_patient is keyed by the ORIGINAL-case State string from id_map.csv
# (its own builder: state_state = str(row.get('State', '')).strip() or 'PreOp')
# - ppblob_key below matches that convention directly (no lowercasing), so
# it both indexes state_per_patient correctly AND keeps this cache's 'state'
# column consistent with cache_prepost_sortedbymjoachange_{FRAME_MODE}.csv's
# own 'state' column (also original-case) that the standalone script reads
# alongside it.
ppblob_adjacency_cache = {}
ppblob_missing_topology = []

for _, row in id_map.iterrows():
    ppblob_state = str(row.get('State', '')).strip() or 'PreOp'
    if ppblob_state.strip().lower() not in (DELTA_TOP_STATE, DELTA_BOTTOM_STATE):
        continue
    ppblob_condition = str(row.get('loading_condition', '')).strip()
    if ppblob_condition.strip().lower() not in ('flexion', 'extension'):
        continue
    ppblob_participant = int(row['participant'])
    ppblob_csv_path = str(row.get('csv_path', '')).strip()
    ppblob_key = (ppblob_participant, ppblob_condition, ppblob_state)
    if ppblob_key not in state_per_patient or not ppblob_csv_path or ppblob_csv_path.lower() == 'nan':
        continue
    ppblob_topology_path = ppblob_csv_path.replace('_mps.csv', '_topology.csv')
    if not os.path.isfile(ppblob_topology_path):
        ppblob_missing_topology.append(ppblob_key)
        continue
    ppblob_adjacency_cache[ppblob_key] = load_adjacency_edges(ppblob_topology_path)

if ppblob_missing_topology:
    print("Skipping {} patient/condition/state combo(s) missing '_topology.csv' for PreOp-NoPreload/PostOp "
          "blob analysis:".format(len(ppblob_missing_topology)))
    for ppblob_participant, ppblob_condition, ppblob_state in ppblob_missing_topology:
        print("  P{} ({}, {})".format(ppblob_participant, ppblob_condition, ppblob_state))


# Flatten (participant, condition, state) -> [(elem_a, elem_b), ...] into a
# CSV: one row per adjacency edge. Written once here; read back fast by
# plot_prepost_blob_distribution.py on every subsequent run.
ppblob_adjacency_rows = []
for (ppblob_cache_p, ppblob_cache_c, ppblob_cache_s), ppblob_cache_edges in ppblob_adjacency_cache.items():
    for ppblob_elem_a, ppblob_elem_b in ppblob_cache_edges:
        ppblob_adjacency_rows.append({
            'participant':       ppblob_cache_p,
            'loading_condition': ppblob_cache_c,
            'state':             ppblob_cache_s,
            'elem_a':            ppblob_elem_a,
            'elem_b':            ppblob_elem_b,
        })

ppblob_adjacency_cache_path = os.path.join(
    OUT_DIR, 'cache_prepost_blob_adjacency_{}.csv'.format(FRAME_MODE))
pd.DataFrame(ppblob_adjacency_rows).to_csv(ppblob_adjacency_cache_path, index=False)
print("Cached blob adjacency edges for standalone replotting: {} ({} edges, {} (participant, condition, state) "
      "combos)".format(ppblob_adjacency_cache_path, len(ppblob_adjacency_rows), len(ppblob_adjacency_cache)))
