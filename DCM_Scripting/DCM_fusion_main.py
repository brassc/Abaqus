"""
DCM_fusion_main.py - PostOp laminectomy-only vs laminectomy + fusion.

Every patient has BOTH surgeries modelled (factual rows in id_map.csv,
counterfactual rows in id_map_counterfactual.csv), for Flexion and Extension,
no preload. Patients are grouped by classification.csv Category
(1 = Normal, 2 = Residual compression, 3 = S-curve).

Stage 1 (this file, current): data assembly + plots only.
  Part A: % cord volume above threshold (0.015 + sensitivity thresholds)
  Part B: % IVD volume above T95 (global and patient-specific)
  Part C: blob-size distributions (one grid per surgery arm)
Part D: statistics. Flexion and Extension are always analysed separately, never pooled.
  Primary: exact Wilcoxon signed-rank (Pratt) on the 12 paired differences (fusion - lami).
  Also: permutation test on the mean, Hodges-Lehmann, sign test, Kenward-Roger LME (+ QQ plots),
  exploratory surgery x category, and a confounded factual-only fusion vs lami comparison.

Category source: CATEGORY_SOURCE = 'PostOp' first; set to 'PreOp' and re-run
for the PreOp categorisation (only the grouping/outputs change - caches are
category-independent).

Run: python DCM_fusion_main.py (first run is slow - builds caches from raw data).
"""

import math
import os
import sys
from collections import defaultdict

# Pin matplotlib's PDF CreationDate so identical data gives identical PDFs
# (same reason as DCM_main.py). Must be set before the matplotlib import.
os.environ.setdefault('SOURCE_DATE_EPOCH', '1700000000')

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import LogLocator, NullFormatter


# ============================================================
# Shared volume helpers (same definitions as DCM_main.py)
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
    total volume lies (Kleiven 2007; see DCM_Scripting/README.md)."""
    s = df.sort_values(mps_col)
    cumvol = s[vol_col].cumsum()
    return s.loc[cumvol >= p * s[vol_col].sum(), mps_col].iloc[0]


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
plt.rcParams.update(PLOT_STYLE)

# ============================================================
# USER SETTINGS
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ID_MAP_PATH = os.path.join(SCRIPT_DIR, 'id_map.csv')
ID_MAP_CF_PATH = os.path.join(SCRIPT_DIR, 'id_map_counterfactual.csv')
CLASSIFICATION_PATH = os.path.join(SCRIPT_DIR, 'classification.csv')

# 'PostOp' first; switch to 'PreOp' to regroup patients by pre-op categorisation.
CATEGORY_SOURCE = 'PostOp'

BASE_DIR = os.path.join(SCRIPT_DIR, 'fusion_vs_lami_results')
CACHE_DIR = os.path.join(BASE_DIR, 'cache')                      # category-independent
OUT_DIR = os.path.join(BASE_DIR, 'category_{}'.format(CATEGORY_SOURCE.lower()))
os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

# Must match across caches. 'peak' = each element's max MPS over all frames.
FRAME_MODE = 'peak'

# 0.015 is the working threshold used for the whole-cord/GM-WM comparisons in
# DCM_main.py (MID_THRESHOLD). No derivation of it is documented in the repo,
# so the other values are plotted as a sensitivity check.
PRIMARY_THRESHOLD = 0.015
CORD_THRESHOLDS = [0.010, 0.015, 0.020, 0.025]

IVD_PERCENTILE = 0.95   # T95

# Blob analysis (same definitions as DCM_main.py)
BLOB_THRESHOLDS = {'t0p01': 0.01, 't0p015': 0.015, 't0p02': 0.02, 't0p025': 0.025}
BLOB_FACE_SHARING_MIN_NODES = 4   # >=4 shared nodes approximates shared face (all C3D8)
BLOB_CONDITION_LINESTYLES = {'Flexion': '-', 'Extension': '--'}

LOADING_CONDITIONS = ['Flexion', 'Extension']   # always analysed separately
ARMS = ['lami', 'fusion']
ARM_LABELS = {'lami': 'Laminectomy Only', 'fusion': 'Laminectomy & Fusion'}   # display only; data keys stay 'lami'/'fusion'
ARM_FILLED = {'lami': True, 'fusion': False}      # solid = lami only, hollow = lami + fusion
CONDITION_MARKERS = {'Flexion': 's', 'Extension': '^'}
ARM_X_OFFSET = {'lami': -0.17, 'fusion': 0.17}
MARKER_COLOR = '#2e75b6'
CATEGORY_NAMES = {1: 'Normal', 2: 'Residual compression', 3: 'S-curve'}
FUSION_SPAN_COLOR = 'orange'    # ground-truth fusion patients (same as DCM_main.py)
GROUP_GAP = 0.8                 # extra x spacing between category groups


# ============================================================
# Blob helpers (same as DCM_main.py)
# ============================================================
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
        uf.find(lbl)
    for a, b in edges:
        if a in exceeding_labels and b in exceeding_labels:
            uf.union(a, b)
    groups = defaultdict(list)
    for lbl in exceeding_labels:
        groups[uf.find(lbl)].append(lbl)
    return list(groups.values())


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


# ============================================================
# Assemble the job table: 12 patients x 2 arms x 2 loadings (PostOp, no preload)
# ============================================================
def _clean(df):
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype(str).str.strip()
    return df


id_map = _clean(pd.read_csv(ID_MAP_PATH, skipinitialspace=True))
id_map_cf = _clean(pd.read_csv(ID_MAP_CF_PATH, skipinitialspace=True))
classification = _clean(pd.read_csv(CLASSIFICATION_PATH, skipinitialspace=True))

factual = id_map[(id_map['State'] == 'PostOp') &
                 (id_map['loading_condition'].isin(LOADING_CONDITIONS))].copy()
counterfactual = id_map_cf[(id_map_cf['State'] == 'PostOp') &
                           (id_map_cf['loading_condition'].isin(LOADING_CONDITIONS))].copy()
jobs = pd.concat([factual, counterfactual], ignore_index=True)
jobs['participant'] = jobs['participant'].astype(int)
jobs['arm'] = jobs['surgery_modelled']                  # 'lami' or 'fusion'
jobs['is_factual'] = jobs['counterfactual'].astype(int) == 0
jobs['p_label'] = 'P' + jobs['participant'].astype(str)

# Ground truth fusion = what the patient actually had (id_map 'surgery' column).
fusion_participants = set(jobs.loc[jobs['surgery'] == 'fusion', 'p_label'])
# classification.csv also carries a Fusion flag on PostOp rows - cross-check.
_cls_fusion = set('P' + classification.loc[
    (classification['State'] == 'PostOp') & (pd.to_numeric(classification['Fusion'], errors='coerce') == 1),
    'participant'].astype(str))
if _cls_fusion != fusion_participants:
    print("WARNING: fusion patients disagree. id_map: {}  classification.csv: {}".format(
        sorted(fusion_participants), sorted(_cls_fusion)))

category_by_participant = {
    'P{}'.format(int(r['participant'])): int(r['Category'])
    for _, r in classification[classification['State'] == CATEGORY_SOURCE].iterrows()}
jobs['category'] = jobs['p_label'].map(category_by_participant)

# Structural checks: exactly one job per (patient, loading, arm), and 1 factual + 1 counterfactual
# arm per patient/loading.
_dups = jobs.duplicated(['participant', 'loading_condition', 'arm'], keep=False)
if _dups.any():
    raise SystemExit("Duplicate (patient, loading, arm) jobs:\n{}".format(
        jobs.loc[_dups, ['participant', 'loading_condition', 'arm', 'is_factual']]))
_expected = len(set(jobs['participant'])) * len(LOADING_CONDITIONS) * len(ARMS)
if len(jobs) != _expected:
    raise SystemExit("Expected {} jobs (patients x loadings x arms) but found {}".format(
        _expected, len(jobs)))
if jobs['category'].isna().any():
    raise SystemExit("No {} category for: {}".format(
        CATEGORY_SOURCE, sorted(jobs.loc[jobs['category'].isna(), 'p_label'].unique())))
jobs['category'] = jobs['category'].astype(int)

print("Jobs: {} ({} patients). Category source: {}.".format(
    len(jobs), jobs['participant'].nunique(), CATEGORY_SOURCE))
print("Ground-truth fusion patients: {}".format(sorted(fusion_participants, key=lambda p: int(p[1:]))))
for cat, name in sorted(CATEGORY_NAMES.items()):
    members = sorted(jobs.loc[jobs['category'] == cat, 'p_label'].unique(), key=lambda p: int(p[1:]))
    print("  Category {} ({}): n={}  {}".format(cat, name, len(members), members))


# ============================================================
# Reduced-data caches (cord + IVD): one row per element per job
# ============================================================
def build_reduced_cache(cache_name, path_fn, label):
    cache_path = os.path.join(CACHE_DIR, cache_name)
    if os.path.isfile(cache_path):
        print("Loading cached {} data: {}".format(label, cache_path))
        return pd.read_csv(cache_path)
    print("No {} cache at {} - building (slow, one-time)...".format(label, cache_path))
    parts, missing = [], []
    for _, row in jobs.iterrows():
        csv_path = path_fn(str(row['csv_path']).strip())
        if not csv_path or csv_path.lower() == 'nan' or not os.path.isfile(csv_path):
            missing.append((row['p_label'], row['loading_condition'], row['arm']))
            continue
        raw = pd.read_csv(csv_path, usecols=['frame_index', 'element_label', 'mps', 'volume'])
        reduced = reduce_to_frame_mode(raw, FRAME_MODE).copy()
        reduced.insert(0, 'arm', row['arm'])
        reduced.insert(0, 'loading_condition', row['loading_condition'])
        reduced.insert(0, 'participant', int(row['participant']))
        parts.append(reduced)
    if missing:
        # Not cached: a partial cache would be trusted on the next run.
        raise SystemExit("Missing {} data for {} job(s) - run the extraction first:\n  {}".format(
            label, len(missing), "\n  ".join(str(m) for m in missing)))
    df = pd.concat(parts, ignore_index=True)
    df.to_csv(cache_path, index=False)
    print("Cached {} data: {}".format(label, cache_path))
    return df


cord_df = build_reduced_cache('cache_fusion_cord_{}.csv'.format(FRAME_MODE), lambda p: p, 'cord')
ivd_df = build_reduced_cache('cache_fusion_ivd_{}.csv'.format(FRAME_MODE),
                             lambda p: p.replace('_mps.csv', '_ivd_mps.csv'), 'IVD')

job_info = jobs.set_index(['participant', 'loading_condition', 'arm'])


def attach_info(summary):
    """Add p_label / category / ground-truth fusion flag to a per-job summary frame."""
    summary['p_label'] = 'P' + summary['participant'].astype(str)
    summary['category'] = summary['p_label'].map(category_by_participant)
    summary['gt_fusion'] = summary['p_label'].isin(fusion_participants)
    return summary


# ============================================================
# Plot helper: Flexion | Extension panels, patients grouped by category,
# lami vs fusion markers joined per patient, orange fill = actually fused.
# ============================================================
def x_layout():
    """x position per patient: grouped by category, then patient number."""
    positions, group_spans = {}, {}
    x = 0.0
    for cat in sorted(CATEGORY_NAMES):
        members = sorted((p for p, c in category_by_participant.items()
                          if c == cat and p in set(jobs['p_label'])), key=lambda p: int(p[1:]))
        if not members:
            continue
        start = x
        for p in members:
            positions[p] = x
            x += 1.0
        group_spans[cat] = (start, x - 1.0, len(members))
        x += GROUP_GAP
    return positions, group_spans


def plot_two_arm(summary, value_col, ylabel, title, filename, ylim_top=None):
    """Single axes. Marker shape = loading condition (square flexion, triangle
    extension); fill = surgery (solid laminectomy only, hollow lami + fusion).
    Each arm sits at its own x offset within the patient's slot; the line joins
    Flexion to Extension for that arm."""
    positions, group_spans = x_layout()
    ymax = summary[value_col].max()
    top = ylim_top if ylim_top is not None else (ymax * 1.15 if ymax > 0 else 1)

    fig, ax = plt.subplots(figsize=(10, 5.5))

    for p in fusion_participants:
        if p in positions:
            ax.axvspan(positions[p] - 0.5, positions[p] + 0.5,
                       color=FUSION_SPAN_COLOR, alpha=0.2, zorder=0)

    for arm in ARMS:
        dx = ARM_X_OFFSET[arm]
        arm_sub = summary[summary['arm'] == arm]
        # Line joins Flexion to Extension for the same patient and surgery.
        for p_label, grp in arm_sub.groupby('p_label'):
            if len(grp) == 2:
                ax.vlines(positions[p_label] + dx, grp[value_col].min(), grp[value_col].max(),
                          color=MARKER_COLOR, linewidth=1.0, alpha=0.5, zorder=2)
        for condition in LOADING_CONDITIONS:
            sub = arm_sub[arm_sub['loading_condition'] == condition]
            xs = [positions[p] + dx for p in sub['p_label']]
            marker = CONDITION_MARKERS[condition]
            if ARM_FILLED[arm]:
                ax.scatter(xs, sub[value_col], color=MARKER_COLOR, marker=marker, s=60, zorder=3)
            else:
                ax.scatter(xs, sub[value_col], facecolors='none', edgecolors=MARKER_COLOR,
                           marker=marker, s=60, linewidths=1.4, zorder=3)

    # Category group separators and labels
    for cat, (x0, x1, n) in group_spans.items():
        ax.text((x0 + x1) / 2, 1.0, '{} (n={})'.format(CATEGORY_NAMES[cat], n),
                transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=10)
        if cat != min(group_spans):
            ax.axvline(x0 - (GROUP_GAP + 1.0) / 2, color='0.7', linestyle=':', linewidth=0.8, zorder=1)

    ordered = sorted(positions, key=lambda p: positions[p])
    ax.set_xticks([positions[p] for p in ordered])
    ax.set_xticklabels(ordered)
    ax.set_xlim(-0.7, max(positions.values()) + 0.7)
    ax.set_ylim(0, top)
    ax.set_xlabel('Participant (grouped by {} category)'.format(CATEGORY_SOURCE))
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=22)

    handles = [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')]
    handles += [Line2D([0], [0], marker=CONDITION_MARKERS[c], linestyle='', color='black',
                       label=c) for c in LOADING_CONDITIONS]
    handles.append(Line2D([0], [0], linestyle='none', marker='None', label=''))
    handles.append(Line2D([0], [0], linestyle='none', marker='None', label='Surgery'))
    handles += [Line2D([0], [0], marker='o', linestyle='', color='black',
                       markerfacecolor='black' if ARM_FILLED[a] else 'none', label=ARM_LABELS[a])
                for a in ARMS]
    handles.append(Line2D([0], [0], linestyle='none', marker='None', label=''))
    handles.append(Patch(facecolor=FUSION_SPAN_COLOR, alpha=0.2, label='Fusion performed (ground truth)'))
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()

    path = os.path.join(OUT_DIR, filename)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    print("Plot saved: {}".format(path))


# ============================================================
# PART A: % cord volume above threshold
# ============================================================
print("")
print("=" * 70)
print("PART A: % cord volume above threshold, laminectomy vs lami+fusion")
print("=" * 70)

cord_records = []
for (participant, condition, arm), grp in cord_df.groupby(['participant', 'loading_condition', 'arm']):
    for threshold in CORD_THRESHOLDS:
        cord_records.append({
            'participant': participant, 'loading_condition': condition, 'arm': arm,
            'threshold': threshold, 'pct_above': pct_volume_above(grp, threshold)})
cord_summary = attach_info(pd.DataFrame(cord_records))
cord_summary['is_factual'] = [
    bool(job_info.loc[(r.participant, r.loading_condition, r.arm), 'is_factual'])
    for r in cord_summary.itertuples()]
cord_summary_path = os.path.join(OUT_DIR, 'fusion_vs_lami_cord_summary.csv')
cord_summary.to_csv(cord_summary_path, index=False)
print("Summary saved: {}".format(cord_summary_path))

# One shared y-range across thresholds would hide the small ones, so each
# threshold is scaled to its own data; the primary threshold is the headline plot.
for threshold in CORD_THRESHOLDS:
    t_sub = cord_summary[cord_summary['threshold'] == threshold]
    tag = 't{}'.format(('%g' % threshold).replace('.', 'p'))
    suffix = '' if threshold == PRIMARY_THRESHOLD else '_sensitivity'
    plot_two_arm(
        t_sub, 'pct_above', '% cord volume\nabove threshold',
        'Cord strain, MPS ≥ {:g}: laminectomy vs laminectomy + fusion (PostOp, no preload)'.format(threshold),
        'fusion_vs_lami_cord_plot_{}{}.pdf'.format(tag, suffix))

print("")
print("Cord % volume above {:g} (paired, per patient):".format(PRIMARY_THRESHOLD))
_piv = (cord_summary[cord_summary['threshold'] == PRIMARY_THRESHOLD]
        .pivot_table(index=['category', 'p_label', 'loading_condition'], columns='arm', values='pct_above')
        .reset_index())
_piv['fusion_minus_lami'] = _piv['fusion'] - _piv['lami']
print(_piv.sort_values(['loading_condition', 'category', 'p_label']).to_string(index=False, float_format='%.3f'))

# ============================================================
# PART B: % IVD volume above T95
# Thresholds are pooled so that the two surgical arms of one patient are
# judged against the SAME cutoff (otherwise the paired comparison is moot):
#   global:           one cutoff pooled over all 48 PostOp jobs
#   patient-specific: one cutoff per patient, pooled over that patient's
#                     4 jobs (lami + fusion x flexion + extension)
# ============================================================
print("")
print("=" * 70)
print("PART B: % IVD volume above T95, laminectomy vs lami+fusion")
print("=" * 70)

ivd_global_threshold = volume_weighted_percentile(ivd_df, p=IVD_PERCENTILE)
ivd_patient_threshold = {
    'P{}'.format(p): volume_weighted_percentile(grp, p=IVD_PERCENTILE)
    for p, grp in ivd_df.groupby('participant')}
print("Global T95 = {:.4f}".format(ivd_global_threshold))

ivd_records = []
for (participant, condition, arm), grp in ivd_df.groupby(['participant', 'loading_condition', 'arm']):
    p_label = 'P{}'.format(participant)
    ivd_records.append({
        'participant': participant, 'loading_condition': condition, 'arm': arm,
        'threshold_global': ivd_global_threshold,
        'threshold_patient': ivd_patient_threshold[p_label],
        'pct_above_global': pct_volume_above(grp, ivd_global_threshold),
        'pct_above_patient': pct_volume_above(grp, ivd_patient_threshold[p_label]),
        'peak_ivd_mps': grp['mps'].max()})
ivd_summary = attach_info(pd.DataFrame(ivd_records))
ivd_summary_path = os.path.join(OUT_DIR, 'fusion_vs_lami_ivd_summary.csv')
ivd_summary.to_csv(ivd_summary_path, index=False)
print("Summary saved: {}".format(ivd_summary_path))

plot_two_arm(
    ivd_summary, 'pct_above_global', '% IVD volume\nabove threshold',
    'IVD strain, global T95 = {:.4f}: laminectomy vs laminectomy + fusion (PostOp, no preload)'.format(
        ivd_global_threshold),
    'fusion_vs_lami_ivd_plot_t95_global.pdf')
plot_two_arm(
    ivd_summary, 'pct_above_patient', '% IVD volume\nabove threshold',
    'IVD strain, patient-specific T95: laminectomy vs laminectomy + fusion (PostOp, no preload)',
    'fusion_vs_lami_ivd_plot_t95_patientwise.pdf')

print("")
print("Patient-specific T95 values (MPS):")
print(pd.Series(ivd_patient_threshold).sort_index(key=lambda s: s.str[1:].astype(int))
      .to_string(float_format='%.4f'))

# ============================================================
# PART C: blob-size distribution (connected components of elements at or
# above threshold; r = equivalent-sphere radius). One grid per surgery arm,
# one panel per threshold; colour = patient, line style = loading condition.
# ============================================================
print("")
print("=" * 70)
print("PART C: blob-size distribution, laminectomy vs lami+fusion")
print("=" * 70)

adjacency_cache_path = os.path.join(CACHE_DIR, 'cache_fusion_blob_adjacency.csv')
if os.path.isfile(adjacency_cache_path):
    print("Loading cached blob adjacency: {}".format(adjacency_cache_path))
    adjacency_df = pd.read_csv(adjacency_cache_path)
else:
    print("No blob adjacency cache - building from '_topology.csv' files (slow, one-time)...")
    adjacency_parts, missing_topology = [], []
    for _, row in jobs.iterrows():
        topology_path = str(row['csv_path']).strip().replace('_mps.csv', '_topology.csv')
        if not os.path.isfile(topology_path):
            missing_topology.append((row['p_label'], row['loading_condition'], row['arm']))
            continue
        edges = load_adjacency_edges(topology_path)
        part = pd.DataFrame(edges, columns=['elem_a', 'elem_b'])
        part.insert(0, 'arm', row['arm'])
        part.insert(0, 'loading_condition', row['loading_condition'])
        part.insert(0, 'participant', int(row['participant']))
        adjacency_parts.append(part)
    if missing_topology:
        raise SystemExit("Missing '_topology.csv' for {} job(s) - re-run Alex_results_extraction.py:\n  {}".format(
            len(missing_topology), "\n  ".join(str(m) for m in missing_topology)))
    adjacency_df = pd.concat(adjacency_parts, ignore_index=True)
    adjacency_df.to_csv(adjacency_cache_path, index=False)
    print("Cached blob adjacency: {}".format(adjacency_cache_path))

job_edges = {key: list(zip(g['elem_a'], g['elem_b']))
             for key, g in adjacency_df.groupby(['participant', 'loading_condition', 'arm'])}
job_elements = {key: g[['element_label', 'mps', 'volume']].reset_index(drop=True)
                for key, g in cord_df.groupby(['participant', 'loading_condition', 'arm'])}
job_total_volume = {key: g['volume'].sum() for key, g in job_elements.items()}

blob_records = []
for threshold_name, threshold_val in BLOB_THRESHOLDS.items():
    for key, edges in job_edges.items():
        participant, condition, arm = key
        elems = job_elements[key]
        vol_lookup = dict(zip(elems['element_label'], elems['volume']))
        exceeding = set(elems.loc[elems['mps'] >= threshold_val, 'element_label'])
        if not exceeding:
            continue
        for blob in find_blobs(exceeding, edges):
            volume = sum(vol_lookup[lbl] for lbl in blob)
            blob_records.append({
                'participant': 'P{}'.format(participant), 'loading_condition': condition,
                'arm': arm, 'threshold': threshold_name, 'n_elements': len(blob),
                'volume': volume, 'r': (3.0 * volume / (4.0 * math.pi)) ** (1.0 / 3.0)})
blob_df = pd.DataFrame(blob_records)
blob_summary_path = os.path.join(OUT_DIR, 'fusion_vs_lami_blob_summary.csv')
blob_df.to_csv(blob_summary_path, index=False)
print("Blob summary saved: {}".format(blob_summary_path))

blob_stats = blob_df.groupby(['participant', 'loading_condition', 'arm', 'threshold']).agg(
    n_blobs=('r', 'size'), median_r=('r', 'median'), max_r=('r', 'max'),
    total_volume=('volume', 'sum'), max_blob_volume=('volume', 'max')).reset_index()
blob_stats['pct_volume_in_largest_blob'] = 100.0 * blob_stats['max_blob_volume'] / blob_stats['total_volume']
blob_stats_path = os.path.join(OUT_DIR, 'fusion_vs_lami_blob_stats.csv')
blob_stats.to_csv(blob_stats_path, index=False)
print("Blob stats saved: {}".format(blob_stats_path))

# Same axes for both arms and all thresholds, so every panel is directly comparable.
blob_xlim = (blob_df['r'].min() * 0.9, blob_df['r'].max() * 1.1)
blob_ylim = (0, 102)

blob_patients = sorted(jobs['p_label'].unique(),
                       key=lambda p: (category_by_participant[p], int(p[1:])))
blob_patient_colors = {p: plt.cm.tab20(i % 20) for i, p in enumerate(blob_patients)}


def blob_plot_threshold(ax, threshold_name, arm):
    sub = blob_df[(blob_df['threshold'] == threshold_name) & (blob_df['arm'] == arm)]
    for (p_label, condition), grp in sub.groupby(['participant', 'loading_condition']):
        grp_sorted = grp.sort_values('r')
        total_vol = job_total_volume[(int(p_label[1:]), condition, arm)]
        cum_pct = list(100.0 * grp_sorted['volume'].cumsum() / total_vol)
        # Extend to the shared axis edges so the curve doesn't float mid-panel.
        rs = [blob_xlim[0]] + list(grp_sorted['r']) + [blob_xlim[1]]
        ys = [0.0] + cum_pct + [cum_pct[-1]]
        ax.plot(rs, ys, color=blob_patient_colors[p_label],
                linestyle=BLOB_CONDITION_LINESTYLES[condition], linewidth=1.2, drawstyle='steps-post')
    ax.set_xscale('log')
    ax.set_xlim(blob_xlim)
    ax.set_ylim(blob_ylim)
    ax.xaxis.set_major_locator(LogLocator(base=10.0))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_title('Threshold = {:g}'.format(BLOB_THRESHOLDS[threshold_name]))
    ax.set_xlabel('MPS concentration effective radius r (mm)')
    ax.set_ylabel('Cumulative % of total cord volume\nabove MPS threshold')


blob_patient_handles = [
    Line2D([0], [0], color=blob_patient_colors[p], linestyle='-',
           label='{} ({}{})'.format(p, CATEGORY_NAMES[category_by_participant[p]],
                                    ', fused' if p in fusion_participants else ''))
    for p in blob_patients]
blob_condition_handles = [Line2D([0], [0], color='black', linestyle=ls, label=c)
                          for c, ls in BLOB_CONDITION_LINESTYLES.items()]

for arm in ARMS:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax, threshold_name in zip(axes.flat, BLOB_THRESHOLDS):
        blob_plot_threshold(ax, threshold_name, arm)
    fig.suptitle(ARM_LABELS[arm])
    fig.legend(handles=blob_patient_handles + blob_condition_handles,
               loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, 'fusion_vs_lami_blob_distribution_{}_grid.pdf'.format(arm))
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    print("Plot saved: {}".format(path))

# ============================================================
# PART D: statistics
#
# H0: adding fusion does not change % cord volume above threshold (the paired
# difference d = fusion - lami is symmetric about 0). Flexion and Extension
# are analysed completely separately - nothing is pooled across loadings.
#
# Primary:    exact Wilcoxon signed-rank (Pratt zero handling), base R
# Effect:     mean difference (exact paired permutation test + CI), Hodges-
#             Lehmann median difference, exact sign test
# Supporting: Kenward-Roger LME (lmerTest) with QQ plots + Shapiro-Wilk,
#             pct ~ surgery (+ category) + (1 | patient); paired t-test
#             fallback if the random intercept is singular (as DCM_main.py)
# Exploratory: surgery x category interaction; factual-only fusion vs lami
#             (confounded - surgeons chose who got fusion)
#
# The signed-rank p-value is from coin::wilcoxsign_test (exact, Pratt) and is
# cross-checked against an independent base-R enumeration of the 2^m sign
# patterns (m = number of non-zero differences; at most 2^12 = 4096). The
# run aborts if the two ever disagree. The mean-difference permutation test,
# its CI and the two-sample tests are base-R enumeration only.
# ============================================================
print("")
print("=" * 70)
print("PART D: statistics, laminectomy vs lami+fusion")
print("=" * 70)

import numpy as np

# R setup - same as DCM_main.py
os.environ.setdefault('R_HOME', r'C:\Program Files\R\R-4.6.1')
os.environ.setdefault('R_LIBS_USER', os.path.join(os.path.expanduser('~'), 'Documents', 'R', 'win-library', '4.6'))
os.environ['PATH'] = os.path.join(os.environ['R_HOME'], 'bin', 'x64') + os.pathsep + os.environ['PATH']
# Git's sh.exe on PATH makes R's "CMD config" shell out to 'make' (not installed) and crashes rpy2.
os.environ['PATH'] = os.pathsep.join(
    p for p in os.environ['PATH'].split(os.pathsep) if 'Git' not in p)

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

importr('lme4')
importr('lmerTest')
importr('coin')

# Jobs known to be invalid, excluded from the stats and listed in the output.
# (patient label, loading) pairs, e.g. {('P5', 'Extension')}. Empty now: P5 PostOp Extension
# originally ran the Flexion loading (Step-2 rotation +0.35 instead of -0.175 in Job-017 /
# Job-517); both jobs were corrected, re-run and re-extracted.
EXCLUDE_PAIRS = set()
TIE_TOL = 0.01          # percentage points; |d| below this is a tie in the sign test
ZERO_TOL = 1e-9         # a pct value at or below this is "no volume above threshold"
NORM_Q75 = 0.6744897501960817   # qnorm(0.75), for the QQ reference line

DIAG_DIR = os.path.join(OUT_DIR, 'diagnostic_plots')
os.makedirs(DIAG_DIR, exist_ok=True)

ro.r('''
    signs_matrix <- function(m) {
        if (m == 0) return(matrix(0, nrow = 1, ncol = 0))
        as.matrix(expand.grid(rep(list(c(-1, 1)), m)))
    }

    # Exact paired sign-flip permutation test on the MEAN difference.
    # Zero differences contribute nothing under any sign flip but stay in the mean (divide by n).
    perm_mean_test <- function(d) {
        n <- length(d); nz <- d[abs(d) > 1e-9]; m <- length(nz)
        if (m == 0) return(c(mean_diff = mean(d), p = NA_real_, n = n, m = 0))
        S <- signs_matrix(m)
        t_null <- as.vector(S %*% nz) / n
        t_obs <- sum(nz) / n
        c(mean_diff = t_obs, p = mean(abs(t_null) >= abs(t_obs) - 1e-9), n = n, m = m)
    }

    # CI for the mean difference by inverting the sign-flip test: all shifts delta for which
    # d - delta is not rejected at alpha. Resolution is the grid spacing.
    perm_mean_ci <- function(d, alpha = 0.05, ngrid = 1501) {
        n <- length(d); S <- signs_matrix(n)
        a <- as.vector(S %*% d) / n; r <- rowSums(S) / n
        rng <- range(d); pad <- (rng[2] - rng[1]) + 1
        grid <- seq(rng[1] - pad, rng[2] + pad, length.out = ngrid)
        keep <- vapply(grid, function(delta) {
            t_null <- a - delta * r; t_obs <- mean(d) - delta
            mean(abs(t_null) >= abs(t_obs) - 1e-9) > alpha
        }, logical(1))
        if (!any(keep)) return(c(NA_real_, NA_real_))
        range(grid[keep])
    }

    # Exact Wilcoxon signed-rank, Pratt zero handling: |d| ranked INCLUDING zeros, zeros then
    # removed; null = all 2^m sign assignments of the remaining ranks (conditional exact).
    signrank_pratt_test <- function(d) {
        n <- length(d); a <- abs(d); zero <- a <= 1e-9
        rk <- rank(a, ties.method = "average")
        rnz <- rk[!zero]; m <- length(rnz)
        if (m == 0) return(c(W_plus = NA_real_, p = NA_real_, p_enum = NA_real_, n = n, m = 0))
        w_obs <- sum(rk[d > 1e-9])
        S <- signs_matrix(m)
        w_null <- as.vector(((S + 1) / 2) %*% rnz)
        mu <- sum(rnz) / 2
        p_enum <- mean(abs(w_null - mu) >= abs(w_obs - mu) - 1e-9)
        # Library implementation (reported); the enumeration above is an independent cross-check.
        p_coin <- tryCatch(as.numeric(coin::pvalue(coin::wilcoxsign_test(
            fusion ~ lami, data = data.frame(fusion = d, lami = 0),
            zero.method = "Pratt", distribution = coin::exact()))), error = function(e) NA_real_)
        c(W_plus = w_obs, p = p_coin, p_enum = p_enum, n = n, m = m)
    }

    # Hodges-Lehmann estimate: median of all Walsh averages (d_i + d_j)/2, i <= j.
    hl_estimate <- function(d) {
        w <- outer(d, d, "+") / 2
        median(w[upper.tri(w, diag = TRUE)])
    }

    # Exact sign test on differences larger than a tolerance (smaller ones are ties).
    sign_test <- function(d, tol) {
        neg <- sum(d < -tol); pos <- sum(d > tol); k <- neg + pos
        c(n_lower = neg, n_higher = pos, n_tied = length(d) - k,
          p = if (k == 0) NA_real_ else binom.test(neg, k, 0.5)$p.value)
    }

    # Two-sample exact permutation (fusion group vs lami group), mean and rank statistics.
    perm_two_sample <- function(x_f, x_l) {
        v <- c(x_f, x_l); nf <- length(x_f); n_all <- length(v)
        idx <- combn(n_all, nf)
        obs <- mean(x_f) - mean(x_l)
        null_m <- apply(idx, 2, function(ix) mean(v[ix]) - mean(v[-ix]))
        rk <- rank(v); obs_r <- sum(rk[seq_len(nf)]); mu <- nf * (n_all + 1) / 2
        null_r <- apply(idx, 2, function(ix) sum(rk[ix]))
        c(diff_mean = obs, p_mean = mean(abs(null_m) >= abs(obs) - 1e-9),
          p_rank = mean(abs(null_r - mu) >= abs(obs_r - mu) - 1e-9),
          n_fusion = nf, n_lami = n_all - nf)
    }

    get_coef_df <- function(model, ddf = "Kenward-Roger") {
        df <- as.data.frame(coef(summary(model, ddf = ddf)))
        df <- cbind(Term = rownames(df), df)
        rownames(df) <- NULL
        df
    }

    # pct_above ~ surgery (+ category) + (1 | patient), Kenward-Roger.
    lme_surgery <- function(dat, with_category) {
        dat$patient <- factor(dat$patient)
        dat$surgery <- factor(dat$surgery, levels = c("Laminectomy only", "Laminectomy + fusion"))
        dat$category <- factor(dat$category)
        f <- if (with_category) pct_above ~ surgery + category + (1 | patient)
             else pct_above ~ surgery + (1 | patient)
        m <- lmerTest::lmer(f, data = dat)
        cf <- get_coef_df(m)
        row <- cf[cf$Term == "surgeryLaminectomy + fusion", ]
        est <- row[["Estimate"]]; se <- row[["Std. Error"]]; dff <- row[["df"]]
        rs <- resid(m)
        sw <- tryCatch(shapiro.test(rs), error = function(e) NULL)
        qq <- qqnorm(rs, plot.it = FALSE)
        list(coef = cf, varcorr = as.data.frame(VarCorr(m)), singular = isSingular(m),
             est = est, se = se, df = dff, p = row[["Pr(>|t|)"]],
             ci_low = est - qt(0.975, dff) * se, ci_high = est + qt(0.975, dff) * se,
             sw_W = if (is.null(sw)) NA_real_ else unname(sw$statistic),
             sw_p = if (is.null(sw)) NA_real_ else sw$p.value,
             qq_x = qq$x, qq_y = qq$y)
    }

    # Exploratory: KR F-test for surgery x category in the full interaction model.
    lme_interaction_p <- function(dat) {
        dat$patient <- factor(dat$patient)
        dat$surgery <- factor(dat$surgery, levels = c("Laminectomy only", "Laminectomy + fusion"))
        dat$category <- factor(dat$category)
        tryCatch({
            m <- lmerTest::lmer(pct_above ~ surgery * category + (1 | patient), data = dat)
            a <- anova(m, ddf = "Kenward-Roger")
            c(p = a["surgery:category", "Pr(>F)"], singular = as.numeric(isSingular(m)))
        }, error = function(e) c(p = NA_real_, singular = NA_real_))
    }

    # Paired t-test fallback (what DCM_main.py does when the random intercept is singular).
    paired_t <- function(d) {
        tt <- t.test(d)
        qq <- qqnorm(d, plot.it = FALSE)
        sw <- tryCatch(shapiro.test(d), error = function(e) NULL)
        list(est = unname(tt$estimate), df = unname(tt$parameter), p = tt$p.value,
             ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
             sw_W = if (is.null(sw)) NA_real_ else unname(sw$statistic),
             sw_p = if (is.null(sw)) NA_real_ else sw$p.value,
             qq_x = qq$x, qq_y = qq$y)
    }
''')


def _named(vec):
    return dict(zip(list(vec.names), [float(v) for v in vec]))


def _to_r_df(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        return ro.conversion.py2rpy(df)


def _from_r_df(r_obj):
    with localconverter(ro.default_converter + pandas2ri.converter):
        return ro.conversion.rpy2py(r_obj)


def _fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return 'n/a'
    return '<0.001' if p < 0.001 else '{:.3f}'.format(p)


# Guard: a patient whose Flexion and Extension results are identical is a copy of one
# loading (this is how the P5 PostOp problem showed up).
_peak = cord_df.groupby(['participant', 'loading_condition', 'arm'])['mps'].max().unstack('loading_condition')
for (_p, _arm), _r in _peak.iterrows():
    if _r['Flexion'] > 0 and abs(_r['Flexion'] - _r['Extension']) < 1e-9:
        _flag = 'excluded' if ('P{}'.format(_p), 'Extension') in EXCLUDE_PAIRS else 'NOT excluded'
        print("WARNING: P{} ({}) has identical peak MPS in Flexion and Extension "
              "({:.5f}) - one loading is probably a copy [{}].".format(_p, _arm, _r['Flexion'], _flag))
print("Excluded from the statistics: {}".format(sorted(EXCLUDE_PAIRS) or 'none'))


def paired_table(condition, threshold):
    """One row per patient: lami, fusion, diff, category, ground-truth fusion."""
    sub = cord_summary[(cord_summary['loading_condition'] == condition) &
                       ((cord_summary['threshold'] - threshold).abs() < 1e-12)]
    wide = sub.pivot_table(index='p_label', columns='arm', values='pct_above')
    wide = wide.loc[[p for p in wide.index if (p, condition) not in EXCLUDE_PAIRS]].copy()
    wide['diff'] = wide['fusion'] - wide['lami']
    wide['category'] = [category_by_participant[p] for p in wide.index]
    wide['gt_fusion'] = [p in fusion_participants for p in wide.index]
    wide['both_zero'] = (wide['fusion'] <= ZERO_TOL) & (wide['lami'] <= ZERO_TOL)
    return wide.sort_values(['category', 'p_label'] if 'p_label' in wide else ['category'], kind='stable')


def run_paired_tests(wide):
    d = wide['diff'].values.astype(float)
    dv = ro.FloatVector(d)
    out = {'n': len(d), 'n_both_zero': int(wide['both_zero'].sum()),
           'n_differ': int((np.abs(d) > ZERO_TOL).sum()), 'diffs': d}
    out['signrank'] = _named(ro.r['signrank_pratt_test'](dv))
    sr = out['signrank']
    if out['n_differ'] > 0 and not (abs(sr['p'] - sr['p_enum']) < 1e-9):
        raise SystemExit("coin and base-R enumeration disagree on the signed-rank p-value "
                         "({} vs {}) - stop and investigate.".format(sr['p'], sr['p_enum']))
    out['perm_mean'] = _named(ro.r['perm_mean_test'](dv))
    out['hl'] = float(ro.r['hl_estimate'](dv)[0])
    out['sign'] = _named(ro.r['sign_test'](dv, ro.FloatVector([TIE_TOL])))
    if out['n_differ'] > 0:
        ci = list(ro.r['perm_mean_ci'](dv))
        out['mean_ci'] = (float(ci[0]), float(ci[1]))
    else:
        out['mean_ci'] = (None, None)
    out['occurrence'] = {
        'both_zero': int(wide['both_zero'].sum()),
        'lami_only': int(((wide['lami'] > ZERO_TOL) & (wide['fusion'] <= ZERO_TOL)).sum()),
        'fusion_only': int(((wide['fusion'] > ZERO_TOL) & (wide['lami'] <= ZERO_TOL)).sum()),
        'both_above': int(((wide['fusion'] > ZERO_TOL) & (wide['lami'] > ZERO_TOL)).sum()),
    }
    return out


def save_qq(xs, ys, path, title, sw_W, sw_p):
    xs, ys = list(xs), list(ys)
    q1, q3 = np.percentile(ys, 25), np.percentile(ys, 75)
    slope = (q3 - q1) / (2 * NORM_Q75)
    intercept = q1 + slope * NORM_Q75
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(xs, ys, color='#003f5c', s=20, alpha=0.8)
    lx = [min(xs), max(xs)]
    ax.plot(lx, [slope * x + intercept for x in lx], color='red')
    ax.set_xlabel('Theoretical Quantiles')
    ax.set_ylabel('Sample Quantiles')
    ax.set_title('{}\nShapiro-Wilk: W={}, p={}'.format(
        title, 'n/a' if np.isnan(sw_W) else '{:.3g}'.format(sw_W), _fmt_p(sw_p)))
    fig.tight_layout()
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    print("  QQ plot saved: {}".format(path))


def fit_surgery_model(wide, condition, with_category):
    """KR LME; paired t-test fallback if the random intercept is singular."""
    tag = 'with_category' if with_category else 'surgery_only'
    d = wide['diff'].values.astype(float)
    if np.all(np.abs(d) <= ZERO_TOL):
        return {'kind': 'degenerate', 'tag': tag, 'note': 'every paired difference is 0 - nothing to fit'}
    long_rows = []
    for p, r in wide.iterrows():
        long_rows.append({'patient': p, 'surgery': 'Laminectomy only',
                          'category': str(int(r['category'])), 'pct_above': float(r['lami'])})
        long_rows.append({'patient': p, 'surgery': 'Laminectomy + fusion',
                          'category': str(int(r['category'])), 'pct_above': float(r['fusion'])})
    r_df = _to_r_df(pd.DataFrame(long_rows))
    res = ro.r['lme_surgery'](r_df, ro.BoolVector([with_category]))
    if not bool(res.rx2('singular')[0]):
        out = {'kind': 'lmm', 'tag': tag, 'est': float(res.rx2('est')[0]), 'se': float(res.rx2('se')[0]),
               'df': float(res.rx2('df')[0]), 'p': float(res.rx2('p')[0]),
               'ci': (float(res.rx2('ci_low')[0]), float(res.rx2('ci_high')[0])),
               'sw': (float(res.rx2('sw_W')[0]), float(res.rx2('sw_p')[0])), 'sw_label': 'LME residuals',
               'qq': (list(res.rx2('qq_x')), list(res.rx2('qq_y'))),
               'coef_df': _from_r_df(res.rx2('coef')), 'varcorr_df': _from_r_df(res.rx2('varcorr'))}
    else:
        tt = ro.r['paired_t'](ro.FloatVector(d))
        out = {'kind': 'ttest', 'tag': tag, 'est': float(tt.rx2('est')[0]), 'se': None,
               'df': float(tt.rx2('df')[0]), 'p': float(tt.rx2('p')[0]),
               'ci': (float(tt.rx2('ci_low')[0]), float(tt.rx2('ci_high')[0])),
               'sw': (float(tt.rx2('sw_W')[0]), float(tt.rx2('sw_p')[0])), 'sw_label': 'paired differences',
               'qq': (list(tt.rx2('qq_x')), list(tt.rx2('qq_y'))), 'coef_df': None, 'varcorr_df': None}
    out['interaction_p'] = None
    save_qq(out['qq'][0], out['qq'][1],
            os.path.join(DIAG_DIR, 'qq_{}_{}.pdf'.format(condition.lower(), tag)),
            '{} ({}): {} Q-Q'.format(condition, tag.replace('_', ' '), out['sw_label']),
            out['sw'][0], out['sw'][1])
    return out


def interaction_p(wide):
    long_rows = []
    for p, r in wide.iterrows():
        for surgery, col in (('Laminectomy only', 'lami'), ('Laminectomy + fusion', 'fusion')):
            long_rows.append({'patient': p, 'surgery': surgery, 'category': str(int(r['category'])),
                              'pct_above': float(r[col])})
    res = _named(ro.r['lme_interaction_p'](_to_r_df(pd.DataFrame(long_rows))))
    return res['p'], res['singular']


# ------------------------------------------------------------
# Run: primary threshold in full, other thresholds as sensitivity
# ------------------------------------------------------------
stats = {}            # (condition, threshold) -> results
for condition in LOADING_CONDITIONS:
    for threshold in CORD_THRESHOLDS:
        wide = paired_table(condition, threshold)
        stats[(condition, threshold)] = {'wide': wide, 'tests': run_paired_tests(wide)}

model_fits = {}       # (condition, tag) -> fit, primary threshold only
interactions = {}
for condition in LOADING_CONDITIONS:
    wide = stats[(condition, PRIMARY_THRESHOLD)]['wide']
    print("")
    print("--- {} at {:g}: LME (Kenward-Roger) ---".format(condition, PRIMARY_THRESHOLD))
    for with_category in (False, True):
        model_fits[(condition, 'with_category' if with_category else 'surgery_only')] = \
            fit_surgery_model(wide, condition, with_category)
    interactions[condition] = interaction_p(wide)

# Secondary, confounded: factual-only fusion (4) vs lami (8)
factual_results = {}
for condition in LOADING_CONDITIONS:
    sub = cord_summary[(cord_summary['loading_condition'] == condition) & cord_summary['is_factual'] &
                       ((cord_summary['threshold'] - PRIMARY_THRESHOLD).abs() < 1e-12)]
    sub = sub[[(p, condition) not in EXCLUDE_PAIRS for p in sub['p_label']]]
    x_f = sub.loc[sub['arm'] == 'fusion', 'pct_above'].values.astype(float)
    x_l = sub.loc[sub['arm'] == 'lami', 'pct_above'].values.astype(float)
    factual_results[condition] = _named(ro.r['perm_two_sample'](ro.FloatVector(x_f), ro.FloatVector(x_l)))

# ------------------------------------------------------------
# Tables / CSV
# ------------------------------------------------------------
diff_rows, summary_rows = [], []
for (condition, threshold), s in stats.items():
    t = s['tests']
    for p_label, r in s['wide'].iterrows():
        diff_rows.append({'participant': p_label, 'loading_condition': condition, 'threshold': threshold,
                          'category': int(r['category']), 'gt_fusion': bool(r['gt_fusion']),
                          'pct_lami': r['lami'], 'pct_fusion': r['fusion'], 'diff_fusion_minus_lami': r['diff']})
    summary_rows.append({
        'loading_condition': condition, 'threshold': threshold, 'n_pairs': t['n'],
        'n_differ': t['n_differ'], 'n_both_zero': t['n_both_zero'],
        'mean_diff': t['perm_mean']['mean_diff'], 'mean_diff_ci_low': t['mean_ci'][0],
        'mean_diff_ci_high': t['mean_ci'][1], 'hodges_lehmann': t['hl'],
        'p_signrank_pratt_exact': t['signrank']['p'], 'p_perm_mean_exact': t['perm_mean']['p'],
        'n_lower': t['sign']['n_lower'], 'n_higher': t['sign']['n_higher'], 'n_tied': t['sign']['n_tied'],
        'p_sign_test': t['sign']['p']})
pd.DataFrame(diff_rows).to_csv(os.path.join(OUT_DIR, 'fusion_vs_lami_paired_differences.csv'), index=False)
stats_summary_df = pd.DataFrame(summary_rows)
stats_summary_df.to_csv(os.path.join(OUT_DIR, 'fusion_vs_lami_stats_summary.csv'), index=False)

# ------------------------------------------------------------
# Forest-style plot of the paired differences (descriptive, grouped by category)
# ------------------------------------------------------------
positions, group_spans = x_layout()
fig, axes = plt.subplots(1, len(LOADING_CONDITIONS), figsize=(13, 5), sharey=True)
all_d = np.concatenate([stats[(c, PRIMARY_THRESHOLD)]['wide']['diff'].values for c in LOADING_CONDITIONS])
ylim = (min(all_d.min() * 1.15, -1), max(all_d.max() * 1.15, 1))
for ax, condition in zip(axes, LOADING_CONDITIONS):
    wide = stats[(condition, PRIMARY_THRESHOLD)]['wide']
    for p in fusion_participants:
        if p in positions:
            ax.axvspan(positions[p] - 0.5, positions[p] + 0.5, color=FUSION_SPAN_COLOR, alpha=0.2, zorder=0)
    ax.axhline(0, color='black', linewidth=0.8, zorder=1)
    for p_label, r in wide.iterrows():
        if r['both_zero']:
            ax.scatter(positions[p_label], 0, facecolors='none', edgecolors='0.4', marker='o', s=55, zorder=3)
        else:
            ax.scatter(positions[p_label], r['diff'], color='#2e75b6', marker='o', s=55, zorder=3)
    for cat, (x0, x1, n) in group_spans.items():
        ax.text((x0 + x1) / 2, 1.0, '{} (n={})'.format(CATEGORY_NAMES[cat], n),
                transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=10)
        if cat != min(group_spans):
            ax.axvline(x0 - (GROUP_GAP + 1.0) / 2, color='0.7', linestyle=':', linewidth=0.8, zorder=1)
    ordered = sorted(positions, key=lambda p: positions[p])
    ax.set_xticks([positions[p] for p in ordered])
    ax.set_xticklabels(ordered)
    ax.set_xlim(-0.7, max(positions.values()) + 0.7)
    ax.set_ylim(ylim)
    ax.set_title(condition, pad=22)
    ax.set_xlabel('Participant (grouped by {} category)'.format(CATEGORY_SOURCE))
    excluded_here = sorted(p for p, c in EXCLUDE_PAIRS if c == condition)
    if excluded_here:
        ax.text(0.01, 0.02, 'Excluded (invalid model): {}'.format(', '.join(excluded_here)),
                transform=ax.transAxes, fontsize=9, color='firebrick')
axes[0].set_ylabel('Fusion \u2212 lami\n(% cord volume above {:g})'.format(PRIMARY_THRESHOLD))
handles = [Line2D([0], [0], marker='o', linestyle='', color='#2e75b6', label='Paired difference'),
           Line2D([0], [0], marker='o', linestyle='', markerfacecolor='none', markeredgecolor='0.4',
                  label='Both arms zero (difference 0)'),
           Patch(facecolor=FUSION_SPAN_COLOR, alpha=0.2, label='Fusion performed (ground truth)')]
fig.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False)
fig.suptitle('Paired difference, fusion minus laminectomy only (MPS \u2265 {:g})'.format(PRIMARY_THRESHOLD), y=1.04)
fig.tight_layout()
forest_path = os.path.join(OUT_DIR, 'fusion_vs_lami_paired_difference_plot.pdf')
fig.savefig(forest_path, bbox_inches='tight')
plt.close(fig)
print("Plot saved: {}".format(forest_path))

# ------------------------------------------------------------
# Paired boxplot, same style as DCM_main.py's whole-cord boxplot: laminectomy only vs
# laminectomy + fusion per patient. Lines join each patient's two values (orange = patient
# actually had fusion); the box median is the summary (no mean marker - the primary test is
# rank-based and the mean is driven by single patients); the p-value is the primary
# exact signed-rank test (each condition is its own hypothesis, tested once, no multiplicity adjustment).
# ------------------------------------------------------------
from matplotlib.colors import to_rgba

BOX_NAVY, BOX_TEAL = '#003f5c', '#58a4b0'


def paired_boxplot(wides, tests, ylabel, title, path, bracket_frac, headroom=1.25):
    """wides[condition]: DataFrame (lami, fusion, gt_fusion) indexed by patient;
    tests[condition]: dict(p, n, n_differ) from the primary exact signed-rank test."""
    ymax = max(w[['lami', 'fusion']].values.max() for w in wides.values()) * headroom
    fig, axes = plt.subplots(1, len(LOADING_CONDITIONS), figsize=(10, 5.5), sharey=True)
    for ax, condition in zip(axes, LOADING_CONDITIONS):
        wide, t = wides[condition], tests[condition]
        lami_vals, fusion_vals = wide['lami'].values, wide['fusion'].values
        box = ax.boxplot([lami_vals, fusion_vals], positions=[0, 1], widths=0.35,
                         showfliers=False, patch_artist=True, zorder=2)
        for patch in box['boxes']:
            patch.set_facecolor(to_rgba(BOX_TEAL, 0.4))
            patch.set_edgecolor('black')
            patch.set_linewidth(0.5)
        for part in ('whiskers', 'caps', 'medians'):
            for line in box[part]:
                line.set_color('black')
                line.set_linewidth(0.5)
        for p_label, r in wide.iterrows():
            ax.plot([0, 1], [r['lami'], r['fusion']],
                    color=FUSION_SPAN_COLOR if r['gt_fusion'] else 'grey',
                    linewidth=1.0 if r['gt_fusion'] else 0.6, alpha=0.9 if r['gt_fusion'] else 0.5, zorder=2.5)
        marker = CONDITION_MARKERS[condition]
        ax.scatter([0] * len(lami_vals), lami_vals, color=BOX_NAVY, s=25, alpha=0.7, marker=marker, zorder=3)
        ax.scatter([1] * len(fusion_vals), fusion_vals, color=BOX_NAVY, s=25, alpha=0.7, marker=marker, zorder=3)

        p = t['p']
        p_text = 'n/a' if (p is None or np.isnan(p)) else ('p < 0.001' if p < 0.001 else 'p = {:.3f}'.format(p))
        bracket_y, tick = bracket_frac[condition] * ymax, 0.03 * ymax
        ax.plot([0, 0, 1, 1], [bracket_y - tick, bracket_y, bracket_y, bracket_y - tick],
                color='black', linewidth=1.2, zorder=6)
        ax.text(0.5, bracket_y + 0.01 * ymax, p_text, ha='center', va='bottom', fontsize=10)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Laminectomy\nOnly', 'Laminectomy\n& Fusion'])
        ax.set_xlabel(condition)
        # Pair counts as sub-text under the condition name (separate annotation so only it is smaller).
        ax.annotate('{} pairs: {} differ, {} tied at zero'.format(t['n'], t['n_differ'], t['n'] - t['n_differ']),
                    xy=(0.5, 0), xycoords=ax.xaxis.label, xytext=(0, -4), textcoords='offset points',
                    ha='center', va='top', fontsize=10)
        ax.set_ylim(-0.05 * ymax, ymax)
        ax.set_xlim(-0.5, 1.5)
    axes[0].set_ylabel(ylabel)
    handles = [
        Line2D([0], [0], marker=CONDITION_MARKERS['Flexion'], linestyle='', markerfacecolor=BOX_NAVY,
               markeredgecolor=BOX_NAVY, alpha=0.7, label='Flexion'),
        Line2D([0], [0], marker=CONDITION_MARKERS['Extension'], linestyle='', markerfacecolor=BOX_NAVY,
               markeredgecolor=BOX_NAVY, alpha=0.7, label='Extension'),
        Line2D([0], [0], color=FUSION_SPAN_COLOR, linewidth=1.2, label='Fusion Performed'),
    ]
    axes[-1].legend(handles=handles, loc='upper right', frameon=False, fontsize=10)
    fig.suptitle(title)
    fig.tight_layout()
    fig.text(0.5, -0.03,
             'p-values: exact Wilcoxon signed-rank test (Pratt zero handling) on the paired patient\n'
             'differences (fusion − lami); flexion and extension tested separately.',
             ha='center', va='top', fontsize=10)
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    print("Plot saved: {}".format(path))


# Cord: % volume above MPS 0.015
paired_boxplot(
    {c: stats[(c, PRIMARY_THRESHOLD)]['wide'] for c in LOADING_CONDITIONS},
    {c: {'p': stats[(c, PRIMARY_THRESHOLD)]['tests']['signrank']['p'],
         'n': stats[(c, PRIMARY_THRESHOLD)]['tests']['n'],
         'n_differ': stats[(c, PRIMARY_THRESHOLD)]['tests']['n_differ']} for c in LOADING_CONDITIONS},
    '% of cord volume above MPS {:g}'.format(PRIMARY_THRESHOLD),
    '% cord volume above threshold: laminectomy only vs laminectomy with fusion (exact signed-rank test)',
    os.path.join(OUT_DIR, 'fusion_vs_lami_boxplot.pdf'),
    {'Flexion': 0.90, 'Extension': 0.50})

# IVD: the same primary test and plot for four thresholds, all from the 48 PostOp models (Part B):
# global T95, patient-specific T95 (per patient, pooled over their 4 jobs), global T90 and global T99.
ivd_groups = dict(list(ivd_df.groupby(['participant', 'loading_condition', 'arm'])))
IVD_VARIANTS = [   # (column in ivd_summary, y-axis threshold label, title threshold label, file tag)
    ('pct_above_global', 'global T95 ({:.4f})'.format(ivd_global_threshold), 'global T95', 't95_global'),
    ('pct_above_patient', 'patient-specific T95', 'patient-specific T95', 't95_patientwise'),
]
for q, tag in ((0.90, 't90'), (0.99, 't99')):
    thr = volume_weighted_percentile(ivd_df, p=q)
    col = 'pct_above_global_' + tag
    ivd_summary[col] = [pct_volume_above(ivd_groups[(r.participant, r.loading_condition, r.arm)], thr)
                        for r in ivd_summary.itertuples()]
    IVD_VARIANTS.append((col, 'global {} ({:.4f})'.format(tag.upper(), thr), 'global ' + tag.upper(),
                         tag + '_global'))

for col, thr_label, title_label, tag in IVD_VARIANTS:
    ivd_wides, ivd_tests = {}, {}
    for condition in LOADING_CONDITIONS:
        sub = ivd_summary[(ivd_summary['loading_condition'] == condition) &
                          ~ivd_summary['p_label'].isin([p for p, c in EXCLUDE_PAIRS if c == condition])]
        w = sub.pivot_table(index='p_label', columns='arm', values=col)
        w['gt_fusion'] = [p in fusion_participants for p in w.index]
        d = (w['fusion'] - w['lami']).values.astype(float)
        sr = _named(ro.r['signrank_pratt_test'](ro.FloatVector(d)))
        if not np.isnan(sr['p']) and abs(sr['p'] - sr['p_enum']) > 1e-9:
            raise SystemExit("coin and base-R enumeration disagree for IVD {} {} ({} vs {})".format(
                title_label, condition, sr['p'], sr['p_enum']))
        ivd_wides[condition] = w
        ivd_tests[condition] = {'p': sr['p'], 'n': len(d), 'n_differ': int((np.abs(d) > ZERO_TOL).sum())}
        print("IVD {} {}: {} pairs, {} differ; signed-rank p = {}".format(
            title_label, condition, len(d), ivd_tests[condition]['n_differ'], _fmt_p(sr['p'])))
    paired_boxplot(
        ivd_wides, ivd_tests,
        '% of IVD volume above {}'.format(thr_label),
        '% IVD volume above {}: laminectomy only vs laminectomy with fusion (exact signed-rank test)'.format(
            title_label),
        os.path.join(OUT_DIR, 'fusion_vs_lami_ivd_boxplot_{}.pdf'.format(tag)),
        {'Flexion': 0.76, 'Extension': 0.76}, headroom=1.5)

# ------------------------------------------------------------
# Markdown write-up
# ------------------------------------------------------------
md = ["### Laminectomy only vs laminectomy + fusion (PostOp, no preload)", ""]
md.append("**H0:** adding fusion does not change % cord volume above MPS = {:g}; the paired difference "
          "d = fusion \u2212 lami is symmetric about 0. Flexion and Extension are analysed separately "
          "(not pooled). Category source: {}.".format(PRIMARY_THRESHOLD, CATEGORY_SOURCE))
md.append("")
md.append("Primary test: exact Wilcoxon signed-rank with Pratt zero handling. Flexion and Extension are separate "
          "hypotheses, each tested once, so no multiplicity adjustment is applied. The other thresholds and all "
          "other tests are supporting or sensitivity analyses.")
md.append("")
if EXCLUDE_PAIRS:
    md.append("**Excluded (invalid model, to be re-run):** {}.".format(
        ', '.join('{} {}'.format(p, c) for p, c in sorted(EXCLUDE_PAIRS))))
    md.append("")

for condition in LOADING_CONDITIONS:
    s = stats[(condition, PRIMARY_THRESHOLD)]
    t = s['tests']
    md.append("#### {} (threshold = {:g})".format(condition, PRIMARY_THRESHOLD))
    md.append("")
    md.append("**{} pairs; {} differ, {} are tied at zero (both arms have no volume above threshold).**".format(
        t['n'], t['n_differ'], t['n_both_zero']))
    md.append("")
    md.append("| Pair type | n |\n|---|---|\n| Both arms zero | {both_zero} |\n| Lami only above threshold | {lami_only} |\n"
              "| Fusion only above threshold | {fusion_only} |\n| Both above threshold | {both_above} |".format(**t['occurrence']))
    md.append("")
    ci = t['mean_ci']
    ci_txt = 'n/a' if ci[0] is None else '{:.2f} to {:.2f}'.format(ci[0], ci[1])
    md.append("| Test | Estimate | p | Note |")
    md.append("|---|---|---|---|")
    md.append("| **Wilcoxon signed-rank, exact, Pratt (primary)** | W+ = {} | **{}** | m = {} non-zero |".format(
        'n/a' if np.isnan(t['signrank']['W_plus']) else '{:.1f}'.format(t['signrank']['W_plus']),
        _fmt_p(t['signrank']['p']), int(t['signrank']['m'])))
    md.append("| Permutation test, mean difference | {:+.2f} pp (95% CI {}) | {} | exact, 2^{} sign patterns |".format(
        t['perm_mean']['mean_diff'], ci_txt, _fmt_p(t['perm_mean']['p']), int(t['perm_mean']['m'])))
    md.append("| Hodges-Lehmann median difference | {:+.2f} pp | | |".format(t['hl']))
    md.append("| Sign test (\\|d\\| > {:g} pp) | lower {:.0f}, higher {:.0f}, tied {:.0f} | {} | |".format(
        TIE_TOL, t['sign']['n_lower'], t['sign']['n_higher'], t['sign']['n_tied'], _fmt_p(t['sign']['p'])))
    for tag, label in (('surgery_only', 'surgery + (1 \\| patient)'),
                       ('with_category', 'surgery + category + (1 \\| patient)')):
        f = model_fits[(condition, tag)]
        if f['kind'] == 'degenerate':
            md.append("| KR LME: {} | n/a | n/a | {} |".format(label, f['note']))
            continue
        kind = 'KR LME' if f['kind'] == 'lmm' else 'paired t-test (LME singular)'
        md.append("| {}: {} | {:+.2f} pp (95% CI {:.2f} to {:.2f}) | {} | Shapiro-Wilk ({}) W = {:.3f}, p = {}{} |".format(
            kind, label, f['est'], f['ci'][0], f['ci'][1], _fmt_p(f['p']), f['sw_label'],
            f['sw'][0], _fmt_p(f['sw'][1]), ' (deviates from normality)' if f['sw'][1] < 0.05 else ''))
    ip, isingular = interactions[condition]
    md.append("| Exploratory: surgery \u00d7 category (KR F-test) | | {} | {} |".format(
        _fmt_p(ip), 'singular fit' if isingular == 1 else ''))
    md.append("")
    fr = factual_results[condition]
    md.append("Secondary, **confounded** (surgeons chose who got fusion): factual-only fusion (n={:.0f}) vs "
              "lami (n={:.0f}), mean difference {:+.2f} pp, exact permutation p = {} (mean), {} (ranks).".format(
                  fr['n_fusion'], fr['n_lami'], fr['diff_mean'], _fmt_p(fr['p_mean']), _fmt_p(fr['p_rank'])))
    md.append("")

md.append("#### Sensitivity to threshold")
md.append("")
md.append("The thresholds are one question measured four ways, not four tests; 0.015 was fixed as primary in advance. "
          "Higher thresholds produce more tied-at-zero pairs.")
md.append("")
md.append("| Condition | Threshold | Differ / pairs | Mean diff (pp) | Signed-rank p | Permutation p | Sign test p |")
md.append("|---|---|---|---|---|---|---|")
for _, r in stats_summary_df.iterrows():
    md.append("| {} | {:g} | {}/{} | {:+.2f} | {} | {} | {} |".format(
        r['loading_condition'], r['threshold'], int(r['n_differ']), int(r['n_pairs']), r['mean_diff'],
        _fmt_p(r['p_signrank_pratt_exact']), _fmt_p(r['p_perm_mean_exact']), _fmt_p(r['p_sign_test'])))
md.append("")
md.append("#### Interpretation notes")
md.append("")
md.append("- Pairs tied at zero stay in the mean but carry no sign information, so the exact tests can only reach "
          "p = 2/2^m (m = number of differing pairs). Subgroups of 3 to 5 patients can never reach p < 0.05; category "
          "effects are shown descriptively (paired-difference plot) and the interaction is exploratory.")
md.append("- These are deterministic simulations: p-values describe this cohort, not new patients.")
md.append("- Exactly-zero strain in many PostOp models is physically plausible (a decompressed cord with no preload "
          "need not touch bone or ligament) but has not been confirmed for every job.")
md_path = os.path.join(OUT_DIR, 'summary_fusion_vs_lami.md')
with open(md_path, 'w', encoding='utf-8') as f:
    f.write("\n".join(md) + "\n")
print("Summary saved: {}".format(md_path))
sys.stdout.reconfigure(encoding='utf-8', errors='replace')   # cp1252 consoles cannot print the minus sign
print("\n".join(md))
