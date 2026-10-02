"""
Does high-strain volume sit in GM or WM, and/or at the GM/WM boundary?

One number per (tissue, region) combo - % of that subset's own volume
above threshold - for GM x Interior, GM x Boundary, WM x Interior,
WM x Boundary. Directly comparable, no combined score.

Runs for PreOp (Flexion/Extension) and Oscillation (cumulative), one plot
per threshold. Needs '_mps_GM_WM.csv' and '_topology.csv' per job; missing
files are skipped, not a crash. Per-element data is cached per dataset -
delete cache_gm_wm_boundary_*.csv to rebuild.

Run: python plot_gm_wm_boundary_enrichment.py
"""

import os
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgba

from mps_common import PLOT_STYLE

# ============================================================
# USER SETTINGS
# ============================================================
ID_MAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'id_map.csv')
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

BLOB_FACE_SHARING_MIN_NODES = 4   # >=4 shared nodes approximates a shared C3D8 hex face

CONDITION_MARKERS = {'flexion': 's', 'extension': '^'}

PREOP_THRESHOLDS = {'t0p10': 0.10, 't0p15': 0.15}
OSC_THRESHOLDS = {'t0p01': 0.01, 't0p02': 0.02, 't0p03': 0.03, 't0p04': 0.04, 't0p05': 0.05,
                   't0p10': 0.10, 't0p15': 0.15}
# ============================================================

plt.rcParams.update(PLOT_STYLE)


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


def reduce_to_peak(df):
    peak_mps = df.groupby('element_label')['mps'].max()
    volume = df.groupby('element_label')['volume'].first()
    return pd.DataFrame({'mps': peak_mps, 'volume': volume}).reset_index()


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
    job_df = reduce_to_peak(raw)
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


id_map = pd.read_csv(ID_MAP_PATH, skipinitialspace=True)

preop_mjoa_by_participant = {}
for _, row in id_map.iterrows():
    if str(row.get('State', '')).strip().lower() != 'preop':
        continue
    p_label = 'P{}'.format(int(row['participant']))
    preop_mjoa_by_participant[p_label] = row.get('mJOA', '')


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


def load_or_build_combined_df(rows, tag, id_fn):
    """Pooled per-element data for this dataset, cached to
    'cache_gm_wm_boundary_<tag>.csv'. Delete the cache to force a rebuild."""
    id_cols = list(id_fn(rows.iloc[0]).keys()) if not rows.empty else []
    cache_path = os.path.join(OUT_DIR, 'cache_gm_wm_boundary_{}.csv'.format(tag))

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


def run_dataset(rows, thresholds, tag, title_prefix, has_condition, id_fn):
    combined_df, id_cols = load_or_build_combined_df(rows, tag, id_fn)

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
        path = os.path.join(OUT_DIR, 'multipatient_gm_wm_tissue_boundary_{}.csv'.format(tag))
        summary.to_csv(path, index=False)
        print(summary.to_string(index=False))
        print("Summary saved: {}".format(path))
        for name, val in thresholds.items():
            plot_path = os.path.join(OUT_DIR, 'multipatient_gm_wm_tissue_boundary_{}_{}_sortedbypreopmJOA.pdf'.format(tag, name))
            plot_tissue_boundary(summary[summary['threshold'] == name], val, plot_path,
                                  '{}: high-strain volume by tissue and region (threshold {:.2f})'.format(title_prefix, val),
                                  has_condition)
    else:
        print("No data available yet.")

    return combined_df, summary


print("=== PreOp (with preload) ===")
preop_rows = id_map[
    (id_map['State'].astype(str).str.strip().str.lower() == 'preop') &
    (id_map['loading_condition'].astype(str).str.strip().str.lower().isin(['flexion', 'extension']))
]
preop_elements, preop_summary = run_dataset(
    preop_rows, PREOP_THRESHOLDS, 'preop', 'PreOp (with preload)',
    has_condition=True,
    id_fn=lambda row: {'participant': 'P{}'.format(int(row['participant'])),
                        'loading_condition': str(row['loading_condition']).strip()})

print("")
print("=== PreOp (no preload) ===")
preop_nopreload_rows = id_map[
    (id_map['State'].astype(str).str.strip().str.lower() == 'preop-nopreload') &
    (id_map['loading_condition'].astype(str).str.strip().str.lower().isin(['flexion', 'extension']))
]
preop_nopreload_elements, preop_nopreload_summary = run_dataset(
    preop_nopreload_rows, PREOP_THRESHOLDS, 'preop_nopreload', 'PreOp (no preload)',
    has_condition=True,
    id_fn=lambda row: {'participant': 'P{}'.format(int(row['participant'])),
                        'loading_condition': str(row['loading_condition']).strip()})

# print("")
# print("=== Oscillation ===")
# osc_rows = id_map[id_map['loading_condition'].astype(str).str.strip().str.lower() == 'oscillation']
# run_dataset(osc_rows, OSC_THRESHOLDS, 'oscillation', 'Oscillation (cumulative)',
#             has_condition=False,
#             id_fn=lambda row: {'participant': 'P{}'.format(int(row['participant']))})


# ============================================================
# Linear mixed-effects tests (R's lme4/lmerTest via rpy2 - Satterthwaite
# t-tests, not the asymptotic z-tests statsmodels gives). Patient is a
# random intercept; loading_condition is averaged out first rather than
# modeled, to keep covariates minimal at N=12. Threshold=0.10 only.
# Run once for PreOp (with preload, Model 1 + Model 2) and once for PreOp
# no preload (Model 1 only).
# ============================================================
LMM_THRESHOLD = 0.10
# Without preload, baseline strain is lower - 0.10 reads as near-zero
# everywhere (no discriminative signal), same reasoning as the lower
# thresholds in plot_oscillation_effect.py.
NOPRELOAD_THRESHOLD = 0.02

os.environ.setdefault('R_HOME', r'C:\Program Files\R\R-4.6.1')
os.environ.setdefault('R_LIBS_USER', os.path.join(os.path.expanduser('~'), 'Documents', 'R', 'win-library', '4.6'))
os.environ['PATH'] = os.path.join(os.environ['R_HOME'], 'bin', 'x64') + os.pathsep + os.environ['PATH']

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
    which would otherwise corrupt the table's column structure. Split out
    from r_table_to_markdown() so callers that need to convert the R object
    to pandas BEFORE it's overwritten (e.g. fitting several models in a
    row that reuse the same R variable name) can do so immediately after
    fitting, then format later."""
    def esc(v):
        s = '{:.4g}'.format(v) if isinstance(v, float) else str(v)
        return s.replace('|', '\\|')

    cols = [esc(c) for c in df.columns]
    lines = ['| ' + ' | '.join(cols) + ' |', '|' + '|'.join(['---'] * len(cols)) + '|']
    for _, row in df.iterrows():
        lines.append('| ' + ' | '.join(esc(v) for v in row) + ' |')
    return '\n'.join(lines)


def r_table_to_markdown(r_expr):
    """Evaluates an R expression returning a data.frame with a 'Term'
    column, and renders it as Markdown (see r_table_to_markdown_from_df)."""
    with localconverter(ro.default_converter + pandas2ri.converter):
        df = ro.conversion.rpy2py(ro.r(r_expr))
    return r_table_to_markdown_from_df(df)


ro.r('''
    get_coef_df <- function(model) {
        df <- as.data.frame(coef(summary(model)))
        df <- cbind(Term = rownames(df), df)
        rownames(df) <- NULL
        df
    }
''')

NAVY = '#003f5c'    # raw data points
TEAL = '#58a4b0'    # box fill


def _fit_one_condition_tissue_model(elements_df, condition, threshold):
    """Fits pct_above ~ tissue + (1 | patient) on ONE condition's data alone
    (Flexion and Extension are different mechanical regimes, tested as
    completely separate models - not pooled with condition as a covariate).
    Returns (df, gm_mean, wm_mean, p_value)."""
    cond_df = elements_df[elements_df['loading_condition'] == condition]
    rows = []
    for (participant, tissue), grp in cond_df.groupby(['participant', 'tissue_type']):
        rows.append({'patient': participant, 'tissue': tissue,
                     'pct_above': pct_above_for_subset(grp, {'t': threshold})['t']})
    df = pd.DataFrame(rows).dropna(subset=['pct_above'])

    print("--- Model 1 data, {} (patient x tissue, N={}) ---".format(condition, df['patient'].nunique()))
    print(df.to_string(index=False))

    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv['m1_data'] = ro.conversion.py2rpy(df)

    print("")
    print("--- Model 1 ({}): pct_above ~ tissue + (1 | patient) ---".format(condition))
    # y_ij = beta_0 + beta_tissueWM * 1[tissue_ij = WM] + u_i + eps_ij
    # $$y_{ij} = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_{ij}=\text{WM}] + u_i + \varepsilon_{ij}$$
    ro.r('''
        m1_data$patient <- factor(m1_data$patient)
        m1_data$tissue  <- factor(m1_data$tissue, levels = c("GM", "WM"))
        model1 <- lmerTest::lmer(pct_above ~ tissue + (1 | patient), data = m1_data)
        print(summary(model1))
    ''')
    ro.r('''
        fe1 <- fixef(model1)
        gm_mean <- as.numeric(fe1['(Intercept)'])
        wm_mean <- as.numeric(fe1['(Intercept)'] + fe1['tissueWM'])
        p1 <- summary(model1)$coefficients['tissueWM', 'Pr(>|t|)']
        coef_df <- get_coef_df(model1)
    ''')
    gm_mean, wm_mean = ro.r('gm_mean')[0], ro.r('wm_mean')[0]
    p1 = ro.r('p1')[0]
    with localconverter(ro.default_converter + pandas2ri.converter):
        coef_df = ro.conversion.rpy2py(ro.r('coef_df'))

    return df, gm_mean, wm_mean, p1, coef_df


def run_model1_lmm(elements_df, tag, title_prefix, threshold=LMM_THRESHOLD, ymax=None):
    """Fits GM vs WM completely separately for Flexion and Extension (two
    independent pct_above ~ tissue + (1 | patient) models, not one pooled
    model with condition as a covariate - the two loading modes differ too
    much in magnitude to assume a shared tissue effect). Saves a 2-panel
    boxplot (Flexion left, Extension right, shared y-axis) and residual QQ
    plots, returns the Model 1 Markdown section covering both conditions."""
    results = {}
    for condition in ('Flexion', 'Extension'):
        df, gm_mean, wm_mean, p_val, coef_df = _fit_one_condition_tissue_model(elements_df, condition, threshold)
        results[condition] = {'df': df, 'gm_mean': gm_mean, 'wm_mean': wm_mean, 'p': p_val, 'coef_df': coef_df}

        qq_path = os.path.join(OUT_DIR, 'lmm_model1_residual_qq_{}_{}.pdf'.format(tag, condition.lower()))
        ro.globalenv['model1_qq_path'] = qq_path
        ro.r('''
            pdf(model1_qq_path, width = 5, height = 5)
            qqnorm(resid(model1), main = "Model 1 ({}) residual Q-Q"); qqline(resid(model1), col = "red")
            dev.off()
            cat("Model 1 ({}) Shapiro-Wilk: "); print(shapiro.test(resid(model1)))
        '''.format(condition, condition))
        print("  Plot saved: {}".format(qq_path))

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
        marker = CONDITION_MARKERS[condition.lower()]
        ax.scatter([0] * len(gm_vals), gm_vals.values, color=NAVY, s=25, alpha=0.7, marker=marker, zorder=3)
        ax.scatter([1] * len(wm_vals), wm_vals.values, color=NAVY, s=25, alpha=0.7, marker=marker, zorder=3)
        ax.scatter([0, 1], [r['gm_mean'], r['wm_mean']], marker='d', s=80, facecolor='red',
                   edgecolor='black', linewidth=1.5, zorder=5)

        p_label = 'p < 0.001' if r['p'] < 0.001 else 'p = {:.3f}'.format(r['p'])
        bracket_y, tick = (58 / 70) * ymax, (1.5 / 70) * ymax
        ax.plot([0, 0, 1, 1], [bracket_y - tick, bracket_y, bracket_y, bracket_y - tick],
                color='black', linewidth=1.2, zorder=6)
        ax.text(0.5, bracket_y + (1 / 70) * ymax, p_label, ha='center', va='bottom', fontsize=10)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['GM', 'WM'])
        ax.set_xlabel(condition)
        ax.set_ylim(0, ymax)

    axes[0].set_ylabel('% of tissue volume above threshold ({:.2f})'.format(threshold))
    lme_handle = [Line2D([0], [0], marker='d', linestyle='', markerfacecolor='red', markeredgecolor='black',
                          label='LME estimate')]
    axes[0].legend(handles=lme_handle, loc='upper left', frameon=False)
    fig.suptitle('{}: Grey Matter vs. White Matter'.format(title_prefix))
    fig.tight_layout()
    box_path = os.path.join(OUT_DIR, 'lmm_model1_boxplot_{}.pdf'.format(tag))
    fig.savefig(box_path, bbox_inches='tight')
    plt.close(fig)
    print("  Plot saved: {}".format(box_path))

    sections = []
    for condition in ('Flexion', 'Extension'):
        r = results[condition]
        sections.append("""\
## Model 1 ({condition}): GM vs WM

$$y_i = \\beta_0 + \\beta_{{\\text{{tissueWM}}}}\\,\\mathbb{{1}}[\\text{{tissue}}_i=\\text{{WM}}] + u_i + \\varepsilon_i$$

- $y_i$: `pct_above` for patient $i$'s tissue observation (GM or WM), {condition} only
- $\\beta_0$: intercept - expected `pct_above` for GM (the reference level)
- $\\beta_{{\\text{{tissueWM}}}}$: fixed effect of WM vs GM
- $u_i \\sim \\mathcal{{N}}(0,\\tau^2)$: random intercept per patient
- $\\varepsilon_i \\sim \\mathcal{{N}}(0,\\sigma^2)$: residual error

`lmer(pct_above ~ tissue + (1 | patient))`, {condition} data only

$$H_0:\\ \\beta_{{\\text{{tissueWM}}}} = 0$$

No difference in % of cord volume above MPS $={threshold:.2f}$ between grey and white \
matter, within {condition} ({pct_lbl}).

{coef_table}
""".format(condition=condition, threshold=threshold,
           pct_lbl='not pooled with Extension' if condition == 'Flexion' else 'not pooled with Flexion',
           coef_table=r_table_to_markdown_from_df(r['coef_df'])))

    return "\n".join(sections)


PREAMBLE = """\
Patient is a random intercept; loading condition (Flexion/Extension) is kept \
as its own main-effect covariate rather than averaged away - they differ \
hugely in magnitude, so averaging would blend two different mechanical \
regimes into one number. Satterthwaite-df t-tests (R `lme4`/`lmerTest`), \
not asymptotic z.
"""

# --- PreOp (with preload): Model 1 + Model 2 ---
print("")
print("=== Linear mixed-effects models: PreOp (with preload), threshold=0.10 ===")
model1_section_preop = run_model1_lmm(preop_elements, 'preop', 'PreOp with Preload', ymax=70)

# Model 2 data: tissue x region, from the summary table already built
# above - fit SEPARATELY per condition (Flexion/Extension are different
# mechanical regimes, not pooled with a covariate - see run_model1_lmm).
model2_sections = []
for condition in ('Flexion', 'Extension'):
    m2_df = (preop_summary[(preop_summary['threshold'] == 't0p10')
                           & (preop_summary['loading_condition'] == condition)]
             .rename(columns={'participant': 'patient'})
             [['patient', 'tissue', 'region', 'pct_above']]
             .dropna(subset=['pct_above']))

    print("")
    print("--- Model 2 data, {} (patient x tissue x region, N={}) ---".format(
        condition, m2_df['patient'].nunique()))
    print(m2_df.to_string(index=False))

    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv['m2_data'] = ro.conversion.py2rpy(m2_df)

    print("")
    print("--- Model 2 ({}): pct_above ~ tissue + region + (1 | patient) ---".format(condition))
    # y_ijk = beta_0 + beta_tissueWM*1[tissue=WM] + beta_regionInterior*1[region=Interior] + u_i + eps_ijk
    # $$y_{ijk} = \beta_0 + \beta_{\text{tissueWM}}\,\mathbb{1}[\text{tissue}_{ijk}=\text{WM}] + \beta_{\text{regionInterior}}\,\mathbb{1}[\text{region}_{ijk}=\text{Interior}] + u_i + \varepsilon_{ijk}$$
    ro.r('''
        m2_data$patient <- factor(m2_data$patient)
        m2_data$tissue  <- factor(m2_data$tissue, levels = c("GM", "WM"))
        m2_data$region  <- factor(m2_data$region, levels = c("Boundary", "Interior"))
        model2 <- lmerTest::lmer(pct_above ~ tissue + region + (1 | patient), data = m2_data)
        print(summary(model2))
        coef_df2 <- get_coef_df(model2)
    ''')
    with localconverter(ro.default_converter + pandas2ri.converter):
        coef_df2 = ro.conversion.rpy2py(ro.r('coef_df2'))

    model2_qq_path = os.path.join(OUT_DIR, 'lmm_model2_residual_qq_preop_{}.pdf'.format(condition.lower()))
    ro.globalenv['model2_qq_path'] = model2_qq_path
    ro.r('''
        pdf(model2_qq_path, width = 5, height = 5)
        qqnorm(resid(model2), main = "Model 2 ({}) residual Q-Q"); qqline(resid(model2), col = "red")
        dev.off()
        cat("Model 2 ({}) Shapiro-Wilk: "); print(shapiro.test(resid(model2)))
    '''.format(condition, condition))
    print("  Plot saved: {}".format(model2_qq_path))

    model2_sections.append("""\
## Model 2 ({condition}): GM vs WM, boundary vs interior

$$y_{{ijk}} = \\beta_0 + \\beta_{{\\text{{tissueWM}}}}\\,\\mathbb{{1}}[\\text{{tissue}}_{{ijk}}=\\text{{WM}}] + \\beta_{{\\text{{regionInterior}}}}\\,\\mathbb{{1}}[\\text{{region}}_{{ijk}}=\\text{{Interior}}] + u_i + \\varepsilon_{{ijk}}$$

- $y_{{ijk}}$: `pct_above` for patient $i$, tissue $j$, region $k$, {condition} only
- $\\beta_0$: intercept - expected `pct_above` for GM x Boundary (the reference levels)
- $\\beta_{{\\text{{tissueWM}}}}$: fixed effect of WM vs GM, controlling for region
- $\\beta_{{\\text{{regionInterior}}}}$: fixed effect of Interior vs Boundary, controlling for tissue
- $u_i \\sim \\mathcal{{N}}(0,\\tau^2)$: random intercept per patient
- $\\varepsilon_{{ijk}} \\sim \\mathcal{{N}}(0,\\sigma^2)$: residual error

`lmer(pct_above ~ tissue + region + (1 | patient))`, {condition} data only

$$H_0:\\ \\beta_{{\\text{{tissueWM}}}} = 0 \\ \\text{{and}} \\ \\beta_{{\\text{{regionInterior}}}} = 0$$

No difference in % of cord volume above MPS $=0.10$ between grey and white \
matter, or between boundary and interior tissue, within {condition}. Both \
factors are 2-level, so an F-test here would just be t^2 - redundant with \
the t-tests below.

{coef_table}
""".format(condition=condition, coef_table=r_table_to_markdown_from_df(coef_df2)))

model2_section = "\n".join(model2_sections)

summary_md_preop = (
    "# GM/WM linear mixed-effects models - PreOp with preload (threshold = 0.10)\n\n"
    + PREAMBLE + "\n" + model1_section_preop + "\n" + model2_section
)
summary_md_preop_path = os.path.join(OUT_DIR, 'lmm_summary_preop.md')
with open(summary_md_preop_path, 'w', encoding='utf-8') as f:
    f.write(summary_md_preop)
print("")
print("Summary saved: {}".format(summary_md_preop_path))

# --- PreOp (no preload): Model 1 only ---
print("")
print("=== Linear mixed-effects models: PreOp (no preload), threshold={:.2f} ===".format(NOPRELOAD_THRESHOLD))
model1_section_nopreload = run_model1_lmm(preop_nopreload_elements, 'preop_nopreload', 'PreOp without Preload',
                                           threshold=NOPRELOAD_THRESHOLD, ymax=70)

summary_md_nopreload = (
    "# GM/WM linear mixed-effects model - PreOp without preload (threshold = {:.2f})\n\n".format(NOPRELOAD_THRESHOLD)
    + PREAMBLE + "\n" + model1_section_nopreload
)
summary_md_nopreload_path = os.path.join(OUT_DIR, 'lmm_summary_preop_nopreload.md')
with open(summary_md_nopreload_path, 'w', encoding='utf-8') as f:
    f.write(summary_md_nopreload)
print("")
print("Summary saved: {}".format(summary_md_nopreload_path))

# TODO: PostOp GM/WM tissue-distinct model goes here, once PostOp's
# '_mps_GM_WM.csv' extraction has run - a pooled
# pct_above ~ tissue * state + (1 | patient) over PreOp-NoPreload + PostOp,
# both at threshold=0.02 (valid now since both states share one threshold,
# unlike the 3-state PreOp-with-preload comparison). Goes above the
# whole-cord analysis below once built.

# ============================================================
# Whole cord (not tissue-split): PreOp (no preload) vs PostOp, threshold=0.02.
# Reuses the existing cache_prepost_sortedbymjoachange_peak.csv (built by
# "Alex_results_multipatient_plot - compare_threshold.py" for
# plot_prepost_sortedbymjoachange.py) rather than the GM/WM-tagged caches
# above - this data has no tissue_type column, just whole-Cord mps/volume.
# ============================================================
WHOLE_CORD_STATE_LABELS = {'preop-nopreload': 'PreOp (no preload)', 'postop': 'PostOp'}


def _fit_one_condition_state_model(df, condition, threshold):
    """Fits pct_above ~ state + (1 | patient) on ONE condition's whole-cord
    data alone (Flexion/Extension fit as completely separate models, not
    pooled with a covariate - see run_model1_lmm's docstring for why).
    Returns (df, pre_mean, post_mean, p_value, coef_df)."""
    cond_df = df[df['loading_condition'] == condition]
    rows = []
    for (participant, state_norm), grp in cond_df.groupby(['participant', 'state_norm']):
        rows.append({'patient': 'P{}'.format(int(participant)), 'state': WHOLE_CORD_STATE_LABELS[state_norm],
                     'pct_above': pct_above_for_subset(grp, {'t': threshold})['t']})
    wc_df = pd.DataFrame(rows).dropna(subset=['pct_above'])

    print("--- Whole-cord data, {} (patient x state, N={}) ---".format(condition, wc_df['patient'].nunique()))
    print(wc_df.to_string(index=False))

    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv['wc_data'] = ro.conversion.py2rpy(wc_df)

    print("")
    print("--- Whole cord ({}): pct_above ~ state + (1 | patient) ---".format(condition))
    # y_i = beta_0 + beta_statePostOp * 1[state_i = PostOp] + u_i + eps_i
    # $$y_i = \beta_0 + \beta_{\text{statePostOp}}\,\mathbb{1}[\text{state}_i=\text{PostOp}] + u_i + \varepsilon_i$$
    ro.r('''
        wc_data$patient <- factor(wc_data$patient)
        wc_data$state   <- factor(wc_data$state, levels = c("PreOp (no preload)", "PostOp"))
        model_wc <- lmerTest::lmer(pct_above ~ state + (1 | patient), data = wc_data)
        print(summary(model_wc))
    ''')
    ro.r('''
        fe_wc <- fixef(model_wc)
        pre_mean  <- as.numeric(fe_wc['(Intercept)'])
        post_mean <- as.numeric(fe_wc['(Intercept)'] + fe_wc['statePostOp'])
        p_wc <- summary(model_wc)$coefficients['statePostOp', 'Pr(>|t|)']
        coef_df_wc <- get_coef_df(model_wc)
    ''')
    pre_mean, post_mean = ro.r('pre_mean')[0], ro.r('post_mean')[0]
    p_wc = ro.r('p_wc')[0]
    with localconverter(ro.default_converter + pandas2ri.converter):
        coef_df = ro.conversion.rpy2py(ro.r('coef_df_wc'))

    return wc_df, pre_mean, post_mean, p_wc, coef_df


def run_wholecord_state_lmm(df, tag, threshold=NOPRELOAD_THRESHOLD):
    """Fits PreOp-NoPreload vs PostOp (whole cord, not tissue-split)
    completely separately for Flexion and Extension (two independent
    pct_above ~ state + (1 | patient) models - see run_model1_lmm's
    docstring for why they aren't pooled with a covariate). Saves a 2-panel
    boxplot (Flexion left, Extension right, shared y-axis) and residual QQ
    plots, returns the Markdown section covering both conditions."""
    df = df.copy()
    df['state_norm'] = df['state'].astype(str).str.strip().str.lower()
    df = df[df['state_norm'].isin(WHOLE_CORD_STATE_LABELS)]

    results = {}
    for condition in ('Flexion', 'Extension'):
        wc_df, pre_mean, post_mean, p_val, coef_df = _fit_one_condition_state_model(df, condition, threshold)
        results[condition] = {'df': wc_df, 'pre_mean': pre_mean, 'post_mean': post_mean,
                               'p': p_val, 'coef_df': coef_df}

        qq_path = os.path.join(OUT_DIR, 'lmm_wholecord_residual_qq_{}_{}.pdf'.format(tag, condition.lower()))
        ro.globalenv['wc_qq_path'] = qq_path
        ro.r('''
            pdf(wc_qq_path, width = 5, height = 5)
            qqnorm(resid(model_wc), main = "Whole-cord ({}) residual Q-Q"); qqline(resid(model_wc), col = "red")
            dev.off()
            cat("Whole-cord ({}) Shapiro-Wilk: "); print(shapiro.test(resid(model_wc)))
        '''.format(condition, condition))
        print("  Plot saved: {}".format(qq_path))

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
        marker = CONDITION_MARKERS[condition.lower()]
        ax.scatter([0] * len(pre_vals), pre_vals.values, color=NAVY, s=25, alpha=0.7, marker=marker, zorder=3)
        ax.scatter([1] * len(post_vals), post_vals.values, color=NAVY, s=25, alpha=0.7, marker=marker, zorder=3)
        ax.scatter([0, 1], [r['pre_mean'], r['post_mean']], marker='d', s=80, facecolor='red',
                   edgecolor='black', linewidth=1.5, zorder=5)

        p_label = 'p < 0.001' if r['p'] < 0.001 else 'p = {:.3f}'.format(r['p'])
        bracket_y, tick = (58 / 70) * ymax, (1.5 / 70) * ymax
        ax.plot([0, 0, 1, 1], [bracket_y - tick, bracket_y, bracket_y, bracket_y - tick],
                color='black', linewidth=1.2, zorder=6)
        ax.text(0.5, bracket_y + (1 / 70) * ymax, p_label, ha='center', va='bottom', fontsize=10)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['PreOp\n(no preload)', 'PostOp'])
        ax.set_xlabel(condition)
        ax.set_ylim(0, ymax)

    axes[0].set_ylabel('% of whole cord volume above threshold ({:.2f})'.format(threshold))
    lme_handle = [Line2D([0], [0], marker='d', linestyle='', markerfacecolor='red', markeredgecolor='black',
                          label='LME estimate')]
    axes[0].legend(handles=lme_handle, loc='upper left', frameon=False)
    fig.suptitle('Whole Cord: PreOp (no preload) vs PostOp')
    fig.tight_layout()
    box_path = os.path.join(OUT_DIR, 'lmm_wholecord_boxplot_{}.pdf'.format(tag))
    fig.savefig(box_path, bbox_inches='tight')
    plt.close(fig)
    print("  Plot saved: {}".format(box_path))

    sections = []
    for condition in ('Flexion', 'Extension'):
        r = results[condition]
        sections.append("""\
## Whole cord ({condition}): PreOp (no preload) vs PostOp

$$y_i = \\beta_0 + \\beta_{{\\text{{statePostOp}}}}\\,\\mathbb{{1}}[\\text{{state}}_i=\\text{{PostOp}}] + u_i + \\varepsilon_i$$

- $y_i$: `pct_above` (whole cord, not split by tissue) for patient $i$, {condition} only
- $\\beta_0$: intercept - expected `pct_above` for PreOp (no preload) (the reference level)
- $\\beta_{{\\text{{statePostOp}}}}$: fixed effect of PostOp vs PreOp (no preload)
- $u_i \\sim \\mathcal{{N}}(0,\\tau^2)$: random intercept per patient
- $\\varepsilon_i \\sim \\mathcal{{N}}(0,\\sigma^2)$: residual error

`lmer(pct_above ~ state + (1 | patient))`, {condition} data only

$$H_0:\\ \\beta_{{\\text{{statePostOp}}}} = 0$$

No difference in % of whole cord volume above MPS $={threshold:.2f}$ between PreOp \
(no preload) and PostOp, within {condition}.

{coef_table}
""".format(condition=condition, threshold=threshold, coef_table=r_table_to_markdown_from_df(r['coef_df'])))

    return "\n".join(sections)


print("")
print("=== Linear mixed-effects model: Whole cord, PreOp (no preload) vs PostOp, threshold={:.2f} ===".format(
    NOPRELOAD_THRESHOLD))
WHOLE_CORD_CACHE_PATH = os.path.join(OUT_DIR, 'cache_prepost_sortedbymjoachange_peak.csv')
whole_cord_df = pd.read_csv(WHOLE_CORD_CACHE_PATH)
wholecord_section = run_wholecord_state_lmm(whole_cord_df, 'prepost')

summary_md_wholecord = (
    "# Whole-cord linear mixed-effects model - PreOp (no preload) vs PostOp (threshold = {:.2f})\n\n".format(
        NOPRELOAD_THRESHOLD)
    + PREAMBLE + "\n" + wholecord_section
)
summary_md_wholecord_path = os.path.join(OUT_DIR, 'lmm_summary_wholecord_prepost.md')
with open(summary_md_wholecord_path, 'w', encoding='utf-8') as f:
    f.write(summary_md_wholecord)
print("")
print("Summary saved: {}".format(summary_md_wholecord_path))
