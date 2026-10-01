"""
plot_preload_effect.py - Compare PreOp WITH the simulated compression-site
preload (id_map.csv State='PreOp') against PreOp WITHOUT it
(State='PreOp-NoPreload'), patient by patient, for the participants that have
both variants (currently P1-P12). Isolates the effect of the preload modeling
assumption itself, independent of surgery - both states are pre-operative, so
patients are ordered by pre-op mJOA (not mJOA change, which needs a post-op
value that doesn't apply here).

Plot A: paired threshold comparison, same visual encoding as the PreOp/PostOp
"Fifth plot" in "Alex_results_multipatient_plot - compare_threshold.py" -
color = threshold, marker shape = loading condition, marker fill = state
(solid = PreOp w/ preload, hollow = PreOp-NoPreload), Flexion/Extension pairs
connected within the same state.

Caches the reduced per-element (mps, volume) data to
cache_preload_effect_{FRAME_MODE}.csv so re-running after changing
thresholds/colors is fast and never re-reads the raw per-frame CSVs. On first
run (no cache yet): reuses PreOp-NoPreload rows from the existing
cache_prepost_sortedbymjoachange_{FRAME_MODE}.csv (built by
"Alex_results_multipatient_plot - compare_threshold.py") if present, instead
of re-reading raw CSVs for that half - that cache has no PreOp-WITH-preload
rows (it was built for the PreOp-NoPreload vs PostOp comparison), so that half
always has to be read from id_map.csv's csv_path column, restricted to just
the participants that have a PreOp-NoPreload row.

Run: python plot_preload_effect.py
"""

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from mps_common import pct_volume_above, PLOT_STYLE

# ============================================================
# USER SETTINGS - edit these and re-run for fast iteration
# ============================================================
ID_MAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'id_map.csv')
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# 'last' - use each element's MPS at the last frame of the step
# 'peak' - use each element's peak (max-ever) MPS across all frames
# Must match the FRAME_MODE used to build cache_prepost_sortedbymjoachange_*.csv
# if that cache is to be reused for the PreOp-NoPreload half.
FRAME_MODE = 'peak'

STATE_PRELOAD = 'PreOp'
STATE_NO_PRELOAD = 'PreOp-NoPreload'
STATE_TITLES = {STATE_PRELOAD: 'PreOp (with preload)', STATE_NO_PRELOAD: 'PreOp (no preload)'}
STATE_FILLED = {STATE_PRELOAD: True, STATE_NO_PRELOAD: False}   # solid vs hollow marker, Plot A

CONDITION_MARKERS = {'flexion': 's', 'extension': '^'}

# Plot A - same 4 manual thresholds used throughout the multipatient scripts.
PLOT_A_THRESHOLDS = {'t0p10': 0.10, 't0p15': 0.15}
PLOT_A_THRESHOLD_COLORS = {'t0p10': '#2e75b6', 't0p15': '#c00000'}
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


# ============================================================
# Cache: reduced (participant, loading_condition, state) -> [element_label, mps, volume]
# ============================================================
own_cache_path = os.path.join(OUT_DIR, 'cache_preload_effect_{}.csv'.format(FRAME_MODE))
shared_cache_path = os.path.join(OUT_DIR, 'cache_prepost_sortedbymjoachange_{}.csv'.format(FRAME_MODE))

id_map = pd.read_csv(ID_MAP_PATH, skipinitialspace=True)

# Participants that have a PreOp-NoPreload row - this comparison is restricted to
# just these (currently P1-P12); reading raw CSVs for PreOp patients outside this
# set would be wasted work since they have nothing to compare against.
target_participants = sorted(set(
    int(row['participant']) for _, row in id_map.iterrows()
    if str(row.get('State', '')).strip() == STATE_NO_PRELOAD
))
if not target_participants:
    raise SystemExit("No '{}' rows found in id_map.csv.".format(STATE_NO_PRELOAD))
print("Participants with both PreOp and PreOp-NoPreload data: {}".format(
    ', '.join('P{}'.format(p) for p in target_participants)))

if os.path.isfile(own_cache_path):
    print("Loading cached reduced data: {}".format(own_cache_path))
    cache_df = pd.read_csv(own_cache_path)
else:
    print("No cache at {} yet - building it.".format(own_cache_path))
    cache_parts = []

    # Reuse PreOp-NoPreload rows from the existing prepost cache if present, instead
    # of re-reading raw per-frame CSVs for data that's already been reduced there.
    reused_nopreload = False
    if os.path.isfile(shared_cache_path):
        shared = pd.read_csv(shared_cache_path)
        shared_nopreload = shared[
            (shared['state'] == STATE_NO_PRELOAD) & (shared['participant'].isin(target_participants))]
        if not shared_nopreload.empty:
            cache_parts.append(shared_nopreload[
                ['participant', 'loading_condition', 'state', 'element_label', 'mps', 'volume']])
            reused_nopreload = True
            print("  Reused {} PreOp-NoPreload row(s) from {}".format(len(shared_nopreload), shared_cache_path))

    # States still needing a raw-CSV read: PreOp-WITH-preload always (never cached
    # anywhere), PreOp-NoPreload only if the shared cache didn't have it.
    states_to_read = {STATE_PRELOAD} if reused_nopreload else {STATE_PRELOAD, STATE_NO_PRELOAD}

    missing = []
    for _, row in id_map.iterrows():
        state = str(row.get('State', '')).strip()
        if state not in states_to_read:
            continue
        participant = int(row['participant'])
        if participant not in target_participants:
            continue
        condition = str(row.get('loading_condition', '')).strip()
        if condition.strip().lower() not in ('flexion', 'extension'):
            continue
        csv_path = str(row.get('csv_path', '')).strip()
        if not csv_path or csv_path.lower() == 'nan' or not os.path.isfile(csv_path):
            missing.append((participant, condition, state))
            continue
        raw = pd.read_csv(csv_path)
        reduced = reduce_to_frame_mode(raw, FRAME_MODE)
        reduced.insert(0, 'state', state)
        reduced.insert(0, 'loading_condition', condition)
        reduced.insert(0, 'participant', participant)
        cache_parts.append(reduced)

    if missing:
        print("Skipping {} row(s) with no csv_path set (or file not found):".format(len(missing)))
        for p, c, s in missing:
            print("  P{} ({}, {})".format(p, c, s))

    if not cache_parts:
        raise SystemExit("No data loaded - check id_map.csv.")

    cache_df = pd.concat(cache_parts, ignore_index=True)
    cache_df.to_csv(own_cache_path, index=False)
    print("Cached reduced data: {}".format(own_cache_path))

state_per_patient = {
    (int(p), c, s): grp[['element_label', 'mps', 'volume']].reset_index(drop=True)
    for (p, c, s), grp in cache_df.groupby(['participant', 'loading_condition', 'state'])
}

# Guard: only keep participants actually present in BOTH states (should already be
# exactly target_participants by construction, but a participant/condition missing
# its csv_path above could narrow this).
participants_with_preload = {p for (p, c, s) in state_per_patient if s == STATE_PRELOAD}
participants_with_nopreload = {p for (p, c, s) in state_per_patient if s == STATE_NO_PRELOAD}
target_participants = sorted(participants_with_preload & participants_with_nopreload)
if not target_participants:
    raise SystemExit("No participants have both PreOp and PreOp-NoPreload reduced data.")

# Pre-op mJOA - from State == 'PreOp' rows only, used to order every plot and label
# each patient (same value applies regardless of which preload variant is plotted).
preop_mjoa_by_participant = {}
for _, row in id_map.iterrows():
    if str(row.get('State', '')).strip() != STATE_PRELOAD:
        continue
    p_label = 'P{}'.format(int(row['participant']))
    preop_mjoa_by_participant[p_label] = row.get('mJOA', '')

participants = sorted(('P{}'.format(p) for p in target_participants),
                       key=lambda p: (preop_mjoa_by_participant.get(p, float('inf')), int(p[1:])))
x_pos = {p: i for i, p in enumerate(participants)}

# ============================================================
# Plot A: paired threshold comparison, PreOp (solid) vs PreOp-NoPreload (hollow)
# ============================================================
a_records = []
for (a_p, a_c, a_s), a_df in sorted(state_per_patient.items()):
    if a_p not in target_participants or a_s not in (STATE_PRELOAD, STATE_NO_PRELOAD):
        continue
    a_record = {'participant': 'P{}'.format(a_p), 'loading_condition': a_c, 'state': a_s}
    for a_name, a_val in PLOT_A_THRESHOLDS.items():
        a_record['pct_above_{}'.format(a_name)] = pct_volume_above(a_df, a_val)
    a_records.append(a_record)
a_summary = pd.DataFrame(a_records)

a_summary_path = os.path.join(OUT_DIR, 'multipatient_mps_summary_preload_effect_sortedbypreopmJOA.csv')
a_summary.to_csv(a_summary_path, index=False)
print(a_summary.to_string(index=False))

a_fig, a_ax = plt.subplots(figsize=(8, 5.5))

for a_name in PLOT_A_THRESHOLDS:
    a_col = 'pct_above_{}'.format(a_name)
    a_color = PLOT_A_THRESHOLD_COLORS[a_name]
    # Connect Flexion<->Extension pairs within the SAME state only.
    for (_, _), grp in a_summary.groupby(['participant', 'state']):
        if len(grp) == 2:
            xp = x_pos[grp['participant'].iloc[0]]
            a_ax.vlines(xp, grp[a_col].min(), grp[a_col].max(), color=a_color, linewidth=1.0, alpha=0.5, zorder=2)
    for (a_cond, a_state), grp in a_summary.groupby(['loading_condition', 'state']):
        marker = CONDITION_MARKERS.get(a_cond.strip().lower(), 'o')
        filled = STATE_FILLED.get(a_state, True)
        xs = [x_pos[p] for p in grp['participant']]
        if filled:
            a_ax.scatter(xs, grp[a_col], color=a_color, marker=marker, s=55, zorder=3)
        else:
            a_ax.scatter(xs, grp[a_col], facecolors='none', edgecolors=a_color, marker=marker,
                         s=55, linewidths=1.3, zorder=3)

a_threshold_handles = [Line2D([0], [0], marker='o', linestyle='', color=PLOT_A_THRESHOLD_COLORS[name],
                               label='{:.2f}'.format(val))
                        for name, val in PLOT_A_THRESHOLDS.items()]
a_condition_handles = [Line2D([0], [0], marker=marker, linestyle='', color='black', label=cond.capitalize())
                        for cond, marker in CONDITION_MARKERS.items()]
a_state_handles = [
    Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='black', label='PreOp'),
    Line2D([0], [0], marker='o', linestyle='', color='black', markerfacecolor='none', label='PreOp-NoPreload'),
]
a_blank = Line2D([0], [0], linestyle='none', marker='None', label='')

a_all_handles = (
    [Line2D([0], [0], linestyle='none', marker='None', label='Threshold')] + a_threshold_handles +
    [a_blank] +
    [Line2D([0], [0], linestyle='none', marker='None', label='Loading condition')] + a_condition_handles +
    [a_blank] +
    [Line2D([0], [0], linestyle='none', marker='None', label='State')] + a_state_handles
)
a_ax.legend(handles=a_all_handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

a_ax.set_xticks(range(len(participants)))
a_ax.set_xticklabels(['{}\n({})'.format(p, preop_mjoa_by_participant.get(p, '?')) for p in participants])
a_ax.set_xlabel('Participant (ordered by pre-op mJOA, ascending)')
a_ax.set_ylabel('% cord volume above threshold')
a_ax.set_title('Effect of simulated compression-site preload (PreOp)')
a_fig.tight_layout()

a_plot_path = os.path.join(OUT_DIR, 'multipatient_mps_plot_preload_effect_sortedbypreopmJOA.pdf')
a_fig.savefig(a_plot_path, bbox_inches='tight')
plt.close(a_fig)

print("Summary saved: {}".format(a_summary_path))
print("Plot saved: {}".format(a_plot_path))

