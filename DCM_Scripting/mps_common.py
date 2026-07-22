"""
mps_common.py - Shared helpers for volume-weighted MPS percentile analysis
(Python 3). Used by Alex_time_history_diagnostic.py and Alex_results_plotting.py.
"""


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


def pct_volume_above(df, threshold, mps_col='mps', vol_col='volume'):
    """% of total volume in df with mps_col >= threshold."""
    total = df[vol_col].sum()
    if total <= 0:
        return 0.0
    above = df.loc[df[mps_col] >= threshold, vol_col].sum()
    return 100.0 * above / total


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
