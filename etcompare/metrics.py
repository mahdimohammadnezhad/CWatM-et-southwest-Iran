"""Agreement statistics between simulated (s) and reference (o) series."""
import numpy as np
import pandas as pd


def pbias(s, o):
    """Percent bias; positive values mean the simulation is higher."""
    return 100.0 * (np.sum(s) - np.sum(o)) / np.sum(o)


def scores(s, o):
    """N, Pearson r, RMSE, PBIAS, NSE and KGE (Gupta et al., 2009)."""
    s, o = np.asarray(s, float), np.asarray(o, float)
    r = np.corrcoef(s, o)[0, 1]
    alpha, beta = s.std() / o.std(), s.mean() / o.mean()
    return {
        "N": len(s),
        "r": r,
        "RMSE": np.sqrt(np.mean((s - o) ** 2)),
        "PBIAS": pbias(s, o),
        "NSE": 1 - np.sum((s - o) ** 2) / np.sum((o - o.mean()) ** 2),
        "KGE": 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2),
        "alpha": alpha,
        "beta": beta,
    }


def by_class(s, o, months):
    """PBIAS and mean difference by reference-quartile class (Q1-Q4) and by calendar quarter."""
    s, o = np.asarray(s, float), np.asarray(o, float)
    cls = np.digitize(o, np.percentile(o, [25, 50, 75]), right=True)
    quarter = (np.asarray(months) - 1) // 3
    rows = []
    for kind, groups, labels in [("class", cls, ["Q1", "Q2", "Q3", "Q4"]), ("quarter", quarter, ["JFM", "AMJ", "JAS", "OND"])]:
        for k, lab in enumerate(labels):
            m = groups == k
            rows.append({"kind": kind, "group": lab, "N": int(m.sum()), "PBIAS": pbias(s[m], o[m]),
                         "ref_mean": o[m].mean(), "sim_mean": s[m].mean(), "diff": (s[m] - o[m]).mean()})
    return pd.DataFrame(rows)


def upper_decile_pbias(s, o):
    s, o = np.asarray(s, float), np.asarray(o, float)
    m = o > np.percentile(o, 90)
    return pbias(s[m], o[m]), int(m.sum())
