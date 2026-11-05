"""Figures 2-9 and supplementary Figure S1.

Usage: python scripts/08_figures.py [fig2 fig3 ... figS1]   (all figures by default)
Fig. 9 needs results/annual_climate.csv from script 07.
"""
import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, LightSource, LinearSegmentedColormap, ListedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from etcompare import cwatm, grid, products, projections  # noqa: E402
from etcompare.config import BASELINE, FIGURES, GCMS, MONTHS, PERIODS, RESULTS, SSPS  # noqa: E402
from etcompare.plotting import (GRID, INK, MUTED, STYLE, WIDTH, basin_outlines, draw_outlines, map_axes,  # noqa: E402
                                map_mask, panel_label, setup)

setup()
FIGURES.mkdir(exist_ok=True)
NAMES = products.PRODUCTS
label = products.label


def fig2():
    env, tr = grid.load("00025")
    h, w = env["basin"].shape
    extent = [tr.c, tr.c + tr.a * w, tr.f + tr.e * h, tr.f]
    whole = env["basin"] == 1
    shown = whole & map_mask(tr, (h, w))
    outlines = basin_outlines()
    e1, tr1 = grid.load("001")
    a1 = grid.area_weights(tr1, e1["basin"].shape)
    b1 = e1["basin"] == 1
    share = {c: float(np.nansum((e1[f"lulc_{c}"] * a1)[b1]) / a1[b1].sum()) for c in (1, 2, 4, 5, 7, 8, 11)}
    fig, axs = plt.subplots(1, 3, figsize=(WIDTH, 2.75), gridspec_kw=dict(wspace=0.08))

    ax = axs[0]
    shade = LightSource(azdeg=315, altdeg=45).hillshade(np.nan_to_num(env["elevation"], nan=0), vert_exag=0.02,
                                                       dx=0.0025 * 92000, dy=0.0025 * 111000)
    hyp = LinearSegmentedColormap.from_list("hyp", ["#4f8a5b", "#a8c686", "#e8dfa6", "#c9a36b", "#8c6a4a", "#f2efe9"])
    ax.imshow(np.where(shown, shade, np.nan), extent=extent, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    im = ax.imshow(np.where(shown, env["elevation"], np.nan), extent=extent, cmap=hyp, vmin=0, vmax=4000, alpha=0.78,
                   interpolation="nearest")
    draw_outlines(ax, outlines)
    map_axes(ax)
    cb = fig.colorbar(im, cax=ax.inset_axes([0.06, -0.2, 0.78, 0.045]), orientation="horizontal")
    cb.set_label("Elevation (m a.s.l.), SRTM 30 m", fontsize=7)
    cb.ax.tick_params(labelsize=7)
    for name, (x, y) in {"Karun": (50.3, 32.0), "Karkheh": (48.0, 33.9), "Marzi-e-Gharb": (46.55, 33.5)}.items():
        ax.text(x, y, name, fontsize=6.5, color=INK, ha="center", style="italic",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
    panel_label(ax, "(a)")

    lc = {1: ("Water", "#419bdf"), 2: ("Trees", "#397d49"), 4: ("Flooded vegetation", "#7a87c6"), 5: ("Crops", "#e49635"),
          7: ("Built area", "#c4281b"), 8: ("Bare ground", "#a59b8f"), 11: ("Rangeland", "#e3e2c3")}
    codes = sorted(lc)
    idx = np.full(env["lulc"].shape, np.nan)
    for i, c in enumerate(codes):
        idx[(env["lulc"] == c) & shown] = i
    ax = axs[1]
    ax.imshow(idx, extent=extent, cmap=ListedColormap([lc[c][1] for c in codes]), vmin=-0.5, vmax=len(codes) - 0.5,
              interpolation="nearest")
    draw_outlines(ax, outlines)
    map_axes(ax)
    ax.set_yticklabels([])
    handles = [Patch(facecolor=lc[c][1], edgecolor="#8c8c8c", linewidth=0.3, label=f"{lc[c][0]} ({share[c] * 100:.1f}%)")
               for c in [5, 11, 8, 7, 1, 2, 4]]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=6.3, handlelength=1,
              columnspacing=0.8)
    panel_label(ax, "(b)")

    kc = {22: ("BWh", "#ff0000"), 27: ("BSh", "#f5a500"), 26: ("BSk", "#ffdc64"), 34: ("Csa", "#ffff00"),
          45: ("Dsa", "#ff00ff"), 46: ("Dsb", "#c800c8")}
    kidx = np.full(env["koppen"].shape, np.nan)
    for i, c in enumerate(kc):
        kidx[(env["koppen"] == c) & shown] = i
    ax = axs[2]
    ax.imshow(kidx, extent=extent, cmap=ListedColormap([v[1] for v in kc.values()]), vmin=-0.5, vmax=len(kc) - 0.5,
              interpolation="nearest")
    draw_outlines(ax, outlines)
    map_axes(ax)
    ax.set_yticklabels([])
    k = env["koppen"][whole]
    handles = [Patch(facecolor=v[1], edgecolor="#8c8c8c", linewidth=0.3, label=f"{v[0]} ({(k == c).mean() * 100:.1f}%)")
               for c, v in kc.items()]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=6.5, handlelength=1,
              columnspacing=0.8)
    panel_label(ax, "(c)")
    fig.savefig(FIGURES / "fig2_environment.png")
    plt.close(fig)


def fig3():
    clim = pd.read_csv(RESULTS / "monthly_climatology.csv")
    t1 = pd.read_csv(RESULTS / "table1_agreement.csv", index_col=0)
    fig, axs = plt.subplots(2, 3, figsize=(WIDTH, 4.2), sharex=True, sharey=True, gridspec_kw=dict(wspace=0.08, hspace=0.32))
    months = np.arange(1, 13)
    for k, (ax, name) in enumerate(zip(axs.ravel(), NAMES)):
        c = clim[clim["dataset"] == name].sort_values("month")
        st = STYLE[name]
        for col, sd, color, lab, marker, ls, z in [("cwatm_mean", "cwatm_sd", "#000000", "CWatM", "o", "-", 3),
                                                   ("ref_mean", "ref_sd", st["color"], label(name), st["marker"], st["ls"], 2)]:
            y, s = c[col].to_numpy(), c[sd].to_numpy()
            ax.fill_between(months, y - s, y + s, color=color, alpha=0.12, linewidth=0, zorder=z - 1)
            ax.plot(months, y, color=color, lw=1.6, marker=marker, ms=3.2, ls=ls, label=lab, zorder=z)
        r = t1.loc[name]
        ax.text(0.97, 0.97, f"r = {r['r']:.2f}   KGE = {r['KGE']:.2f}\nPBIAS = {r['PBIAS']:+.1f}%   N = {int(r['N'])}",
                transform=ax.transAxes, va="top", ha="right", fontsize=6.8, color=INK, linespacing=1.4)
        ax.set_title(f"{label(name)} ({r['start'][:4]}–{r['end'][:4]})", loc="left", fontsize=8, color=INK)
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_xticks(months)
        ax.set_xticklabels(list("JFMAMJJASOND"))
        ax.set_ylim(0, 2.6)
        ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.78), fontsize=6.6, handlelength=1.6)
        panel_label(ax, f"({'abcdef'[k]})")
        if k % 3 == 0:
            ax.set_ylabel("ET (mm d$^{-1}$)")
    fig.savefig(FIGURES / "fig3_seasonal_cycle.png")
    plt.close(fig)


def fig4():
    pb = pd.read_csv(RESULTS / "pbias_by_class.csv")
    t1 = pd.read_csv(RESULTS / "table1_agreement.csv", index_col=0)
    cl, qu = pb[pb["kind"] == "class"], pb[pb["kind"] == "quarter"]
    A = cl.pivot(index="dataset", columns="group", values="PBIAS").loc[NAMES, ["Q1", "Q2", "Q3", "Q4"]]
    A["Top 10%"] = t1.loc[NAMES, "PBIAS_upper_decile"]
    Ad = cl.pivot(index="dataset", columns="group", values="diff").loc[NAMES, ["Q1", "Q2", "Q3", "Q4"]]
    B = qu.pivot(index="dataset", columns="group", values="PBIAS").loc[NAMES, ["JFM", "AMJ", "JAS", "OND"]]
    Bd = qu.pivot(index="dataset", columns="group", values="diff").loc[NAMES, ["JFM", "AMJ", "JAS", "OND"]]
    cmap = LinearSegmentedColormap.from_list("div", ["#2166ac", "#92c5de", "#f2f2f2", "#f4a582", "#b2182b"])
    norm = TwoSlopeNorm(vmin=-100, vcenter=0, vmax=100)
    fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 2.9), gridspec_kw=dict(width_ratios=[5, 4], wspace=0.05))
    panels = [(axs[0], A, Ad, "By reference-ET magnitude class",
               ["Q1\n(≤p25)", "Q2\n(p25–p50)", "Q3\n(p50–p75)", "Q4\n(>p75)", "Top 10%\n(>p90)"]),
              (axs[1], B, Bd, "By calendar quarter", ["JFM", "AMJ", "JAS", "OND"])]
    for ax, M, Dm, title, xl in panels:
        ax.imshow(np.clip(M.to_numpy(float), -100, 100), cmap=cmap, norm=norm, aspect="auto")
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                v = M.iat[i, j]
                txt = f"{v:+.0f}%" + (f"\n{Dm.iat[i, j]:+.2f}" if j < Dm.shape[1] else "")
                ax.text(j, i, txt, ha="center", va="center", fontsize=6.6, color="white" if abs(v) > 70 else INK,
                        linespacing=1.15)
        ax.set_xticks(range(M.shape[1]))
        ax.set_xticklabels(xl, fontsize=7)
        ax.set_yticks(range(M.shape[0]))
        ax.set_yticklabels([label(n) for n in M.index] if ax is axs[0] else [])
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xticks(np.arange(-0.5, M.shape[1]), minor=True)
        ax.set_yticks(np.arange(-0.5, M.shape[0]), minor=True)
        ax.grid(which="minor", color="white", lw=1.5)
        ax.tick_params(which="minor", length=0)
        ax.set_title(title, loc="left", fontsize=8, color=INK)
    panel_label(axs[0], "(a)", x=-0.25)
    panel_label(axs[1], "(b)", x=-0.02)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs, orientation="horizontal", fraction=0.05,
                      pad=0.14, aspect=50, extend="both")
    cb.set_label("PBIAS of CWatM (%); colour saturates at ±100%. Second line: mean difference CWatM − reference (mm d$^{-1}$)",
                 fontsize=7)
    fig.savefig(FIGURES / "fig4_pbias.png")
    plt.close(fig)


def fig5():
    env, tr = grid.load("001")
    shape = env["basin"].shape
    cw, _, _ = grid.cwatm_to_grid(grid.cwatm_mean_map(cwatm.load_cells_historical()), tr, shape)
    inside = (env["basin"] == 1) & np.isfinite(cw)
    data = {"CWatM": cw, **products.load_maps(tr, shape)}
    data = {k: np.where(inside, v, np.nan) for k, v in data.items()}
    shown = map_mask(tr, shape)
    extent = [tr.c, tr.c + tr.a * shape[1], tr.f + tr.e * shape[0], tr.f]
    area = grid.area_weights(tr, shape)
    outlines = basin_outlines()
    names = list(data)
    fig, axs = plt.subplots(2, 4, figsize=(WIDTH, 2.35 * 2 + 0.5), squeeze=False, gridspec_kw=dict(wspace=0.06, hspace=0.18))
    norm = BoundaryNorm([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0], 256)
    for k, ax in enumerate(axs.ravel()):
        if k >= len(names):
            ax.axis("off")
            continue
        n = names[k]
        im = ax.imshow(np.where(shown, data[n], np.nan), extent=extent, cmap="YlGnBu", norm=norm, interpolation="nearest")
        draw_outlines(ax, outlines, lw=0.45)
        map_axes(ax)
        if k % 4:
            ax.set_yticklabels([])
        if k < len(names) - 4:
            ax.set_xticklabels([])
        m = grid.weighted_mean(data[n], area)
        ax.set_title(f"{label(n)}\nbasin mean {m:.2f} mm d$^{{-1}}$", fontsize=7.2, color=INK, loc="left", linespacing=1.2)
        ax.text(0.97, 0.97, f"({'abcdefgh'[k]})", transform=ax.transAxes, ha="right", va="top", fontsize=9,
                fontweight="bold", color=INK)
    cb = fig.colorbar(im, ax=axs.ravel().tolist(), orientation="horizontal", fraction=0.04, pad=0.1, aspect=45)
    cb.set_label("Mean actual evapotranspiration (mm d$^{-1}$)")
    fig.savefig(FIGURES / "fig5_mean_et_maps.png")
    plt.close(fig)


def fig6():
    g = json.load(open(RESULTS / "gradients.json"))
    cells = pd.read_csv(RESULTS / "cwatm_cells_environment.csv")
    names = [n for n in ["CWatM"] + NAMES if n in g["elevation"]]
    fig = plt.figure(figsize=(WIDTH, 5.3))
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.22)
    sub = gs[0, 0].subgridspec(2, 1, height_ratios=[4, 1], hspace=0.08)
    ax, axb = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
    bands = list(g["share"]["elevation"])
    x = np.arange(len(bands))
    for n in names:
        st = STYLE[n]
        ax.plot(x, [g["elevation"][n][b] for b in bands], color=st["color"], marker=st["marker"], ls=st["ls"], lw=1.4,
                ms=3.5, label=label(n))
    ax.set_xticks(x)
    ax.set_xticklabels([])
    ax.set_ylabel("Mean ET (mm d$^{-1}$)")
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_title("Elevation band (m a.s.l.)", loc="left", fontsize=8)
    axb.bar(x, [g["share"]["elevation"][b] * 100 for b in bands], color="#bdbdbd", width=0.7)
    axb.set_xlim(ax.get_xlim())
    axb.set_xticks(x)
    axb.set_xticklabels(bands, rotation=35, ha="right", fontsize=6.8)
    axb.set_ylabel("Area (%)", fontsize=7)
    axb.tick_params(axis="y", labelsize=6.5)
    panel_label(ax, "(a)")
    for pos, key, title, letter in [(gs[0, 1], "koppen", "Köppen–Geiger class", "(b)"),
                                    (gs[1, 0], "lulc", "Land-cover class (Esri 10 m, 2017)", "(c)")]:
        ax = fig.add_subplot(pos)
        cls = [c for c, s in g["share"][key].items() if s >= 0.005 and np.isfinite(g[key]["CWatM"].get(c, np.nan))]
        xx = np.arange(len(cls))
        for o, n in zip(np.linspace(-0.25, 0.25, len(names)), names):
            st = STYLE[n]
            ax.plot(xx + o, [g[key][n][c] for c in cls], ls="none", marker=st["marker"], color=st["color"], ms=4.5)
        ax.set_xticks(xx)
        ax.set_xticklabels([f"{c.replace(' ', chr(10))}\n({g['share'][key][c] * 100:.0f}%)" for c in cls], fontsize=6.8)
        ax.set_ylabel("Mean ET (mm d$^{-1}$)")
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_title(title, loc="left", fontsize=8)
        ax.set_xlim(-0.6, len(cls) - 0.4)
        panel_label(ax, letter)
    ax = fig.add_subplot(gs[1, 1])
    kcol = {"BWh": "#d7301f", "BSh": "#fc8d59", "Csa": "#7a9a01", "Dsa": "#7b3294"}
    for k, col in kcol.items():
        c = cells[cells["koppen"] == k]
        ax.scatter(c["elevation"], c["ET"], s=12 + 260 * c["f_crops"], color=col, alpha=0.8, edgecolor="white", linewidth=0.5)
    rr = np.corrcoef(cells["elevation"], cells["ET"])[0, 1]
    ax.text(0.03, 0.97, f"CWatM cells (n = {len(cells)})\nr(ET, elevation) = {rr:.2f}", transform=ax.transAxes, va="top",
            fontsize=6.8)
    ax.set_xlabel("Mean cell elevation (m a.s.l.)")
    ax.set_ylabel("CWatM mean ET (mm d$^{-1}$)")
    ax.grid(color=GRID, lw=0.5)
    h1 = [Line2D([], [], ls="none", marker="o", color=c, ms=5, label=k) for k, c in kcol.items()]
    h2 = [Line2D([], [], ls="none", marker="o", color="#9e9e9e", ms=np.sqrt(12 + 260 * f) * 0.95, label=f"{int(f * 100)}% crops")
          for f in (0.05, 0.3, 0.6)]
    ax.legend(handles=h1 + h2, loc="lower right", fontsize=6.3, ncol=2, handletextpad=0.2, columnspacing=0.6, labelspacing=0.6)
    panel_label(ax, "(d)")
    fig.legend(*fig.axes[0].get_legend_handles_labels(), loc="upper center", ncol=len(names), bbox_to_anchor=(0.5, 0.995),
               fontsize=7.2)
    fig.savefig(FIGURES / "fig6_gradients.png")
    plt.close(fig)


def fig7():
    a = pd.read_csv(RESULTS / "annual_series.csv", index_col=0)
    gcm_col = {"GFDL": "#56B4E9", "IPSL": "#E69F00", "MPI": "#CC79A7", "MRI": "#9467bd", "UKESM": "#8c564b"}
    fig, axs = plt.subplots(3, 1, figsize=(WIDTH, 7.4), sharex=True, gridspec_kw=dict(hspace=0.16))
    for k, (ax, s) in enumerate(zip(axs, SSPS)):
        for x0 in (2025, 2050, 2075):
            ax.axvline(x0, color="#9a9a9a", ls=(0, (4, 3)), lw=0.7)
        for lab, (x0, x1) in {"Historical": (2000, 2019), **{p.replace("-", "–"): v for p, v in PERIODS.items()}}.items():
            ax.text((x0 + x1 + 1) / 2, 1.47, lab, ha="center", va="top", fontsize=7.2, color=MUTED)
        for g, col in gcm_col.items():
            v = a[f"{g}_{s}"].dropna()
            ax.plot(v.index, v, color=col, lw=0.8, alpha=0.75, label=f"{g} (scaled)")
        e = a[f"ensemble_{s}"].dropna()
        ax.plot(e.index, e, color="#d62728", lw=2.0, label="Ensemble mean")
        h = a["CWatM"].loc[2000:2019]
        ax.plot(h.index, h, color="#000000", lw=1.6, marker="o", ms=2.6, label="CWatM (retrospective)")
        for p, col in [("GLDAS", "#0072B2"), ("SMAP L4", "#009E73")]:
            v, va = a[p].dropna(), a[f"{p} adjusted"].dropna()
            ax.plot(v.index, v, color=col, lw=1.2, marker="s", ms=2.4, label=f"{label(p)} (original)")
            ax.plot(va.index, va, color=col, lw=1.2, ls=(0, (3, 2)), marker="s", ms=2.4, mfc="white",
                    label=f"{label(p)} (mean-adjusted)")
        ax.set_ylim(0.35, 1.5)
        ax.set_xlim(2000, 2100)
        ax.grid(axis="y", color=GRID, lw=0.5)
        ax.set_ylabel("Annual mean ET (mm d$^{-1}$)")
        ax.set_title(SSPS[s], loc="left", fontsize=8.5, fontweight="bold")
        panel_label(ax, f"({'abc'[k]})")
    h, l = axs[0].get_legend_handles_labels()
    order = ["CWatM (retrospective)", "GLDAS-Noah (original)", "GLDAS-Noah (mean-adjusted)", "SMAP L4 (original)",
             "SMAP L4 (mean-adjusted)", *[f"{g} (scaled)" for g in GCMS], "Ensemble mean"]
    fig.legend([h[l.index(x)] for x in order], order, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.005), fontsize=7)
    axs[-1].set_xlabel("Year")
    fig.subplots_adjust(bottom=0.09)
    fig.savefig(FIGURES / "fig7_annual_et.png")
    plt.close(fig)


def fig8():
    sig = pd.read_csv(RESULTS / "monthly_signals.csv")
    fig = plt.figure(figsize=(WIDTH, 5.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 2.55], hspace=0.28)
    late = list(PERIODS)[-1]
    ax = fig.add_subplot(gs[0])
    x = np.arange(12)
    first = sig[(sig["period"] == late) & (sig["ssp"] == "SSP1-2.6")].set_index("unit").loc[MONTHS]
    ax.bar(x, first["baseline"], width=0.62, color="#d9d9d9", edgecolor="none",
           label=f"Baseline {BASELINE[0]}–{BASELINE[1]}", zorder=1)
    marks = {"SSP1-2.6": ("o", "#0072B2"), "SSP3-7.0": ("s", "#E69F00"), "SSP5-8.5": ("^", "#D55E00")}
    for k, (s, (mk, col)) in enumerate(marks.items()):
        d = sig[(sig["period"] == late) & (sig["ssp"] == s)].set_index("unit").loc[MONTHS]
        ax.plot(x + (k - 1) * 0.16, d["baseline"] * (1 + d["mean"] / 100), ls="none", marker=mk, ms=4.2, mfc=col,
                mec="white", mew=0.5, label=f"{late.replace('-', '–')} {s}", zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(MONTHS)
    ax.set_xlim(-0.5, 13.1)
    ax.set_ylim(0, 2.6)
    ax.set_ylabel("ET (mm d$^{-1}$)")
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", length=0)
    ax.legend(loc="upper right", ncol=2, fontsize=7, handletextpad=0.3, columnspacing=1.0)
    ax.set_title("Monthly mean ET: baseline and late-century ensemble means", loc="left", fontsize=8, color=INK)
    ax.text(-0.075, 1.02, "(a)", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right", color=INK)

    ax = fig.add_subplot(gs[1])
    cols = MONTHS + ["Annual"]
    rows = [(p, s) for p in PERIODS for s in SSPS.values()]
    cmap = LinearSegmentedColormap.from_list("brbg", ["#8c510a", "#d8b365", "#f2f2f2", "#5ab4ac", "#01665e"])
    norm = TwoSlopeNorm(vmin=-0.3, vcenter=0, vmax=0.3)
    xpos = np.r_[np.arange(12), 12.6]
    for i, (p, s) in enumerate(rows):
        d = sig[(sig["period"] == p) & (sig["ssp"] == s)].set_index("unit")
        for j, c in enumerate(cols):
            v, pct, agree = d.loc[c, "change_mm_d"], d.loc[c, "mean"], d.loc[c, "n_agree"]
            ax.add_patch(plt.Rectangle((xpos[j] - 0.5, i - 0.5), 1, 1, facecolor=cmap(norm(np.clip(v, -0.3, 0.3))),
                                       edgecolor="white", lw=1.2))
            txt = f"{pct:+.1f}" if (j == 12 or abs(pct) < 0.95) else f"{pct:+.0f}"
            ax.text(xpos[j], i + 0.04, txt.replace("-", "−"), ha="center", va="center", fontsize=6.6,
                    color="white" if abs(v) > 0.2 else INK)
            if agree >= 4:
                ax.plot(xpos[j] + 0.33, i - 0.30, marker="o", ms=1.9, color=INK, mec="none")
    ax.set_xlim(-0.5, 13.1)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xticks(xpos)
    ax.set_xticklabels(cols)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([s for _, s in rows], fontsize=7)
    for k, p in enumerate(PERIODS):
        ax.text(-2.55, 3 * k + 1, p.replace("-", "–"), ha="center", va="center", rotation=90, fontsize=7.5, color=INK)
        if k:
            ax.axhline(3 * k - 0.5, color="white", lw=4)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(f"Ensemble-mean change relative to {BASELINE[0]}–{BASELINE[1]} (colour: mm d$^{{-1}}$; numbers: %)",
                 loc="left", fontsize=8, color=INK)
    ax.text(-0.075, 1.02, "(b)", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right", color=INK)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation="horizontal", fraction=0.045,
                      pad=0.07, aspect=45, extend="both")
    cb.set_label("Change in ET (mm d$^{-1}$); dot: at least four of the five GCMs agree on the sign of the change", fontsize=7)
    cb.ax.tick_params(labelsize=7)
    fig.savefig(FIGURES / "fig8_monthly_change.png")
    plt.close(fig)


def figS1():
    runs = cwatm.load_runs()
    hist = cwatm.run_series(runs, "historical")
    base = projections.window(hist, *BASELINE).to_numpy()
    fig, axs = plt.subplots(3, 3, figsize=(WIDTH, 5.8), sharey=True, gridspec_kw=dict(hspace=0.45, wspace=0.08))
    for i, (p, (y0, y1)) in enumerate(PERIODS.items()):
        for j, (s, sname) in enumerate(SSPS.items()):
            ax = axs[i, j]
            data = [base]
            for g in GCMS:
                x, _ = projections.scale_to_overlap(hist, cwatm.run_series(runs, f"{g}_{s}"))
                data.append(projections.window(x, y0, y1).to_numpy())
            bp = ax.boxplot(data, widths=0.6, patch_artist=True, showmeans=True, meanline=True,
                            medianprops=dict(color="#e66101", lw=1.4), meanprops=dict(color=INK, ls="--", lw=0.9),
                            flierprops=dict(marker="o", ms=2.2, mfc="none", mec="#7f7f7f", mew=0.5),
                            whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7), boxprops=dict(lw=0.7))
            for k, b in enumerate(bp["boxes"]):
                b.set_facecolor("#bdbdbd" if k == 0 else "#6baed6")
            ax.set_xticks(range(1, 7))
            ax.set_xticklabels(["Baseline"] + GCMS, fontsize=6.3, rotation=30)
            ax.grid(axis="y", color=GRID, lw=0.5)
            ax.set_title(f"({'abcdefghi'[3 * i + j]}) {sname} | {p.replace('-', '–')}", loc="left", fontsize=7.5)
            if j == 0:
                ax.set_ylabel("Monthly ET (mm d$^{-1}$)")
    fig.savefig(FIGURES / "figS1_monthly_distributions.png")
    plt.close(fig)


def nice_limits(values, n=7):
    """Axis limits on a round step (1, 2, 2.5 or 5 times a power of ten) enclosing the values."""
    lo, hi = np.min(values), np.max(values)
    raw = (hi - lo) / n
    exp = np.floor(np.log10(raw))
    step = next(f for f in (1, 2, 2.5, 5, 10) if raw / 10 ** exp <= f) * 10 ** exp
    return np.floor(lo / step) * step, np.ceil(hi / step) * step, step


def clear_boxes(ax, x, y, boxes, step, pad=0.02):
    """Widen the y-axis until no point lies under the given boxes (axes fractions); limits end on major ticks."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xf = (x - x0) / (x1 - x0)
    for _ in range(30):
        yf = (y - y0) / (y1 - y0)
        hit = [b for b in boxes if np.any((xf > b.x0 - pad) & (xf < b.x1 + pad) & (yf > b.y0 - pad) & (yf < b.y1 + pad))]
        if not hit:
            ticks = ax.yaxis.get_major_locator().tick_values(y0, y1)
            tick = ticks[1] - ticks[0]
            snapped = np.floor(y0 / tick + 1e-9) * tick, np.ceil(y1 / tick - 1e-9) * tick
            if np.allclose(snapped, (y0, y1)):
                break
            y0, y1 = snapped
            continue
        for b in hit:
            if b.y0 > 0.5:
                y1 += step
            else:
                y0 -= step
    ax.set_ylim(y0, y1)


def fig9():
    f = RESULTS / "annual_climate.csv"
    if not f.exists():
        print("fig9: results/annual_climate.csv not found (run script 07)")
        return
    a = pd.read_csv(f)
    drivers = [("P_mm", "Precipitation", "Precipitation (mm)"), ("T_C", "Temperature", "Temperature (°C)")]
    fig, axs = plt.subplots(2, 3, figsize=(WIDTH, 5.4), gridspec_kw=dict(hspace=0.45, wspace=0.42))
    panels = []
    for j, (s, sname) in enumerate(SSPS.items()):
        ens = a[a["run"].isin([f"{g}_{s}" for g in GCMS])].groupby("year")[["ET_mm", "P_mm", "T_C"]].mean()
        for i, (var, name, xlabel) in enumerate(drivers):
            ax = axs[i, j]
            x, y = ens[var].to_numpy(), ens["ET_mm"].to_numpy()
            fit = stats.linregress(x, y)
            x0, x1, _ = nice_limits(x)
            y0, y1, dy = nice_limits(y)
            ax.scatter(x, y, s=14, color="#1f77b4", alpha=0.85, edgecolor="black", linewidth=0.4,
                       label="Annual data", zorder=3)
            xx = np.linspace(x.min(), x.max(), 100)
            ax.plot(xx, fit.intercept + fit.slope * xx, color="#1f77b4", lw=1.3, label="Linear regression", zorder=4)
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
            p = "p < 0.001" if fit.pvalue < 0.001 else f"p = {fit.pvalue:.3f}"
            sign = "+" if fit.intercept >= 0 else "−"
            text = (f"Pearson r = {fit.rvalue:.3f}\nR² = {fit.rvalue ** 2:.3f}\n{p}\n"
                    f"y = {fit.slope:.3f}x {sign} {abs(fit.intercept):.3f}").replace("-", "−")
            box = ax.text(0.04, 0.96, text, transform=ax.transAxes,
                          ha="left", va="top", fontsize=5.6, linespacing=1.2, zorder=5,
                          bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="black", lw=0.5))
            leg = ax.legend(loc="lower right", fontsize=5.4, frameon=True, framealpha=1, edgecolor="#bdbdbd",
                            handlelength=1.5, borderpad=0.35, labelspacing=0.3)
            ax.grid(True, ls="--", color="#d9d9d9", lw=0.5, zorder=0)
            for sp in ax.spines.values():
                sp.set_visible(True)
                sp.set_linewidth(0.6)
            ax.tick_params(labelsize=6.5)
            ax.set_title(f"Actual ET vs {name} - {sname}", fontsize=6.8, fontweight="bold", pad=4)
            ax.set_xlabel(xlabel, fontsize=7)
            ax.set_ylabel("Actual evapotranspiration (mm)", fontsize=7)
            ax.text(-0.24, 1.08, f"({'abcdef'[3 * i + j]})", transform=ax.transAxes, fontsize=9,
                    fontweight="bold", ha="left", va="bottom")
            panels.append((ax, x, y, box, leg, dy))
    fig.canvas.draw()
    for ax, x, y, box, leg, dy in panels:
        to_axes = ax.transAxes.inverted()
        extents = [box.get_bbox_patch().get_window_extent(), leg.get_window_extent()]
        clear_boxes(ax, x, y, [e.transformed(to_axes) for e in extents], dy)
    fig.savefig(FIGURES / "fig9_et_climate.png")
    plt.close(fig)


FIGS = {"fig2": fig2, "fig3": fig3, "fig4": fig4, "fig5": fig5, "fig6": fig6, "fig7": fig7, "fig8": fig8,
        "fig9": fig9, "figS1": figS1}
for name in sys.argv[1:] or FIGS:
    FIGS[name]()
    print("done", name)
