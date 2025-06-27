# %%
import pickle

import koolstof as ks
import numpy as np
from calkulate.density import seawater_1atm_MP81
from current_version import current_version
from matplotlib import dates as mdates
from matplotlib import pyplot as plt
from scipy import stats

# Import nutrient subsamples
with open("pickles/dbs_tdata_pdata_{}.pkl".format(current_version), "rb") as f:
    dbs = pickle.load(f)[0]
dbs["nuts"] = dbs.bottle.str.lower().str.startswith("nuts")
dbs["junk"] = dbs.bottle.str.lower().str.startswith("junk")
dbs["blank_diff"] = dbs.blank - dbs.blank_here

# Get session stats for blanks
dbs_sessions = (
    dbs[["dic_cell_id", "blank_diff"]][dbs.blank_good]
    .groupby("dic_cell_id")
    .apply(
        lambda g: np.sqrt((g.blank_diff.abs() ** 2).mean()),
        include_groups=False,
    )
)
dbs["blank_rmse"] = dbs_sessions.loc[dbs.dic_cell_id].values
dbs["counts_u_blank"] = dbs.run_time * dbs.blank_rmse
dbs["counts_u_blank_pct"] = 100 * dbs.counts_u_blank / dbs.counts_corrected

nuts = dbs[dbs.nuts].copy()

# %% Add salinity value and convert counts to DIC
# TODO: use this salinity value in the main data1.py script
nuts["salinity"] = 32
nuts["density_analysis_dic"] = seawater_1atm_MP81(
    temperature=nuts.temperature_analysis_dic,
    salinity=nuts.salinity,
)
nuts["dic"] = nuts.counts_corrected * nuts.k_dic / nuts.density_analysis_dic

# Recalculate DIC with a constant blank (of 40) and with blank_here
nuts["dic_cb"] = (
    (nuts.counts - nuts.run_time * 40) * nuts.k_dic / nuts.density_analysis_dic
)
nuts["dic_bh"] = (
    (nuts.counts - nuts.run_time * nuts.blank_here)
    * nuts.k_dic
    / nuts.density_analysis_dic
)

# QC
nuts["nuts_good"] = ~(
    ((nuts.bottle == "Nutslab1") & (nuts.dic_cell_id == "C_Nov20-18_0811"))
    | ((nuts.bottle == "NUTSLAB3") & (nuts.dic_cell_id == "C_May02-19_0705"))
    | ((nuts.bottle == "NUTSLAB8") & (nuts.dic_cell_id == "C_Nov25-21_0711"))
    | ((nuts.bottle == "NUTSLAB2") & (nuts.dic_cell_id == "C_Aug01-22_0808"))
    | ((nuts.bottle == "NUTSLAB3") & (nuts.dic_cell_id == "C_Aug01-22_0808"))
    | ((nuts.bottle == "NUTS1") & (nuts.dic_cell_id == "C_Apr15-24_0804"))
    | ((nuts.bottle == "NUTS7") & (nuts.dic_cell_id == "C_Nov12-24_0811"))
    | ((nuts.bottle == "NUTS14") & (nuts.dic_cell_id == "C_Nov15-24_0811"))
    | ((nuts.bottle == "NUTS2") & (nuts.dic_cell_id == "C_Jun12-25_0906"))
    | ((nuts.bottle == "NUTS12") & (nuts.dic_cell_id == "C_Jun17-25_0806"))
)

# Get sessions and plot
nuts["dic_cell_number"] = (nuts.dic_cell_id != nuts.dic_cell_id.shift(1)).cumsum()
sessions = nuts[nuts.nuts_good].groupby("dic_cell_number")
fig, ax = plt.subplots(dpi=300)
dt = []
for dcn, session in sessions:
    nuts.loc[session.index, "dic_norm"] = session.dic - session.dic.mean()
    nuts.loc[session.index, "dic_norm_cb"] = session.dic_cb - session.dic_cb.mean()
    nuts.loc[session.index, "dic_norm_bh"] = session.dic_bh - session.dic_bh.mean()
    ax.scatter(
        session.index,
        session.dic - session.dic.mean(),
        s=10,
    )
    # L = session.dic - session.dic.mean() > 10
    # if L.any():
    #     print(session.dic_cell_id.iloc[0])
    #     print(session[L].bottle)
    if session.nuts_good.any():
        dt.append(session.datetime_analysis.max() - session.datetime_analysis.min())

ax.axhline(0, c="k", lw=0.8)
ax.grid(alpha=0.2)
# ax.set_ylim(-15, 15)


# %% Histograms
def std_Sn(a):
    # https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/diffsn.htm
    #
    # Peter J. Rousseuw and Christophe Croux (1993), "Alternatives to the Median
    # Absolute Deviation," Journal of the American Statistical Association, Vol. 88,
    # No. 424, pp. 1273-1283.
    a = a[~np.isnan(a)].ravel()
    sn0 = np.full((np.size(a), np.size(a)), np.nan)
    for i in range(len(a)):
        sn0[i] = np.abs(a[i] - a)
    np.fill_diagonal(sn0, np.nan)
    sn1 = np.nanmedian(sn0, axis=1)
    sn = 1.1926 * np.median(sn1)
    return sn


fnuts = nuts[nuts.nuts_good & nuts.dic.notnull()]
bins = np.arange(-15, 15.1, 1)
fig, ax = plt.subplots(dpi=300, figsize=(8, 4))
ax.hist(
    fnuts.dic_norm_cb,
    bins=bins,
    alpha=0.7,
    facecolor="xkcd:cerulean",
    label="Constant blank",
    zorder=1,
)
ax.hist(
    fnuts.dic_norm_bh,
    bins=bins,
    alpha=0.8,
    facecolor="xkcd:royal purple",
    label="Per-sample blank",
    zorder=0,
)
ax.hist(
    fnuts.dic_norm,
    bins=bins,
    alpha=0.6,
    facecolor="xkcd:tangerine",
    label="Fitted blank",
    zorder=2,
)
bin_width = (bins.max() - bins.min()) / bins.size
fx = np.linspace(-12, 12, num=1000)
ax.plot(
    fx,
    stats.norm.pdf(
        fx,
        # scale=fnuts.dic_norm_cb.std(),
        scale=std_Sn(fnuts.dic_norm_cb.values),
    )
    * nuts.nuts_good.sum()
    * bin_width,
    c="xkcd:cerulean",
    label=r"$\sigma(S_n)$ = "
    + f"{std_Sn(fnuts.dic_norm_cb.values):.2f}"
    + r" µmol kg$^{-1}$",
    alpha=0.8,
    zorder=4,
    lw=2,
)
ax.plot(
    fx,
    stats.norm.pdf(
        fx,
        # scale=fnuts.dic_norm.std(),
        scale=std_Sn(fnuts.dic_norm_bh.values),
    )
    * nuts.nuts_good.sum()
    * bin_width,
    c="xkcd:royal purple",
    label=r"$\sigma(S_n)$ = "
    + f"{std_Sn(fnuts.dic_norm_bh.values):.2f}"
    + r" µmol kg$^{-1}$",
    alpha=0.8,
    zorder=3,
    lw=2,
)
ax.plot(
    fx,
    stats.norm.pdf(
        fx,
        # scale=fnuts.dic_norm.std(),
        scale=std_Sn(fnuts.dic_norm.values),
    )
    * nuts.nuts_good.sum()
    * bin_width,
    c="xkcd:tangerine",
    label=r"$\sigma(S_n)$ = "
    + f"{std_Sn(fnuts.dic_norm.values):.2f}"
    + r" µmol kg$^{-1}$",
    alpha=0.8,
    zorder=5,
    lw=2,
)
ax.axvline(0, c="k", lw=0.8)
ax.set_xlabel("∆DIC / µmol kg$^{-1}$")
ax.set_ylabel("Frequency")
ax.legend()
ax.grid(alpha=0.2)
fig.tight_layout()
fig.savefig("figures/test-koolstof-histogram.png")

# %%
l1 = stats.levene(fnuts.dic_norm_cb.values, fnuts.dic_norm_bh.values)
l2 = stats.levene(fnuts.dic_norm_cb.values, fnuts.dic_norm.values)
l3 = stats.levene(fnuts.dic_norm_bh.values, fnuts.dic_norm.values)
print("BF test")
print(f"Const vs per-m: stat={l1.statistic:.1f}, p={l1.pvalue:.0e}")
print(f"Const vs fittd: stat={l2.statistic:.1f}, p={l2.pvalue:.0e}")
print(f"Per-m vs fittd: stat={l3.statistic:.2f}, p={l3.pvalue:.2f}")
print("")
print("Kurtosis")
print(f"Const: {stats.kurtosis(fnuts.dic_norm_cb.values):.1f}")
print(f"Per-m: {stats.kurtosis(fnuts.dic_norm_bh.values):.1f}")
print(f"Fittd: {stats.kurtosis(fnuts.dic_norm.values):.1f}")
print("")
print("Range")
print(f"Const: {fnuts.dic_norm_cb.max() - fnuts.dic_norm_cb.min():.1f}")
print(f"Const: {fnuts.dic_norm_bh.max() - fnuts.dic_norm_bh.min():.1f}")
print(f"Const: {fnuts.dic_norm.max() - fnuts.dic_norm.min():.1f}")

# %%
fyl = 15  # 4 or 15
fig, ax = plt.subplots(dpi=300)
fx = range(fnuts.dic.size)
fy = fnuts.dic_norm.sort_values()
fy_cb = fnuts.dic_norm_cb.sort_values()
fy_bh = fnuts.dic_norm_bh.sort_values()
kws = dict()
ax.plot(
    fx,
    fy_cb,
    **kws,
    c="xkcd:cerulean",
    alpha=0.8,
    lw=2,
    label="Constant blank",
    zorder=1,
)
ax.plot(
    fx,
    fy_bh,
    **kws,
    c="xkcd:royal purple",
    alpha=0.8,
    lw=2,
    label="Per-sample blank",
    zorder=0,
)
ax.plot(
    fx,
    fy,
    **kws,
    c="xkcd:tangerine",
    alpha=0.8,
    lw=2,
    label="Fitted blank",
    zorder=2,
)
# kws = dict(s=10, edgecolor="none")
# ax.scatter(fx, fy, **kws, c="xkcd:tangerine", alpha=0.7)
# ax.scatter(fx, fy_cb, **kws, c="xkcd:cerulean", alpha=0.5)
# ax.scatter(fx, fy_bh, **kws, c="xkcd:kelly green", alpha=0.5)
ax.axhline(0, c="k", lw=0.8)
ax.set_ylim(-fyl, fyl)
ax.grid(alpha=0.2)
ax.set_ylabel("∆DIC / µmol kg$^{-1}$")
ax.legend()
ax.set_xlabel("Sorted measurement number")
fig.tight_layout()
fig.savefig(f"figures/test-koolstof-lines-{fyl}.png")


# %%
methods = [  # list of VINDTA methods that are considered as measurements
    "3C standard",
    "3C standardRWS",
]
vindta_path = "data/vindta15_data/"  # path to the logfile and .dbs files
logfile_fname = "logfile_20241115.bak"  # the logfile's filename
logfile = ks.vindta.read_logfile(
    vindta_path + logfile_fname,
    methods=methods,
    ignore_lines=[31691, 42776, 47717, 76543, 76560],
)
use_from = 6
# TODO: it appears that the figures don't change perceptibly when adjusting
# use_from to a different value?? Or maybe it's too small to notice...
sessions = ks.vindta.blank_correction(dbs, logfile, use_from=use_from)

# %%
marker = "o"
c = "xkcd:navy"
ax = None
# Prepare to draw the figure
sample_types = {
    "nuts": "xkcd:strawberry",
    "junk": "xkcd:grass",
    "crm": "xkcd:azure",
}
# Get preliminary blank correction before QCing and re-doing
session = sessions.index[80]
s = sessions.loc[session]
L = dbs[sessions.index.name] == session
if ax is None:
    fig, ax = plt.subplots(dpi=300)
# Create and draw fitted line
fx = np.linspace(
    dbs[L].datenum_analysis_scaled.min(), dbs[L].datenum_analysis_scaled.max(), 500
)
fy = ks.vindta.get._blank_progression(s.blank_progression, fx)
fx = mdates.num2date(
    ks.vindta.get._de_centre_and_scale(
        fx, s.datenum_analysis_std, s.datenum_analysis_mean
    )
)
ax.plot(fx, fy, c=c, label="Best fit")
# Draw errorbars
ax.errorbar(
    "datetime_analysis",
    "blank_here",
    yerr="blank_here_std",
    data=dbs[L & dbs.blank_good & ~dbs.nuts & ~dbs.crm & ~dbs.junk],
    alpha=0.3,
    ecolor=c,
    label=None,
    linestyle="none",
)
for st, stc in sample_types.items():
    ax.errorbar(
        "datetime_analysis",
        "blank_here",
        yerr="blank_here_std",
        data=dbs[L & dbs.blank_good & dbs[st]],
        alpha=0.3,
        ecolor=stc,
        label=None,
        linestyle="none",
    )
# Draw the rest of the figure
dbs[L & dbs.blank_good & ~dbs.nuts & ~dbs.crm & ~dbs.junk].plot.scatter(
    "datetime_analysis",
    "blank_here",
    ax=ax,
    c=c,
    marker=marker,
    label="Samples",
)
for st, stc in sample_types.items():
    dbs[L & dbs.blank_good & dbs[st]].plot.scatter(
        "datetime_analysis",
        "blank_here",
        ax=ax,
        c=stc,
        marker=marker,
        label=st,
    )
y_max = np.max([dbs[L & dbs.blank_good].blank_here.max(), np.max(fy)]) * 1.2
l_ignored = L & ~dbs.blank_good & (dbs.blank_here <= y_max)
if l_ignored.any():
    dbs[l_ignored].plot.scatter(
        "datetime_analysis",
        "blank_here",
        ax=ax,
        c="none",
        edgecolor=c,
        marker=marker,
        label="Ignored",
    )
l_offscale = L & (dbs.blank_here > y_max)
if l_offscale.any():
    off_x = dbs[l_offscale].datetime_analysis.values
    ax.scatter(
        off_x,
        np.full(np.size(off_x), y_max * 0.99999),
        c="none",
        edgecolor=c,
        marker="^",
        label="Off scale (ignored)",
        clip_on=False,
    )
ax.set_ylim([0, y_max])
ax.legend(edgecolor="k")
ax.xaxis.set_major_locator(mdates.HourLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H"))
ax.set_xlabel("Time of day of analysis")
ax.set_ylabel(r"Coulometer blank / cts min$^{-1}$")
ax.set_title(session)
ax.grid(alpha=0.2)
# add_credit(ax)
fig.tight_layout()

# %%
fig, ax = plt.subplots(dpi=300)
ax.scatter(
    "datetime_analysis",
    "dic",
    data=nuts[nuts.nuts_good],
    c="xkcd:dark",
    s=10,
    alpha=0.5,
)
ax.grid(alpha=0.2)
ax.set_xlabel("Analysis date")
ax.set_ylabel("Calibrated DIC / µmol kg$^{-1}$")

# %%
fig, ax = plt.subplots(dpi=300)
ax.hist(
    nuts[nuts.nuts_good].counts_u_blank_pct,
    bins=np.arange(0, 0.62, 0.01),
    facecolor="xkcd:dark",
    alpha=0.8,
)
ax.grid(alpha=0.2)
ax.set_ylabel("Frequency")
ax.set_xlabel("$σ(C)$ due to blank / %")
fig.tight_layout()
fig.savefig("figures/test-koolstof-uncert.png")

# %%
