# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from read_datasets import dbs, std_Sn

bins = np.arange(-12, 12.1, 1)
fig, axs = plt.subplots(nrows=2, dpi=300, figsize=(5, 7))

ax = axs[0]
ax.hist(
    dbs.dic_norm_cb,
    bins=bins,
    alpha=0.8,
    facecolor="xkcd:cerulean",
    edgecolor="xkcd:cerulean",
    label="Constant",
    zorder=1,
    histtype="stepfilled",
)
ax.hist(
    dbs.dic_norm_cb,
    bins=bins,
    alpha=0.8,
    facecolor="none",
    edgecolor="xkcd:cerulean",
    zorder=10,
    histtype="stepfilled",
)
ax.hist(
    dbs.dic_norm_bh,
    bins=bins,
    alpha=0.8,
    facecolor="xkcd:royal purple",
    edgecolor="xkcd:royal purple",
    label="Per-measurement",
    zorder=0,
    histtype="stepfilled",
)
ax.hist(
    dbs.dic_norm,
    bins=bins,
    alpha=0.5,
    facecolor="xkcd:tangerine",
    edgecolor="xkcd:tangerine",
    label="Fitted",
    zorder=2,
    histtype="stepfilled",
)
bin_width = (bins.max() - bins.min()) / bins.size
fx = np.linspace(-12, 12, num=1000)
ax.plot(
    fx,
    stats.norm.pdf(
        fx,
        # scale=dbs.dic_norm_cb.std(),
        scale=std_Sn(dbs.dic_norm_cb.values),
    )
    * dbs.nuts_usable.sum()
    * bin_width,
    c="xkcd:cerulean",
    # label=r"$\sigma(S_n)$ = "
    # + f"{std_Sn(dbs.dic_norm_cb.values):.2f}"
    # + r" µmol kg$^{-1}$",
    alpha=0.8,
    zorder=100,
    lw=2,
)
ax.plot(
    fx,
    stats.norm.pdf(
        fx,
        # scale=dbs.dic_norm_bh.std(),
        scale=std_Sn(dbs.dic_norm_bh.values),
    )
    * dbs.nuts_usable.sum()
    * bin_width,
    c="xkcd:royal purple",
    # label=r"$\sigma(S_n)$ = "
    # + f"{std_Sn(dbs.dic_norm_bh.values):.2f}"
    # + r" µmol kg$^{-1}$",
    alpha=0.8,
    zorder=99,
    lw=2,
)
ax.plot(
    fx,
    stats.norm.pdf(
        fx,
        # scale=dbs.dic_norm.std(),
        scale=std_Sn(dbs.dic_norm.values),
    )
    * dbs.nuts_usable.sum()
    * bin_width,
    c="xkcd:tangerine",
    # label=r"$\sigma(S_n)$ = "
    # + f"{std_Sn(dbs.dic_norm.values):.2f}"
    # + r" µmol kg$^{-1}$",
    alpha=0.7,
    zorder=101,
    lw=2,
)
ax.axvline(0, c="k", lw=0.8)
ax.set_xlabel("∆DIC / µmol kg$^{-1}$")
ax.set_ylabel("Number of measurements")
ax.legend(fontsize=9)
ax.grid(alpha=0.2)
ax.text(0, 1.03, "(a)", transform=ax.transAxes)

ax = axs[1]
fx = range(dbs[dbs.nuts_usable].dic.size)
fy = dbs[dbs.nuts_usable].dic_norm.sort_values()
fy_cb = dbs[dbs.nuts_usable].dic_norm_cb.sort_values()
fy_bh = dbs[dbs.nuts_usable].dic_norm_bh.sort_values()
kws = dict()
ax.plot(
    fx,
    fy_cb,
    **kws,
    c="xkcd:cerulean",
    alpha=0.8,
    lw=2,
    label="Constant",
    zorder=1,
)
ax.plot(
    fx,
    fy_bh,
    **kws,
    c="xkcd:royal purple",
    alpha=0.8,
    lw=2,
    label="Per-measurement",
    zorder=0,
)
ax.plot(
    fx,
    fy,
    **kws,
    c="xkcd:tangerine",
    alpha=0.8,
    lw=2,
    label="Fitted",
    zorder=2,
)
fxl = ax.get_xlim()
lw = 1
ax.fill_between(
    fxl, fy_bh.iloc[-1], fy_cb.iloc[-1], facecolor="xkcd:cerulean", alpha=0.2
)
ax.fill_between(
    fxl, fy.iloc[-1], fy_bh.iloc[-1], facecolor="xkcd:royal purple", alpha=0.2
)
ax.fill_between(fxl, fy.iloc[0], fy.iloc[-1], facecolor="xkcd:tangerine", alpha=0.2)
ax.fill_between(fxl, fy.iloc[0], fy_cb.iloc[0], facecolor="xkcd:cerulean", alpha=0.2)
ax.fill_between(
    fxl, fy_cb.iloc[0], fy_bh.iloc[0], facecolor="xkcd:royal purple", alpha=0.2
)
# ax.axhline(fy_cb.iloc[0], ls=":", lw=lw, c="xkcd:cerulean")
# ax.axhline(fy_cb.iloc[-1], ls=":", lw=lw, c="xkcd:cerulean")
# ax.axhline(fy_bh.iloc[0], ls="--", lw=lw, c="xkcd:royal purple")
# ax.axhline(fy_bh.iloc[-1], ls="--", lw=lw, c="xkcd:royal purple")
# ax.axhline(fy.iloc[0], ls="-", lw=lw, c="xkcd:tangerine", alpha=0.5)
# ax.axhline(fy.iloc[-1], ls="-", lw=lw, c="xkcd:tangerine", alpha=0.5)
ax.set_xlim(fxl)
# kws = dict(s=10, edgecolor="none")
# ax.scatter(fx, fy, **kws, c="xkcd:tangerine", alpha=0.7)
# ax.scatter(fx, fy_cb, **kws, c="xkcd:cerulean", alpha=0.5)
# ax.scatter(fx, fy_bh, **kws, c="xkcd:kelly green", alpha=0.5)
ax.axhline(0, c="k", lw=0.8)
ax.set_ylim(-15, 15)
ax.grid(alpha=0.2)
ax.set_ylabel("∆DIC / µmol kg$^{-1}$")
ax.legend(fontsize=9)
ax.set_xlabel("Sorted measurement number")
ax.text(0, 1.03, "(b)", transform=ax.transAxes)

fig.tight_layout()
fig.savefig("figures/plot_histograms.png")
