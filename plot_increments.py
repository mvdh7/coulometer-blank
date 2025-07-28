# %%
import numpy as np
from matplotlib import pyplot as plt

from read_datasets import dbs, logfile

lines = dict(c="xkcd:dark", alpha=0.05, lw=0.5)

increments = {k: [] for k in range(13)}
for i, row in logfile.iterrows():
    # dbsrow = dbs[dbs.name_anon == row.name_anon]
    # if dbsrow.shape[0] dbsrow.name_anon != "":
    if row.name_anon in dbs.name_anon.values:
        for k in increments:
            L = row.table["time"] == k
            if L.sum() == 1:
                incr = row.table["increments"][L][0]
                if incr < 10000000:  # some dodgy really high ones out there
                    increments[k].append(incr)
increments = {k: np.array(v) for k, v in increments.items()}

fig, axs = plt.subplots(dpi=300, nrows=2, figsize=(5, 6.5))
# vstyle = dict(facecolor="xkcd:dark")
for k, v in increments.items():
    parts0 = axs[0].violinplot(
        v * 1e-3,
        [k],
        showextrema=False,
        points=100,
        widths=0.8,
    )
    parts1 = axs[1].violinplot(
        v,
        [k],
        showextrema=False,
        points=5000,
        widths=0.8,
    )
    for parts in [parts0, parts1]:
        parts["bodies"][0].set_facecolor("xkcd:dark")
        parts["bodies"][0].set_alpha(0.4)
    axs[0].scatter(k, np.median(v) * 1e-3, marker="+", c="xkcd:dark")
    axs[1].scatter(k, np.median(v), marker="+", c="xkcd:dark")
axs[0].set_ylim([0, 200])
axs[1].set_ylim([0, 500])
for ax in axs:
    ax.set_xlim(0, 12.5)
    ax.set_xlabel("Titration time / min")
    ax.set_xticks(range(0, 13, 2))
axs[0].set_ylabel("Increments / 10$^3$ counts min$^{-1}$")
axs[1].set_ylabel("Increments / counts min$^{-1}$")
axs[0].text(0, 1.03, "(a)", transform=axs[0].transAxes)
axs[1].text(0, 1.03, "(b)", transform=axs[1].transAxes)
fig.tight_layout()
fig.savefig("figures/plot_increments.png")
