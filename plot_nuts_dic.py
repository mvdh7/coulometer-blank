# %%
from matplotlib import pyplot as plt

from read_datasets import dbs

fig, ax = plt.subplots(dpi=300, figsize=(5, 6))
ax.scatter(
    "datetime_analysis",
    "dic",
    data=dbs[dbs.nuts & ~dbs.nuts_good],
    s=15,
    alpha=0.8,
    c="xkcd:strawberry",
    marker="x",
    label="Removed in QC",
)
ax.scatter(
    "datetime_analysis",
    "dic",
    data=dbs[dbs.nuts & dbs.nuts_good & ~dbs.nuts_usable],
    s=15,
    alpha=0.9,
    c="none",
    edgecolor="xkcd:sea blue",
    marker="s",
    label="Not used (one in session)",
)
ax.scatter(
    "datetime_analysis",
    "dic",
    data=dbs[dbs.nuts_usable],
    s=20,
    alpha=0.5,
    c="xkcd:dark",
    edgecolor="none",
    label="Used for analysis",
)
ax.grid(alpha=0.2)
ax.set_xlabel("Analysis date")
ax.set_ylabel("DIC / µmol kg$^{-1}$")
ax.legend(loc="lower left")
fig.tight_layout()
fig.savefig("figures/plot_nuts_dic.png")
