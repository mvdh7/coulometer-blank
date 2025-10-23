# %%
import matplotlib.pyplot as plt
import numpy as np

from read_datasets import dbs

fig, ax = plt.subplots(dpi=300, figsize=(5, 4))
ax.hist(
    dbs[dbs.nuts_usable].counts_u_blank_pct,
    # dbs[dbs.nuts_usable].blank_u_full_pct,
    bins=np.arange(0, 4, 0.01),
    facecolor="xkcd:dark",
    alpha=0.8,
)
ax.grid(alpha=0.2)
ax.set_ylabel("Number of measurements")
ax.set_xlabel("$σ(C)$ due to blank / %")
ax.set_xlim(0, 0.8)
fig.tight_layout()
fig.savefig("figures/plot_uncertainty.png")
fig.savefig("figures/figures-final/fig04.pdf")
