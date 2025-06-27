# %%
import matplotlib.pyplot as plt
import numpy as np

from read_datasets import dbs

fig, ax = plt.subplots(dpi=300, figsize=(5, 4))
ax.hist(
    dbs[dbs.nuts_good].counts_u_blank_pct,
    bins=np.arange(0, 0.62, 0.01),
    facecolor="xkcd:dark",
    alpha=0.8,
)
ax.grid(alpha=0.2)
ax.set_ylabel("Frequency")
ax.set_xlabel("$σ(C)$ due to blank / %")
fig.tight_layout()
fig.savefig("figures/plot_uncertainty.png")
