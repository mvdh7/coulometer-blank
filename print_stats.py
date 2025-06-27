# %%
from scipy import stats

from read_datasets import dbs

l1 = stats.levene(
    dbs[dbs.nuts_good].dic_norm_cb.values, dbs[dbs.nuts_good].dic_norm_bh.values
)
l2 = stats.levene(
    dbs[dbs.nuts_good].dic_norm_cb.values, dbs[dbs.nuts_good].dic_norm.values
)
l3 = stats.levene(
    dbs[dbs.nuts_good].dic_norm_bh.values, dbs[dbs.nuts_good].dic_norm.values
)
print("BF test")
print(f"Const vs per-m: stat={l1.statistic:.1f}, p={l1.pvalue:.0e}")
print(f"Const vs fittd: stat={l2.statistic:.1f}, p={l2.pvalue:.0e}")
print(f"Per-m vs fittd: stat={l3.statistic:.2f}, p={l3.pvalue:.2f}")
print("")
print("Kurtosis")
print(f"Const: {stats.kurtosis(dbs[dbs.nuts_good].dic_norm_cb.values):.1f}")
print(f"Per-m: {stats.kurtosis(dbs[dbs.nuts_good].dic_norm_bh.values):.1f}")
print(f"Fittd: {stats.kurtosis(dbs[dbs.nuts_good].dic_norm.values):.1f}")
print("")
print("Range")
print(
    f"Const: {dbs[dbs.nuts_good].dic_norm_cb.max() - dbs[dbs.nuts_good].dic_norm_cb.min():.1f}"
    + f" ({dbs[dbs.nuts_good].dic_norm_cb.min():.1f} to {dbs[dbs.nuts_good].dic_norm_cb.max():.1f})"
)
print(
    f"Per-m: {dbs[dbs.nuts_good].dic_norm_bh.max() - dbs[dbs.nuts_good].dic_norm_bh.min():.1f}"
    + f" ({dbs[dbs.nuts_good].dic_norm_bh.min():.1f} to {dbs[dbs.nuts_good].dic_norm_bh.max():.1f})"
)
print(
    f"Fittd: {dbs[dbs.nuts_good].dic_norm.max() - dbs[dbs.nuts_good].dic_norm.min():.1f}"
    + f" ({dbs[dbs.nuts_good].dic_norm.min():.1f} to {dbs[dbs.nuts_good].dic_norm.max():.1f})"
)
