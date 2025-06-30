# %%
import numpy as np
from scipy import stats

from read_datasets import dbs, sessions, std_Sn

print("Counts")
print(f"{dbs.nuts.sum()} nuts measurements total")
print(f"{dbs.nuts_good.sum()} nuts measurements after QC")
print(f" (so {dbs.nuts.sum() - dbs.nuts_good.sum()} not used)")
print(f"{dbs.nuts_usable.sum()} nuts measurements usable for analysis")
print(" (i.e. where there was more than one nuts measurement per session)")
print(
    f"There were {(sessions.n_nuts > 0).sum()} analysis sessions over {dbs[dbs.nuts_usable].datetime_analysis.max() - dbs[dbs.nuts_usable].datetime_analysis.min()} from"
)
print(
    f"{dbs[dbs.nuts_usable].datetime_analysis.min()} to {dbs[dbs.nuts_usable].datetime_analysis.max()}"
)
print(f"but only {sessions.usable.sum()} sessions had usable data in the end.")
l1 = stats.levene(
    dbs[dbs.nuts_usable].dic_norm_cb.values, dbs[dbs.nuts_usable].dic_norm_bh.values
)
l2 = stats.levene(
    dbs[dbs.nuts_usable].dic_norm_cb.values, dbs[dbs.nuts_usable].dic_norm.values
)
l3 = stats.levene(
    dbs[dbs.nuts_usable].dic_norm_bh.values, dbs[dbs.nuts_usable].dic_norm.values
)
print("")
print("BF test")
print(f"Const vs per-m: stat={l1.statistic:.1f}, p={l1.pvalue:.0e}")
print(f"Const vs fittd: stat={l2.statistic:.1f}, p={l2.pvalue:.0e}")
print(f"Per-m vs fittd: stat={l3.statistic:.2f}, p={l3.pvalue:.2f}")
print("")
print("Kurtosis")
print(f"Const: {stats.kurtosis(dbs[dbs.nuts_usable].dic_norm_cb.values):.1f}")
print(f"Per-m: {stats.kurtosis(dbs[dbs.nuts_usable].dic_norm_bh.values):.1f}")
print(f"Fittd: {stats.kurtosis(dbs[dbs.nuts_usable].dic_norm.values):.1f}")
print("")
print("Range")
print(
    f"Const: {dbs[dbs.nuts_usable].dic_norm_cb.max() - dbs[dbs.nuts_usable].dic_norm_cb.min():.1f}"
    + f" ({dbs[dbs.nuts_usable].dic_norm_cb.min():.1f} to {dbs[dbs.nuts_usable].dic_norm_cb.max():.1f})"
)
print(
    f"Per-m: {dbs[dbs.nuts_usable].dic_norm_bh.max() - dbs[dbs.nuts_usable].dic_norm_bh.min():.1f}"
    + f" ({dbs[dbs.nuts_usable].dic_norm_bh.min():.1f} to {dbs[dbs.nuts_usable].dic_norm_bh.max():.1f})"
)
print(
    f"Fittd: {dbs[dbs.nuts_usable].dic_norm.max() - dbs[dbs.nuts_usable].dic_norm.min():.1f}"
    + f" ({dbs[dbs.nuts_usable].dic_norm.min():.1f} to {dbs[dbs.nuts_usable].dic_norm.max():.1f})"
)
print("")
print("Sn stdev")
print(f"Const: {std_Sn(dbs.dic_norm_cb.values):.2f}")
print(f"Per-m: {std_Sn(dbs.dic_norm_bh.values):.2f}")
print(f"Fittd: {std_Sn(dbs.dic_norm.values):.2f}")
print("Normal stdev")
print(f"Const: {np.nanstd(dbs.dic_norm_cb.values):.2f}")
print(f"Per-m: {np.nanstd(dbs.dic_norm_bh.values):.2f}")
print(f"Fittd: {np.nanstd(dbs.dic_norm.values):.2f}")
print("")
print("Uncertainty")
print(
    f"{(dbs[dbs.nuts_usable].counts_u_blank_pct < 0.1).sum() * 100 / dbs.nuts_usable.sum():.0f} % of samples have blank uncertainty < 0.1 %,"
)
print(f" with a median of {dbs[dbs.nuts_usable].counts_u_blank_pct.median():.2f} %")
