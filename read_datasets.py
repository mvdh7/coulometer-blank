# %%
import numpy as np
import pandas as pd
from calkulate.density import seawater_1atm_MP81

dbs = pd.read_parquet("data/ecb-dbs.parquet")
sessions = pd.read_parquet("data/ecb-sessions.parquet")

# Get session stats for blanks
dbs["blank_diff"] = dbs.blank - dbs.blank_here
sessions["blank_rmsd"] = (
    dbs[["dic_cell_id", "blank_diff"]][dbs.blank_good]
    .groupby("dic_cell_id")
    .apply(
        lambda g: np.sqrt((g.blank_diff.abs() ** 2).mean()),
        include_groups=False,
    )
)
dbs["blank_rmsd"] = sessions.blank_rmsd.loc[dbs.dic_cell_id].values
dbs["counts_u_blank"] = dbs.run_time * dbs.blank_rmsd
dbs["counts_u_blank_pct"] = 100 * dbs.counts_u_blank / dbs.counts_corrected
dbs["dic_cell_number"] = (dbs.dic_cell_id != dbs.dic_cell_id.shift(1)).cumsum()

# Fill in missing salinity and DIC values with assumed salinity
dbs.loc[dbs.salinity.isnull(), "salinity"] = 35
dbs["density_analysis_dic"] = seawater_1atm_MP81(
    temperature=dbs.temperature_analysis_dic,
    salinity=dbs.salinity,
)
dbs["dic"] = dbs.counts_corrected * dbs.k_dic / dbs.density_analysis_dic

# Recalculate DIC with a constant blank (of 40) and with blank_here
dbs["dic_cb"] = (dbs.counts - dbs.run_time * 40) * dbs.k_dic / dbs.density_analysis_dic
dbs["dic_bh"] = (
    (dbs.counts - dbs.run_time * dbs.blank_here) * dbs.k_dic / dbs.density_analysis_dic
)

# Remove large outliers
dbs["nuts_good"] = (
    dbs.nuts
    & dbs.dic.notnull()
    & ~dbs.name_anon.isin(
        [
            "nuts-0012",
            "nuts-0018",
            "nuts-0083",
            "nuts-0107",
            "nuts-0108",
            "nuts-0217",
            "nuts-0251",
            "nuts-0258",
            "nuts-0270",
            "nuts-0280",
        ]
    )
)

# Get session stats and normalised DIC for the nuts samples
sessions["hours_range"] = np.nan
for session in sessions.index:
    L = (dbs.dic_cell_id == session) & dbs.nuts_good
    if L.sum() > 0:
        dbs.loc[L, "dic_norm"] = dbs[L].dic - dbs[L].dic.mean()
        dbs.loc[L, "dic_norm_cb"] = dbs[L].dic_cb - dbs[L].dic_cb.mean()
        dbs.loc[L, "dic_norm_bh"] = dbs[L].dic_bh - dbs[L].dic_bh.mean()
        sessions.loc[session, "hours_range"] = (
            dbs[L].datenum_analysis.max() - dbs[L].datenum_analysis.min()
        ) * 24


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
