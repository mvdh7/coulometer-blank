# %%
import koolstof as ks
import numpy as np
from matplotlib import dates as mdates
from matplotlib import pyplot as plt

from read_datasets import dbs, sessions

marker = "o"
c = "xkcd:navy"
# Prepare to draw the figure
sample_types = {
    "nuts": "xkcd:strawberry",
    "junk": "xkcd:grass",
    "crm": "xkcd:azure",
}
# Get preliminary blank correction before QCing and re-doing
# for session in [sessions.index[0]]:
for session in sessions.index:
    s = sessions.loc[session]
    L = dbs[sessions.index.name] == session
    if dbs[L].nuts_good.any():
        fig, ax = plt.subplots(dpi=300, figsize=(5, 4))
        # Create and draw fitted line
        fx = np.linspace(
            dbs[L].datenum_analysis_scaled.min(),
            dbs[L].datenum_analysis_scaled.max(),
            500,
        )
        fy = ks.blank._blank_progression(s.blank_progression, fx)
        fx = mdates.num2date(
            ks.blank._de_centre_and_scale(
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
        fig.savefig(f"figures/plot_sessions_all/psa_{session}.png")
        plt.close()
