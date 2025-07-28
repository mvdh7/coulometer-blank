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
labels = {
    "nuts": "Internal standard",
    "junk": "Junk",
    "crm": "CRM",
}
markers = {"nuts": "s", "crm": "d", "junk": "v"}
fig, axs = plt.subplots(dpi=300, figsize=(7.5, 6), nrows=2, ncols=2)
fsessions = {
    "a": "C_Dec07-22_0912",
    "b": "C_Feb01-22_0802",
    "c": "C_Nov15-24_0811",
    "d": "C_Jun22-21_0806",
}
for i, (lcl, session) in enumerate(fsessions.items()):
    # do_legend = i == 1
    do_legend = False
    ax = axs.ravel()[i]
    s = sessions.loc[session]
    L = dbs[sessions.index.name] == session
    if dbs[L].nuts_good.any():
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
        ax.plot(fx, fy, c=c, label="Fitted blank", alpha=0.7)
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
            legend=do_legend,
        )
        for st, stc in sample_types.items():
            dbs[L & dbs.blank_good & dbs[st]].plot.scatter(
                "datetime_analysis",
                "blank_here",
                ax=ax,
                c=stc,
                marker=markers[st],
                label=labels[st],
                legend=do_legend,
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
                legend=do_legend,
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
        if do_legend:
            ax.legend(
                loc=3,
                bbox_to_anchor=(-0.5, 1.1),
                ncol=3,
                edgecolor="k",
                # mode="expand",
            )
        ax.xaxis.set_major_locator(mdates.HourLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H"))
        ax.set_xlabel("Analysis time of day")
        ax.set_ylabel(r"Coulometer blank / counts min$^{-1}$")
        ax.text(0, 1.03, f"({lcl}) {session}", transform=ax.transAxes)
        ax.grid(alpha=0.2)
fig.tight_layout()
fig.savefig("figures/plot_sessions_tidy4.png")
# fig.savefig("figures/plot_sessions_tidy4_legend.png")
