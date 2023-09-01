import argparse
import numpy as np

import matplotlib.pyplot as plt

from astropy.modeling import models, fitting

from calc_calfactors import get_calfactors

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
               "F1500W", "F1800W", "F2100W", "F2550W"]

    # make plot
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))

    startday = 59720.
    for k, cfilter in enumerate(filters):
        if cfilter == "F2550W":
            bkgsub = True
        else:
            bkgsub = False

        cfacs = get_calfactors(
            "ADwarfs",
            cfilter,
            xaxisval="timemid",
            repeat=True,
            startday=startday,
            bkgsub=bkgsub,
        )
        yvals = cfacs[0]
        yvals_unc = cfacs[1]
        xvals = cfacs[2]

        # fit an exponential
        # ignore the bad data point for F770W
        gvals = abs(xvals - (60070. - startday)) > 20.

        fit = fitting.LevMarLSQFitter()
        mod_init = (models.Exponential1D(tau=-200., amplitude=-0.2) 
                    + models.Const1D(amplitude=0.70))
        mod_init[0].amplitude.bounds = [None, 0.0]
        mod_init[0].tau.fixed = True
        # mod_init[1].amplitude.fixed = True
        fitx = xvals[gvals]
        fity = yvals[gvals]
        sindxs = np.argsort(fitx)
        mod_fit = fit(mod_init, fitx[sindxs], fity[sindxs])

        per_dev = (mod_fit(fitx) - fity) / mod_fit(fitx)
        per_dev = 100.0 * np.sqrt(np.sum(np.square(per_dev) / (len(fitx) - 2)))

        mod_dev = (mod_fit(fitx) - fity)
        mod_dev = np.sqrt(np.sum(np.square(mod_dev) / (len(fitx) - 2)))

        # pxvals = np.arange(min(fitx), max(fitx))
        pxvals = np.arange(0.0, 500.)

        per_amp = 100. * (mod_fit[0].amplitude.value / mod_fit[1].amplitude.value)

        meanval = mod_fit[1].amplitude.value
        yvals = meanval / yvals
        # yvals_unc /= np.nanmean(yvals)

        yoff = k * 0.25
        ax.plot(xvals, yvals+yoff, "ko")
        ax.plot([0., 500.], [1. + yoff, 1. + yoff], "k:", alpha=0.5)

        modvals = (meanval / mod_fit(pxvals))
        ax.plot(pxvals, modvals + yoff, "m-")

        shifty = 0.05
        ax.text(450., 1. + yoff + shifty, cfilter)
        ax.text(0.0, yoff + shifty + modvals[0],
                f"A={-1.*per_amp:.1f}% / $\sigma$={per_dev:.1f}%",
                backgroundcolor="w",
                fontsize=0.8*fontsize)

    ax.set_ylim(0.9, 3.5)
    ax.set_xlabel(f"Time [MJD] - {startday}")
    ax.set_ylabel("Fractional change (+ const)")

    plt.tight_layout()

    fname = "all_repeatability"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
