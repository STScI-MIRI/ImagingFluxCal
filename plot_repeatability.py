import argparse
import numpy as np

import matplotlib.pyplot as plt

from astropy.modeling import models, fitting
from astropy.table import QTable

from calc_calfactors import get_calfactors

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
               "F1500W", "F1800W", "F2100W", "F2550W"]

    # make plot
    fontsize = 16
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))

    ax = axs[0]
    startday = 59720.
    for k, cfilter in enumerate(filters):
        if cfilter == "F2550W":
            bkgsub = True
            rstr = "_bkgsub"
        else:
            bkgsub = False
            rstr = ""

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

        # now get the two stars that we repeated twice to fill in the time gap

        for stype, sname in zip(["SolarAnalogs", "ADwarfs"], ["HD 37962", "del UMi"]):
            cfacs2 = get_calfactors(
                stype,
                cfilter,
                xaxisval="timemid",
                startday=startday,
                bkgsub=bkgsub,
            )       
            nvals = [pname == sname for pname in cfacs2[4]]
            if (np.sum(nvals) > 1):
                nxvals = cfacs2[2][nvals]
                nyvals = cfacs2[0][nvals]
                nyvals_unc = cfacs2[1][nvals]
                #print(sname, cfilter)
                #print(nyvals)
                #exit()
                # find the value of BD+60 1753 that is closest to the last value
                # use this value to adjust the new star to be on the same scale
                #    later interpolate between values
                sindxs = np.argsort(np.absolute(nxvals[-1] - xvals))
                nyvals *= yvals[sindxs[0]] / nyvals[-1]

                yvals = np.append(yvals, nyvals)
                yvals_unc = np.append(yvals_unc, nyvals_unc)
                xvals = np.append(xvals, nxvals)

        # fit an exponential
        # ignore the bad data point for F770W
        gvals = abs(xvals - (60070. - startday)) > 20.

        fit = fitting.LevMarLSQFitter()
        mod_init = (models.Exponential1D(tau=-150., amplitude=-0.2) 
                    + models.Const1D(amplitude=0.70))
        mod_init[0].amplitude.bounds = [None, 0.0]
        if cfilter in ["F560W", "F770W", "F1000W", "F1130W", "F1280W", "F1500W"]:
            mod_init[0].tau.fixed = True
        else:
            mod_init[0].tau.bounds = [-400., -100.]
        fitx = xvals[gvals]
        fity = yvals[gvals]
        sindxs = np.argsort(fitx)
        mod_fit = fit(mod_init, fitx[sindxs], fity[sindxs])
        print(cfilter, mod_fit[0].tau.value)

        per_dev = (mod_fit(fitx) - fity) / mod_fit(fitx)
        per_dev = 100.0 * np.sqrt(np.sum(np.square(per_dev) / (len(fitx) - 2)))

        mod_dev = (mod_fit(fitx) - fity)
        mod_dev = np.sqrt(np.sum(np.square(mod_dev) / (len(fitx) - 2)))

        # save the fit results
        atab = QTable()
        atab[f"fit_exp_amp_{cfilter}"] = [mod_fit[0].amplitude.value]
        atab[f"fit_exp_tau_{cfilter}"] = [mod_fit[0].tau.value]
        atab[f"fit_exp_const_{cfilter}"] = [mod_fit[1].amplitude.value]
        atab[f"fit_exp_startday_{cfilter}"] = [startday]
        atab[f"fit_exp_std_{cfilter}"] = [mod_dev]
        sext = "_fit.dat"
        atab.write(
            f"CalFacs/miri_calfactors{rstr}_repeat_{cfilter}_fit.dat",
            format="ascii.commented_header",
            overwrite=True,
        )

        pxvals = np.arange(0, max(fitx))

        per_amp = 100. * (mod_fit[0].amplitude.value / mod_fit[1].amplitude.value)

        meanval = mod_fit[1].amplitude.value
        yvals = meanval / yvals
        # yvals_unc /= np.nanmean(yvals)
        modvals = (meanval / mod_fit(pxvals))

        yoff = k * 0.25
        yoff2 = k * 0.05
        ax.errorbar(xvals, yvals+yoff, yerr=yvals_unc, fmt="ko")
        ax.plot([0., 700.], [1. + yoff, 1. + yoff], "k:", alpha=0.5)
        ax.plot(pxvals, modvals + yoff, "m-")

        axs[1].errorbar(xvals, yvals+yoff2, yerr=yvals_unc, fmt="ko")
        axs[1].plot([0., 700.], [1. + yoff2, 1. + yoff2], "k:", alpha=0.5)
        axs[1].plot(pxvals, modvals + yoff2, "m-")

        shifty = 0.05
        ax.text(425., 1. + yoff + shifty, cfilter)
        ax.text(0.0, yoff + shifty + modvals[0],
                rf"A={-1.*per_amp:.1f}% / $\tau$={-1.*mod_fit[0].tau.value:.0f} days / $\sigma$={per_dev:.1f}%",
                backgroundcolor="w",
                fontsize=0.8*fontsize)

        shifty2 = 0.01
        axs[1].text(550., 1. + yoff2 + shifty2, cfilter)

    ax.set_ylim(0.9, 3.5)
    ax.set_xlabel(f"Time [MJD] - {startday} [days]")
    ax.set_ylabel("Fractional change (+ const)")

    axs[1].set_xlim(400., 700.)
    axs[1].set_ylim(0.98, 1.50)
    axs[1].set_xlabel(f"Time [MJD] - {startday} [days]")
    # axs[1].set_ylabel("Fractional change (+ const)")

    plt.tight_layout()

    fname = "all_repeatability"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
