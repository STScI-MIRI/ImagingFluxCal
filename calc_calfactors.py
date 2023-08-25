from os.path import exists
import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np

from astropy.table import QTable
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.modeling import models, fitting


def get_calfactors(dir, filter, xaxisval="mflux", bkgsub=False, indivmos=False, indivcals=False, eefraction=0.7,
                   repeat=False, subtrans=False, startday=59720.):
    """
    Read in the observed and mdoel fluxes and computer the calibration factors
    """
    if bkgsub:
        extstr = "_bkgsub"
    else:
        extstr = ""
    if indivmos:
        extstr = "_indivmos"
    if indivcals:
        extstr = f"{extstr}_indivcals"
    # read in observed fluxes
    obstab = QTable.read(f"{dir}/{filter}{extstr}_eefrac{eefraction}_phot.fits")
    # read in model fluxes
    modtab = QTable.read("Models/model_phot.fits")

    if repeat:  # only use observations of BD+60 1753 and HD 2811
        # gvals = (obstab["name"] == "BD+60 1753") | (obstab["name"] == "HD 2811")
        gvals = (obstab["name"] == "BD+60 1753")
        obstab = obstab[gvals]
    elif subtrans:  # only use 2MASS J17571324+6703409 HJD around 59818.2693130625
        gvals = (obstab["name"] == "2MASS J17571324+6703409") & (abs(obstab["timemid"].value - 59818.3) < 1.)
        obstab = obstab[gvals]

    names = []
    xvals = []
    cfactors = []
    cfactors_unc = []
    subarrs = []
    for k, cname in enumerate(obstab["name"]):
        cfilter = obstab["filter"][k]
        oflux = obstab["aperture_sum_bkgsub"][k]
        oflux_unc = obstab["aperture_sum_bkgsub_err"][k]
        apcorr = obstab["apcorr"][k]
        pixarea = obstab["pixarea"][k]

        (mindx,) = np.where(modtab["name"] == cname)
        if len(mindx) < 1:
            print(f"Model fluxes for {cname} not present")
            exit()
        mflux = modtab[cfilter][mindx[0]]

        if xaxisval == "timemid":
            xval = obstab["timemid"][k]
        elif xaxisval == "rate":
            xval = obstab["pix_max"][k]
        elif xaxisval == "welldepth":
            xval = obstab["tgroup"][k] * obstab["ngroups"][k] * obstab["pix_max"][k]
        elif xaxisval == "bkg":
            xval = obstab["mean_bkg"][k]
        else:
            xval = mflux * 1e3

        if np.isfinite(oflux.value):
            cfactor = 1e-6 * mflux.value / (oflux.value * apcorr * pixarea.value)
            cfactor_unc = (oflux_unc / oflux) * cfactor
            # if obstab["name"][k] not in ["HD 167060", "16 Cyg B", "HD 37962", "del UMi"]:
            # print(obstab["name"][k], cfactor)
            names.append(obstab["name"][k])
            xvals.append(xval.value)
            cfactors.append(cfactor)
            cfactors_unc.append(cfactor_unc)
            subarrs.append(obstab["subarray"][k])

    if xaxisval == "timemid":
        xvals = np.array(xvals) - startday

    res = (np.array(cfactors), np.array(cfactors_unc), xvals, subarrs, names)
    return res


def plot_calfactors(
    ax,
    filter,
    xaxisval,
    showleg=True,
    savefile=None,
    applysubarrcor=True,
    showcurval=True,
    bkgsub=False,
    indivmos=False,
    indivcals=False,
    eefraction=0.7,
    repeat=False,
    subtrans=False,
):
    """
    Plot the calibration factors versus the requested xaxis.
    """
    dirs = ["HotStars", "ADwarfs", "SolarAnalogs"]
    if repeat or subtrans:
        dirs = ["ADwarfs"]
    pcols = ["b", "g", "r"]

    psubsym = {
        "FULL": "o",
        "BRIGHTSKY": "s",
        "SUB256": "p",
        "SUB128": "P",
        "SUB64": "*",
        "MASK1065": "^",
        "MASK1140": ">",
        "MASK1550": "<",
        "MASKLYOT": "v",
    }

    # efac = 1.04
    efac = 1.0
    # updated based on array-bkg subtraction reductions - better centroids (9 Mar 2023)       
    subarr_cor = {
        "FULL": 1.0,
        "BRIGHTSKY": 1.0005 * efac,
        "SUB256": 1.005 * efac,
        "SUB128": 1.005 * efac,
        "SUB64": 1.0111 * efac,
        "MASK1065": 1.0,
        "MASK1140": 1.0,
        "MASK1550": 1.0,
        "MASKLYOT": 1.0,
    }

    # print(subarr_cor)
    ignore_names = ["HD 167060", "16 Cyg B", "HD 37962", "del UMi", "HD 106252", "HD 142331"]
    modfac = {"HD 167060": 1.0/1.09,
              "16 Cyg B": 1.0/1.08,
              "HD 37962": 1.0/1.06,
              "del UMi": 1.0/1.03,
              "HD 106252": 1.0/1.06,
              "HD 142331": 1.0/1.05,
              "BD+60 1753": 1.0}
    # modfac = {"HD 167060": 1.0,
    #           "16 Cyg B": 1.0,
    #           "HD 37962": 1.0,
    #           "del UMi": 1.0,
    #           "HD 180609": 1.0}

    startday = 59720.

    allfacs = []
    allfacuncs = []
    allnames = []
    xvals = []
    for k, dir in enumerate(dirs):
        print(f"{dir}/{filter}_eefrac{eefraction}_phot.fits")
        if exists(f"{dir}/{filter}_eefrac{eefraction}_phot.fits"):
            cfacs = get_calfactors(
                dir, filter, xaxisval=xaxisval, bkgsub=bkgsub,
                indivmos=indivmos, indivcals=indivcals,
                eefraction=eefraction, repeat=repeat, subtrans=subtrans,
                startday=startday,
            )
            # allfacs.append(cfacs[0])
            allfacuncs.append(cfacs[1])
            xvals.append(cfacs[2])
            allnames.append(cfacs[4])
            for cfactor, cfactor_unc, xval, subarray, cname in zip(
                cfacs[0], cfacs[1], cfacs[2], cfacs[3], cfacs[4]
            ):
                # print(cname, xval, cfactor)
                ax.errorbar(
                    [xval],
                    [cfactor],
                    yerr=[cfactor_unc],
                    fmt=f"{pcols[k]}{psubsym[subarray]}",
                    alpha=0.1,
                    markersize=10,
                )
                if applysubarrcor:
                    cfactor = cfactor / subarr_cor[subarray]
                allfacs.append(cfactor)
                ax.errorbar(
                    [xval],
                    [cfactor],
                    yerr=[cfactor_unc],
                    fmt=f"{pcols[k]}{psubsym[subarray]}",
                    alpha=0.5,
                    markersize=10,
                )
                # plot a red circle around those not used in the average
                if ((cname in ignore_names) or
                    ((cname == "BD+60 1753") and (abs(xval - 60070.) < 20.))):
                    ax.scatter(
                        [xval], [cfactor], s=150, facecolor="none", edgecolor="m",
                    )
                    #print(modfac[cname])
                    ax.scatter(
                        [xval], [cfactor * modfac[cname]], s=150, facecolor="k", edgecolor="m",
                    )                    
                if subarray == "FULL":
                    meanfull = cfactor
            # special code to give the differneces between the subarrays
            if args.nosubarrcor:
                print(cfacs[3])
                print(cfacs[0] / meanfull)
            # exit()
    # allfacs = np.concatenate(allfacs)
    allfacs = np.array(allfacs)
    allfacuncs = np.concatenate(allfacuncs)
    medval = np.nanmedian(allfacs)
    allnames = np.concatenate(allnames)
    xvals = np.concatenate(xvals)

    # print the top 4 calibration factors with names
    # aindxs = np.flip(np.argsort(allfacs))
    # print(allfacs[aindxs[0:4]])
    # print(allnames[aindxs[0:4]])

    # ignore the 4 "high" stars
    gvals = []
    for k, cname in enumerate(allnames):
        if ((cname in ignore_names) or
                    ((cname == "BD+60 1753") and (abs(xvals[k] - 60070.) < 20.))):
            gvals.append(False)
        else:
            gvals.append(True)

    meanvals = sigma_clipped_stats(allfacs[gvals], sigma=3, maxiters=5)
    ax.axhline(y=meanvals[0], color="k", linestyle="-", alpha=0.5)
    medval = meanvals[0]
    filtered_data = sigma_clip(allfacs[gvals], sigma=3, maxiters=5)
    ax.plot((xvals[gvals])[filtered_data.mask], (allfacs[gvals])[filtered_data.mask], "x", color='#d62728')
    # ax.axhline(y=meanvals[0] + meanvals[2], color="k", linestyle=":", alpha=0.5)
    # ax.axhline(y=meanvals[0] - meanvals[2], color="k", linestyle=":", alpha=0.5)
    perstd = 100.0 * meanvals[2] / meanvals[0]
    # print(meanvals[0], meanvals[2], perstd)

    if xaxisval == "timemid":
        fit = fitting.LevMarLSQFitter()
        mod_init = models.Exponential1D(tau=-200., amplitude=-0.2) + models.Const1D(amplitude=0.70)
        mod_init[0].amplitude.bounds = [None, 0.0]
        mod_init[0].tau.fixed = True
        # mod_init[1].amplitude.fixed = True
        fitx = ((xvals[gvals])[~filtered_data.mask])
        fity = (allfacs[gvals])[~filtered_data.mask]
        sindxs = np.argsort(fitx)
        mod_fit = fit(mod_init, fitx[sindxs], fity[sindxs])
        per_dev = (mod_fit(fitx) - fity) / mod_fit(fitx)
        per_dev = 100.0 * np.sqrt(np.sum(np.square(per_dev) / (len(fitx) - 2)))

        mod_dev = (mod_fit(fitx) - fity)
        mod_dev = np.sqrt(np.sum(np.square(mod_dev) / (len(fitx) - 2)))

        pxvals = np.arange(min(fitx), max(fitx))
        ax.plot(pxvals, mod_fit(pxvals), "m-")

        per_amp = 100. * (mod_fit[0].amplitude.value / mod_fit[1].amplitude.value)
        ax.text(
            0.95,
            0.02,
            fr"fit: {mod_fit[0].amplitude.value:.3f} exp(x/{mod_fit[0].tau.value:.1f}) + {mod_fit[1].amplitude.value:.3f} ({per_dev:.1f}%) [amp: {per_amp:.1f}%]",
            transform=ax.transAxes,
            ha="right",
        )

    ax.text(
        0.95,
        0.08,
        fr"average: {meanvals[0]:.3f} $\pm$ {meanvals[2]:.3f} ({perstd:.1f}%)",
        transform=ax.transAxes,
        ha="right",
    )

    if savefile is not None:
        atab = QTable()
        atab[f"avecalfac_{filter}"] = [meanvals[0]]
        atab[f"avecalfac_unc_{filter}"] = [meanvals[2] / np.sqrt(len(allfacs[gvals]))]
        atab[f"avecalfac_std_{filter}"] = [meanvals[2]]
        if xaxisval == "timemid":
            atab[f"fit_exp_amp_{filter}"] = [mod_fit[0].amplitude.value]
            atab[f"fit_exp_tau_{filter}"] = [mod_fit[0].tau.value]
            atab[f"fit_exp_const_{filter}"] = [mod_fit[1].amplitude.value]
            atab[f"fit_exp_std_{filter}"] = [mod_dev]
            sext = "_fit.dat"
        else:
            sext = "_ave.dat"
        atab.write(
            savefile.replace(".fits", sext),
            format="ascii.commented_header",
            overwrite=True,
        )

        otab = QTable()
        otab["name"] = allnames
        otab[f"calfac_{filter}"] = allfacs
        otab[f"calfac_{filter}_mean_dev"] = allfacs / meanvals[0]
        otab[f"calfac_{filter}_med_dev"] = allfacs / meanvals[1]
        otab[f"modflux_{filter}"] = xvals
        otab.write(savefile, overwrite=True)

    # get the current pipeline calibration factor and plot
    if showcurval:
        cftab = QTable.read("CalFactors/jwst_miri_photom_0079.fits")
        pipe_cfactor = cftab["photmjsr"][cftab["filter"] == filter.split("_")[0]][0]
        ax.axhline(y=pipe_cfactor, color="b", linestyle="--", alpha=0.5)

    # now make the plot nice
    if xaxisval == "timemid":
        ax.set_xlabel(f"Time [MJD] - {startday}")
    elif xaxisval == "rate":
        ax.set_xscale("log")
        ax.set_xlabel("Central Pixel Rate [DN/s]")
    elif xaxisval == "welldepth":
        ax.set_xlabel("Central Pixel Well Depth [DN]")
    elif xaxisval == "bkg":
        ax.set_xscale("log")
        ax.set_xlabel("Background [DN/s]")
    else:
        ax.set_xscale("log")
        ax.set_xlabel("Flux [mJy]")
    ax.set_ylabel("Calibration Factors [(MJy/sr) / (DN/s)]")
    ax.set_title(f"{filter} / EEFRAC {args.eefrac}")

    def val2per(val):
        return (val / medval) * 100.0 - 100.0

    def per2val(per):
        return ((per + 100) / 100.0) * medval

    secax = ax.secondary_yaxis("right", functions=(val2per, per2val))
    secax.set_ylabel("percentage")

    if showleg:

        # make space for the legend
        ylim = ax.get_ylim()
        ax.set_ylim(ylim[0], ylim[1] + 0.4 * (ylim[1] - ylim[0]))

        first_legend = [
            Patch(facecolor=ccol, edgecolor=ccol, label=cdir, alpha=0.5)
            for cdir, ccol in zip(dirs, pcols)
        ]
        first_legend.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label="Not in average or fit",
                markerfacecolor="none",
                markeredgecolor="m",
                markersize=13,
                alpha=0.5,
            )
        )
        leg1 = ax.legend(handles=first_legend, loc="upper center")
        ax.add_artist(leg1)

        second_legend = []
        for ckey in psubsym.keys():
            if ckey[0:4] != "MASK":
                second_legend.append(
                    Line2D(
                        [0],
                        [0],
                        marker=psubsym[ckey],
                        color="w",
                        label=ckey,
                        markerfacecolor="k",
                        markersize=10,
                        alpha=0.5,
                    )
                )
        ax.legend(handles=second_legend, loc="upper left")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter",
        help="filter to process",
        default="F770W",
        # fmt: off
        choices=["F560W", "F770W", "F1000W", "F1130W", "F1280W",
                 "F1500W", "F1800W", "F2100W", "F2550W",
                 "F1065C", "F1140C", "F1550C", "F2300C",
                 "F770W_subtrans", "F770W_repeat"]
        # fmt: on
    )
    parser.add_argument(
        "--bkgsub", help="use results from background subtraction run", action="store_true"
    )
    parser.add_argument(
        "--indivmos",
        help="use results from individual mosaics (1 per cal image) instead of combined mosaics",
        action="store_true",
    )
    parser.add_argument(
        "--indivcals",
        help="use results from individual cal images instead of combined mosaics",
        action="store_true",
    )
    parser.add_argument(
        "--eefrac", default=0.7, help="Enclosed energy fraction to use", type=float,
    )
    parser.add_argument(
        "--xaxisval",
        help="x-axis values",
        default="mflux",
        choices=["mflux", "timemid", "rate", "welldepth", "bkg"],
    )
    parser.add_argument(
        "--nosubarrcor",
        help="do not apply subarray correction factors",
        action="store_true",
    )
    parser.add_argument(
        "--nocurval", help="do not plot the current calfactor", action="store_true",
    )
    parser.add_argument(
        "--repeat", help="plot the repeatability observations", action="store_true",
    )
    parser.add_argument(
        "--subtrans", help="plot the subarray transfer observations", action="store_true",
    )
    parser.add_argument("--multiplot", help="4 panel plot", action="store_true")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # make plot
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    if args.bkgsub:
        extstr = "_bkgsub"
    else:
        extstr = ""
    if args.indivmos:
        extstr = "_indivmos"
    if args.indivcals:
        extstr = "_indivcals"
    if args.repeat:
        extstr = f"{extstr}_repeat"


    savefacs = f"CalFacs/miri_calfactors{extstr}_{args.filter}.fits"
    if args.multiplot:
        fontsize = 10
        font = {"size": fontsize}
        plt.rc("font", **font)

        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
        plot_calfactors(
            ax[0, 0],
            args.filter,
            "mflux",
            showleg=True,
            applysubarrcor=(not args.nosubarrcor),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans,           
        )
        plot_calfactors(
            ax[0, 1],
            args.filter,
            "timemid",
            savefile=savefacs,
            showleg=False,
            applysubarrcor=(not args.nosubarrcor),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans, 
        )
        plot_calfactors(
            ax[1, 0],
            args.filter,
            "welldepth",
            showleg=False,
            applysubarrcor=(not args.nosubarrcor),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans, 
        )
        plot_calfactors(
            ax[1, 1],
            args.filter,
            "bkg",
            showleg=False,
            applysubarrcor=(not args.nosubarrcor),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans, 
        )
        fname = f"miri_calfactors_{args.filter}_many"
    else:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
        plot_calfactors(
            ax,
            args.filter,
            args.xaxisval,
            savefile=savefacs,
            applysubarrcor=(not args.nosubarrcor),
            showcurval=(not args.nocurval),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans, 
        )
        fname = f"miri_calfactors_{args.filter}_{args.xaxisval}"

    plt.tight_layout()

    fname = f"{fname}_eefrac{args.eefrac}"
    if args.bkgsub:
        fname = f"{fname}_bkgsub"
    if args.indivmos:
        fname = f"{fname}_indivmos"
    if args.indivcals:
        fname = f"{fname}_indivcals"
    if args.repeat:
        fname = f"{fname}_repeat"
    if args.subtrans:
        fname = f"{fname}_subtrans"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
