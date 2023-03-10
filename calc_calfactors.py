from os.path import exists
import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np

from astropy.table import QTable
from astropy.stats import sigma_clipped_stats


def get_calfactors(dir, filter, xaxisval="mflux", bkgsub=False, eefraction=0.8):
    """
    Read in the observed and mdoel fluxes and computer the calibration factors
    """
    if bkgsub:
        extstr = "_bkgsub"
    else:
        extstr = ""

    # read in observed fluxes
    obstab = QTable.read(f"{dir}/{filter}{extstr}_eefrac{eefraction}_phot.fits")
    # read in model fluxes
    modtab = QTable.read("Models/model_phot.fits")

    names = []
    xvals = []
    cfactors = np.zeros(len(obstab["name"]))
    cfactors_unc = np.zeros(len(obstab["name"]))
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

        cfactor = 1e-6 * mflux.value / (oflux.value * apcorr * pixarea.value)
        cfactor_unc = (oflux_unc / oflux) * cfactor
        # if obstab["name"][k] not in ["HD 167060", "16 Cyg B", "HD 37962", "del UMi"]:
        # print(obstab["name"][k], cfactor)
        names.append(obstab["name"][k])
        xvals.append(xval.value)
        cfactors[k] = cfactor
        cfactors_unc[k] = cfactor_unc
        subarrs.append(obstab["subarray"][k])

    res = (cfactors, cfactors_unc, xvals, subarrs, names)
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
    eefraction=0.8,
):
    """
    Plot the calibration factors versus the requested xaxis.
    """
    dirs = ["HotStars", "ADwarfs", "SolarAnalogs"]
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

    efac = 1.04
    # efac = 1.0
    # updated based on array-bkg subtraction reductions - better centroids (9 Mar 2023)
    subarr_cor = {
        "FULL": 1.0,
        "BRIGHTSKY": 0.93862651 * efac,
        "SUB256": 0.98265799 * efac,
        "SUB128": 0.96159878 * efac,
        "SUB64": 0.97713754 * efac,
        "MASK1065": 1.0,
        "MASK1140": 1.0,
        "MASK1550": 1.0,
        "MASKLYOT": 1.0,
    }

    # print(subarr_cor)
    ignore_names = ["HD 167060", "16 Cyg B", "HD 37962", "del UMi", "HD 180609"]

    allfacs = []
    allfacuncs = []
    allnames = []
    xvals = []
    for k, dir in enumerate(dirs):
        print(f"{dir}/{filter}_eefrac{eefraction}_phot.fits")
        if exists(f"{dir}/{filter}_eefrac{eefraction}_phot.fits"):
            cfacs = get_calfactors(
                dir, filter, xaxisval=xaxisval, bkgsub=bkgsub, eefraction=eefraction
            )
            # allfacs.append(cfacs[0])
            allfacuncs.append(cfacs[1])
            xvals.append(cfacs[2])
            allnames.append(cfacs[4])
            for cfactor, cfactor_unc, xval, subarray, cname in zip(
                cfacs[0], cfacs[1], cfacs[2], cfacs[3], cfacs[4]
            ):
                print(cname, xval, cfactor)
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
                if cname in ignore_names:
                    ax.scatter(
                        [xval], [cfactor], s=150, facecolor="none", edgecolor="m",
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
    for cname in allnames:
        if cname in ignore_names:
            gvals.append(False)
        else:
            gvals.append(True)

    meanvals = sigma_clipped_stats(allfacs[gvals], sigma=3, maxiters=5)
    ax.axhline(y=meanvals[0], color="k", linestyle="-", alpha=0.5)
    # ax.axhline(y=meanvals[0] + meanvals[2], color="k", linestyle=":", alpha=0.5)
    # ax.axhline(y=meanvals[0] - meanvals[2], color="k", linestyle=":", alpha=0.5)
    perstd = 100.0 * meanvals[2] / meanvals[0]
    # print(meanvals[0], meanvals[2], perstd)

    ax.text(
        0.7,
        0.02,
        fr"{meanvals[0]:.3f} $\pm$ {meanvals[2]:.3f} ({perstd:.1f}%)",
        transform=ax.transAxes,
    )

    if savefile is not None:
        atab = QTable()
        atab[f"avecalfac_{filter}"] = [meanvals[0]]
        atab[f"avecalfac_unc_{filter}"] = [meanvals[2] / np.sqrt(len(allfacs[gvals]))]
        atab[f"avecalfac_std_{filter}"] = [meanvals[2]]
        atab.write(
            savefile.replace(".fits", "_ave.dat"),
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
        ax.set_xlabel("Time [MJD]")
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

        # make space for th legend
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
                label="Not in average",
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
        choices=["F560W", "F770W", "F770W_subarray", "F770W_repeat", "F1000W",
                 "F1130W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W",
                 "F1065C", "F1140C", "F1550C", "F2300C"],
        # fmt: on
    )
    parser.add_argument(
        "--bkgsub", help="compute and subtract background image", action="store_true"
    )
    parser.add_argument(
        "--eefrac", default=0.8, help="Enclosed energy fraction to use", type=float,
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
            savefile=savefacs,
            applysubarrcor=(not args.nosubarrcor),
            bkgsub=args.bkgsub,
            eefraction=args.eefrac,
        )
        plot_calfactors(
            ax[0, 1],
            args.filter,
            "timemid",
            showleg=False,
            applysubarrcor=(not args.nosubarrcor),
            bkgsub=args.bkgsub,
            eefraction=args.eefrac,
        )
        plot_calfactors(
            ax[1, 0],
            args.filter,
            "welldepth",
            showleg=False,
            applysubarrcor=(not args.nosubarrcor),
            bkgsub=args.bkgsub,
            eefraction=args.eefrac,
        )
        plot_calfactors(
            ax[1, 1],
            args.filter,
            "bkg",
            showleg=False,
            applysubarrcor=(not args.nosubarrcor),
            bkgsub=args.bkgsub,
            eefraction=args.eefrac,
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
            eefraction=args.eefrac,
        )
        fname = f"miri_calfactors_{args.filter}_{args.xaxisval}"

    plt.tight_layout()

    fname = f"{fname}_eefrac{args.eefrac}"
    if args.bkgsub:
        fname = f"{fname}_bkgsub"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
