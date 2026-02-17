import argparse
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from astropy.table import QTable
from statsmodels.stats.weightstats import DescrStatsW
from astropy.stats import sigma_clip


def compute_stats(allfacs, weights, sigcut):
    "compute the weighted mean, weighted std, and weighted std of the mean"

    if len(allfacs) > 0:
        filtered_data = sigma_clip(allfacs, sigma=sigcut, maxiters=5)

        # compute the weighted mean
        newvals = allfacs[~filtered_data.mask]
        newweights = weights[~filtered_data.mask]
        meanval = np.average(newvals, weights=newweights)
        # compute weighted standard deviation
        if len(allfacs) > 3:
            meanstd = DescrStatsW(newvals, weights=newweights, ddof=1).std
            meanstdmean = meanstd / np.sqrt(np.sum(newweights))
        else:
            meanstd = 0.0
            meanstdmean = 0.0
        return (meanval, meanstd, meanstdmean, filtered_data, np.sum(newweights))
    else:
        return (None, None, None, None, None)


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
                 "FND"]
        # fmt: on
    )
    parser.add_argument(
        "--indivcals",
        help="use results from individual cal images instead of the individual mosaics",
        action="store_true",
    )
    parser.add_argument("--bkg", help="include background measurements", action="store_true")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    pcols = {"HotStars": "b",
             "ADwarfs": "g", 
             "SolarAnalogs": "r",
             "Asteroids": "k"}

    psubsym = {
        "FULL": "o",
        "BRIGHTSKY": "s",
        "SUB256": "p",
        "SUB128": "P",
        "SUB64": "*",
    }

    if args.indivcals:
        stype = "cals"
    else:
        stype = "mos"

    # make plot
    fontsize = 22
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 8))

    mu, sigma = 0, 0.12  # mean and standard deviation
    rng = np.random.default_rng()

    all_ratios = []
    all_uncs = []
    all_dithers = []
    all_dithers_rng = []
    for cdir in pcols.keys():
        filter = args.filter
        obstab = QTable.read(f"{cdir}/{filter}_indiv{stype}_eefrac0.7_phot.fits")

        subfnames = np.array([(cfile.split("/")[3].split("_0000"))[0] for cfile in obstab["filename"]])
        usub = np.unique(subfnames)
        for csub in usub:
            gvals = np.array([True if csub == tsub else False for tsub in subfnames])
            dithvals = np.array([int(cdith.decode()) for cdith in (obstab["exposure"].value)[gvals]])
            subarrs = obstab["subarray"][gvals]
            obsflux = obstab["aperture_sum_bkgsub"][gvals]
            obsfluxunc = obstab["aperture_sum_bkgsub_err"][gvals]
            aveobs = np.average(obsflux)
            xvals = dithvals + rng.normal(mu, sigma)
            ax.errorbar(xvals, obsflux / aveobs, yerr=obsfluxunc / aveobs,
                        fmt=f"{pcols[cdir]}{psubsym[subarrs[0]]}",
                        alpha=0.5,
                        markersize=10)

            if args.bkg:
                obsbkg = obstab["mean_bkg"][gvals]
                avebkg = np.average(obsbkg)
                ax.errorbar(xvals, 0.05 + obsbkg / avebkg,
                            fmt=f"k{psubsym[subarrs[0]]}",
                            alpha=0.5,
                            markersize=10)

            all_ratios = np.concatenate([all_ratios, obsflux / aveobs])
            all_uncs = np.concatenate([all_uncs, obsfluxunc / aveobs])
            all_dithers = np.concatenate([all_dithers, dithvals])
            all_dithers_rng = np.concatenate([all_dithers_rng, xvals])

    all_dithers = np.array(all_dithers)

    # compute the averages by dither position
    for k in range(4):
        gvals = ((k + 1) == all_dithers) & (np.isfinite(all_ratios))
        res = compute_stats(all_ratios[gvals], 1.0/np.square(all_uncs[gvals]), 3.5)
        # aveval = np.average(all_ratios[gvals], weights=1.0/np.square(all_uncs[gvals]))
        ax.errorbar([k+1], [res[0]], yerr=[res[1]], fmt="k*", markersize=20)
        ax.text(k+1, 0.97, rf"{res[0].value:.4f} $\pm$ {res[1]:.4f}", fontsize=0.7*fontsize, alpha=0.75,
                ha="center")
        print(k+1, np.absolute(res[0] - 1.0))

        filtered_data = res[3]
        ax.scatter(
            (all_dithers_rng[gvals])[filtered_data.mask],
            (all_ratios[gvals])[filtered_data.mask],
            s=200,
            facecolor="none",
            edgecolor="m",
        )

    ax.axhline(1.0, color="k", linestyle="-", alpha=0.5)

    ax.set_title(filter)
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xticklabels(["1", "2", "3", "4"])
    ax.set_xlabel("Dither #")
    ax.set_ylabel("flux/(ave flux)")
    if args.bkg:
        ax.set_ylim(0.965, 1.075)
    else:
        ax.set_ylim(0.965, 1.025)

    # make space for the legend
    ylim = ax.get_ylim()
    ax.set_ylim(ylim[0], ylim[1] + 0.4 * (ylim[1] - ylim[0]))

    first_legend = [
        Patch(facecolor=ccol, edgecolor=ccol, label=cdir, alpha=0.5)
        for cdir, ccol in zip(pcols.keys(), pcols.values())
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
    loc = "upper right"
    leg1 = ax.legend(handles=first_legend, fontsize=0.7*fontsize, loc=loc)
    ax.add_artist(leg1)

    second_legend = []
    for ckey in psubsym.keys():
        if ckey[0:4] != "xMASK":
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
    ncol = 1
    ax.legend(handles=second_legend, fontsize=0.7*fontsize, ncol=ncol, loc="upper left")


    plt.tight_layout()

    fname = f"ditherpos_{filter}"
    if args.bkg:
        fname = f"{fname}_wbkg"
    if args.indivcals:
        fname = f"{fname}_cals"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()