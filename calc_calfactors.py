from os.path import exists
import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np

from astropy.table import QTable


def get_calfactors(dir, filter, xaxisval="mflux"):
    """
    Read in the observed and mdoel fluxes and computer the calibration factors
    """
    # read in observed fluxes
    obstab = QTable.read(f"{dir}/{filter}_phot.fits")
    # read in model fluxes
    modtab = QTable.read("Models/model_phot.fits")

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
        else:
            xval = mflux * 1e3

        cfactor = 1e-6 * mflux.value / (oflux.value * apcorr * pixarea.value)
        cfactor_unc = (oflux_unc / oflux) * cfactor
        xvals.append(xval.value)
        cfactors.append(cfactor)
        cfactors_unc.append(cfactor_unc)
        subarrs.append(obstab["subarray"][k])

    res = (cfactors, cfactors_unc, xvals, subarrs)

    return res


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
        "--xaxisval",
        help="x-axis values",
        default="mflux",
        choices=["mflux", "timemid"],
    )
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
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

    dirs = ["HotStars", "ADwarfs", "SolarAnalogs"]
    pcols = ["b", "g", "r"]
    filter = args.filter

    psubsym = {
        "FULL": "o",
        "BRIGHTSKY": "s",
        "SUB256": "^",
        "SUB128": "P",
        "SUB64": "*",
    }

    allfacs = []
    for k, dir in enumerate(dirs):
        if exists(f"{dir}/{filter}_phot.fits"):
            cfacs = get_calfactors(dir, filter, xaxisval=args.xaxisval)
            allfacs.append(cfacs[0])
            for cfactor, cfactor_unc, xval, subarray in zip(
                cfacs[0], cfacs[1], cfacs[2], cfacs[3]
            ):
                ax.errorbar(
                    [xval],
                    [cfactor],
                    yerr=[cfactor_unc],
                    fmt=f"{pcols[k]}{psubsym[subarray]}",
                    alpha=0.5,
                    markersize=10,
                )
    allfacs = np.concatenate(allfacs)
    medval = np.nanmedian(allfacs)

    # get the current pipeline calibration factor and plot
    cftab = QTable.read("CalFactors/jwst_miri_photom_0079.fits")
    pipe_cfactor = cftab["photmjsr"][cftab["filter"] == filter.split("_")[0]][0]
    ax.axhline(y=pipe_cfactor, color="b", linestyle="--", alpha=0.5)

    # now make the plot nice
    if args.xaxisval == "timemid":
        ax.set_xlabel("Time [MJD]")
    else:
        ax.set_xscale("log")
        ax.set_xlabel("Flux [mJy]")
    ax.set_ylabel("Calibration Factors [(MJy/sr) / (DN/s)]")
    ax.set_title(f"{filter} (fixed aperture, no aperture correction)")

    def val2per(val):
        return (val / medval) * 100.0 - 100.

    def per2val(per):
        return ((per + 100) / 100.0) * medval

    secax = ax.secondary_yaxis("right", functions=(val2per, per2val))
    secax.set_ylabel("percentage")

    first_legend = [
        Patch(facecolor=ccol, edgecolor=ccol, label=cdir, alpha=0.5)
        for cdir, ccol in zip(dirs, pcols)
    ]
    leg1 = ax.legend(handles=first_legend, loc="upper center")
    ax.add_artist(leg1)

    second_legend = []
    for ckey in psubsym.keys():
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

    plt.tight_layout()

    fname = f"miri_calfactors_{filter}"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
