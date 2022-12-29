import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np

from astropy.table import QTable


def get_calfactors(dir, filter):
    """
    Read in the observed and mdoel fluxes and computer the calibration factors
    """
    # read in observed fluxes
    obstab = QTable.read(f"{dir}/{filter}_phot.fits")
    # read in model fluxes
    modtab = QTable.read("Models/model_phot.fits")

    mfluxes = []
    cfactors = []
    cfactors_unc = []
    subarrs = []
    for k, cname in enumerate(obstab["name"]):
        oflux = obstab["aperture_sum_bkgsub"][k]
        oflux_unc = obstab["aperture_sum_bkgsub_err"][k]

        (mindx,) = np.where(modtab["name"] == cname)
        if len(mindx) < 1:
            print(f"Model fluxes for {cname} not present")
            exit()
        mflux = modtab[filter][mindx[0]]

        cfactor = mflux.value / oflux.value
        cfactor_unc = (oflux_unc / oflux) * cfactor
        mfluxes.append(mflux.value)
        cfactors.append(cfactor)
        cfactors_unc.append(cfactor_unc)
        subarrs.append(obstab["subarray"][k])

    res = zip(cfactors, cfactors_unc, mfluxes, subarrs)

    return res


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
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
    filter = "F770W"

    psubsym = {
        "FULL": "o",
        "BRIGHTSKY": "s",
        "SUB256": "^",
        "SUB128": "P",
        "SUB64": "*",
    }

    for k, dir in enumerate(dirs):
        cfacs = get_calfactors(dir, filter)
        for cfactor, cfactor_unc, mflux, subarray in cfacs:

            ax.errorbar(
                [mflux],
                [cfactor],
                yerr=[cfactor_unc],
                fmt=f"{pcols[k]}{psubsym[subarray]}",
                alpha=0.5,
                markersize=10,
            )

    ax.set_xscale("log")
    ax.set_xlabel("Flux [Jy]")
    ax.set_ylabel("Calibration Factors [Jy / (DN/s)]")
    ax.set_title(f"{filter} (fixed aperture, no aperture correction)")

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
