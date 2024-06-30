import argparse
import numpy as np

import matplotlib.pyplot as plt

from astropy.table import QTable

from calc_calfactors import get_calfactors

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # make plot
    fontsize = 16
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

    # subarrs = np.array(["FULL", "BRIGHTSKY", "SUB256", "SUB128", "SUB64"])

    delt = 0.05
    subarrs = {"FULL": 0.0,
               "BRIGHTSKY": 1.0,
               "SUB256": 2.0,
               "SUB128": 3.0,
               "SUB64": 4.0}

    adopted_vals = [1.02, 1.03, 1.0, 1.01, 0.985]
    ax.plot(list(subarrs.keys()), adopted_vals, "*", markersize=20, label="Adopted")

    atab = QTable.read("CalFacs/subarray_transfer_F770W.dat",
                      format="ascii.commented_header")
    xvals = np.array([subarrs[csub] for csub in atab["name"].data])
    ax.plot(xvals + delt, atab["FracChange"] / atab["FracChange"][2], "ks", alpha=0.7, label="dedicated F770W")

    atab = QTable.read("CalFacs/subarray_transfer_F1280W.dat",
                      format="ascii.commented_header")
    xvals = np.array([subarrs[csub] for csub in atab["name"].data])
    ax.plot(xvals + 2 * delt, atab["FracChange"] / atab["FracChange"][2], "ko", alpha=0.7, label="dedicated F1280W")

    # add in measurements from each of the bands
    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
               "F1500W", "F1800W", "F2100W", "F2550W"]
    subarrs_vals = np.array((list(subarrs.values())))
    for k, cfilter in enumerate(filters):
        tab = QTable.read(f"CalFacs/miri_calfactors_timecor_{cfilter}_subarr.dat",
                          format="ascii.commented_header")
        # relative to SUB256 as it is the one always present
        relvals = (tab["calfacs"][2] / tab["calfacs"])  # * atab["FracChange"][2] 
        relvals_unc = relvals * (tab["calfacs_uncmean"] / tab["calfacs"])
        gvals = relvals > 0.0
        ax.errorbar(subarrs_vals[gvals] + (k+3)*delt, relvals[gvals], 
                    yerr=relvals_unc[gvals],
                    fmt="o", label=cfilter, alpha=0.5)
    
    ax.set_ylabel("SUB256 Fractional change")
    ax.set_ylim(0.95, 1.06)

    ax.legend(ncol=3, fontsize=0.7*fontsize)

    plt.tight_layout()

    fname = f"subtrans_allobs"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
