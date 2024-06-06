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

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 5))

    subarrs = np.array(["FULL", "BRIGHTSKY", "SUB256", "SUB128", "SUB64"])

    atab = QTable.read("CalFacs/subarray_transfer_F1280W.dat",
                      format="ascii.commented_header")
    ax.plot(atab["name"], atab["FracChange"], "ko", alpha=0.7, label="dedicated F1280W")

    # add in measurements from each of the bands
    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
               "F1500W", "F1800W", "F2100W", "F2550W"]
    for cfilter in filters:
        tab = QTable.read(f"CalFacs/miri_calfactors_timecor_{cfilter}_subarr.dat",
                          format="ascii.commented_header")
        # relative to SUB256 as it is the one always present
        relvals = atab["FracChange"][2] * tab["calfacs"][2] / tab["calfacs"]
        gvals = relvals > 0.0
        ax.scatter(subarrs[gvals], relvals[gvals], label=cfilter, alpha=0.5)

    ax.set_ylabel("Fractional change")
    ax.set_ylim(0.95, 1.05)

    ax.legend(ncol=2)

    plt.tight_layout()

    fname = f"subtrans_allobs"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
