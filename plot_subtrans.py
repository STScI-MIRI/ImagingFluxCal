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
    fontsize = 16
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5))

    startday = 59720.

    cfacs = get_calfactors(
        "ADwarfs",
        "F770W",
        xaxisval="timemid",
        subtrans=True,
        startday=startday,
        bkgsub=False,
    )

    yvals = cfacs[0]
    yvals_unc = cfacs[1]
    print(cfacs[3])
    xvals = [2, 5, 4, 1, 3]

    mval = np.average(yvals)
    yvals /= mval
    yvals_unc /= mval

    mstd = np.std(yvals) * 100.

    ax.errorbar(xvals, yvals, yerr=yvals_unc, fmt="ko")

    ax.set_xticks([1, 2, 3, 4, 5])
    ax.set_xticklabels(["FULL", "BRIGHTSKY", "SUB256", "SUB128", "SUB64"])

    ax.plot([1.0, 5.0], [1.0, 1.0], "k:", alpha=0.5)
    ax.text(3, 0.9925, rf"$\sigma$ = {mstd:.2f}%")

    ax.set_ylim(0.99, 1.01)
    ax.set_ylabel("Fractional change")

    plt.tight_layout()

    fname = "subtrans"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
