import argparse
import numpy as np

import matplotlib.pyplot as plt

from astropy.modeling import models, fitting

from calc_calfactors import get_calfactors

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter",
        help="filter to process",
        default="F770W",
        choices=["F770W", "F1280W"]
    )
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

    startday = 59720.

    cfacs = get_calfactors(
        "ADwarfs",
        args.filter,
        xaxisval="timemid",
        subtrans=True,
        startday=startday,
        bkgsub=False,
    )

    yvals = cfacs[0]
    yvals_unc = cfacs[1]
    #if args.filter == "F1280W":
    xvals = np.arange(len(yvals)) + 1
    #else:
    #    xvals = [2, 5, 4, 1, 3]

    mval = np.average(yvals)
    yvals /= mval
    yvals_unc /= mval

    yvals = 1 / yvals

    sindxs = np.argsort(cfacs[2])

    mstd = np.std(yvals) * 100.

    ax.errorbar(xvals, yvals[sindxs], yerr=yvals_unc, fmt="ko")

    ax.set_xticks(xvals)
    ax.set_xticklabels(np.array(cfacs[3])[sindxs])

    ax.plot([1.0, len(yvals)], [1.0, 1.0], "k:", alpha=0.5)
    ax.text(3, 0.9925, rf"$\sigma$ = {mstd:.2f}%")

    if args.filter == "F770W":
        ax.set_ylim(0.99, 1.01)
    ax.set_ylabel("Fractional change")

    plt.tight_layout()

    # compute the factors compared to FULL
    aves = []
    for csub in ["FULL", "BRIGHTSKY", "SUB256", "SUB128", "SUB64"]:
        fvals = [True if tsub == csub else False for tsub in cfacs[3]]
        aves.append(np.average(yvals[fvals]))
    print(aves / aves[0])


    fname = f"subtrans_{args.filter}"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
