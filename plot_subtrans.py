import argparse
import numpy as np
import os

import matplotlib.pyplot as plt

from astropy.table import QTable

from calc_calfactors import get_calfactors

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter",
        help="filter to process",
        default="F770W",
        choices=["F770W", "F1280W"]
    )

    # HDL: I'm not sure this is actually the best way to designate data
    #      could be worth it to have some kind of control file with a list of PIDs and Obs
    #      to use in the analysis? Would be kind of a big overhaul. Think more about this.
    parser.add_argument(
        "--subtrans_pid",
        help="select the subarray transfer calibrator by pid(s)",
        default=None,
        nargs="+",
        type=int,
    )

    parser.add_argument(
        "--startday",
        help="start day for the observations",
        default=59720.,
        type=float,
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

    #startday = 59720.

    # HDL: not sure if this makes sense or meaningfully catches any logical errors
    assert not (args.subtrans_pid is None and args.startday is None), "Give either PID or startday"

    cfacs = get_calfactors(
        "ADwarfs",
        args.filter,
        xaxisval="timemid",
        subtrans=True,
        subtrans_pid=args.subtrans_pid,
        startday=args.startday,
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

    # compute the factors compared to FULL
    aves = []
    subarrs = np.array(["FULL", "BRIGHTSKY", "SUB256", "SUB128", "SUB128_IP", "SUB64", "SUB64_IP"])
    for csub in subarrs:
        fvals = [True if tsub == csub else False for tsub in cfacs[3]]
        aves.append(np.average(yvals[fvals]))
    yvals /= aves[0]
    aves = aves / aves[0]

    sindxs = np.argsort(cfacs[2])
    mstd = np.std(yvals) * 100.

    ax.errorbar(xvals, yvals[sindxs], yerr=yvals_unc, fmt="ko")
    #if args.filter == "F1280W":
    #    new_sub128 = np.average([aves[2], aves[4]])
    #    aves[3] = new_sub128
    #    ax.text(xvals[4], new_sub128 + 0.003, "Adopted", ha="center")
    #    ax.errorbar(xvals[4], new_sub128, fmt="ko", alpha=0.5)

        # based on visualizing all bands
    #    new_bs = 1.015
    #    aves[1] = new_bs
    #    ax.text(xvals[2], new_bs + 0.003, "Adopted", ha="center")
    #    ax.errorbar(xvals[2], new_bs, fmt="ko", alpha=0.5)

    ax.set_xticks(xvals)
    ax.set_xticklabels(np.array(cfacs[3])[sindxs])

    ax.plot([1.0, len(yvals)], [1.0, 1.0], "k:", alpha=0.5)
    # ax.text(3, 0.9925, rf"$\sigma$ = {mstd:.2f}%")

    for cave, csub in zip(aves, subarrs):
        print(f"{csub} & {cave:.3f} \\\\")
    otab = QTable()
    otab["name"] = subarrs
    otab["FracChange"] = aves

    
    if not os.path.isdir("./CalFacs/"): os.makedirs("./CalFacs/")
    otab.write(f"CalFacs/subarray_transfer_{args.filter}.dat",
               overwrite=True,
               format="ascii.commented_header")

    if args.filter == "F770W":
        ax.set_ylim(0.99, 1.01)
    ax.set_ylabel("Fractional change")

    plt.tight_layout()

    fname = f"subtrans_{args.filter}"
    if args.png:
        if not os.path.isdir("./Figs/"): os.makedirs("./Figs/")
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        if not os.path.isdir("./Figs/"): os.makedirs("./Figs/")
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
