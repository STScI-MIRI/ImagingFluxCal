import argparse
import numpy as np

import matplotlib.pyplot as plt

from astropy.modeling import models, fitting

from calc_calfactors import plot_calfactors

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--grieke", help="use GRieke models for the 10 G-stars, CALSPEC models for the rest", action="store_true",
    )
    parser.add_argument(
        "--xaxisval",
        help="x-axis values",
        default="mflux",
        choices=["mflux", "timemid", "rate", "welldepth", "bkg", "inttime", "srctype", "subarr"],
    )
    parser.add_argument(
        "--subarrcor",
        help="Apply subarray correction factors",
        action="store_true",
    )
    parser.add_argument(
        "--noignore",
        help="Do not ignore any stars (useful for measuring the subarray differences)",
        action="store_true",
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
               "FND", "F1500W", "F1800W", "F2100W", "F2550W",
               "F1065C", "F1140C", "F1550C", "F2300C"]

    # make plot
    fontsize = 12
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=8, ncols=2, figsize=(14, 16), sharex=True)

    # setup the legend
    plot_calfactors(
        ax[0, 0],
        "F560W",
        args.xaxisval,
        showleg=True,
        showcurval=False, 
        x2ndaxis=False,
        notext=True,
        legonly=True,
        grieke=args.grieke,
        nofiltername=True,
        applytime=True,
    )
    ax[0,0].get_xaxis().set_visible(False)
    ax[0,0].get_yaxis().set_visible(False)
    ax[0,0].set_ylim(1.0, 2.)
    ax[0,0].set_title("")
    ax[0,0].axis('off')

    ax[0,1].get_xaxis().set_visible(False)
    ax[0,1].get_yaxis().set_visible(False)
    ax[0,1].set_ylim(1.0, 2.)
    ax[0,1].set_title("")
    ax[0,1].axis('off')

    startday = 59720.
    for k, cfilter in enumerate(filters):

        if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
            bkgsub = True
            extstr = "_bkgsub"
        else:
            bkgsub = False
            extstr = ""
        if args.grieke:
            extstr = f"{extstr}_grieke"
        if args.subarrcor:
            extstr = f"{extstr}_subarracor"

        applytime = True

        savefacs = f"CalFacs/miri_calfactors{extstr}_timecor_{cfilter}.fits"

        # if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
        #     noignore = False
        # else:
        #     noignore = False

        m = k + 2
        px = m // 2
        py = m % 2
        plot_calfactors(
            ax[px, py],
            cfilter,
            args.xaxisval,
            applysubarrcor=args.subarrcor,
            savefile=savefacs,
            showleg=False,
            showcurval=False,
            bkgsub=bkgsub, 
            applytime=applytime,
            grieke=args.grieke,
            noignore=args.noignore,
            fontsize=fontsize,
            nofiltername=True,
            # fitline=True,
        )
        ax[px, py].set_ylabel("C")
        ax[px, py].set_title("")
        if px < 7:
            ax[px, py].set_xlabel("")

        ax[px, py].text(0.1, 0.9, cfilter,
                        transform=ax[px, py].transAxes)

    #ax[0, 1].set_ylim(0.42, 0.48)
    #ax[2, 0].set_ylim(1.0, 1.2)
    #ax[4, 1].set_ylim(0.65, 0.8)

    if args.xaxisval == "welldepth":
        ax[0, 0].set_xlim(1e2, 65000.)
    elif args.xaxisval == "rate":
        ax[0, 0].set_xlim(1e1, 1e5)
    elif args.xaxisval == "inttime":
        ax[0, 0].set_xlim(0.1, 400)
    elif args.xaxisval == "subarr":
        ax[0, 0].set_xlim(-0.5, 8.5)

    plt.tight_layout()

    fname = f"multi_calfacs_{args.xaxisval}"
    if args.grieke:
        fname = f"{fname}_grieke"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()