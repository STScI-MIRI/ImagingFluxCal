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
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
               "F1500W", "F1800W", "F2100W", "F2550W",
               "F1065C", "F1140C", "F1550C", "F2300C"]

    # make plot
    fontsize = 12
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=7, ncols=2, figsize=(12, 18))

    startday = 59720.
    for k, cfilter in enumerate(filters):

        if cfilter in ["F2550W", "F1065C", "F1140C", "F1550C", "F2300C"]:
            bkgsub = True
            extstr = "_bkgsub"
        else:
            bkgsub = False
            extstr = ""

        applytime = True

        savefacs = f"CalFacs/miri_calfactors{extstr}_timecor_{cfilter}.fits"

        m = k + 1
        px = m // 2
        py = m % 2
        plot_calfactors(
            ax[px, py],
            cfilter,
            "mflux",
            # "inttime",
            # "welldepth",
            # "timemid",
            savefile=savefacs,
            showleg=False,
            showcurval=False,
            bkgsub=bkgsub, 
            applytime=applytime,
            grieke=args.grieke,
        )
        ax[px, py].set_ylabel("CalFactor")
        ax[px, py].set_title("")

        ax[px, py].text(0.1, 0.9, cfilter,
                        transform=ax[px, py].transAxes)

    ax[0, 1].set_ylim(0.4, 0.5)
    ax[2, 0].set_ylim(1.0, 1.2)
    ax[4, 1].set_ylim(0.65, 0.8)

    plt.tight_layout()

    fname = "multi_calfacs"
    if args.grieke:
        fname = f"{fname}_grieke"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()