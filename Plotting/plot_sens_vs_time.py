import argparse
import matplotlib.pyplot as plt
import numpy as np

from astropy.table import QTable
import astropy.units as u

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

    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(14, 10), sharex=True, sharey=True)

    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
               "F1500W", "F1800W", "F2100W", "F2550W"]

    stype = {"BD+60 1753": "go",
             "HD 2811": "bo"}
    mfctype = [None, "none"]
    bkgtype = ["annulus", "image+annulus"]
    for k, cfilter in enumerate(filters):
        for j, ptype in enumerate(["_", "_indivcals_"]):
            ptab = QTable.read(f"ADwarfs/{cfilter}{ptype}eefrac0.7_phot.fits")
            for cname in ["BD+60 1753", "HD 2811"]:
                gvals = ptab["name"] == cname
                if np.sum(gvals) > 0:
                    mtime = ptab["timemid"][gvals] - 59700. * u.day
                    oflux = ptab["aperture_sum_bkgsub"][gvals]
                    print(cname, cfilter)
                    print(oflux)
                    oflux_unc = ptab["aperture_sum_bkgsub_err"][gvals]
                    normval = np.nanmean(oflux[mtime < 100. * u.day])
                    ax[k//3, k%3].errorbar(mtime, oflux / normval, yerr=oflux_unc / normval,
                                        fmt=stype[cname], label=f"{cname} / bkg {bkgtype[j]}", alpha=0.5, mfc=mfctype[j])        
        ax[k//3, k%3].plot([0, 500], [1.0, 1.0], "k--", alpha=0.5)
        ax[k//3, k%3].set_title(cfilter)
        ax[k//3, k%3].set_xlim(left=0.)

    for k in range(3):
        ax[k, 0].set_ylabel("Rel Flux")
        ax[2, k].set_xlabel("MJD-57900")

    ax[0, 1].legend(fontsize=0.7*fontsize)
    plt.tight_layout()

    fname = "miri_imaging_sens_vs_time"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()