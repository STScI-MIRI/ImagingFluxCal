import argparse
import matplotlib.pyplot as plt
import numpy as np

from astropy.table import QTable

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="name of star")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    eefraction = 0.8
    extstr = "_bkgsub"
    dir = "SolarAnalogs"

    # make plot
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 8))

    # read in the calfactors for each filter
    waves = np.array([5.6, 7.7, 10.0, 11.3, 12.8, 15.0, 18.0, 21.0, 25.5])
    filters = [
        "F560W",
        "F770W",
        "F1000W",
        "F1130W",
        "F1280W",
        "F1500W",
        "F1800W",
        "F2100W",
        "F2550W",
    ]
    atab = None
    flux = []
    flux_unc = []
    ofilters = []
    for cfilter in filters:
        ctab = QTable.read(
            f"CalFacs/miri_calfactors{extstr}_{cfilter}_ave.dat",
            format="ascii.commented_header",
        )
        obstab = QTable.read(f"{dir}/{cfilter}{extstr}_eefrac{eefraction}_phot.fits")
        (mindx,) = np.where(obstab["name"] == args.name)
        if len(mindx) > 0:
            ofilters.append(cfilter)
            apcor = obstab["apcorr"][mindx[0]]
            pixarea = obstab["pixarea"][mindx[0]]
            tcfac = ctab[f"avecalfac_{cfilter}"] * apcor * pixarea * 1e6
            flux.append(obstab["aperture_sum_bkgsub"][mindx[0]] * tcfac)
            flux_unc.append(obstab["aperture_sum_bkgsub_err"][mindx[0]] * tcfac)

    outtab = QTable()
    outtab["filter"] = ofilters
    outtab["flux"] = flux
    outtab["unc"] = flux_unc
    print(outtab)

    ax.set_xlabel(r"$\lambda$ [$\mu$m]")
    ax.set_ylabel("calfac / (median calfac)")

    ax.legend(fontsize=0.8 * fontsize, ncol=2)

    plt.tight_layout()

    fname = "Figs/miri_calfactors_allwaves"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
