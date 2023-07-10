import argparse
import matplotlib.pyplot as plt
import numpy as np

from astropy.table import QTable
import astropy.units as u

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bkgsub", help="compute and subtract background image", action="store_true"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    eefraction = 0.7
    if args.bkgsub:
        extstr = "_bkgsub"
    else:
        extstr = ""
    extstr = "_repeat"
    dir = "ADwarfs"
    sname = "BD+60 1753"

    # make plot
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))

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
    flux1 = []
    flux1_unc = []
    flux2 = []
    flux2_unc = []
    bkg1 = []
    bkg2 = []
    for cfilter in filters:
        obstab = QTable.read(f"{dir}/{cfilter}{extstr}_eefrac{eefraction}_phot.fits")

        (mindx,) = np.where((obstab["name"] == sname) & (obstab["timemid"] < 59800. * u.day))
        # get the earliest one
        sindxs = np.argsort(obstab["timemid"][mindx])
        mindx = mindx[sindxs]
        print(cfilter, obstab["timemid"][mindx[0]], obstab["subarray"][mindx[0]])
        flux1.append(obstab["aperture_sum_bkgsub"][mindx[0]].value)
        flux1_unc.append(obstab["aperture_sum_bkgsub_err"][mindx[0]].value)
        bkg1.append(obstab["mean_bkg"][mindx[0]].value)
        # print(cfilter, obstab["timemid"][mindx[0]])

        (mindx,) = np.where((obstab["name"] == sname) & (obstab["timemid"] > 60000. * u.day))
        print(cfilter, obstab["timemid"][mindx[0]], obstab["subarray"][mindx[0]])
        flux2.append(obstab["aperture_sum_bkgsub"][mindx[0]].value)
        flux2_unc.append(obstab["aperture_sum_bkgsub_err"][mindx[0]].value)
        bkg2.append(obstab["mean_bkg"][mindx[0]].value)
        # print(cfilter, obstab["timemid"][mindx[0]])

    flux1 = np.array(flux1)
    flux1_unc = np.array(flux1_unc)
    flux2 = np.array(flux2)
    flux2_unc = np.array(flux2_unc)
    bkg1 = np.array(bkg1)
    bkg2 = np.array(bkg2)

    ratio = flux2 / flux1
    ratio_unc = np.square(flux1_unc / flux1) + np.square(flux2_unc / flux2)
    ratio_unc = ratio * np.sqrt(ratio_unc)

    print(filters)
    print(ratio)
    print(ratio_unc)

    ax.set_xlabel(r"$\lambda$ [$\mu$m]")
    ax.set_ylabel("flux (MJD=60131 d)/flux (MJD=59724 d/59762 d) ")

    ax.plot([7.0, 26.0], [1., 1.], "k:", alpha=0.5)

    ax.errorbar(waves, ratio, fmt="bo", yerr=ratio_unc, label="flux")
    # ax.errorbar(waves, bkg2 / bkg1, fmt="ko", label="bkg", alpha=0.5)

    # ax.set_ylim(0.9, 1.05)

    # ax.legend(fontsize=0.8 * fontsize)

    ax.set_title("BD+60 1753 repeat test")

    plt.tight_layout()

    fname = "hd+60d1753_repeat"
    if args.bkgsub:
        fname = f"{fname}_bkgsub"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
