import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable, join

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

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
    for cfilter in filters:
        ctab = QTable.read(f"CalFacs/miri_calfactors_{cfilter}.fits")

        # merge multiple measurements of the same star
        gvals = np.full(len(ctab), False)
        for cname in np.unique(ctab["name"]):
            sindxs, = np.where(ctab["name"] == cname)
            if len(sindxs) > 1:
                ctab[f"calfac_{cfilter}"][sindxs[0]] = np.average(ctab[f"calfac_{cfilter}"][sindxs])
                ctab[f"calfac_{cfilter}_mean_dev"][sindxs[0]] = np.average(ctab[f"calfac_{cfilter}_mean_dev"][sindxs])
                ctab[f"calfac_{cfilter}_med_dev"][sindxs[0]] = np.average(ctab[f"calfac_{cfilter}_med_dev"][sindxs])
            gvals[sindxs[0]] = True
        ctab = ctab[gvals]

        if atab is not None:
            atab = join(atab, ctab, join_type="outer")
        else:
            atab = ctab
    nstars = len(atab)

    atab.write("miri_calfactors_all.fits", overwrite=True)

    # these two do not have F1500W obs, so nothing in the model column
    mindx, = np.where(atab["name"] == "GD 71")
    atab["modflux_F1500W"][mindx[0]] = 0.0
    mindx, = np.where(atab["name"] == "2MASS J17571324+6703409")
    atab["modflux_F1500W"][mindx[0]] = 0.0

    # sort by F1150W model flux as all are measured in this band
    aindxs = np.argsort(atab["modflux_F1500W"])

    atab = atab[aindxs]

    # print(atab)
    # atab.write("test.dat", format="ipac", overwrite=True)
    # exit()

    # make plot
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 8))

    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, nstars))))

    lstyle = ["solid", "dotted", "dashed", "dashdot"]

    for k, cname in enumerate(atab["name"]):
        pfacs = []
        pwaves = []
        for j, cfilter in enumerate(filters):
            pwaves.append(waves[j])
            pfacs.append(atab[k][f"calfac_{cfilter}_med_dev"])
        pfacs = np.array(pfacs)
        if np.sum(np.isfinite(pfacs)) > 1:
            ax.plot(pwaves, pfacs, linestyle=lstyle[k % 4], label=cname)
        else:
            ax.plot(pwaves, pfacs, "ko", label=cname)

    # ax.set_xlim(5.0, 35.0)
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
