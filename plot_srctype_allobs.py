import argparse
import numpy as np

import matplotlib.pyplot as plt
from statsmodels.stats.weightstats import DescrStatsW

from astropy.table import QTable

from calc_calfactors import get_calfactors

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # make plot
    fontsize = 20
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=3)
    plt.rc("axes", linewidth=3)
    plt.rc("xtick.major", width=3)
    plt.rc("ytick.major", width=3)

    cftab = QTable.read("Photom/jwst_miri_photom_flight_31jul24.fits", hdu=1)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 7))

    delt = 0.05

    dirs = {"HotStars": 0.0,
            "ADwarfs": 1.0,
            "SolarAnalogs": 2.0}

    ax.plot(list(dirs.keys()), [10.0, 10.0, 10.0], "*", markersize=20)

    # add in measurements from each of the bands
    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
               "F1500W", "F1800W", "F2100W", "F2550W"]
    syms = ["v", "^", ">", "<", "p", "X", "D", "P", "d"]
    subarrs_vals = np.array((list(dirs.values())))
    arelvals = np.zeros((9, 3))
    arelvals_unc = np.zeros((9, 3))
    for k, cfilter in enumerate(filters):
        pipe_cfactor = cftab["photmjsr"][cftab["filter"] == cfilter][0]
        pipe_unc = cftab["uncertainty"][cftab["filter"] == cfilter][0]

        cfname = f"CalFacs/miri_calfactors_grieke_subarracor_timecor_{cfilter}_srctype.dat"
        # cfname = f"CalFacs/miri_calfactors_subarracor_timecor_{cfilter}_srctype.dat"
        tab = QTable.read(cfname,
                          format="ascii.commented_header")
        # relative to SUB256 as it is the one always present
        relvals = (pipe_cfactor / tab["calfacs"])
        relvals_unc = np.square(tab["calfacs_uncmean"] / tab["calfacs"])
        gvals = relvals_unc > 0.0
        relvals_unc += np.square(pipe_unc / pipe_cfactor)
        relvals_unc = relvals * np.sqrt(relvals_unc)
        ax.errorbar(subarrs_vals[gvals] + (k+3)*delt, relvals[gvals], 
                    yerr=relvals_unc[gvals],
                    fmt=syms[k], label=f"{cfilter} $\pm 1\sigma$", alpha=0.5)

        arelvals[k, :] = relvals
        arelvals_unc[k, gvals] = relvals_unc[gvals]

    # determine the mean ratios
    ameans = []
    auncs = []
    for i, csrc in enumerate(list(dirs.keys())):
        gvals = arelvals_unc[:, i] > 0
        weights = 1.0 / np.square(arelvals_unc[gvals, i])
        meanval = np.average(arelvals[gvals, i], weights=weights)
        # compute weighted standard deviation
        if np.sum(gvals) > 3:
            meanstd = DescrStatsW(arelvals[gvals, i], weights=weights, ddof=1).std
            meanstdmean = meanstd / np.sqrt(np.sum(gvals))
        else:
            meanstd = 0.0
            meanstdmean = 0.0
        print(csrc, meanval, meanstd, meanstdmean)
        ameans.append(meanval)
        auncs.append(meanstdmean)

    ax.errorbar(list(dirs.keys()), ameans, yerr=np.array(auncs) * 3, fmt="k*", markersize=20, label=r"Average $\pm 3\sigma$")

    ax.plot([-1.0, 3.5], [1.0, 1.0], "k:")

    ax.set_ylabel("C(srctype) / C")
    ax.set_xlim(-0.25, 2.75)
    ax.set_ylim(0.97, 1.03)

    ax.legend(ncol=3, fontsize=0.7*fontsize)

    plt.tight_layout()

    fname = f"srctype_allobs"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
