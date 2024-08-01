import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from astropy.table import QTable


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

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(14, 7), width_ratios=[9, 4])

    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
               "F1500W", "F1800W", "F2100W", "F2550W",
               "F1065C", "F1140C", "F1550C", "F2300C"]
    fvals = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 0.0, 1.0, 2.0, 3.0]
    # syms = ["v", "^", ">", "<", "p", "X", "D", "P", "d", "o", "o", "o", "o"]

    files = ["CalFactors/jwst_miri_photom_0079.fits",
             "CalFactors/jwst_miri_photom_0201.fits",
             "Photom/jwst_miri_photom_flight_31jul24.fits"]
    dates = ["2022-06", "2023-09", "2024-07"]
    cols = ["b", "c", "g"]

    xtick_num = []
    xtick_txt = []
    delt = 0.25
    for ll, cfile in enumerate(files):
        cftab = QTable.read(cfile, hdu=1)
        if ll > 0:
            cftab_time = QTable.read(cfile, hdu=2)

        for k, cfilter in enumerate(filters):
            pipe_cfactor = cftab["photmjsr"][cftab["filter"] == cfilter][0]
            pipe_unc = cftab["uncertainty"][cftab["filter"] == cfilter][0]

            if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
                ax = axs[1]
            else:
                ax = axs[0]
            ax.errorbar(fvals[k] + ll * delt, pipe_cfactor, yerr=pipe_unc, fmt=f"{cols[ll]}o")

            if ll > 0:
                amp = cftab_time["amplitude"][cftab["filter"] == cfilter][0]
                ax.plot([fvals[k] + ll * delt], [pipe_cfactor - amp], f"{cols[ll]}o", fillstyle="none")
                ax.plot(np.array([1.0, 1.0]) * (fvals[k] + ll * delt), [pipe_cfactor, pipe_cfactor - amp],
                        f"{cols[ll]}-", alpha=0.5)

            if ll == 0:
                ax.plot(cfilter, [10.0])
                ax.plot(np.array([1.0, 1.0])*(fvals[k] - delt), [0.0, 5.0], "k:", alpha=0.7)

    axs[0].tick_params(axis='x', labelrotation=60)
    axs[1].tick_params(axis='x', labelrotation=60)

    axs[0].set_ylabel("C")
    axs[1].set_ylabel("C")
    axs[0].set_ylim(0.0, 1.6)
    axs[1].set_ylim(0.0, 5.0)
    axs[1].yaxis.tick_right()
    axs[1].yaxis.set_label_position("right")

    axs[0].text(0.0, 1.35, "Imaging", bbox=dict(facecolor='white', edgecolor='none'))
    axs[1].text(0.0, 4.5, "Coronagraphy", bbox=dict(facecolor='white', edgecolor='none'))

    first_legend = [
        Patch(facecolor=ccol, edgecolor=ccol, label=cdir, alpha=0.5)
        for cdir, ccol in zip(dates, cols)
    ]
    leg1 = axs[0].legend(handles=first_legend, loc="upper right", fontsize=0.8*fontsize)
    axs[0].add_artist(leg1)

    #ax.legend(ncol=3, fontsize=0.7*fontsize)

    plt.tight_layout()

    fname = f"delivered_photom_vs_time"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
