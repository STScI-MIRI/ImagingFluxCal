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

    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(14, 7), width_ratios=[9, 1, 4])

    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
               "F1500W", "F1800W", "F2100W", "F2550W", "FND",
               "F1065C", "F1140C", "F1550C", "F2300C"]
    fvals = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 0.0, 0.0, 1.0, 2.0, 3.0]
    # syms = ["v", "^", ">", "<", "p", "X", "D", "P", "d", "o", "o", "o", "o"]

    files = ["CalFactors/jwst_miri_photom_0079.fits",
             "CalFactors/jwst_miri_photom_0201.fits",
             "CalFactors/jwst_miri_photom_0218.fits",
             "Photom/jwst_miri_photom_flight_03sep25.fits",
            ]
    dates = ["2022-06", "2023-09", "2024-08", "2025-09"]
    cols = ["b", "c", "g", "r"]

    xtick_num = []
    xtick_txt = []
    delt = 0.2
    for ll, cfile in enumerate(files):
        cftab = QTable.read(cfile, hdu=1)
        if ll > 0:
            cftab_time = QTable.read(cfile, hdu=2)
            if ll > 2:
                cftab_time2 = QTable.read(cfile, hdu=3)

        for k, cfilter in enumerate(filters):
            if ((cfilter == "FND") & (ll < 2)):
                doband = False
            else:
                doband = True

            if doband:

                pipe_cfactor = cftab["photmjsr"][cftab["filter"] == cfilter][0]
                pipe_unc = cftab["uncertainty"][cftab["filter"] == cfilter][0]

                if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
                    ax = axs[2]
                elif cfilter == "FND":
                    ax = axs[1]
                else:
                    ax = axs[0]

                if ll > 0:
                    if ll > 2:
                        amp = cftab_time2["amplitude"][cftab["filter"] == cfilter][0]
                        const = cftab_time2["const"][cftab["filter"] == cfilter][0]
                        tau = cftab_time2["tau"][cftab["filter"] == cfilter][0]
                        t0 = cftab_time2["t0"][cftab["filter"] == cfilter][0]
                        lossperyear = cftab_time["lossperyear"][cftab["filter"] == cfilter][0]
                        t = 1200.0  # days
                        timefac = (1 - lossperyear*t/365.) * (amp * np.exp(-1. * t / tau) + const)
                        # print(cfilter, timefac)
                        startval = pipe_cfactor / (amp + const)
                        endval = pipe_cfactor / timefac
                    else:
                        amp = cftab_time["amplitude"][cftab["filter"] == cfilter][0]
                        startval = pipe_cfactor + amp
                        endval = pipe_cfactor
                    ax.plot([fvals[k] + ll * delt], [startval], f"{cols[ll]}o")
                    ax.errorbar(fvals[k] + ll * delt, endval, yerr=pipe_unc, fmt=f"{cols[ll]}o", fillstyle="none")

                    ax.plot(np.array([1.0, 1.0]) * (fvals[k] + ll * delt), [endval, startval],
                            f"{cols[ll]}-", alpha=0.5)
                else:
                    ax.errorbar(fvals[k] + ll * delt, pipe_cfactor, yerr=pipe_unc, fmt=f"{cols[ll]}o")

                if ll == 2:
                    ax.plot(cfilter, [100.0])
                    ax.plot(np.array([1.0, 1.0])*(fvals[k] - delt), [0.0, 100.0], "k:", alpha=0.7)

    axs[0].tick_params(axis='x', labelrotation=60)
    axs[1].tick_params(axis='x', labelrotation=60)
    axs[2].tick_params(axis='x', labelrotation=60)

    axs[1].set_xlim(-0.3, 1.0)

    axs[0].set_ylabel("C [(MJy/sr) / (DN/s/pix)]")
    axs[2].set_ylabel("C [(MJy/sr) / (DN/s/pix)]")
    axs[0].set_ylim(0.0, 1.6)
    axs[1].set_ylim(0.0, 52.0)
    axs[2].set_ylim(0.0, 5.0)
    axs[2].yaxis.tick_right()
    axs[2].yaxis.set_label_position("right")

    axs[0].text(0.0, 1.35, "Imaging", bbox=dict(facecolor='white', edgecolor='none'))
    axs[2].text(0.0, 4.5, "Coronagraphy", bbox=dict(facecolor='white', edgecolor='none'))

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
