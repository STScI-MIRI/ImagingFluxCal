import argparse
import matplotlib.pyplot as plt

from astropy.table import QTable

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--example", help="show two filters with model", action="store_true"
    )
    parser.add_argument(
        "--coron", help="plot coronagraphic filters", action="store_true"
    )
    parser.add_argument(
        "--filter",
        help="filter to process",
        default="all",
        # fmt: off
        choices=["F560W", "F770W", "F1000W", "F1130W", "F1280W",
                 "F1500W", "F1800W", "F2100W", "F2550W",
                 "F1065C", "F1140C", "F1550C", "F2300C"]
        # fmt: on
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # make plot
    fontsize = 22
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

    pixscale = 0.110

    files = {}

    dsets = {"F560W": ["1", "2", "3", "4", "6"],
             "F770W": ["1", "2", "3", "4", "5", "6", "7", "8", "10", "11", "12", "13", "14", "15"],
             "F1000W": ["1", "2", "3", "4", "5"],
             "F1130W": ["1", "2", "3"],
             "F1280W": ["1", "2", "3"],
             "F1500W": ["2", "3", "4", "5", "6"],
             "F1800W": ["1", "2", "3"],
             "F2100W": ["1", "2", "3"],
             "F2550W": ["1", "2", "3"],
             }

    dfwhmfac = {"F560W": "25.0",
                "F770W": "25.0",
                "F1000W": "15.0",
                "F1130W": "10.0",
                "F1280W": "10.0",
                "F1500W": "10.0",
                "F1800W": "10.0",
                "F2100W": "5.0",
                "F2550W": "5.0",
                "F1065C": "15.0",
                "F1140C": "15.0",
                "F1550C": "10.0",
                "F2300C": "10.0",
               }

    if args.filter != "all":
        star = "BD+60 1753"
        cfilter = args.filter
        files[cfilter] = []
        if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
            files[cfilter] = [
                # f"ADwarfs/{cfilter}/BD+60 1753_set1/miri_BD+60 1753_set1_stage3_asn_i2d_ee_fwhmfac{dfwhmfac[cfilter]}.dat",
                f"ADwarfs/{cfilter}/del UMi_set1/miri_del UMi_set1_stage3_asn_i2d_ee_fwhmfac{dfwhmfac[cfilter]}.dat",
                f"ADwarfs/{cfilter}/HD 2811_set1/miri_HD 2811_set1_stage3_asn_i2d_ee_fwhmfac{dfwhmfac[cfilter]}.dat",
                f"SolarAnalogs/{cfilter}/HD 167060_set1/miri_HD 167060_set1_stage3_asn_i2d_ee_fwhmfac{dfwhmfac[cfilter]}.dat",
            ]
        else:
            for cset in dsets[cfilter]:
                files[cfilter].append(
                    f"ADwarfs/{cfilter}/{star}_set{cset}/miri_{star}_set{cset}_stage3_asn_i2d_ee_fwhmfac{dfwhmfac[cfilter]}.dat"
                )
    elif args.coron:
        for cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
            files[cfilter] = [
                f"ADwarfs/{cfilter}/del UMi_set1/miri_del UMi_set1_stage3_asn_i2d_ee_fwhmfac{dfwhmfac[cfilter]}.dat",
                f"ADwarfs/{cfilter}/HD 2811_set1/miri_HD 2811_set1_stage3_asn_i2d_ee_fwhmfac{dfwhmfac[cfilter]}.dat",
                f"SolarAnalogs/{cfilter}/HD 167060_set1/miri_HD 167060_set1_stage3_asn_i2d_ee_fwhmfac{dfwhmfac[cfilter]}.dat",
            ]
        # remove HD 2811 as the jwst 1.13.4 reductions look odd
        cfilter = "F1550C"
        files[cfilter] = [
            f"ADwarfs/{cfilter}/del UMi_set1/miri_del UMi_set1_stage3_asn_i2d_ee_fwhmfac{dfwhmfac[cfilter]}.dat",
            f"SolarAnalogs/{cfilter}/HD 167060_set1/miri_HD 167060_set1_stage3_asn_i2d_ee_fwhmfac{dfwhmfac[cfilter]}.dat",
        ]
    else:
        star = "BD+60 1753"
        if args.example:
            dsets = {}
            dsets["F770W"] = ["1"]
            dsets["F2100W"] = ["1"]

        for cfilter in dsets.keys():
            files[cfilter] = []
            for cset in dsets[cfilter]:
                files[cfilter].append(
                    f"ADwarfs/{cfilter}/{star}_set{cset}/miri_{star}_set{cset}_stage3_asn_i2d_ee_fwhmfac{dfwhmfac[cfilter]}.dat"
                )

    filters = list(files.keys())

    pcol = {
        "F560W": ("tab:blue", "solid"),
        "F770W": ("tab:orange", "dashed"),
        "F1000W": ("tab:green", "dotted"),
        "F1130W": ("tab:purple", "dashdot"),
        "F1280W": ("tab:brown", "solid"),
        "F1500W": ("tab:pink", "dashed"),
        "F1800W": ("tab:gray", "dotted"),
        "F2100W": ("tab:olive", "dashdot"),
        "F2550W": ("tab:cyan", "solid"),
        "F1065C": ("tab:blue", "solid"),
        "F1140C": ("tab:orange", "dashed"),
        "F1550C": ("tab:green", "dotted"),
        "F2300C": ("tab:purple", "dashdot"),
    }

    for k, cfilter in enumerate(filters):

        if args.filter == "all":
            if args.coron:
                offval = (3 - k) * 0.1
            else:
                offval = (8 - k) * 0.1
        else:
            offval = 0.0

        # get aperture corrections and uncertainties
        tab = QTable.read("ApCor/jwst_miri_apcorr_flight_2jul24_full.fits")
        gtab = tab[(tab["subarray"] == "FULL") & (tab["filter"] == cfilter)]

        ax.errorbar(gtab["radius"] * pixscale, gtab["eefraction"] + offval, xerr=gtab["radius_unc"] * pixscale, fmt="k.")
        # plot the sky outer
        ax.errorbar([gtab["skyout"][0] * pixscale], [0.85 + offval], xerr=[gtab["skyout_unc"][0] * pixscale], fmt="k.")

        for j, cfile in enumerate(files[cfilter]):
            print(cfilter, cfile)
            ctab = QTable.read(cfile, format="ascii.commented_header")
            aptab = QTable.read(cfile.replace("_ee_", "").replace(".dat", "_apcor.dat"), format="ascii.commented_header")

            ccol, cline = pcol[cfilter]
            if j == 0:
                clabel = cfilter
            else:
                clabel = None

            ax.plot(
                ctab["radius"] * pixscale,
                ctab["ee"] + offval,
                color=ccol,
                linestyle=cline,
                alpha=0.5,
                label=clabel,
            )

            if args.example:
                ax.plot(
                    ctab["radius"] * pixscale,
                    ctab["ee_model"] + offval,
                    color="k",
                    linestyle="solid",
                    alpha=0.5,
                    label=f"{clabel} model",
                )

    if args.coron:
        xmax_fac = 50
    else:
        xmax_fac = 30

    ax.set_ylim(bottom=0.0)
    ax.set_xlim(-1 * pixscale, xmax_fac * pixscale)
    ax.set_xlabel("radius [arcsec]")
    ax.set_ylabel("Encircled Energy + const")

    ax.legend(ncol=2)

    plt.tight_layout()

    fname = "encirciled_energy"
    if args.example:
        fname = f"{fname}_example"
    if args.coron:
        fname = f"{fname}_coron"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
