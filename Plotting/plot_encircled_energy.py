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

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

    pixscale = 0.110

    files = {}
    dsets = {}

    if args.coron:
        for cfilter in ["F1065C", "F1140C"]:
            files[cfilter] = [
                f"ADwarfs/{cfilter}/del UMi_set1/miri_del UMi_set1_stage3_asn_i2d_ee.dat",
                f"ADwarfs/{cfilter}/HD 2811_set1/miri_HD 2811_set1_stage3_asn_i2d_ee.dat",
                f"SolarAnalogs/{cfilter}/HD 167060_set1/miri_HD 167060_set1_stage3_asn_i2d_ee.dat",
            ]
        # for copmarison
        # files["F1130W"] = ["ADwarfs/F1130W/BD+60 1753_set1/miri_BD+60 1753_set1_stage3_asn_i2d_ee.dat"]
    else:
        star = "BD+60 1753"
        if args.example:
            dsets["F770W"] = ["1"]
            dsets["F2100W"] = ["1"]
        else:
            dsets["F560W"] = ["1", "2", "3", "4"]
            dsets["F770W"] = ["1", "2", "3"]
            dsets["F1000W"] = ["1", "2", "3"]
            dsets["F1130W"] = ["1"]
            dsets["F1280W"] = ["1"]
            dsets["F1500W"] = ["1", "2", "3", "4"]
            dsets["F1800W"] = ["1", "2"]
            dsets["F2100W"] = ["1"]
            dsets["F2550W"] = ["1"]

            for cfilter in dsets.keys():
                files[cfilter] = []
                for cset in dsets[cfilter]:
                    files[cfilter].append(
                        f"ADwarfs/{cfilter}/{star}_set{cset}/miri_{star}_set{cset}_stage3_asn_i2d_ee.dat"
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
        "F2550C": ("tab:purple", "dashdot"),
    }

    for k, cfilter in enumerate(filters):
        for j, cfile in enumerate(files[cfilter]):
            ctab = QTable.read(cfile, format="ascii.commented_header")
            ccol, cline = pcol[cfilter]
            if j == 0:
                clabel = cfilter
            else:
                clabel = None
            ax.plot(
                ctab["radius"] * pixscale,
                ctab["ee"] + (len(dsets.keys()) - k) * 0.05,
                color=ccol,
                linestyle=cline,
                alpha=0.5,
                label=clabel,
            )

            if args.example:
                ax.plot(
                    ctab["radius"] * pixscale,
                    ctab["ee_model"] + (len(dsets.keys()) - k) * 0.05,
                    color="k",
                    linestyle="solid",
                    alpha=0.5,
                    label=f"{clabel} model",
                )

    ax.set_xlim(-1 * pixscale, 30 * pixscale)
    ax.set_xlabel("radius [arcsec]")
    ax.set_ylabel("Encircled Energy + const")

    ax.legend()

    plt.tight_layout()

    fname = "encirciled_energy"
    if args.example:
        fname = f"{fname}_example"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
