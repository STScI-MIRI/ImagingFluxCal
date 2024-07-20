import argparse
from astropy.io import fits
import matplotlib.pyplot as plt
from poppy import radial_profile
import webbpsf

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter",
        help="filter to process",
        default="F770W",
        # fmt: off
        choices=["F560W", "F770W", "F1000W",
                 "F1130W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W",
                 "F1065C", "F1140C", "F1550C", "F2300C"],
        # fmt: on
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

    cfilter = args.filter
    psf_fname = f"PSFs/miri_{cfilter}_psf.fits"
    psf = fits.open(psf_fname)

    ee2 = webbpsf.measure_ee(psf, ext=3)
    ee1 = webbpsf.measure_ee(psf)

    radius, profile, ee = radial_profile(psf, ee=True, ext=0)
    radius2, profile2, ee2 = radial_profile(psf, ee=True, ext=3)

    ax.plot(radius * 0.11 * 4, ee, label="optical only")
    ax.plot(radius2 * 0.11 * 4, ee2, label="w/ detector")
    ax.set_xlim(0.0, 3.5)

    #webbpsf.display_ee(psf)
    #webbpsf.display_ee(psf, ext=3, overplot=True, ax=ax)

    psf.close()

    ax.legend()

    plt.tight_layout()

    fname = f"Figs/{filter}_webbpsf_ee_comp"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()