import argparse
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import webbpsf

from aper_one_filter import aper_image

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter",
        help="filter to process",
        default="F770W",
        # fmt: off
        choices=["F560W", "F770W", "F770W_subarray", "F770W_repeat", "F1000W",
                 "F1130W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W",
                 "F1065C", "F1140C", "F1550C", "F2300C"],
        # fmt: on
    )
    parser.add_argument(
        "--saveimg", help="save images for each aperture", action="store_true"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    cfilter = args.filter
    filename = (
        f"ADwarfs/{cfilter}/BD+60 1753_set1/miri_BD+60 1753_set1_stage3_asn_i2d.fits"
        # "ADwarfs/F1000W/BD+60 1753_set1/miri_BD+60 1753_set1_stage3_bkgsub_asn_i2d.fits"
    )

    filter_fwhm = {
        "F560W": 1.0,
        "F770W": 2.27,
        "F1000W": 2.91,
        "F1130W": 3.27,
        "F1280W": 3.73,
        "F1500W": 4.36,
        "F1800W": 5.27,
        "F2100W": 6.09,
        "F2550W": 7.45,
        "F2550WR": 7.45,
        "F1065C": 2.91,
        "F1140C": 3.27,
        "F1550C": 4.36,
        "F2300C": 6.50,
    }

    # Calculate encircled energy as a function of distance for the PSF
    psf = fits.open(f"PSFs/miri_{cfilter}_psf.fits")
    mod_pixscale = psf[0].header["PIXELSCL"] * psf[0].header["DET_SAMP"]
    ee = webbpsf.measure_ee(psf)
    cradii = np.logspace(np.log10(1.0), np.log10(20.0 * filter_fwhm[cfilter]), 50)
    model_eenergy = ee(cradii * mod_pixscale)
    psf.close()

    annrad = np.array([max(cradii) * 1.1, max(cradii) * 1.2])
    apsum = np.zeros(len(cradii))
    for k, crad in enumerate(cradii):
        if args.saveimg:
            imgfile = filename.replace(".fits", f"_manyap_{crad}.png")
        else:
            imgfile = None

        cphot = aper_image(
            filename,
            crad,
            annrad,
            1.0,
            imgfile=imgfile,
        )
        apsum[k] = cphot["aperture_sum_bkgsub"][0].value
        print(crad, apsum[k])

    # calculate enclosed energy
    eenergy = apsum / np.max(apsum)

    # find the values at a fixed radius and adjust the empirical
    # to match webbpsf
    pix_rad = 30.0
    obs_val = np.interp([pix_rad], cradii, eenergy)
    mod_val = np.interp([pix_rad], cradii, model_eenergy)

    eenergy_orig = np.array(eenergy)
    eenergy += mod_val - obs_val

    print("obs:", np.interp([0.8, 0.85], eenergy, cradii))
    print("model:", np.interp([0.8, 0.85], model_eenergy, cradii))

    # make plot
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

    ax.plot(cradii, eenergy_orig, "bo", label="Observed", alpha=0.25)
    ax.plot(cradii, eenergy, "bo", label="Observed (corrected)")
    ax.plot(cradii, model_eenergy, "gs", label="WebbPSF")
    ax.plot([np.min(cradii), np.max(cradii)], [1.0, 1.0], "k:")

    ax.set_xlabel("radius [pixels]")
    ax.set_ylabel("Fractional enclosed energy")
    ax.set_title(cfilter)

    ax.legend()

    plt.tight_layout()

    fname = f"Figs/miri_many_apertures_{cfilter}"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
