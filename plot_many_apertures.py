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
    )
    filename_bkg = f"ADwarfs/{cfilter}/BD+60 1753_set1/miri_BD+60 1753_set1_stage3_bkgsub_asn_i2d.fits"

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
    cradii = np.logspace(
        np.log10(0.1 * filter_fwhm[cfilter]), np.log10(20.0 * filter_fwhm[cfilter]), 50
    )
    model_eenergy = ee(cradii * mod_pixscale)
    psf.close()

    if cfilter == "F2550W":
        minbkg = 1.0
        maxbkg = 1.1
    else:
        minbkg = 1.1
        maxbkg = 1.4
    annrad = np.array([max(cradii) * minbkg, max(cradii) * maxbkg])

    # get the
    tmp, ncenter = aper_image(
        filename_bkg,
        filter_fwhm[cfilter] * 5.0,
        annrad,
        1.0,
        imgfile=filename.replace(".fits", "_manyap_centerap.png"),
        return_center=True,
    )

    apsum = np.zeros(len(cradii))
    apsum_bkg = np.zeros(len(cradii))
    for k, crad in enumerate(cradii):
        if args.saveimg:
            imgfile = filename.replace(".fits", f"_manyap_{crad}.png")
            imgfile_bkg = filename_bkg.replace(".fits", f"_manyap_{crad}.png")
        else:
            imgfile = None
            imgfile_bkg = None

        cphot = aper_image(
            filename,
            crad,
            annrad,
            1.0,
            imgfile=imgfile,
            override_center=ncenter,
        )

        cphot_bkg = aper_image(
            filename_bkg,
            crad,
            annrad,
            1.0,
            imgfile=imgfile_bkg,
            override_center=ncenter,
        )
        apsum[k] = cphot["aperture_sum_bkgsub"][0].value
        apsum_bkg[k] = cphot_bkg["aperture_sum_bkgsub"][0].value
        print(crad, apsum[k], apsum_bkg[k])

    # calculate enclosed energy
    eenergy = apsum / np.max(apsum)
    eenergy_bkg = apsum_bkg / np.max(apsum)

    # find the values at a fixed radius and adjust the empirical
    # to match webbpsf
    pix_rad = 10.0 * filter_fwhm[cfilter]
    obs_val = np.interp([pix_rad], cradii, eenergy)
    obs_val_bkg = np.interp([pix_rad], cradii, eenergy_bkg)
    mod_val = np.interp([pix_rad], cradii, model_eenergy)
    print("normalizing at (rad, obs, mod)", pix_rad, obs_val, mod_val)

    eenergy_orig = np.array(eenergy)
    eenergy += mod_val - obs_val
    eenergy_bkg += mod_val - obs_val_bkg

    ee_vals = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85])
    # get the radii at fixed enclosed energy
    obs_rad_ee = np.interp(ee_vals, eenergy, cradii)
    mod_rad_ee = np.interp(ee_vals, model_eenergy, cradii)

    apcor_vals = ee_vals[0:-1] * 0.0
    print(apcor_vals)
    bkg_pix_val = (ee_vals[-1] - ee_vals[-2]) / (np.pi * (obs_rad_ee[-1] ** 2 - obs_rad_ee[-2] ** 2))
    for k, cee in enumerate(ee_vals[0:-1]):
        ee_w_bkg = ee_vals[k] - bkg_pix_val * np.pi * obs_rad_ee[k] ** 2
        apcor_vals[k] = 1.0 / ee_w_bkg

    print("apeture corrections for bkg from 0.8 to 0.85 EE")
    print(ee_vals[0:-1])
    print(apcor_vals)

    # make plot
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

    # ax.plot(cradii, eenergy_orig, "bo", label="Observed", alpha=0.25)
    ax.plot(cradii, eenergy, "b-", alpha=0.5, label="Observed (corrected)")
    ax.plot(
        cradii, eenergy_bkg, "r-", alpha=0.5, label="Observed w/ bkgsub (corrected)"
    )
    ax.plot(cradii, model_eenergy, "g-", alpha=0.5, label="WebbPSF")
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
