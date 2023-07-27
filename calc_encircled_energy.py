import argparse
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.table import QTable

import webbpsf

from aper_one_filter import aper_image

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset",
        default="BD+60 1753_set1",
        help="use this MOS file instead of BD+60 1753_set1",
    )
    parser.add_argument(
        "--dir",
        choices=["HotStars", "ADwarfs", "SolarAnalogs", "all"],
        default="ADwarfs",
        help="directory to process",
    )
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
    parser.add_argument(
        "--fwhmfac",
        default=20.0,
        help="FWHM factor for normlization of empirical PSF to WebbPSF PSF",
        type=float,
    )
    parser.add_argument(
        "--saveimg", help="save images for each aperture", action="store_true"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    cfilter = args.filter
    filename = (
        f"{args.dir}/{cfilter}/{args.dataset}/miri_{args.dataset}_stage3_asn_i2d.fits"
    )
    print(f"mosaic filename = {filename}")
    filename_bkg = filename.replace("stage3", "stage3_bkgsub")

    # in 0.11 arcsec pixels
    filter_fwhm = {
        "F560W": 1.636,
        "F770W": 2.187,
        "F1000W": 2.888,
        "F1065C": 2.910,
        "F1130W": 3.318,
        "F1140C": 3.270,
        "F1280W": 3.713,
        "F1500W": 4.354,
        "F1550C": 4.360,
        "F1800W": 5.224,
        "F2100W": 5.989,
        "F2300C": 6.090,
        "F2550W": 7.312,
    }

    # Calculate encircled energy as a function of distance for the PSF
    if cfilter in ["F560W", "F770W"]:
        psf_fname = f"PSFs/miri_{cfilter}_psf_wcurciform.fits"
    else:
        psf_fname = f"PSFs/miri_{cfilter}_psf.fits"
    psf = fits.open(psf_fname)
    mod_pixscale = psf[0].header["PIXELSCL"] * psf[0].header["DET_SAMP"]
    if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
        wcenter = (896, 315)
    else:
        wcenter = None
    ee = webbpsf.measure_ee(psf, center=wcenter)
    psf.close()

    norm_factor = args.fwhmfac
    minbkg = 1.0
    maxbkg = 1.2

    cradii = np.logspace(
        np.log10(0.1 * filter_fwhm[cfilter]),
        np.log10(2 * norm_factor * filter_fwhm[cfilter]),
        100
    )
    model_eenergy = ee(cradii * mod_pixscale)

    annrad = np.array([max(cradii) * minbkg, max(cradii) * maxbkg])

    # get the center for all the photometry
    tmp, ncenter = aper_image(
        filename_bkg,
        filter_fwhm[cfilter] * 5.0,
        [filter_fwhm[cfilter] * 5.0, filter_fwhm[cfilter] * 6.0],
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
            filename, crad, annrad, 1.0, imgfile=imgfile, override_center=ncenter,
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
    eenergy = apsum / np.nanmax(apsum)
    eenergy_bkg = apsum_bkg / np.nanmax(apsum)

    # find the values at a fixed radius and adjust the empirical
    # to match webbpsf
    pix_rad = norm_factor * filter_fwhm[cfilter]
    obs_val = np.interp([pix_rad], cradii, eenergy)
    obs_val_bkg = np.interp([pix_rad], cradii, eenergy_bkg)
    mod_val = np.interp([pix_rad], cradii, model_eenergy)
    print("normalizing at (rad, obs, mod)", pix_rad, obs_val, mod_val)

    eenergy_orig = np.array(eenergy)
    eenergy *= mod_val / obs_val
    eenergy_bkg *= mod_val / obs_val_bkg

    # create the final ee profile
    # observed to the pix_rad radius and model for the rest
    if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
        fin_eenergy = np.array(eenergy_bkg)
    else:
        fin_eenergy = np.array(eenergy)
    gvals = cradii > pix_rad
    fin_eenergy[gvals] = model_eenergy[gvals]
    if cfilter == "F2550W":
        fin_eenergy = np.array(model_eenergy)

    # create the table to save the ee values
    atab = QTable()
    atab["radius"] = cradii
    atab["ee"] = fin_eenergy
    atab["ee_obs"] = eenergy
    atab["ee_obs_bkg"] = eenergy_bkg
    atab["ee_model"] = model_eenergy
    atab.write(
        filename.replace(".fits", f"_ee_fwhmfac{args.fwhmfac}.dat"),
        format="ascii.commented_header",
        overwrite=True,
    )

    ee_vals = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85])
    # get the radii at fixed enclosed energy
    obs_rad_ee = np.interp(ee_vals, fin_eenergy, cradii)
    mod_rad_ee = np.interp(ee_vals, model_eenergy, cradii)

    apcor_vals = ee_vals[0:-1] * 0.0
    bkg_pix_val = (ee_vals[-1] - ee_vals[-2]) / (
        np.pi * (obs_rad_ee[-1] ** 2 - obs_rad_ee[-2] ** 2)
    )
    for k, cee in enumerate(ee_vals[0:-1]):
        ee_w_bkg = ee_vals[k] - bkg_pix_val * np.pi * obs_rad_ee[k] ** 2
        apcor_vals[k] = 1.0 / ee_w_bkg

    print("apeture corrections for bkg from 0.8 to 0.85 EE")
    print(obs_rad_ee)
    print(ee_vals[0:-1])
    print(apcor_vals)

    # create the table to save the values in
    atab = QTable()
    atab["ee"] = ee_vals
    atab["radii"] = obs_rad_ee
    atab["apcor"] = np.concatenate([apcor_vals, [1.0]])
    atab.write(
        filename.replace(".fits", f"fwhmfac{args.fwhmfac}_apcor.dat"),
        format="ascii.commented_header",
        overwrite=True,
    )

    # make plot
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

    ax.plot(cradii, eenergy_orig, "b:", label="Observed", alpha=0.25)
    ax.plot(
        cradii, eenergy_bkg, "r-", alpha=0.5, label="Observed w/ bkgsub (corrected)"
    )
    ax.plot(cradii, model_eenergy, "g-", alpha=0.5, label="WebbPSF")
    ax.plot(cradii, fin_eenergy, "b-", alpha=0.5, label="Final (Obs+WebbPSF)")
    ax.plot([np.min(cradii), np.max(cradii)], [1.0, 1.0], "k:")

    ax.set_xlabel("radius [pixels]")
    ax.set_ylabel("Fractional enclosed energy")
    ax.set_title(f"{cfilter} / FWHMFAC = {args.fwhmfac}")

    ax.set_ylim(0.0, 1.1)

    ax.legend()

    plt.tight_layout()

    fname = filename.replace(".fits", "_ee")
    fname = f"{fname}_fwhmfac{args.fwhmfac}"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
