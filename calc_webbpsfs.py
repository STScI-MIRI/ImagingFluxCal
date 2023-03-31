import argparse
import webbpsf

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
    args = parser.parse_args()
    cfilter = args.filter

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

    # Set input parameters:
    samp = 4  # psf oversample factor
    parity = "even"  # set psf model parity to 'odd' or 'even'
    fov = 30 * 0.11 * filter_fwhm[cfilter]  # fov for the PSF model in arcsec

    # Create a MIRI instance and calculate PSF
    miri = webbpsf.MIRI()
    miri.options["parity"] = parity
    miri.filter = cfilter

    shiftx = 74.0
    shifty = 72.0
    if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
        miri.options["parity"] = "odd"
        miri.options["source_offset_x"] = shiftx * 0.11
        miri.options["source_offset_y"] = -shifty * 0.11

        if cfilter in ["F1065C", "F1140C", "F1550C"]:
            miri.pupil_mask = "MASKFQPM"
            miri.image_mask = f"FQPM{cfilter[1:-1]}"
        elif cfilter == "F2300C":
            miri.pupil_mask = "MASKLYOT"
            miri.image_mask = "LYOT2300"

        psf = miri.calc_psf(
            fov_pixels=300,
            oversample=samp,
            add_distortion=False,
            normalize="exit_pupil",
        )
    else:
        psf = miri.calc_psf(fov_arcsec=fov, oversample=samp, add_distortion=False)
    psf.writeto(f"PSFs/miri_{cfilter}_psf.fits", overwrite=True)

    # Calculate encircled energy as a function of distance for the PSF
    # ee = webbpsf.measure_ee(psf)
    # aradii = np.arange(1.0, 10.0, 0.5)
    # ee_frac = ee(aradii)
    # print(ee_frac)
