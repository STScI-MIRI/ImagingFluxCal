import argparse
import webbpsf

if __name__ == '__main__':
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

    # Set input parameters:
    samp = 4  # psf oversample factor
    parity = "even"  # set psf model parity to 'odd' or 'even'
    fov = 30 * 0.11 * filter_fwhm[cfilter]  # fov for the PSF model in arcsec

    # Create a MIRI instance and calculate PSF
    miri = webbpsf.MIRI()
    miri.options['parity'] = parity
    miri.filter = cfilter
    psf = miri.calc_psf(fov_arcsec=fov, oversample=samp, add_distortion=False)
    psf.writeto(f"PSFs/miri_{cfilter}_psf.fits", overwrite=True)

    # Calculate encircled energy as a function of distance for the PSF
    # ee = webbpsf.measure_ee(psf)
    # aradii = np.arange(1.0, 10.0, 0.5)
    # ee_frac = ee(aradii)
    # print(ee_frac)
