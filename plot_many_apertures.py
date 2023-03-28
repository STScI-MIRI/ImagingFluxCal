import argparse
import matplotlib.pyplot as plt
import numpy as np

import webbpsf

from aper_one_filter import aper_image

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filename = (
        "ADwarfs/F1000W/BD+60 1753_set1/miri_BD+60 1753_set1_stage3_bkgsub_asn_i2d.fits"
    )
    cfilter = "F1000W"

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
    samp = 5  # psf oversample factor
    parity = 'odd'   # set psf model parity to 'odd' or 'even'
    fov = 12  # fov for the PSF model in arcsec

    # Create a MIRI instance and calculate PSF
    # miri = webbpsf.MIRI()
    # miri.options['parity'] = parity
    # miri.filter = cfilter
    # psf = miri.calc_psf(fov_arcsec=fov, oversample=samp, add_distortion=True)

    cradii = np.logspace(np.log10(1.0), np.log10(20.0 * filter_fwhm[cfilter]), 50)
    annrad = np.array([max(cradii)*1.1, max(cradii)*1.2])
    apsum = np.zeros(len(cradii))
    for k, crad in enumerate(cradii):
        cphot = aper_image(filename, crad, annrad, 1.0,
                           imgfile=filename.replace(".fits", f"_manyap_{crad}.png"))
        apsum[k] = cphot["aperture_sum_bkgsub"][0].value
        print(crad, apsum[k])

    # calculate enclosed energy
    eenergy = apsum / np.max(apsum)

    print(np.interp([0.8], eenergy, cradii))

    # make plot
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14, 8))

    ax.plot(cradii, eenergy, "ko")
    ax.plot([np.min(cradii), np.max(cradii)], [1.0, 1.0], "k:")

    ax.set_xlabel("radius [arcsec]")
    ax.set_ylabel("Fractional enclosed energy")
    ax.set_title(cfilter)

    plt.tight_layout()

    fname = f"Figs/miri_many_apertures_{cfilter}"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
