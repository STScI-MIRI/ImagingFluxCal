import argparse
import copy
import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from astropy.convolution import convolve_fft  # , Box2DKernel, convolve
import webbpsf

# from poppy import utils
from matplotlib.gridspec import GridSpec


def add_cruciform(psf_webbpsf, f):
    """
    Code adapted from CAP202 notebook by Mattia L.
    """
    # Amplitude and scale for the cruciform modelling of commissioning data
    ampl = {"F560W": 0.00773238, "F770W": 0.00877061}
    scal = {"F560W": 0.01895400, "F770W": 0.02492824}
    jump = {"F560W": 3.0, "F770W": 3.0}

    # Taken from the WebbPSF source code
    ext = 0
    test_psf = copy.deepcopy(psf_webbpsf)

    oversample = test_psf[ext].header["DET_SAMP"]

    x = np.arange(test_psf[ext].data.shape[1], dtype=float)
    x -= (test_psf[ext].data.shape[1] - 1) / 2
    x /= oversample

    kernel_x = ampl[f] * np.exp(-np.abs(x) * scal[f]) / oversample
    kernel_x[np.abs(x) < jump[f]] = 0.0
    kernel_x.shape = (1, test_psf[ext].data.shape[1])
    im_conv_x = convolve_fft(
        test_psf[ext].data,
        kernel_x,
        boundary="fill",
        fill_value=0.0,
        normalize_kernel=False,
        nan_treatment="fill",
        allow_huge=True,
    )
    kernel_y = ampl[f] * np.exp(-np.abs(x) * scal[f]) / oversample
    kernel_y[np.abs(x) < jump[f]] = 0.0
    kernel_y.shape = (1, test_psf[ext].data.shape[1])
    im_conv_y = convolve_fft(
        test_psf[ext].data,
        kernel_y.T,
        boundary="fill",
        fill_value=0.0,
        normalize_kernel=False,
        nan_treatment="fill",
        allow_huge=True,
    )

    # compute and save convolution kernel
    im_conv_both = (im_conv_x + im_conv_y) / (oversample ** 2)
    fits.writeto(f"PSFs/cruciform_kernel_{f}.fits", im_conv_both, overwrite=True)

    psf_new = test_psf[ext].data + im_conv_both

    # keep the the total energy in the PSF as originally computed the same
    psf_new *= test_psf[ext].data.sum() / psf_new.sum()

    test_psf[ext].data = copy.deepcopy(psf_new)

    return test_psf

#    psf_data_wcruciform = copy.deepcopy(psf_new)

#    psf_webbpsf_all[f][ext + 1].data = utils.rebin_array(
#        psf_data_wcruciform, rc=(oversample, oversample)
#    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter",
        help="filter to process",
        default="F560W",
        choices=["F560W", "F770W"],
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    f = args.filter

    # read in the webbpsf file
    psf_webbpsf = fits.open(f"PSFs/miri_{f}_psf.fits")
    psf_webbpsf_wcurciform = add_cruciform(psf_webbpsf, f)
    psf_webbpsf_wcurciform.writeto(f"PSFs/miri_{f}_psf_wcurciform.fits", overwrite=True)

    fig = plt.figure(constrained_layout=True, figsize=(20, 5))
    gs = GridSpec(1, 6, figure=fig)

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])
    ax3 = fig.add_subplot(gs[0, 3:5])

    ext = 0
    if (ext == 0) | (ext == 2):
        over = 4.0

    psf_data = psf_webbpsf[ext].data
    psf_data_wcruciform = psf_webbpsf_wcurciform[ext].data

    ax0.imshow(
        psf_data, vmin=-0.00001, vmax=0.00001, origin="lower",
    )
    ax1.imshow(psf_data_wcruciform, vmin=-0.00001, vmax=0.00001, origin="lower")
    ax2.imshow(
        psf_data_wcruciform - psf_data, vmin=-0.00001, vmax=0.00001, origin="lower",
    )
    # ax2.axhline(500)
    # ax2.axvline(500)
    ax0.set_xlabel("X")
    ax0.set_ylabel("Y")
    ax1.set_xlabel("X")
    ax2.set_xlabel("X")
    ax0.set_title(f + "\nORIGINAL i2d WEBBPSF\n OVER SAMPLED")
    ax1.set_title(f + "\nORIGINAL i2d WEBBPSF + CRUCIFORM")
    ax2.set_title(f + " CRUCIFOR")

    new = webbpsf.measure_ee(psf_webbpsf, ext=ext)
    i2d_rp = webbpsf.measure_ee(psf_webbpsf_wcurciform, ext=ext)
    # cal_rp = webbpsf.measure_ee(psf_webbpsf_wcurciform, ext=ext + 2)

    r = np.arange(0, 100, 0.001)
    ax3.plot(r, i2d_rp(r), label="i2d")
    # ax3.plot(r, cal_rp(r), label="cal")
    ax3.plot(r, new(r), label="new")
    ax3.legend(loc="lower right")
    ax3.set_xlabel("Radius")
    ax3.set_ylabel("Encircled energy")

    plt.tight_layout()

    fname = f"PSFs/webbpsf_w_wo_cruciform_{args.filter}"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
