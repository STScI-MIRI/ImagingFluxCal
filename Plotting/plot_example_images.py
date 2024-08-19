import argparse
import warnings
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import numpy as np

import astropy.units as u
from astropy.visualization import simple_norm
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.io import fits
from astropy.wcs import FITSFixedWarning
from astropy.time import Time
from astropy.stats import sigma_clipped_stats

from photutils.detection import find_peaks

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--coron", help="show coronagraphic PSFs", action="store_true")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    if args.coron:
        fontsize = 18
        nrows = 2
        ncols = 2
        imsize = 150.
        fsize = (10, 10)
    else:
        fontsize = 16
        nrows = 2
        ncols = 5
        imsize = 60.
        fsize = (14, 6)

    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=fsize)

    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W", "FND", "F1500W",
               "F1800W", "F2100W", "F2550W"]
    if args.coron:
        filters = ["F1065C", "F1140C", "F1550C", "F2300C"]
    for k, cfilter in enumerate(filters):
        if cfilter in ["F1500W", "F1800W", "F2100W", "F2550W"]:
            fname = f"ADwarfs/{cfilter}/del UMi_set1/miri_del UMi_set1_stage3_bkgsub_asn_i2d.fits"
        elif cfilter in ["F560W", "F770W", "F1000W", "F1130W", "F1280W"]:
            fname = f"ADwarfs/{cfilter}/BD+60 1753_set1/miri_BD+60 1753_set1_stage3_bkgsub_asn_i2d.fits"
        elif cfilter == "FND":
            fname = f"ADwarfs/{cfilter}/del UMi_set1/miri_del UMi_set1_stage3_asn_i2d.fits"
        else:
            fname = f"ADwarfs/{cfilter}/HD 2811_set1/miri_HD 2811_set1_stage3_bkgsub_asn_i2d.fits"
        hdul = fits.open(fname)
        orig_data = hdul[1].data
        targra = hdul[0].header["TARG_RA"]
        targdec = hdul[0].header["TARG_DEC"]

        # suppress warning given *every* time
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            w = WCS(hdul[1].header)
        orig_coord = SkyCoord(
            targra,
            targdec,
            unit="deg",
            pm_ra_cosdec=hdul[0].header["MU_RA"] * u.arcsec / u.yr,
            pm_dec=hdul[0].header["MU_DEC"] * u.arcsec / u.yr,
        )
        new_obstime = Time(hdul[0].header["DATE-BEG"])
        orig_coord.obstime = Time(hdul[0].header["MU_EPOCH"])
        hdul.close()

        coord = orig_coord
        cutout = Cutout2D(orig_data, coord, (imsize, imsize), wcs=w, fill_value=0.0, mode="partial")
        data = cutout.data
        data_wcs = cutout.wcs

        # find the star
        mean, median, std = sigma_clipped_stats(data, sigma=3.0)
        threshold = median + (5.0 * std)
        tbl = find_peaks(data, threshold, box_size=11)

        # get the new coordinates of the star in the original image
        # use the brightest source for the new center
        sindx = np.flip(np.argsort(tbl["peak_value"]))
        ncoord = data_wcs.pixel_to_world(tbl["x_peak"][sindx[0]], tbl["y_peak"][sindx[0]])

        cutout = Cutout2D(orig_data, ncoord, (imsize, imsize), wcs=w)
        data = cutout.data
        data_wcs = cutout.wcs

        data[np.isnan(data)] = np.median(data[data > 0.0])

        norm = simple_norm(data, "sqrt", percent=99)
        ll = k // ncols
        mm = k % ncols
        ax[ll, mm].imshow(data, norm=norm, interpolation="nearest", origin="lower")
        ax[ll, mm].text(0.05, 0.95, cfilter, horizontalalignment='left', color="white",
                        verticalalignment='top', transform=ax[ll, mm].transAxes,
                        fontsize=fontsize)

        scalebar = AnchoredSizeBar(ax[ll, mm].transData,
                                1./0.11, '1"', 'lower right', 
                                pad=0.1,
                                color='white',
                                frameon=False,
                                size_vertical=1)

        ax[ll, mm].add_artist(scalebar)

        ax[ll, mm].set_yticks([])
        ax[ll, mm].set_xticks([])

    plt.tight_layout()

    fname = "Figs/example_images"
    if args.coron:
        fname = f"{fname}_coron"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
