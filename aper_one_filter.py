import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clipped_stats, SigmaClip
import astropy.units as u
from astropy.table import QTable, vstack, hstack

from photutils.detection import find_peaks
from photutils.centroids import centroid_com
from photutils.aperture import (
    CircularAnnulus,
    CircularAperture,
    RectangularAperture,
    ApertureStats,
    aperture_photometry,
)


def aper_image(filename, aprad, annrad, imgfile=None):
    """
    Measure aperture photometry on one image for target source
    """
    hdul = fits.open(filename)
    targname = hdul[0].header["TARGNAME"]
    if targname == "J1757132":
        targname = "2MASS J17571324+6703409"
    filter = hdul[0].header["FILTER"]
    photmjysr = hdul[1].header["PHOTMJSR"]
    print(targname, filter)
    targra = hdul[0].header["TARG_RA"]
    targdec = hdul[0].header["TARG_DEC"]

    orig_data = hdul[1].data
    orig_err = hdul["ERR"].data
    orig_data[orig_data == 0.0] = np.NaN
    orig_data /= photmjysr
    orig_err /= photmjysr

    w = WCS(hdul[1].header)
    coord = SkyCoord(targra, targdec, unit="deg")
    pix_coord = w.world_to_pixel(coord)
    # check if outside the image
    if ((pix_coord[0] < 0) | (pix_coord[0] > orig_data.shape[0])) | (
        (pix_coord[1] < 0) | (pix_coord[1] > orig_data.shape[1])
    ):
        # this indicates something bad happened, currently only seen with SUB64
        # just use the full image in this case
        data = orig_data
        data_wcs = w
    else:
        imsize = annrad[1] * 6.0
        cutout = Cutout2D(orig_data, coord, (imsize, imsize), wcs=w)
        data = cutout.data
        data_wcs = cutout.wcs

    fits.writeto("test.fits", data, overwrite=True)

    # find the star
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    threshold = median + (5.0 * std)
    tbl = find_peaks(data, threshold, box_size=11)
    # tbl["peak_value"].info.format = "%.8g"  # for consistent table output
    # print(tbl[:10])  # print only the first 10 peaks

    # get the new coordinates of the star in the original image
    # use the brightest source for the new center
    sindx = np.flip(np.argsort(tbl["peak_value"]))
    ncoord = data_wcs.pixel_to_world(tbl["x_peak"][sindx[0]], tbl["y_peak"][sindx[0]])

    # offset from expected position
    npix_coord = w.world_to_pixel(ncoord)
    xoff = pix_coord[0] - npix_coord[0]
    yoff = pix_coord[1] - npix_coord[1]
    # print("pixel offsets from expected position: ", xoff, yoff)

    # recutout the region around the star
    imsize = annrad[1] * 3.0  # use a smaller size for the refined cutout
    cutout = Cutout2D(orig_data, ncoord, (imsize, imsize), wcs=w)
    cutout_err = Cutout2D(orig_err, ncoord, (imsize, imsize), wcs=w)
    data = cutout.data
    data_err = cutout_err.data

    # define for plotting
    extract_aper = RectangularAperture(npix_coord, imsize, imsize)

    # find the "exact peak" to center the apertures
    # peaksize = 4
    # xy1 = int(0.5 * imsize - peaksize)
    # xy2 = int(0.5 * imsize + peaksize)
    # print(centroid_2dg(data[xy1:xy2, xy1:xy2]))
    # print(npix_coord)
    # pix_coord = centroid_2dg(data)
    # pix_coord = centroid_1dg(data)
    pix_coord = centroid_com(data)
    # pix_coord = npix_coord

    # convert pix_coord in cutout to coordinates in the original image
    tcoord = cutout.wcs.pixel_to_world(pix_coord[0], pix_coord[1])
    full_coord = w.world_to_pixel(tcoord)

    # define apertures
    aper = CircularAperture(pix_coord, r=aprad)
    annulus_aperture = CircularAnnulus(pix_coord, r_in=annrad[0], r_out=annrad[1])

    # do the aperture photometry
    phot = aperture_photometry(data, aper, error=data_err)
    phot_stats = ApertureStats(data, aper, sigma_clip=None)

    # modify the properites of the output table
    phot.remove_column("id")

    # start with separate table to ensure specific columns are in the first columns
    tphot = QTable()
    tphot["name"] = [targname]
    tphot["filter"] = filter.upper()
    tphot["subarray"] = hdul[0].header["SUBARRAY"]
    phot = hstack([tphot, phot])

    # now add more info
    phot["aperture_sum"] *= u.DN / u.s
    phot["xcenter_full"] = full_coord[0] * u.pixel
    phot["ycenter_full"] = full_coord[1] * u.pixel

    # do background subtraction
    sigclip = SigmaClip(sigma=3.0, maxiters=10)
    bkg = ApertureStats(data, annulus_aperture, sigma_clip=sigclip)
    tot_bkg = bkg.mean * aper.area
    tot_bkg_err = bkg.std * np.sqrt(aper.area)
    phot["pix_max"] = phot_stats.max * u.DN / u.s
    phot["mean_bkg"] = bkg.mean * u.DN / u.s
    phot["aperture_area"] = aper.area
    phot["total_bkg"] = tot_bkg * u.DN / u.s
    phot["aperture_sum_bkgsub"] = phot["aperture_sum"] - phot["total_bkg"]
    phot["aperture_sum_bkgsub_err"] = (
        np.sqrt((phot["aperture_sum_err"] ** 2) + (tot_bkg_err**2)) * u.DN / u.s
    )
    phot["x_offset_from_expected"] = xoff * u.pixel
    phot["y_offset_from_expected"] = yoff * u.pixel
    phot["ngroups"] = hdul[0].header["NGROUPS"]
    phot["tgroup"] = hdul[0].header["TGROUP"]
    phot["timemid"] = hdul[0].header["EXPMID"]

    if imgfile is not None:
        # show an image of the source and apertures used
        fontsize = 14
        font = {"size": fontsize}
        plt.rc("font", **font)
        plt.rc("lines", linewidth=2)
        plt.rc("axes", linewidth=2)
        plt.rc("xtick.major", width=2)
        plt.rc("ytick.major", width=2)
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))

        norm = simple_norm(orig_data, "sqrt", percent=99)
        ax[0].imshow(orig_data, norm=norm, interpolation="nearest", origin="lower")
        extract_aper.plot(ax=ax[0], color="white")
        norm = simple_norm(data, "sqrt", percent=99.9)
        ax[1].imshow(data, norm=norm, interpolation="nearest", origin="lower")
        aper.plot(ax=ax[1], color="white", lw=2, label="Photometry aperture")
        annulus_aperture.plot(ax=ax[1], color="blue", lw=2, label="Background annulus")
        ax[1].plot(phot["xcenter"][0], phot["ycenter"][0], "b+")
        ax[1].plot([phot_stats.centroid[0]], [phot_stats.centroid[1]], "k+")

        fig.suptitle(f"{targname} / {filter}")
        plt.tight_layout()
        plt.savefig(imgfile)
        # plt.show()
    hdul.close()

    return phot


def aper_one_filter(subdir, filter):
    """
    Do aperture photometry on all mosaic files for one filter and one class
    of stars.
    """
    mosfiles = glob.glob(f"{subdir}/{filter}/*/miri*_i2d.fits")
    print(f"{subdir}/{filter}/*/miri*_i2d.fits")

    aprad = 5.0
    annrad = [10.0, 15.0]
    mres = None
    for cfile in mosfiles:
        one_res = aper_image(
            cfile, aprad, annrad, imgfile=cfile.replace(".fits", "_absfluxapers.png")
        )
        if mres is None:
            mres = one_res
        else:
            mres = vstack([mres, one_res])

    # save table
    mres.write(f"{subdir}/{filter}_phot.fits", overwrite=True)
    print(mres)


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
        "--dir",
        choices=["HotStars", "ADwarfs", "SolarAnalogs"],
        default="ADwarfs",
        help="directory to process",
    )
    args = parser.parse_args()

    aper_one_filter(args.dir, args.filter)
