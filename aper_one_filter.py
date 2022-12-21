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

from photutils.detection import find_peaks
from photutils.centroids import centroid_2dg
from photutils.aperture import (
    CircularAnnulus,
    CircularAperture,
    ApertureStats,
    aperture_photometry,
)


def aper_image(filename, aprad, annrad):
    """
    Measure aperture photometry on one image for target source
    """
    hdul = fits.open(filename)
    targname = hdul[0].header["TARGNAME"]
    filter = hdul[0].header["FILTER"]
    photmjysr = hdul[1].header["PHOTMJSR"]
    print(targname, filter)
    targra = hdul[0].header["TARG_RA"]
    targdec = hdul[0].header["TARG_DEC"]

    orig_data = hdul[1].data
    orig_data[orig_data == 0.0] = np.NaN
    orig_data /= photmjysr

    w = WCS(hdul[1].header)
    coord = SkyCoord(targra, targdec, unit="deg")
    pix_coord = w.world_to_pixel(coord)

    imsize = annrad[1] * 3.0
    cutout = Cutout2D(orig_data, coord, (imsize, imsize), wcs=w)
    data = cutout.data

    # find the star
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    threshold = median + (5.0 * std)
    tbl = find_peaks(data, threshold, box_size=11)
    # tbl["peak_value"].info.format = "%.8g"  # for consistent table output
    # print(tbl[:10])  # print only the first 10 peaks

    # get the new coordinates of the star in the original image
    ncoord = cutout.wcs.pixel_to_world(tbl["x_peak"][0], tbl["y_peak"][0])

    # offset from expected position
    npix_coord = w.world_to_pixel(ncoord)
    xoff = pix_coord[0] - npix_coord[0]
    yoff = pix_coord[1] - npix_coord[1]
    print("pixel offsets from expected position: ", xoff, yoff)

    # recutout the region around the star
    cutout = Cutout2D(orig_data, ncoord, (imsize, imsize), wcs=w)
    data = cutout.data

    # find the "exact peak" to center the apertures
    # peaksize = 4
    # xy1 = int(0.5 * imsize - peaksize)
    # xy2 = int(0.5 * imsize + peaksize)
    # print(centroid_2dg(data[xy1:xy2, xy1:xy2]))
    pix_coord = centroid_2dg(data)

    aper = CircularAperture(pix_coord, r=aprad)
    annulus_aperture = CircularAnnulus(pix_coord, r_in=annrad[0], r_out=annrad[1])

    # do the aperture photometry
    phot = aperture_photometry(data, aper)
    sigclip = SigmaClip(sigma=3.0, maxiters=10)
    bkg = ApertureStats(data, annulus_aperture, sigma_clip=sigclip)
    tot_bkg = bkg.mean * aper.area
    phot["total bkg"] = tot_bkg
    phot["aperture_sum_bkgsub"] = phot["aperture_sum"] - tot_bkg
    print(phot)

    # show an image of the source and apertures used
    norm = simple_norm(data, "sqrt", percent=99.9)
    plt.imshow(data, norm=norm, interpolation="nearest")
    plt.title(f"{targname} / {filter}")
    aper.plot(color="white", lw=2, label="Photometry aperture")
    annulus_aperture.plot(color="blue", lw=2, label="Background annulus")
    plt.show()

    hdul.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter",
        help="filter to process",
        default="F770W",
        choices=["F560W", "F770W"],
        nargs=1,
    )
    parser.add_argument(
        "--dir",
        help="directory to process",
    )
    args = parser.parse_args()

    mosfiles = glob.glob(f"{args.dir}/{args.filter}/*/miri*_i2d.fits")

    aprad = 5.0
    annrad = [10.0, 15.0]
    for cfile in mosfiles:
        one_res = aper_image(cfile, aprad, annrad)
