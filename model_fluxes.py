import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.table import QTable, vstack, hstack

import warnings
from astropy.units import UnitsWarning

from jwstabsfluxcal.Webb.read_webb import read_miri


model_names = {"1743045": "2MASS J17430448+6655015",
               "1802271": "2MASS J18022716+6043356",
               "bd60d1753": "BD+60 1753",
               "hd180609": "HD 180609",
               "hd2811": "HD 2811"}


def compute_bandflux(wave, flux_source, bwave, bandpass):
    """
    Compute the band flux given the bandpass, reference spectrum,
    and source spectrum.  Assumes a flat reference spectrum
    (for motivation see Gordon et al. 2022).

    Parameters
    ----------
    wave : nd float array
       the wavelengths of flux_source
    flux_source : nd float array
        source flux density F(lambda) as a function of wave
    bwave : nd float array
        the wavelengths of bandpass
    bandpass : nd float array
        end-to-end, total throughput bandpass of filter in fractional units
    """
    flux_source_bp = np.interp(bwave, wave, flux_source)

    # compute the the integrals
    inttop = np.trapz(bwave * bandpass * flux_source_bp)
    intbot = np.trapz(bwave * bandpass)

    return inttop / intbot


def get_band_fluxes(cfile, bandpasses, imgfile=None):
    """
    Calculated the band fluxes for each possible filter
    """
    # surpress the annoying units warning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UnitsWarning)
        modspec = QTable.read(cfile)

    mwave = (modspec["WAVELENGTH"].value * u.angstrom).to(u.micron)
    mflux = modspec["FLUX"].value * u.erg / (u.cm * u.cm * u.s * u.angstrom)
    mflux = mflux.to(u.Jy, equivalencies=u.spectral_density(mwave))

    # get the bandpasses
    bfluxes = QTable()
    rwaves = {}
    for cband in bandpasses.keys():
        rwave, cwave, ceff = bandpasses[cband]
        rwaves[cband] = rwave
        bfluxes[cband] = [compute_bandflux(mwave, mflux, cwave, ceff)]

    if imgfile is not None:
        # show an image of the source and apertures used
        fontsize = 14
        font = {"size": fontsize}
        plt.rc("font", **font)
        plt.rc("lines", linewidth=2)
        plt.rc("axes", linewidth=2)
        plt.rc("xtick.major", width=2)
        plt.rc("ytick.major", width=2)
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

        gvals = (mwave > 0.6 * u.micron) & (mwave < 32.0 * u.micron)
        ax.plot(mwave[gvals], mflux[gvals] * (mwave[gvals] ** 2), "k-", alpha=0.5)

        for cband in bfluxes.keys():
            ax.plot(
                rwaves[cband], bfluxes[cband] * (rwaves[cband] ** 2), "go", alpha=0.5
            )

        ax.set_xscale("log")
        ax.set_xlabel(r"wavelength [$\mu$m]")
        ax.set_ylabel(r"Flux [Jy $\mu$m$^2$]")

        # plt.show()
        fig.suptitle(f"{cfile}")
        plt.tight_layout()
        plt.savefig(imgfile)

    return bfluxes


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    # models to use
    modfiles = glob.glob("Models/*_mod_*.fits")

    # bandpasses to use
    bandpasses = read_miri()

    mmods = []
    for cfile in modfiles:
        ntab = QTable()
        ntab["name"] = [model_names[(cfile.split("/")[1]).split("_")[0]]]
        print(f"working on {cfile}")
        onemod = get_band_fluxes(cfile, bandpasses, cfile.replace(".fits", "_absfluxbands.png"))
        onemod = hstack([ntab, onemod])

        onemod["modfile"] = cfile

        if mmods is None:
            mmods = onemod
        else:
            mmods = vstack([mmods, onemod])

    # remove "col0" column added by vstack
    mmods.remove_column("col0")

    # save table
    mmods.write("Models/model_phot.fits", overwrite=True)
    print(mmods)
