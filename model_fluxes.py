import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.table import QTable, vstack, hstack

import warnings
from astropy.units import UnitsWarning

from jwstabsfluxcal.Webb.read_webb import read_miri


model_names = {
    "g191b2b": "G 191-B2B",
    "gd71": "GD 71",
    "1743045": "2MASS J17430448+6655015",
    "1757132": "2MASS J17571324+6703409",
    "1802271": "2MASS J18022716+6043356",
    "bd60d1753": "BD+60 1753",
    "hd163466": "HD 163466",
    "hd167060": "HD 167060",
    "hd180609": "HD 180609",
    "hd106252": "HD 106252",
    "hd142331": "HD 142331",
    "hd37962": "HD 37962",
    "hd2811": "HD 2811",
    "p177d": "GSPC P177-D",
    "p330e": "GSPC P330-E",
    "16cygb": "16 Cyg B",
    "delumi": "del UMi",
    "10lac": "10 Lac",
}


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

    # replace with G. Rieke's models when available
    cname = cfile.split("/")[1].split("_")[0]
    if cname in ["xxp330e", "xx16cygb"]:
        if cname == "16cygb":
            gcfile = "Models/grieke_16cygb_final.csv"
            mfac = 1e-3
        else:
            gcfile = "Models/grieke_p330e_final.csv"
            mfac = 1.0
        print(f"using grieke model for {cname}")
        modspec = QTable.read(gcfile, format="ascii.csv")
        mwave = modspec["wave_um"].value * u.micron
        mflux = modspec["flux_W_cm-2_um-1"].value * u.W / (u.cm * u.cm * u.micron)
        mflux = mflux.to(u.Jy, equivalencies=u.spectral_density(mwave))
        mflux *= mfac

        # save the file with Jy units
        otab = QTable()
        otab["wave"] = mwave
        otab["flux"] = mflux
        otab.write(gcfile.replace(".csv", ".fits"), overwrite=True)
        otab.write(gcfile.replace(".csv", ".dat"), format="ascii.ipac", overwrite=True)
    else:
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
        rwaves[cband.upper()] = rwave
        bfluxes[cband.upper()] = [compute_bandflux(mwave, mflux, cwave, ceff)]

    if imgfile is not None:
        fontsize = 14
        font = {"size": fontsize}
        plt.rc("font", **font)
        plt.rc("lines", linewidth=2)
        plt.rc("axes", linewidth=2)
        plt.rc("xtick.major", width=2)
        plt.rc("ytick.major", width=2)
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

        gvals = (mwave > 0.01 * u.micron) & (mwave < 32.0 * u.micron)
        ax.plot(mwave[gvals], mflux[gvals] * (mwave[gvals] ** 2), "k-", alpha=0.5)

        for cband in bfluxes.keys():
            ax.plot(
                rwaves[cband], bfluxes[cband] * (rwaves[cband] ** 2), "go", alpha=0.5
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
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
    # modfiles = glob.glob("Models/*_mod_*.fits")
    modfiles = glob.glob("Models/*_stis*.fits")

    # bandpasses to use
    bandpasses = read_miri()

    # save the band response functions
    for cband in bandpasses.keys():
        rwave, cwave, ceff = bandpasses[cband]
        otab = QTable()
        otab["wave"] = cwave
        otab["bandpass"] = ceff
        otab.write(
            f"Models/bandpass_{cband}.dat",
            format="ascii.commented_header",
            overwrite=True,
        )

    mmods = []
    for cfile in modfiles:
        ntab = QTable()
        ntab["name"] = [model_names[(cfile.split("/")[1]).split("_")[0]]]
        print(f"working on {cfile}")
        onemod = get_band_fluxes(
            cfile, bandpasses, cfile.replace(".fits", "_absfluxbands.png")
        )
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
