import argparse
import numpy as np
import matplotlib.pyplot as plt

from astropy.table import QTable
import astropy.units as u

import warnings
from astropy.units import UnitsWarning

from jwstabsfluxcal.Webb.read_webb import read_miri
from model_fluxes import compute_bandflux

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # models to use
    modfiles = ["Models/g191b2b_stiswfcnic_004.fits",
                "Models/bd60d1753_stiswfc_004.fits",
                "Models/p330e_stiswfcnic_005.fits"]

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

    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

    for k, cfile in enumerate(modfiles):

        # surpress the annoying units warning
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UnitsWarning)
            modspec = QTable.read(cfile)

        mwave = (modspec["WAVELENGTH"].value * u.angstrom).to(u.micron)
        mflux = modspec["FLUX"].value * u.erg / (u.cm * u.cm * u.s * u.angstrom)
        mflux = mflux.to(u.Jy, equivalencies=u.spectral_density(mwave))

        # get the band fluxes
        bfluxes = QTable()
        rwaves = {}
        for cband in bandpasses.keys():
            rwave, cwave, ceff = bandpasses[cband]
            rwaves[cband.upper()] = rwave
            bfluxes[cband.upper()] = [compute_bandflux(mwave, mflux, cwave, ceff)]
            print(cfile, cband, bfluxes[cband.upper()])

        gvals = (mwave > 4.5 * u.micron) & (mwave < 32.0 * u.micron)

        ave_mflux = (np.average(mflux[gvals] * (mwave[gvals] ** 2))).value
        offval = k * 0.2

        ax.plot(mwave[gvals], mflux[gvals].value * (mwave[gvals].value ** 2) / ave_mflux + offval, "k-", alpha=0.5)

        for cband in bfluxes.keys():
            ax.plot(
                rwaves[cband], bfluxes[cband].value * (rwaves[cband] ** 2) / ave_mflux + offval, "go", alpha=0.5
            )

    ax.set_xscale("log")
    ax.set_xlabel(r"wavelength [$\mu$m]")
    ax.set_ylabel(r"Normalized RJ Flux [Jy $\mu$m$^2$] + const")

    plt.tight_layout()

    fname = "Figs/miri_example_model_fluxes"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()