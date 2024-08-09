import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

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
                "Models/bd60d1753_stiswfc_005.fits",
                "Models/p330e_stiswfcnic_007.fits"]

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
        mflux_Jy = mflux.to(u.Jy, equivalencies=u.spectral_density(mwave))

        # get the band fluxes
        bfluxes = QTable()
        rwaves = {}
        bwidths = {}
        for cband in bandpasses.keys():
            rwave, cwave, ceff = bandpasses[cband]
            rwaves[cband.upper()] = rwave
            bflux_lambda = compute_bandflux(mwave, mflux, cwave, ceff)
            bflux_nu = bflux_lambda.to(u.Jy, equivalencies=u.spectral_density(rwave * u.micron))
            bfluxes[cband.upper()] = [bflux_lambda]

            # bandwidth
            intresp = np.trapz(ceff, cwave)
            bwidth = intresp / max(ceff)
            bwidths[cband.upper()] = bwidth

        gvals = (mwave > 4.5 * u.micron) & (mwave < 32.0 * u.micron)

        ave_mflux = (np.average(mflux[gvals] * (mwave[gvals] ** 0))).value
        offval = 1. + k * 0.2

        ax.plot(mwave[gvals], (mflux[gvals].value * (mwave[gvals].value ** 0) / ave_mflux) * offval, "k-", alpha=0.5)

        for cband in bfluxes.keys():
            if cband in ["F1065C", "F1140C", "F1550C", "F2300C"]:
                ptype = "bo"
            else:
                ptype = "go"
            if cband != "FND":
                ax.errorbar(
                    rwaves[cband], (bfluxes[cband].value * (rwaves[cband] ** 0) / ave_mflux) * offval,
                    xerr=0.5 * bwidths[cband].value, fmt=ptype, alpha=0.5
                )

    ax.set_xscale("log")
    ax.set_xlabel(r"wavelength [$\mu$m]")
    ax.set_yscale("log")
    ax.set_ylabel(r"Normalized RJ Flux [Jy $\mu$m$^2$] + const")

    ax.text(20., 1.0, "Hot Star (G191B2B)", alpha=0.7)
    ax.text(8., 1.25, "A Dwarf (BD+60 1753)", alpha=0.7)
    ax.text(4.5, 1.42, "Solar Analog (GSPC P330-E)", alpha=0.7)

    leg = []
    leg.append(Line2D([0], [0], marker="o", color="w", label="Imaging", markerfacecolor="g", alpha=0.5, markersize=8))
    leg.append(Line2D([0], [0], marker="o", color="w", label="Coronagraphy", markerfacecolor="b", alpha=0.5, markersize=8))
    ax.legend(handles=leg, loc=(0.7, 0.7))

    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    ax.xaxis.set_minor_formatter(ScalarFormatter())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xticks([10.0])
    ax.set_xticks([5, 6, 7, 8, 9, 12., 15.0, 20.0, 25., 30.], minor=True)

    plt.tight_layout()

    fname = "Figs/miri_example_model_fluxes"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()