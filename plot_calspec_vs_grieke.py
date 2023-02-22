import warnings
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.units import UnitsWarning
from astropy.table import QTable
import astropy.units as u
from scipy import interpolate


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--star",
        help="star to plot",
        default="p330e",
        choices=["p330e", "16cygb"],
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    if args.star == "16cygb":
        cfile = "Models/16cygb_stis_003.fits"
        gcfile = "Models/grieke_16cygb_final.csv"
        mfac = 1e-3
    else:
        cfile = "Models/p330e_stiswfcnic_005.fits"
        gcfile = "Models/grieke_p330e_final.csv"
        mfac = 1.0

    print(cfile, gcfile)

    # get CALSPEC model
    # surpress the annoying units warning
    # cfile = "Models/p330e_mod_005.fits"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UnitsWarning)
        modspec = QTable.read(cfile)

    mwave = (modspec["WAVELENGTH"].value * u.angstrom).to(u.micron)
    mflux = modspec["FLUX"].value * u.erg / (u.cm * u.cm * u.s * u.angstrom)
    mflux = mflux.to(u.Jy, equivalencies=u.spectral_density(mwave))

    # get grieke model
    gmodspec = QTable.read(gcfile, format="ascii.csv")
    gmwave = gmodspec["wave_um"].value * u.micron
    gmflux = gmodspec["flux_W_cm-2_um-1"].value * u.W / (u.cm * u.cm * u.micron)
    gmflux = gmflux.to(u.Jy, equivalencies=u.spectral_density(gmwave))
    gmflux *= mfac

    # make plot
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

    # get ratios
    f = interpolate.interp1d(gmwave, gmflux)

    gvals = (mwave > 0.6 * u.micron) & (mwave < 31.0 * u.micron)
    ax.plot(mwave[gvals], mflux[gvals] * mwave[gvals] ** 2, "b-", label="CALSPEC")

    gvals = (gmwave > 0.6 * u.micron) & (gmwave < 31.0 * u.micron)
    ax.plot(gmwave[gvals], gmflux[gvals] * gmwave[gvals] ** 2, "g-", label="GRIEKE")

    gvals = (mwave > 0.6 * u.micron) & (mwave < 5.0 * u.micron)
    igmflux = f(mwave)
    ratio = mflux[gvals].value / igmflux[gvals]
    averatio = np.average(ratio)
    ax.text(0.4, 0.2, fr"CALSPEC/GRIEKE (0.6-5 $\mu$m) = {averatio:.3f}", transform=ax.transAxes)

    gvals = (mwave > 5.0 * u.micron) & (mwave < 28 * u.micron)
    ratio = mflux[gvals].value / igmflux[gvals]
    averatio = np.average(ratio)
    ax.text(0.4, 0.15, fr"CALSPEC/GRIEKE (5-28 $\mu$m) = {averatio:.3f}", transform=ax.transAxes)

    ax.set_xscale("log")
    ax.set_xlabel(r"wavelength [$\mu$m]")
    ax.set_ylabel(r"$\lambda^2 F(\nu)$ [Jy $\mu$m$^2$]")

    ax.legend(fontsize=0.8 * fontsize, ncol=2)

    plt.tight_layout()

    fname = f"comp_calspec_vs_grieke_{args.star}"
    if args.png:
        fig.savefig(f"{fname}.png")
    elif args.pdf:
        fig.savefig(f"{fname}.pdf")
    else:
        plt.show()
