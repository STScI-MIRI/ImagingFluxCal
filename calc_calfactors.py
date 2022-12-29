import matplotlib.pyplot as plt
import numpy as np

from astropy.table import QTable


def get_calfactors(dir, filter):
    """
    Read in the observed and mdoel fluxes and computer the calibration factors
    """
    # read in observed fluxes
    obstab = QTable.read(f"{dir}/{filter}_phot.fits")
    # read in model fluxes
    modtab = QTable.read("Models/model_phot.fits")

    mfluxes = []
    cfactors = []
    cfactors_unc = []
    for k, cname in enumerate(obstab["name"]):
        oflux = obstab["aperture_sum_bkgsub"][k]
        oflux_unc = obstab["aperture_sum_bkgsub_err"][k]

        mindx, = np.where(modtab["name"] == cname)
        mflux = modtab[filter][mindx[0]]

        cfactor = mflux.value / oflux.value
        cfactor_unc = (oflux_unc / oflux) * cfactor
        mfluxes.append(mflux.value)
        cfactors.append(cfactor)
        cfactors_unc.append(cfactor_unc)

    return (cfactors, cfactors_unc, mfluxes)


if __name__ == '__main__':

    dir = "ADwarfs"
    filter = "F770W"

    cfactors, cfactors_unc, mfluxes = get_calfactors(dir, filter)

    # make plot
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

    ax.errorbar(mfluxes, cfactors, yerr=cfactors_unc, fmt="bo", alpha=0.5)

    ax.set_xscale("log")
    ax.set_xlabel("Flux [Jy]")
    ax.set_ylabel("Calibration Factors [Jy / (DN/s)]")
    ax.set_title(filter)

    plt.tight_layout()
    plt.show()
