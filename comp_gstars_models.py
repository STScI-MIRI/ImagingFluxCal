import numpy as np
from astropy.table import QTable

if __name__ == "__main__":

    caltab = QTable.read("Models/model_phot.fits")
    gcaltab = QTable.read("Models/model_w_grieke_phot.fits")

    # get the ratio of fluxes between the models in two bands
    print("****")
    print("ratio between models in same band")
    compbands = ["F1280W", "F1800W"]
    compstars = ["16 Cyg B", "GSPC P330-E"]
    for cstar in compstars:
        (mindx,) = np.where(caltab["name"] == cstar)
        (gmindx,) = np.where(gcaltab["name"] == cstar)
        for cband in compbands:
            cflux = caltab[cband][mindx[0]]
            gflux = gcaltab[cband][gmindx[0]]
            print(f"{cstar}, {cband}: grieke/calspec, {gflux/cflux}")

    print("****")
    print("ratio between stars for same band")
    # get the ratio of the models between the two bands
    print("-----")
    print("observations")
    # read in observed fluxes
    cstar1 = compstars[0]
    cstar2 = compstars[1]
    obsratio = []
    for cband in compbands:
        obstab = QTable.read(f"SolarAnalogs/{cband}_phot.fits")
        (mindx1,) = np.where(obstab["name"] == cstar1)
        (mindx2,) = np.where(obstab["name"] == cstar2)
        cflux1 = obstab["aperture_sum_bkgsub"][mindx1[0]]
        cflux2 = obstab["aperture_sum_bkgsub"][mindx2[0]]
        obsratio.append(cflux1 / cflux2)
        print(f"{cband}: calspec {cstar1}/{cstar2}, {cflux1/cflux2}")

    print("-----")
    print("models")
    for k, cband in enumerate(compbands):
        (mindx,) = np.where(caltab["name"] == cstar1)
        (gmindx,) = np.where(gcaltab["name"] == cstar1)
        (mindx2,) = np.where(caltab["name"] == cstar2)
        (gmindx2,) = np.where(gcaltab["name"] == cstar2)

        cflux1 = caltab[cband][mindx[0]]
        cflux2 = caltab[cband][mindx2[0]]
        cratio = cflux1 / cflux2
        print(f"{cband}: calspec {cstar1}/{cstar2}, {cratio}, mod/obs: {cratio/obsratio[k]}")

        gcflux1 = gcaltab[cband][gmindx[0]]
        gcflux2 = gcaltab[cband][gmindx2[0]]
        gcratio = gcflux1 / gcflux2
        print(f"{cband}: grieke {cstar1}/{cstar2}, {gcratio}, mod/obs: {gcratio/obsratio[k]}")

    print("****")
