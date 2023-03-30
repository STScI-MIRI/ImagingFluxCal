import glob
import numpy as np
from astropy.table import QTable


if __name__ == '__main__':

    filters = ["F770W", "F1000W"]

    for cfilter in filters:
        afiles = glob.glob(f"ADwarfs/{cfilter}/BD+60 1753_set?/miri_BD+60 1753_set?_stage3_asn_i2d_apcor.dat")

        nfiles = len(afiles)
        radiis = np.zeros((7, nfiles))
        apcors = np.zeros((7, nfiles))
        for k, cfile in enumerate(afiles):
            itab = QTable.read(cfile, format="ascii.commented_header")
            ees = itab["ee"]
            radiis[:, k] = itab["radii"][0:-1]
            apcors[:, k] = itab["apcor"][0:-1]

        aradii = np.average(radiis, axis=1)
        aradii_std = np.std(radiis, axis=1)
        aapcor = np.average(apcors, axis=1)
        aapcor_std = np.std(apcors, axis=1)

        print(aradii)
        print(aradii_std)
        print(aapcor)
        print(aapcor_std)
