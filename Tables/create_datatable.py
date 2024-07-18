import os
import numpy as np
from astropy.table import QTable
import astropy.units as u


if __name__ == "__main__":

    # fmt: off
    filters = ["F560W", "F770W", "F1000W",
               "F1130W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W",
               "F1065C", "F1140C", "F1550C", "F2300C"]
    # fmt: on

    dirs = ["HotStars", "ADwarfs", "SolarAnalogs"]

    otab = QTable(
        names=(
            "name",
            "filter",
            "subarray",
            "timemid",
            "inst flux",
            "inst flux unc",
            "background"
        ),
        units=("", "", "", u.d, "DN/s", "DN/s", "DN/s"),
        dtype=("str", "str", "str", "float64", "float64", "float64", "float64"),
    )

    for cdir in dirs:
        for cfilter in filters:

            fname = f"{cdir}/{cfilter}_eefrac0.7_phot.fits"
            if os.path.isfile(fname): 
                tab = QTable.read(fname)
                # print(tab.colnames)

                for cline in tab:
                    otab.add_row(
                        (
                            cline["name"],
                            cline["filter"],
                            cline["subarray"],
                            cline["timemid"],
                            cline["aperture_sum_bkgsub"],
                            cline["aperture_sum_bkgsub_err"],
                            cline["mean_bkg"]
                        )
                    )

        # sort by name
        sindxs = np.argsort(otab["name"])
        otab = otab[sindxs]
        print(otab)
        exit()
