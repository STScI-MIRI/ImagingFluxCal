import argparse

from astropy.table import QTable, hstack


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter",
        help="filter to process",
        default="F770W",
        # fmt: off
        choices=["F560W", "F770W", "F1000W",
                 "F1130W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W",
                 "F1065C", "F1140C", "F1550C", "F2300C",
                 "FND"],
        # fmt: on
    )
    args = parser.parse_args()

    ofile = "CalFactors/jwst_miri_photom_0201.fits"
    nfile = "Photom/jwst_miri_photom_flight_30aug24.fits"

    print("Old")
    tab = QTable.read(ofile, hdu=1)
    tab2 = QTable.read(ofile, hdu=2)
    gvals = tab["filter"] == args.filter
    tab3 = QTable()
    tab3["p amp"] = 100.0 * tab2[gvals]["amplitude"] / tab[gvals]["photmjsr"]
    print(hstack([tab[gvals], tab2[gvals], tab3]))

    print("New")
    tab = QTable.read(nfile, hdu=1)
    tab2 = QTable.read(nfile, hdu=2)
    gvals = tab["filter"] == args.filter
    tab3 = QTable()
    tab3["p amp"] = 100.0 * tab2[gvals]["amplitude"] / tab[gvals]["photmjsr"]
    print(hstack([tab[gvals], tab2[gvals], tab3]))
