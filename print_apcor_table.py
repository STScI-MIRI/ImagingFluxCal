import argparse

from astropy.table import QTable


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter",
        help="filter to process",
        default="F770W",
        # fmt: off
        choices=["F560W", "F770W", "F770W_subarray", "F770W_repeat", "F1000W",
                 "F1130W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W",
                 "F1065C", "F1140C", "F1550C", "F2300C"],
        # fmt: on
    )
    args = parser.parse_args()

    print("Old")
    tab = QTable.read("ApCor/jwst_miri_apcorr_0008.fits")
    print(tab[(tab["subarray"] == "FULL") & (tab["filter"] == args.filter)])

    print("New")
    tab = QTable.read("jwst_miri_apcorr_flight_8may8.fits")
    print(tab[(tab["subarray"] == "FULL") & (tab["filter"] == args.filter)])
