from astropy.table import QTable


if __name__ == "__main__":

    files = ["CalFactors/jwst_miri_photom_0218.fits",
             "Photom/jwst_miri_photom_flight_8oct25.fits",
            ]
    

    cftab1 = QTable.read(files[0], hdu=1)
    cftab2 = QTable.read(files[1], hdu=1)

    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
               "F1500W", "F1800W", "F2100W", "F2550W", "FND",
               "F1065C", "F1140C", "F1550C", "F2300C"]

    for cfilter in filters:
        gvals1 = cftab1["filter"] == cfilter
        subarrs = cftab1[gvals1]["subarray"]
        for csubarr in subarrs:
            gvals1 = (cftab1["filter"] == cfilter) & (cftab1["subarray"] == csubarr)
            gvals2 = (cftab2["filter"] == cfilter) & (cftab2["subarray"] == csubarr)
            fac1 = cftab1[gvals1]["photmjsr"].value[0]
            fac2 = cftab2[gvals2]["photmjsr"].value[0]
            print(fac1/fac2, cfilter, csubarr)
