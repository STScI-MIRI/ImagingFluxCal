import numpy as np
import datetime
from astropy.table import QTable

from jwst.datamodels import MirImgApcorrModel


if __name__ == "__main__":

    files = {}
    star = "BD+60 1753"

    dsets = {"F560W": ["1", "2", "3", "4"],
             "F770W": ["1", "2", "3", "4", "5", "6", "7", "8", "10", "11", "12", "13", "14"],
             "F1000W": ["1", "2", "3", "4"],
             "F1130W": ["1", "2"],
             "F1280W": ["1", "2"],
             "F1500W": ["1", "2", "3", "4", "5", "6"],
             "F1800W": ["1", "2", "3"],
             "F2100W": ["1", "2", "3"],
             "F2550W": ["1", "2", "3"],
             }

    dfwhmfac = {"F560W": "25.0",
                "F770W": "25.0",
                "F1000W": "15.0",
                "F1130W": "10.0",
                "F1280W": "10.0",
                "F1500W": "10.0",
                "F1800W": "10.0",
                "F2100W": "5.0",
                "F2550W": "5.0",
                "F1065C": "15.0",
                "F1140C": "15.0",
                "F1550C": "15.0",
                "F2300C": "10.0",
               }

    for cfilter in dsets.keys():
        files[cfilter] = []
        for cset in dsets[cfilter]:
            files[cfilter].append(
                f"ADwarfs/{cfilter}/{star}_set{cset}/miri_{star}_set{cset}_stage3_asn_i2dfwhmfac{dfwhmfac[cfilter]}_apcor.dat"
            )
    filters = list(dsets.keys())

    # add in corongraphic obs
    csubarray = {"F1065C": "MASK1065",
                 "F1140C": "MASK1140",
                 "F1550C": "MASK1550",
                 "F2300C": "MASKLYOT"}
    for cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
        files[cfilter] = [f"ADwarfs/{cfilter}/del UMi_set1/miri_del UMi_set1_stage3_asn_i2dfwhmfac{dfwhmfac[cfilter]}_apcor.dat",
                          f"ADwarfs/{cfilter}/HD 2811_set1/miri_HD 2811_set1_stage3_asn_i2dfwhmfac{dfwhmfac[cfilter]}_apcor.dat",
                          f"SolarAnalogs/{cfilter}/HD 167060_set1/miri_HD 167060_set1_stage3_asn_i2dfwhmfac{dfwhmfac[cfilter]}_apcor.dat"]
        filters.append(cfilter)

    fulltab = QTable(names=("filter", "subarray", "eefraction", "radius", "radius_unc", "apcor", "apcor_unc",
                            "skyin", "skyin_unc", "skyout", "skyout_unc"),
                     dtype=("str", "str", "f", "f", "f", "f", "f", "f", "f", "f", "f"))

    data_list = []
    for cfilter in filters:
        afiles = files[cfilter]

        nfiles = len(afiles)
        radiis = np.zeros((8, nfiles))
        apcors = np.zeros((7, nfiles))
        for k, cfile in enumerate(afiles):
            itab = QTable.read(cfile, format="ascii.commented_header")
            ees = itab["ee"]
            radiis[:, k] = itab["radii"]
            apcors[:, k] = itab["apcor"][0:-1]

        aradii = np.average(radiis, axis=1)
        aradii_std = np.std(radiis, axis=1)
        aapcor = np.average(apcors, axis=1)
        aapcor_std = np.std(apcors, axis=1)

        # print(cfilter)
        # print(ees)
        # print(aradii)
        # print(aradii_std)
        # print(aapcor)
        # print(aapcor_std)

        if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
            subarray_values = ["FULL", csubarray[cfilter]]
        else:
            subarray_values = ["FULL", "BRIGHTSKY", "SUB256", "SUB128", "SUB64"]

        for csub in subarray_values:
            for k, cee in enumerate(ees[:-1]):
                data_list.append(
                    (cfilter, csub, cee, aradii[k], aapcor[k], aradii[-2], aradii[-1])
                )
                fulltab.add_row([cfilter, csub, cee, aradii[k], aradii_std[k], aapcor[k], aapcor_std[k], 
                                 aradii[-2], aradii_std[-2], aradii[-1], aradii_std[-1]])

    data = np.array(
        data_list,
        dtype=[
            ("filter", "S6"),
            ("subarray", "S9"),
            ("eefraction", "<f4"),
            ("radius", "<f4"),
            ("apcorr", "<f4"),
            ("skyin", "<f4"),
            ("skyout", "<f4"),
        ],
    )

    new_model = MirImgApcorrModel(apcorr_table=data)
    d1 = datetime.datetime
    new_model.meta.date = d1.isoformat(d1.today())
    new_model.meta.filename = "jwst_miri_apcorr.fits"
    new_model.meta.telescope = "JWST"
    new_model.meta.instrument.name = "MIRI"
    new_model.meta.instrument.detector = "MIRIMAGE"
    new_model.meta.exposure.type = "MIR_IMAGE"
    new_model.meta.exposure.p_exptype = "MIR_IMAGE|MIR_4QPM|MIR_LYOT|MIR_TACQ|MIR_TACONFIRM|"
    new_model.meta.subarray = "GENERIC"
    new_model.meta.reftype = "APCORR"
    new_model.meta.author = "Karl Gordon"
    new_model.meta.pedigree = "INFLIGHT 2022-05-25 2023-05-25"
    new_model.meta.useafter = "2022-04-01T00:00:00"
    new_model.meta.description = "Apcorr reference file."
    entry = "The aperture corrections were computed starting from encircled-energy"
    new_model.history.append(entry)
    entry = "profiles measured with flight LVL-3 data and normalized to infinity"
    new_model.history.append(entry)
    entry = "using WebbPSF."
    new_model.history.append(entry)
    entry = " "
    new_model.history.append(entry)
    entry = "The inner and outer radii for the sky-background annulus are given by"
    new_model.history.append(entry)
    entry = "the SKYIN and SKYOUT columns."
    new_model.history.append(entry)
    new_model.save("ApCor/jwst_miri_apcorr_flight_27jul23.fits")

    print(fulltab)
    fulltab.write("ApCor/jwst_miri_apcorr_flight_27jul23_full.fits", overwrite=True)