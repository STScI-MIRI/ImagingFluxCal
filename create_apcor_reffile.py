import numpy as np
import datetime
from astropy.table import QTable

from jwst.datamodels import MirImgApcorrModel


if __name__ == "__main__":

    subarray_values = ["FULL", "BRIGHTSKY", "SUB256", "SUB128", "SUB64"]

    files = {}
    star = "BD+60 1753"
    dsets = {}
    dsets["F560W"] = ["1", "2", "3", "4"]
    dsets["F770W"] = ["1", "2", "3"]
    dsets["F1000W"] = ["1", "2", "3"]
    dsets["F1130W"] = ["1"]
    dsets["F1280W"] = ["1"]
    dsets["F1500W"] = ["1", "2", "3", "4"]
    dsets["F1800W"] = ["1", "2"]
    dsets["F2100W"] = ["1"]
    dsets["F2550W"] = ["1"]

    for cfilter in dsets.keys():
        files[cfilter] = []
        for cset in dsets[cfilter]:
            files[cfilter].append(
                f"ADwarfs/{cfilter}/{star}_set{cset}/miri_{star}_set{cset}_stage3_asn_i2d_apcor.dat"
            )

    # add in corongraphic obs
    files["F1140C"] = ["ADwarfs/F1140C/del UMi_set1/miri_del UMi_set1_stage3_asn_i2d_apcor.dat",
                       "ADwarfs/F1140C/HD 2811_set1/miri_HD 2811_set1_stage3_asn_i2d_apcor.dat",
                       "SolarAnalogs/F1140C/HD 167060_set1/miri_HD 167060_set1_stage3_asn_i2d_apcor.dat"]
    filters = list(dsets.keys())
    filters.append("F1140C")

    data_list = []
    for cfilter in filters:
        afiles = files[cfilter]
        # afiles = glob.glob(
        #     f"ADwarfs/{cfilter}/BD+60 1753_set?/miri_BD+60 1753_set?_stage3_asn_i2d_apcor.dat"
        # )

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

        for csub in subarray_values:
            for k, cee in enumerate(ees[:-1]):
                data_list.append(
                    (cfilter, csub, cee, aradii[k], aapcor[k], aradii[-2], aradii[-1])
                )

    print(data_list[0])

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
    new_model.meta.exposure.p_exptype = "MIR_IMAGE|MIR_TACQ|MIR_TACONFIRM|"
    new_model.meta.subarray = "GENERIC"
    new_model.meta.reftype = "APCORR"
    new_model.meta.author = "Karl Gordon"
    new_model.meta.pedigree = "INFLIGHT 2022-05-25 2023-05-08"
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
    new_model.save("jwst_miri_apcorr_flight_8may8.fits")
