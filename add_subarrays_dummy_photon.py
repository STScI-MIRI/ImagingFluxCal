import numpy as np
import datetime
from jwst.datamodels import MirImgPhotomModel, open
from astropy.table import QTable


if __name__ == "__main__":
    cur_model = open("Photom/jwst_miri_photom_flight_8oct25.fits")
    cur_phot = cur_model.phot_table
    cur_exp = cur_model.timecoeff_exponential
    cur_line = cur_model.timecoeff_linear

    data_list = []
    data_list_time_exp = []
    data_list_time_linear = []
    for cline, cline2, cline3 in zip(cur_phot, cur_exp, cur_line):
        data_list.append((cline[0], cline[1], cline[2], cline[3]))
        data_list_time_exp.append((cline2[0], cline2[1], cline2[2], cline2[3], cline2[4], cline2[5]))
        data_list_time_linear.append((cline3[0], cline3[1], cline3[2], cline3[3]))
        if cline["subarray"] == "SUB128":
            data_list.append((cline[0], "SUB128_IP", cline[2], cline[3]))
            data_list_time_exp.append((cline2[0], "SUB128_IP", cline2[2], cline2[3], cline2[4], cline2[5]))
            data_list_time_linear.append((cline3[0], "SUB128_IP", cline3[2], cline3[3]))
        if cline["subarray"] == "SUB64":
            data_list.append((cline[0], "SUB64_IP", cline[2], cline[3]))
            data_list_time_exp.append((cline2[0], "SUB64_IP", cline2[2], cline2[3], cline2[4], cline2[5]))
            data_list_time_linear.append((cline3[0], "SUB64_IP", cline3[2], cline3[3]))

    # create the photom reference file
    data = np.array(
        data_list,
        dtype=[
            ("filter", "S12"),
            ("subarray", "S15"),
            ("photmjsr", "<f4"),
            ("uncertainty", "<f4")
        ],
    )

    # create the exponentital time dependence
    data_time_exp = np.array(
        data_list_time_exp,
        dtype=[
                ("filter", "S12"),
                ("subarray", "S15"),
                ("t0", "<f4"),
                ("amplitude", "<f4"),
                ("tau", "<f4"),
                ("const", "<f4"),
        ],
    )

    # create the linear time dependence
    data_time_linear = np.array(
        data_list_time_linear,
        dtype=[
                ("filter", "S12"),
                ("subarray", "S15"),
                ("t0", "<f4"),
                ("lossperyear", "<f4"),
        ],
    )

    print(data)
    print("---")
    print(data_time_exp)
    print("----")
    print(data_time_linear)

    new_model = MirImgPhotomModel(phot_table=data,
                                  timecoeff_exponential=data_time_exp,
                                  timecoeff_linear=data_time_linear)
    d1 = datetime.datetime
    new_model.meta.date = d1.isoformat(d1.today())
    new_model.meta.filename = f"jwst_miri_photom_27may26.fits"
    new_model.meta.telescope = "JWST"
    new_model.meta.instrument.name = "MIRI"
    new_model.meta.instrument.detector = "MIRIMAGE"
    new_model.meta.exposure.type = "MIR_IMAGE"
    new_model.meta.photometry.pixelarea_steradians = 2.8606325654256E-13
    new_model.meta.photometry.pixelarea_arcsecsq = 0.01217199
    new_model.meta.instrument.band = "N/A"
    new_model.meta.exposure.p_exptype = "MIR_IMAGE|MIR_4QPM|MIR_LYOT|MIR_TACQ|MIR_TACONFIRM|MIR_CORONCAL|"
    new_model.meta.subarray = "GENERIC"
    new_model.meta.reftype = "PHOTOM"
    new_model.meta.author = "Karl Gordon"
    # updates to next 2 lines needed
    new_model.meta.pedigree = "INFLIGHT 2022-05-21 2025-08-27"
    new_model.meta.useafter = "2022-04-01T00:00:00"
    new_model.meta.description = "Photom reference file."
    entry = "The flux calibration factors calculated from exponential + linear"
    new_model.history.append(entry)
    entry = "fits to the time dependent flux calibration factors.  "
    new_model.history.append(entry)
    new_model.save("Photom/jwst_miri_photom_flight_27may26.fits")

    ttab = QTable.read("Photom/jwst_miri_photom_flight_27may26.fits", hdu=3)
    print(ttab)
