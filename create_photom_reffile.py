import numpy as np
import datetime
from astropy.table import QTable

from jwst.datamodels import MirImgPhotomModel

if __name__ == "__main__":

    # fmt: off
    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
               "F1500W", "F1800W", "F2100W", "F2550W"]

    coron_filters = ["F1065C", "F1140C", "F1550C", "F2300C"]
    # fmt: on
    csubarray = {"F1065C": "MASK1065",
                 "F1140C": "MASK1140",
                 "F1550C": "MASK1550",
                 "F2300C": "MASKLYOT"}

    # current calibration factors
    cftab = QTable.read("CalFactors/jwst_miri_photom_0079.fits")

    startday = 59720.
    days = np.arange(0.0, 1000.0, 1.0)
    comvals = days < 100.

    print("filter, nfac, comfac, oldfac / comfac")
    for cfilter in np.flip(filters):
        if cfilter in ["F2550W", "F1065C", "F1140C", "F1550C", "F2300C"]:
            rstr = "_bkgsub"
        else:
            rstr = ""
        ntab = QTable.read(f"CalFacs/miri_calfactors{rstr}_{cfilter}_fit.dat", format="ascii.commented_header")
        # calculate the calibration factor versus time
        amp = ntab[f"fit_exp_amp_{cfilter}"][0]
        tau = ntab[f"fit_exp_tau_{cfilter}"][0]
        c = ntab[f"fit_exp_const_{cfilter}"][0]
        unc = ntab[f"fit_exp_std_{cfilter}"][0]
        ncfacs = (amp * np.exp(days/tau)) + c

        # calculated the value for the first 100 days
        #  approximates Commissioning so we can compare to the previous value
        #  not used otherwise
        new_cfactor = np.average(ncfacs[comvals])

        pipe_cfactor = cftab["photmjsr"][cftab["filter"] == cfilter.split("_")[0]][0]

        print(cfilter, c, new_cfactor, pipe_cfactor / new_cfactor)

        # use the time dependent factors for F2550W to define the ranges for the multiple
        # photom reference files
        if cfilter == "F2550W":
            val_allowed = 0.05

            tval = ncfacs[0]
            begday = []
            begval = []
            begday.append(days[0])
            begval.append(tval)
            for cval, cday in zip(ncfacs, days):
                if ((cval - tval) / tval) > val_allowed:
                    print(begday[-1], cday, begval[-1], cval)
                    begday.append(cday)
                    begval.append(cval)
                    tval = cval
            begday.append(cday)
            begval.append(cval)

            n_photom = len(begday) - 1
            print(f"# photom files needed: {n_photom}")

        # allowed subarrays
        if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
            subarray_values = ["FULL", csubarray[cfilter]]
        else:
            subarray_values = ["FULL", "BRIGHTSKY", "SUB256", "SUB128", "SUB64"]

        data_list = []
        for k in range(n_photom):
            gvals = (days >= begday[k]) & (days < begday[k+1])
            cfac = np.average(ncfacs[gvals])

            for csub in subarray_values:
                data_list.append(
                    (cfilter, csub, cfac, unc)
                )

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

            new_model = MirImgPhotomModel(apcorr_table=data)
            d1 = datetime.datetime
            new_model.meta.date = d1.isoformat(d1.today())
            new_model.meta.filename = f"jwst_miri_photom_{k+1}.fits"
            new_model.meta.telescope = "JWST"
            new_model.meta.instrument.name = "MIRI"
            new_model.meta.instrument.detector = "MIRIMAGE"
            new_model.meta.exposure.type = "MIR_IMAGE"
            new_model.meta.exposure.p_exptype = "MIR_IMAGE|MIR_4QPM|MIR_LYOT|MIR_TACQ|MIR_TACONFIRM|"
            new_model.meta.subarray = "GENERIC"
            new_model.meta.reftype = "PHOTOM"
            new_model.meta.author = "Karl Gordon"
            # updates to next 2 lines needed
            new_model.meta.pedigree = "INFLIGHT 2022-05-25 2023-05-25"
            new_model.meta.useafter = "2022-04-01T00:00:00"
            #####
            new_model.meta.description = "Photom reference file."
            entry = "The flux calibration factors calculated from exponential fits to the"
            new_model.history.append(entry)
            entry = "time dependent flux calibration factors.  "
            new_model.history.append(entry)
            new_model.save(f"Photom/jwst_miri_photom_flight_{k+1}.fits")
