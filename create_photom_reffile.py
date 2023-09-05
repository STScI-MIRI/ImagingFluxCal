import numpy as np
import datetime
from astropy.table import QTable
from astropy.time import Time

from jwst.datamodels import MirImgPhotomModel

if __name__ == "__main__":

    # fmt: off
    image_filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
                     "F1500W", "F1800W", "F2100W", "F2550W"]
    image_filters = list(np.flip(image_filters))

    coron_filters = ["F1065C", "F1140C", "F1550C", "F2300C"]

    filters = image_filters + coron_filters
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

    data_list = {}

    print("filter, nfac, comfac, oldfac / comfac")
    for cfilter in filters:
        if cfilter in ["F2550W", "F1065C", "F1140C", "F1550C", "F2300C"]:
            rstr = "_bkgsub"
        else:
            rstr = ""

        # repeatability measurements for time dependence
        ntab_repeat = QTable.read(f"CalFacs/miri_calfactors{rstr}_repeat_{cfilter}_fit.dat",
                                  format="ascii.commented_header")
        amp = ntab_repeat[f"fit_exp_amp_{cfilter}"][0]
        tau = ntab_repeat[f"fit_exp_tau_{cfilter}"][0]
        if cfilter not in coron_filters:  # amplitude is reported in percentage, not absolute
            c = ntab_repeat[f"fit_exp_const_{cfilter}"][0]
            per_amp = amp / c
        else:
            per_amp = amp

        # average of all stars after correcting for time dependence
        ntab = QTable.read(f"CalFacs/miri_calfactors{rstr}_timecor_{cfilter}_ave.dat",
                           format="ascii.commented_header")
        cfac_ave = ntab[f"avecalfac_{cfilter}"][0]
        cfac_unc = ntab[f"avecalfac_unc_{cfilter}"][0]
        cfac_std = ntab[f"avecalfac_std_{cfilter}"][0]

        # account the sensitivity loss to the startday
        perfac = (per_amp * np.exp(days/tau)) + 1.0
        ncfacs = (perfac / (per_amp + 1)) * cfac_ave

        # calculated the value for the first 100 days
        #  approximates Commissioning so we can compare to the previous value
        #  not used otherwise
        new_cfactor = np.average(ncfacs[comvals])

        pipe_cfactor = cftab["photmjsr"][cftab["filter"] == cfilter.split("_")[0]][0]

        print(cfilter, cfac_ave, new_cfactor, pipe_cfactor / new_cfactor)

        # use the time dependent factors for F2550W to define the ranges for the multiple
        # photom reference files
        if cfilter == "F2550W":
            val_allowed = 0.01

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

        for k in range(n_photom):
            gvals = (days >= begday[k]) & (days < begday[k+1])
            cfac = np.average(ncfacs[gvals])

            if cfilter == "F2550W":
                data_list[f"{k}"] = []

            for csub in subarray_values:
                data_list[f"{k}"].append(
                    (cfilter, csub, cfac, cfac_unc)
                )

    for k in range(n_photom):
        # create the photom reference file
        data = np.array(
            data_list[f"{k}"],
            dtype=[
                ("filter", "S12"),
                ("subarray", "S15"),
                ("photmjsr", "<f4"),
                ("uncertainty", "<f4")
            ],
        )

        new_model = MirImgPhotomModel(phot_table=data)
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
        new_model.meta.pedigree = "INFLIGHT 2022-05-21 2023-09-10"
        tm = Time(begday[k] + startday, format="mjd")
        new_model.meta.useafter = tm.to_value(format="fits")
        new_model.meta.description = "Photom reference file."
        entry = "The flux calibration factors calculated from exponential fits to the"
        new_model.history.append(entry)
        entry = "time dependent flux calibration factors.  "
        new_model.history.append(entry)
        new_model.save(f"Photom/jwst_miri_photom_flight_{k+1}.fits")
