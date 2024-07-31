import numpy as np
import datetime
from astropy.table import QTable
from astropy.time import Time

from jwst.datamodels import MirImgPhotomModel

if __name__ == "__main__":

    # fmt: off
    image_filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W",
                     "F1500W", "F1800W", "F2100W", "F2550W"]
    # image_filters = list(np.flip(image_filters))

    coron_filters = ["F1065C", "F1140C", "F1550C", "F2300C", "FND"]

    filters = image_filters + coron_filters
    # fmt: on
    csubarray = {"F1065C": "MASK1065",
                 "F1140C": "MASK1140",
                 "F1550C": "MASK1550",
                 "F2300C": "MASKLYOT"}
    fndsubarray = ["FULL", "MASK1065", "MASK1140", "MASK1550", "MASKLYOT", "SLITLESSPRISM"]

    # based on calibration factor ratios and dedicated subarray transfer observations
    subarr_cor = {
        "FULL": 1.0,
        "BRIGHTSKY": 1.005, 
        "SUB256": 0.98,
        "SUB128": 1.00, 
        "SUB64": 0.966,
        "MASK1065": 1.0,
        "MASK1140": 1.0,
        "MASK1550": 1.0,
        "MASKLYOT": 1.0,
        "SLITLESSPRISM": 1.0,
    }

    # current calibration factors
    cftab = QTable.read("CalFactors/jwst_miri_photom_0079.fits")

    startday = 59720.
    days = np.arange(0.0, 1000.0, 1.0)
    comvals = days < 50.

    data_list = []
    data_list_time = []

    fulltab = QTable(names=("filter", "amplitude", "tau", "photmjysr", "startday", "uncertainty"),
                     dtype=("str", "f", "f", "f", "f", "f"))

    # print("filter, nfac, comfac, oldfac / comfac, comfac / oldfac")
    print("filter   CF      amp     amp_per   tau   CF_unc  CF_unc_per repeat_per")
    for cfilter in filters:
        if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
            rstr = "_bkgsub"
        else:
            rstr = ""
        if cfilter == "F2550W":
            rstr2 = "_bkgsub"
        else:
            rstr2 = rstr

        # repeatability measurements for time dependence
        ntab_repeat = QTable.read(f"CalFacs/miri_calfactors{rstr2}_repeat_{cfilter}_fit.dat",
                                  format="ascii.commented_header")
        amp = ntab_repeat[f"fit_exp_amp_{cfilter}"][0]
        tau = ntab_repeat[f"fit_exp_tau_{cfilter}"][0]
        if cfilter not in coron_filters:  # amplitude is reported in percentage, not absolute
            c = ntab_repeat[f"fit_exp_const_{cfilter}"][0]
            per_amp = amp / c
        else:
            per_amp = amp
        # repeatability as a percentage for paper table
        if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C", "FND"]:
            repeat_per = 0.0
        else:
            repeat_per = ntab_repeat[f"fit_exp_std_per_{cfilter}"][0]

        # average of all stars after correcting for time dependence
        ntab = QTable.read(f"CalFacs/miri_calfactors{rstr}_grieke_subarracor_timecor_{cfilter}_ave.dat",
                           format="ascii.commented_header")
        cfac_ave = ntab[f"avecalfac_{cfilter}"][0]
        cfac_unc = ntab[f"avecalfac_unc_{cfilter}"][0]
        cfac_std = ntab[f"avecalfac_std_{cfilter}"][0]
        cfac_unc_per = 100.0 * cfac_unc / cfac_ave
        cfac_npts = ntab[f"avecalfac_npts_{cfilter}"][0]

        # account the sensitivity loss to the startday
        #perfac = (per_amp * np.exp(days/tau)) + 1.0
        #ncfacs = (perfac / (per_amp + 1)) * cfac_ave

        amp = (per_amp / (per_amp + 1)) * cfac_ave
        const = (1.0 / (per_amp + 1)) * cfac_ave
        # ncfacs = (amp * np.exp(days/tau)) + const

        fulltab.add_row([cfilter, amp, -1.*tau, const, startday, cfac_std])

        frac_change = (const + amp) / const
        amp_per = (np.absolute(amp) / const) * 100.0

        # calculated the value for the first 100 days
        #  approximates Commissioning so we can compare to the previous value
        #  not used otherwise
        # new_cfactor = np.average(ncfacs[comvals])
        # pipe_cfactor = cftab["photmjsr"][cftab["filter"] == cfilter.split("_")[0]][0]
        
        # build the data structure needed
        # allowed subarrays
        if cfilter in ["F1065C", "F1140C", "F1550C", "F2300C"]:
            subarray_values = ["FULL", csubarray[cfilter]]
        elif cfilter in ["FND"]:
            subarray_values = fndsubarray
        else:
            subarray_values = ["FULL", "BRIGHTSKY", "SUB256", "SUB128", "SUB64"]

        print(cfilter, cfac_unc)

        print(f"{cfilter} & {cfac_ave:.4f} & {amp:.4f} & {amp_per:.1f} & {-1.*tau:.1f} & {cfac_unc:.5f} & {cfac_unc_per:.2f} & {cfac_npts:.2f} & {repeat_per:.2f} \\\\ ")

        for csub in subarray_values:
            data_list.append((cfilter, csub, cfac_ave / subarr_cor[csub], cfac_unc / subarr_cor[csub]))
            data_list_time.append((amp / subarr_cor[csub], tau, startday))

    # temp fix for FND - remove once photom file includes this filter
    #cfilter = "FND"
    #csub = "FULL"
    #data_list.append((cfilter, csub, 1.0, 0.1))
    #data_list_time.append((0.0, -200., startday))

    # save time dependent coefficients
    fulltab.write("CalFacs/jwst_miri_photom_coeff.dat", format="ipac", overwrite=True)

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

    # create the photom reference file
    data_time = np.array(
        data_list_time,
        dtype=[
            ("amplitude", "<f4"),
            ("tau", "<f4"),
            ("t0", "<f4")
        ],
    )

    new_model = MirImgPhotomModel(phot_table=data, timecoeff=data_time)
    d1 = datetime.datetime
    new_model.meta.date = d1.isoformat(d1.today())
    new_model.meta.filename = f"jwst_miri_photom_31jul24.fits"
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
    new_model.meta.pedigree = "INFLIGHT 2022-05-21 2024-07-02"
    new_model.meta.useafter = "2022-04-01T00:00:00"
    new_model.meta.description = "Photom reference file."
    entry = "The flux calibration factors calculated from exponential fits to the"
    new_model.history.append(entry)
    entry = "time dependent flux calibration factors.  "
    new_model.history.append(entry)
    new_model.save(f"Photom/jwst_miri_photom_flight_31jul24.fits")
