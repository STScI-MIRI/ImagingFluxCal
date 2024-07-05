from astropy.table import QTable

if __name__ == "__main__":

    filter_pairs = {"F2300C": (["F2100W", "F2550W"], [20.80, 25.23], 22.66),
                    "F1550C": (["F1500W", "F1800W"], [15.05, 17.97], 15.51),
                    "F1140C": (["F1000W", "F1130W"], [9.95, 11.31], 11.30),
                    "F1065C": (["F1000W", "F1130W"], [9.95, 11.31], 10.60)}
    for cfilter in filter_pairs.keys():
        fpair = filter_pairs[cfilter]
        ifilter = fpair[0]
        iwave = fpair[1]
        cwave = fpair[2]

        nfile1 = f"CalFacs/miri_calfactors_repeat_{ifilter[0]}_fit.dat"
        nfile2 = f"CalFacs/miri_calfactors_repeat_{ifilter[1]}_fit.dat"
        if ifilter[1] == "F2550W":
            nfile2 = nfile2.replace("_repeat", "_bkgsub_repeat")
        ntab1 = QTable.read(nfile1, format="ascii.commented_header")
        ntab2 = QTable.read(nfile2, format="ascii.commented_header")

        # calculate the calibration factor versus time
        amp1 = ntab1[f"fit_exp_amp_{ifilter[0]}"][0]
        tau1 = ntab1[f"fit_exp_tau_{ifilter[0]}"][0]
        c1 = ntab1[f"fit_exp_const_{ifilter[0]}"][0]
        startday1 = ntab1[f"fit_exp_startday_{ifilter[0]}"][0]

        amp2 = ntab2[f"fit_exp_amp_{ifilter[1]}"][0]
        tau2 = ntab2[f"fit_exp_tau_{ifilter[1]}"][0]
        c2 = ntab2[f"fit_exp_const_{ifilter[1]}"][0]
        startday2 = ntab2[f"fit_exp_startday_{ifilter[1]}"][0]

        weight = (cwave - iwave[0])/(iwave[1] - iwave[0])
        
        amp = (amp1 / c1) * (1 - weight) + (amp2 / c2) * weight
        tau = tau1 * (1 - weight) + tau2 * weight
        startday = startday1 * (1 - weight) + startday2 * weight
        print(cfilter, amp, tau, startday)

        otab = QTable()
        otab[f"fit_exp_amp_{cfilter}"] = [amp]
        otab[f"fit_exp_tau_{cfilter}"] = [tau]
        otab[f"fit_exp_startday_{cfilter}"] = [startday]
        otab.write(f"CalFacs/miri_calfactors_bkgsub_repeat_{cfilter}_fit.dat",
                   format="ascii.commented_header",
                   overwrite=True)