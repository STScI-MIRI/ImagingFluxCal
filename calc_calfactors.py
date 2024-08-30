from os.path import exists
import argparse
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np

from statsmodels.stats.weightstats import DescrStatsW

from astropy.table import QTable
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.modeling import models, fitting
import astropy.units as u


def compute_stats(allfacs, weights, sigcut):
    "compute the weighted mean, weighted std, and weighted std of the mean"

    if len(allfacs) > 0:
        filtered_data = sigma_clip(allfacs, sigma=sigcut, maxiters=5)

        # compute the weighted mean
        newvals = allfacs[~filtered_data.mask]
        newweights = weights[~filtered_data.mask]
        meanval = np.average(newvals, weights=newweights)
        # compute weighted standard deviation
        if len(allfacs) > 3:
            meanstd = DescrStatsW(newvals, weights=newweights, ddof=1).std
            meanstdmean = meanstd / np.sqrt(np.sum(newweights))
        else:
            meanstd = 0.0
            meanstdmean = 0.0
        return (meanval, meanstd, meanstdmean, filtered_data, np.sum(newweights))
    else:
        return (None, None, None, None, None)


def get_calfactors(
    dir,
    filter,
    xaxisval="mflux",
    bkgsub=False,
    indivmos=False,
    indivcals=False,
    eefraction=0.7,
    repeat=False,
    subtrans=False,
    startday=59720.0,
    applytime=False,
    grieke=False,
):
    """
    Read in the observed and model fluxes and computer the calibration factors
    """
    if bkgsub:
        extstr = "_bkgsub"
    else:
        extstr = ""
    if indivmos:
        extstr = "_indivmos"
    if indivcals:
        extstr = f"{extstr}_indivcals"
    # read in observed fluxes
    obstab = QTable.read(f"{dir}/{filter}{extstr}_eefrac{eefraction}_phot.fits")
    # print(f"{dir}/{filter}{extstr}_eefrac{eefraction}_phot.fits")
    # read in model fluxes
    if grieke:
        mfilename = "Models/model_phot_grieke.fits"
    else:
        mfilename = "Models/model_phot.fits"
    modtab = QTable.read(mfilename)

    # get the info to remove the time dependent variation using the repeatability fit
    if applytime:
        if filter in ["F2550W", "F1065C", "F1140C", "F1550C", "F2300C"]:
            repstr = "_bkgsub"
        else:
            repstr = ""
        ffilename = f"CalFacs/miri_calfactors{repstr}_repeat_{filter}_fit.dat"
        ntab = QTable.read(ffilename, format="ascii.commented_header")
        # calculate the calibration factor versus time
        amp = ntab[f"fit_exp_amp_{filter}"][0]
        tau = ntab[f"fit_exp_tau_{filter}"][0]
        if filter not in ["F1065C", "F1140C", "F1550C", "F2300C", "FND"]:
            c = ntab[f"fit_exp_const_{filter}"][0]
            amp = amp / c  # put in percentage terms like the coronagraphs
        # print(filter, amp)

    if repeat:  # only use observations of BD+60 1753 and HD 2811
        # gvals = (obstab["name"] == "BD+60 1753") | (obstab["name"] == "HD 2811")
        gvals = obstab["name"] == "BD+60 1753"
        obstab = obstab[gvals]
    elif subtrans: 
        if filter == "F1280W":
            gvals = (obstab["name"] == "2MASS J18022716+6043356") & (
                abs(obstab["timemid"].value - 60177.3) < 1.0)
            obstab = obstab[gvals]
        else: # only use 2MASS J17571324+6703409 HJD around 59818.2693130625
            gvals = (obstab["name"] == "2MASS J17571324+6703409") & (
                abs(obstab["timemid"].value - 59818.3) < 1.0
            )
            obstab = obstab[gvals]

    names = []
    xvals = []
    cfactors = []
    cfactors_unc = []
    subarrs = []
    for k, cname in enumerate(obstab["name"]):
        if cname == "J1802271":
            cname = "2MASS J18022716+6043356"
        if cname == "GD153":
            cname = "GD 153"

        cfilter = obstab["filter"][k]
        oflux = obstab["aperture_sum_bkgsub"][k]
        oflux_unc = obstab["aperture_sum_bkgsub_err"][k]
        apcorr = obstab["apcorr"][k]
        pixarea = obstab["pixarea"][k]

        if applytime:  # fix the time dependancies
            ncfac = (amp * np.exp((obstab["timemid"][k].value - startday) / tau)) + 1.0
            # correct the sensitivity loss to the startday
            # oflux *= ncfac  #  / (amp + 1.)
            # correct the sensitivity loss to the value at infinite time
            oflux *= ncfac

        (mindx,) = np.where(modtab["name"] == cname)
        if len(mindx) < 1:
            print(f"Model fluxes for {cname} not present")
            exit()
        mflux = modtab[cfilter][mindx[0]]

        if xaxisval == "timemid":
            xval = obstab["timemid"][k]
        elif xaxisval == "rate":
            xval = obstab["pix_max"][k]
        elif xaxisval == "welldepth":
            xval = obstab["tgroup"][k] * obstab["ngroups"][k] * obstab["pix_max"][k]
        elif xaxisval == "bkg":
            xval = obstab["mean_bkg"][k]
        elif xaxisval == "inttime":
            xval = obstab["tgroup"][k] * obstab["ngroups"][k] * u.s
        elif xaxisval == "srctype":
            xval = dir
        elif xaxisval == "subarr":
            xval = obstab["subarray"][k]
        else:
            xval = mflux * 1e3

        if np.isfinite(oflux.value):
            cfactor = 1e-6 * mflux.value / (oflux.value * apcorr * pixarea.value)
            cfactor_unc = (oflux_unc / oflux) * cfactor
            # if obstab["name"][k] not in ["HD 167060", "16 Cyg B", "HD 37962", "del UMi"]:
            # print(obstab["name"][k], cfactor)
            names.append(obstab["name"][k])
            if isinstance(xval, u.Quantity):
                xvals.append(xval.value)
            else:
                xvals.append(xval)
            cfactors.append(cfactor)
            cfactors_unc.append(cfactor_unc)
            subarrs.append(obstab["subarray"][k])

    if xaxisval == "timemid":
        xvals = np.array(xvals) - startday

    cfactors = np.array(cfactors)
    cfactors_unc = np.array(cfactors_unc)
    xvals = np.array(xvals)
    subarrs = np.array(subarrs)
    names = np.array(names)

    # sort the subarrays
    if xaxisval == "subarr":
        sindxs = np.argsort(xvals)
        cfactors = cfactors[sindxs]
        cfactors_unc = cfactors_unc[sindxs]
        xvals = xvals[sindxs]
        subarrs = subarrs[sindxs]
        names = names[sindxs]

    res = (cfactors, cfactors_unc, xvals, subarrs, names)
    return res


def plot_calfactors(
    ax,
    filter,
    xaxisval,
    showleg=True,
    savefile=None,
    applysubarrcor=False,
    showcurval=True,
    bkgsub=False,
    indivmos=False,
    indivcals=False,
    eefraction=0.7,
    repeat=False,
    subtrans=False,
    applytime=False,
    grieke=False,
    shownames=False,
    x2ndaxis=True,
    notext=False,
    noignore=False,
    legonly=False,
    fitline=False,
):
    """
    Plot the calibration factors versus the requested xaxis.
    """
    dirs = ["HotStars", "ADwarfs", "SolarAnalogs"]
    if repeat or subtrans:
        dirs = ["ADwarfs"]
    pcols = ["b", "g", "r"]

    if legonly:
        psubsym = {
            "FULL": "o",
            "BRIGHTSKY": "s",
            "SUB256": "p",
            "SUB128": "P",
            "SUB64": "*",
            "MASK1065": "^",
            "MASK1140": ">",
            "MASK1550": "<",
            "MASKLYOT": "v",
        }
    elif "C" in filter:
        psubsym = {
            "MASK1065": "^",
            "MASK1140": ">",
            "MASK1550": "<",
            "MASKLYOT": "v",
        }
    elif filter == "FND":
        psubsym = {"FULL": "o"}
    else:
        psubsym = {
            "FULL": "o",
            "BRIGHTSKY": "s",
            "SUB256": "p",
            "SUB128": "P",
            "SUB64": "*",
        }

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
    }

    # ignore_names = ["HD 167060", "16 Cyg B", "HD 37962", "del UMi", "HD 106252", "HD 142331"]
    if noignore:
        ignore_names = []
    else:
        # ignore_names = ["GSPC P177-D", "16 Cyg B", "del UMi"]
        # ignore_names = ["HD 180609", "mu Col"]
        ignore_names = ["HD 180609"]
    # numbers in comments are from Bohlin et al. 2022 between IRAC 3.6/CALSPEC
    #  (is ratio "measured" from MIRI)
    # not used
    # modfac = {"HD 167060": 1.0/1.09,
    #           "16 Cyg B": 1.0/1.08,  # (0.93) 0.95 +/- 0.03
    #           "HD 37962": 1.0/1.06,  # (0.94) 0.96 +/- 0.01
    #           "del UMi": 1.0/1.03,  # (0.97) 0.98 +/- 0.01
    #           "HD 106252": 1.0/1.06, # (0.94) 0.98 +/- 0.02
    #           "HD 142331": 1.0/1.05,
    #           "BD+60 1753": 1.0}  # 0.99 +/- 0.02
    # modfac = {"HD 167060": 1.0,
    #           "16 Cyg B": 1.0,
    #           "HD 37962": 1.0,
    #           "del UMi": 1.0,
    #           "HD 180609": 1.0}

    startday = 59720.0

    # if xaxisval == "subarr":
    #     ax.scatter(
    #         psubsym.keys(),
    #         np.full(len(psubsym.keys()), 1.0),
    #         facecolors="none",
    #         edgecolors="none",
    #     )

    if bkgsub:
        extstr = "_bkgsub"
    else:
        extstr = ""

    # used for making srctype and subarr plots have better x distributions
    srctype_vals = {"HotStars": 0.0,
                    "ADwarfs": 1.0,
                    "SolarAnalogs": 2.0}
    subarr_vals = {"FULL": 0.0,
                    "BRIGHTSKY": 1.0,
                    "SUB256": 2.0,
                    "SUB128": 3.0,
                    "SUB64": 4.0,
                    "MASK1065": 5.0,
                    "MASK1140": 6.0,
                    "MASK1550": 7.0,
                    "MASKLYOT": 8.0}
    mu, sigma = 0, 0.12 # mean and standard deviation
    rng = np.random.default_rng()

    allfacs = []
    allfacuncs = []
    allnames = []
    allsubarr = []
    xvals = []
    meanfull = None
    for k, dir in enumerate(dirs):
        # print(f"{dir}/{filter}{extstr}_eefrac{eefraction}_phot.fits")
        if exists(f"{dir}/{filter}{extstr}_eefrac{eefraction}_phot.fits"):
            cfacs = get_calfactors(
                dir,
                filter,
                xaxisval=xaxisval,
                bkgsub=bkgsub,
                indivmos=indivmos,
                indivcals=indivcals,
                eefraction=eefraction,
                repeat=repeat,
                subtrans=subtrans,
                startday=startday,
                applytime=applytime,
                grieke=grieke,
            )
            # allfacs.append(cfacs[0])
            allfacuncs.append(cfacs[1])
            xvals.append(cfacs[2])
            allnames.append(cfacs[4])
            allsubarr.append(cfacs[3])
            for cfactor, cfactor_unc, xval, subarray, cname in zip(
                cfacs[0], cfacs[1], cfacs[2], cfacs[3], cfacs[4]
            ):
                # print(cname, xval, cfactor)
                # if not applysubarrcor:
                #     ax.errorbar(
                #         [xval],
                #         [cfactor],
                #         yerr=[cfactor_unc],
                #         fmt=f"{pcols[k]}{psubsym[subarray]}",
                #         alpha=0.1,
                #         markersize=10,
                #     )
                if applysubarrcor:
                    cfactor = cfactor * subarr_cor[subarray]
                if xaxisval == "welldepth":  # set the maximum well depth to the saturation level
                    xval = min([xval, 60000.0])

                if xaxisval == "srctype":
                    pxval = srctype_vals[xval] + rng.normal(mu, sigma)
                elif xaxisval == "subarr":
                    pxval = subarr_vals[xval] + rng.normal(mu, sigma)
                else:
                    pxval = xval

                allfacs.append(cfactor)
                ax.errorbar(
                    [pxval],
                    [cfactor],
                    yerr=[cfactor_unc],
                    fmt=f"{pcols[k]}{psubsym[subarray]}",
                    alpha=0.5,
                    markersize=10,
                )
                if shownames:
                    ax.text(pxval, cfactor, cname, rotation=45.0)
                # plot a red circle around those not used in the average
                if (cname in ignore_names) or (
                    (cname == "BD+60 1753")
                    and (xaxisval == "timmid")
                    and (abs(pxval - 60070.0) < 20.0)
                ):
                    ax.scatter(
                        [pxval],
                        [cfactor],
                        s=150,
                        facecolor="none",
                        edgecolor="m",
                    )
                    # print(modfac[cname])
                    # ax.scatter(
                    #    [pxval], [cfactor * modfac[cname]], s=200, facecolor="k", edgecolor="m",
                    # )
                if subarray == "FULL":
                    meanfull = cfactor
            # special code to give the differneces between the subarrays
            # if not applysubarrcor:
            #     print(cfacs[3])
            #     if meanfull is not None:
            #         print(cfacs[0] / meanfull)
            # exit()
    # allfacs = np.concatenate(allfacs)
    allfacs = np.array(allfacs)
    allfacuncs = np.concatenate(allfacuncs)
    allnames = np.concatenate(allnames)
    allsubarr = np.concatenate(allsubarr)
    xvals = np.concatenate(xvals)

    # compute the number of times each star has been observed
    # usful for reducing the weight of often observed stars so they do not dominate the average
    uniq_names = np.unique(allnames, return_counts=True)
    weights = np.zeros(len(xvals))
    for cname, ccount in zip(uniq_names[0], uniq_names[1]):
        weights[cname == allnames] = 1.0 / ccount

    # print the top 4 calibration factors with names
    # aindxs = np.flip(np.argsort(allfacs))
    # print(allfacs[atindxs[0:4]])
    # print(allnames[aindxs[0:4]])

    # ignore the 4 "high" stars
    gvals = []
    for k, cname in enumerate(allnames):
        if (cname in ignore_names) or (
            (cname == "BD+60 1753")
            and (xaxisval == "timmid")
            and (abs(xval - 60070.0) < 20.0)
        ):
            gvals.append(False)
        else:
            gvals.append(True)

    # use sigma clipping to remove the extreme outliers
    if filter == "FND":
        sigcut = 3.5
    else:
        sigcut = 3.5

    meanval, meanstd, meanstdmean, filtered_data, npts = compute_stats(
        allfacs[gvals], weights[gvals], sigcut
    )

    ax.axhline(y=meanval, color="k", linestyle="-", alpha=0.5)
    ax.axhline(y=meanval + meanstd, color="k", linestyle=":", alpha=0.5)
    ax.axhline(y=meanval - meanstd, color="k", linestyle=":", alpha=0.5)
    ax.scatter(
        (xvals[gvals])[filtered_data.mask],
        (allfacs[gvals])[filtered_data.mask],
        s=200,
        facecolor="none",
        edgecolor="m",
    )
    if np.sum(filtered_data.mask) > 0:
        print(f"{filter} sigma-clipped names", (allnames[gvals])[filtered_data.mask])
    perstd = 100.0 * meanstd / meanval
    perstdmean = 100.0 * meanstdmean / meanval

    # compute averages in bins
    if (xaxisval == "subarr") and "SUB256" in psubsym.keys():
        txvals = xvals[gvals][~filtered_data.mask]
        tallfacs = allfacs[gvals][~filtered_data.mask]
        tweights = weights[gvals][~filtered_data.mask]

        gvals2 = txvals == "SUB256"
        refres = compute_stats(tallfacs[gvals2], tweights[gvals2], sigcut)
        outvals = np.zeros((len(psubsym.keys()), 3))
        for k, csub in enumerate(psubsym.keys()):
            gvals2 = txvals == csub
            res = compute_stats(tallfacs[gvals2], tweights[gvals2], sigcut)
            if res[0] is not None:
                print(csub, res[0] / refres[0])
                outvals[k, :] = res[0:3]
        if savefile:
            otab = QTable()
            otab["name"] =  list(psubsym.keys())
            otab["calfacs"] = outvals[:, 0]
            otab["calfacs_unc"] = outvals[:, 1]
            otab["calfacs_uncmean"] = outvals[:, 2]
            otab.write(
                savefile.replace(".fits", "_subarr.dat"),
                overwrite=True,
                format="ascii.commented_header",
            )

    # compute averages in bins
    if xaxisval == "srctype":
        txvals = xvals[gvals][~filtered_data.mask]
        tallfacs = allfacs[gvals][~filtered_data.mask]
        tweights = weights[gvals][~filtered_data.mask]

        gvals2 = txvals == "ADwarfs"
        refres = compute_stats(tallfacs[gvals2], tweights[gvals2], sigcut)
        outvals = np.zeros((len(dirs), 3))
        for k, csub in enumerate(dirs):
            gvals2 = txvals == csub
            res = compute_stats(tallfacs[gvals2], tweights[gvals2], sigcut)
            if res[0] is not None:
                print(csub, res[0] / refres[0])
                outvals[k, :] = res[0:3]
        if savefile:
            otab = QTable()
            otab["name"] =  list(dirs)
            otab["calfacs"] = outvals[:, 0]
            otab["calfacs_unc"] = outvals[:, 1]
            otab["calfacs_uncmean"] = outvals[:, 2]
            otab.write(
                savefile.replace(".fits", "_srctype.dat"),
                overwrite=True,
                format="ascii.commented_header",
            )

    if (xaxisval == "timemid") and (not applytime):
        fit = fitting.LevMarLSQFitter()
        mod_init = models.Exponential1D(tau=-200.0, amplitude=-0.2) + models.Const1D(
            amplitude=0.70
        )
        mod_init[0].amplitude.bounds = [None, 0.0]
        mod_init[0].tau.fixed = True
        # mod_init[1].amplitude.fixed = True
        fitx = (xvals[gvals])[~filtered_data.mask]
        fity = (allfacs[gvals])[~filtered_data.mask]
        sindxs = np.argsort(fitx)
        mod_fit = fit(mod_init, fitx[sindxs], fity[sindxs])

        per_dev = (mod_fit(fitx) - fity) / mod_fit(fitx)
        per_dev = 100.0 * np.sqrt(np.sum(np.square(per_dev) / (len(fitx) - 2)))

        mod_dev = mod_fit(fitx) - fity
        mod_dev = np.sqrt(np.sum(np.square(mod_dev) / (len(fitx) - 2)))

        pxvals = np.arange(min(fitx), max(fitx))
        ax.plot(pxvals, mod_fit(pxvals), "m-")

        per_amp = 100.0 * (mod_fit[0].amplitude.value / mod_fit[1].amplitude.value)
        ax.text(
            0.95,
            0.02,
            rf"fit: {mod_fit[0].amplitude.value:.3f} exp(x/{mod_fit[0].tau.value:.1f}) + {mod_fit[1].amplitude.value:.3f} ({per_dev:.1f}%) [amp: {per_amp:.1f}%]",
            transform=ax.transAxes,
            ha="right",
        )

        # now see if we can derive the function to remove the trend
        # mod_div = (mod_fit[0].amplitude.value - mod_fit(pxvals))
        # print(mod_div)

    # fit a line
    if fitline:
        fit = fitting.LinearLSQFitter()
        mod_init = models.Linear1D()
        fitx = (xvals[gvals])[~filtered_data.mask]
        fity = (allfacs[gvals])[~filtered_data.mask]
        fitweights = weights[gvals][~filtered_data.mask]

        sindxs = np.argsort(fitx)
        mod_fit = fit(
            mod_init, fitx[sindxs], fity[sindxs], weights=1.0 / fitweights[sindxs]
        )
        pxvals = np.arange(min(fitx), max(fitx))
        ax.plot(pxvals, mod_fit(pxvals), "m-")

    if not notext:
        ax.text(
            0.95,
            0.08,
            rf"average: {meanval:.3f} $\pm$ {meanstd:.3f} ({perstd:.1f}% / {perstdmean:.1f}%)",
            transform=ax.transAxes,
            ha="right",
        )

    if savefile is not None:
        atab = QTable()
        atab[f"avecalfac_{filter}"] = [meanval]
        atab[f"avecalfac_unc_{filter}"] = [meanstdmean]
        atab[f"avecalfac_std_{filter}"] = [meanstd]
        atab[f"avecalfac_npts_{filter}"] = [npts]
        if (xaxisval == "timemid") and (not applytime):
            atab[f"fit_exp_amp_{filter}"] = [mod_fit[0].amplitude.value]
            atab[f"fit_exp_tau_{filter}"] = [mod_fit[0].tau.value]
            atab[f"fit_exp_const_{filter}"] = [mod_fit[1].amplitude.value]
            atab[f"fit_exp_std_{filter}"] = [mod_dev]
            sext = "_fit.dat"
        else:
            sext = "_ave.dat"
        atab.write(
            savefile.replace(".fits", sext),
            format="ascii.commented_header",
            overwrite=True,
        )

        otab = QTable()
        otab["name"] = allnames
        otab[f"calfac_{filter}"] = allfacs
        otab[f"calfac_{filter}_mean_dev"] = allfacs / meanval
        otab[f"calfac_{filter}_med_dev"] = allfacs / meanstd
        otab[f"modflux_{filter}"] = xvals
        otab.write(savefile, overwrite=True)

    # get the current pipeline calibration factor and plot
    if showcurval:
        # cftab = QTable.read("CalFactors/jwst_miri_photom_0079.fits")
        cftab = QTable.read("CalFactors/jwst_miri_photom_0201.fits", hdu=1)
        pipe_cfactor = cftab["photmjsr"][cftab["filter"] == filter.split("_")[0]][0]
        ax.axhline(y=pipe_cfactor, color="b", linestyle="--", alpha=0.5)

    # now make the plot nice
    if xaxisval == "timemid":
        ax.set_xlabel(f"Time [MJD] - {startday}")
    elif xaxisval == "rate":
        ax.set_xscale("log")
        ax.set_xlabel("Central Pixel Rate [DN/s]")
    elif xaxisval == "welldepth":
        ax.set_xlabel("Central Pixel Well Depth [DN]")
    elif xaxisval == "bkg":
        ax.set_xscale("log")
        ax.set_xlabel("Background [DN/s]")
    elif xaxisval == "inttime":
        ax.set_xscale("log")
        ax.set_xlabel("Integration Time [s]")
    elif xaxisval == "srctype":
        ax.set_xticks(list(srctype_vals.values()))
        ax.set_xticklabels(srctype_vals.keys())
        ax.tick_params(axis="x", labelrotation=60)
    elif xaxisval == "subarr":
        ax.set_xticks(list(subarr_vals.values()))
        ax.set_xticklabels(subarr_vals.keys())
        ax.tick_params(axis="x", labelrotation=60)
        if "C" in filter:
            ax.set_xlim(4.5, 8.5)
        else:
            ax.set_xlim(-0.5, 4.5)
    else:
        ax.set_xscale("log")
        ax.set_xlabel("Model Flux Density [mJy]")

    ax.set_ylabel("C [(MJy/sr) / (DN/s/pix)]")
    ax.set_title(f"{filter} / EEFRAC {eefraction}")

    ax.set_ylim(0.90 * meanval, 1.10 * meanval)

    def val2per(val):
        return (val / meanval) * 100.0 - 100.0

    def per2val(per):
        return ((per + 100) / 100.0) * meanval

    if x2ndaxis:
        secax = ax.secondary_yaxis("right", functions=(val2per, per2val))
        secax.set_ylabel("percentage", rotation=270)

    if showleg:

        # make space for the legend
        ylim = ax.get_ylim()
        ax.set_ylim(ylim[0], ylim[1] + 0.4 * (ylim[1] - ylim[0]))

        first_legend = [
            Patch(facecolor=ccol, edgecolor=ccol, label=cdir, alpha=0.5)
            for cdir, ccol in zip(dirs, pcols)
        ]
        first_legend.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label="Not in average",
                markerfacecolor="none",
                markeredgecolor="m",
                markersize=13,
                alpha=0.5,
            )
        )
        if legonly:
            loc = "upper right"
        else:
            loc = "upper center"
        leg1 = ax.legend(handles=first_legend, loc=loc)
        ax.add_artist(leg1)

        second_legend = []
        for ckey in psubsym.keys():
            if ckey[0:4] != "xMASK":
                second_legend.append(
                    Line2D(
                        [0],
                        [0],
                        marker=psubsym[ckey],
                        color="w",
                        label=ckey,
                        markerfacecolor="k",
                        markersize=10,
                        alpha=0.5,
                    )
                )
        if legonly:
            ncol = 2
        else:
            ncol = 1
        ax.legend(handles=second_legend, fontsize=9, ncol=ncol, loc="upper left")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter",
        help="filter to process",
        default="F770W",
        # fmt: off
        choices=["F560W", "F770W", "F1000W", "F1130W", "F1280W",
                 "F1500W", "F1800W", "F2100W", "F2550W",
                 "F1065C", "F1140C", "F1550C", "F2300C",
                 "FND"]
        # fmt: on
    )
    parser.add_argument(
        "--bkgsub",
        help="use results from background subtraction run",
        action="store_true",
    )
    parser.add_argument(
        "--indivmos",
        help="use results from individual mosaics (1 per cal image) instead of combined mosaics",
        action="store_true",
    )
    parser.add_argument(
        "--indivcals",
        help="use results from individual cal images instead of combined mosaics",
        action="store_true",
    )
    parser.add_argument(
        "--eefrac",
        default=0.7,
        help="Enclosed energy fraction to use",
        type=float,
    )
    parser.add_argument(
        "--xaxisval",
        help="x-axis values",
        default="mflux",
        choices=[
            "mflux",
            "timemid",
            "rate",
            "welldepth",
            "bkg",
            "inttime",
            "srctype",
            "subarr",
        ],
    )
    parser.add_argument(
        "--subarrcor",
        help="Apply subarray correction factors",
        action="store_true",
    )
    parser.add_argument(
        "--nocurval",
        help="do not plot the current calfactor",
        action="store_true",
    )
    parser.add_argument(
        "--noignore",
        help="do not ignore any stars",
        action="store_true",
    )
    parser.add_argument(
        "--repeat",
        help="plot the repeatability observations",
        action="store_true",
    )
    parser.add_argument(
        "--subtrans",
        help="plot the subarray transfer observations",
        action="store_true",
    )
    parser.add_argument(
        "--applytime",
        help="remove the time dependent variation",
        action="store_true",
    )
    parser.add_argument(
        "--grieke",
        help="use GRieke models for the 10 G-stars, CALSPEC models for the rest",
        action="store_true",
    )
    parser.add_argument(
        "--shownames",
        help="show the names for each point",
        action="store_true",
    )
    parser.add_argument(
        "--detmulti",
        help="4 panel plot of detector characteristics",
        action="store_true",
    )
    parser.add_argument(
        "--sourcemulti",
        help="2 panel plot of source characteristics",
        action="store_true",
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # make plot
    fontsize = 14
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    if args.bkgsub:
        extstr = "_bkgsub"
    else:
        extstr = ""
    if args.indivmos:
        extstr = "_indivmos"
    if args.indivcals:
        extstr = "_indivcals"
    if args.repeat:
        extstr = f"{extstr}_repeat"
    if args.grieke:
        extstr = f"{extstr}_grieke"
    if args.subarrcor:
        extstr = f"{extstr}_subarracor"
    if args.applytime:
        extstr = f"{extstr}_timecor"

    savefacs = f"CalFacs/miri_calfactors{extstr}_{args.filter}.fits"
    if args.detmulti:
        fontsize = 10
        font = {"size": fontsize}
        plt.rc("font", **font)

        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
        plot_calfactors(
            ax[0, 0],
            args.filter,
            "rate",
            showleg=True,
            applysubarrcor=args.subarrcor,
            showcurval=(not args.nocurval),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans,
            applytime=args.applytime,
            grieke=args.grieke,
            noignore=args.noignore,
        )
        plot_calfactors(
            ax[0, 1],
            args.filter,
            "inttime",
            savefile=savefacs,
            showleg=False,
            applysubarrcor=args.subarrcor,
            showcurval=(not args.nocurval),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans,
            applytime=args.applytime,
            grieke=args.grieke,
            noignore=args.noignore,
        )
        plot_calfactors(
            ax[1, 0],
            args.filter,
            "welldepth",
            showleg=False,
            applysubarrcor=args.subarrcor,
            showcurval=(not args.nocurval),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans,
            applytime=args.applytime,
            grieke=args.grieke,
            noignore=args.noignore,
        )
        plot_calfactors(
            ax[1, 1],
            args.filter,
            "subarr",
            showleg=False,
            applysubarrcor=args.subarrcor,
            showcurval=(not args.nocurval),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans,
            applytime=args.applytime,
            grieke=args.grieke,
            noignore=args.noignore,
        )
        fname = f"miri_calfactors_{args.filter}_detmulti"
    elif args.sourcemulti:
        fontsize = 10
        font = {"size": fontsize}
        plt.rc("font", **font)

        fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
        plot_calfactors(
            ax[0,0],
            args.filter,
            "mflux",
            showleg=True,
            applysubarrcor=args.subarrcor,
            showcurval=(not args.nocurval),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans,
            applytime=args.applytime,
            grieke=args.grieke,
            noignore=args.noignore,
        )
        plot_calfactors(
            ax[0,1],
            args.filter,
            "bkg",
            showleg=False,
            applysubarrcor=args.subarrcor,
            showcurval=(not args.nocurval),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans,
            applytime=args.applytime,
            grieke=args.grieke,
            noignore=args.noignore,
        )
        plot_calfactors(
            ax[1,0],
            args.filter,
            "timemid",
            showleg=False,
            applysubarrcor=args.subarrcor,
            showcurval=(not args.nocurval),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans,
            applytime=args.applytime,
            grieke=args.grieke,
            noignore=args.noignore,
        )
        plot_calfactors(
            ax[1,1],
            args.filter,
            "srctype",
            showleg=False,
            applysubarrcor=args.subarrcor,
            showcurval=(not args.nocurval),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans,
            applytime=args.applytime,
            grieke=args.grieke,
            noignore=args.noignore,
        )
        fname = f"miri_calfactors_{args.filter}_sourcemulti"
    else:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
        plot_calfactors(
            ax,
            args.filter,
            args.xaxisval,
            savefile=savefacs,
            applysubarrcor=args.subarrcor,
            showcurval=(not args.nocurval),
            bkgsub=args.bkgsub,
            indivmos=args.indivmos,
            indivcals=args.indivcals,
            eefraction=args.eefrac,
            repeat=args.repeat,
            subtrans=args.subtrans,
            applytime=args.applytime,
            grieke=args.grieke,
            shownames=args.shownames,
            noignore=args.noignore,
        )
        fname = f"miri_calfactors_{args.filter}_{args.xaxisval}"

    plt.tight_layout()

    fname = f"{fname}_eefrac{args.eefrac}"
    if args.bkgsub:
        fname = f"{fname}_bkgsub"
    if args.indivmos:
        fname = f"{fname}_indivmos"
    if args.indivcals:
        fname = f"{fname}_indivcals"
    if args.repeat:
        fname = f"{fname}_repeat"
    if args.subtrans:
        fname = f"{fname}_subtrans"
    if args.grieke:
        fname = f"{fname}_grieke"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
