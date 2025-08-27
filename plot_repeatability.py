import argparse
import numpy as np

import matplotlib.pyplot as plt

from astropy.modeling import models, fitting
from astropy.modeling import Fittable1DModel, Parameter
from astropy.table import QTable
from astropy.time import Time

from calc_calfactors import get_calfactors


class PowerLaw1D_Shift(Fittable1DModel):
    """
    One dimensional power law model.

    Parameters
    ----------
    amplitude : float
        Model amplitude at the reference point
    x_0 : float
        Reference point
    alpha : float
        Power law index

    See Also
    --------
    BrokenPowerLaw1D, ExponentialCutoffPowerLaw1D, LogParabola1D

    Notes
    -----
    Model formula (with :math:`A` for ``amplitude`` and :math:`\\alpha` for ``alpha``):

        .. math:: f(x) = A (x / x_0) ^ {-\\alpha}

    """

    amplitude = Parameter(default=1, description="Peak value at the reference point")
    x_0 = Parameter(default=1, description="Reference point")
    alpha = Parameter(default=1, description="Power law index")

    @staticmethod
    def evaluate(x, amplitude, x_0, alpha):
        """One dimensional power law model function."""
        xx = x + x_0
        return amplitude * xx ** (-alpha)

    # @staticmethod
    # def fit_deriv(x, amplitude, x_0, alpha):
    #     """One dimensional power law derivative with respect to parameters."""
    #     xx = x / x_0

    #     d_amplitude = xx ** (-alpha)
    #     d_x_0 = amplitude * alpha * d_amplitude / x_0
    #     d_alpha = -amplitude * d_amplitude * np.log(xx)

    #     return [d_amplitude, d_x_0, d_alpha]

    @property
    def input_units(self):
        if self.x_0.input_unit is None:
            return None
        return {self.inputs[0]: self.x_0.input_unit}

    def _parameter_units_for_data_units(self, inputs_unit, outputs_unit):
        return {
            "x_0": inputs_unit[self.inputs[0]],
            "amplitude": outputs_unit[self.outputs[0]],
        }


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--show_prev", help="show previous time dependence", action="store_true"
    )
    parser.add_argument("--dexp", help="include exp+exp model", action="store_true")
    parser.add_argument("--docs", help="include only exp+line model", action="store_true")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filters = [
        "F560W",
        "F770W",
        "F1000W",
        "F1130W",
        "F1280W",
        "F1500W",
        "F1800W",
        "F2100W",
        "F2550W",
    ]
    # filters = np.flip(filters)

    # make plot
    fontsize = 16
    font = {"size": fontsize}
    plt.rc("font", **font)
    plt.rc("lines", linewidth=2)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("ytick.major", width=2)

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(16, 10))

    if args.show_prev:
        # cftab = QTable.read("CalFactors/jwst_miri_photom_0201.fits", hdu=1)
        # cftab_time = QTable.read("CalFactors/jwst_miri_photom_0201.fits", hdu=2)
        cftab = QTable.read("Photom/jwst_miri_photom_flight_30aug24.fits", hdu=1)
        cftab_time = QTable.read("Photom/jwst_miri_photom_flight_30aug24.fits", hdu=2)

    ax = axs[0]
    startday = 59720
    for k, cfilter in enumerate(filters):
        if cfilter == "F2550W":
            bkgsub = True
            rstr = "_bkgsub"
        else:
            bkgsub = False
            rstr = ""

        cfacs = get_calfactors(
            "ADwarfs",
            cfilter,
            xaxisval="timemid",
            repeat=True,
            startday=startday,
            bkgsub=bkgsub,
            grieke=True,
        )

        # remove point that is near 350 - excess stripping in images
        gvals = np.absolute(cfacs[2] - 350.0) < 5.0

        yvals = cfacs[0][~gvals]
        yvals_unc = cfacs[1][~gvals]
        xvals = cfacs[2][~gvals]

        # now get the two stars that we repeated twice to fill in the time gap
        for stype, sname in zip(["SolarAnalogs", "ADwarfs"], ["HD 37962", "del UMi"]):
            cfacs2 = get_calfactors(
                stype,
                cfilter,
                xaxisval="timemid",
                startday=startday,
                bkgsub=bkgsub,
                grieke=True,
            )
            nvals = [pname == sname for pname in cfacs2[4]]
            if np.sum(nvals) > 1:
                nxvals = cfacs2[2][nvals]
                nyvals = cfacs2[0][nvals]
                nyvals_unc = cfacs2[1][nvals]

                # find the value of BD+60 1753 that is closest to the last value
                # use this value to adjust the new star to be on the same scale
                #    later interpolate between values
                sindxs = np.argsort(np.absolute(nxvals[-1] - xvals))
                nyvals *= yvals[sindxs[0]] / nyvals[-1]

                yvals = np.append(yvals, nyvals)
                yvals_unc = np.append(yvals_unc, nyvals_unc)
                xvals = np.append(xvals, nxvals)

        # fit an exponential
        # ignore the bad data point for F770W
        gvals = abs(xvals - (60070.0 - startday)) > 20.0

        fit = fitting.LevMarLSQFitter()
        mod_init = models.Exponential1D(tau=-200.0, amplitude=-0.2) + models.Const1D(
            amplitude=0.70
        )
        # mod_init2 = models.PowerLaw1D(amplitude=0.70, x_0=1.0, alpha=-1.0)
        # mod_init = (models.Exponential1D(tau=-150., amplitude=-0.2)
        #             + models.Linear1D(intercept=0.70, slope=0.0))
        mod_init[0].amplitude.bounds = [None, 0.0]
        if cfilter in ["F560W", "F770W", "F1000W", "F1130W", "F1280W", "F1500W"]:
            mod_init[0].tau.fixed = True
        else:
            mod_init[0].tau.bounds = [-400.0, -100.0]

        mod_init2 = PowerLaw1D_Shift(amplitude=0.70, x_0=60.0, alpha=0.5)+ models.Const1D(
            amplitude=0.70
        )
        mod_init2[0].x_0.bounds = [50.0, 200.0]
        mod_init2[1].amplitude = 0.0
        mod_init2[1].amplitude.fixed = True

        mod_init3 = models.Exponential1D(tau=-200.0, amplitude=-0.2) + models.Linear1D(
            slope=-0.5, intercept=1.0
        )

        mod_init4 = (
            models.Exponential1D(tau=-100.0, amplitude=-0.2)
            + models.Exponential1D(tau=-500.0, amplitude=-0.1)
            + models.Const1D(amplitude=0.70)
        )
        # mod_init4[0].tau.bounds = (-150.0, -50.0)
        mod_init4[0].amplitude.bounds = (-0.2, 0.0)
        #mod_init4[0].tau = mod_fit3[0].tau.value
        #mod_init4[0].tau.fixed = True
        # mod_init4[1].tau.bounds = (-1000.0, -300.0)
        # mod_init4[1].amplitude.bounds = (-1.0, 0.0)

        fitx = xvals[gvals]
        fity = yvals[gvals]
        sindxs = np.argsort(fitx)

        # plot the data
        meanval = np.average(yvals)
        yvals = meanval / yvals
        sindxs = np.argsort(xvals)
        yoff0 = k * 0.25
        ydiff0 = np.average(yvals) - np.average(yvals[sindxs[-5:]])
        yoff = yoff0 + ydiff0
        yoff2 = k * 0.12
        if cfilter == "F2550W":
            lname = "data"
        else:
            lname = None
        ax.errorbar(
            xvals, yvals + yoff, yerr=yvals_unc, fmt="ko", alpha=0.5, label=lname
        )
        ax.plot([0.0, max(fitx)], [1.0 + yoff0, 1.0 + yoff0], "k:", alpha=0.5)
        axs[1].plot([0.0, max(fitx)], [0.0 + yoff2, 0.0 + yoff2], "k:", alpha=0.5)

        # print(
        #     "powerlaw amp/shift/alpha:",
        #     mod_fit2[0].amplitude.value,
        #     mod_fit2[0].x_0.value,
        #     mod_fit2[0].alpha.value,
        # )
        # print("exp tau/amp", mod_fit[0].tau.value, mod_fit[0].amplitude.value)
        # print("exp+line tau/amp", mod_fit3[0].tau.value, mod_fit3[0].amplitude.value)
        # if args.dexp:
        #     print("exp+exp taus", mod_fit4[0].tau.value, mod_fit4[1].tau.value)
        #     print(
        #         "exp+exp amps", mod_fit4[0].amplitude.value, mod_fit4[1].amplitude.value
        #     )

        modnames = ["exp", "powerlaw", "exp+line"]
        allmods = [mod_init, mod_init2, mod_init3]
        allnparam = [3, 3, 4]
        pcol = ["m", "r", "g"]
        if args.dexp:
            modnames += ["exp+exp"]
            allmods += [mod_init4]
            allnparam += [5]
            pcol += ["c"]

        sigtext = ""
        for cname, cmod, cparam, ccol in zip(modnames, allmods, allnparam, pcol):

            mod_fit = fit(cmod, fitx[sindxs], fity[sindxs])

            per_dev = (mod_fit(fitx) - fity) / mod_fit(fitx)
            per_dev = 100.0 * np.sqrt(np.sum(np.square(per_dev) / (len(fitx) - cparam)))

            mod_dev = mod_fit(fitx) - fity
            mod_dev = np.sqrt(np.sum(np.square(mod_dev) / (len(fitx) - cparam)))

            # save the fit results - update for exp+line 26 Aug 2025
            if cname == "exp+line":
                atab = QTable()
                atab[f"fit_exp_amp_{cfilter}"] = [mod_fit[0].amplitude.value]
                atab[f"fit_exp_tau_{cfilter}"] = [mod_fit[0].tau.value]
                atab[f"fit_exp_intercept_{cfilter}"] = [mod_fit[1].intercept.value]
                atab[f"fit_exp_slope_{cfilter}"] = [mod_fit[1].slope.value]
                atab[f"fit_exp_startday_{cfilter}"] = [startday]
                atab[f"fit_exp_std_{cfilter}"] = [mod_dev]
                atab[f"fit_exp_std_per_{cfilter}"] = [per_dev]
                sext = "_fit.dat"
                atab.write(
                    f"CalFacs/miri_calfactors{rstr}_repeat_{cfilter}_fit.dat",
                    format="ascii.commented_header",
                    overwrite=True,
                )

            pxvals = np.arange(0, max(fitx))
            modvals = np.average(mod_fit(xvals)) / mod_fit(pxvals)

            show_plot = False
            if args.docs:
                if cname == "exp+line":
                    show_plot = True
            else:
                show_plot = True

            if cfilter == "F2550W":
                lname = cname
            else:
                lname = None

            if show_plot:
                ax.plot(pxvals, modvals + yoff, f"{ccol}-", label=lname)
                modxvals = meanval / mod_fit(xvals)
                axs[1].errorbar(
                    xvals, (yvals - modxvals) + yoff2, yerr=yvals_unc, fmt=f"{ccol}o", alpha=0.5
                )
                sigtext = rf"{sigtext}$\sigma$({cname})={per_dev:.2f}% "

            # show the delta change
            if cname == "exp+line":
                bpredx = [min(fitx), max(fitx)]
                bvals = 1.0 / mod_fit(bpredx)
                per_decrease = 100.0 * (bvals[1] - bvals[0]) / bvals[0]

                shifty = 0.05
                ax.text(425.0, 1.0 + yoff + shifty, cfilter)
                ax.text(
                    0.0,
                    yoff + shifty + 0.03 + modvals[0],
                    rf"$\Delta$={-1.*per_decrease:.1f}%",
                    backgroundcolor="w",
                    fontsize=0.8 * fontsize,
                )

            if args.show_prev & (cname == "exp"):
                amp = cftab_time["amplitude"][cftab["filter"] == cfilter][0]
                tau = cftab_time["tau"][cftab["filter"] == cfilter][0]
                startday = cftab_time["t0"][cftab["filter"] == cfilter][0]
                mod_fit[0].amplitude = amp
                mod_fit[0].tau = -1.0 * tau
                modvals = meanval / mod_fit(pxvals)

                # get the offset to match the range used for the current photom up to 800 days 
                # a better visual comparison of the change between the delivered
                # and current time dependence
                tmodvals = meanval / mod_fit(xvals)
                tvals = xvals < 800.0
                extoff = np.average(yvals[tvals] - tmodvals[tvals])

                if cfilter == "F2550W":
                    tlab = "current photom"
                else:
                    tlab = None

                ax.plot(pxvals, modvals + yoff + extoff, "b:", label=tlab)

        # show percentage sigma text on diff plot
        shifty2 = 0.04
        axs[1].text(
            100.0,
            yoff2 + shifty2,
            sigtext,
            backgroundcolor="w",
            fontsize=0.6 * fontsize,
        )

        # # predict the throughput in 3 and 6 years
        # predx = 0.0 + np.array([0.0, 3 * 365, 6 * 365])
        # print(cfilter)
        # print(
        #     Time(predx + startday, format="mjd").to_value(format="iso", subfmt="date")
        # )
        # all_fits = [mod_fit, mod_fit2, mod_fit3]
        # all_names = ["     exp:", "powerlaw:", "exp+line:"]
        # if args.dexp:
        #     all_fits.append(mod_fit4)
        #     all_names.append(" exp+exp:")
        # for cfit, cfname in zip(all_fits, all_names):
        #     predy = meanval / cfit(predx)
        #     print(
        #         cfname,
        #         predy / predy[0],
        #         f"{100.0 * (predy[2] - predy[1]) / 3.0:.3f}%/year",
        #     )

    ax.legend(fontsize=0.7 * fontsize, ncol=2)

    ax.set_ylim(0.9, 3.6)
    ntvals = np.arange(0, max(fitx) + 50, 100)
    ax.set_xticks(ntvals)
    ax.set_xticklabels(
        Time(ntvals + startday, format="mjd").to_value(format="iso", subfmt="date")
    )
    ax.tick_params(axis="x", labelrotation=60)
    # ax.set_xlabel(f"Date")
    ax.set_ylabel("Fractional change (+ const)")

    axs[1].yaxis.tick_right()
    axs[1].yaxis.set_label_position("right")
    # axs[1].set_xlim(400., 800.)
    axs[1].set_ylim(-0.04, 1.25)
    # axs[1].set_xlabel(f"Date")

    axs[1].set_xticks(ntvals)
    axs[1].set_xticklabels(
        Time(ntvals + startday, format="mjd").to_value(format="iso", subfmt="date")
    )
    axs[1].tick_params(axis="x", labelrotation=60)

    secax = ax.secondary_xaxis("top")
    secax.set_xlabel(f"MJD Time - {startday} [days]")
    secax = axs[1].secondary_xaxis("top")
    secax.set_xlabel(f"MJD Time - {startday} [days]")

    axs[1].set_ylabel("Data - Model Residual (+ const)")
    # axs[1].set_ylabel("Fractional change (+ const)")

    plt.tight_layout()

    fname = "all_repeatability"
    if args.docs:
        fname = f"{fname}_docs"
    if args.png:
        fig.savefig(f"Figs/{fname}.png")
    elif args.pdf:
        fig.savefig(f"Figs/{fname}.pdf")
    else:
        plt.show()
