import os
import numpy as np
from astropy.table import QTable
import astropy.units as u


if __name__ == "__main__":

    pfile = "Photom/jwst_miri_photom_flight_2jul24.fits"
    cftab = QTable.read(pfile, hdu=1)
    cftab_time = QTable.read(pfile, hdu=2)
    pixarea = 2.8606325654256e-13

    # fmt: off
    filters = ["F2100W", "F560W", "F770W", "F1000W",
               "F1130W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W",
               "F1065C", "F1140C", "F1550C", "F2300C"]
    # fmt: on

    dirs = ["HotStars", "ADwarfs", "SolarAnalogs"]

    colformats = {'time': '.1f', 
                  "pixrate": ".2e",
                  "pixwelldepth": ".2e",
                  "inttime": ".2f",
                  "iflux": ".2e",
                  "ifluxunc": ".2e",
                  "ibkg": ".2e",
                  "flux": ".2f",
                  "fluxunc": ".2f",
                  "bkg": ".2e"}

    otab = QTable(
        names=(
            "name",
            "PID",
            "srctype",
            "filter",
            "subarray",
            "time",
            "pixrate",
            "pixwelldepth",
            "inttime",
            "iflux",
            "ifluxunc",
            "ibkg",
            "flux",
            "fluxunc",
            "bkg",
        ),
        units=(
            "",
            "",
            "",
            "",
            "",
            u.d,
            u.DN / u.s,
            u.DN,
            u.s,
            u.DN / u.s,
            u.DN / u.s,
            u.DN / u.s,
            u.mJy,
            u.mJy,
            u.MJy / u.sr,
        ),
        dtype=(
            "str",
            "int",
            "str",
            "str",
            "str",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
            "float32",
        ),
    )

    for cdir in dirs:
        for cfilter in filters:

            fname = f"{cdir}/{cfilter}_eefrac0.7_phot.fits"
            if os.path.isfile(fname):

                tab = QTable.read(fname)
                # print(tab.colnames)
                # exit()

                for cline in tab:

                    # compute calibrated flux
                    gval = (cftab["filter"] == cfilter) & (
                        cftab["subarray"] == cline["subarray"]
                    )
                    pipe_cfactor = cftab["photmjsr"][gval][0]
                    pipe_amp = cftab_time["amplitude"][gval][0]
                    pipe_tau = cftab_time["tau"][gval][0]
                    pipe_t0 = cftab_time["t0"][gval][0]
                    cur_cf = pipe_cfactor + pipe_amp * np.exp(
                        (cline["timemid"].value - pipe_t0) / pipe_tau
                    )
                    phys_flux = (
                        1e9
                        * cline["aperture_sum_bkgsub"].value
                        * cline["apcorr"]
                        * pixarea
                        * cur_cf
                    ) * u.mJy
                    phys_flux_unc = (
                        1e9
                        * cline["aperture_sum_bkgsub_err"].value
                        * cline["apcorr"]
                        * pixarea
                        * cur_cf
                    ) * u.mJy
                    phys_bkg = (cline["mean_bkg"].value * cur_cf) * u.MJy / u.sr

                    otab.add_row(
                        (
                            cline["name"],
                            cline["program"],
                            cdir,
                            cline["filter"],
                            cline["subarray"],
                            cline["timemid"],
                            cline["pix_max"],
                            cline["tgroup"] * u.s * cline["ngroups"] * cline["pix_max"],
                            cline["tgroup"] * cline["ngroups"] * u.s,
                            cline["aperture_sum_bkgsub"],
                            cline["aperture_sum_bkgsub_err"],
                            cline["mean_bkg"],
                            phys_flux,
                            phys_flux_unc,
                            phys_bkg,
                        )
                    )

    # sort by name
    sindxs = np.argsort(otab["name"])
    otab = otab[sindxs]

    #otab.write("miri_absflux_program_data.dat", format="ipac", overwrite=True)
    # otab.write('ecliptic_cols.dat', format='ascii.mrt', overwrite=True, formats=colformats)

    otab.write("miri_absflux_program_data.tex", overwrite=True, format="aastex",
               formats=colformats)
