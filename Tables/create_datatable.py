import os
import numpy as np
from astropy.table import QTable, Table
import astropy.units as u

import cdspyreadme


if __name__ == "__main__":

    pfile = "Photom/jwst_miri_photom_flight_2jul24.fits"
    cftab = QTable.read(pfile, hdu=1)
    cftab_time = QTable.read(pfile, hdu=2)
    pixarea = 2.8606325654256e-13

    # fmt: off
    filters = ["F560W", "F770W", "F1000W",
               "F1130W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W",
               "F1065C", "F1140C", "F1550C", "F2300C"]
    # fmt: on

    dirs = ["HotStars", "ADwarfs", "SolarAnalogs"]

    colformats = {
        "name": "s",
        "time": ".1f",
        "pixrate": ".2e",
        "pixwelldepth": ".2e",
        "inttime": ".2f",
        "iflux": ".2e",
        "ifluxunc": ".2e",
        "ibkg": ".2e",
        "flux": ".2f",
        "fluxunc": ".2f",
        "bkg": ".2e",
    }

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
                        -1.0 * (cline["timemid"].value - pipe_t0) / pipe_tau
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

    # sort each name by filter
    fvals = {"F560W": 5.6, 
             "F770W": 7.7,
             "F1000W" : 10.0,
             "F1130W": 11.3,
             "F1280W": 12.8,
             "FND": 13.0,
             "F1500W": 15.0,
             "F1800W": 18.0,
             "F2100W": 21.0,
             "F2550W": 25.5,
             "F1065C": 106.5,
             "F1140C": 114.0,
             "F1550C": 155.0,
             "F2300C": 230.0}
    unames = np.unique(otab["name"])
    for cname in unames:
        gvals = otab["name"] == cname
        targfilters = otab["filter"][gvals].data
        svals = [fvals[cfilter] for cfilter in targfilters]
        sindxs = np.argsort(svals)
        otab[gvals] = otab[gvals][sindxs]

    otab.write("miri_absflux_program_data.dat", format="ipac", overwrite=True)
    # otab.write('miri_absflux_program_data.dat', format='ascii.mrt', overwrite=True, formats=colformats)

    # write latex table to provide example lines for paper
    otab.write(
        "miri_absflux_program_data.tex",
        overwrite=True,
        format="aastex",
        formats=colformats,
    )

    # write MRT table
    tablemaker = cdspyreadme.CDSTablesMaker()
    tablemaker.addTable(Table(otab), name="table1")

    # Customize ReadMe output
    tablemaker.title = (
        "MIRI Imaging and Coronagraphic Absolute Flux Calibration Measurements"
    )
    tablemaker.author = "Karl D. Gordon"
    tablemaker.date = 2020
    tablemaker.abstract = "This is my abstract..."
    tablemaker.more_description = "Additional information of the data context."
    #tablemaker.putRef("II/246", "2mass catalogue")
    #tablemaker.putRef("http://...", "external link")

    tablemaker.writeCDSTables()
    # Save ReadMe into a file
    with open("ReadMe", "w") as fd:
        tablemaker.makeReadMe(out=fd)
