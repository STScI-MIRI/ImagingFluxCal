import os
import glob
import argparse

import numpy as np
from astropy.io import fits

from jwst.associations import asn_from_list
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from helpers.miri_helpers import miri_detector1, miri_image2, miri_image3


def sort_uncal(subdir, filter):
    """
    Sort through the files in a directory into sets per star.  Checks that
    all are for the same filter.
    """
    rfiles = glob.glob(f"{subdir}/{filter}/*_uncal.fits")
    filenames = []
    objects = []
    filters = []
    obssets = []
    for cfile in rfiles:
        filenames.append(cfile)
        chdr = fits.getheader(cfile, 0)
        ctargname = chdr["TARGNAME"]
        objects.append(ctargname)
        cfilter = chdr["FILTER"]
        filters.append(cfilter)
        obssets.append((cfile.split("/")[2]).split("_")[0])
        if cfilter != filter.split("_")[0]:
            print("filter does not match")
            print(f"expected: {filter}, actual: {cfilter}")
            exit()
    filenames = np.array(filenames)
    objects = np.array(objects)
    filters = np.array(filters)
    obssets = np.array(obssets)

    # get the filenames for each unique object/obsset combination
    uobj, rind = np.unique(objects, return_inverse=True)
    setobj = {}
    for k, cobj in enumerate(uobj):
        gvals = cobj == objects
        # create a seperate set for each unique set of observations
        usets, sind = np.unique(obssets[gvals], return_inverse=True)
        for m, cset in enumerate(usets):
            gvals = (cobj == objects) & (cset == obssets)
            setobj[f"{cobj}_set{m+1}"] = filenames[gvals]

    return setobj


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--filter",
        help="filter to process",
        default="F770W",
        # fmt: off
        choices=["F560W", "F770W", "F770W_subarray", "F770W_repeat", "F1000W",
                 "F1130W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W",
                 "F1065C", "F1140C", "F1550C", "F2300C"],
        # fmt: on
    )
    parser.add_argument(
        "--dir",
        choices=["HotStars", "ADwarfs", "SolarAnalogs"],
        default="ADwarfs",
        help="directory to process",
    )
    args = parser.parse_args()

    # get the files for each object for the specified filter
    objsets = sort_uncal(args.dir, args.filter)

    # for each object, run the pipeline
    for ckey in objsets.keys():
        print(f"pipelining {ckey}")

        # directory for pipeline results
        ndir = f"{args.dir}/{args.filter}/{ckey}"
        if not os.path.exists(ndir):
            print(f"making dir {ndir}")
            os.mkdir(ndir)

        # calwebb_detector1
        cbase = f"{ndir}/{ckey}_stage1"
        with open(f"{cbase}.cfg", "w") as f:
            f.write("[*]\n")
            f.write(f"handler = file:{cbase}.log\n")
            f.write("level = INFO\n")

        print(f"detector1 for {ckey}")
        miri_detector1(objsets[ckey], ndir, logfile=f"{cbase}.cfg")

        # calwebb_image2
        cbase = f"{ndir}/{ckey}_stage2"
        with open(f"{cbase}.cfg", "w") as f:
            f.write("[*]\n")
            f.write(f"handler = file:{cbase}.log\n")
            f.write("level = INFO\n")

        ratefiles = glob.glob(f"{ndir}/*_rate.fits")
        print(f"image2 for {ckey}")
        miri_image2(ratefiles, ndir, logfile=f"{cbase}.cfg")

        # calwebb_image3
        # generate association file
        calfiles = glob.glob(f"{ndir}/*_cal.fits")
        miri_asn_name = f'miri_{ckey}_stage3_asn'
        miri_asn = asn_from_list.asn_from_list(calfiles, rule=DMS_Level3_Base,
                                               product_name=miri_asn_name)
        miri_asn_file = f'{miri_asn_name}.json'
        with open(miri_asn_file, 'w') as outfile:
            name, serialized = miri_asn.dump(format='json')
            outfile.write(serialized)

        # log file setup
        cbase = f"{ndir}/{ckey}_stage3"
        with open(f"{cbase}.cfg", "w") as f:
            f.write("[*]\n")
            f.write(f"handler = file:{cbase}.log\n")
            f.write("level = INFO\n")

        print(f"image3 for {ckey}")
        miri_image3(miri_asn_file, ndir, logfile=f"{cbase}.cfg")
