import glob
from astropy.io import fits
import numpy as np

if __name__ == '__main__':

    filters = ["F560W", "F770W", "F1000W",
               "F1130W", "F1280W", "F1500W", "F1800W", "F2100W", "F2550W", "FND",
               "F1065C", "F1140C", "F1550C", "F2300C"]

    types = ["HotStars", "ADwarfs", "SolarAnalogs"]
    # types = [types[0]]

    # get the info
    stype = []
    name = []
    filter = []
    pid = []
    subarray = []
    for ctype in types:
        for cfilter in filters:
            files = glob.glob(f"{ctype}/{cfilter}/*uncal.fits")
            for cfile in files:
                hdr = fits.getheader(cfile)
                stype.append(ctype)
                targname = hdr["TARGNAME"]
                if targname == "J1757132":
                    targname = "2MASS J17571324+6703409"
                if targname == "J1802271":
                    targname = "2MASS J18022716+6043356"
                if targname == "GD153":
                    targname == "GD 153"
                if targname == "PG 1057+719":
                    targname == "WD 1057+719"
                name.append(targname)
                filter.append(hdr["FILTER"])
                pid.append(hdr["PROGRAM"][1:])
                subarray.append(hdr["SUBARRAY"])

    stype = np.array(stype)
    name = np.array(name)
    filter = np.array(filter)
    pid = np.array(pid)
    subarray = np.array(subarray)

    # condense to a nice table format and print as the needed latex table data contents
    for ctype in types:
        print("\\hline \\multicolumn{4}{c}{" + ctype + "} \\\\ \\hline")
        gtype = ctype == stype
        ustars = np.unique(name[gtype])
        for cname in ustars:
            pname = cname
            gname = name[gtype] == cname
            upids = np.unique(pid[gtype][gname])
            for cpid in upids:
                ppid = cpid
                gpid = cpid == pid[gtype][gname]
                usubarrays = np.unique(subarray[gtype][gname][gpid])
                for csubarray in usubarrays:
                    gsubarrays = csubarray == subarray[gtype][gname][gpid]
                    ufilters = np.unique(filter[gtype][gname][gpid][gsubarrays])
                    for tfilter in ["F770W", "F560W"]:
                        if tfilter in ufilters:
                            ufilters = np.concatenate(([tfilter], ufilters[0:-1]))
                    print(f"{pname} & {ppid} & {csubarray} & {', '.join(ufilters)} \\\\")
                    pname = ""
                    ppid = ""
