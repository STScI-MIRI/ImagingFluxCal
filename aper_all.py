from os.path import exists

from aper_one_filter import aper_one_filter


if __name__ == '__main__':

    types = ["HotStars", "ADwarfs", "SolarAnalogs"]
    filters = ["F2300C", "F1550C", "F1140C", "F1065C",
               "F560W", "F770W", "F1000W", "F1130W", "F1280W", "F1500W",
               "F1800W", "F2100W", "F2550W"]
    # filters = ["F770W"]

    for ctype in types:
        for cfilter in filters:
            if exists(f"{ctype}/{cfilter}/"):
                print(f"{ctype}/{cfilter}/")
                aper_one_filter(ctype, cfilter)
                aper_one_filter(ctype, cfilter, bkgsub=True)

#    aper_one_filter("ADwarfs", "F770W_subarray")
#    aper_one_filter("ADwarfs", "F770W_repeat")
