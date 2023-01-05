from os.path import exists

from aper_one_filter import aper_one_filter


if __name__ == '__main__':

    types = ["HotStars", "ADwarfs", "SolarAnalogs"]
    filters = ["F560W", "F770W", "F1000W", "F1130W", "F1500W", "F1880W",
               "F2100W", "F2550W"]

    for ctype in types:
        for cfilter in filters:
            print(f"{ctype}/{cfilter}/")
            print(exists(f"{ctype}/{cfilter}/"))
            if exists(f"{ctype}/{cfilter}/"):
                aper_one_filter(ctype, cfilter)
