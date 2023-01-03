from aper_one_filter import aper_one_filter


if __name__ == '__main__':

    types = ["HotStars", "ADwarfs", "SolarAnalogs"]
    filters = ["F560W", "F770W"]

    for ctype in types:
        for cfilter in filters:
            aper_one_filter(ctype, cfilter)
