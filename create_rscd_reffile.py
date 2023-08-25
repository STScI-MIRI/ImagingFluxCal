from jwst.datamodels import RSCDModel

if __name__ == "__main__":
    rscd_model = RSCDModel("RSCD/jwst_miri_rscd_0017.fits")

    print('before')
    print(rscd_model.rscd_group_skip_table)

    for k, tabdata in enumerate(rscd_model.rscd_group_skip_table):
        subarray_table = tabdata['subarray']
        readpatt_table = tabdata['readpatt']
        group_skip_table = tabdata['group_skip']
        if readpatt_table == "FASTR1":
            if subarray_table == "BRIGHTSKY":
                rscd_model.rscd_group_skip_table[k] = ("BRIGHTSKY", "FASTR1", 4)
            elif subarray_table == "SUB256":
                rscd_model.rscd_group_skip_table[k] = ("SUB256", "FASTR1", 12)
            elif subarray_table == "SUB128":
                rscd_model.rscd_group_skip_table[k] = ("SUB128", "FASTR1", 16)
            elif subarray_table == "SUB64":
                rscd_model.rscd_group_skip_table[k] = ("SUB64", "FASTR1", 18)

    print('after')
    print(rscd_model.rscd_group_skip_table)

    rscd_model.save("RSCD/jwst_miri_rscd_kdg_updated.fits")