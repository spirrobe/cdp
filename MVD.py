#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: spirrobe -> github.com/spirrobe/
"""

def MVD(inputdata,
        data_is_lwc=True,
        data_is_conc=False,
        binsizes=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20,
                  22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50],
        windspeed=[1],  # in m/s
        samplearea=0.298,  # as area in square millimeters
        samplefrequency=0.1,
        ):
    import numpy as np
    # if the data is a concentration it cannot be a lwc
    if data_is_conc:
        data_is_lwc = False

    inputdata = inputdata.copy()

    if data_is_lwc:
        # data is as we expect it usually, and we can go on
        _lwc = inputdata
    elif data_is_conc:
        # data is given as concentration, not liquid water content, hence
        # we have to first calculate the lwc
        _lwc = LWC(inputdata,
                    windspeed=windspeed,
                    samplearea=samplearea,
                    samplefrequency=samplefrequency,
                    data_is_conc=data_is_conc,
                    binsizes=binsizes,
                    combined=False,
                    )
    else:
        # data_is_lwc and data_is_conc are set to false and we have to
        # calculate everything from scratch, starting with concentration
        _conc = ConcPerCCM(inputdata,
                            windspeed=windspeed,
                            samplearea=samplearea,
                            samplefrequency=samplefrequency,
                            combined=False,
                            )
        # calculate liquid water concentration from concentration
        _lwc = LWC(_conc,
                    samplearea=samplearea,
                    samplefrequency=samplefrequency,
                    data_is_conc=data_is_conc,
                    binsizes=binsizes,
                    combined=False)

    _lwc = np.nan_to_num(_lwc)

    if len(_lwc.shape) == 1:
        _lwc = _lwc[np.newaxis, :]

    # to be used as lookup table
    binsizes = np.asarray(binsizes)

    #adjust for when there is no lwc, and the first value is the index found
    _lwc = np.ma.masked_array(_lwc, mask=(_lwc <= 0))

    #divide the array by the respective sum to achieve a value in percent
    pro = (_lwc.T / np.nansum(_lwc, axis=1)).T

    #create an array with the time-wise sum of LWC
    mvd = np.nancumsum(pro, axis=1)

    mvd_check = mvd.copy()
    #change values so nanargmin works in finding the first one that is above 0.5
    mvd_check[mvd_check < 0.5] = 1

    # gives us the indices where the first bin is above 1
    mvd_check = np.nanargmin(mvd_check, axis=1)

    ix = np.arange(mvd.shape[0])

    # follows the formula in PADS Vers. 3.6.3
    # b_i + (b_i+1 - b_i) * (0.5 - cum_i-1) / pro_i
    mvd = binsizes[mvd_check] + (binsizes[mvd_check+1] - binsizes[mvd_check]) * (0.5 - mvd[ix, mvd_check-1])/ pro[ix, mvd_check]
    mvd[mvd.mask] = 0
    mvd = mvd.data


    # # #oldway to do calculation, way slower as loop goes over all records
    # # initiate empty array with the second dimensions
    # mvd = np.zeros(_lwc.shape[0])

    # # # divide the single bins by the sum of the bins
    # for i in range(_lwc.shape[0]):

    #     # optimization: if sum is zero, no droplets are found  -> skip record
    #     if np.nansum(_lwc[i, :]) == 0:
    #         continue

    #     # calculate the terms
    #     _lwc[i, :] /= np.nansum(_lwc[i, :])
    #     _cumsum = np.nancumsum(_lwc[i, :])

    #     for j in range(_cumsum.shape[0]):
    #         if _cumsum[j] > 0.5:
    #             mvd[i] = ((0.5 - _cumsum[j-1])/_lwc[i, j])
    #             mvd[i] *= (binsizes[j+1]-binsizes[j])
    #             mvd[i] += binsizes[j]
    #             # optimization: break if we found it
    #             break



    return mvd
