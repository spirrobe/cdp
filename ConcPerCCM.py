#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: spirrobe -> github.com/spirrobe/
"""


def ConcPerCCM(bins,
               windspeed=[1],  # in m/s
               samplearea=0.298,  # as area in square millimeters
               samplefrequency=10,  # as Hz
               combined=True):
    import numpy as np
    # windspeed defaults to 1 meter per seconds windspeed
    # samplearea defaults to 0.298 as given by CDP specs
    # make cm/s out of it
    windspeed = (
        windspeed.copy()
        if type(windspeed) == np.ndarray
        else np.asarray(windspeed, dtype=np.float)
    )

    # convert windspeed first to cm per s and then to cm adjusted to
    # sampling frequency, as we look for cm and we still have cm/s
    windspeed = windspeed * (100 / samplefrequency)

    if np.nanmin(windspeed) < 0.0:
        print('*' * 30 + '\nWarning: Wind speed contains negative values.\n' +
              ' Concentration and derived parameters may be ' +
              'negative therefore\n' +
              'Use np.abs(x) on result to ammend this if needed\n' + '*' * 30)

    elif np.nanmax(windspeed) > 1000:
        print('*' * 30,
              'Warning: Wind speed max seems high with >1000 cm ',
              'per samplinginterval', samplefrequency,
              '\n Make sure all parameters are proper',
              '*' * 30)

    samplearea /= 10**2  # make cm instead of millimeters out of it
    volume = windspeed * samplearea  # make a volume out of it

    # set volume to nan if its zero because then its not physically correct
    volume[volume == 0] = np.nan

    if combined:
        concperccm = np.nansum(bins, axis=1) / volume
    else:
        concperccm = bins / np.repeat(volume[:, np.newaxis],
                                      bins.shape[1], axis=1)

    concperccm = np.nan_to_num(concperccm)

    return concperccm
