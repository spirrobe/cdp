#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: spirrobe -> github.com/spirrobe/
"""



def LWC(inputdata, # needs to be in # / m^3
        data_is_conc=True,
        binsizes=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20,
                  22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50],
        windspeed=[1],  # in m/s
        samplearea=0.298,  # as area in square millimeters
        samplefrequency=10,
        combined=True,
        quiet=True,
        ):
    import numpy as np
    inputdata = inputdata.copy()
    # binsizes are NOT the treshold values that are sent via setup cmd
    if not data_is_conc:
        from ConcPerCCM import ConcPerCCM
        _lwc = ConcPerCCM(inputdata,
                          windspeed=windspeed,
                          samplearea=samplearea,
                          samplefrequency=samplefrequency,
                          combined=False,
                          )
        # if you have given bins directly, ConcPerCCM will give you back #/cm^3
        # which we need to scale for the LWC

        # convert the concentration from cm^3 (default) to m^3
        # as we get per cm^3 from ConcPerCCM (as the name implies)
        _lwc *= 10**6

    else:
        _lwc = np.asarray(inputdata, np.float)

    if len(binsizes) < _lwc.shape[1]:
        print('binsizes length must be one longer than length of bins')
        return False

    # midpoints are in micrometers
    midpoints = [sum(_)/2 for _ in zip(binsizes[1:], binsizes[:-1])]

    dropsize = [(midpoint)**3 * np.pi/6 * 10**-12 for midpoint in midpoints]

    dropsize = np.asarray(dropsize)

    if np.nanmax(_lwc) > 200000 and not quiet:
        print('Warning from LWC():')
        print('Are you sure the concentration is in # / m^3?')
        print('It seems somewhat high with a max of', np.nanmax(_lwc))

    if type(_lwc) != np.ndarray:
        _lwc = np.asarray(_lwc, dtype=np.float)

    _lwc = (_lwc * dropsize)

    # choose either the sum of bins for a record or the single bins
    if combined:
        # summarize the concenctrations of the bins toghether
        return np.nansum(_lwc, axis=1)
    else:
        # do not summarize, each bin is already the lwc
        return _lwc

