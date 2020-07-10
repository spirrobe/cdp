#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: spirrobe -> github.com/spirrobe/
"""


def ED(bins,
       # the 2 as first binsize is owed to the ADC treshold, which denotes
       # the lower binborder
       binsizes=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20,
                 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50],
       binsize_as_diameter=True):
    import numpy as np
    # binsizes are NOT the treshold values that are sent via setup cmd
    if len(binsizes) < bins.shape[1]:
        print('binsizes length must be at least one longer than bins length')
        return False

    midpoints = [(binsizes[_+1]+binsizes[_])/2 for _ in range(len(binsizes)-1)]
    if binsize_as_diameter:
        midpoints = [midpoint/2 for midpoint in midpoints]

    _r2 = np.zeros(bins.shape)
    _r3 = np.zeros(bins.shape)

    # midpoints divided by 1000 as a security precaution of float overflow
    for mno, midp in enumerate(midpoints):
        _r2[:, mno] = bins[:, mno].squeeze() * midp**2
        _r3[:, mno] = bins[:, mno].squeeze() * midp**3

    # set stuff with zero to -1 so we know that negative numbers
    # which are not possible here were once zero and can be set to that after
    # the division

    _r2 = np.nansum(_r2, axis=1)
    _r3 = np.nansum(_r3, axis=1)
    _r2[_r2 == 0] = -1.0
    ed = _r3/_r2 * 2
    ed[ed < 0] = 0

    return ed
