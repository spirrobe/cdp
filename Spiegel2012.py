#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: spirrobe -> github.com/spirrobe/
"""

# From Spiegel et al 2012, more details in
# https://www.atmos-meas-tech.net/5/2237/2012/

def Spiegel2012(bincounts,
                returnnewsizes=True):
    import numpy as np
    if bincounts.shape[1] != 30:
        bincounts = bincounts.copy().transpose()

    # directly from paper applied onto the 30 CDP bins
    # take note that this will also mean, that the passed in bins to
    # all derived quantities like LWC/Conc need to be adjusted for
    # the proper calculations thereof
    sca = {
        2.66: [[0, 0.66]],
        4.66: [[0, 0.34], [1, 1.0], [2, 0.66]],
        7.41: [[2, 0.34], [3, 1.0], [4, 1], [5, 0.41]],
        9.66: [[5, 0.59], [6, 1.0], [7, 0.66]],
        11.36: [[7, 0.34], [8, 1.0], [9, 0.36]],
        13.21: [[9, 0.64], [10, 1.0], [11, 0.21]],
        15.26: [[11, 0.79], [12, 1.26 / 2]],
        18.01: [[12, 0.74 / 2], [13, 1.0], [14, 0.01 / 2]],
        19.96: [[14, 1.95 / 2]],
        21.96: [[14, 0.04 / 2], [15, 1.96 / 2]],
        23.61: [[15, 0.04 / 2], [16, 1.61 / 2]],
        25.01: [[16, 0.39 / 2], [17, 1.01 / 2]],
        26.51: [[17, 0.99 / 2], [18, 0.51 / 2]],
        28.46: [[18, 1.49 / 2], [19, 0.46 / 2]],
        30.51: [[19, 1.54 / 2], [20, 0.51 / 2]],
        33.36: [[20, 1.49 / 2], [21, 1.36 / 2]],
        35.26: [[21, 0.64 / 2], [22, 1.26 / 2]],
        37.06: [[22, 0.74 / 2], [23, 1.06 / 2]],
        39.01: [[23, 0.94 / 2], [24, 1.01 / 2]],
        41.11: [[24, 0.99 / 2], [25, 1.11 / 2]],
        43.76: [[25, 0.89 / 2], [26, 1.76 / 2]],
        46.41: [[26, 0.24 / 2], [27, 1.0], [28, 0.41 / 2]],
        50: [[28, 1.59 / 2], [29, 1.0]],
    }


    binsizes = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20,
                22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50]
    _ = [0] * len(binsizes)
    output = np.zeros((bincounts.shape[0], len(sca.keys())))
    for newbinno, newbinsize in enumerate(sca.keys()):
        for scalepair in sca[newbinsize]:
            output[:, newbinno] += bincounts[:, scalepair[0]] * scalepair[1]

    return (output, [2.0] + list(sca.keys())) if returnnewsizes else output


