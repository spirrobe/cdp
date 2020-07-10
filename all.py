# Some utilites for CDP Data, depends on numpy for the array calculations

# Dependencies between functions:
# MVD depends on LWC
# LWC depends on ConcPerCCM
# thus MVD depends also on ConcPerCCM
# ConcPerCCM depends on nothing
# ED depends on nothing

# The calling should be as follows
# first get the concentration by ConcPerCCM in # per cm^3 for a timestep
# get the liquid water concentration by LWC(ConcPerCCM(conc,combined=False))
# or by LWC(bindata, data_is_conc=False)
# then get the median volume diameter by
# MVD(LWC(ConcPerCCM(bindata, combined=False),combined=False))
# or by
# MVD(LWC(concentration,combined=False))
# or by
# MVD(lwc)
# Finally ED can be called at any time since it does not depend on anything


# Note 1:
# For MVD and LWC the windspeed, samplearea and frequency
# are optional and default to 1, 0.298, 0.1 Hz as given in ConcPerCCM as these
# are forwared to ConcPerCCM in case the concentration needs calculation

# Note 2:
# The windspeed has no bearing on the MVD as the windspeed will be the same
# for one record and while lwc depends on it, the formula used uses relative
# values and a mulitiplication of all values does not impact this
import numpy as np


def ConcPerCCM(bins,
               windspeed=[1],  # in m/s
               samplearea=0.298,  # as area in square millimeters
               samplefrequency=10,  # as Hz
               combined=True):

    # windspeed defaults to 1 meter per seconds windspeed
    # samplearea defaults to 0.298 as given by CDP specs
    # make cm/s out of it
    if not type(windspeed) == np.ndarray:
        windspeed = np.asarray(windspeed, dtype=np.float)
    else:
        windspeed = windspeed.copy()

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


def ED(bins,
       # the 2 as first binsize is owed to the ADC treshold, which denotes
       # the lower binborder
       binsizes=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20,
                 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50],
       binsize_as_diameter=True):

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
    inputdata = inputdata.copy()
    # binsizes are NOT the treshold values that are sent via setup cmd
    if not data_is_conc:
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




def MVD(inputdata,
        data_is_lwc=True,
        data_is_conc=False,
        binsizes=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20,
                  22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50],
        windspeed=[1],  # in m/s
        samplearea=0.298,  # as area in square millimeters
        samplefrequency=0.1,
        ):

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

# EXAMPLE DATA FOR MVD FROM PADS MANUAL V3.6.3
# col = [i.strip() for i in "binindex, binsize, conc, midpnt, LWC, pro, cum".split(',')]
# ix = [[1, 10, 100000, 15, 0.00018, 0.00058, 0.00058],
# [2, 20, 200000, 25, 0.00164, 0.00538, 0.00596],
# [3, 30, 300000, 35, 0.00673, 0.02214, 0.02810],
# [4, 40, 400000, 45, 0.01909, 0.06274, 0.09084],
# [5, 50, 500000, 55, 0.04356, 0.14320, 0.23404],
# # [5, 50, 500000, 55, 0.04356, 0.14320, 0.23404],
# [6, 60, 400000, 65, 0.05752, 0.18909, 0.42313],
# [7, 70, 300000, 75, 0.06627, 0.21786, 0.64099],
# [8, 80, 200000, 85, 0.06431, 0.21143, 0.85242],
# [9, 90, 100000, 95, 0.04489, 0.14758, 1.00000],
# # [10, 100, 0, 105, 0.00000, 0.00000, 1.00000],
# ]
# dx = {j: [k[i] for k in ix] for i,j in enumerate(col)}
# import pandas as pd
# z = pd.DataFrame(dx)
# zz = MVD(z['LWC'].copy(), data_is_lwc=True, binsizes=z['binsize'].copy().values.T, )
# print(zz)


# follows the Koschmeider equation, 3.91/sigma_e
def conc2visibility():
    pass


def FluxGrav(lwc,
             T=[20],
             binsizes=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20,
                       22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50],
             combined=True,
             lat=None,
             asl=0,
             latisrad=False,
             ):

    import numpy as np
    from mcr.met.constants.gravity import gravity
    from mcr.met.density.airdensity import airdensity
    from mcr.met.density.waterdensity import waterdensity
    from mcr.met.constants.viscosityair import viscosityair

    # get gravity, possibly adjusted to latitude and elevation asl
    g = gravity(lat=lat, asl=asl, latisrad=latisrad)

    if g >= 0:
        g *= -1

    # midpoints are in micrometers, i.e. to scale to meters we need 10**-6 here
    midpoints = [(binsizes[_+1]+binsizes[_])/2 * 10**-6
                 for _ in range(len(binsizes)-1)]

    midpoints = np.asarray(midpoints)
    # diameter is double the radius, but we already have diameter.
    # dropletdiameter = np.asarray(midpoints)*2

    # just make it an array.. ffs
    if type(T) in [list, range]:
        T = np.asarray(T)

    # ccalculate density difference
    rho_water = waterdensity(T)
    rho_air = airdensity(T)
    muair = viscosityair(T)

    # follows equation (2) from
    # 1. Westbeld, A. et al.
    # Fog deposition to a Tillandsia carpet in the Atacama Desert.
    # Ann. Geophys. 27, 3571–3576 (2009).
    # the original formula in Beswick 1991) employs the radius and this is
    # the diameter derived version
    sedspeed = -g * midpoints**2 * (rho_water - rho_air) / 18 / muair

    gravflux = lwc * sedspeed

    # choose either the sum of bins for a record or the single bins
    if combined:
        # summarize the gravfluxed of each bin toghether
        # summarize because this will be the complete water flux down to
        # the surface
        return np.nansum(gravflux, axis=1)
    else:
        # do not summarize, each bin is already the
        # graviational flux of that bin
        return gravflux


def Gonser2011(bincounts=np.ones(30),
               returnnewsizes=True,
               getbins=False):

    # directly from paper applied onto the 30 CDP bins
    # take note that this will also mean, that the passed in bins to
    # all derived quantities like LWC/Conc need to be adjusted for
    # the proper calculations thereof
    sca = {2.66: [[0, 0.66]],
           4.66: [[0, 0.34], [1, 1.0], [2, 0.66]],
           7.41: [[2, 0.34], [3, 1.0], [4, 1], [5, 0.41]],
           9.66: [[5, 0.59], [6, 1.0], [7, 0.66]],
           11.36: [[7, 0.34], [8, 1.0], [9, 0.36]],
           13.21: [[9, 0.64], [10, 1.0], [11, 0.21]],
           15.26: [[11, 0.79], [12, 1.26/2]],
           18.01: [[12, 0.74/2], [13, 2.0/2.0], [14, 0.01/2]],
           19.96: [[14, 1.95/2]],
           21.96: [[14, 0.04/2], [15, 1.96/2]],
           23.61: [[15, 0.04/2], [16, 1.61/2]],
           25.01: [[16, 0.39/2], [17, 1.01/2]],
           26.51: [[17, 0.99/2], [18, 0.51/2]],
           28.46: [[18, 1.49/2], [19, 0.46/2]],
           30.51: [[19, 1.54/2], [20, 0.51/2]],
           33.36: [[20, 1.49/2], [21, 1.36/2]],
           35.26: [[21, 0.64/2], [22, 1.26/2]],
           37.06: [[22, 0.74/2], [23, 1.06/2]],
           39.01: [[23, 0.94/2], [24, 1.01/2]],
           41.11: [[24, 0.99/2], [25, 1.11/2]],
           43.76: [[25, 0.89/2], [26, 1.76/2]],
           46.41: [[26, 0.24/2], [27, 2.0/2.0], [28, 0.41/2]],
           50: [[28, 1.59/2], [29, 2.0/2.0]],
           }

    if getbins:
        return [2.0] + list(sca.keys())

    binsizes = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20,
                22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50]

    if bincounts.shape[1] != 30:
        bincounts = bincounts.copy().transpose()


    _ = [0] * len(binsizes)
    output = np.zeros((bincounts.shape[0], len(sca.keys())))
    for newbinno, newbinsize in enumerate(sca.keys()):
        for scalepair in sca[newbinsize]:
            output[:, newbinno] += bincounts[:, scalepair[0]] * scalepair[1]

    if returnnewsizes:
        newsizes = [2.00]
        newsizes.extend(sorted(list(sca.keys())))
        return output,  newsizes
    else:
        return output


def Spiegel2012(bincounts,
                returnnewsizes=True):

    if bincounts.shape[1] != 30:
        bincounts = bincounts.copy().transpose()

    # directly from paper applied onto the 30 CDP bins
    # take note that this will also mean, that the passed in bins to
    # all derived quantities like LWC/Conc need to be adjusted for
    # the proper calculations thereof
    sca = {2.66: [[0, 0.66]],
           4.66: [[0, 0.34], [1, 1.0], [2, 0.66]],
           7.41: [[2, 0.34], [3, 1.0], [4, 1], [5, 0.41]],
           9.66: [[5, 0.59], [6, 1.0], [7, 0.66]],
           11.36: [[7, 0.34], [8, 1.0], [9, 0.36]],
           13.21: [[9, 0.64], [10, 1.0], [11, 0.21]],
           15.26: [[11, 0.79], [12, 1.26/2]],
           18.01: [[12, 0.74/2], [13, 2.0/2.0], [14, 0.01/2]],
           19.96: [[14, 1.95/2]],
           21.96: [[14, 0.04/2], [15, 1.96/2]],
           23.61: [[15, 0.04/2], [16, 1.61/2]],
           25.01: [[16, 0.39/2], [17, 1.01/2]],
           26.51: [[17, 0.99/2], [18, 0.51/2]],
           28.46: [[18, 1.49/2], [19, 0.46/2]],
           30.51: [[19, 1.54/2], [20, 0.51/2]],
           33.36: [[20, 1.49/2], [21, 1.36/2]],
           35.26: [[21, 0.64/2], [22, 1.26/2]],
           37.06: [[22, 0.74/2], [23, 1.06/2]],
           39.01: [[23, 0.94/2], [24, 1.01/2]],
           41.11: [[24, 0.99/2], [25, 1.11/2]],
           43.76: [[25, 0.89/2], [26, 1.76/2]],
           46.41: [[26, 0.24/2], [27, 2.0/2.0], [28, 0.41/2]],
           50: [[28, 1.59/2], [29, 2.0/2.0]],
           }

    binsizes = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20,
                22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50]
    _ = [0] * len(binsizes)
    output = np.zeros((bincounts.shape[0], len(sca.keys())))
    for newbinno, newbinsize in enumerate(sca.keys()):
        for scalepair in sca[newbinsize]:
            output[:, newbinno] += bincounts[:, scalepair[0]] * scalepair[1]

    if returnnewsizes:
        return output,  [2.0] + list(sca.keys())
    else:
        return output







