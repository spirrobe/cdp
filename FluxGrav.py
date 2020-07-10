#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: spirrobe -> github.com/spirrobe/
"""

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
