# Cloud Droplet Probe / Fogg Monitor Postprocessing functions
Code related to the **CDP-2** (Cloud Droplet Probe 2) of DMT (Droplet Measurement Technologies) that (should) also work for the **FM** (Fog Monitor).

Basic starting point is either direct sampling of CDP/FM when the DMT software is not available, for example when a laptop would need too much power.
The data therefore is the initial count of droplets per bin. This can also work with reprocessing of original data when e.g. the pitot tube on the fog monitor broke

This repository is linked to the publication [Droplet size distribution, liquid water content and water input of the seasonally variable, nocturnal fog in the Central Namib Desert by Spirig et al. 2021](https://doi.org/10.1016/j.atmosres.2021.105765). 

This repository also includes functions for rebinning based on the two publications:
- [Evaluating the capabilities and uncertainties of droplet measurements for the fog droplet spectrometer (FM-100) by Spiegel et al 2012](https://doi.org/10.5194/amt-5-2237-2012)
- [The relation between humidity and liquid water content in fog: an experimental approach by Gonser et al. 2011](https://doi.org/10.1007/s00024-011-0270-x)

## Requirements and extra
- numpy

The functions were used in conjunction with pandas time series of the bins (1 bin -> 1 column) as the numpy functions do broadcasting.

## Function names
MVD is the mean volume diameter
ED calculates the effective diameter
LWC calculates the liquid water content
ConcPerCCM calculates the concentration per cubic centimeter 

## internal dependencies
Dependencies between functions:
MVD depends on LWC (either calculated before or during use)
LWC depends on ConcPerCCM (either calculated before or during use)
meaning MVD depends also on ConcPerCCM
ConcPerCCM depends on nothing 
ED depends on nothing

## Basic usage and examples
As above the starting point are the bincounts (for the CDP per default 30 bins).
The first decision is whether rebinning should take place via Gonser2011 or Spiegel2012 as rebinning affects all other steps.
_Note that both Gonser2011 and Spiegel2012 have a keyword to return the bin sizes they create when rebinning without making a rebinning as this is required for all other functions (LWC, ED, MVD) except ConcPerCCM._

```
# dependencies
import pandas as pd
import numpy as np

# load functions from all.py
from all import *

# these variables you need to get from somewhere else
# windspeed
ws = 10 # in m s-1, will be converted in ConcPerCCM to cm s-1
freq = 10 # in Hz, 1 would be one sample per second, this is needed as the volume that passes by with the wind depends on how long we measured

# example loading csv file
data = pd.read_csv(FILENAME, index_col=0)

# rebin data according to Gonser2011
data, binsizes = Gonser2011(data)

# combined = False means ConcPerCCM is given for each bin instead of the sum of bins per sampling volume
# also needed for calculating lwd
conc = ConcPerCCM(data, windspeed=ws, samplefrequency=freq, combined=False)

# the factor 10**6 is to convert the # cm-3 to # m-3
lwc = LWC(conc * 10**6, combined=False, binsizes=binsizes) # the keywords windspeed and samplefrequency are ignored when data_is_conc = True (the default)
total_lwc = LWC(conc * 10**6, combined=True, binsizes=binsizes) # the keywords windspeed and samplefrequency are ignored when data_is_conc = True (the default)
# note also that LWC takes binsizes as a keyword in case Gonser2011 or Spiegel2011 have been applied; without these the default bin sizes are taken

# mvd is the diameter where half of the particle size distribution is reached
mvd = MVD(lwc, binsizes=binsizes) 
```

## Processing pipeline
The order of calling should be as follows:
0. Load the bincounts from your file somehow, depending on how they are saved
1. first get the concentration by ```conc = ConcPerCCM(bincounts)``` in per cm^3 for a timestep
2. get the liquid water concentration (per bin or total) with either
    -  ```lwc = LWC(bindata, data_is_conc=False)```
    -  ```lwc = LWC(bindata, data_is_conc=False, combined=False)```
    - ```lwc = LWC(ConcPerCCM(bincounts,combined=False))```
3. then get the median volume diameter by one of the following
    - ```MVD(LWC(ConcPerCCM(bindata, combined=False),combined=False))```
    - ```MVD(LWC(concentration,combined=False))```
    - ```MVD(lwc)```
4. Finally ```ED(bincounts)``` can be called at any time since it does not depend on anything

## Notes
### Note 1
For MVD and LWC the windspeed, samplearea and frequency
are optional and default to 1, 0.298, 0.1 Hz as given in ConcPerCCM as these
are forwarded to ConcPerCCM in case the concentration needs calculation

### Note 2
The windspeed has no bearing on the MVD as the windspeed will be the same
for one record and while the lwc calculation depends on it, the formula used for MVD uses relative
values and a mulitiplication of all values does not impact this
