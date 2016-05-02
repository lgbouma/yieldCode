# see 16-05-02-moonEarthOpaque
# we throw out roughly 1 in 8 (really, 1 in 8.66) fields near the ecliptic
# because of earth/moon crossings.

import numpy as np
import pandas as pd

dat = np.genfromtxt('ecliptic_coord_transparentEarthMoon.dat')
df = pd.DataFrame(dat, columns=['camera', 'subCamera', 'elat', 'elon'])

toKeep = np.random.permutation(np.arange(len(dat)))[:len(dat)-6] 

dfKeep = df.ix[toKeep]

dfKeep = dfKeep.sort(['camera', 'subCamera'])

keepOut = np.array(dfKeep)
np.savetxt('ecliptic_coord.dat', keepOut, fmt=['%d','%d','%.10f','%.10f'])
