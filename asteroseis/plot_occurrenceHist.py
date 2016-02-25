'''
The matrix of occurrence rates contains only 1 bin of planetary radius and period. 
    The occurrence rate for this single bin is 10%.
Planetary radii are drawn from a log-normal distribution. The underlying normal 
    distribution has mean 1 (i.e., mean is e^1=2.7R_earth and standard deviation 0.5.)
Orbital periods are drawn from a log-normal distribution. The underlying normal 
    distribution has mean 3 (i.e., mean is e^3=20.1day) and standard deviation 1.31. 

A variable _x_ has log-normal distribution if _log(x)_ is normally distributed
(note this is log_e = ln). Thus the PDF for the log-normal distribution is
$$ p(x) = \frac{1}{\sigma x \sqrt{2\pi}} \exp{-\frac{(\ln(x) - \mu)**2}{2\sigma**2}}

Doing this in IDL is literally behind another paywall (...).

So we do it via a janky I/O procedure.
'''

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(color_codes=True, style='whitegrid')
sns.set_context('talk', font_scale=1.5, rc={'lines.linewidth': 2.5})

nDraws = 1e6
nBins = 30
radii = np.random.lognormal(1., 0.5, nDraws)
periods = np.random.lognormal(3., 1.31, nDraws)

plt.close('all')
plt.ion()

plt.figure(1)
sns.distplot(radii, bins=nBins, kde=True, hist=True, norm_hist=False)
plt.xlabel('Planet radius [R_earth]')

plt.figure(2)
sns.distplot(periods, bins=nBins, kde=True, hist=True, norm_hist=False)
plt.xlabel('Planet period [day]')

np.savetxt('radiiVals.dat', radii, delimiter=',')
np.savetxt('periodVals.dat', periods, delimiter=',')
