# renaming columns so idl's read_csv is happy
import pandas as pd
import numpy as np

dat = pd.read_csv('10000brightest_BC_SGBstars.csv')


dat.columns = [u'glon', u'glat', u'elon',
                u'elat', u'maxNpt', u'dart', u'rad',
                u'mass', u'teff', u'cosi', u'age', u'mini', u'feh', u'logg', u'lum',
                u'av', u'dm', u'mv', u'mic', u'mj', u'v', u'g', u'r', u'i', u'imag',
                u'z', u'kp', u't', u'j', u'icsys', u'jsys', u'ksys', u'tsys', u'kpsys',
                u'mvsys', u'micsys', u'mjsys', u'h', u'k', u'pri', u'sec', u'spl',
                u'ffi', u'gc']

dat.to_csv('asteroseisStars.csv')
