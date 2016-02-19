'''
Make many lists (5e5 long) of random coordinates for each tile.
'''

import numpy as np
from os import listdir
from os.path import isfile, join

trilegalPath = '/home/luke/Dropbox/tessSimulation/trilegal'
onlyFiles = np.array([f for f in listdir(trilegalPath) if isfile(join(trilegalPath, f))])
onlyFiles.sort()
nTiles = len(onlyFiles)/3.

fourDigitNumList = onlyFiles[:nTiles]
numList = np.zeros(len(fourDigitNumList)).astype('str')
# Not sure if there's easy vectorization here
for i in range(len(fourDigitNumList)):
    fourDigitNumList[i]  = fourDigitNumList[i][2:-4] # based on PS's numbering scheme. keep as string.
    numList[i] = fourDigitNumList[i]
    print i
    while numList[i][0] == '0':
        if numList[i] == '0000':
            numList[i] = 0
            break
        numList[i] = numList[i][1:]
        
# Now, you have array of 4 digit strings, & their 1/2/3/4 digit counterparts
