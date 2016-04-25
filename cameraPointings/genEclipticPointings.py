''' 
Generate array of geocentric ecliptic coordinates that boresights of
the four TESS cameras will point to in an ecliptic focused extended
mission (1yr)
'''

import numpy as np

fov = 24
missionLen = 1 # year
nCams = 4
nPointingsPerHemi = 13 # per year
nPointings = nPointingsPerHemi * missionLen # total number of camera fields from mission

# Pointing number, camera number, ecl lat/ ecl long
pointingCoords = np.zeros([nPointings, nCams, 2]) 

firstNorthPoint = [[0., 0.], [0., 1.*fov], [0., 2.*fov], [0., 3.*fov]]
pointingCoords[0, :, :] = firstNorthPoint

longDiff = 360./float(nPointingsPerHemi)
singleLongRot = [[0, longDiff], [0, longDiff], [0, longDiff], [0, longDiff]]

j = 0
while j < 12:
    pointingCoords[j+1] = np.mod(pointingCoords[j] + singleLongRot, 360.)
    j += 1

'''
firstSouthPoint = [[-(90.-1.5*fov), 0.], [-(90.-0.5*fov), 0.], [-(90.-0.5*fov), 180.], 
        [-(90.-1.5*fov), 180.]]
pointingCoords[j+1] = firstSouthPoint
j += 1

while j < 25:
    pointingCoords[j+1] = pointingCoords[j] + singleLongRot
    pointingCoords[j+1,:,1] = np.mod(pointingCoords[j+1,:,1], 360.)
    j += 1

nCoords = len(pointingCoords) * nCams
assert (nCoords == nPointings * nCams)
'''

# output format:
# pointingNum / ecl lat / ecl long
k = 0
for i in range(len(pointingCoords)):
    for j in range(nCams):
        print k//nCams, j, pointingCoords[i,j,0], pointingCoords[i,j,1]
        k += 1

# Now, being too lazy to learn proper file I/O with python, run in terminal with python,
# bash pipe it into nominalCamPointings.txt
