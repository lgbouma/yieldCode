''' 
Generate array of geocentric ecliptic coordinates that boresights of
the four TESS cameras will point to in the nominal (2017-19) mission. 
'''

import numpy as np

radPerDeg = 2.*np.pi/360.
fov = 24*radPerDeg

missionLen = 24 # months

nCams = 4
nPointingsPerHemi = 13
nPointings = nPointingsPerHemi * 2

# Pointing number, camera number, ecl lat/ ecl long
pointingCoords = np.zeros([nPointings, nCams, 2]) 

firstNorthPoint = [[90., 0.], [90.-24., 0.], [90.-48., 0.], [90.-72., 0.]]
pointingCoords[0, :, :] = firstNorthPoint

longDiff = 360./float(nPointingsPerHemi)
singleLongRot = [[0, longDiff], [0, longDiff], [0, longDiff], [0, longDiff]]

j = 0
while j < 12:
    pointingCoords[j+1] = pointingCoords[j] + singleLongRot
    j += 1

firstSouthPoint = [[-90., 0.], [-90.+24., 0.], [-90.+48., 0.], [-90.+72., 0.]]
pointingCoords[j+1] = firstSouthPoint
j += 1

while j < 25:
    pointingCoords[j+1] = pointingCoords[j] + singleLongRot
    j += 1


nCoords = len(pointingCoords) * nCams
assert (nCoords == nPointings * nCams)

# output format:
# pointingNum / ecl lat / ecl long

k = 0
for i in range(len(pointingCoords)):
    for j in range(nCams):
        print k//nCams, j, pointingCoords[i,j,0], pointingCoords[i,j,1]
        k += 1

# Now, being too lazy to learn proper file I/O with python, run in terminal with python,
# bash pipe it into nominalCamPointings.txt
