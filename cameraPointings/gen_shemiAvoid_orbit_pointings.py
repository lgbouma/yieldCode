''' 
Generate array of geocentric ecliptic coordinates that boresights of
the four TESS cameras will point to in the nominal (2017-19) mission. 
Resolve them over each __orbit__ (not e.g., every 2 orbits).
Nomenclature:
    SECTOR := there are 13 of these per year
    POINTING := there are 26 per year (spacecraft pointing)
'''

import numpy as np

fov = 24. # degrees
missionLen = 1 # years
nCams = 4
nSectorsPerHemi = 13 # per year
nSectors = nSectorsPerHemi * missionLen # 26 over entire mission
nPointings = nSectors * 2 # 52 over both years

# Pointing number, camera number, ecl lat/ ecl long
pointingCoords = np.zeros([nPointings, nCams, 2]) 

longDiff = 360./float(nSectorsPerHemi)
singleLongRot = [[0, longDiff], [0, longDiff], [0, longDiff], [0, longDiff]]

pntg_ind = 0
firstSouthPoint = [[-90+fov/2., 0.], [-90.+3*fov/2., 0.], \
                   [-90.+5*fov/2., 0.], [-90.+7*fov/2., 0.]]
pointingCoords[pntg_ind] = firstSouthPoint
while pntg_ind < 26:
    if pntg_ind == 0:
        pass
    else:
        pointingCoords[pntg_ind] = pointingCoords[pntg_ind-1] + singleLongRot \
                                   if pntg_ind % 2 == 0 else pointingCoords[pntg_ind-1]
    pntg_ind += 1

'''
pntg_ind += 1

while pntg_ind < 52:
    pointingCoords[pntg_ind] = pointingCoords[pntg_ind-1] + singleLongRot \
                               if pntg_ind % 2 == 0 else pointingCoords[pntg_ind-1]
    pntg_ind += 1
'''

# OUTPUT:
k = 0
for i in range(nPointings):
    for j in range(nCams):
        print k//nCams, j, pointingCoords[i,j,0], pointingCoords[i,j,1]
        k += 1

# Now, being too lazy to learn proper file I/O with python, run in terminal with python,
# bash pipe it into nominalCamSectors.txt


'''
# nominal 2yr scenario
k = 0
for i in range(len(pointingCoords)):
    for j in range(nCams):
        print k//nCams, j, pointingCoords[i,j,0], pointingCoords[i,j,1]
        k += 1

# Now, being too lazy to learn proper file I/O with python, run in terminal with python,
# bash pipe it into nominalCamSectors.txt
'''

'''
# one 4 yr scenario: N, S, N, S hemis
pCoord = np.zeros([2*nSectors, nCams, 2])
pCoord[:nSectors] = pointingCoords
pCoord[nSectors:] = pointingCoords

k = 0
for i in range(len(pCoord)):
    for j in range(nCams):
        print k//nCams, j, pCoord[i,j,0], pCoord[i,j,1]
        k += 1
'''

'''
# another: N, S, S, N hemis
pCoord = np.zeros([2*nSectors, nCams, 2])
pCoord[:nSectors] = pointingCoords
pCoord[nSectors:nSectors*1.5] = pointingCoords[nSectors/2.:]
pCoord[nSectors*1.5:] = pointingCoords[:nSectors/2.]

k = 0
for i in range(len(pCoord)):
    for j in range(nCams):
        print k//nCams, j, pCoord[i,j,0], pCoord[i,j,1]
        k += 1
'''
