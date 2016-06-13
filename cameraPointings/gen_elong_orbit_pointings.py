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
nPointings = nSectors * 2 # 52 over both years, 26 over one year

# Pointing number, camera number, ecl lat/ ecl long
pointingCoords = np.zeros([nPointings, nCams, 2]) 

longDiff = 360./float(nSectorsPerHemi)
singleLongRot = [[0, longDiff], [0, longDiff], [0, longDiff], [0, longDiff]]

pntg_ind = 0
firstNorthPoint = [[90.-3*fov/2., 180.], [90.-fov/2., 180.], \
                   [90.-fov/2., 0.], [90.-3*fov/2., 0.]]
pointingCoords[pntg_ind] = firstNorthPoint
# First do 6 npole sectors (get longest term continuity)
while pntg_ind < 12:
    if pntg_ind == 0:
        pass
    else:
        pointingCoords[pntg_ind] = pointingCoords[pntg_ind-1] + singleLongRot \
                                   if pntg_ind % 2 == 0 else pointingCoords[pntg_ind-1]
    pntg_ind += 1

# COLLECT RELEVANT ELAT / ELON INFO
orbit_number, camera_number, elat, elon = [], [], [], []
k = 0
for i in range(12):
    for j in range(nCams):
        orbit_number.append(k//nCams)
        camera_number.append(j)
        elat.append(pointingCoords[i,j,0])
        elon.append(pointingCoords[i,j,1])
        k += 1

elon = [e%360. for e in elon]

# Do 7 elong sectors
center_elon = elon[-1]+longDiff # to be advanced on as "central point" for elong section
firstEclipPoint = [[0, center_elon-1.5*fov], [0, center_elon-0.5*fov], \
                   [0, center_elon+0.5*fov], [0, center_elon+1.5*fov]]
pointingCoords[pntg_ind] = firstEclipPoint
pntg_ind += 1
while pntg_ind < 26:
    pointingCoords[pntg_ind] = pointingCoords[pntg_ind-1] + singleLongRot \
                               if pntg_ind % 2 == 0 else pointingCoords[pntg_ind-1]
    pntg_ind += 1


'''
# then do 7 sectors of elong:
firstNorthPoint = [[0, 180.], [90.-fov/2., 180.], \
                   [90.-fov/2., 0.], [90.-3*fov/2., 0.]]

while pntg_ind < 12:
    if pntg_ind == 0:
        pass
    else:
        pointingCoords[pntg_ind] = pointingCoords[pntg_ind-1] + singleLongRot \
                                   if pntg_ind % 2 == 0 else pointingCoords[pntg_ind-1]
    pntg_ind += 1
    '''

'''
pntg_ind += 1

while pntg_ind < 52:
    pointingCoords[pntg_ind] = pointingCoords[pntg_ind-1] + singleLongRot \
                               if pntg_ind % 2 == 0 else pointingCoords[pntg_ind-1]
    pntg_ind += 1
'''

# OUTPUT:
orbit_number, camera_number, elat, elon = [], [], [], []
k = 0
for i in range(nPointings):
    for j in range(nCams):
        #print k//nCams, j, pointingCoords[i,j,0], pointingCoords[i,j,1]
        orbit_number.append(k//nCams)
        camera_number.append(j)
        elat.append(pointingCoords[i,j,0])
        elon.append(pointingCoords[i,j,1])
        k += 1

elon = [e%360. for e in elon]

k = 0
for ix in range(len(elon)):
    print orbit_number[ix], camera_number[ix], elat[ix], elon[ix]


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
