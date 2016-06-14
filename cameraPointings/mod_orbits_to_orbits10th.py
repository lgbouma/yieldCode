# We want to run these at 1/10th the regular orbit frequency (i.e. to count
# pointings as 1/20th of what we did on a per-sector basis)

import os
import numpy as np
f_names = [n for n in os.listdir('.') if '_orbits.dat' in n and '10th_orbits' not in n]

for f_name in f_names:
    dat = np.genfromtxt(f_name)
    out_name = f_name[:-10]+'10th_orbits.dat'
    o = open(out_name, 'w')
    for row in dat:
        print row[0],row[1],row[2],row[3]
        for i in range(10):
            o.write(str(int(row[0]))+' '+str(int(row[1]))+' '+str(row[2])+' '+str(row[3])+'\n')
    o.close()
