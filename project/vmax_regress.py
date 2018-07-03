from __future__ import print_function
from scipy import stats
import numpy as np

n = 4
nfirstlines = []

with open("VMAXDATA.txt") as f, open("VMAXDATAtemp.txt", "w") as out:
    for x in xrange(n):
        nfirstlines.append(next(f))
    for line in f:
        out.write(line)

i = 0

# It will throw an error at the end when it reaches EOF, but I cba to deal with
# it because by then I will have already had all my stat values
with open("VMAXDATAtemp.txt") as fobj:
    while True:
        print(fobj.readline().strip('\n') + ', ', end='')
        i += 1
        x_temp = fobj.readline().strip('[]\n').split(', ')
        x = [float(j) for j in x_temp]
        i += 1
        y_temp = fobj.readline().strip('[]\n').split(', ')
        y = [float(j) for j in y_temp   ]
        i += 1

        slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        print(str(r_value**2) + ', ' + str(p_value))
