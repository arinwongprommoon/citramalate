#!/usr/bin/env python
# Gets the mins and maxes from CSV
# Used to create minmax_citra.csv

import csv

with open('1db_ana_minmax.csv', 'wb') as f:
    pass

with open('1db_ana_combined.csv', 'rb') as csvfile:

    fluxreader = csv.reader(csvfile, skipinitialspace=True, delimiter=',', quoting=csv.QUOTE_NONE)
    fluxlist = list(fluxreader)
    #print fluxlist
 
    for row in fluxlist:
        r = []
        for el in row:
            try:
                r.append(float(el))
            except ValueError:
                pass
        print(min(r), max(r))
        
        with open('1db_ana_minmax.csv', 'ab') as outputfile:
            outputwriter = csv.writer(outputfile, delimiter=',')
            outputwriter.writerow([min(r), max(r)])
