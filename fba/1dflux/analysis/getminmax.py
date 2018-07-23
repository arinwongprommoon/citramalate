import csv

with open('MINMAX_WITHOUT_ATP_MAINTENANCE.csv', 'wb') as f:
    pass

with open('ALLFLUXCOMBINED_WITHOUT_ATP_MAINTENANCE.csv', 'rb') as csvfile:

    fluxreader = csv.reader(csvfile, skipinitialspace=True, delimiter=',', quoting=csv.QUOTE_NONE)
    fluxlist = list(fluxreader)
    print fluxlist
 
    for row in fluxlist:
        r = []
        for el in row:
            try:
                r.append(float(el))
            except ValueError:
                pass
        print(min(r), max(r))
        
        with open('MINMAX_WITHOUT_ATP_MAINTENANCE.csv', 'ab') as outputfile:
            outputwriter = csv.writer(outputfile, delimiter=',')
            outputwriter.writerow([min(r), max(r)])
