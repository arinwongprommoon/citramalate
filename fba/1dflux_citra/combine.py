import csv
import itertools as IT
import glob, os

filenames = glob.glob("*.csv")
#filenames = ['1.csv', '2.csv']
handles = [open(filename, 'rb') for filename in filenames]    
readers = [csv.reader(f, delimiter=',') for f in handles]

with  open('combined.csv', 'wb') as h:
    writer = csv.writer(h, delimiter=',', lineterminator='\n', )
    for rows in IT.izip_longest(*readers, fillvalue=['']*2):
        combined_row = []
        for row in rows:
            row = row[::-1]
            row = row[:2] # select the columns you want
            if len(row) == 2:
                combined_row.extend(row)
            else:
                combined.extend(['']*2)
        writer.writerow(combined_row)

for f in handles:
    f.close()
