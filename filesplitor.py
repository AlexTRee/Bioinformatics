#!/usr/bin/env python

chunksize = 330
fid = 1
with open('aaaaa.txt') as infile:
    f = open('output.%d.txt' %fid, 'w')
    for i,line in enumerate(infile):
        f.write(line)
        if not (i+1)%chunksize:
            f.close()
            fid += 1
            f = open('output.%d.txt' %fid, 'w')
    f.close()