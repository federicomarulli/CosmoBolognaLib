#!/usr/bin/python
"""inside.py
Ben Granett ben.granett@brera.inaf.it
Time-stamp: <2011-05-09 16:32:29 ben>

Wrapper routine for Mangle polyid to find points which are
inside and outside the mask.
"""
import sys,os
import argparse
import numpy as N
import time


parser = argparse.ArgumentParser(
    description='Separate points that are inside and outside the mask')

parser.add_argument('--in', dest='inside', metavar='interior_file',nargs=1,
                    help='output file where interior points will be written')
parser.add_argument('--out',  metavar='exterior_file',nargs=1,
                    help='output file where exterior points will be written')
parser.add_argument('mask_file', nargs=1, help='mask file in mangle format')
parser.add_argument('points_file', nargs=1, help='input file to be passed to polyid')
args = parser.parse_args()

inpath = args.inside[0]
outpath = args.out[0]
points = args.points_file[0]
mask = args.mask_file[0]


tmpfile = os.tempnam(".")
cmd = 'polyid %s %s %s'%(mask, points, tmpfile)
r = os.system(cmd)
if not r==0:
    print >>sys.error, "polyid failed",r
    sys.exit(r)

# how many numbers per line are in the input file...
# probably 2, but here we look at a few lines to make sure
count = -1
for line in file(points):
    try:
        w = line.split()
        float(w[0])
        float(w[1])
        count = len(w)
        break
    except:
        continue
if count<=0:
    print >>sys.stderr, "Nothing was found in the input file"
    sys.exit(1)

IN = file(inpath,'w')
OUT = file(outpath,'w')

cin = 0
cout = 0

lcount = 0
for line in file(tmpfile):
    lcount += 1
    if lcount <= 1: continue  # skip the first line, which just says az, alt
    w = line.split()
        
    if len(w)==count:
        OUT.write(line)
        cout += 1
    elif len(w)==count+1:
        print >>IN, " ".join(w[:count])  # this leaves off the polygon id at the end of the line
        cin += 1
    else:
        pass
IN.close()
OUT.close()

print "number of lines written to %s: %i"%(inpath,cin)
print "number of lines written to %s: %i"%(outpath,cout)
