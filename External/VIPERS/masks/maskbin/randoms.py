#!/usr/bin/python
"""randoms.py
Ben Granett ben.granett@brera.inaf.it
Time-stamp: <2011-05-09 16:21:41 ben>

Wrapper routines for Mangle ransack and polyid to generate random
points and test whether they fall within a survey mask.
"""

import sys,os
import argparse
import numpy as N
import random
import time


def randoms(n=1e3, fields=[], holes=None, verbose=False):
    """ Generate random points inside the survey area using Mangle routines.

    About hole files: The mangle hole files should be pixelized for
    efficient lookup.  To combine many hole masks together into one
    file, you may use the Mangle snap command with the -S flag:
    snap -S -m1e-5s file1 file2 ... outputmask

    Timing info: the mangle code is fast and scales as ~n for large n.
    The scaling with number of polygons has not been investigated.
    The time to generate 1e6 points in the W1 area including the
    photometric mask is 16 seconds on Ben's desktop.

    inputs:
      n - number of points to generate.
      fields -  a list of paths to survey field mask files in mangle format.
      holes - optionally set the mask file describing holes in mangle format.
      verbose - print messages or not.

    outputs: a numpy array in the form array([x,y]) where x and y are
      length-n vectors of coordinates.
    """
    t0 = time.time()

    masks = " ".join(fields)
    tmppath = os.tempnam(".")


    cmd = "poly2poly -q -m1e-5s %s %s"%(masks,tmppath)
    r=os.system(cmd)
    if not r==0:
        print >>sys.stderr, cmd
        print >>sys.stderr,"poly2poly failed"
        sys.exit(r)

    seed = int(random.uniform(0,1e10))
    out = os.tempnam(".")

    if not holes==None:        
        cmd = 'ransack -q -m1e-5s -r%i -c%i %s - | polyid -q %s - %s'%(n,seed,tmppath, holes, out)
        howmany=2
    else:
        cmd = 'ransack -q -m1e-5s -r%i -c%i %s %s'%(n,seed,tmppath, out)
        howmany=3

    if verbose: print cmd
    r = os.system(cmd)
    if not r==0:
        print >>sys.stderr, cmd
        print >>sys.stderr, "!! ransack failed with error",r
        sys.exit(r)

    x = []
    y = []
    for i,line in enumerate(file(out)):
        if i==0: continue
        w = line.split()
        if len(w)==howmany:
            x.append(float(w[0]))
            y.append(float(w[1]))

    xy = N.array([x,y])
    
    os.unlink(out)
    if os.path.exists(tmppath): os.unlink(tmppath)

    dt = time.time()-t0
    if verbose: print "Done with %i randoms, time %g sec"%(n,dt)

    return xy

def timeit(fields,holes,loops=3, runs = [1e3,1e4,1e5,2e5,4e5,1e6,2e6,4e6,8e6]):
    """ Execute the code with a range of n and compile the timing statistics.
    """
    import pylab
    runs = N.array(runs)
    data = []
    for n in runs:
        tt = []
        for i in range(loops):
            t0 = time.time()
            randoms(n=n, fields=fields, holes=holes, verbose=False)
            dt = time.time()-t0
            print n,dt
            tt.append(dt)
        data.append(N.mean(tt))

    data = N.array(data)
    
    pylab.loglog(runs,data,"o")

    ii = N.where(runs>=2e5)
    
    fit = N.polyfit(N.log10(runs[ii]),N.log10(data[ii]),1)
    x = N.linspace(1e5,max(runs),10)
    y = N.polyval(fit,N.log10(x))
    pylab.plot(x,10**y,"--")

    print "fit",fit
    
    pylab.show()
    

if __name__=="__main__":
    parser = argparse.ArgumentParser(description=
                                     'Generate random points inside a mask using Mangle tools.')

    parser.add_argument('fields',metavar='field_file', nargs='+',
                        help='path to mask file')
    parser.add_argument('--holes',metavar='hole_file', nargs=1,default=None,
                        help='path to mask file representing holes. Accepts only one file.')
    parser.add_argument('-n',nargs=1,default=[1e4],type=float,
                        help='number of points to generate')

    parser.add_argument('--out',nargs=1,default=[None],
                        help='output file in which to write the random coordinates.  If this flag is not specified, the output will be written to STDOUT.')

    parser.add_argument('--time',action='store_true',help='collect timing statistics')
    parser.add_argument('--plot',action='store_true',help='show a plot with matplotlib')
    
    args = parser.parse_args()

    if not args.holes == None:
        holes = args.holes[0]
    else:
        holes = None
        
    outfile = args.out[0]
    n = args.n[0]
    
    verbose=True
    if outfile==None:
        verbose=False

    if args.time:
        timeit(args.fields,holes)
        sys.exit()

    # call the routine
    x,y = randoms(n=n, fields=args.fields, holes=holes, verbose=verbose)
    
    if outfile == None:
        out = sys.stdout
    else:
        out = file(outfile,'w')

    for i in range(len(x)):
        print >>out, "%f %f"%(x[i],y[i])

    if not outfile==None:
        out.close()
        print "> randoms written to %s"%outfile

    if args.plot:
        import pylab
        pylab.plot(x,y,",",alpha=0.2)
        pylab.show()
