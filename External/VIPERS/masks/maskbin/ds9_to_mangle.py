#!/usr/bin/python
""" ds9_to_mangle
Ben Granett ben.granett@brera.inaf.it
Time-stamp: <2011-05-06 16:40:20 ben>

A tool to convert a ds9 region file that defines polygons to a mask
that can be input into mangle2.  The generated masks are not pretty,
often with long skinny triangles, but they work.

requirements:
 numpy
 matplotlib for plots
"""

import sys,os
import argparse
import numpy as N


def to_xyz(ra, dec):
    """convert longitude and latitude (ra,dec) angles to an xyz vector"""
    x = N.cos(dec*N.pi/180)*N.cos(ra*N.pi/180)
    y = N.cos(dec*N.pi/180)*N.sin(ra*N.pi/180)
    z = N.sin(dec*N.pi/180)
    return N.array([x,y,z])

def to_radec(x,y,z):
    """ convert an xyz vector to longitude and latitude angles."""
    ra = N.arctan2(y,x)*180/N.pi
    dec = 90-N.arctan2(N.sqrt(x**2+y**2), z)*180/N.pi
    return N.array([ra,dec])


def findcenter(poly):
    """ The geometric center of the polygon. """
    return N.mean(poly,axis=0)

def cross(a,b,c):
    """ cross product of two vectors defined by 3 points. """
    xyz = N.array([a,b,c])

    d1 = xyz[0]-xyz[1]
    d2 = xyz[2]-xyz[1]

    return N.cross(d1,d2)


def outsideTriangle(tri, points):
    """Determine if any points are within the triangle"""
    v = N.zeros((3,3))
    for i in range(3):
        j = tri[i]
        k = tri[(i+1)%3]
        v[i] = N.cross(points[j],points[k])

    for l in range(len(points)):
        if l in tri: continue
        p = N.array(points[l])
        r = N.dot(p,v.transpose())
        if N.all(r>0):
            return False
    return True

def reorder(poly):
    """ Test if a polygon is ordered properly.
    """
    c = findcenter(poly)
    out = N.zeros(len(poly))
    for i in range(len(poly)):
        j = (i+1)%len(poly)
        out[i] = N.dot(c,cross(c,poly[i],poly[j]))
    if N.all(out < 0):
        #print "--ok"
        return False
    if N.all(out > 0):
        print "reordering!",poly
        return True
    raise("crappy ordering")

def isTiny(tri):
    """ test if a triangle is really small """
    x = cross(*tri)
    m = N.sqrt(N.dot(x,x))

    if m < 1e-13:
        print "tiny tri yes", m
        return True
    return False

def chop(poly, direction=1):
    """ Divide a polygon into triangles. """
    n = len(poly)
    if n < 3: return None

    center = findcenter(poly)
    c = to_radec(*center)

    stuckinaloop = 0
    j = 0
    out = []
    while len(poly)>3:
        stuckinaloop += 1
        if stuckinaloop > 1000:
            print poly
            print "stuckinaloop! bailing on this poly"
            return None
        
        n = len(poly)
        j %= n
        i = (j-1)%n
        k = (j+1)%n
                
        p = [poly[i],poly[j],poly[k]]

        if isTiny(p):
            poly.pop(j)
            continue

        if not outsideTriangle((i,j,k), poly):
            j += 1
            j %= n
            continue

        s = N.dot(center,cross(poly[i],poly[j],poly[k]))

        if s > 0:
            j += 1
            j %= n
            continue

        if reorder(p):
            print "shouldnt need a reorder here unless something is wrong"
            p = p[::-1]
            
        out.append(p)
        poly.pop(j)

        
    assert(len(poly)==3)
    if isTiny(poly): return out

    if reorder(poly):
        out.append(poly[::-1])
    else:
        out.append(poly)
    return out


def plot(poly, ls="ro-"):
    import pylab
    n = len(poly)
    xx = []
    yy = []
    for i in range(len(poly)):
        j = (i+1)%n
        x,y = to_radec(*poly[i])
        x2,y2=to_radec(*poly[j])
        xx.append(x)
        xx.append(x2)
        yy.append(y)
        yy.append(y2)        
    pylab.fill(xx,yy,lw=0,alpha=.2)

def plotpoint(p):
    """ """
    import pylab
    c = to_radec(*p)
    pylab.plot([c[0]],[c[1]],"s")
    


def printPoly(poly, out=sys.stdout):
    for i in range(len(poly)):
        x,y = to_radec(*poly[i])
        out.write("%13.10f %13.10f "%(x,y))
    out.write("\n")

def mangle(vertfile, outfile=None, pixscheme='s0,11',unifyflag=False, pixelize=True, balkanize=True, unify=True, cleanup=True):
    """ Convert the mask to Mangle format and carry out mangle processing
    steps including pixelizing and snapping."""
    if outfile == None:
        outfile = '%s.mangle'%vertfile
    
    if pixelize:
        os.system('pixelize -m1e-5s -iv -P%s %s %s.pix'%(pixscheme,vertfile,vertfile))
        next = '%s.pix'%vertfile
        fmt = ''
    else:
        next = vertfile
        fmt = '-iv'

    os.system('snap  -m1e-5s  %s %s %s.snap'%(fmt, next, vertfile))

    if balkanize:
        os.system('balkanize -m1e-5s %s.snap %s.balk'%(vertfile,vertfile))
    else:
        os.system('cp %s.snap %s.balk'%(vertfile,vertfile))

    if unify:
        U=''
        if unifyflag: U='-U'
        os.system('unify %s -m1e-5s  %s.balk %s.uni'%(U, vertfile,vertfile))
    else:
        os.system('cp %s.balk %s.uni'%(vertfile,vertfile))
        
    os.system('poly2poly -m1e-5s -k1e-10 %s.uni %s'%(vertfile, outfile))

    if os.path.exists('%s.pix'%vertfile):
        os.unlink('%s.pix'%vertfile)
    os.unlink('%s.snap'%vertfile)
    os.unlink('%s.balk'%vertfile)
    os.unlink('%s.uni'%vertfile)
    os.unlink('%s'%vertfile)


    return outfile



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=
                                     """A tool to convert a ds9 region
file that defines polygons to a mask that can be input into mangle2.
The generated masks are not pretty, often with long skinny triangles,
but they work.""")
    
    parser.add_argument('--multi', action='store_true',
                        help='Overlapping polygons will not be combined (skips balkanize and unify steps).  Use this to speed up processing.')
    parser.add_argument('--notpix', action='store_true',
                        help='The output should not be pixelized (sets unify -U flag).')
    parser.add_argument('--pixscheme', nargs=1, default=['s0,11'],
                        help='Specify the pixelization scheme for mangle pixelize (default s0,11)')

    parser.add_argument('maskfile',nargs='+',
                        help='path to region file')
    parser.add_argument('--overwrite',action='store_true',help='overwrite previous output files')
    parser.add_argument('--plot',action='store_true',help='Save each mask to a plot with matplotlib')

    args = parser.parse_args()

    if args.plot:
        import plotmask


    unify = True
    balkanize = True
    if args.multi:
        unify = False
        balkanize = False

    unifyflag = args.notpix
    pixscheme = args.pixscheme[0]
    
    for filename in args.maskfile:

        mangleout = '%s.mangle'%filename
        pngname = filename+".png"

        if not args.overwrite and os.path.exists(mangleout):
            print "skipping",filename
            if args.plot: plotmask.saveimage(mangleout, pngname,overwrite=True)
            continue
        
        outfile = filename+".vert"

        poly = []
        for line in file(filename):
            """ parse the ds9 region file.  look for lines that have
            polygon(a,b,c,d,...)
            """
            line = line.strip()
            if line.startswith("#"): continue

            i = line.find("(")
            j = line.find(")")
            
            if not line[:i].endswith("polygon"): continue
            line = line[i+1:j]

            xx = []
            for v in line.split(","):
                try:
                    xx.append(float(v))
                except ValueError:
                    break
            
            ra = N.array(xx[::2])
            dec = N.array(xx[1::2])
        
            xyz = N.transpose(to_xyz(ra,dec)).tolist()
        
            poly.append(xyz)

        out = file(outfile,"w")
        i = 0
        for pp in poly:
            b = chop(pp)
            if b == None: continue
            for p in b:
                plot(p)

                printPoly(p, out)
                i +=1
        out.close()

        # free memory used by poly array
        del poly

        outfilev = mangle(outfile, mangleout, unifyflag=unifyflag, pixscheme=pixscheme,
                          pixelize=True, unify=unify, balkanize=balkanize)

        if args.plot:
            plotmask.saveimage(mangleout, pngname, overwrite=args.overwrite)

