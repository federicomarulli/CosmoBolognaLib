#!/usr/bin/python
"""plotmask.py
Time-stamp: <2011-05-09 11:20:11 ben>
Ben Granett ben.granett@brera.inaf.it

"""
import sys,os
import argparse
import pylab
import numpy as N
from matplotlib.patches import Polygon

def plotgraphics(fname, alpha=.5, linestyle='solid', edgecolor="b", linewidth=.1, facecolor="b", margin=.05, npoly=None, dontplot=False):
    fid = file(fname)
    line1 = fid.readline()
    np = int(line1.split()[0])
    print "n polygons in file:",np
    if npoly==None:
        npoly = np
    print "plotting this many polygons",npoly
    line = fid.readline()

    allx,ally = [],[]
    ax = pylab.gca()
    for i in range(npoly):
        line = fid.readline().strip()
        line2 = fid.readline().strip()
        if line2.startswith('graphics'): continue
        vert = [float(v) for v in line2.split()]
        if len(vert)==0: continue
        x,y = vert[:-1:2],vert[1::2]
        if dontplot:
            continue
        xy = zip(x,y)
        p = Polygon(xy)
        a = ax.add_artist(p)
        pylab.setp(a, alpha=alpha)
        pylab.setp(a, linestyle=linestyle)
        pylab.setp(a, edgecolor=edgecolor)
        pylab.setp(a, facecolor=facecolor)
        pylab.setp(a, linewidth=linewidth)
        
        allx.append(min(x))
        allx.append(max(x))
        ally.append(min(y))
        ally.append(max(y))
    bounds = (min(allx),max(allx),min(ally),max(ally))

    mx = (bounds[1]-bounds[0])*margin
    my = (bounds[3]-bounds[2])*margin

    pylab.xlim(bounds[0]-mx,bounds[1]+mx)
    pylab.ylim(bounds[2]-my,bounds[3]+my)
    return bounds

def plotds9reg(path, label=False,alpha=.5, linestyle="solid", linewidth=.5,edgecolor="k", facecolor="red", margin=.05,
               npoly=None, dontplot=False):
    """ """
    if label:
        label = path.split("_")[-1]
        label = label.split(".")[0][3:]
        label = "\n".join(label.split("Q"))
        print "label",label
    else:
        label = None
    
    allx,ally = [],[]
    ax = pylab.gca()
    for line in file(path):
        line = line.strip()
        #if line.startswith("#"):
        #    label = line.split('#')[1].strip()
        #    if label.startswith('W'):
        #        label = label[3:]
        #        label = label.replace(' ','\n')
        #    print label
        #    continue
        if not line.startswith('polygon'): continue

        line = line.replace("("," ")
        line = line.replace(")"," ")
        line = line.replace(","," ")
        
        
        i = line.find('#')
        if i > -1:
            w = line[:i].split()
        else:
            w = line.split()

        c = [float(v) for v in w[1:]]
        x = N.take(c,N.arange(0,len(c),2))
        y = N.take(c,N.arange(1,len(c),2))
        xy = N.transpose([x,y])

        p = Polygon(xy,facecolor=facecolor)
        if label != None:
            pylab.text((x.min()+x.max())/2.,(y.min()+y.max())/2.,label,fontsize=6,
                       horizontalalignment='center',
                       verticalalignment='center')
        a = ax.add_artist(p)
        pylab.setp(a, alpha=alpha)
        pylab.setp(a, linestyle=linestyle)
        pylab.setp(a, edgecolor=edgecolor)
        pylab.setp(a, facecolor=facecolor)
        pylab.setp(a, linewidth=linewidth)
        
        allx.append(min(x))
        allx.append(max(x))
        ally.append(min(y))
        ally.append(max(y))
        
    bounds = (min(allx),max(allx),min(ally),max(ally))

    mx = (bounds[1]-bounds[0])*margin
    my = (bounds[3]-bounds[2])*margin

    pylab.xlim(bounds[0]-mx,bounds[1]+mx)
    pylab.ylim(bounds[2]-my,bounds[3]+my)
    return bounds

def plotvert(path,label=None,alpha=.5, linestyle="solid", edgecolor="r", facecolor="red", margin=.05,
               npoly=None, dontplot=False):

    allx,ally = [],[]
    ax = pylab.gca()

    stuff = N.loadtxt(path)
    print "npoly in vert",len(stuff)
    if len(stuff.shape)==1:
        stuff = [stuff]
    
    for k,line in enumerate(stuff):
        x = line[::2]
        y = line[1::2]
        xy = N.transpose([x,y])

        p = Polygon(xy,facecolor=facecolor)
        ax = pylab.gca()

        a = ax.add_artist(p)
        pylab.setp(a, alpha=alpha)
        pylab.setp(a, linestyle=linestyle)
        pylab.setp(a, edgecolor=edgecolor)
        pylab.setp(a, facecolor=facecolor)

        #if k==0:
        for i in range(len(x)):
            pylab.text(x[i],y[i],"%i"%(i+1),fontsize=8)

        

        pylab.xlim(min(x),max(x))
        pylab.ylim(min(y),max(y))
        #pylab.figure()

        
        allx.append(min(x))
        allx.append(max(x))
        ally.append(min(y))
        ally.append(max(y))
        
    bounds = (min(allx),max(allx),min(ally),max(ally))

    mx = (bounds[1]-bounds[0])*margin
    my = (bounds[3]-bounds[2])*margin

    pylab.xlim(bounds[0]-mx,bounds[1]+mx)
    pylab.ylim(bounds[2]-my,bounds[3]+my)
    return bounds


def convert(fname, npoly=None, label=False):
    """ Load a mask file for plotting.
    If the file name ends with .reg it is assumed to be a ds9 region file
    otherwise it should be a file that mangle understands.
    """
    if fname.endswith(".reg"):
        bounds = plotds9reg(fname, npoly=npoly, dontplot=False, label=label)
    elif fname.endswith(".vert"):
        bounds = plotvert(fname)
    else:
        # assume it is a mangle file
        tmp = os.tempnam('.')
        cmd = 'poly2poly -m1e-5s -k1e-10 -og12 %s %s'%(fname,tmp)
        r=os.system(cmd)
        if not r==0:
            print "!! mangle failed (code %i)"%r
            print "dead"
            sys.exit(r)
        bounds=plotgraphics(tmp, npoly=npoly, dontplot=False)
        os.unlink(tmp)
    #else:
    #    return None
    return bounds


def saveimage(manglefile, pngfile, dpi=100, overwrite=False):
    """ """
    if os.path.exists(pngfile):
        if not overwrite:
            return
        
    bounds = convert(manglefile)
    aspect = N.cos((bounds[2]+bounds[3])/2.*N.pi/180)
    print "aspect",aspect
    pylab.figure(aspect=aspect)

    pylab.xlim(pylab.xlim[1],pylab.xlim[0])
    pylab.xlabel("RA")
    pylab.ylabel("Dec")

    title = manglefile.split("/")[-1]
    title = title.split(".")[0]
    pylab.title(title)

    print "wrote file",pngfile
    pylab.savefig(pngfile,dpi=dpi)
    pylab.close()

if __name__=="__main__":

    margin = .05

    parser = argparse.ArgumentParser(description='Plot a Mangle mask file')
    parser.add_argument('--npoly',metavar='N',type=int,default=None,
                        help='only plot the first N polygons in the file')
    parser.add_argument('--save',metavar='maskplot.png',default='maskplot.png',
                        help='The file to save the plot to.')

    parser.add_argument('--label',action='store_true',help='Label polygons from the file name')
    parser.add_argument('maskfile',nargs='*',
                        help='path to mask file.  If no path is given, read paths from STDIN')

    parser.add_argument('--show',action='store_true',
                        help='')

    args = parser.parse_args()

    print args.maskfile
    if args.maskfile == []:
        args.maskfile = []
        for line in sys.stdin.read().split("\n"):
            line = line.strip()
            if line=="": continue
            if line.startswith("#"): continue
            args.maskfile.append(line.strip())

    outname = args.save
    npoly = args.npoly
    label = args.label

    pylab.figure(figsize=(20,20))
    pylab.subplot(111, aspect=1)
    
    allbounds = []
    for f in args.maskfile:
        bounds=convert(f, npoly=npoly, label=label)
        if bounds==None: continue
        allbounds.append(bounds)

    allbounds = N.array(allbounds)
    print allbounds
    low = allbounds.min(axis=0)
    high = allbounds.max(axis=0)
    bounds = [low[0],high[1],low[2],high[3]]
    print bounds
    
    aspect = N.cos((bounds[2]+bounds[3])/2.*N.pi/180)
    print "aspect",aspect
    #pylab.subplot(111,aspect=aspect)

    
    mx = (bounds[1]-bounds[0])*margin
    my = (bounds[3]-bounds[2])*margin


    pylab.xlim(bounds[1]+mx,bounds[0]-mx)  # axis is reversed
    pylab.ylim(bounds[2]-my,bounds[3]+my)


    pylab.xlabel("RA")
    pylab.ylabel("Dec")

    pylab.savefig(outname, dpi=300)
    tmpname = os.tempnam(".")
    os.system('convert -trim %s %s'%(outname,tmpname))
    os.system('mv %s %s'%(tmpname, outname))
    if args.show: pylab.show()
        
