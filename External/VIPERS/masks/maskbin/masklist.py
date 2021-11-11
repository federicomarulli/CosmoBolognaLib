#!/usr/bin/python
"""masklist.py
Time-stamp: <2011-05-06 16:43:00 ben>
Ben Granett ben.granett@brera.inaf.it
"""
import sys
import argparse


parser = argparse.ArgumentParser(description='List fields')
parser.add_argument('-f',metavar='field',type=type('string'),default='W1',
                    choices=['W1','W4'],
                    help='either W1 or W4')

parser.add_argument('-b',metavar='base',default='regions/mask_%s.reg',
                    help='')

parser.add_argument('-1',dest='one',action='store_true',
                    help='')

args=parser.parse_args()
pointinglist = '/home/ben/vipersmask/data_1.1/pointings_1.1.txt'


base = args.b
field = args.f
one = args.one

quads = ['Q1','Q2','Q3','Q4']

def listmasks(field,base, pointinglist=pointinglist):
    """ """
    paths = []
    for pointing in file(pointinglist):
        pointing = pointing.strip()
        if pointing.startswith("#"): continue
        if not pointing.startswith(field): continue

        for q in quads:
            print base%(pointing+q),
            if not one:
                print ""
    #print " ".join(paths)
            

listmasks(field=field, base=base)
