#!/usr/bin/python3

import sys

CRED = '\033[91m'
CEND = '\033[0m'

if (len(sys.argv)!=2 and len(sys.argv)!=3) :
    print(CRED, 'Please, specify either the redshift range [2 floats] or the reference work [1 string]', CEND)
    
else:
    file = 'fsigma8.dat'

    with open(file) as fin:
        data = list(zip(*[line.split() for line in fin]))

    if (len(sys.argv)==2):
        check = False
        for ll in range(len(data[0]))[1:]:
            if (data[0][ll]==sys.argv[1]):
                check = True
                print(data[0][ll], 'f*sigma8( z=', data[1][ll], ') =', data[2][ll], 'with errors -', data[3][ll], ' +', data[4][ll])
        if (check==False):
            print(CRED, '\nNo fsigma8 measurements are stored in the table related to the chosen reference. Please, specify one of the following:\n', CEND)
            print(data[0][1:])
        
    if (len(sys.argv)==3):
        check = False
        for ll in range(len(data[0]))[1:]:
            if (float(sys.argv[1])<float(data[1][ll]) and float(data[1][ll])<float(sys.argv[2])):
                check = True
                print(data[0][ll], 'f*sigma8( z=', data[1][ll], ') =', data[2][ll], 'with errors -', data[3][ll], ' +', data[4][ll])
        if (check==False):
            print(CRED, '\nNo fsigma8 measurements are stored in the table in the choosen redshift range\n', CEND)
