#!/usr/bin/python

import os
import glob

os.remove("CBL.h")

fileCBL = open("CBL", "w")

for file in glob.glob("*.h"):
    fileCBL.write('#include "'+file+'"\n')

fileCBL.close()

os.rename("CBL", "CBL.h")
