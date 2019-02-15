#!/bin/bash
# Author: Johannes Buchner (C) 2015
# For creating a shared library (libcuba.so)

echo $1
echo "rebuilding libcuba.a archive"
./configure
sed 's/CFLAGS = -O3 -fomit-frame-pointer/CFLAGS = -O3 -fPIC -fomit-frame-pointer/g' --in-place makefile
make -B libcuba.a
echo "unpacking libcuba.a"

FILES=$(ar xv libcuba.a |sed 's/x - //g')
ES="so"
FLAGS_LINK=-shared

LIB=libcuba.$ES

echo "making libcuba."$ES
gcc $FLAGS_LINK $FILES -o $LIB
rm $FILES
