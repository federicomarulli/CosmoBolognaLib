## Makefile for FFTLog-f90 (https://github.com/coccoinomane/fftlog-f90)

# Settings for gfortran
F90 = gfortran -g
F90FLAGS = -O3 -ftree-vectorize -ffast-math -fimplicit-none 

# Settings for ifort by Intel
# F90 = ifort -g
# F90FLAGS = -fast -par-report0 -vec-report0 -parallel -xsse4.1 -w -arch i386

COMPILE = $(F90) $(F90FLAGS)

SOURCES = fftlog_driver.f90 drffti.f drfftb.f drfftf.f fftlog.f cdgamma.f spline.f90 

OBJECTS = $(SOURCES:.f=.o)

EXECUTABLE = fftlog-f90

all : $(SOURCES) $(EXECUTABLE)

%.o : %.f
	$(COMPILE) -c $< -o $@

$(EXECUTABLE) : $(OBJECTS)
	$(COMPILE) $(OBJECTS) -o $@

clean :
	rm -f *.o
