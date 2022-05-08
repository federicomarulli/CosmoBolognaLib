F90C = gfortran
switch = -O3

OBJ1  = read_pk_sum.o
OBJ2  = fit_pkred.o
OBJ3  = read_pk2.o
OBJ4  = stdPT2.o
OBJ5  = fit_pkred3.o
OBJ6  = LPTred.o
OBJ7  = LPTreal.o
OBJ8  = read_pk_real.o
OBJ9  = stdPTred.o
OBJ10  = stdPTreal.o
OBJ11  = fit_pkred4.o
OBJ12  = calc_pk_correction.o
OBJ13  = calc_pk_correction2.o

all: read_pk_sum fit_pkred read_pk2 stdPT2 fit_pkred3 LPTred LPTreal read_pk_real stdPTred stdPTreal fit_pkred4 calc_pk_correction calc_pk_correction2

read_pk_sum: $(OBJ1) 
	$(F90C) $(OBJ1) -o read_pk_sum $(FLAGS_LIB)

fit_pkred: $(OBJ2) 
	$(F90C) $(OBJ2) -o fit_pkred $(FLAGS_LIB) 

read_pk2: $(OBJ3) 
	$(F90C) $(OBJ3) -o read_pk2 $(FLAGS_LIB)

stdPT2: $(OBJ4) 
	$(F90C) $(OBJ4) -o stdPT2 $(FLAGS_LIB) 

fit_pkred3: $(OBJ5) 
	$(F90C) $(OBJ5) -o fit_pkred3 $(FLAGS_LIB) 

LPTred: $(OBJ6) 
	$(F90C) $(OBJ6) -o LPTred $(FLAGS_LIB) 

LPTreal: $(OBJ7) 
	$(F90C) $(OBJ7) -o LPTreal $(FLAGS_LIB) 

read_pk_real: $(OBJ8) 
	$(F90C) $(OBJ8) -o read_pk_real $(FLAGS_LIB) 

stdPTred: $(OBJ9) 
	$(F90C) $(OBJ9) -o stdPTred $(FLAGS_LIB) 

stdPTreal: $(OBJ10) 
	$(F90C) $(OBJ10) -o stdPTreal $(FLAGS_LIB) 

fit_pkred4: $(OBJ11) 
	$(F90C) $(OBJ11) -o fit_pkred4 $(FLAGS_LIB) 

calc_pk_correction: $(OBJ12) 
	$(F90C) $(OBJ12) -o calc_pk_correction $(FLAGS_LIB) 

calc_pk_correction2: $(OBJ13) 
	$(F90C) $(OBJ13) -o calc_pk_correction2 $(FLAGS_LIB) 

clean:
	rm -f *.o read_pk_sum fit_pkred read_pk2 stdPT2 fit_pkred3 LPTred LPTreal read_pk_real stdPTred stdPTreal fit_pkred4 calc_pk_correction calc_pk_correction2 *~ \#* temp* core*

read_pk_sum.o: read_pk_sum.f makefile 
	$(F90C) -c $(switch) read_pk_sum.f

fit_pkred.o: fit_pkred.f makefile 
	$(F90C) -c $(switch) fit_pkred.f

read_pk2.o: read_pk2.f makefile 
	$(F90C) -c $(switch) read_pk2.f

stdPT2.o: stdPT2.f makefile 
	$(F90C) -c $(switch) stdPT2.f

fit_pkred3.o: fit_pkred3.f makefile 
	$(F90C) -c $(switch) fit_pkred3.f

LPTred.o: LPTred.f makefile 
	$(F90C) -c $(switch) LPTred.f

LPTreal.o: LPTreal.f makefile 
	$(F90C) -c $(switch) LPTreal.f

read_pk_real.o: read_pk_real.f makefile 
	$(F90C) -c $(switch) read_pk_real.f

stdPTred.o: stdPTred.f makefile 
	$(F90C) -c $(switch) stdPTred.f

stdPTreal.o: stdPTreal.f makefile 
	$(F90C) -c $(switch) stdPTreal.f

fit_pkred4.o: fit_pkred4.f makefile 
	$(F90C) -c $(switch) fit_pkred4.f

calc_pk_correction.o: calc_pk_correction.f makefile 
	$(F90C) -c $(switch) calc_pk_correction.f

calc_pk_correction2.o: calc_pk_correction2.f makefile 
	$(F90C) -c $(switch) calc_pk_correction2.f

