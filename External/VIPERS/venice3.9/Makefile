#------------------------------------------------#
#Makefile for venice                             #	
#------------------------------------------------#

CC             = icc -use-asm
FITS           = yes
SRC            = main.c
LDFLAGS        = -L$(HOME)/local/lib -lgsl -lgslcblas -lm
#LDFLAGS        = -L/sw64/lib/ -lgsl -lgslcblas -lm
CFLAGS_PYTHON  = -I/usr/include/python2.6
LDFLAGS_PYTHON = -ldl -lpython2.6

ifeq ($(FITS),yes)
    LDFLAGS += -lcfitsio
    SRC     += fits.c
else
   SRC     += withoutFits.c
endif

CFLAGS	= -I$(HOME)/local/include -fPIC #-Wall -Wuninitialized -O3 
EXEC	= venice
OBJ	= $(SRC:.c=.o)

all: $(EXEC)

$(EXEC) : $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	rm -rf *.o $(EXEC).tar *.so $(EXEC)py.py $(EXEC)py.pyc $(SRC:.c=_wrap.c)

mrproper: clean
	rm -rf $(EXEC)

tar: 
	tar cvf $(EXEC).tar Makefile *.c *.h README

python:
	swig -python main.i	
	$(CC) -c $(SRC) main_wrap.c $(CFLAGS) $(CFLAGS_PYTHON)
	$(CC) -bundle $(SRC:.c=.o) main_wrap.o -o _$(EXEC)py.so $(LDFLAGS) $(LDFLAGS_PYTHON)
