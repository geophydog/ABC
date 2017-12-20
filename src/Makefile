OBJ = abc_egf.o sacio.o

# You should know where the FFTW3 exists

LDLIBS = -L/home/feng_xuping/MY_LIB/lib -lfftw3 -lm
mycorr : $(OBJ)
	cc -o abc_egf $(OBJ) $(LDLIBS)

$(OBJ) : sacio.h

clean : 
	rm -f abc_egf $(OBJ)
