COMPILER=gfortran
LIBS = -llapack -lfftw3

CFLAGS = -std=gnu -fbounds-check -Wall -Wuninitialized

FILES = myKinds.o fftw3constants.o fft.o pa.o hi.o

####################################################################################################
# Pattern rules
%.o : %.f95
	$(COMPILER) $(CFLAGS) -c $<

%.o : %.f90
	$(COMPILER) $(CFLAGS) -c $<

%.o : %.f
	$(COMPILER) $(CFLAGS) -c $<

####################################################################################################
# targets

prog: $(FILES)
	$(COMPILER) $(CFLAGS) $(FILES) $(LIBS) -o $@

new:	clean prog

# directly run program after compilation
run:	prog
	./prog

# directly run program after compilation
all:	run

####################################################################################################
# special targets
.PHONY: clean new cleanout

####################################################################################################
# Dependencies
myKinds.o:
fftw3constants.o:
fft.o: myKinds.o fftw3constants.o
pa.o: myKinds.o
prog.o: myKinds.o fft.o pa.o

####################################################################################################
# clean up rules
cleanout: clean
	rm -f out.*

clean:
	rm -f *.o *.mod

####################################################################################################
