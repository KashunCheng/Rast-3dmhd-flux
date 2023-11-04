#FC     =       ifort
FC	=	mpif90
#FC	=	ftn
# gcc:
FFLAGS  =       -O2 -g -acc -Minfo
FFLAGS_EXTRA =
# intel:
#FFLAGS =       -O2 -g -traceback
#FFLAGS =       -O2 -pad -ip -unroll -align -w -i-static -opt-report
#LMPI   =       -lmpichf90 -lmpichfarg -lmpich -L/lib64 -libt -lpublic -lmpicm -lmtl_common -lvapi -lmpga -lmosal -lpthread
LMPI    =       #-lmpichf90
IDIR    =       #-I/coral/local/mpich64/include
LDIR    =       #-L/coral/local/mpich64/lib64

BIN 	=	 3dmhd.exe
 
SRC 	=	 3dmhd.f 3dmhdset.f 3dmhdsub.f 3dmhd-flux.f90

OBJ 	=	 3dmhd.o 3dmhdset.o 3dmhdsub.o 3dmhd-flux.o

$(BIN): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LDIR) $(LMPI) -o $(BIN)

3dmhd.o: 3dmhdparam.f
3dmhdset.o: 3dmhdparam.f
3dmhdsub.o: 3dmhdparam.f
3dmhd-flux.o: 3dmhdparam.f 3dmhd-flux.f90
	$(FC) $(FFLAGS) -c $(IDIR) 3dmhd-flux.f90

clean:
	rm -rf *.o $(BIN)

.f.o:
	$(FC) $(FFLAGS) $(FFLAGS_EXTRA) -c $(IDIR) $*.f
