
#default settings for ifort

F90C     = mpiifort
#F90C    = ifort

HEALPIX =  /home/roger/Downloads/Healpix_3.31
healpix = $(HEALPIX)
LAPACKL = -mkl=sequential -lmpi -lhealpix -openmp

#remove -xHost if cluster is not homogeneous
#add -DHEALPIXI4B if using older healpix and get errors about arguments not matching
FFLAGS = -O3 -xHost -ip -fpp -error-limit 500 -DMPIPIX -DMPI -heap-arrays -g -traceback

#cfitsio = /usr/local/Cluster-Apps/cfitsio/intel/3.300
cfitsio = /usr/lib/x86_64-linux-gnu
cfitsio ?= $(CFITSIO)


F90FLAGS = $(FFLAGS) -I$(healpix)/include -L$(cfitsio) -I$(healpix)/lib -L$(healpix)/lib $(LAPACKL) -lcfitsio 

OBJFILES= toms760.o inifile.o utils.o spin_alm_tools.o \
   HealpixObj.o HealpixVis.o 

LENSPIX = $(OBJFILES) SimLens.o

PLOT = $(OBJFILES) plot.o


default: simlens
all: simlens recon plot

spin_alm_tools.o:  utils.o toms760.o
HealpixObj.o: spin_alm_tools.o
HealpixVis.o: HealpixObj.o
SimLens.o: HealpixVis.o inifile.o
plot.o: plot.o

.f.o:
	f77 $(F90FLAGS) -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $*.F90


simlens: $(LENSPIX) 	
	$(F90C) -o simlens $(LENSPIX) $(F90FLAGS) $(LINKFLAGS)

recon: $(OBJFILES) LensReconExample.o
	$(F90C) -o recon $(OBJFILES) LensReconExample.o $(F90FLAGS) $(LINKFLAGS)

plot: $(PLOT) 
	$(F90C) -o plot $(PLOT) $(F90FLAGS) $(LINKFLAGS)



clean:
	rm -f *.o* *.e* *.mod *.d *.pc *.obj core* *.il
