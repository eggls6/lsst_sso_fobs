#Makefile for p-type hierarchical triple simualtion in 2d 


#compiler
cf=gfortran
#compilder flags
cflags=-O  -fbounds-check #-fopenmp 
#program output name
pg=lsst_mba.exe

modk=kind_m.o
mod0=const_m.o
mod1=trafo3_m.o
mod2=photometry_m.o
mod3=camera_m.o
mod4=orbdyn2_m.o
mod5=fileio_m.o
mod6=earth_m.o
mod7=spicelib.a
mod8=support.a
mod9=spice_m.o
mod10=lsst_mba.o


obj= $(modk) $(mod0) $(mod1) $(mod2) $(mod3) $(mod4) $(mod5) $(mod6) $(mod7) $(mod8) $(mod9) $(mod10)
targ:$(obj)
	$(cf) ${cflags} -o $(pg) $(obj)
	
%.o: %.f90
	$(cf) ${cflags} -c $<

%.o: %.f
	$(cf) ${cflags} -c $<
    
# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

