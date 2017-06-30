#$Id: Makefile.stable_netcdf,v 1.7 2003/03/28 09:27:56 kbk Exp $
#
# Master Makefile for making the GOTM executable.
#

.SUFFIXES:
.SUFFIXES: .f90

SHELL	= /bin/sh

## BEGIN CUSTOM BUILD (FOR GOTM SOSSTA PROJECT, WT 2017-06-29)

# Use a locally built netcdf-3.6.2 static library. 
# Prefer Intel Fortran over GNU Fortran compiler (if available)

# Use ifort if available.
ifneq (,$(shell which ifort)) # if `which ifort` does not return nothing, assume that's the right one to use
FC=ifort
else
# Fallback to gfortran if ifort is unavailable.
FC=gfortran
endif
NETCDFLIB = libnetcdf.a # Build it and copy it to our working folder

# The following build options are all irrelevant.

## END CUSTOM BUILD

# Set up the path to the NetCDF library - very important.

#NETCDFINCDIR	= ./netcdf-3.6.2/libsrc/
#/users/zxb/netcdf/netcdf-350/src/f90/netcdf.mod
#NETCDFINCDIR = /usr/include

#NAG F95 compiler for Linux
#FC     = f95nag
#FFLAGS = -f77 -I$(NETCDFINCDIR)

#Fujitsu F95 compiler for Linux
#FC     = f95
#FFLAGS = -Am -fw -I$(NETCDFINCDIR)

#DECFOR F95 compiler for OSF1 - alpha
#FC    = f95
#FFLAGS = -I$(NETCDFINCDIR)

#Intel Fortran  compiler for 
#FC=ifort
#FFLAGS=-g
#FFLAGS=-g -fbacktrace -ffpe-trap=zero,overflow,underflow
#FFLAGS=-g -fbacktrace -ffpe-trap=overflow,underflow

#FC=gfortran
#FFLAGS=-g -fbacktrace
## No more need.
#FFLAGS = -I$(NETCDFINCDIR)

#Mac Fortran  compiler for 
#FC    = gfortran
#FFLAGS = -I$(NETCDFINCDIR)

#sun fortran compiler
#FC = f90
#FFLAGS = -I$(NETCDFINCDIR)

MODULES	= \
libutil.a(time.o)			\
libutil.a(tridiagonal.o)		\
libutil.a(eqstate.o)			\
libobservations.a(observations.o)	\
libmeanflow.a(meanflow.o)		\
libairsea.a(airsea.o)			\
libturbulence.a(turbulence.o)		\
liboutput.a(asciiout.o)			\
liboutput.a(ncdfout.o)			\
liboutput.a(output.o)			\
libturbulence.a(sediment.o)		\
libturbulence.a(seagrass.o)


UTIL	= \
libutil.a(advection.o)		\
libutil.a(w_split_it_adv.o)	\
libutil.a(gridinterpol.o)	\
libutil.a(yevol.o)

MEANFLOW	= \
libmeanflow.a(updategrid.o)		\
libmeanflow.a(adaptivegrid.o)		\
libmeanflow.a(coriolis.o)		\
libmeanflow.a(uequation.o)		\
libmeanflow.a(vequation.o)		\
libmeanflow.a(extpressure.o)		\
libmeanflow.a(intpressure.o)		\
libmeanflow.a(friction.o)		\
libmeanflow.a(temperature.o)		\
libmeanflow.a(salinity.o)		\
libmeanflow.a(stratification.o)		\
libmeanflow.a(buoyancy.o)		\
libmeanflow.a(convectiveadjustment.o)	\
libmeanflow.a(production.o)

TURBULENCE   = \
libturbulence.a(tkeeq.o)		\
libturbulence.a(q2over2eq.o)		\
libturbulence.a(lengthscaleeq.o)	\
libturbulence.a(dissipationeq.o)	\
libturbulence.a(genericeq.o)		\
libturbulence.a(tkealgebraic.o)		\
libturbulence.a(algebraiclength.o)	\
libturbulence.a(ispralength.o)		\
libturbulence.a(potentialml.o)		\
libturbulence.a(cmue_bb.o)		\
libturbulence.a(cmue_bbqe.o)		\
libturbulence.a(cmue_ca.o)		\
libturbulence.a(cmue_caqe.o)		\
libturbulence.a(cmue_cb.o)		\
libturbulence.a(cmue_cbqe.o)		\
libturbulence.a(cmue_kc.o)		\
libturbulence.a(cmue_kcqe.o)		\
libturbulence.a(cmue_my.o)		\
libturbulence.a(cmue_gpqe.o)		\
libturbulence.a(cmue_ma.o)		\
libturbulence.a(cmue_sg.o)		\
libturbulence.a(cmue_rf.o)		\
libturbulence.a(fk_craig.o)		\
libturbulence.a(turbulence_adv.o)	\
libturbulence.a(gotm.o)
#libturbulence.a(gotm_lib_version.o)

OBSERVATIONS   = \
libobservations.a(analytical_profile.o)	\
libobservations.a(get_eps_profile.o)	\
libobservations.a(get_ext_pressure.o)	\
libobservations.a(get_int_pressure.o)	\
libobservations.a(get_s_profile.o)	\
libobservations.a(get_t_profile.o)	\
libobservations.a(get_vel_profile.o)	\
libobservations.a(get_w_adv.o)	\
libobservations.a(get_zeta.o)	\
libobservations.a(read_extinction.o)	\
libobservations.a(read_chlo.o)	\


LIBS	=	libairsea.a		\
		libturbulence.a 	\
		libmeanflow.a 		\
		libobservations.a	\
		liboutput.a		\
		libutil.a		

all: gotm 

gotm: $(NETCDFLIB) $(MODULES) $(LIBS) main.o
	$(FC) main.o -o $@ gotm.o $(LIBS) $(NETCDFLIB)
	-rm main.o

main.o: gotm.o

modules: $(MODULES)

libairsea.a: $(AIRSEA)

libturbulence.a: $(TURBULENCE)

libutil.a: $(UTIL)

libmeanflow.a: $(MEANFLOW)

libobservations.a: $(OBSERVATIONS)

# The following commands needs to be chained together, and must end with
# backslash for MAKE to know it is one single command.
#
# STEPS:
# 1. Remove the existing source folder if any.
# 2. Extract the archive
# 3. Go inside the subfolder
# 4. Make only for running Fortran code
# 5. Make with unlimited number of threads.
# 6. Copy our static library
# 7. Go back to parent folder
# 8. Remove the source
libnetcdf.a: 
	@rm -Rf netcdf-3.6.2 && \
	tar -xzf netcdf-3.6.2.tar.gz && \
	cd ./netcdf-3.6.2 && \
	echo "Configuring netCDF-3.6.2..." && \
	./configure --silent --disable-utilities --disable-v2 --disable-examples --disable-cxx --disable-f90 && \
	echo "Running GNU make to create libnetcdf.a..." && \
	make -j --silent 2>&1 > make_netCDF-3.6.2.log && \
	cp ./libsrc/.libs/libnetcdf.a .. && \
	cd .. && \
	rm -Rf netcdf-3.6.2

# Do not remove libnetcdf.a which takes a long time to build.
clean:	
	-mv libnetcdf.a libnetcdf.a.backup && \
	rm -f lib*.a  *.mod *.o && \
	mv libnetcdf.a.backup libnetcdf.a

realclean: clean
	-rm libnetcdf.a
	-rm -f gotm 

%.o: %.f90
	echo "Using $(FC)..."
	$(FC) $(FFLAGS) -c $< -o $@