#$Id: Makefile.stable_netcdf,v 1.7 2003/03/28 09:27:56 kbk Exp $
#
# Master Makefile for making the GOTM executable.
#

.SUFFIXES:
.SUFFIXES: .f90

SHELL	= /bin/sh

## BEGIN CUSTOM BUILD (FOR GOTM SOSSTA PROJECT, WT 2017-06-29)

# Remember the current directry
CWD = $(CURDIR)
# Remember today's date
TODAY = $(shell date +'%Y%m%d')

# Set a temporary directory to speed up, but defaults to the current directory.
TMP = /dev/shm
ifneq (,$(wildcard $(TMP)))
TMP = .
endif

# Prevent the object archives to be created in deterministic mode.
ARFLAGS=rU

# Use a locally built netcdf-3.6.2 static library. 
# Prefer Intel Fortran over GNU Fortran compiler (if available)

# Use ifort if available.
ifneq (,$(shell which ifort)) # if `which ifort` does not return nothing, assume that's the right one to use
FC=ifort
FFLAGS=-assume bscc # To use '\n', '\r' etc...
else
# Fallback to gfortran if ifort is unavailable.
FC=gfortran
FFLAGS=-fbackslash # To use '\n', '\r' etc...
endif
NETCDFLIB=libnetcdf.a # Build it and copy it to our working folder

# The build options below are all irrelevant for usual PC, Westgrid and ComputeCanada servers.

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
libutil.a(util/time.o)			\
libutil.a(util/tridiagonal.o)		\
libutil.a(util/eqstate.o)			\
libobservations.a(observations/observations.o)	\
libmeanflow.a(meanflow/meanflow.o)		\
libturbulence.a(turbulence/turbulence.o)		\
liboutput.a(output/asciiout.o)			\
liboutput.a(output/ncdfout.o)			\
liboutput.a(output/output.o)			\
libturbulence.a(turbulence/sediment.o)		\
libturbulence.a(turbulence/seagrass.o)

#WT Other than the key module defining source code above,
# The rest of the code under each library is to be grouped
# into subdirectories.

AIRSEA	= \
libairsea.a(airsea/airsea.o)		\
libairsea.a(airsea/short_wave_radiation.o)	\

UTIL	= \
libutil.a(util/advection.o)		\
libutil.a(util/w_split_it_adv.o)	\
libutil.a(util/gridinterpol.o)	\
libutil.a(util/yevol.o)

MEANFLOW	= \
libmeanflow.a(meanflow/updategrid.o)		\
libmeanflow.a(meanflow/adaptivegrid.o)		\
libmeanflow.a(meanflow/coriolis.o)		\
libmeanflow.a(meanflow/uequation.o)		\
libmeanflow.a(meanflow/vequation.o)		\
libmeanflow.a(meanflow/extpressure.o)		\
libmeanflow.a(meanflow/intpressure.o)		\
libmeanflow.a(meanflow/friction.o)		\
libmeanflow.a(meanflow/temperature.o)		\
libmeanflow.a(meanflow/salinity.o)		\
libmeanflow.a(meanflow/stratification.o)		\
libmeanflow.a(meanflow/buoyancy.o)		\
libmeanflow.a(meanflow/convectiveadjustment.o)	\
libmeanflow.a(meanflow/production.o)

TURBULENCE   = \
libturbulence.a(turbulence/tkeeq.o)		\
libturbulence.a(turbulence/q2over2eq.o)		\
libturbulence.a(turbulence/lengthscaleeq.o)	\
libturbulence.a(turbulence/dissipationeq.o)	\
libturbulence.a(turbulence/genericeq.o)		\
libturbulence.a(turbulence/tkealgebraic.o)		\
libturbulence.a(turbulence/algebraiclength.o)	\
libturbulence.a(turbulence/ispralength.o)		\
libturbulence.a(turbulence/potentialml.o)		\
libturbulence.a(turbulence/cmue_bb.o)		\
libturbulence.a(turbulence/cmue_bbqe.o)		\
libturbulence.a(turbulence/cmue_ca.o)		\
libturbulence.a(turbulence/cmue_caqe.o)		\
libturbulence.a(turbulence/cmue_cb.o)		\
libturbulence.a(turbulence/cmue_cbqe.o)		\
libturbulence.a(turbulence/cmue_kc.o)		\
libturbulence.a(turbulence/cmue_kcqe.o)		\
libturbulence.a(turbulence/cmue_my.o)		\
libturbulence.a(turbulence/cmue_gpqe.o)		\
libturbulence.a(turbulence/cmue_ma.o)		\
libturbulence.a(turbulence/cmue_sg.o)		\
libturbulence.a(turbulence/cmue_rf.o)		\
libturbulence.a(turbulence/fk_craig.o)		\
libturbulence.a(turbulence/turbulence_adv.o)	\
libturbulence.a(gotm.o) #WT Curious. A particular module depends on main loop gotm.f90!?
#libturbulence.a(gotm_lib_version.o)

OBSERVATIONS   = \
libobservations.a(observations/analytical_profile.o)	\
libobservations.a(observations/get_eps_profile.o)	\
libobservations.a(observations/get_ext_pressure.o)	\
libobservations.a(observations/get_int_pressure.o)	\
libobservations.a(observations/get_s_profile.o)	\
libobservations.a(observations/get_t_profile.o)	\
libobservations.a(observations/get_vel_profile.o)	\
libobservations.a(observations/get_w_adv.o)	\
libobservations.a(observations/get_zeta.o)	\
libobservations.a(observations/read_extinction.o)	\
#libobservations.a(observations/read_chlo.o)	\


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

backup: gotm
	tar -czvf .backup/gotm-src-$(TODAY).tar.gz *.f90 Makefile

backup-all: gotm
	tar -czvf .backup/gotm-$(TODAY).tar.gz *

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
	@cd $(TMP) && \
	rm -Rf netcdf-3.6.2 && \
	tar -xzf $(CWD)/.backup/netcdf-3.6.2.tar.gz && \
	cd ./netcdf-3.6.2 && \
	echo "Configuring netCDF-3.6.2..." && \
	./configure --silent --disable-utilities --disable-v2 --disable-examples --disable-cxx && \
	echo "Running GNU make to create libnetcdf.a..." && \
	make -j --silent 2>&1 > make_netCDF-3.6.2.log && \
	cp ./libsrc/.libs/libnetcdf.a $(CWD) && \
	cd $(CWD) && \
	rm -Rf $(TMP)/netcdf-3.6.2

# Do not remove libnetcdf.a which takes a long time to build.
clean:	
	-mv libnetcdf.a libnetcdf.a.backup
	-rm -f lib*.a  *.mod *.o
	-mv libnetcdf.a.backup libnetcdf.a

# Recompile but only leave the GOTM executable in the folder
gotm-only: gotm clean

# Remove also libnetcdf.a and force a recompilation of the big library. Necessary if you're switching compiler, for instance.
realclean: clean
	-rm libnetcdf.a
	-rm -f gotm 

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
