ifeq ($(OS),Windows_NT)
    CCFLAGS += -D WIN32
    OSFLAG= WIN32
    ifeq ($(PROCESSOR_ARCHITEW6432),AMD64)
        CCFLAGS += -D AMD64
	   ARCH=amd64
    else
        ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
            CCFLAGS += -D AMD64
		  ARCH=amd64
        endif
        ifeq ($(PROCESSOR_ARCHITECTURE),x86)
            CCFLAGS += -D IA32
		  ARCH=i386
        endif
    endif
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        CCFLAGS += -D LINUX
	   OSFLAG= LINUX
    endif
    ifeq ($(UNAME_S),Darwin)
        CCFLAGS += -D OSX
	   OSFLAG= OSX
    endif
    UNAME_P := $(shell uname -p)
    ifeq ($(UNAME_P),x86_64)
        CCFLAGS += -D AMD64
    endif
    ifneq ($(filter %86,$(UNAME_P)),)
        CCFLAGS += -D IA32
    endif
    ifneq ($(filter arm%,$(UNAME_P)),)
        CCFLAGS += -D ARM
    endif
endif

PG_COMPILER=
PG_COMPILER=pgfortran
ifeq ($(OSFLAG),WIN32)
INC_ADD= '$(subst \,/,$(MSMPI_INC));$(subst \,/,$(MSMPI_INC))$(ARCH)'
else
INC_ADD= '/opt/pgi/linux86-64/2017/mpi/openmpi/include/'
PG_COMPILER=mpifort
endif





PG_COPT  =
#PG_COPT += -acc=autopar
#PG_COPT += -Minfo
#PG_COPT += -Mprof
#PG_COPT += -Mconcur
#PG_COPT += -Mpfi
#PG_COPT += -Mpfo
#PG_COPT += -fast
PG_COPT += -module obj
PG_COPT += -I$(INC_ADD)
PG_COPT += -cpp
PG_COPT += -D _USE_MPI_

PG_COPT += -Mbounds
PG_COPT += -Minfo=all
PG_COPT += -traceback
PG_COPT += -Mchkfpstk
PG_COPT += -Mchkstk
PG_COPT += -Mdalign
##PG_COPT +=-Mdclchk 
PG_COPT +=-Mdepchk 
PG_COPT +=-Miomutex 
PG_COPT +=-Mrecursive 
PG_COPT +=-Msave 
PG_COPT +=-Ktrap=align
PG_COPT +=-Ktrap=denorm
PG_COPT +=-Ktrap=divz
PG_COPT +=-Ktrap=fp
PG_COPT +=-Ktrap=inexact
PG_COPT +=-Ktrap=inv
PG_COPT +=-Ktrap=none
#PG_COPT +=-Ktrap=unf
#PG_COPT +=-Ktrap=ovf
PG_COPT +=-g 
PG_COPT +=-byteswapio
PG_COPT +=-Kieee

PG_LOPT  =
PG_LOPT += -acclibs
PG_LOPT += -Mprof
PG_LOPT += -Mconcur
PG_LOPT += -Mpfi
PG_LOPT += -Mpfo
PG_LOPT += -module obj

ifeq ($(OSFLAG),WIN32)
PG_LOPT += -Mmpi=msmpi
else
#PG_LOPT += -Mmpi
endif






FC=$(PG_COMPILER)


FCCOPT = $(PG_COPT)
FCLOPT = $(PG_LOPT)


default:

.PHONY:default
default: sami2.x
.PHONY: clean
clean:
	rm -rf obj
	rm -rf build
.PHONY:sami2.x
sami2.x: build/sami2.x


build/sami2.x:|build
build/sami2.x:|build/input
build/sami2.x: obj/sami2-1.00.o obj/grid-1.00.o obj/chapman.o obj/nrlmsise00.o obj/hwm93.o obj/com-1.00.o obj/com-subroutines.o obj/vdrift_model.o obj/mpi_client.o obj/param-1.00.o obj/com-1.00.o
	$(FC) $(FCLOPT) -o $@ $^
	cp $@ build/sami2.exe

obj build:
	mkdir $@
build/input:input
	cp -r $^ $@
#input:deni-init.inp euvflux.inp ichem.inp phabsdt.inp phiondt.inp
build/input/sami2-1.00.namelist:input/sami2-1.00.namelist
	cp -f $^ $@
obj/sami2-1.00.o:|obj obj/parameters.mod obj/commons.mod obj/inputfiles.mod obj/commonsubroutines.mod obj/vdrift_model.mod obj/mpi_client.mod
obj/sami2-1.00.o:sami2-1.00.f90 com-1.00.inc param-1.00.inc gonamelist.inc
	$(FC) $(FCCOPT) -c -o $@ $<

obj/parameters.mod:obj/param-1.00.o
obj/param-1.00.o:|obj
obj/param-1.00.o:param-1.00.f90
	$(FC) $(FCCOPT) -c -o $@ $<



obj/mpi_client.mod:obj/mpi_client.o
obj/mpi_client.o:|obj
obj/mpi_client.o:mpi_client.f90
	$(FC) $(FCCOPT) -c -o $@ $<

obj/vdrift_model.mod:obj/vdrift_model.o
obj/vdrift_model.o:|obj
obj/vdrift_model.o:vdrift_model.f90
	$(FC) $(FCCOPT) -c -o $@ $<

obj/commonsubroutines.mod:obj/com-subroutines.o
obj/com-subroutines.o:|obj
obj/com-subroutines.o:com-subroutines.f90
	$(FC) $(FCCOPT) -c -o $@ $<

obj/inputfiles.mod:obj/com-1.00.o
obj/commons.mod:obj/com-1.00.o
obj/com-1.00.o:|obj
obj/com-1.00.o:com-1.00.f90
	$(FC) $(FCCOPT) -c -o $@ $<

obj/grid-1.00.o:|obj
obj/grid-1.00.o:grid-1.00.f90 com-1.00.inc param-1.00.inc
	$(FC) $(FCCOPT) -c -o $@ $<


obj/chapman.o:|obj
obj/chapman.o:chapman.f90
	$(FC) $(FCCOPT) -c -o $@ $^


obj/nrlmsise00.o:|obj
obj/nrlmsise00.o:nrlmsise00.f90
	$(FC) $(FCCOPT) -c -o $@ $^

obj/hwm93.o:|obj
obj/hwm93.o:hwm93.f90
	$(FC) $(FCCOPT) -c -o $@ $^

