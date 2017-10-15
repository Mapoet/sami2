ARCH=i386
ARCH=amd64
INC_ADD= '$(subst \,/,$(MSMPI_INC));$(subst \,/,$(MSMPI_INC))$(ARCH)'
INC_ADD= '/opt/pgi/linux86-64/2017/mpi/openmpi/include/'


PG_COMPILER=pgfortran
PG_COMPILER=mpifort

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

PG_LOPT  =
#PG_LOPT += -acclibs
#PG_LOPT += -Mprof
#PG_LOPT += -Mconcur
#PG_LOPT += -Mpfi
#PG_LOPT += -Mpfo
PG_LOPT += -module obj
#PG_LOPT += -Mmpi=msmpi






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

