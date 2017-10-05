
PG_COMPILER=pgfortran

PG_COPT  =
PG_COPT += -acc=autopar
#PG_COPT += -Minfo
PG_COPT += -Mprof
PG_COPT += -Mconcur
PG_COPT += -Mpfi
PG_COPT += -Mpfo
#PG_COPT += -fast
PG_COPT += -module obj 

PG_LOPT  =
PG_LOPT += -acclibs
PG_LOPT += -Mprof
PG_LOPT += -Mconcur
PG_LOPT += -Mpfi
PG_LOPT += -Mpfo






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
build/sami2.x: obj/sami2-1.00.o obj/grid-1.00.o obj/chapman.o obj/nrlmsise00.o obj/hwm93.o
	$(FC) $(FCLOPT) -o $@ $^
	cp $@ build/sami2.exe

obj build:
	mkdir $@
build/input:input
	cp -r $^ $@
#input:deni-init.inp euvflux.inp ichem.inp phabsdt.inp phiondt.inp
build/sami2-1.00.namelist:sami2-1.00.namelist
	cp $^ $@
obj/sami2-1.00.o:|obj obj/parameters.mod obj/commons.mod
obj/sami2-1.00.o:sami2-1.00.f90 com-1.00.inc param-1.00.inc gonamelist.inc
	$(FC) $(FCCOPT) -c -o $@ $<
#sami2-1.00.f90:sami2-1.00.f
#	cp $^ $@

obj/parameters.mod:obj/param-1.00.o
obj/param-1.00.o:|obj
obj/param-1.00.o:param-1.00.f90
	$(FC) $(FCCOPT) -c -o $@ $<

obj/commons.mod:obj/com-1.00.o
obj/com-1.00.o:|obj
obj/com-1.00.o:com-1.00.f90
	$(FC) $(FCCOPT) -c -o $@ $<

obj/grid-1.00.o:|obj
obj/grid-1.00.o:grid-1.00.f90 com-1.00.inc param-1.00.inc
	$(FC) $(FCCOPT) -c -o $@ $<
#grid-1.00.f90:grid-1.00.f
#	cp $^ $@


obj/chapman.o:|obj
obj/chapman.o:chapman.f90
	$(FC) $(FCCOPT) -c -o $@ $^
#chapman.f90:chapman.f
#	cp $^ $@


obj/nrlmsise00.o:|obj
obj/nrlmsise00.o:nrlmsise00.f90
	$(FC) $(FCCOPT) -c -o $@ $^
#nrlmsise00.f90:nrlmsise00.f
#	cp $^ $@

obj/hwm93.o:|obj
obj/hwm93.o:hwm93.f90
	$(FC) $(FCCOPT) -c -o $@ $^
#hwm93.f90:hwm93.f
#	cp $^ $@

