
PG_COMPILER=pgfortran

PG_COPT  =
PG_COPT += -Minfo -fast

PG_LOPT  =
PG_LOPT +=-acclibs





FC=$(PG_COMPILER)


FCCOPT = $(PG_COPT)
FCLOPT = $(PG_LOPT)


default:

.PHONY:default
default: obj/chapman.o


obj:
	mkdir $@


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


OBJ= hwm93.o nrlmsise00.o grid-1.00.o sami2-1.00.o

# Uncomment the compiler you want to use

#   1. Absoft

#  f77 = f77 -s -O3 -g

#    2. Lahey

# f77 = lf95  --sav -O

#    3. Portland Group

# f77 = pgf77 -Msave -fast

#    4. intel

  f77 = ifort -save -O2 -vec_report0

.f.o:
	$(f77) -c $*.f
#
sami2.x:    $(OBJ)
	$(f77)  -o sami2.x   $(OBJ)
#
clean:
	rm *.x *.o
#
$(OBJ): com-1.00.inc
$(OBJ): param-1.00.inc


