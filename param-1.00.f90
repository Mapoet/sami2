!*******************************************
!*******************************************

!            PARAM-1.00.f90

!*******************************************
!*******************************************
module parameters
implicit none
!      number of altitudes
      INTEGER::nf=52
      INTEGER::nfp1 
      INTEGER::nfm1 
      INTEGER::nfm2 

!      number of grid cells
      INTEGER,PARAMETER::nz = 101
      INTEGER,PARAMETER::nzp1 = nz + 1
      INTEGER,PARAMETER::nzm1 = nz - 1

!      ion densities

       INTEGER,PARAMETER::nion  = 7    ! number of ions
       INTEGER,PARAMETER::pthp  = 1    ! h+
       INTEGER,PARAMETER::pthep = 5    ! he+
       INTEGER,PARAMETER::ptnp  = 7    ! n+
       INTEGER,PARAMETER::ptop  = 2    ! o+
       INTEGER,PARAMETER::ptn2p = 6    ! n2+
       INTEGER,PARAMETER::ptnop = 3    ! no+
       INTEGER,PARAMETER::pto2p = 4    ! o2+

!      neutrals 


       INTEGER,PARAMETER::nneut = 7    ! number of neutrals
       INTEGER,PARAMETER::pth   = 1    ! h
       INTEGER,PARAMETER::pthe  = 5    ! he
       INTEGER,PARAMETER::ptn   = 7    ! n
       INTEGER,PARAMETER::pto   = 2    ! o
       INTEGER,PARAMETER::ptn2  = 6    ! n2
       INTEGER,PARAMETER::ptno  = 3    ! no
       INTEGER,PARAMETER::pto2  = 4    ! o2

!      number of chemical reactions

       INTEGER,PARAMETER::nchem = 21

!      various constants

!      ftnchek is giving some meaningless errors about precision, 
!      but i am going to lower the precision of some of these
!      variables to keep down the error messages


       REAL,PARAMETER::pie    = 3.1415927   
       REAL,PARAMETER::po180  = 1.745329e-02
       REAL,PARAMETER::rtod   = 57.295780   
       REAL,PARAMETER::tm18   = 1.e-18      

       REAL,PARAMETER::spday  = 86400.
       REAL,PARAMETER::sphr = 3600.

       REAL,PARAMETER::gzero  = 980.665
       REAL,PARAMETER::re = 6370.
       REAL,PARAMETER::bmag = 0.25

       REAL,PARAMETER::bolt   = 1.38044e-16 
       REAL,PARAMETER::amu    = 1.67252e-24 
       REAL,PARAMETER::evtok  = 1.1604e+04   

       INTEGER,PARAMETER::linesuv = 37
       INTEGER,PARAMETER::linesnt = 4

       REAL,PARAMETER::dayve = 80.
       REAL,PARAMETER::sidyr = 365.4
       REAL,PARAMETER::solinc = 23.5

!      these are for the error function

       REAL,PARAMETER::pas =   .3275911
       REAL,PARAMETER::z1 =  .2548295
       REAL,PARAMETER::z2  = - .284496
       REAL,PARAMETER::z3 = 1.421413
       REAL,PARAMETER::z4  = -1.453152
       REAL,PARAMETER::z5 = 1.0614054

contains
      subroutine init_param
      implicit none
      nfp1 = nf + 1
      nfm1 = nf - 1
      nfm2 = nf - 2
      end subroutine init_param
end module parameters