module commons
use parameters
implicit none
!*******************************************
!*******************************************

!            COM-1.00.INC

!*******************************************
!*******************************************

!     namelist data
      logical::fejer
      logical::fmtout
      real,dimension(nneut)::snn


!     s grid data
      real,dimension(nz,nf)::alts

      real,dimension(nz,nf)::grs
      real,dimension(nz,nf)::glats
      real,dimension(nz,nf)::glons
      real,dimension(nz,nf)::bms
      real,dimension(nz)::gs
      real,dimension(nz,nf)::ps
      real,dimension(nz,nf)::blats
      real,dimension(nz,nf)::coschicrit
      real,dimension(nz,nf)::ds
      real,dimension(nz,nf)::d2s
      real,dimension(nz,nf)::d22s
      real,dimension(nz,nf)::dels
      real,dimension(nf)::grad_inp

      real,dimension(nzp1,nf)::xnorms
      real,dimension(nzp1,nf)::ynorms
      real,dimension(nzp1,nf)::znorms
      real,dimension(nz,nfp1)::xnormp
      real,dimension(nz,nfp1)::ynormp
      real,dimension(nz,nfp1)::znormp
      real,dimension(nz,nf)::arg
      real,dimension(nz,nf)::athg
      real,dimension(nz,nf)::aphig


!       alts     altitude  (in km) on s mesh
!       grs      radial geographic distance to field line on s mesh
!       glats    geographic latitude on s mesh
!       glons    geographic longitude on s mesh
!       bms      normalized magnetic field on s mesh (b/b0)
!       ds,d2,
!       d22s     differential `distances' used in diff eqs
!       dels     actual arc length of grid in s direction
 
!     p grid data

      real,dimension(nz,nfp1)::delsp
      real,dimension(nz,nf)::vol
      real,dimension(nz,nfp1)::areap
      real,dimension(nz,nfp1)::areas
      real,dimension(nzp1,nfp1)::vnx
      real,dimension(nzp1,nfp1)::vny
      real,dimension(nzp1,nfp1)::vnz
      real,dimension(nz,nfp1,2)::xdels
      real,dimension(nzp1,nf,2)::xdelp
      real,dimension(nzp1,nf)::vexbs
      real,dimension(nz,nfp1)::vexbp
      real,dimension(nzp1,nfp1)::vexb

!     delsp      actual arc length of grid in s direction on p grid
!     vol        volume (i.e., area) of cell

!     chemical reaction data

      integer,dimension(nchem,3):: ichem
      real,dimension(nion,nneut,nchem)::ireact


!     variables

      real,dimension(nz,nf,nion)::deni
      real,dimension(nz,nf,nneut)::denn
      real,dimension(nz,nf)::ne
      real,dimension(nz,nf,nion)::vsi
      real,dimension(nz,nf,nion)::vsid
      real,dimension(nz,nf,nion)::sumvsi
      real,dimension(nz,nf,nion)::vsic
      real,dimension(nz,nf)::te
      real,dimension(nz,nf,nion)::ti
      real,dimension(nz,nf)::tn
      real,dimension(nz,nf)::u
      real,dimension(nz,nf)::v
      real,dimension(nz,nf)::vpi


!     velocity in radial (vor) and theta (vot) directions 

      real,dimension(nz,nf,nion)::vot
      real,dimension(nz,nf,nion)::vor


!     atomic masses

      real,dimension(nion)::ami
      real,dimension(nneut)::amn
      real,dimension(nneut)::alpha0
      real,dimension(7)::aap

!     zenith datt

      real,dimension(nz,nf)::cx

!     photodeposition rates
!     used 3 (absorption) and 7 (nion) explicitly
!     used 4 (number of angles in nighttime deposition)

      real,dimension(linesuv,3):: sigabsdt
      real,dimension(linesuv)::flux
      real,dimension(linesuv,7)::sigidt
      real,dimension(linesnt,7)::sigint
      real,dimension(nz,nf,91,linesnt)::fluxnt
      real,dimension(linesnt,4)::thetant
      real,dimension(linesnt,2)::zaltnt


!     diagnostic variables

      real,dimension(linesnt,2)::t1(nz,nf,nion)
      real,dimension(linesnt,2)::t2(nz,nf,nion)
      real,dimension(linesnt,2)::t3(nz,nf,nion)
      real,dimension(linesnt,2)::u1(nz,nf)
      real,dimension(linesnt,2)::u2(nz,nf)
      real,dimension(linesnt,2)::u3(nz,nf)
      real,dimension(linesnt,2)::u4(nz,nf)
      real,dimension(linesnt,2)::u5(nz,nf)

      real::x0
      real::y0
      real::z0
      real::plat
      real::plon
      real::bb0
 
      



end module commons



module inputfiles
      integer,parameter:: inputf=1
      char(len=20)::inputpath='input\/'

      integer,parameter::sami2_1_00_namelist=10
      integer,parameter::deni_init_inp=20
      integer,parameter::ichem_inp=30
      integer,parameter::phabsdt_inp=50
      integer,parameter::phiondt_inp=60
      integer,parameter::phionnt_inp=61
      integer,parameter::euvflux_inp=65
      integer,parameter::thetant_inp=66
      integer,parameter::zaltnt_inp=67
contains
      subroutine open_input_files
      open ( unit=sami2_1_00_namelist, file=trim(inputpath)//'sami2-1.00.namelist')
      open ( unit=deni_init_inp, file='input\/deni-init.inp')
      open ( unit=ichem_inp, file='input\/ichem.inp')
      open ( unit=phabsdt_inp, file='input\/phabsdt.inp')
      open ( unit=phiondt_inp, file='input\/phiondt.inp')
      open ( unit=phionnt_inp, file='input\/phionnt.inp')
      open ( unit=euvflux_inp, file='input\/euvflux.inp')
      open ( unit=thetant_inp, file='input\/thetant.inp')
      open ( unit=zaltnt_inp, file='input\/zaltnt.inp')

      end subroutine open_input_files

      subroutine close_input_files
      close ( sami2_1_00_namelist)
      close ( deni_init_inp)
      close ( unit=ichem_inp)
      close ( unit=phabsdt_inp)
      close ( unit=phiondt_inp)
      close ( unit=phionnt_inp)
      close ( unit=euvflux_inp)
      close ( unit=thetant_inp)
      close ( unit=zaltnt_inp)


      end subroutine close_input_files
end module inputfiles