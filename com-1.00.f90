module commons
use parameters
implicit none
!*******************************************
!*******************************************

!            COM-1.00.INC

!*******************************************
!*******************************************

!     namelist data
      logical,save::fejer
      logical,save::fmtout
      real,dimension(nneut)::snn
      INTEGER,save::maxstep
      REAL,save::hrmax
      REAL,save::dt0
      REAL,save::dthr
      REAL,save::hrpr
      REAL,save::grad_in
      REAL,save::glat_in
      REAL,save::glon_in
      REAL,save::rmin
      REAL,save::rmax
      REAL,save::altmin
      REAL,save::fbar
      REAL,save::f10p7
      REAL,save::ap
      INTEGER,save::year
      INTEGER,save::day
      INTEGER,save::mmass
      REAL,save::nion1
      REAL,save::nion2
      REAL,save::hrinit
      REAL,save::tvn0
      REAL,save::tvexb0
      REAL,save::ve01
      REAL,save::gams
      REAL,save::gamp
      REAL,save::stn
      REAL,save::denmin
      REAL,save::alt_crit
      REAL,save::cqe


!     s grid data
      real,save,dimension(nz,nf)::alts

      real,save,dimension(nz,nf)::grs
      real,save,dimension(nz,nf)::glats
      real,save,dimension(nz,nf)::glons
      real,save,dimension(nz,nf)::bms
      real,save,dimension(nz)::gs
      real,save,dimension(nz,nf)::ps
      real,save,dimension(nz,nf)::blats
      real,save,dimension(nz,nf)::coschicrit
      real,save,dimension(nz,nf)::ds
      real,save,dimension(nz,nf)::d2s
      real,save,dimension(nz,nf)::d22s
      real,save,dimension(nz,nf)::dels
      real,save,dimension(nf)::grad_inp

      real,save,dimension(nzp1,nf)::xnorms
      real,save,dimension(nzp1,nf)::ynorms
      real,save,dimension(nzp1,nf)::znorms
      real,save,dimension(nz,nfp1)::xnormp
      real,save,dimension(nz,nfp1)::ynormp
      real,save,dimension(nz,nfp1)::znormp
      real,save,dimension(nz,nf)::arg
      real,save,dimension(nz,nf)::athg
      real,save,dimension(nz,nf)::aphig


!       alts     altitude  (in km) on s mesh
!       grs      radial geographic distance to field line on s mesh
!       glats    geographic latitude on s mesh
!       glons    geographic longitude on s mesh
!       bms      normalized magnetic field on s mesh (b/b0)
!       ds,d2,
!       d22s     differential `distances' used in diff eqs
!       dels     actual arc length of grid in s direction
 
!     p grid data
      real,save:: dt
      real,save,dimension(nz,nfp1)::delsp
      real,save,dimension(nz,nf)::vol
      real,save,dimension(nz,nfp1)::areap
      real,save,dimension(nz,nfp1)::areas
      real,save,dimension(nzp1,nfp1)::vnx
      real,save,dimension(nzp1,nfp1)::vny
      real,save,dimension(nzp1,nfp1)::vnz
      real,save,dimension(nz,nfp1,2)::xdels
      real,save,dimension(nzp1,nf,2)::xdelp
      real,save,dimension(nzp1,nf)::vexbs
      real,save,dimension(nz,nfp1)::vexbp
      real,save,dimension(nzp1,nfp1)::vexb

!     delsp      actual arc length of grid in s direction on p grid
!     vol        volume (i.e., area) of cell

!     chemical reaction data

      integer,dimension(nchem,3):: ichem
      real,dimension(nion,nneut,nchem)::ireact


!     variables

      real,save,dimension(nz,nf,nion)::deni
      real,save,dimension(nz,nf,nneut)::denn
      real,save,dimension(nz,nf)::ne
      real,save,dimension(nz,nf,nion)::vsi
      real,save,dimension(nz,nf,nion)::vsid
      real,save,dimension(nz,nf,nion)::sumvsi
      real,save,dimension(nz,nf,nion)::vsic
      real,save,dimension(nz,nf)::te
      real,save,dimension(nz,nf,nion)::ti
      real,save,dimension(nz,nf)::tn
      real,save,dimension(nz,nf)::u
      real,save,dimension(nz,nf)::v
      real,save,dimension(nz,nf)::vpi


!     velocity in radial (vor) and theta (vot) directions 

      real,save,dimension(nz,nf,nion)::vot
      real,save,dimension(nz,nf,nion)::vor


!     atomic masses

      real,save,dimension(nion)::ami
      real,save,dimension(nneut)::amn
      real,save,dimension(nneut)::alpha0
      real,save,dimension(7)::aap

!     zenith datt

      real,save,dimension(nz,nf)::cx

!     photodeposition rates
!     used 3 (absorption) and 7 (nion) explicitly
!     used 4 (number of angles in nighttime deposition)

      real,save,dimension(linesuv,3):: sigabsdt
      real,save,dimension(linesuv)::flux
      real,save,dimension(linesuv,7)::sigidt
      real,save,dimension(linesnt,7)::sigint
      real,save,dimension(nz,nf,91,linesnt)::fluxnt
      real,save,dimension(linesnt,4)::thetant
      real,save,dimension(linesnt,2)::zaltnt


!     diagnostic variables

      real,save,dimension(nz,nf,nion)::t1
      real,save,dimension(nz,nf,nion)::t2
      real,save,dimension(nz,nf,nion)::t3
      real,save,dimension(nz,nf)::u1
      real,save,dimension(nz,nf)::u2
      real,save,dimension(nz,nf)::u3
      real,save,dimension(nz,nf)::u4
      real,save,dimension(nz,nf)::u5

      real,save::x0
      real,save::y0
      real,save::z0
      real,save::plat
      real,save::plon
      real,save::bb0
 


      








end module commons



module inputfiles
      integer,parameter:: inputf=1

      character(6),parameter::inputpath="input"
      character(7),parameter::outputpath="output"
      character::delimeter='\/'
      integer,parameter::sami2_1_00_namelist=10
      integer,parameter::deni_init_inp=20
      integer,parameter::ichem_inp=30
      integer,parameter::phabsdt_inp=50
      integer,parameter::phiondt_inp=60
      integer,parameter::phionnt_inp=61
      integer,parameter::euvflux_inp=65
      integer,parameter::thetant_inp=66
      integer,parameter::zaltnt_inp=67

      !grid output files
      integer,parameter::zaltf_dat=69
      integer,parameter::glatf_dat=76
      integer,parameter::glonf_dat=77

      integer,parameter::time_dat=70
      integer,parameter::denif_dat=71
      integer,parameter::tif_dat=72
      integer,parameter::vsif_dat=73
      integer,parameter::tef_dat=75
      integer,parameter::vnf_dat=78

      integer,parameter::vtf_dat=90
      integer,parameter::vrf_dat=91
      integer,parameter::dennf_dat=92
      integer,parameter::vexbf_dat=93

      integer,parameter::t1f_dat=81
      integer,parameter::t2f_dat=82
      integer,parameter::t3f_dat=83
      integer,parameter::u1f_dat=84
      integer,parameter::u2f_dat=85
      integer,parameter::u3f_dat=86
      integer,parameter::u4f_dat=87
      integer,parameter::u5f_dat=88
      

contains


      subroutine open_file(fileunit,filename)
            implicit none
      INTEGER::ios=0
      CHARACTER(20)::filename
      INTEGER::fileunit
      
      open ( STATUS='OLD',IOSTAT=ios, unit=fileunit,  file=trim(inputpath)//delimeter//trim(filename))
      IF( ios.ne.0) THEN
          WRITE(6,*) 'Error opening file: '//trim(inputpath)//delimeter//trim(filename)
          STOP
      ENDIF
      end subroutine open_file

      subroutine open_input_files
      implicit none
      INTEGER::ios=0
      CHARACTER(20)::filename
!      CHARACTER(11)::fileform
      INTEGER::fileunit=0

      filename='sami2-1.00.namelist'
      fileunit=sami2_1_00_namelist
      call open_file(fileunit,filename)


      filename='deni-init.inp'
      fileunit=deni_init_inp
      call open_file(fileunit,filename)

      filename='ichem.inp'
      fileunit=ichem_inp
      call open_file(fileunit,filename)

      filename='phabsdt.inp'
      fileunit=phabsdt_inp
      call open_file(fileunit,filename)

      filename='phiondt.inp'
      fileunit=phiondt_inp
      call open_file(fileunit,filename)

      filename='phionnt.inp'
      fileunit=phionnt_inp
      call open_file(fileunit,filename)

      filename='euvflux.inp'
      fileunit=euvflux_inp
      call open_file(fileunit,filename)

      filename='thetant.inp'
      fileunit=thetant_inp
      call open_file(fileunit,filename)

      filename='zaltnt.inp'
      fileunit=zaltnt_inp
      call open_file(fileunit,filename)

      end subroutine open_input_files

      subroutine close_input_files
      implicit none
      close ( unit=sami2_1_00_namelist)
      close ( unit=deni_init_inp)
      close ( unit=ichem_inp)
      close ( unit=phabsdt_inp)
      close ( unit=phiondt_inp)
      close ( unit=phionnt_inp)
      close ( unit=euvflux_inp)
      close ( unit=thetant_inp)
      close ( unit=zaltnt_inp)


      end subroutine close_input_files

      subroutine open_output_grid_files
      use commons
      implicit none
      INTEGER::ios=0
      CHARACTER(20)::filename
      CHARACTER(11)::fileform
      INTEGER::fileunit=0

      call system('mkdir ' // adjustl(trim( outputpath ) ) )
!      call system('mkdir -p ' // adjustl(trim( outputpath ) ) )
      ios=0
      fileform=merge('formatted  ','unformatted',fmtout)

      filename='zalt'//merge('f','u',fmtout)//'.dat'
      fileunit=zaltf_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='glat'//merge('f','u',fmtout)//'.dat'
      fileunit=glatf_dat
      call open_file_with_replace(fileunit,filename,fileform)
 
      filename='glon'//merge('f','u',fmtout)//'.dat'
      fileunit=glonf_dat
      call open_file_with_replace(fileunit,filename,fileform)

      end subroutine open_output_grid_files

      subroutine close_output_grid_files
      implicit none

      close ( unit=zaltf_dat)
      close ( unit=glatf_dat)
      close ( unit=glonf_dat)


      end subroutine close_output_grid_files

      subroutine write_output_grid_files
      use commons
      implicit none
      if ( fmtout ) then
        write(zaltf_dat,100) alts
        write(glatf_dat,100) glats
        write(glonf_dat,100) glons
      else
        write(zaltf_dat) alts
        write(glatf_dat) glats
        write(glonf_dat) glons
      endif
 100  format (1x,1p10e16.6)

      end subroutine write_output_grid_files
      



      subroutine open_file_with_replace(fileunit,filename,fileform)
      implicit none
      INTEGER::ios=0
      CHARACTER(20)::filename
      CHARACTER(11)::fileform
      INTEGER::fileunit

      open ( STATUS='REPLACE',IOSTAT=ios, unit=fileunit,  file=trim(outputpath)//delimeter//trim(filename))
      IF( ios.ne.0) THEN
          WRITE(6,*) 'Error replacing file: '//trim(outputpath)//delimeter//trim(filename)
          STOP
      ENDIF
      close(unit=fileunit)
      ios=0
      open (ACCESS='APPEND',IOSTAT=ios,form=trim(fileform), unit=fileunit,  file=trim(outputpath)//delimeter//trim(filename))
      IF( ios.ne.0) THEN
          WRITE(6,*) 'Error opening file: '//trim(outputpath)//delimeter//trim(filename),ios
          STOP
      ENDIF

      end subroutine open_file_with_replace

      subroutine open_uf
      use commons
      implicit none
      INTEGER::ios=0
      CHARACTER(20)::filename
      CHARACTER(11)::fileform
      INTEGER::fileunit=0


      filename='time'//'.dat'
      fileunit=time_dat
      fileform=merge('formatted  ','unformatted',.true.)
      call open_file_with_replace(fileunit,filename,fileform)

      fileform=merge('formatted  ','unformatted',fmtout)

      filename='deni'//merge('f','u',fmtout)//'.dat'
      fileunit=denif_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='ti'//merge('f','u',fmtout)//'.dat'
      fileunit=tif_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='vsi'//merge('f','u',fmtout)//'.dat'
      fileunit=vsif_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='te'//merge('f','u',fmtout)//'.dat'
      fileunit=tef_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='vn'//merge('f','u',fmtout)//'.dat'
      fileunit=vnf_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='vt'//merge('f','u',fmtout)//'.dat'
      fileunit=vtf_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='vr'//merge('f','u',fmtout)//'.dat'
      fileunit=vrf_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='denn'//merge('f','u',fmtout)//'.dat'
      fileunit=dennf_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='vexb'//merge('f','u',fmtout)//'.dat'
      fileunit=vexbf_dat
      call open_file_with_replace(fileunit,filename,fileform)

!     diagnostic files (formatted)
      filename='t1'//merge('f','u',fmtout)//'.dat'
      fileunit=t1f_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='t2'//merge('f','u',fmtout)//'.dat'
      fileunit=t2f_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='t3'//merge('f','u',fmtout)//'.dat'
      fileunit=t3f_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='u1'//merge('f','u',fmtout)//'.dat'
      fileunit=u1f_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='u2'//merge('f','u',fmtout)//'.dat'
      fileunit=u2f_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='u3'//merge('f','u',fmtout)//'.dat'
      fileunit=u3f_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='u4'//merge('f','u',fmtout)//'.dat'
      fileunit=u4f_dat
      call open_file_with_replace(fileunit,filename,fileform)

      filename='u5'//merge('f','u',fmtout)//'.dat'
      fileunit=u5f_dat
      call open_file_with_replace(fileunit,filename,fileform)

      end subroutine open_uf

      subroutine open_u
      call open_uf
      return
      end

      subroutine open_f
      call open_uf
      return
      end

      subroutine close_uf
      close ( unit=zaltf_dat)
      close ( unit=time_dat)
      close ( unit=denif_dat)
      close ( unit=tif_dat)
      close ( unit=vsif_dat)
      close ( unit=tef_dat)
      close ( unit=vnf_dat)
      close ( unit=vtf_dat)
      close ( unit=vrf_dat)
      close ( unit=dennf_dat)
      close ( unit=vexbf_dat)
      close ( unit=t1f_dat)
      close ( unit=t2f_dat)
      close ( unit=t3f_dat)
      close ( unit=u1f_dat)
      close ( unit=u2f_dat)
      close ( unit=u3f_dat)
      close ( unit=u4f_dat)
      close ( unit=u5f_dat)
      
      end subroutine close_uf




end module inputfiles