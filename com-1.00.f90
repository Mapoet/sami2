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
      real,save,dimension(:,:),allocatable::alts

      real,save,dimension(:,:),allocatable::grs
      real,save,dimension(:,:),allocatable::glats
      real,save,dimension(:,:),allocatable::glons
      real,save,dimension(:,:),allocatable::bms
      real,save,dimension(:),allocatable::gs
      real,save,dimension(:,:),allocatable::ps
      real,save,dimension(:,:),allocatable::blats
      real,save,dimension(:,:),allocatable::coschicrit
      real,save,dimension(:,:),allocatable::ds
      real,save,dimension(:,:),allocatable::d2s
      real,save,dimension(:,:),allocatable::d22s
      real,save,dimension(:,:),allocatable::dels
      real,save,dimension(:),allocatable::grad_inp

      real,save,dimension(:,:),allocatable::xnorms
      real,save,dimension(:,:),allocatable::ynorms
      real,save,dimension(:,:),allocatable::znorms
      real,save,dimension(:,:),allocatable::xnormp
      real,save,dimension(:,:),allocatable::ynormp
      real,save,dimension(:,:),allocatable::znormp
      real,save,dimension(:,:),allocatable::arg
      real,save,dimension(:,:),allocatable::athg
      real,save,dimension(:,:),allocatable::aphig


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
      real,save,dimension(:,:),allocatable::delsp
      real,save,dimension(:,:),allocatable::vol
      real,save,dimension(:,:),allocatable::areap
      real,save,dimension(:,:),allocatable::areas
      real,save,dimension(:,:),allocatable::vnx
      real,save,dimension(:,:),allocatable::vny
      real,save,dimension(:,:),allocatable::vnz
      real,save,dimension(:,:,:),allocatable::xdels
      real,save,dimension(:,:,:),allocatable::xdelp
      real,save,dimension(:,:),allocatable::vexbs
      real,save,dimension(:,:),allocatable::vexbp
      real,save,dimension(:,:),allocatable::vexb

!     delsp      actual arc length of grid in s direction on p grid
!     vol        volume (i.e., area) of cell

!     chemical reaction data

      integer,dimension(nchem,3):: ichem
      real,dimension(nion,nneut,nchem)::ireact


!     variables

      real,save,dimension(:,:,:),allocatable::deni
      real,save,dimension(:,:,:),allocatable::denn
      real,save,dimension(:,:),allocatable::ne
      real,save,dimension(:,:,:),allocatable::vsi
      real,save,dimension(:,:,:),allocatable::vsid
      real,save,dimension(:,:,:),allocatable::sumvsi
      real,save,dimension(:,:,:),allocatable::vsic
      real,save,dimension(:,:),allocatable::te
      real,save,dimension(:,:,:),allocatable::ti
      real,save,dimension(:,:),allocatable::tn
      real,save,dimension(:,:),allocatable::u
      real,save,dimension(:,:),allocatable::v
      real,save,dimension(:,:),allocatable::vpi


!     velocity in radial (vor) and theta (vot) directions 

      real,save,dimension(:,:,:),allocatable::vot
      real,save,dimension(:,:,:),allocatable::vor


!     atomic masses

      real,save,dimension(nion)::ami
      real,save,dimension(nneut)::amn
      real,save,dimension(nneut)::alpha0
      real,save,dimension(7)::aap

!     zenith datt

      real,save,dimension(:,:),allocatable::cx

!     photodeposition rates
!     used 3 (absorption) and 7 (nion) explicitly
!     used 4 (number of angles in nighttime deposition)

      real,save,dimension(linesuv,3):: sigabsdt
      real,save,dimension(linesuv)::flux
      real,save,dimension(linesuv,7)::sigidt
      real,save,dimension(linesnt,7)::sigint
      real,save,dimension(:,:,:,:),allocatable::fluxnt
      real,save,dimension(linesnt,4)::thetant
      real,save,dimension(linesnt,2)::zaltnt


!     diagnostic variables

      real,save,dimension(:,:,:),allocatable::t1
      real,save,dimension(:,:,:),allocatable::t2
      real,save,dimension(:,:,:),allocatable::t3
      real,save,dimension(:,:),allocatable::u1
      real,save,dimension(:,:),allocatable::u2
      real,save,dimension(:,:),allocatable::u3
      real,save,dimension(:,:),allocatable::u4
      real,save,dimension(:,:),allocatable::u5

      real,save::x0
      real,save::y0
      real,save::z0
      real,save::plat
      real,save::plon
      real,save::bb0
 

contains
      subroutine init_memory
      include 'param-1.00.inc'
      implicit none

      ALLOCATE(alts(nz,nf))
      ALLOCATE(grs(nz,nf))
      ALLOCATE(glats(nz,nf))
      ALLOCATE(glons(nz,nf))
      ALLOCATE(bms(nz,nf))
      ALLOCATE(gs(nz))
      ALLOCATE(ps(nz,nf))
      ALLOCATE(blats(nz,nf))
      ALLOCATE(coschicrit(nz,nf))
      ALLOCATE(ds(nz,nf))
      ALLOCATE(d2s(nz,nf))
      ALLOCATE(d22s(nz,nf))
      ALLOCATE(dels(nz,nf))
      ALLOCATE(grad_inp(nf))

      ALLOCATE(xnorms(nzp1,nf))
      ALLOCATE(ynorms(nzp1,nf))
      ALLOCATE(znorms(nzp1,nf))
      ALLOCATE(xnormp(nz,nfp1))
      ALLOCATE(ynormp(nz,nfp1))
      ALLOCATE(znormp(nz,nfp1))
      ALLOCATE(arg(nz,nf))
      ALLOCATE(athg(nz,nf))
      ALLOCATE(aphig(nz,nf))

      ALLOCATE(delsp(nz,nfp1))
      ALLOCATE(vol(nz,nf))
      ALLOCATE(areap(nz,nfp1))
      ALLOCATE(areas(nz,nfp1))
      ALLOCATE(vnx(nzp1,nfp1))
      ALLOCATE(vny(nzp1,nfp1))
      ALLOCATE(vnz(nzp1,nfp1))
      ALLOCATE(xdels(nz,nfp1,2))
      ALLOCATE(xdelp(nzp1,nf,2))
      ALLOCATE(vexbs(nzp1,nf))
      ALLOCATE(vexbp(nz,nfp1))
      ALLOCATE(vexb(nzp1,nfp1))

      ALLOCATE(deni(nz,nf,nion))
      ALLOCATE(denn(nz,nf,nneut))
      ALLOCATE(ne(nz,nf))
      ALLOCATE(vsi(nz,nf,nion))
      ALLOCATE(vsid(nz,nf,nion))
      ALLOCATE(sumvsi(nz,nf,nion))
      ALLOCATE(vsic(nz,nf,nion))
      ALLOCATE(te(nz,nf))
      ALLOCATE(ti(nz,nf,nion))
      ALLOCATE(tn(nz,nf))
      ALLOCATE(u(nz,nf))
      ALLOCATE(v(nz,nf))
      ALLOCATE(vpi(nz,nf))


      ALLOCATE(vot(nz,nf,nion))
      ALLOCATE(vor(nz,nf,nion))



      ALLOCATE(cx(nz,nf))

      ALLOCATE(fluxnt(nz,nf,91,linesnt))

      ALLOCATE(t1(nz,nf,nion))
      ALLOCATE(t2(nz,nf,nion))
      ALLOCATE(t3(nz,nf,nion))
      ALLOCATE(u1(nz,nf))
      ALLOCATE(u2(nz,nf))
      ALLOCATE(u3(nz,nf))
      ALLOCATE(u4(nz,nf))
      ALLOCATE(u5(nz,nf))

      end subroutine init_memory

      subroutine deinit_memory
      implicit none
      DEALLOCATE(alts(nz,nf))
      DEALLOCATE(grs(nz,nf))
      DEALLOCATE(glats(nz,nf))
      DEALLOCATE(glons(nz,nf))
      DEALLOCATE(bms(nz,nf))
      DEALLOCATE(gs(nz))
      DEALLOCATE(ps(nz,nf))
      DEALLOCATE(blats(nz,nf))
      DEALLOCATE(coschicrit(nz,nf))
      DEALLOCATE(ds(nz,nf))
      DEALLOCATE(d2s(nz,nf))
      DEALLOCATE(d22s(nz,nf))
      DEALLOCATE(dels(nz,nf))
      DEALLOCATE(grad_inp(nf))

      DEALLOCATE(xnorms(nzp1,nf))
      DEALLOCATE(ynorms(nzp1,nf))
      DEALLOCATE(znorms(nzp1,nf))
      DEALLOCATE(xnormp(nz,nfp1))
      DEALLOCATE(ynormp(nz,nfp1))
      DEALLOCATE(znormp(nz,nfp1))
      DEALLOCATE(arg(nz,nf))
      DEALLOCATE(athg(nz,nf))
      DEALLOCATE(aphig(nz,nf))

      DEALLOCATE(delsp(nz,nfp1))
      DEALLOCATE(vol(nz,nf))
      DEALLOCATE(areap(nz,nfp1))
      DEALLOCATE(areas(nz,nfp1))
      DEALLOCATE(vnx(nzp1,nfp1))
      DEALLOCATE(vny(nzp1,nfp1))
      DEALLOCATE(vnz(nzp1,nfp1))
      DEALLOCATE(xdels(nz,nfp1,2))
      DEALLOCATE(xdelp(nzp1,nf,2))
      DEALLOCATE(vexbs(nzp1,nf))
      DEALLOCATE(vexbp(nz,nfp1))
      DEALLOCATE(vexb(nzp1,nfp1))

      DEALLOCATE(deni(nz,nf,nion))
      DEALLOCATE(denn(nz,nf,nneut))
      DEALLOCATE(ne(nz,nf))
      DEALLOCATE(vsi(nz,nf,nion))
      DEALLOCATE(vsid(nz,nf,nion))
      DEALLOCATE(sumvsi(nz,nf,nion))
      DEALLOCATE(vsic(nz,nf,nion))
      DEALLOCATE(te(nz,nf))
      DEALLOCATE(ti(nz,nf,nion))
      DEALLOCATE(tn(nz,nf))
      DEALLOCATE(u(nz,nf))
      DEALLOCATE(v(nz,nf))
      DEALLOCATE(vpi(nz,nf))


      DEALLOCATE(vot(nz,nf,nion))
      DEALLOCATE(vor(nz,nf,nion))



      DEALLOCATE(cx(nz,nf))

      DEALLOCATE(fluxnt(nz,nf,91,linesnt))

      DEALLOCATE(t1(nz,nf,nion))
      DEALLOCATE(t2(nz,nf,nion))
      DEALLOCATE(t3(nz,nf,nion))
      DEALLOCATE(u1(nz,nf))
      DEALLOCATE(u2(nz,nf))
      DEALLOCATE(u3(nz,nf))
      DEALLOCATE(u4(nz,nf))
      DEALLOCATE(u5(nz,nf))
      
      end subroutine deinit_memory






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

      call system('mkdir ' // adjustl(trim( outputpath ) ) )
      
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



!*******************************************
!*******************************************

!             output

!*******************************************
!*******************************************

       subroutine output ( hrut,ntm,istep )

       include 'param-1.00.inc'
      use commons
      implicit none
      REAL::hrut
      INTEGER::ntm
      INTEGER::istep

      !LOCAL VARIABLES
      REAL::hr24,totsec,thr,tmin,tsec
      INTEGER::nthr,ntmin,ntsec

       hr24   = mod (hrut,24.)
       totsec = hr24 * 3600.
       thr    = totsec / 3600.
       nthr   = int(thr)
       tmin   = ( thr - nthr ) * 60.
       ntmin  = int(mod(tmin,60.))
       tsec   = ( tmin - ntmin ) * 60.
       ntsec  = int(tsec)

!      put data into output variables

       write (*,*) 'istep = ',istep,'ntm = ',ntm,&
                   'time step = ',dt,' hour = ',hrut
       write (70,100) ntm,nthr,ntmin,ntsec

       if ( fmtout ) then
         write(71,101) deni
         write(72,101) ti
         write(73,101) vsi
         write(75,101) te
!         write(78,101) vn
!         write(81,101) t1
!         write(82,101) t2
!         write(83,101) t3
!         write(84,101) u1
!         write(85,101) u2
!         write(86,101) u3
!         write(87,101) u4
!         write(88,101) u5
!         write(90,101) vot
!         write(91,101) vor
!         write(92,101) denn
       endif 

       if ( .not. fmtout ) then
         write(71) deni
         write(72) ti
         write(73) vsi
         write(75) te
!         write(78) vn
!         write(81) t1
!         write(82) t2
!         write(83) t3
         write(84) u1
         write(85) u2
         write(86) u3
         write(87) u4
         write(88) u5
!         write(90) vot
!         write(91) vor
!         write(92) denn
!         write(93) vexbp
       endif

 100   format(1x,4i6)
 101   format(1x,1p10e16.6)

       return
       end




end module inputfiles