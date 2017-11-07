!     *******************************************
!     *******************************************
!
!                  SAMI2_MPI-1.00_HEATING
!
!     *******************************************
!     *******************************************


!                   LICENSE

!     *******************************************
!     *******************************************

!     I hereby agree to the following terms governing the use and 
!     redistribution of the SAMI2 software release written and 
!     developed by J.D. Huba, G. Joyce and M. Swisdak.

!     Redistribution and use in source and binary forms, with or 
!     without modification, are permitted provided that (1) source code 
!     distributions retain this paragraph in its entirety, (2) distributions 
!     including binary code include this paragraph in its entirety in 
!     the documentation or other materials provided with the distribution, 
!     (3) improvements, additions and upgrades to the software will be 
!     provided to NRL Authors in computer readable form, with an unlimited, 
!     royalty-free license to use these improvements, additions and upgrades, 
!     and the authority to grant unlimited royalty-free sublicenses to these 
!     improvements and (4) all published research 
!     using this software display the following acknowledgment ``This 
!     work uses the SAMI2 ionosphere model written and developed 
!     by the Naval Research Laboratory.'' 

!     Neither the name of NRL or its contributors, nor any entity of the 
!     United States Government may be used to endorse or promote products 
!     derived from this software, nor does the inclusion of the NRL written 
!     and developed software directly or indirectly suggest NRL's or the 
!     United States Government's endorsement of this product.


!     THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR IMPLIED 
!     WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF 
!     MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.


!     *******************************************
!     *******************************************
!
      program main
      

      include 'param-1.00.inc'
      include 'com-1.00.inc' 
      use commonsubroutines
#ifdef _USE_MPI_
      use mpi_client
#endif
      implicit none
#ifdef _USE_MPI_
      include "mpif.h"
      integer::status(MPI_STATUS_SIZE)
      !MPI_Request request
      INTEGER,DIMENSION(:),allocatable::nfsize,nfstindex
      INTEGER::itask
      REAL::dtnew
#endif

      integer::ntm,istep,nfl
      real::hrut,timemax,time,tprnt,tneut
      INTEGER::i

#ifdef _USE_MPI_
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
      if(taskid .EQ. 0) then
#endif
      call init_param
      call initial
  
!     open output files

      call open_uf

      ntm   = 0
      istep = 0
!      call output ( hrinit,ntm,istep )

      close (68)

#ifdef _USE_MPI_      
      endif
      if((taskid .EQ. 0).and.(numtasks.gt.1)) then
      ALLOCATE(nfsize(numtasks-1),nfstindex(numtasks-1))
      nfstindex=1
      i=1
      do i = 1,numtasks-1,1
          call MPI_SEND(merge(i-1,numtasks-1,i .ne. 1),1,MPI_INT,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(merge(i+1,1,i .ne. numtasks-1),1,MPI_INT,i,0,MPI_COMM_WORLD,status,ierr)

          nfsize(i)=merge((nf-2)/(numtasks-1),nf-((nf-2)/(numtasks-1))*(numtasks-2)-1,i.ne.numtasks-1)+merge(1,0,(i .eq.1).and.(i .ne.(numtasks-1)))
          if(i.gt.1) then 
          nfstindex(i)=merge(nfstindex(i-1)+nfsize(i-1),1,i .ne. 1)
          endif

          call MPI_SEND(nfsize(i)+merge(1,0,i.ne.1)+merge(1,0,i.ne.numtasks-1),1,MPI_INT,i,0,MPI_COMM_WORLD,ierr)
          call share_data_server(nfstindex(i)-merge(1,0,i.ne.1),nfstindex(i)+nfsize(i)-1+merge(1,0,i.ne.numtasks-1),i)
          
      enddo 
      endif

      if(taskid .NE. 0) then
          call MPI_RECV(left,1,MPI_INT,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(right,1,MPI_INT,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(nf,1,MPI_INT,0,0,MPI_COMM_WORLD,status,ierr)

          call init_param
          call init_memory
          call share_data_client!share chanks
          

      endif

      !if(taskid .EQ. 0) then
      !if((taskid .NE. 0).or.(numtasks.eq.1)) then
#endif

!     time loop
      hrut    = hrinit
      timemax = hrmax * sphr
      istep   = 0
      tprnt   = 0.
      tneut   = 0.
      time    = 0.
      call neutambt (hrinit) 
      do while (      istep .le. maxstep &
                .and. time  .lt. timemax  )

#ifdef _USE_MPI_      
      if((taskid .NE. 0).or.(numtasks.eq.1)) then
#endif


!       parallel transport
        do nfl = nf,1,-1
          call zenith (hrut,nfl)
          call transprt (nfl)
        enddo         

!       perpendicular transport
        call exb(hrut)         


#ifdef _USE_MPI_
        endif
#endif


#ifdef _USE_MPI_      
      if((taskid .NE. 0).or.(numtasks.eq.1)) then
#endif


#ifdef _USE_MPI_
        endif
#endif

!       time/step advancement

        istep  = istep + 1
        time   = time  + dt
        hrut   = time / sphr + hrinit
        tprnt  = tprnt + dt / sphr
        tneut  = tneut + dt / sphr

#ifdef _USE_MPI_      
      if((taskid .NE. 0).or.(numtasks.eq.1)) then
#endif

        call courant ( hrut )

!       update neutrals

#ifdef _USE_MPI_
        if(numtasks.gt.1) then
          call MPI_SEND(dt,1,MPI_REAL,0,MPI_TAG_DT_SYNC,MPI_COMM_WORLD,ierr)
        endif
        endif

        if((taskid.eq.0).and.(numtasks.gt.1)) then
          dt = dt0
          do itask=1,numtasks-1
            call MPI_Probe(MPI_ANY_SOURCE, MPI_TAG_DT_SYNC, MPI_COMM_WORLD, status, ierr)
            CALL MPI_RECV(dtnew,1,MPI_REAL,status(MPI_SOURCE),MPI_TAG_DT_SYNC,MPI_COMM_WORLD, status, ierr)
            dt=min(dt,dtnew)
          enddo

        endif
        call mpi_bcast(dt,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

#ifdef _USE_MPI_      
      if((taskid .NE. 0).or.(numtasks.eq.1)) then
#endif
        if ( tneut .ge. 0.25 ) then
          !call neutambt (hrut) 
          tneut = 0.
        endif


#ifdef _USE_MPI_
        endif
#endif


!       output data

        if ( tprnt .ge. dthr .and. hrut .gt. hrpr+hrinit) then
          ntm      = ntm + 1
#ifdef _USE_MPI_
        if((taskid .ne.0)) then  
          !send data to master
          call share_output_data_client()
          else
        if(numtasks .gt.1)then
          !Gather data from workers
          do itask=1,numtasks-1
            call share_output_data_server(nfsize,nfstindex)
          enddo
         call flush(6)
        endif   
#endif
          call output ( hrut,ntm,istep )
#ifdef _USE_MPI_
        endif
#endif
          tprnt   = 0.
          print *,'ouput - hour = ',hrut
          
        elseif ( tprnt .ge. dthr ) then
#ifdef _USE_MPI_
        if(taskid.EQ.0)then      
#endif
          print *,'no ouput yet - hour = ',hrut
#ifdef _USE_MPI_
        endif      
#endif
          tprnt   = 0.
        endif
         call flush(6)
#ifdef _USE_MPI_
        if(taskid.ne.0) then
              
              !send data to neighbrs 
        endif      

#endif

      enddo    ! end time loop

!     close files

      call close_uf
      
#ifdef _USE_MPI_      
!      endif


      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if((taskid.eq.0).and.(numtasks.gt.1))then
      DEALLOCATE(nfsize(numtasks-1),nfstindex(numtasks-1))
      endif
#endif

      call deinit_memory

#ifdef _USE_MPI_
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
#endif
!      stop
      end


!*******************************************
!*******************************************

!            initial

!*******************************************
!*******************************************

      subroutine initial

      include 'param-1.00.inc'
      include 'com-1.00.inc' 
      use commonsubroutines
      implicit none
      include "gonamelist.inc"
      !LOCAL VARIABLES
      real,dimension(nz,nf,91):: f1026
      real,dimension(nz,nf,91):: f584
      real,dimension(nz,nf,91):: f304 
      real,dimension(nz,nf,91):: f1216
      real,dimension(29):: zi
      real,dimension(29,7):: denii
      real,dimension(9):: d
      real,dimension(2):: temp
      real,dimension(2):: app
      real,dimension(2):: whm93
      real,dimension(linesuv,5):: phionr
      real,dimension(linesuv,2):: fluxdat
      REAL::slope,hrl,f74,ai,xflux,sec,p
      INTEGER::i,j,k,k0,j0,n,l,nn,iyd


!     open input files
      call open_input_files

!     read in parameters and initial ion density data 

      read(sami2_1_00_namelist,go)
      call init_memory





      dt = dt0

      ami(pthp)  = 1.
      ami(pthep) = 4.
      ami(ptnp)  = 14.
      ami(ptop)  = 16.
      ami(ptn2p) = 28.
      ami(ptnop) = 30.
      ami(pto2p) = 32.

      amn(pth)  = 1.
      amn(pthe) = 4.
      amn(ptn)  = 14.
      amn(pto)  = 16.
      amn(ptn2) = 28.
      amn(ptno) = 30.
      amn(pto2) = 32.

      alpha0(pth)  = 0.67
      alpha0(pthe) = 0.21
      alpha0(ptn)  = 1.10
      alpha0(pto)  = 0.79
      alpha0(ptn2) = 1.76
      alpha0(ptno) = 1.74
      alpha0(pto2) = 1.59

      aap = ap

!     read in initial density data

      do i = 1,29
        read(deni_init_inp,102) zi(i),(denii(i,j),j=1,7)
 102    format(1x,f7.1,1p7e8.1)
      enddo


!     read in chemistry data
!     in format statement 104 need to 'hardwire' nneut (= 7)

      do k = 1,nchem
        read(ichem_inp,103) (ichem(k,j),j=1,3)
 103    format(3i3)
      enddo

!     build reaction matrix

      do k = 1,nchem
        do j = 1,nneut
          do i = nion1,nion2
            ireact(i,j,k) = 0
            if (      i .eq. ichem(k,1)&
               .and. j .eq. ichem(k,2) ) ireact(i,j,k) = 1.
          enddo
        enddo
      enddo

!     generate the mesh data

      call grid

!     output grid data
      call open_output_grid_files
      call write_output_grid_files
      call close_output_grid_files


! MS: chicrit is the zenith angle below which the Sun is visible.
! For points on the surface this is just pi/2, but at higher
! altitudes it is bigger.

        do j = 1,nf
          do i = 1,nz
            coschicrit(i,j) = cos(pie -&
                       asin( 1./ (1. + alts(i,j)/re) ))
          enddo
        enddo


!     put deni on mesh via linear interpolation
!     and put on lower limit

!     initialize all ion densities
      j0 = 1
      do n = 1,nion
        do l = 1,nf
          do i = 1,nz
            j = 1
            do while (  alts(i,l) .ge. zi(j) .and. j .le. 28 )
              j0 = j
              j  = j + 1
            enddo
            if ( n .eq. 1 ) k = pthp
            if ( n .eq. 2 ) k = pthep
            if ( n .eq. 3 ) k = ptnp
            if ( n .eq. 4 ) k = ptop
            if ( n .eq. 5 ) k = ptn2p
            if ( n .eq. 6 ) k = ptnop
            if ( n .eq. 7 ) k = pto2p
            slope   = ( denii(j0+1,n) - denii(j0,n) ) &
                      / ( zi   (j0+1)   - zi   (j0) )
            deni(i,l,k) = denii(j0,n) + ( alts(i,l) - zi(j0) ) * slope
            deni(i,l,k) = amax1 ( deni(i,l,k) , denmin )
            if ( alts(i,l) .gt. zi(29) ) then
              if ( n .eq. 1 )  then
                nn = pthp
                deni(i,l,k) = denii(29,n)
              else
                deni(i,l,k) = denmin
              endif
            endif
          enddo
        enddo
      enddo


!     initialize neutrals

!     neutral density and temperature
!      stop
      do j = 1,nf
        do i = 1,nz
          hrl = hrinit + glons(i,j) / 15.
          do while ( hrl .ge. 24. ) 
               hrl = hrl - 24.
         enddo
          
          call msistim ( int(year),int(day),hrl,glons(i,j),iyd,sec )
          call gtd7 ( iyd,sec,alts(i,j),glats(i,j),glons(i,j),&
                     hrl,fbar,f10p7,aap,mmass,d,temp )
          denn(i,j,pth )  = snn(pth)  * d(7)
          denn(i,j,pthe)  = snn(pthe) * d(1)
          denn(i,j,ptn )  = snn(ptn)  * d(8)
          denn(i,j,pto )  = snn(pto)  * d(2)
          denn(i,j,ptn2)  = snn(ptn2) * d(3) + 1.e-30
          denn(i,j,pto2)  = snn(pto2) * d(4) + 1.e-30
          tn(i,j)         = stn * temp(2)
          denn(i,j,ptno)  = 0.4 * exp( -3700. / tn(i,j) )&
                           * denn(i,j,pto2)&
                           + 5.0e-7 * denn(i,j,pto)
        enddo
      enddo
!     electron and ion temperature initialization

      do k = nion1,nion2
        do j = 1,nf
          do i = 1,nz
            te(i,j)   = tn(i,j)
            ti(i,j,k) = tn(i,j)
          enddo
        enddo
      enddo
!     initialize ion velocity to zero


      do k = nion1,nion2
        do j = 1,nf
          do i = 1,nz
            vsi(i,j,k)     = 0.
            sumvsi(i,j,k)  = 0.
          enddo
        enddo
      enddo

!     neutral winds (convert from m/s to cm/s)

      do j = 1,nf
        do i = 1,nz
          app(1)   = ap
          app(2)   = ap
          hrl = hrinit + glons(i,j) / 15.
          do while ( hrl .ge. 24. ) 
               hrl = hrl - 24.
         enddo
          call msistim ( int(year),int(day),hrl,glons(i,j),iyd,sec )
          call gws5 ( iyd,sec,alts(i,j),glats(i,j),glons(i,j),&
                     hrl,fbar,f10p7,app,whm93                )
          v(i,j)   = 100. * whm93(1)
          u(i,j)   = 100. * whm93(2)
        enddo
      enddo

!     read in photoabsorption rates

      do i = 1,linesuv
        read (phabsdt_inp,105) (sigabsdt(i,j), j=1,3)
 105    format (3f7.2) 
      enddo 

      do j = 1,3
        do i = 1,linesuv
          sigabsdt(i,j) = tm18 * sigabsdt(i,j) 
        enddo 
      enddo 
 
!     initialize photoionization rates to zero

      do j = 1,nneut
        do i = 1,linesuv
          sigidt(i,j)  = 0.
        enddo
        do i = 1,linesnt
          sigint(i,j)  = 0.
        enddo
      enddo

!     read in daytime photoionization line data
!     (only n, o, he, n_2, o_2)

      do i = 1,linesuv
        read(phiondt_inp,106) (phionr(i,j), j=1,5)
        sigidt(i,ptn ) = phionr(i,1)
        sigidt(i,pto ) = phionr(i,2)
        sigidt(i,pthe) = phionr(i,3)
        sigidt(i,ptn2) = phionr(i,4)
        sigidt(i,pto2) = phionr(i,5)
      enddo
 106  format(5f7.2)

      do j = 1,nion
        do i = 1,linesuv
          sigidt(i,j) = tm18 * sigidt(i,j) 
        enddo 
      enddo 

!     read in nighttime photoionization line data
!     (only o, n_2, n0, o_2)

      do i = 1,linesnt
        read(phionnt_inp,106) (phionr(i,j), j=1,4)
        sigint(i,pto ) = phionr(i,1)
        sigint(i,ptn2) = phionr(i,2)
        sigint(i,ptno) = phionr(i,3)
        sigint(i,pto2) = phionr(i,4)
      enddo

      do j = 1,nion
        do i = 1,linesnt
          sigint(i,j) = tm18 * sigint(i,j) 
        enddo 
      enddo 

!     read in f74113, ai data and set euv flux
!     (from richards et al., jgr 99, 8981, 1994)

      p  = 0.5 * ( f10p7 + fbar )

      do i = 1,linesuv
        read (euvflux_inp,107) (fluxdat(i,j),j=1,2)
        f74   = fluxdat(i,1)
        ai    = fluxdat(i,2)
        xflux = 1. + ai * ( p - 80.)
        if ( xflux .lt. 0.8 ) xflux = 0.8
        flux(i) = f74 * xflux * 1.e9
!        if ( flux(i) .lt. 0 ) flux(i) = 0.
!        print *,'i,flux',i,flux(i)
      enddo 
 107  format (f6.3,1pe11.4)

!      stop

!     read in angles for nighttime deposition fluxes

      do i = 1,linesnt
        read(thetant_inp,108) (thetant(i,j), j=1,4)
      enddo
 108  format (4f7.1)

!     read in min/max altitude for nighttime deposition fluxes
!       zaltnt(i,1): zmin(i)
!       zaltnt(i,2): zmax(i)

      do i = 1,linesnt
        read(zaltnt_inp,108) (zaltnt(i,j), j=1,2)
      enddo
 109  format (2f7.1)

!     call nighttime euv flux subroutines
!     (lyman beta 1026, he i 584, he ii 304, lyman alpha 1216)
      do j = 1,nf
        call sf1026 ( f1026,1,j )

        call sf584  ( f584 ,2,j )
        call sf304  ( f304 ,3,j )
        call sf1216 ( f1216,4,j )
        do k = 1,91
          do i = 1,nz
            fluxnt(i,j,k,1) = f1026(i,j,k)
            fluxnt(i,j,k,2) = f584 (i,j,k)
            fluxnt(i,j,k,3) = f304 (i,j,k)
            fluxnt(i,j,k,4) = f1216(i,j,k)
          enddo
        enddo
      enddo
!     intialize diagnostic variables to 0

      u1=0.
      u2=0.
      u3=0.
      u4=0.
      u5=0.
      
      
      t1 = 0.
      t2 = 0.
      t3 = 0.




!         close opened files
      call close_input_files

      print *,' finished initialization'


      return
      end

!*******************************************
!*******************************************

!            transprt

!*******************************************
!*******************************************

      subroutine transprt (nfl)

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      use commonsubroutines
      implicit none
      real prod(nz,nion),loss(nz,nion),lossr,&
          phprodr(nz,nion),chrate(nz,nchem),&
          chloss(nz,nion),chprod(nz,nion),relossr(nz,nion)
      real deni_old(nz,nion),te_old(nz),ti_old(nz,nion),vsi_old(nz,nion)
      real tvn(nz)
      real nuin(nz,nion,nneut),&
          nuij(nz,nion,nion),sumnuj(nz,nion)
      real vsin(nz,nion),vsidn(nz,nion),denin(nz,nion)
      real ten(nz),tin(nz,nion)
      !LOCAL VARIABLES
      INTEGER::i,j,nfl,nni,nj,ni

!     calculation of production and loss
!       phprodr: photo production rates
!       chrate:  chemical rates (ichem)
!       chloss:  chemical loss term
!       chprod:  chemical production term
!       relossr: recombination loss rates

!     initialize tvn and gs

!      do i = 1,nz
!        tvn(i) = 0.
!        gs(i)  = 0.
!      enddo

       
      tvn = 0.
      gs = 0.

      do i = 1,nz

        ne(i,nfl)   = 1.
        te_old(i)   = te(i,nfl)
        do j = nion1,nion2
          deni_old(i,j) = deni(i,nfl,j)
          ne(i,nfl)     = ne(i,nfl) + deni(i,nfl,j)
          ti_old(i,j)   = ti(i,nfl,j)
          vsi_old(i,j)  = vsi(i,nfl,j)
        enddo

      enddo
 
       call photprod ( cx,phprodr,nfl   )         ! calculates phprodr
       call chemrate ( chrate,nfl               ) ! calculates chrate
       call chempl   ( chrate,chloss,chprod,nfl ) ! calcualtes chloss,chprod
       call recorate ( relossr,nfl              ) ! calculates relossr

       do i = 1,nz

        do j = nion1,nion2
          prod  (i,j) =  phprodr(i,j) * denn(i,nfl,j)&
                        + chprod(i,j)
          lossr       =  relossr(i,j) * deni(i,nfl,j) * ne(i,nfl)&
                        + chloss(i,j)
          loss (i,j)  =  lossr / deni(i,nfl,j)
        enddo

!       gravity and neutral wind 
!       modified 9/19/05 (MS)

        gs(i)   =  gzero * arg(i,nfl)&
                  * ( re / (re + alts(i,nfl)) ) ** 2

        tvn(i)  = (  v(i,nfl) * athg(i,nfl)&
                  - u(i,nfl) * aphig(i,nfl) )

        tvn(i)    = tvn0 * tvn(i) ! tvn0 used to modify tvn

        u3(i,nfl) = u(i,nfl)
        u4(i,nfl) = v(i,nfl)
        u5(i,nfl) = tvn(i)

      enddo

      call update ( tvn,nuin,sumnuj,nuij,nfl )

      do i = 1,nz
        do nni = nion1,nion2
          sumvsi(i,nfl,nni) = 0.
          do nj = nion1,nion2
          sumvsi(i,nfl,nni) =   sumvsi(i,nfl,nni)&
                          + nuij(i,nni,nj)*vsi(i,nfl,nj)
          enddo
        enddo
      enddo

!     define new arrays for velocity and density

      do ni = nion1,nion2
        do i = 1,nz
          vsin (i,ni) = vsi(i,nfl,ni)
          vsidn(i,ni) = vsid(i,nfl,ni)
          denin(i,ni) = deni(i,nfl,ni)
        enddo
      enddo

!     update variables

      do ni = nion1,nion2

        call vsisolv ( vsin(1,ni),vsidn(1,ni),vsi_old(1,ni)&
                     ,sumnuj(1,ni),nfl )

! compensating filter

        call smoothz ( vsin(1,ni), 1 )

!       put stuff back into velocity array

        do i = 1,nz
          vsi(i,nfl,ni)  = vsin(i,ni)
          vsid(i,nfl,ni) = vsidn(i,ni)
        enddo

        if ( nfl .eq. 1) then
          do i = 1,nz
            vsi(i,1,ni) = vsi(i,2,ni)
          enddo
        endif

        if ( nfl .eq. nf ) then
          do i = 1,nz
            vsi(i,nf,ni) = vsi(i,nf-1,ni)
          enddo
        endif

        call densolv2 ( ni,denin(1,ni)&
            ,prod(1,ni),loss(1,ni),deni_old(1,ni),nfl )

!       put stuff back into density array

        do i = 1,nz
          deni(i,nfl,ni) = denin(i,ni)
        enddo

!       put floor on density

        do i = 1,nz
          deni(i,nfl,ni) = amax1 ( deni(i,nfl,ni), denmin )
        enddo

      enddo

!     define new arrays for temperature

      do ni = nion1,nion2
        do i = 1,nz
          tin(i,ni)  = ti(i,nfl,ni)
        enddo
      enddo

      do i = 1,nz
        ten(i)  = te(i,nfl)
      enddo

!     temperatures (with floors and warnings)
      call etemp  (ten,te_old,phprodr,nfl)
      do i = 1,nz
        te(i,nfl)  = amax1(tn(i,nfl),ten(i))
        te(i,nfl)  = amin1(te(i,nfl),1.e4)
        if ( te(i,nfl) .lt. 0 ) then
          print *,' T(e) negative: i,nfl',i,nfl
          stop
        endif
      enddo

      call htemp  (tin(1,pthp) ,ti_old(1,pthp) ,tvn,nuin,nfl)
      do i = 1,nz
        ti(i,nfl,pthp)  = amax1(tn(i,nfl),tin(i,pthp))
        ti(i,nfl,pthp)  = amin1(ti(i,nfl,pthp),1.e4)
        if ( ti(i,nfl,pthp) .lt. 0 ) then
          print *,' T(H) negative: i,nfl',i,nfl
          stop
        endif
      enddo

      call hetemp (tin(1,pthep),ti_old(1,pthep),tvn,nuin,nfl)
      do i = 1,nz
        ti(i,nfl,pthep)  = amax1(tn(i,nfl),tin(i,pthep))
        ti(i,nfl,pthep)  = amin1(ti(i,nfl,pthep),1.e4)
        if ( ti(i,nfl,pthep) .lt. 0 ) then
          print *,' T(He) negative: i,nfl',i,nfl
          stop
        endif
      enddo

      call otemp  (tin(1,ptop) ,ti_old(1,ptop) ,tvn,nuin,nfl)
      do i = 1,nz
        ti(i,nfl,ptop)  = amax1(tn(i,nfl),tin(i,ptop))
        ti(i,nfl,ptop)  = amin1(ti(i,nfl,ptop),1.e4)
        if ( ti(i,nfl,ptop) .lt. 0 ) then
          print *,' T(O) negative: i,nfl',i,nfl
          stop
        endif
      enddo

      do i = 1,nz
        ti(i,nfl,ptnp )    = ti(i,nfl,ptop)
        ti(i,nfl,ptn2p)    = ti(i,nfl,ptop)
        ti(i,nfl,ptnop)    = ti(i,nfl,ptop)
        ti(i,nfl,pto2p)    = ti(i,nfl,ptop)
      enddo

      return
      end





!*******************************************
!*******************************************

!             xerfcexp

!*******************************************
!*******************************************

        function xerfcexp(x)

        include 'param-1.00.inc'
      use vdrift_model
          t          = 1. / (1 + pas * x)
          xerfcexp   = (   z1 * t&
                        + z2 * t ** 2&
                        + z3 * t ** 3&
                        + z4 * t ** 4&
                        + z5 * t ** 5  )

        return
        end


!*******************************************
!*******************************************

!             EXB

!*******************************************
!*******************************************

       subroutine exb(hrut)

       include 'param-1.00.inc'
       include 'com-1.00.inc'
      use vdrift_model
       real denic(nz,nf,nion)
       real tic(nz,nf,nion)
       real tec(nz,nf)
       real fluxnp(nz,nf,nion),fluxtp(nz,nf,nion)
       real fluxtep(nz,nf)
       real fluxns(nz,nf,nion),fluxts(nz,nf,nion)
       real fluxtes(nz,nf)
       real param(2)

!     define the e x b drift
     flush(6)

      param(1) = day
      param(2) = f10p7
      nzh      = ( nz - 1 ) / 2
!     note: modification of vexb because of field line variation
!           uses cos^3/sqrt(1.+3sin^2) instead of
!           uses sin^3/sqrt(1.+3cos^2) because
!           blats = 0 at the magnetic equator 
!           (as opposed to pi/2 in spherical coordinates)


      hrl = hrut + glons(nzh,1) / 15.
      do while ( hrl .ge. 24. ) 
               hrl = hrl - 24.
         enddo
      call vdrift_model(hrl,glon_in,param,vd,fejer,ve01)

      do j = 1,nf
!        altfac = ( alts(nzh,j) + re ) / re ! L^2 dependence on E x B drift
        altfac = 1.
        vexb0 = 100. * vd * tvexb0 * altfac * altfac ! convert to cm/s 
        do i = 1,nz
          vexb(i,j) = 0.
          blat = blats(i,j) * po180
          vexb(i,j) = vexb0 * &
                      cos(blat) * cos(blat) * cos(blat)/&
                      sqrt( 1. + 3. * sin(blat)* sin(blat) )
        enddo
      enddo

      do j = 1,nf
        do i = 1,nz
          if ( alts(i,j) .lt. alt_crit ) then
            dela      = 20.
            farg      = ( alt_crit - alts(i,j) ) / dela
            vexb(i,j) = vexb(i,j) * exp (-farg*farg)
          endif
        enddo
      enddo

      do j = 1,nf
        vexb(nzp1,j) = vexb(nz,j)
      enddo
      do i = 1,nzp1
        vexb(i,nfp1) = vexb(i,nf) 
      enddo

!     sweep in p-direction

      do j = 1,nfp1
        do i = 1,nz
          vexbp(i,j) = vexb(i,j) * ( vnx(i,j) * xnormp(i,j) +&
                                     vny(i,j) * ynormp(i,j) +&
                                     vnz(i,j) * znormp(i,j)   )
        enddo
      enddo

!     sweep in s-direction
      
      do j = 1,nf
        do i = 1,nzp1
          vexbs(i,j) = vexb(i,j) * ( vnx(i,j) * xnorms(i,j) +&
                                     vny(i,j) * ynorms(i,j) +&
                                     vnz(i,j) * znorms(i,j)   )
        enddo
      enddo

!     output e x b velocities

      do j = 1,nf
        do i = 1,nz
          u1(i,j) = vexbp(i,j)
          u2(i,j) = vexbs(i,j)
        enddo
      enddo

!      dump velocity
!        vot: velocity in theta direction
!        vor: velocity in radial direction
!      needs to be redone ....


!       do ni = nion1,nion2
!         do j = 1,nf
!           do i = 1,nz
!             vot(i,j,ni)   =  vsi(i,j,ni) * cosdips(i,j) + 
!     .                        vexbp(i,j)  * sindips(i,j)
!             vor(i,j,ni)   = -vsi(i,j,ni) * sindips(i,j) + 
!     .                        vexbp(i,j)  * cosdips(i,j)
!           enddo
!         enddo
!       enddo

!      calculate conserved particle number: denic
!      and 'conserved' temperature: tic,tec

       do ni = nion1,nion2
         do j = 1,nf
           do i = 1,nz
             denic(i,j,ni) = deni(i,j,ni) * vol(i,j)
             tic(i,j,ni)   = ti(i,j,ni) * vol(i,j)
           enddo
         enddo
       enddo

       do j = 1,nf
         do i = 1,nz
           tec(i,j)   = te(i,j) * vol(i,j)
         enddo
       enddo

!      calculate flux in p-direction at interface

       do ni = nion1,nion2
         do j = 2,nf
           do i = 1,nz
             if ( vexbp(i,j) .ge. 0 ) then
               fluxnp(i,j,ni) = deni(i,j-1,ni) * vexbp(i,j)
               fluxtp(i,j,ni) = ti(i,j-1,ni)   * vexbp(i,j)
             else
               fluxnp(i,j,ni) = deni(i,j,ni) * vexbp(i,j)
               fluxtp(i,j,ni) = ti(i,j,ni)   * vexbp(i,j)
               if ( j .eq. nf ) then
                 fluxnp(i,j,ni) = fluxnp(i,j-1,ni)&
                                  * areap(i,j-1)/areap(i,j)       
                 fluxtp(i,j,ni) = fluxtp(i,j-1,ni)&
                                  * areap(i,j-1)/areap(i,j)       

                 if ( ni .eq. 1 ) then
                 fluxnp(i,j,ni) = (deni(i,j-1,ni) -&
                                  (areas(i,j-1) / areas(i,j-2)) *&
                                  (deni(i,j-2,ni) - deni(i,j-1,ni)))*&
                                  vexbp(i,j)

                 fluxtp(i,j,ni) = fluxtp(i,j-1,ni) -&
                                  (areas(i,j-1) / areas(i,j-2)) *&
                                  (fluxtp(i,j-2,ni) - fluxtp(i,j-1,ni))

                  endif

!                 fluxnp(i,j,ni) = fluxnp(i,j-1,ni)
!                 fluxtp(i,j,ni) = fluxtp(i,j-1,ni)
               endif
             endif
           enddo
         enddo
       enddo


       do j = 2,nf
         do i = 1,nz
           if ( vexbp(i,j) .ge. 0 ) then
             fluxtep(i,j) = te(i,j-1)   * vexbp(i,j)
           else
             fluxtep(i,j) = te(i,j)   * vexbp(i,j)
             if ( j .eq. nf ) then
                fluxtep(i,j) = fluxtep(i,j-1)&
                              * areap(i,j-1)/areap(i,j)       
                
                if (ni .eq. 1 ) &
                fluxtep(i,j) = fluxtep(i,j-1) -&
                               (areas(i,j-1) / areas(i,j-2)) *&
                               (fluxtep(i,j-2) - fluxtep(i,j-1))

!              fluxtep(i,j) = fluxtep(i,j-1)
             endif
           endif
         enddo
       enddo

!      calculate flux in s-direction at interface

       do ni = nion1,nion2
         do j = 1,nf
           do i = 2,nz
             if ( vexbs(i,j) .ge. 0 ) then
               fluxns(i,j,ni) = deni(i-1,j,ni) * vexbs(i,j)
               fluxts(i,j,ni) = ti(i-1,j,ni)   * vexbs(i,j)
             else
               fluxns(i,j,ni) = deni(i,j,ni) * vexbs(i,j)
               fluxts(i,j,ni) = ti(i,j,ni)   * vexbs(i,j)
             endif
           enddo
         enddo
       enddo

       do j = 1,nf
         do i = 2,nz
           if ( vexbs(i,j) .ge. 0 ) then
             fluxtes(i,j) = te(i-1,j)   * vexbs(i,j)
           else
             fluxtes(i,j) = te(i,j)   * vexbs(i,j)
           endif
         enddo
       enddo

!      update total particle number and density
!      and temperatures

       do ni = nion1,nion2
         do j = 2,nfm1
           do i = 2,nzm1
             denic(i,j,ni) = denic(i,j,ni) &
                             + dt * ( areap(i,j)   * fluxnp(i,j,ni) -&
                                      areap(i,j+1) * fluxnp(i,j+1,ni) )&
                             + dt * ( areas(i,j)   * fluxns(i,j,ni) -&
                                      areas(i+1,j) * fluxns(i+1,j,ni) )
             deni(i,j,ni)  = denic(i,j,ni) / vol(i,j)
             tic(i,j,ni) = tic(i,j,ni) &
                             + dt * ( areap(i,j)   * fluxtp(i,j,ni) -&
                                      areap(i,j+1) * fluxtp(i,j+1,ni) )&
                             + dt * ( areas(i,j)   * fluxts(i,j,ni) -&
                                      areas(i+1,j) * fluxts(i+1,j,ni) )
             ti(i,j,ni)  = tic(i,j,ni) / vol(i,j)
           enddo
         enddo
       enddo

       do j = 2,nfm1
         do i = 2,nzm1
           tec(i,j) = tec(i,j) &
                           + dt * ( areap(i,j)   * fluxtep(i,j) -&
                                    areap(i,j+1) * fluxtep(i,j+1) )&
                           + dt * ( areas(i,j)   * fluxtes(i,j) -&
                                    areas(i+1,j) * fluxtes(i+1,j) )
           te(i,j)  = tec(i,j) / vol(i,j)
         enddo
       enddo

!      fill cells at j = 1 and nf with j = 2 and nfm1

       do ni = nion1,nion2
         do i = 2,nzm1
           deni(i,1,ni)  = deni(i,2,ni)
           deni(i,nf,ni) = deni(i,nfm1,ni)
           ti(i,1,ni)    = ti(i,2,ni)
           ti(i,nf,ni)   = ti(i,nfm1,ni)
         enddo
       enddo

       do i = 2,nzm1
         te(i,1)    = te(i,2)
         te(i,nf)   = te(i,nfm1)
       enddo

!      fill cells at i = 1 and nz with i = 2 and nzm1

       do ni = nion1,nion2
         do j = 1,nf
           deni(1,j,ni)  = deni(2,j,ni)
           deni(nz,j,ni) = deni(nzm1,j,ni)
           ti(1,j,ni)    = ti(2,j,ni)
           ti(nz,j,ni)   = ti(nzm1,j,ni)
         enddo
       enddo

       do j = 1,nf
         te(1,j)    = te(2,j)
         te(nz,j)   = te(nzm1,j)
       enddo

!      E x B drift can cause excessive density
!      if the first field line peak is below alt_crit
!      and the second field line peak is above alt_crit
!      using zero gradient boundary condition -
!      this is a fix

       if ( alts(nz/2,2) .ge. alt_crit .and.&
            alts(nz/2,1) .lt. alt_crit       ) then
         do ni = nion1,nion2
           do i = 1,nz
             deni(i,1,ni)  = denmin
             ti(i,1,ni)    = tn(i,1)
           enddo
         enddo
         do ie = 1,nz
           te(ie,1)  = tn(ie,1)
         enddo
       endif
       

       return
       end


!*******************************************
!*******************************************

!             courant

!*******************************************
!*******************************************

       subroutine courant (hrut)

       include 'param-1.00.inc'
       include 'com-1.00.inc'

       hr24ut = mod(hrut,24.)
       dt00   = dt0
       hr24l = hr24ut + glon_in / 15.
       if ( hr24l .ge. 24. ) hr24l = hr24l - 24.
       if ( hr24l .gt. 18. .and. hr24l .le. 24. ) dt00 = .4 * dt0
       if ( hr24l .gt.  0. .and. hr24l .le.  6. ) dt00 = .4 * dt0

!      parallel motion

       dtnew = 1.e6
       do k = nion1,nion2
         do j = 1,nf
           do i = 1,nz
             dt1 = dels(i,j) / amax1(1.,abs(vsi(i,j,k)))
             if ( dt1 .lt. dtnew ) dtnew = dt1
           enddo
         enddo
       enddo

!      perpendicular motion

       do j = 1,nf
         do i = 1,nz
           dts = xdels(i,j,1) / amax1(1.,abs(vexbs(i,j)))
           dtp = xdelp(i,j,1) / amax1(1.,abs(vexbp(i,j)))
           dt1 = amin1 ( dts,dtp )
           if ( dt1 .lt. dtnew ) dtnew = dt1
         enddo
       enddo

       if ( dtnew .le. .01 ) then
         print *,' Time step too small: dtnew',dtnew
         stop
       elseif ( dtnew .ge. 5e4 ) then
         print *,' Time step too big: dtnew',dtnew
         stop
       endif
       dt = .25 * dtnew
       if ( dtnew/dt .le. 1.0  ) dt = amin1(dt00,dtnew   )
       if ( dtnew/dt .ge. 1.2  ) dt = amin1(dt00,dt * 1.2)

       return
       end




