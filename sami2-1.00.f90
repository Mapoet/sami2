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
      implicit none
      integer::ntm,istep,nfl
      real::hrut,timemax,time,tprnt,tneut

      call initial
  
!     open output files

      call open_uf

      ntm   = 0
      istep = 0
!      call output ( hrinit,ntm,istep )

      close (68)

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

!       parallel transport

        do nfl = nf,1,-1
          call zenith (hrut,nfl)
          call transprt (nfl)
        enddo         

!       perpendicular transport

        call exb(hrut)         

!       time/step advancement

        istep  = istep + 1
        time   = time  + dt
        hrut   = time / sphr + hrinit
        tprnt  = tprnt + dt / sphr
        tneut  = tneut + dt / sphr

        call courant ( hrut )

!       update neutrals

        if ( tneut .ge. 0.25 ) then
          call neutambt (hrut) 
          tneut = 0.
        endif

!       output data

        if ( tprnt .ge. dthr .and. hrut .gt. hrpr+hrinit) then
          ntm      = ntm + 1
          call output ( hrut,ntm,istep )
          tprnt   = 0.
          print *,'ouput - hour = ',hrut
        elseif ( tprnt .ge. dthr ) then
          print *,'no ouput yet - hour = ',hrut
          tprnt   = 0.
        endif
      enddo    ! end time loop

!     close files

      call close_uf

      stop
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

!            zenith

!*******************************************
!*******************************************

       subroutine zenith (hrut,nfl)

       include 'param-1.00.inc'
       include 'com-1.00.inc' 

!      geometric variables

!      sdec: solar zenith angle
!      cx:  cos of the zenith angle

       do i = 1,nz
         hrl = hrut + glons(i,nfl) / 15.
         do while ( hrl .ge. 24. ) 
               hrl = hrl - 24.
         enddo
         sdec         = rtod * asin (  sin (2.*pie*(day-dayve)/sidyr)&
                               * sin (solinc/rtod)             )
         cossdec      = cos ( po180 * sdec )
         sinsdec      = sin ( po180 * sdec )
         clat         = cos ( po180 * glats(i,nfl) )
         slat         = sin ( po180 * glats(i,nfl) )
         cx(i,nfl)    =   slat * sinsdec&
                  - clat * cossdec * cos ( 15.0*po180*hrl )
!         u3(i,nfl)    = cx(i,nfl)
! MS: Since we will be taking acos of this value in photprod, make
! sure that the absolute value does not minutely exceed 1 because of
! round-off error.
         if (abs(abs(cx(i,nfl))-1.) .lt. 1.e-6)&
      cx(i,nfl) = sign(1.,cx(i,nfl))
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

!             output

!*******************************************
!*******************************************

       subroutine output ( hrut,ntm,istep )

       include 'param-1.00.inc'
       include 'com-1.00.inc'

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


!*******************************************
!*******************************************

!             EXB

!*******************************************
!*******************************************

       subroutine exb(hrut)

       include 'param-1.00.inc'
       include 'com-1.00.inc'

       real denic(nz,nf,nion)
       real tic(nz,nf,nion)
       real tec(nz,nf)
       real fluxnp(nz,nf,nion),fluxtp(nz,nf,nion)
       real fluxtep(nz,nf)
       real fluxns(nz,nf,nion),fluxts(nz,nf,nion)
       real fluxtes(nz,nf)
       real param(2)

!     define the e x b drift

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
             if ( dt1 .le. dtnew ) dtnew = dt1
           enddo
         enddo
       enddo

!      perpendicular motion

       do j = 1,nf
         do i = 1,nz
           dts = xdels(i,j,1) / amax1(1.,abs(vexbs(i,j)))
           dtp = xdelp(i,j,1) / amax1(1.,abs(vexbp(i,j)))
           dt1 = amin1 ( dts,dtp )
           if ( dt1 .le. dtnew ) dtnew = dt1
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


!*******************************************
!*******************************************

!             vdrift_model

!*******************************************
!*******************************************

!       ************************************************************
!       ************************************************************

      subroutine vdrift_model(xt,xl,param,y,fejer,ve01)

!       ************************************************************

!       ************************************************************
!       SUBROUTINE CALCULATES EQUATORIAL VERTICAL DRIFT AS DESCRIBED 
!       IN SCHERLIESS AND FEJER, JGR, 104, 6829-6842, 1999
!       ************************************************************

!       INPUT:   XT: SOLAR LOCAL TIME
!                XL: GEOGRAPHIC LONGITUDE (+ EAST)
!               
!             PARAM: 2-DIM ARRAY (DOY,F10.7CM)
!                    DOY     :Day of Year has to run from 1 to 365 (366)
!                    F10.7cm : F10.7cm solar flux
!             
!       OUTPUT:   Y: EQUATORIAL VERTICAL DRIFT

!       ************************************************************

        include 'param-1.00.inc'

!	implicit none

        real ve01
        logical fejer

        real param(2),coeff(624),funct(6)
        real coeff1(312),coeff2(312)
      real xt,xl,y
      real bspl4,bspl4_time,bspl4_long

      integer i,j,ind,il,kk
      integer index_t/13/,dim_t/78/
      integer index_l/8/,dim_l/48/
      integer index/104/,dim/624/
      integer nfunc/6/

        data coeff1/&
        -10.80592, -9.63722,-11.52666,  -0.05716,-0.06288,  0.03564,    &
         -5.80962, -7.86988, -8.50888, -0.05194, -0.05798, -0.00138,&
          2.09876,-19.99896, -5.11393, -0.05370, -0.06585,  0.03171,&
        -10.22653, -3.62499,-14.85924, -0.04023, -0.01190, -0.09656,&
         -4.85180,-26.26264, -6.20501, -0.05342, -0.05174,  0.02419,&
        -13.98936,-18.10416, -9.30503, -0.01969, -0.03132, -0.01984,&
        -18.36633,-24.44898,-16.69001,  0.02033, -0.03414, -0.02062,&
        -20.27621,-16.95623,-36.58234,  0.01445, -0.02044, -0.08297,&
          1.44450,  5.53004,  4.55166, -0.02356, -0.04267,  0.05023,&
          5.50589,  7.05381,  1.94387, -0.03147, -0.03548,  0.01166,&
          3.24165, 10.05002,  4.26218, -0.03419, -0.02651,  0.07456,&
          7.02218,  0.06708,-11.31012, -0.03252, -0.01021, -0.09008,&
         -3.47588, -2.82534, -4.17668, -0.03719, -0.01519,  0.06507,&
         -4.02607,-11.19563,-10.52923, -0.00592, -0.01286, -0.00477,&
        -11.47478, -9.57758,-10.36887,  0.04555, -0.02249,  0.00528,&
        -14.19283,  7.86422, -8.76821,  0.05758, -0.02398, -0.04075,&
         14.58890, 36.63322, 27.57497,  0.01358, -0.02316,  0.04723,&
         12.53122, 29.38367, 21.40356, -0.00071, -0.00553,  0.01484,&
         18.64421, 26.27327, 18.32704,  0.00578,  0.03349,  0.11249,&
          4.53014,  6.15099,  7.41935, -0.02860, -0.00395, -0.08394,&
         14.29422,  9.77569,  2.85689, -0.00107,  0.04263,  0.10739,&
          7.17246,  4.40242, -1.00794,  0.00089,  0.01436,  0.00626,&
          7.75487,  5.01928,  4.36908,  0.03952, -0.00614,  0.03039,&
         10.25556,  8.82631, 24.21745,  0.05492, -0.02968,  0.00177,&
         21.86648, 24.03218, 39.82008,  0.00490, -0.01281, -0.01715,&
         19.18547, 23.97403, 34.44242,  0.01978,  0.01564, -0.02434,&
         26.30614, 14.22662, 31.16844,  0.06495,  0.19590,  0.05631,&
         21.09354, 25.56253, 29.91629, -0.04397, -0.08079, -0.07903,&
         28.30202, 16.80567, 38.63945,  0.05864,  0.16407,  0.07622,&
         22.68528, 25.91119, 40.45979, -0.03185, -0.01039, -0.01206,&
         31.98703, 24.46271, 38.13028, -0.08738, -0.00280,  0.01322,&
         46.67387, 16.80171, 22.77190, -0.13643, -0.05277, -0.01982,&
         13.87476, 20.52521,  5.22899,  0.00485, -0.04357,  0.09970,&
         21.46928, 13.55871, 10.23772, -0.04457,  0.01307,  0.06589,&
         16.18181, 16.02960,  9.28661, -0.01225,  0.14623, -0.01570,&
         18.16289, -1.58230, 14.54986, -0.00375, -0.00087,  0.04991,&
         10.00292, 11.82653,  0.44417, -0.00768,  0.15940, -0.01775,&
         12.15362,  5.65843, -1.94855, -0.00689,  0.03851,  0.04851,&
         -1.25167,  9.05439,  0.74164,  0.01065,  0.03153,  0.02433,&
        -15.46799, 18.23132, 27.45320,  0.00899, -0.00017,  0.03385,&
          2.70396, -0.87077,  6.11476, -0.00081,  0.05167, -0.08932,&
          3.21321, -1.06622,  5.43623,  0.01942,  0.05449, -0.03084,&
         17.79267, -3.44694,  7.10702,  0.04734, -0.00945,  0.11516,&
          0.46435,  6.78467,  4.27231, -0.02122,  0.10922, -0.03331,&
         15.31708,  1.70927,  7.99584,  0.07462,  0.07515,  0.08934,&
          4.19893,  6.01231,  8.04861,  0.04023,  0.14767, -0.04308,&
          9.97541,  5.99412,  5.93588,  0.06611,  0.12144, -0.02124,&
         13.02837, 10.29950, -4.86200,  0.04521,  0.10715, -0.05465,&
          5.26779,  7.09019,  1.76617,  0.09339,  0.22256,  0.09222,&
          9.17810,  5.27558,  5.45022,  0.14749,  0.11616,  0.10418,&
          9.26391,  4.19982, 12.66250,  0.11334,  0.02532,  0.18919,&
         13.18695,  6.06564, 11.87835,  0.26347,  0.02858,  0.14801/

        data coeff2/&
         10.08476,  6.14899, 17.62618,  0.09331,  0.08832,  0.28208,&
         10.75302,  7.09244, 13.90643,  0.09556,  0.16652,  0.22751,&
          6.70338, 11.97698, 18.51413,  0.15873,  0.18936,  0.15705,&
          5.68102, 23.81606, 20.65174,  0.19930,  0.15645,  0.08151,&
         29.61644,  5.49433, 48.90934,  0.70710,  0.40791,  0.26325,&
         17.11994, 19.65380, 44.88810,  0.45510,  0.41689,  0.22398,&
          8.45700, 34.54442, 27.25364,  0.40867,  0.37223,  0.22374,&
         -2.30305, 32.00660, 47.75799,  0.02178,  0.43626,  0.30187,&
          8.98134, 33.01820, 33.09674,  0.33703,  0.33242,  0.41156,&
         14.27619, 20.70858, 50.10005,  0.30115,  0.32570,  0.45061,&
         14.44685, 16.14272, 45.40065,  0.37552,  0.31419,  0.30129,&
          6.19718, 18.89559, 28.24927,  0.08864,  0.41627,  0.19993,&
          7.70847, -2.36281,-21.41381,  0.13766,  0.05113, -0.11631,&
         -9.07236,  3.76797,-20.49962,  0.03343,  0.08630,  0.00188,&
         -8.58113,  5.06009, -6.23262,  0.04967,  0.03334,  0.24214,&
        -27.85742,  8.34615,-27.72532, -0.08935,  0.15905, -0.03655,&
          2.77234,  0.14626, -4.01786,  0.22338, -0.04478,  0.18650,&
          5.61364, -3.82235,-16.72282,  0.26456, -0.03119, -0.08376,&
         13.35847, -6.11518,-16.50327,  0.28957, -0.01345, -0.19223,&
         -5.37290, -0.09562,-27.27889,  0.00266,  0.22823, -0.35585,&
        -15.29676,-18.36622,-24.62948, -0.31299, -0.23832, -0.08463,&
        -23.37099,-13.69954,-26.71177, -0.19654, -0.18522, -0.20679,&
        -26.33762,-15.96657,-42.51953, -0.13575, -0.00329, -0.28355,&
        -25.42140,-14.14291,-21.91748, -0.20960, -0.19176, -0.32593,&
        -23.36042,-23.89895,-46.05270, -0.10336,  0.03030, -0.21839,&
        -19.46259,-21.27918,-32.38143, -0.17673, -0.15484, -0.11226,&
        -19.06169,-21.13240,-34.01677, -0.25497, -0.16878, -0.11004,&
        -18.39463,-16.11516,-19.55804, -0.19834, -0.23271, -0.25699,&
        -19.93482,-17.56433,-18.58818,  0.06508, -0.18075,  0.02796,&
        -23.64078,-18.77269,-22.77715, -0.02456, -0.12238,  0.02959,&
        -12.44508,-21.06941,-19.36011,  0.02746, -0.16329,  0.19792,&
        -26.34187,-19.78854,-24.06651, -0.07299, -0.03082, -0.03535,&
        -10.71667,-26.04401,-16.59048,  0.02850, -0.09680,  0.15143,&
        -18.40481,-23.37770,-16.31450, -0.03989, -0.00729, -0.01688,&
         -9.68886,-20.59304,-18.46657,  0.01092, -0.07901,  0.03422,&
         -0.06685,-19.24590,-29.35494,  0.12265, -0.24792,  0.05978,&
        -15.32341, -9.07320,-13.76101, -0.17018, -0.15122, -0.06144,&
        -14.68939,-14.82251,-13.65846, -0.11173, -0.14410, -0.07133,&
        -18.38628,-18.94631,-19.00893, -0.08062, -0.14481, -0.12949,&
        -16.15328,-17.40999,-14.08705, -0.08485, -0.06896, -0.11583,&
        -14.50295,-16.91671,-25.25793, -0.06814, -0.13727, -0.12213,&
        -10.92188,-14.10852,-24.43877, -0.09375, -0.11638, -0.09053,&
        -11.64716,-14.92020,-19.99063, -0.14792, -0.08681, -0.12085,&
        -24.09766,-16.14519, -8.05683, -0.24065, -0.05877, -0.23726,&
        -25.18396,-15.02034,-15.50531, -0.12236, -0.09610, -0.00529,&
        -15.27905,-19.36708,-12.94046, -0.08571, -0.09560, -0.03544,&
         -7.48927,-16.00753,-13.02842, -0.07862, -0.10110, -0.05807,&
        -13.06383,-27.98698,-18.80004, -0.05875, -0.03737, -0.11214,&
        -13.67370,-16.44925,-16.12632, -0.07228, -0.09322, -0.05652,&
        -22.61245,-21.24717,-18.09933, -0.05197, -0.07477, -0.05235,&
        -27.09189,-21.85181,-20.34676, -0.05123, -0.05683, -0.07214,&
        -27.09561,-22.76383,-25.41151, -0.10272, -0.02058, -0.16720/

        do i = 1,312
          coeff(i)     = coeff1(i)
          coeff(i+312) = coeff2(i)
        enddo

        xt = mod(xt,24.)

!       sinusoid e x b model

        if ( .not. fejer ) then
          y = ve01 * sin ( 2 * pie * ( xt - 7. ) / 24. )
          return
        endif

!       fejer-scherliess e x b model

      call g(param,funct,xl)

!       **********************************
      y=0.
!       **********************************
      do i=1,index_t
        do il=1,index_l
          kk=index_l*(i-1)+il
          do j=1,nfunc
             ind=nfunc*(kk-1)+j
             bspl4=bspl4_time(i,xt)*bspl4_long(il,xl)
               y=y+bspl4*funct(j)*coeff(ind)
          end do
          end do
         end do

      end
!------------------------------------------------------------------



!       *************************************************
!       *************************************************
        real function bspl4_time(i,x1)
!       *************************************************
      implicit none 

      integer i,order/4/,j,k
      real t_t(0:39)
      real x,b(20,20),x1

        data t_t/&
                0.00,2.75,4.75,5.50,6.25,&
                7.25,10.00,14.00,17.25,18.00,&
                18.75,19.75,21.00,24.00,26.75,&
                28.75,29.50,30.25,31.25,34.00,&
                38.00,41.25,42.00,42.75,43.75,&
                45.00,48.00,50.75,52.75,53.50,&
                54.25,55.25,58.00,62.00,65.25,&
                66.00,66.75,67.75,69.00,72.00/
       
      x=x1
        if(i.ge.0) then
          if (x.lt.t_t(i-0)) then
           x=x+24
        end if
      end if
      do j=i,i+order-1
         if(x.ge.t_t(j).and.x.lt.t_t(j+1)) then
             b(j,1)=1
         else
             b(j,1)=0
         end if
      end do

      do j=2,order
           do k=i,i+order-j
      	b(k,j) = ( x - t_t(k) ) / ( t_t(k+j-1) - t_t(k) ) * &
                         b(k,j-1)
      	b(k,j) = b(k,j) + &
                         ( t_t(k+j)-x ) / ( t_t(k+j) - t_t(k+1) ) *&
                          b(k+1,j-1)
           end do
      end do

      bspl4_time=b(i,order)
      end
!------------------------------------------------------------------



!       *************************************************
!       *************************************************
        real function bspl4_long(i,x1)
!       *************************************************
      implicit none 

      integer i,order/4/,j,k
      real t_l(0:24)
      real x,b(20,20),x1

        data t_l/&
                0,10,100,190,200,250,280,310,&
                360,370,460,550,560,610,640,670,&
                720,730,820,910,920,970,1000,1030,1080/
       
      x=x1
        if(i.ge.0) then
          if (x.lt.t_l(i-0)) then
           x=x+360
        end if
      end if
      do j=i,i+order-1
         if(x.ge.t_l(j).and.x.lt.t_l(j+1)) then
             b(j,1)=1
         else
             b(j,1)=0
         end if
      end do

      do j=2,order
           do k=i,i+order-j
      	b(k,j)=(x-t_l(k))/(t_l(k+j-1)-t_l(k))*b(k,j-1)
      	b(k,j)=b(k,j)+(t_l(k+j)-x)/(t_l(k+j)-t_l(k+1))*&
                       b(k+1,j-1)
           end do
      end do

      bspl4_long=b(i,order)
      end
!------------------------------------------------------------------



!       *************************************************
!       *************************************************
        subroutine g(param,funct,x)
!       *************************************************
      implicit none

        integer i
      real param(2),funct(6)
      real x,a,sigma,gauss,flux,cflux

!       *************************************************
      flux=param(2)
        if(param(2).le.75)  flux=75.
        if(param(2).ge.230) flux=230.
      cflux=flux

      a=0.
        if((param(1).ge.120).and.(param(1).le.240)) a=170.
        if((param(1).ge.120).and.(param(1).le.240)) sigma=60
        if((param(1).le.60).or.(param(1).ge.300)) a=170.
        if((param(1).le.60).or.(param(1).ge.300)) sigma=40

      if((flux.le.95).and.(a.ne.0)) then
       gauss=exp(-0.5*((x-a)**2)/sigma**2)
         cflux=gauss*95.+(1-gauss)*flux
        end if
!       *************************************************

!       *************************************************
        do i=1,6
         funct(i)=0.
        end do
!       *************************************************

!       *************************************************
        if((param(1).ge.135).and.(param(1).le.230)) funct(1)=1
        if((param(1).le.45).or.(param(1).ge.320)) funct(2)=1
        if((param(1).gt.75).and.(param(1).lt.105)) funct(3)=1
        if((param(1).gt.260).and.(param(1).lt.290)) funct(3)=1
!       *************************************************

        if((param(1).ge.45).and.(param(1).le.75)) then  ! W-E
       funct(2)=1.-(param(1)-45.)/30.
       funct(3)=1-funct(2)
        end if
        if((param(1).ge.105).and.(param(1).le.135)) then  ! E-S
       funct(3)=1.-(param(1)-105.)/30.
       funct(1)=1-funct(3)
        end if
        if((param(1).ge.230).and.(param(1).le.260)) then  ! S-E
       funct(1)=1.-(param(1)-230.)/30.
       funct(3)=1-funct(1)
        end if
        if((param(1).ge.290).and.(param(1).le.320)) then  ! E-W
       funct(3)=1.-(param(1)-290.)/30.
       funct(2)=1-funct(3)
        end if

!       *************************************************
        funct(4)=(cflux-140)*funct(1)
        funct(5)=(cflux-140)*funct(2)
        funct(6)=(flux-140)*funct(3)
!       *************************************************

      end
!------------------------------------------------------------------



