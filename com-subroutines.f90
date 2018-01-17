module commonsubroutines
contains

!*******************************************
!*******************************************

!            msistim

!*******************************************
!*******************************************


      subroutine msistim ( iyr,iday,hrl,glong,iyd,secut )

!     msistim calculates time parameters for the
!     nrlmsise00 neutral atmosphere model.

!     the arguments are defined as follows:

!       iyr    the julian year
!       iday   the day of the year
!       hr     the local time in hours
!       glong  the geocentric longitude in degrees east
!       iyd    the year and day in the form yydd
!       secut  the universal time in seconds
      implicit none
      INTEGER::iyr
      INTEGER::iday
      REAL::hrl
      REAL::glong
      INTEGER::iyd
      REAL::secut

      REAL::hrut
      iyd    = 1000 * mod(iyr,100) + iday
      hrut   = hrl - glong /15.

      do while ( hrut .lt. 0.  )
        hrut = hrut + 24.
      enddo

      do while ( hrut.ge. 24. )
        hrut = hrut - 24.
      enddo

      secut  = hrut * 3600.

      return
      end


!*******************************************
!*******************************************

!            neutambt

!*******************************************
!*******************************************


!     calculate neutral densities and temperature
!     from nrlmsise00

      subroutine neutambt (hrut)


      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit  none
      REAL,intent(in)::hrut

     !LOCAL VARIABLES
      INTEGER::i
      INTEGER::j
      INTEGER::iyd
      REAL::hrl
      REAL::SEC
      REAL,DIMENSION(9):: d
      REAL,DIMENSION(2):: temp
      REAL,DIMENSION(2):: whm93
      REAL,DIMENSION(2):: app
 
!     no obtained from eq. (128) - bailey and balan (red book)

!     neutral density and temperature

!     input:
!        iyd - year and day as yyddd
!        sec - ut(sec)
!        alt - altitude(km) (greater than 85 km)
!        glat - geodetic latitude(deg)
!        glong - geodetic longitude(deg)
!        stl - local apparent solar time(hrs)
!        f107a - 3 month average of f10.7 flux
!        f107 - daily f10.7 flux for previous day
!        ap - magnetic index(daily) or when sw(9)=-1. :
!           - array containing:
!             (1) daily ap
!             (2) 3 hr ap index for current time
!             (3) 3 hr ap index for 3 hrs before current time
!             (4) 3 hr ap index for 6 hrs before current time
!             (5) 3 hr ap index for 9 hrs before current time
!             (6) average of eight 3 hr ap indicies from 12 to 33 hrs prior
!                    to current time
!             (7) average of eight 3 hr ap indicies from 36 to 59 hrs prior
!                    to current time
!        mass - mass number (only density for selected gas is
!                 calculated.  mass 0 is temperature.  mass 48 for all.
!     output:
!        d(1) - he number density(cm-3)
!        d(2) - o number density(cm-3)
!        d(3) - n2 number density(cm-3)
!        d(4) - o2 number density(cm-3)
!        d(5) - ar number density(cm-3)
!        d(6) - total mass density(gm/cm3)
!        d(7) - h number density(cm-3)
!        d(8) - n number density(cm-3)
!        d(9) - anomalous O (see msis)
!        t(1) - exospheric temperature
!        t(2) - temperature at alt

      do j = 1,nf
        do i = 1,nz
          hrl = hrut + glons(i,j) / 15.
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

!     neutral winds

!        iyd - year and day as yyddd
!        sec - ut(sec)  (not important in lower atmosphere)
!        alt - altitude(km)
!        glat - geodetic latitude(deg)
!        glong - geodetic longitude(deg)
!        stl - local apparent solar time(hrs)
!        f107a - 3 month average of f10.7 flux (use 150 in lower atmos.)
!        f107 - daily f10.7 flux for previous day ( " )
!        ap - two element array with
!             ap(1) = magnetic index(daily) (use 4 in lower atmos.)
!             ap(2)=current 3hr ap index (used only when sw(9)=-1.)
!     note:  ut, local time, and longitude are used independently in the
!            model and are not of equal importance for every situation.
!            for the most physically realistic calculation these three
!            variables should be consistent.
!      output
!        w(1) = meridional (m/sec + northward)
!        w(2) = zonal (m/sec + eastward)

      do j = 1,nf
        do i = 1,nz
          app(1)   = ap
          app(2)   = ap
          hrl = hrut + glons(i,j) / 15.
          do while ( hrl .ge. 24. ) 
               hrl = hrl - 24.
         enddo
          call msistim ( int(year),int(day),hrl,glons(i,j),iyd,sec )
          call gws5 ( iyd,sec,alts(i,j),glats(i,j),glons(i,j),&
                     hrl,fbar,f10p7,app,whm93                )
          v(i,j)   = 100. * whm93(1) ! convert to cm/sec
          u(i,j)   = 100. * whm93(2) ! convert to cm/sec
        enddo
      enddo

      return
      end

!*******************************************
!*******************************************

!            f1026

!*******************************************
!*******************************************

!     subroutine to calculate the nighttime flux of
!     lyman beta (1026) (note: line = 1)

      subroutine sf1026 ( f,line,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none
      REAL,DIMENSION(nz,nf,91):: f
      INTEGER::line
      INTEGER::nfl

      !LOCAL VARIABLES
      INTEGER::imax
      INTEGER::i,k,k90,ji,ki,j,jip1,kip1
      REAL::delk,flog

      imax = 1

!     determine f for the 4 known values of theta

      do i = 1,nz
        if ( alts(i,nfl) .lt. zaltnt(line,1) ) then
          do k = 1,4
            f( i,nfl,int(thetant(line,k))+1-90 ) = 1.
          enddo
        elseif ( zaltnt(line,1) .le. alts(i,nfl) .and.&
                alts(i,nfl) .le. zaltnt(line,2)       ) then
          f( i,nfl,int(thetant(line,1))+1-90 ) =&
            1.4e8 * tanh ( (alts(i,nfl) - 90.) / 50. )
          f( i,nfl,int(thetant(line,2))+1-90 ) =&
            3.8e7 * tanh ( (alts(i,nfl) - 90.) / 50. )
          f( i,nfl,int(thetant(line,3))+1-90 ) =&
            1.4e7 * tanh ( (alts(i,nfl) - 93.) / 55. )
          f( i,nfl,int(thetant(line,4))+1-90 ) =&
            9.2e6 * tanh ( (alts(i,nfl) - 94.) / 55. )
          imax = i
        else
          do k = 1,4
            f( i,nfl,   int(thetant(line,k))+1-90 ) =&
           f( imax,nfl,int(thetant(line,k))+1-90 )
          enddo
        endif
      enddo

      do k = 1,4
        do i = 1,nz
          f( i,nfl,int(thetant(line,k))+1-90 ) =&
         amax1 ( 1., f( i,nfl,int(thetant(line,k))+1-90 ) )
        enddo
      enddo

!     now interpolate to all valuse of theta (90 - 180)

      do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
          if ( k90 .gt. int(thetant(line,j)) ) then
            ji = j
            ki = int(thetant(line,ji))
          endif
        enddo
        jip1 = ji + 1
        kip1 = int(thetant(line,jip1))
        delk = float (   int(thetant(line,jip1))&
                      - int(thetant(line,ji  )) )
        do i = 1,nz
          flog =   alog10(f(i,nfl,ki+1-90))&
                + (k90 - ki) / delk&
                             * (  alog10(f(i,nfl,kip1+1-90))&
                                - alog10(f(i,nfl,ki  +1-90)) )
          f(i,nfl,k) = 10 ** flog
        enddo
      enddo

      return
      end

!*******************************************
!*******************************************

!            f584

!*******************************************
!*******************************************

!     subroutine to calculate the nighttime flux of
!     he i (584) (note: line = 2)

      subroutine sf584 ( f,line,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none
      REAL,DIMENSION(nz,nf,91):: f
      INTEGER::line
      INTEGER::nfl

      !LOCAL VARIABLES
      INTEGER::imax
      INTEGER::i,k,k90,ji,ki,j,jip1,kip1
      REAL::delk,flog

      imax = 1

!     determine f for the 4 known values of theta

      do i = 1,nz
        if ( alts(i,nfl) .lt. zaltnt(line,1) ) then
          do k = 1,4
            f( i,nfl,int(thetant(line,k))+1-90 ) = 1.
          enddo
        elseif ( zaltnt(line,1) .le. alts(i,nfl) .and.&
                alts(i,nfl) .le. zaltnt(line,2)       ) then
          f( i,nfl,int(thetant(line,1))+1-90 ) =&
            1.85e5 * ( alts(i,nfl) - 170. ) ** 1.20
          f( i,nfl,int(thetant(line,2))+1-90 ) =&
            2.60e4 * ( alts(i,nfl) - 170. ) ** 1.25
          f( i,nfl,int(thetant(line,3))+1-90 ) =&
            2.60e3 * ( alts(i,nfl) - 170. ) ** 1.20
          f( i,nfl,int(thetant(line,4))+1-90 ) =&
            2.60e2 * ( alts(i,nfl) - 170. ) ** 1.20
          imax = i
        else
          do k = 1,4
            f( i   ,nfl,int(thetant(line,k))+1-90 ) =&
           f( imax,nfl,int(thetant(line,k))+1-90 )
          enddo
        endif
      enddo

      do k = 1,4
        do i = 1,nz
          f( i,nfl,int(thetant(line,k))+1-90 ) =&
         amax1 ( 1., f( i,nfl,int(thetant(line,k))+1-90 ) )
        enddo
      enddo

!     now interpolate to all valuse of theta (90 - 180)
!     set f(i,nfl,theta=180) = 1.

      do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
          if ( k90 .gt. int(thetant(line,j)) ) then
            ji = j
            ki = int(thetant(line,ji))
          endif
        enddo
        if ( ji .ne. 4 ) then
          jip1 = ji + 1
          kip1 = int(thetant(line,jip1))
          delk = float (   int(thetant(line,jip1))&
                        - int(thetant(line,ji  )) )
          do i = 1,nz
            flog =   alog10(f(i,nfl,ki+1-90))&
                  + (k90 - ki) / delk&
                               * (  alog10(f(i,nfl,kip1+1-90))&
                                  - alog10(f(i,nfl,ki  +1-90)) )
            f(i,nfl,k) = 10 ** flog
          enddo
        else
          delk = float (   180&
                        - int(thetant(line,ji  )) )
          do i = 1,nz
            flog =   alog10(f(i,nfl,ki+1-90))&
                  + (k90 - ki) / delk&
                               * (  alog10(1.)&
                                  - alog10(f(i,nfl,ki  +1-90)) )
            f(i,nfl,k) = 10 ** flog
          enddo
        endif
      enddo

      return
      end

!*******************************************
!*******************************************

!            f304

!*******************************************
!*******************************************

!     subroutine to calculate the nighttime flux of
!     he ii (304) (note: line = 3)

      subroutine sf304 ( f,line,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none

      REAL,DIMENSION(nz,nf,91):: f
      INTEGER::line
      INTEGER::nfl

      !LOCAL VARIABLES
      INTEGER::imax
      INTEGER::i,k,k90,ji,ki,j,jip1,kip1
      REAL::delk,flog

      imax = 1

!     determine f for the 4 known values of theta

      do i = 1,nz
        if ( alts(i,nfl) .lt. zaltnt(line,1) ) then
          do k = 1,4
            f( i,nfl,int(thetant(line,k))+1-90 ) = 1.
          enddo
        elseif ( zaltnt(line,1) .le. alts(i,nfl) .and.&
                alts(i,nfl) .le. zaltnt(line,2)       ) then
          f( i,nfl,int(thetant(line,1))+1-90 ) =&
            3.8e6 * tanh ( (alts(i,nfl) - 138.) / 80. )
          f( i,nfl,int(thetant(line,2))+1-90 ) =&
            3.0e6 * tanh ( (alts(i,nfl) - 138.) / 80. )
          f( i,nfl,int(thetant(line,3))+1-90 ) =&
            2.5e6 * tanh ( (alts(i,nfl) - 138.) / 80. )
          f( i,nfl,int(thetant(line,4))+1-90 ) =&
            2.5e6 * tanh ( (alts(i,nfl) - 138.) / 80. )
          imax = i
        else
          do k = 1,4
            f( i,   nfl,int(thetant(line,k))+1-90 ) =&
           f( imax,nfl,int(thetant(line,k))+1-90 )
          enddo
        endif
      enddo

      do k = 1,4
        do i = 1,nz
          f( i,nfl,int(thetant(line,k))+1-90 ) =&
         amax1 ( 1., f( i,nfl,int(thetant(line,k))+1-90 ) )
        enddo
      enddo

!     now interpolate to all valuse of theta (90 - 180)
!     set f(i,nfl,theta=180) = 1.

      do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
          if ( k90 .gt. int(thetant(line,j)) ) then
            ji = j
            ki = int(thetant(line,ji))
          endif
        enddo
        if ( ji .ne. 4 ) then
          jip1 = ji + 1
          kip1 = int(thetant(line,jip1))
          delk = float (   int(thetant(line,jip1))&
                        - int(thetant(line,ji  )) )
          do i = 1,nz
            flog =   alog10(f(i,nfl,ki+1-90))&
                  + (k90 - ki) / delk&
                               * (  alog10(f(i,nfl,kip1+1-90))&
                                  - alog10(f(i,nfl,ki  +1-90)) )
            f(i,nfl,k) = 10 ** flog
          enddo
        else
          delk = float (   180&
                        - int(thetant(line,ji  )) )
          do i = 1,nz
            flog =   alog10(f(i,nfl,ki+1-90))&
                  + (k90 - ki) / delk&
                               * (  alog10(1.)&
                                  - alog10(f(i,nfl,ki  +1-90)) )
            f(i,nfl,k) = 10 ** flog
          enddo
        endif
      enddo

      return
      end

!*******************************************
!*******************************************

!            f1216

!*******************************************
!*******************************************

!     subroutine to calculate the nighttime flux of
!     lyman alpha (1216) (note: line = 4)

      subroutine sf1216 ( f,line,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none

      REAL,DIMENSION(nz,nf,91):: f
      INTEGER::line
      INTEGER::nfl

      !LOCAL VARIABLES
      INTEGER::imax
      INTEGER::i,k,k90,ji,ki,j,jip1,kip1
      REAL::delk,flog

      imax = 1

!     determine f for the 4 known values of theta
      do i = 1,nz
        if ( alts(i,nfl) .lt. zaltnt(line,1) ) then
          do k = 1,4
            f( i,nfl,int(thetant(line,k))+1-90 ) = 1.
          enddo
        elseif ( zaltnt(line,1) .le. alts(i,nfl) .and.&
                alts(i,nfl) .le. zaltnt(line,2)       ) then
          f( i,nfl,int(thetant(line,1))+1-90 ) =&
            1.2e10 * tanh ( (alts(i,nfl) - 80.) / 50. ) + 3.e9
          f( i,nfl,int(thetant(line,2))+1-90 ) =&
            4.0e9  * tanh ( (alts(i,nfl) - 80.) / 50. ) + 1.e9
          f( i,nfl,int(thetant(line,3))+1-90 ) =&
            2.0e9  * tanh ( (alts(i,nfl) - 65.) / 50. ) + 1.e8
          f( i,nfl,int(thetant(line,4))+1-90 ) =&
            1.5e9  * tanh ( (alts(i,nfl) - 75.) / 50. ) + 1.e8
          imax = i
        else
          do k = 1,4
            f( i,   nfl,int(thetant(line,k))+1-90 ) =&
           f( imax,nfl,int(thetant(line,k))+1-90 )
          enddo
        endif
      enddo

      do k = 1,4
        do i = 1,nz
          f( i,nfl,int(thetant(line,k))+1-90 ) =&
         amax1 ( 1., f( i,nfl,int(thetant(line,k))+1-90 ) )
        enddo
      enddo
!     now interpolate to all valuse of theta (90 - 180)

      do k = 1,91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1,4
          if ( k90 .gt. int(thetant(line,j)) ) then
            ji = j
            ki = int(thetant(line,ji))
          endif
        enddo
        jip1 = ji + 1
        kip1 = int(thetant(line,jip1))
        delk = float (   int(thetant(line,jip1))&
                      - int(thetant(line,ji  )) )
        do i = 1,nz
          flog =   alog10(f(i,nfl,ki+1-90))&
                + (k90 - ki) / delk&
                             * (  alog10(f(i,nfl,kip1+1-90))&
                                - alog10(f(i,nfl,ki  +1-90)) )
          f(i,nfl,k) = 10 ** flog
        enddo
      enddo

      return
      end


!*******************************************
!*******************************************

!            photprod

!*******************************************
!*******************************************

!     photoproduction rates

      subroutine photprod ( cxl,phprodr,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none

      REAL atm_chapman

      REAL,DIMENSION(nz,nf)::cxl
      REAL,DIMENSION(nz,nion):: phprodr(nz,nion)
      INTEGER::nfl

      !LOCAL VARIABLES
      real xmass(3)
      integer idx(3)

      REAL::hcof,rp,rp2,xscale,hscale,exa,coschi,y1,ch1,flx,ang
      INTEGER::iz,i,j,l,itheta
!     scale height of neutral atmosphere

      hcof = 1.e-5 * bolt / ( gzero * amu * re ** 2 )

      do iz = 1,nz

         coschi = cxl(iz,nfl)

         do j = nion1,nion2
            phprodr ( iz,j ) = 0.
         enddo

!     only consider o, n2, o2 for absorption

         idx(1) = pto
         idx(2) = ptn2
         idx(3) = pto2

         rp    = alts(iz,nfl) + re
         rp2   = rp * rp

!         if ( coschi .ge. 0. ) then ! sun is up
         if ( coschi .ge. coschicrit(iz,nfl) ) then ! sun is up

!     daytime deposition

            do i = 1,3
               hscale   = hcof * tn(iz,nfl) * rp2 / amn(idx(i))
               xscale   = rp / hscale
               y1       = sqrt ( .5 * xscale ) * abs(coschi)
               ch1      = atm_chapman(xscale,rtod*acos(coschi))
              if (ch1 .gt. 1.e22) ch1 = 1.e22
              xmass(i) = denn(iz,nfl,idx(i)) * hscale * ch1 * 1.e5
            enddo

            do l=1,linesuv
               exa =   xmass(1) * sigabsdt(l,1)&
                   + xmass(2) * sigabsdt(l,2)&
                   + xmass(3) * sigabsdt(l,3)
               if(exa .gt. 35.) exa = 35.
               flx = flux(l) * exp(-exa)
               do j=nion1,nion2
                  phprodr(iz,j) = phprodr(iz,j) + sigidt(l,j) * flx
               enddo
            enddo

!     nighttime deposition

         else                   ! sun is dowm

            ang    = acos ( coschi )
            itheta = nint ( ang / po180 ) - 90
            itheta = int ( amax1 ( float(itheta), 1. ) )
            do l = 1,linesnt
               do j=nion1,nion2
                  phprodr(iz,j) =   phprodr(iz,j)&
                      + sigint(l,j) * fluxnt(iz,nfl,itheta,l)
               enddo
            enddo

         endif
      enddo

      return
      end


!*******************************************
!*******************************************

!            chemrate

!*******************************************
!*******************************************

!     chemical producation and loss rates
!     bb: bailley and balan (red book, 1996)

      subroutine chemrate ( chrate,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none

      REAL,DIMENSION(nz,nchem)::chrate
      INTEGER::nfl

      !LOCAL VARIABLES
      INTEGER::iz
      REAL::ti300o


      do iz = 1,nz

      ti300o = ti(iz,nfl,ptop) / 300.

      chrate (iz,1) = 2.2e-11&
                 * sqrt( ti(iz,nfl,pthp) )       ! h+ + o --> o+ + h (bb)

      chrate (iz,2) = 3.5e-10                        ! he+ + n2 --> n2+ + he (bb)

      chrate (iz,3) = 8.5e-10                        ! he+ + n2 --> n+ + n + he (schunk)

      chrate (iz,4) = 8.0e-10                        ! he+ + o2 --> o+ + o + he (bb)

      chrate (iz,5) = 2.0e-10                        ! he+ + o2 --> o2+ + he

      chrate (iz,6) = 2.0e-10                        ! n+ + o2 --> no+ + o  (schunk)

      chrate (iz,7) = 4.0e-10                        ! n+ + o2 --> o2+ + n(2d) (schunk)

      chrate (iz,8) = 1.0e-12                        ! n+ + 0 --> o+ + n

      chrate (iz,9) = 2.0e-11                        ! n+ + no --> no+ + o (schunk)

      chrate(iz,10) = 2.5e-11&
                  * sqrt( tn(iz,nfl) )           ! o+ + h --> h+ + o   (bb)

      chrate(iz,11) = 1.533e-12 -&
                   5.920e-13 * ti300o +&
                   8.600e-14 * ti300o ** 2                    ! o+ + n2 --> no+ + n (bb)
      if ( ti(iz,nfl,ptop) .gt. 1700 )&
        chrate(iz,11) = 2.730e-12 -&
                      1.155e-12 * ti300o +&
                      1.483e-13 * ti300o ** 2

      chrate(iz,12) = 2.820e-11 -&
                   7.740e-12 * ti300o +&
                   1.073e-12 * ti300o ** 2 -&
                   5.170e-14 * ti300o ** 3 +&
                   9.650e-16 * ti300o ** 4                    ! o+ + o2 --> o2+ + o

      chrate(iz,13) = 1.0e-12                        ! o+ + no --> no+ + o

      chrate(iz,14) = 1.4e-10 / ti300o ** .44        ! n2+ + o --> no+ + n(2d) (bb)

      chrate(iz,15) = 5.0e-11 / sqrt( ti300o )       ! n2+ + o2 --> o2+ + n2 (schunk)

      chrate(iz,16) = 1.0e-14                        ! n2+ + o2 --> no+ + no

      chrate(iz,17) = 3.3e-10                        ! n2+ + no --> no+ + n2 (schunk)

      chrate(iz,18) = 1.2e-10                        ! o2+ + n --> no+ + o (schunk)

      chrate(iz,19) = 2.5e-10                        ! o2+ + n(2d) --> n+ + o2

      chrate(iz,20) = 4.4e-10                        ! o2+ + no --> no+ + o2 (bb)

      chrate(iz,21) = 5.0e-16                        ! o2+ + n2 --> no+ + no (schunk)


      enddo

      return
      end


!*******************************************
!*******************************************

!          chempl

!*******************************************
!*******************************************

!     chemical loss (chloss) and production (chprod)

!     chrate: chemical reaction rates calculated in chemrate
!     ichem: input data file showing loss, neutral, production
!            species for each reaction

      subroutine chempl ( chrate,chloss,chprod,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none

      REAL,DIMENSION(nz,nchem)::chrate
      REAL,DIMENSION(nz,nion)::chloss
      REAL,DIMENSION(nz,nion)::chprod
      INTEGER::nfl

      !LOCAL VARIABLES
      INTEGER::i,iz,k,il,in,ip
      REAL::chem,tdeni
      do i = nion1,nion2
        do iz = 1,nz
          chloss(iz,i)   = 0.
          chprod(iz,i)   = 0.
        enddo
      enddo

      do k = 1,nchem
        il = ichem(k,1) ! ion species (reacting) loss
        in = ichem(k,2) ! neutral species reacting
        ip = ichem(k,3) ! ion species produced
        do iz = 1,nz
           chem  = denn(iz,nfl,in) * chrate(iz,k)
           tdeni = deni(iz,nfl,il) * chem
           chloss(iz,il) = tdeni + chloss(iz,il)
           chprod(iz,ip) = tdeni + chprod(iz,ip)
        enddo
      enddo

      return
      end


!*******************************************
!*******************************************

!            recorate

!*******************************************
!*******************************************

!     recombination rates
!     bb: bailley and balan (red book, 1996)

      subroutine recorate ( relossr,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none

      REAL,DIMENSION(nz,nion):: relossr
      INTEGER::nfl

      !LOCAL VARIABLES
      INTEGER::iz
      REAL::te300

      do iz = 1,nz

        te300 = te(iz,nfl) / 300.

        relossr(iz,pthp)  = 4.43e-12 / te300 ** .7
        relossr(iz,pthep) = relossr(iz,pthp)
        relossr(iz,ptnp)  = relossr(iz,pthp)
        relossr(iz,ptop)  = relossr(iz,pthp)
        relossr(iz,ptn2p) = 1.8e-7 / te300 ** 0.39     !   (schunk)
        relossr(iz,ptnop) = 4.2e-7 / te300 ** 0.85     !   (bb)
        relossr(iz,pto2p) = 1.6e-7 / te300 ** 0.55     !   (schunk)

      enddo

      return
      end


!*******************************************
!*******************************************

!            update

!*******************************************
!*******************************************

      subroutine update ( tvn,nuin,sumnuj,nuij,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none

      REAL,DIMENSION(nz):: tvn
      REAL,DIMENSION(nz,nion,nneut):: nuin
      REAL,DIMENSION(nz,nion)::sumnuj
      REAL,DIMENSION(nz,nion,nion)::nuij
      INTEGER::nfl

      !LOCAL VARIABLES
      real nuint(nz,nion)
      real nufacij,nufacin
      real k0,mi
      INTEGER::ni,nn,i,nj
      REAL::teff,fac,tfactor,amuf,amimn,alame1,alame2,alame,alam,amufac,term1,term2,term3,pip,pim,pep,pem,denid,dened
!     ion-neutral collision frequency

!     initialize everything to 0

      nuin = 0.
      nuint = 0.

!     collision frequencies/factors

!     hydrogen (H)

      ni = pthp
      do nn = 1,nneut
         do i = 1,nz
          if ( nn .eq. pto ) then
!     MS: According to both the SAMI2 paper and the red book the
!     temperature used here should be the H+ temperature, and not the
!     hybrid used in the other terms. I've changed that.

!            teff    = 0.5 * ( ti(i,nfl,ni) + tn(i,nfl) ) ! original
            teff    = ti(i,nfl,ni)
            fac     = ( 1.00 - .047 * alog10(teff) ) ** 2
            tfactor = sqrt(teff) * fac
            nuin(i,ni,nn)  = 6.61e-11 * denn(i,nfl,nn) * tfactor
          else
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          endif
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     helium (He)

      ni = pthep
      do nn = 1,nneut
         do i = 1,nz
          amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
          amimn   = amn(nn) / ( ami(ni) + amn(nn) )
          nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
          nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     nitrogen (N)

      ni = ptnp
      do nn = 1,nneut
         do i = 1,nz
          amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
          amimn   = amn(nn) / ( ami(ni) + amn(nn) )
          nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
          nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     oxygen (O)

      ni = ptop
      do nn = 1,nneut
         do i = 1,nz
          if ( nn .eq. pto ) then
            teff    = 0.5 * ( ti(i,nfl,ni) + tn(i,nfl) )
            fac     = ( 1.04 - .067 * alog10(teff) ) ** 2
            tfactor = sqrt(teff) * fac
            nuin(i,ni,nn)  = 4.45e-11 * denn(i,nfl,nn) * tfactor
          else
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          endif
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     nitrogen 2(N2)

      ni = ptn2p
      do nn = 1,nneut
         do i = 1,nz
          if ( nn .eq. ptn2 ) then
            teff    = 0.5 * ( ti(i,nfl,ni) + tn(i,nfl) )
!       MS: According to the SAMI2 paper the first coefficient in
!       fac should be 1.04, not 1.00.
!            fac     = ( 1.00 - .069 * alog10(teff) ) ** 2 ! original
            fac     = ( 1.04 - .069 * alog10(teff) ) ** 2
            tfactor = sqrt(teff) * fac
            nuin(i,ni,nn) = 5.14e-11 * denn(i,nfl,nn) * tfactor
          else
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          endif
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     nitrous oxide (N0)

      ni = ptnop
      do nn = 1,nneut
         do i = 1,nz
          amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
          amimn   = amn(nn) / ( ami(ni) + amn(nn) )
          nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
          nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     oxygen 2(O2)

      ni = pto2p
      do nn = 1,nneut
         do i = 1,nz
          if ( nn .eq. pto2 ) then
            teff    = 0.5 * ( ti(i,nfl,ni) + tn(i,nfl) )
            fac     = ( 1.00 - .073 * alog10(teff) ) ** 2
            tfactor = sqrt(teff) * fac
            nuin(i,ni,nn) = 2.59e-11 * denn(i,nfl,nn) * tfactor
          else
            amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
            amimn   = amn(nn) / ( ami(ni) + amn(nn) )
            nufacin = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
            nuin(i,ni,nn) = nufacin * denn(i,nfl,nn)
          endif
          nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
      enddo

!     ion-ion collision frequency

      nuij = 0.

      do nj = nion1,nion2
        do ni = nion1,nion2
          if ( ni .ne. nj ) then
             do i = 1,nz
              alame1  = ( ami(ni) + ami(nj) ) * evtok /&
                     ( ami(ni)*ti(i,nfl,nj) + ami(nj)*ti(i,nfl,ni) )
              alame2  = deni(i,nfl,ni) * evtok / ti(i,nfl,ni) +&
                       deni(i,nfl,nj) * evtok / ti(i,nfl,nj)
              if ( alame2 .lt. 0 ) then
                print *,'ni,i,nj,nfl,tii,tij,alame1,alame2,nii,nij',&
                         ni,i,nj,nfl,ti(i,nfl,ni),ti(i,nfl,nj),&
                         alame1,alame2,&
                         deni(i,nfl,ni),deni(i,nfl,nj)
                print *,' code has gone bad ...'
                stop
              endif
              alame   = alame1 * sqrt(alame2)
              alam    = 23. - alog(alame)
              amufac  = (ami(nj)/ami(ni))/(ami(ni) +ami(nj))
              nufacij = 9.2e-2*alam*sqrt(amufac)
              nuij(i,ni,nj) =  nufacij * deni(i,nfl,nj)&
                              / sqrt( ti(i,nfl,ni)**3 )
          enddo
          endif
        enddo
      enddo

 100  format(1x,2e12.2)

!     sumnuj: sum of ion-ion coll freq and nuin

      do ni = nion1,nion2
        do i = 1,nz
          sumnuj(i,ni) = 0.
          do nj = nion1,nion2
            sumnuj(i,ni) = sumnuj(i,ni) + nuij(i,ni,nj)
          enddo
          sumnuj(i,ni) = sumnuj(i,ni) + nuint(i,ni)
        enddo
      enddo

!     update ne

      ne = 1.

      do ni = nion1,nion2
        do i = 1,nz
          ne(i,nfl) = ne(i,nfl) + deni(i,nfl,ni)
        enddo
      enddo

!     get a new value for vsid

      do ni = nion1,nion2
        do i = 2,nz-1
          mi    = amu * ami(ni)
          k0    = bolt / mi
          term1 = nuint(i,ni) * tvn(i) + sumvsi(i,nfl,ni) + gs(i)
          pip   = 0.5 * (   deni(i+1,nfl,ni) * ti(i+1,nfl,ni)&
                         + deni(i,nfl,ni)   * ti(i,nfl,ni)   )
          pim   = 0.5 * (   deni(i,nfl,ni)   * ti(i,nfl,ni)&
                         + deni(i-1,nfl,ni) * ti(i-1,nfl,ni) )
          denid =&
                  (        deni(i-1,nfl,ni)&
                    + 4. * deni(i,nfl,ni)&
                    +      deni(i+1,nfl,ni)  ) / 6.
          term2 =  - bms(i,nfl) * k0 /  denid&
                  * ( pip - pim ) / d22s(i,nfl)
          pep   = 0.5 * (   ne(i+1,nfl) * te(i+1,nfl)&
                         + ne(i,nfl)   * te(i,nfl)   )
          pem   = 0.5 * (   ne(i,nfl)   * te(i,nfl)&
                         + ne(i-1,nfl) * te(i-1,nfl) )
          dened =&
       ( ne(i-1,nfl) + 4. * ne(i,nfl) + ne(i+1,nfl) ) / 6.
          term3 =  - bms(i,nfl) * k0 /  dened&
                  * ( pep - pem ) / d22s(i,nfl)

          vsid(i,nfl,ni)  =  term1 + term2 + term3

          if ( deni(i,nfl,ni) .le. .0001*ne(i,nfl) )&
           vsid(i,nfl,ni) =   vsid(i,nfl,ni)&
                             * exp ( -.0001*ne(i,nfl)/deni(i,nfl,ni) )

        enddo
      enddo

!     fix up end points for vsid

      do ni = nion1,nion2
        vsid (1,nfl,ni)    = vsid (2,nfl,ni)
        vsid (nz,nfl,ni)   = vsid (nz-1,nfl,ni)
      enddo

!     calculate collisional ion velocity
!     not used; simply a diagnostic

!      do i = 1,nz
!        do ni = nion1,nion2
!          vsic(i,nfl,ni) = vsid(i,nfl,ni) / sumnuj(i,ni)
!        enddo
!      enddo

      return
      end


!*******************************************
!*******************************************

!            rtryds

!*******************************************
!*******************************************

      subroutine rtryds(a,b,c,d,x,n)

      include 'param-1.00.inc'
      include 'com-1.00.inc'

!     arrays a,b,c, and d may be used for storage of alfa, beta and x
!     in the actual call of this routine, but remember, whatever you
!     use will be lost by the definition of of alfa and beta here.
!     form,  a(k)*x(k-1) + b(k)*x(k) + c(k)*x(k+1) = d(k)

!     i have modified the input sequence to the routine, but have left it
!     otherwise intact.  we may  want to eventually change this (gj)
      implicit none

      real,DIMENSION(:):: a,b,c,d,x
      INTEGER::n
      real:: alfa(nz),beta(nz)

      !LOCAL VARIABLES
      INTEGER::nm1,k,i
      REAL::dst,rb,ast,z
      nm1=n-1

!     apply the boundary condition at x(1)
!     alfa(1) and beta(1) determined from b(1),c(1),d(1),a(1)=0.

      dst     = d(1)
      rb      = 1. / b(1)
      alfa(1) = -c(1) * rb
      beta(1) =   dst * rb

!     calculate the alfas and betas of k on forward sweep

      do k=2,nm1
        ast     =  a(k)
        z       =  1. / ( b(k) + ast * alfa(k-1) )
        dst     =  d(k)
        alfa(k) = -c(k) * z
        beta(k) =  ( dst - ast * beta(k-1) ) * z
      enddo

!     apply the boundary condition at x(n)
!     x(n) determined from a(n),b(n),d(n),c(n)=0.

      x(n) = ( d(n) - a(n) *beta(nm1) ) / ( b(n) + a(n) * alfa(nm1) )

!     calculate x of k from the alfas and betas on backward sweep

      do i=2,n
        k    = n + 1 - i
        x(k) = x(k+1) * alfa(k) + beta(k)
      enddo

      return
      end


!*******************************************
!*******************************************

!            vsisolv

!*******************************************
!*******************************************

      subroutine vsisolv ( vi,vid,viold,snuj,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none

      REAL,DIMENSION(nz)::vi
      REAL,DIMENSION(nz)::vid
      REAL,DIMENSION(nz)::viold
      REAL,DIMENSION(nz)::snuj
      INTEGER::nfl

      !LOCAL VARIABLES
      REAL,DIMENSION(nz):: a, b, c, d
      INTEGER::j
      REAL::ujm1,uj,ujp1,ur,ul,a0,b0,c0
      

!     initialize

      a = 0.
      b = 0.
      c = 0.
      d = 0.

      do j = 2,nz-1

        ujm1 = vi(j-1)
        uj   = vi(j)
        ujp1 = vi(j+1)
        ur = .25*( uj +ujp1)
        ul = .25*( uj +ujm1)

        if (ur .ge. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 =  ur
          c0 =  0.
        endif
        if (ur .le. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = -ul
          c0 = ur
        endif
        if (ur .ge. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = ur - ul
          c0 = 0.
        endif
        if (ur .le. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 = 0.
          c0 = ur
        endif

        a(j) = a0 / d22s(j,nfl) * bms(j,nfl)

        b(j) = 1/dt + snuj(j) + b0 / d22s(j,nfl) * bms(j,nfl)

        c(j) = c0 / d22s(j,nfl) * bms(j,nfl)

        d(j) = viold(j)/dt + vid(j)

      enddo

!     we will assume that the bc's are the neutral temperature
!     at both ends of the field line

!     lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) = 0.

!     upper bc

      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) = 0.

      call rtryds(a,b,c,d,vi,nz)

      return
      end


! *********************
!
!     smoothz
!
! *********************

      subroutine smoothz(finout,ncomp)

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none


      REAL,DIMENSION(nz)::finout
      INTEGER::ncomp

      !LOCAL VARIABLES
      REAL:: const
      INTEGER::i,ip1,im1
      REAL,DIMENSION(nz)::tempz
! 
! This is the binomial filter (in x space) as described in  
! Birdsall appendix C. 
! We have the choice of a compensating filter or not. 
! if ncomp=0, no compensation, else compensation 
! 
 
! do smoothz in the z direction 

       do i = 1,nz
          ip1 = i +1
          if(i .eq. nz) ip1 = 1
          im1 = i -1
          if(i .eq. 1) im1 = nz
          tempz(i) = .25*(finout(im1) +2.*finout(i)&
                         +finout(ip1))
       enddo
       do i = 1,nz
          finout(i) = tempz(i)
       enddo

       if ( ncomp .ne. 0 ) then

! put in compensator
! the alogrithm below is equivalent to
! fftmp(i)=(1./16.)*(-ff0(i-2)+4.*ff0(i-1)+10.*ff0(i)+4.*ff0(i+1)-ff0(i+2))

! do compensation in the z direction

       const = sqrt(1.4571072)
       do i = 1,nz
          ip1 = i +1
          if(i .eq. nz) ip1 = 1
          finout(i) = const*(finout(i) -.171573*finout(ip1))
       enddo
       do i = nz,1,-1
          im1 = i -1
          if(i .eq. 1) im1 = nz
          finout(i) = const*(finout(i) -.171573*finout(im1))
       enddo

      endif

      return
      end

!*******************************************
!*******************************************

!            densolv2

!*******************************************
!*******************************************

      subroutine densolv2( ni,tdeni,prod,loss,oldion,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none

      INTEGER::ni
      REAL,DIMENSION(nz)::tdeni
      REAL,DIMENSION(nz)::prod
      REAL,DIMENSION(nz)::loss
      REAL,DIMENSION(nz)::oldion
      INTEGER::nfl

      !LOCAL VARIABLES
      REAL,DIMENSION(nz)::a
      REAL,DIMENSION(nz)::b
      REAL,DIMENSION(nz)::c
      REAL,DIMENSION(nz)::d
      INTEGER::j
      REAL::ujm1,uj,ujp1,ur,ul,a0,b0,c0

!     initialize

      a = 0.
      b = 0.
      c = 0.
      d = 0.

      do j = 2,nz-1

      ujm1  = vsi(j-1,nfl,ni)/bms(j-1,nfl)
      uj    = vsi(j,nfl,ni)  /bms(j,nfl)
      ujp1  = vsi(j+1,nfl,ni)/bms(j+1,nfl)
      ur = .5*( uj +ujp1)
      ul = .5*( uj +ujm1)

      if (ur .ge. 0. .and. ul .ge. 0.) then
        a0 = -ul
        b0 =  ur
        c0 =  0.
      endif
      if (ur .le. 0. .and. ul .le. 0.) then
        a0 = 0.
        b0 = -ul
        c0 = ur
      endif
      if (ur .ge. 0. .and. ul .le. 0.) then
        a0 = 0.
        b0 = ur - ul
        c0 = 0.
      endif
      if (ur .le. 0. .and. ul .ge. 0.) then
        a0 = -ul
        b0 = 0.
        c0 = ur
      endif

      a(j) =  a0 * bms(j,nfl) ** 2 / d22s(j,nfl)

      b(j) = 1. / dt + loss(j) + b0 * bms(j,nfl) ** 2 / d22s(j,nfl)

      c(j) = c0 * bms(j,nfl) ** 2 / d22s(j,nfl)

      d(j) = oldion(j) / dt + prod(j)

      enddo

!     we will assume that they are determined by the production and loss
!     at both ends of the field line

!     lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) =&
  sqrt ( tdeni(1) * prod(1) / loss(1) ) + denmin

!     upper bc

      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) =&
  sqrt ( tdeni(nz) * prod(nz) / loss(nz) ) + denmin


      call rtryds ( a,b,c,d,tdeni,nz )

      return
      end


!*******************************************
!*******************************************

!            etemp

!*******************************************
!*******************************************

      subroutine etemp ( tte,te_old,phprodr,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none

      REAL,DIMENSION(nz)::tte
      REAL,DIMENSION(nz)::te_old
      REAL,DIMENSION(nz,nion)::phprodr
      INTEGER::nfl

      !LOPCAL VARIABLES
      real kape(nz)
      real s1e(nz),s2e(nz),s3e(nz),s4e(nz)
      real s5e(nz),qphe(nz),phprod(nz)
      real qen(nz,nneut)
      real ratio(nz)
      real divvexb(nz)
      integer iz300s(nf),iz300n(nf)

      INTEGER::i,ni,nn,izs,izn,nzh
      REAL::fac1,fac2,fac3,akpefac,xs3e,xs4e
      REAL::xarg,x,earg,epsi,facts
      REAL::ne300s,o2300,n2300,o300,phprod300
      REAL::xarg300,x300,earg300,epsi300
      REAL::q0s,factn,ne300n,q0n,xbms,xbmn
      REAL::dels300s,dels300n,xn,xs,xints,xintn,xqs
      REAL::xqn,vexbeq

      s1e = 0.
      s2e = 0.
      s3e = 0.
      s4e = 0.
      kape = 0.
      qen = 0.

      do i = 1,nz

        fac1 = denn(i,nfl,pto)  * 1.1e-16&
                * ( 1. + 5.7e-4 * te(i,nfl) )
        fac2 = denn(i,nfl,ptn2) * 2.82e-17 * sqrt(te(i,nfl))&
                * ( 1  - 1.2e-4 * te(i,nfl) )
        fac3 = denn(i,nfl,pto2) * 2.2e-16&
               * ( 1. + 3.6e-2  * sqrt(te(i,nfl)) )
        akpefac = fac1 + fac2 + fac3

        kape(i) = 7.7e5 * sqrt ( te(i,nfl)**5 ) * 0.6667 * evtok&
           / ( 1. + 3.22e4 * ( te(i,nfl)**2 / ne(i,nfl) * akpefac) )


!       neutrals (Tn - Te) term

!       N2

!       vibrational state from red book (p. 269) milward et al.
!       removed (2/16/01)

        qen(i,ptn2) = .6667 *  evtok * denn(i,nfl,ptn2) *&
                       ( 1.2e-19 * ( 1. - 1.2e-4 * te(i,nfl) )&
                                 * te(i,nfl) +&
                         2.e-14 / sqrt(te(i,nfl))&
                         + 6.5e-22 * ( tn(i,nfl) - 310 ) ** 2 *&
                           exp(.0023*(te(i,nfl) - tn(i,nfl))) )

!       O2

        qen(i,pto2) = .6667 * evtok * denn(i,nfl,pto2) *&
                      ( 7.9e-19 * ( 1. + 3.6e-2 * sqrt(te(i,nfl)) ) *&
                          sqrt(te(i,nfl)) +&
                        7.e-14 / sqrt(te(i,nfl)) )

!       O

        qen(i,pto) = .6667 * 7.2e-18 * evtok * denn(i,nfl,pto) *&
                       sqrt(te(i,nfl))

!       H

        qen(i,pth) = .6667 * 6.3e-16 * evtok * denn(i,nfl,pth) *&
                       ( 1. - 1.35e-4 * te(i,nfl) ) *&
                       sqrt(te(i,nfl))

        do nn = 1,nneut
          s2e(i) = s2e(i) + qen(i,nn)
        enddo

        s1e(i) = s2e(i) * tn(i,nfl)

!       ions (Ti - Te) term

        do ni = nion1,nion2
          xs3e    = 7.7e-6 * deni(i,nfl,ni) / ami(ni)&
                          / te(i,nfl) / sqrt(te(i,nfl))&
                          * .66667 * evtok
          xs4e    = xs3e * ti(i,nfl,ni)
          s3e(i) = s3e(i) + xs3e
          s4e(i) = s4e(i) + xs4e
        enddo

      enddo

!     photoelectron heating
!     red book (millward et al. p. 269)

!     calculate total ion photoproduction (= photoelectron)

      do i = 1,nz
        phprod(i)   = 0.
        do ni = nion1,nion2
          phprod(i) = phprodr(i,ni) * denn(i,nfl,ni) + phprod(i)
        enddo
      enddo

! comment out: use iz300s/n calculated in grid

! iz300s/iz300n are redefined here

      do i = 1,nz
        ratio(i) = ne(i,nfl) / &
                   (0.1*denn(i,nfl,pto)+&
                    denn(i,nfl,pto2)+denn(i,nfl,ptn2)) 
      enddo
      
      iz300s=1
      i = 1 
      do while ( ratio(i) .le. 3.e-3 .and. i .lt. nz ) 
         iz300s(nfl) = i 
         i         = i + 1 
      enddo 
 
      iz300n=nz
      i = nz 
      do while ( ratio(i) .le. 3.e-3 .and. i .gt. 1 )  
         iz300n(nfl) = i 
         i         = i - 1 
      enddo 

      if ( iz300s(nfl) .gt. iz300n(nfl) ) then

        do i = 1,nz
            xarg =   ne(i,nfl)&
             / (        denn(i,nfl,pto2)&
                 +      denn(i,nfl,ptn2)&
                 + .1 * denn(i,nfl,pto)   )
            x    = alog ( xarg )
            earg =     12.75&
               + 6.941 * x&
               + 1.166 * x ** 2&
               + 0.08034 * x ** 3&
               + 0.001996 * x ** 4
            epsi = exp ( -earg )
            qphe(i) = epsi * phprod(i)
          enddo
      else
          do i = 1,iz300s(nfl)
            xarg =   ne(i,nfl)&
                   / (        denn(i,nfl,pto2)&
                       +      denn(i,nfl,ptn2)&
                       + .1 * denn(i,nfl,pto)   )
            x    = alog ( xarg )
            earg =     12.75&
                     + 6.941 * x&
                     + 1.166 * x ** 2&
                     + 0.08034 * x ** 3&
                     + 0.001996 * x ** 4
            epsi = exp ( -earg )
            qphe(i) = epsi * phprod(i)
          enddo

!       smooth things at 300 km 

        izs       = iz300s(nfl)
!        facts = (250.-alts(izs,nfl)) / 
!     .             (alts(izs+1,nfl)-alts(izs,nfl))
        facts = (3.e-3-ratio(izs)) / &
                (ratio(izs+1)-ratio(izs))
        ne300s = ne(izs,nfl) + (ne(izs+1,nfl)-ne(izs,nfl)) * facts
        o2300 = denn(izs,nfl,pto2) + &
               (denn(izs+1,nfl,pto2)-denn(izs,nfl,pto2)) * facts
        n2300 = denn(izs,nfl,ptn2) + &
               (denn(izs+1,nfl,ptn2)-denn(izs,nfl,ptn2)) * facts
        o300 = denn(izs,nfl,pto) + &
               (denn(izs+1,nfl,pto)-denn(izs,nfl,pto)) * facts
        phprod300 = phprod(izs) + &
              (phprod(izs+1)-phprod(izs)) * facts
        xarg300 = ne300s / ( o2300 + n2300 + 0.1*o300 )
        x300 = alog( xarg300)
        earg300 =     12.75 +&
             6.941 * x300 +&
             1.166 * x300 ** 2 +&
             0.08034 * x300 ** 3 +&
             0.001996 * x300 ** 4
        epsi300 = exp ( -earg300 )
        q0s = epsi300 * phprod300 / ne300s

        do i = iz300n(nfl),nz
          xarg =   ne(i,nfl)&
                 / (       denn(i,nfl,pto2)&
                    +      denn(i,nfl,ptn2)&
                    + .1 * denn(i,nfl,pto) )
          x    = alog ( xarg )
          earg =     12.75&
                  + 6.941 * x&
                  + 1.166 * x ** 2&
                  + 0.08034 * x ** 3&
                  + 0.001996 * x ** 4
          epsi = exp ( -earg )
          qphe(i) = epsi * phprod(i)
        enddo

        izn      = iz300n(nfl)
!        factn = (250.-alts(izn,nfl)) / 
!     .             (alts(izn-1,nfl)-alts(izn,nfl))
        factn = (3.e-3-ratio(izn)) / &
                 (ratio(izn-1)-ratio(izn))
        ne300n = ne(izn,nfl) + &
              (ne(izn-1,nfl)-ne(izn,nfl)) * factn
        o2300 = denn(izn,nfl,pto2) + &
              (denn(izn-1,nfl,pto2)-denn(izn,nfl,pto2)) * factn
        n2300 = denn(izn,nfl,ptn2) + &
              (denn(izn-1,nfl,ptn2)-denn(izn,nfl,ptn2)) * factn
        o300 = denn(izn,nfl,pto) + &
              (denn(izn-1,nfl,pto)-denn(izn,nfl,pto)) * factn
        phprod300 = phprod(izn) + &
              (phprod(izn-1)-phprod(izn)) * factn
        xarg300 = ne300n / ( o2300 + n2300 + 0.1*o300 )
        x300 = alog( xarg300)
        earg300 =     12.75 +&
             6.941 * x300 +&
             1.166 * x300 ** 2 +&
             0.08034 * x300 ** 3 +&
             0.001996 * x300 ** 4
        epsi300 = exp ( -earg300 )
        q0n = epsi300 * phprod300 / ne300n

        xbms = bms(izs,nfl) + (bms(izs+1,nfl)-bms(izs,nfl)) * facts
        xbmn = bms(izn,nfl) + (bms(izn-1,nfl)-bms(izn,nfl)) * factn


        dels300s = dels(iz300s(nfl),nfl) * facts
        dels300n = dels(iz300n(nfl)-1,nfl) * factn

        ! MS: Old code used a wasteful way to calculate xn. 
        ! Cleaner version here. 
        xn = 0. 
        ! Set bottom integration bound to 300 km. 
        xn =   xn + 0.5 * ( ne(iz300n(nfl)-1,nfl) + ne300n ) * &
              (dels(iz300n(nfl)-1,nfl) - dels300n ) 
        do i =iz300n(nfl)-2,iz300s(nfl)+1,-1 
           xn = xn + 0.5 * ( ne(i,nfl) + ne(i+1,nfl) ) * dels(i,nfl) 
        enddo

!        cqe   = 6.e-14                     ! constant (now in namelist)

        if ( q0s .lt. 0 .or. q0n .lt. 0 ) then
          print *,' q0s = ',q0s,' q0n = ',q0n,' nfl = ',nfl
        endif

!       1/22/00

!       put in dels (arc length along field line)

        xs    = 0.

        do i = iz300s(nfl)+1,iz300n(nfl)-1 
           if (i .eq. iz300s(nfl)+1) then 
             xs = xs + 0.5*( ne300s + ne(i,nfl) ) * &
                 (dels(iz300s(nfl),nfl) - dels300s) 
           else 
              xs = xs + 0.5 * ( ne(i,nfl) + ne(i-1,nfl) ) &
                               * dels(i-1,nfl) 
              xn = xn - 0.5 * ( ne(i,nfl) + ne(i-1,nfl) ) &
                               * dels(i-1,nfl) 
           endif 
 
           xints = cqe*xs
           xintn = cqe*xn
           xqs    = ne(i,nfl) * q0s * bms(i,nfl) / xbms * exp(-xints) 
           xqn    = ne(i,nfl) * q0n * bms(i,nfl) / xbmn * exp(-xintn) 
           qphe(i) = xqs + xqn 
        enddo 

      endif

      do i = 1,nz
        s5e(i) = 0.66667 * evtok * qphe(i) / ne(i,nfl) ! * .15
      enddo 

! MS: Neglected term, divergence of ExB drift
! Divergence of the ExB drift; requires equatorial drift

      nzh    = (nz+1)/2
      vexbeq = vexb(nzh,nfl)
      do i = 1,nz
        divvexb(i) = 6.*vexbeq /&
                       (ps(i,nfl)*re*1.e5) *&
                       cos(blats(i,nfl)*po180)**2 *&
                       (1.+sin(blats(i,nfl)*po180)**2) /&
                       (1.+3.*sin(blats(i,nfl)*po180)**2)**2
        s2e(i) = s2e(i) - 0.6667 * divvexb(i)
      enddo

      call tesolv(tte,te_old,kape,s1e,s2e,s3e,s4e,s5e,nfl)

      return
      end


!*******************************************
!*******************************************

!            htemp

!*******************************************
!*******************************************

      subroutine htemp ( tti,tiold,tvn,nuin,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none
      REAL,DIMENSION(nz)::tti
      REAL,DIMENSION(nz)::tiold
      REAL,DIMENSION(nz)::tvn
      REAL,DIMENSION(nz,nion,nneut)::nuin
      INTEGER::nfl

      !LOCAL VARIABLES

      real kapi(nz),s1i(nz),s2i(nz),s3i(nz),s4i(nz),s5i(nz)
      real s6i(nz),s7i(nz)
      real divvexb(nz)
      REAL::convfac,redmass,tfac,xs6i,xs7i,vexbeq
      INTEGER::i,nn,ni,nzh

      convfac = amu / bolt / 3.

      s1i = 0.
      s2i = 0.
      s3i = 0.
      s4i = 0.
      s5i = 0.
      s6i = 0.
      s7i = 0.
      kapi = 0.

      do i = 1,nz

        kapi(i) = 4.6e+4 * sqrt ( ti(i,nfl,pthp)**5 ) / ne(i,nfl) *&
                 deni(i,nfl,pthp) / sqrt(ami(pthp))

        kapi(i)  = 0.6667 * kapi(i) * evtok

!       neutrals

        do nn = 1,nneut
          redmass =&
          ami(pthp) * amn(nn) / ( ami(pthp) + amn(nn) ) ** 2
          s2i(i) = s2i(i) + 2. * nuin(i,pthp,nn) * redmass
          s3i(i) = s3i(i)&
           + convfac * amn(nn)&
                     * abs ( vsi(i,nfl,pthp) - tvn(i) ) ** 2&
           * 2. * nuin(i,pthp,nn) * redmass
        enddo

        s1i(i) = s2i(i) * tn(i,nfl)

!       electrons

        s4i(i) = 7.7e-6 * ne(i,nfl) / ami(pthp)&
                        / te(i,nfl) / sqrt(te(i,nfl))&
                        * .66667 * evtok
        s5i(i) = s4i(i) * te(i,nfl)

!       other ions

        do ni = nion1,nion2
!          if ( ni .ne. ptop ) then
          if ( ni .ne. pthp ) then
            tfac    =    ti(i,nfl,pthp) / ami(pthp)&
                     +  ti(i,nfl,ni) / ami(ni)
            xs6i    = 3.3e-4 * deni(i,nfl,ni) / ami(pthp) / ami(ni)&
                     / tfac / sqrt(tfac) * .66667 * evtok
            xs7i    = xs6i * ti(i,nfl,ni)
            s6i(i) = s6i(i) + xs6i
            s7i(i) = s7i(i) + xs7i
          endif
        enddo

      enddo

! MS: Neglected term, divergence of ExB drift
! Divergence of the ExB drift; requires equatorial drift

      nzh = (nz+1)/2
      vexbeq = vexb(nzh,nfl)
      do i = 1,nz
        divvexb(i) = 6.*vexbeq /&
                    (ps(i,nfl)*re*1.e5) *&
                    cos(blats(i,nfl)*po180)**2 *&
                    (1.+sin(blats(i,nfl)*po180)**2) /&
                    (1.+3.*sin(blats(i,nfl)*po180)**2)**2
        s2i(i) = s2i(i) - 0.6667 * divvexb(i)
      enddo

      call tisolv(tti,tiold,kapi,s1i,s2i,s3i,s4i,s5i,s6i,s7i,pthp,nfl)

      return
      end


!*******************************************
!*******************************************

!            hetemp

!*******************************************
!*******************************************

      subroutine hetemp ( tti,tiold,tvn,nuin,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none
      REAL,DIMENSION(nz)::tti
      REAL,DIMENSION(nz)::tiold
      REAL,DIMENSION(nz)::tvn
      REAL,DIMENSION(nz,nion,nneut)::nuin
      INTEGER::nfl

      !LOCAL VARIABLES
      real kapi(nz),s1i(nz),s2i(nz),s3i(nz),s4i(nz),s5i(nz)
      real s6i(nz),s7i(nz)
      real divvexb(nz)
      REAL::convfac,redmass,tfac,xs6i,xs7i,vexbeq
      INTEGER::i,nn,ni,nzh


      convfac = amu / bolt / 3.

      s1i = 0.
      s2i = 0.
      s3i = 0.
      s4i = 0.
      s5i = 0.
      s6i = 0.
      s7i = 0.
      kapi = 0.

      do i = 1,nz

        kapi(i) = 4.6e+4 * sqrt ( ti(i,nfl,pthep)**5 ) / ne(i,nfl) *&
                 deni(i,nfl,pthep) / sqrt(ami(pthep))
        kapi(i)  = 0.6667 * kapi(i) * evtok

!       neutrals

        do nn = 1,nneut
          redmass =&
          ami(pthep) * amn(nn) / ( ami(pthep) + amn(nn) ) ** 2
          s2i(i) = s2i(i) + 2. * nuin(i,pthep,nn) * redmass
          s3i(i) = s3i(i)&
           + convfac * amn(nn)&
                     * abs ( vsi(i,nfl,pthep) - tvn(i) ) ** 2&
           * 2. * nuin(i,pthep,nn) * redmass
        enddo

        s1i(i) = s2i(i) * tn(i,nfl)

!       electrons

        s4i(i) = 7.7e-6 * ne(i,nfl) / ami(pthep)&
                        / te(i,nfl) / sqrt(te(i,nfl))&
                        * .66667 * evtok
        s5i(i) = s4i(i) * te(i,nfl)

!       other ions

        do ni = nion1,nion2
!          if ( ni .ne. ptop ) then
          if ( ni .ne. pthep ) then
            tfac    =   ti(i,nfl,pthep) / ami(pthep)&
                     + ti(i,nfl,ni) / ami(ni)
            xs6i    = 3.3e-4 * deni(i,nfl,ni) / ami(pthep) / ami(ni)&
                     / tfac / sqrt(tfac) * .66667 * evtok
            xs7i    = xs6i * ti(i,nfl,ni)
            s6i(i) = s6i(i) + xs6i
            s7i(i) = s7i(i) + xs7i
          endif
        enddo

      enddo

! MS: Neglected term, divergence of ExB drift
! Divergence of the ExB drift; requires equatorial drift

      nzh = (nz+1)/2
      vexbeq = vexb(nzh,nfl)
      do i = 1,nz 
        divvexb(i) = 6.*vexbeq / &
                     (ps(i,nfl)*re*1.e5) * &
                     cos(blats(i,nfl)*po180)**2 * &
                     (1.+sin(blats(i,nfl)*po180)**2) / &
                     (1.+3.*sin(blats(i,nfl)*po180)**2)**2
        s2i(i) = s2i(i) - 0.6667 * divvexb(i)
      enddo

      call tisolv(tti,tiold,kapi,s1i,s2i,s3i,s4i,s5i,s6i,s7i,pthep,nfl)

      return
      end


!*******************************************
!*******************************************

!            otemp

!*******************************************
!*******************************************

      subroutine otemp ( tti,tiold,tvn,nuin,nfl )

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none
      REAL,DIMENSION(nz)::tti
      REAL,DIMENSION(nz)::tiold
      REAL,DIMENSION(nz)::tvn
      REAL,DIMENSION(nz,nion,nneut)::nuin
      INTEGER::nfl

      real kapi(nz),s1i(nz),s2i(nz),s3i(nz),s4i(nz),s5i(nz)
      real s6i(nz),s7i(nz)
      real divvexb(nz)
      REAL::convfac,redmass,tfac,xs6i,xs7i,vexbeq
      INTEGER::i,nn,ni,nzh


      convfac = amu / bolt / 3.

      s1i = 0.
      s2i = 0.
      s3i = 0.
      s4i = 0.
      s5i = 0.
      s6i = 0.
      s7i = 0.
      kapi = 0.

      do i = 1,nz

        kapi(i) = 4.6e+4 * sqrt ( ti(i,nfl,ptop)**5 ) / ne(i,nfl) *&
                 deni(i,nfl,ptop) / sqrt(ami(ptop))
        kapi(i)  = 0.6667 * kapi(i) * evtok

!       neutrals

        do nn = 1,nneut
          redmass =&
          ami(ptop) * amn(nn) / ( ami(ptop) + amn(nn) ) ** 2
          s2i(i) = s2i(i) + 2. * nuin(i,ptop,nn) * redmass
          s3i(i) = s3i(i)&
           + convfac * amn(nn)&
                     * abs ( vsi(i,nfl,ptop) - tvn(i) ) ** 2&
           * 2. * nuin(i,ptop,nn) * redmass
        enddo

        s1i(i) = s2i(i) * tn(i,nfl)

!       electrons

        s4i(i) = 7.7e-6 * ne(i,nfl) / ami(ptop)&
                        / te(i,nfl) / sqrt(te(i,nfl))&
                        * .66667 * evtok
        s5i(i) = s4i(i) * te(i,nfl)

!       other ions

        do ni = nion1,nion2
          if ( ni .ne. ptop ) then
            tfac    =    ti(i,nfl,ptop) / ami(ptop)&
                      + ti(i,nfl,ni) / ami(ni)
            xs6i    = 3.3e-4 * deni(i,nfl,ni) / ami(ptop) / ami(ni)&
                     / tfac / sqrt(tfac) * .66667 * evtok
            xs7i    = xs6i * ti(i,nfl,ni)
            s6i(i) = s6i(i) + xs6i
            s7i(i) = s7i(i) + xs7i
          endif
        enddo

      enddo

! MS: Neglected term, divergence of ExB drift
! Divergence of the ExB drift; requires equatorial drift

      nzh = (nz+1)/2
      vexbeq = vexb(nzh,nfl)
      do i = 1,nz 
        divvexb(i) = 6.*vexbeq / &
                     (ps(i,nfl)*re*1.e5) * &
                     cos(blats(i,nfl)*po180)**2 * &
                     (1.+sin(blats(i,nfl)*po180)**2) / &
                     (1.+3.*sin(blats(i,nfl)*po180)**2)**2
        s2i(i) = s2i(i) - 0.6667 * divvexb(i)
      enddo

      call tisolv(tti,tiold,kapi,s1i,s2i,s3i,s4i,s5i,s6i,s7i,ptop,nfl)

      return
      end



!*******************************************
!*******************************************

!            tisolv

!*******************************************
!*******************************************

      subroutine tisolv(tti,tio,kap,s1,s2,s3,s4,s5,s6,s7,npt,nfl)

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none
      REAL,dimension(nz)::tti
      REAL,dimension(nz)::tio
      REAL,dimension(nz)::kap
      REAL,dimension(nz)::s1,s2,s3,s4,s5,s6,s7
      INTEGER::npt,nfl

      !LOCAL VARIABLES
      real a(nz),b(nz),c(nz),d(nz)
      INTEGER::j
      REAL::ujm1,uj,ujp1,ur,ul,a0,b0,c0
      

!     initialize

      a = 0.
      b = 0.
      c = 0.
      d = 0.

      do j = 2,nz-1
        ujm1 = bms(j-1,nfl)*vsi(j-1,nfl,npt)
        uj   = bms(j,nfl)  *vsi(j,nfl,npt)
        ujp1 = bms(j+1,nfl)*vsi(j+1,nfl,npt)
        ur = .5*( uj +ujp1)
        ul = .5*( uj +ujm1)

        if (ur .ge. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 =  ur
          c0 =  0.
        endif
        if (ur .le. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = -ul
          c0 = ur
        endif
        if (ur .ge. 0. .and. ul .le. 0.) then
          a0 = 0.
          b0 = ur - ul
          c0 = 0.
        endif
        if (ur .le. 0. .and. ul .ge. 0.) then
          a0 = -ul
          b0 = 0.
          c0 = ur
        endif

        a(j) =     a0 / d22s(j,nfl)&
         - ( bms(j,nfl)**2 / deni(j,nfl,npt) ) / d22s(j,nfl)&
           *.5 * ( kap(j) + kap(j-1) ) / ds(j,nfl)

        b(j) = 1. / dt + b0 / d22s(j,nfl)&
         -.333333 * ( bms(j,nfl)&
                     * (vsi(j+1,nfl,npt) - vsi(j-1,nfl,npt) )&
                     + 5. * vsi(j,nfl,npt)&
                          * (bms(j+1,nfl) - bms(j-1,nfl) ) )&
         / d2s(j,nfl)&
         +  ( bms(j,nfl)**2 / deni(j,nfl,npt) ) / d22s(j,nfl)&
           *(.5* (kap(j+1) + kap(j) ) / ds(j+1,nfl)&
         +.5 * (kap(j) + kap(j-1) ) / ds(j,nfl))&
         + s2(j) + s4(j) + s6(j)

        c(j) =     c0 / d22s(j,nfl)&
         - ( bms(j,nfl)**2 / deni(j,nfl,npt) ) /d22s(j,nfl)&
           *.5 * (kap(j+1) + kap(j) ) / ds(j+1,nfl)

        d(j) = tio(j)/dt + s1(j) + s3(j) + s5(j) + s7(j)

      enddo

!     we will assume that the bc's are the neutral temperature
!     at both ends of the field line

!     lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) = tn(1,nfl)

!     upper bc

      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) = tn(nz,nfl)

      call rtryds ( a,b,c,d,tti,nz )


      return
      end

!*******************************************
!*******************************************

!            tesolv

!*******************************************
!*******************************************

      subroutine tesolv(tte,te_old,kap,s1,s2,s3,s4,s5,nfl)

      include 'param-1.00.inc'
      include 'com-1.00.inc'
      implicit none
      REAL,dimension(nz)::tte
      REAL,dimension(nz)::te_old
      REAL,dimension(nz)::kap
      REAL,dimension(nz)::s1,s2,s3,s4,s5
      INTEGER::nfl

      REAL:: a(nz),b(nz),c(nz),d(nz)
      INTEGER::j
!     initialize

      a = 0.
      b = 0.
      c = 0.
      d = 0.

!     note: ne used here is in a common block

      do j = 2,nz-1

        a(j) = - bms(j,nfl)**2 / ne(j,nfl) / d22s(j,nfl)&
         *.5 * ( kap(j) + kap(j-1) ) / ds(j,nfl)

        b(j) = 1. / dt + bms(j,nfl)**2 / ne(j,nfl) / d22s(j,nfl)&
        *(  .5 * (kap(j+1) + kap(j)   ) /ds(j+1,nfl)&
           +.5 * (kap(j)   + kap(j-1) ) /ds(j,nfl)   )&
        + s2(j) + s3(j)

        c(j) = - bms(j,nfl)**2 / ne(j,nfl) /d22s(j,nfl)&
         *.5 * ( kap(j+1) + kap(j) )/ ds(j+1,nfl)

        d(j) = te_old(j)/dt + s1(j) + s4(j) + s5(j)

       enddo

!     we will assume that the bc's are the neutral temperature
!     at both ends of the field line

!     lower bc

      a(1) = 0.
      b(1) = 1.
      c(1) = 0.
      d(1) = tn(1,nfl)

!     upper bc

      a(nz) = 0.
      b(nz) = 1.
      c(nz) = 0.
      d(nz) = tn(nz,nfl)

      call rtryds(a,b,c,d,tte,nz)

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
      implicit none

      REAL::hrut,sdec,cossdec,sinsdec,clat,slat
      INTEGER::nfl
!      geometric variables

!      sdec: solar zenith angle
!      cx:  cos of the zenith angle

      !LOCAL VARIABLES
      INTEGER::i
      REAL::hrl
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

!             EXB

!*******************************************
!*******************************************

       subroutine exb(hrut)

       include 'param-1.00.inc'
       include 'com-1.00.inc'
      use vdrift_module
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
          vexb(i,j) = vexb0 *&
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
                 fluxnp(i,j,ni) = fluxnp(i,j-1,ni)*&
                                  areap(i,j-1)/areap(i,j)       
                 fluxtp(i,j,ni) = fluxtp(i,j-1,ni)*&
                                  areap(i,j-1)/areap(i,j)       

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
                fluxtep(i,j) = fluxtep(i,j-1)*&
                              areap(i,j-1)/areap(i,j)       
                
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
             denic(i,j,ni) = denic(i,j,ni)&
                             + dt * ( areap(i,j)   * fluxnp(i,j,ni) -&
                                      areap(i,j+1) * fluxnp(i,j+1,ni) )&
                             + dt * ( areas(i,j)   * fluxns(i,j,ni) -&
                                      areas(i+1,j) * fluxns(i+1,j,ni) )
             deni(i,j,ni)  = denic(i,j,ni) / vol(i,j)
             tic(i,j,ni) = tic(i,j,ni)&
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
           tec(i,j) = tec(i,j)&
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
!          implicit none
          REAL::hrut
          REAL::hr24ut,dt00,hr24l,dtnew,dt1
          INTEGER::k,j,i
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
!$OMP PARALLEL PRIVATE(I,J,dts,dtp,dt1) SHARED(xdels,xdelp,vexbs,vexbp,dtnew)
!!$OMP DO
       do j = 1,nf
!$OMP DO
         do i = 1,nz
           dts = xdels(i,j,1) / amax1(1.,abs(vexbs(i,j)))
           dtp = xdelp(i,j,1) / amax1(1.,abs(vexbp(i,j)))
           dt1 = amin1 ( dts,dtp )
!$OMP CRITICAL
           if ( dt1 .lt. dtnew ) dtnew = dt1
!$OMP END CRITICAL

         enddo
!$OMP END DO
       enddo 
!!$OMP END DO
!$OMP END PARALLEL

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

!             xerfcexp

!*******************************************
!*******************************************

        function xerfcexp(x)

        include 'param-1.00.inc'
      use vdrift_module
          t          = 1. / (1 + pas * x)
          xerfcexp   = (   z1 * t&
                        + z2 * t ** 2&
                        + z3 * t ** 3&
                        + z4 * t ** 4&
                        + z5 * t ** 5  )

        return
        end




end module commonsubroutines
