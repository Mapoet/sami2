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

      
      do while ( hrut .ge. 24. )
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
      REAL::hrut

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

end module commonsubroutines