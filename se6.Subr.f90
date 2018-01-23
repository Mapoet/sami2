        subroutine s6esubr(nfl,s6e)
        include 'param-1.00.inc'
        include 'com-1.00.inc'
        integer nfl
        real,dimension(nz)::s6e

          do i = 1,nz
               s6e(i) = 0.
          enddo

        return
        end
