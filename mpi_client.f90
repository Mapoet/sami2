module mpi_client
      integer,save::left,right,taskid,numtasks,ierr
      
contains
      subroutine share_data
      use parameters
      use commons
      implicit none
      include "mpif.h"
      integer::status(MPI_STATUS_SIZE)

         call MPI_RECV(grs,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(glats,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(glons,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(bms,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(gs,nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(ps,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(blats,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(coschicrit,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(ds,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(d2s,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(d22s,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(dels,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(grad_inp,nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)


          call MPI_RECV(xnorms,nzp1*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(ynorms,nzp1*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(znorms,nzp1*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(xnormp,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(ynormp,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(znormp,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(arg,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(athg,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(aphig,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

          call MPI_RECV(delsp,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vol,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(areap,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(areas,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vnx,nzp1*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vny,nzp1*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vnz,nzp1*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(xdels,nz*nfp1*2,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(xdelp,nzp1*nf*2,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vexbs,nzp1*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vexbp,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vexb,nzp1*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

          call MPI_RECV(deni,nz*nf*nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(denn,nz*nf*nneut,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(ne,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vsi,nz*nf*nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vsid,nz*nf*nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(sumvsi,nz*nf*nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vsic,nz*nf*nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(te,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(ti,nz*nf*nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(tn,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(u,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(v,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vpi,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

          call MPI_RECV(vot,nz*nf*nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vor,nz*nf*nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

          call MPI_RECV(cx,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

          call MPI_RECV(fluxnt,nz*nf*91*linesnt,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

          call MPI_RECV(t1,nz*nf*nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(t2,nz*nf*nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(t3,nz*nf*nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(u1,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(u2,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(u3,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(u4,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(u5,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

      end subroutine share_data
end module mpi_client