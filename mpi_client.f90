module mpi_client
      integer,save::left,right,taskid,numtasks,ierr
      
contains
      subroutine share_data_server(nfl,nfr,i)
      use parameters
      use commons
      implicit none
      include "mpif.h"
      INTEGER,INTENT(IN)::nfl,nfr,i
      integer::status(MPI_STATUS_SIZE)
      INTEGER::size,sizep1



          call MPI_SEND(hrinit,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(maxstep,1,MPI_INT,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(hrmax,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

          call MPI_SEND(year,1,MPI_INT,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(day,1,MPI_INT,i,0,MPI_COMM_WORLD,status,ierr)

          call MPI_SEND(dt,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)



      size=(nfr-nfl+1)
      sizep1=size+1
          call MPI_SEND(alts(:,nfl:nfr),size*nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(glats(:,nfl:nfr),size*nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(grs(:,nfl:nfr),size*nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

          call MPI_SEND(glons(:,nfl:nfr),size*nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(bms(:,nfl:nfr),size*nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(gs(:),nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(ps(:,nfl:nfr),size*nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(blats(:,nfl:nfr),size*nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(coschicrit(:,nfl:nfr),size*nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(ds(:,nfl:nfr),size*nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(d2s(:,nfl:nfr),size*nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(d22s(:,nfl:nfr),size*nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(dels(:,nfl:nfr),size*nz,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(grad_inp(nfl:nfr),size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)


          call MPI_SEND(xnorms(:,nfl:nfr),nzp1*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(ynorms(:,nfl:nfr),nzp1*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(znorms(:,nfl:nfr),nzp1*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(xnormp,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(ynormp,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(znormp,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(arg(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(athg(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(aphig(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

          call MPI_SEND(delsp,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vol(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(areap,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(areas,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vnx,nzp1*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vny,nzp1*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vnz,nzp1*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(xdels,nz*sizep1*2,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(xdelp(:,nfl:nfr,:),nzp1*size*2,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vexbs(:,nfl:nfr),nzp1*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vexbp,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vexb,nzp1*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

          call MPI_SEND(deni(:,nfl:nfr,:),nz*size*nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(denn(:,nfl:nfr,:),nz*size*nneut,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(ne(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vsi(:,nfl:nfr,:),nz*size*nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vsid(:,nfl:nfr,:),nz*size*nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(sumvsi(:,nfl:nfr,:),nz*size*nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vsic(:,nfl:nfr,:),nz*size*nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(te(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(ti(:,nfl:nfr,:),nz*size*nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(tn(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(u(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(v(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vpi(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

          call MPI_SEND(vot(:,nfl:nfr,:),nz*size*nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vor(:,nfl:nfr,:),nz*size*nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

          call MPI_SEND(cx(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

          call MPI_SEND(fluxnt(:,nfl:nfr,:,:),nz*size*91*linesnt,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

          call MPI_SEND(t1(:,nfl:nfr,:),nz*size*nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(t2(:,nfl:nfr,:),nz*size*nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(t3(:,nfl:nfr,:),nz*size*nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(u1(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(u2(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(u3(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(u4(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(u5(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

      end subroutine share_data_server
      

      subroutine share_data_client
      use parameters
      use commons
      implicit none
      include "mpif.h"
      integer::status(MPI_STATUS_SIZE)
      INTEGER::size

          call MPI_RECV(hrinit,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(maxstep,1,MPI_INT,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(hrmax,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

          call MPI_RECV(year,1,MPI_INT,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(day,1,MPI_INT,0,0,MPI_COMM_WORLD,status,ierr)

          call MPI_RECV(dt,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)


          call MPI_RECV(alts,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(glats,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(grs,nf*nz,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          
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

      end subroutine share_data_client
end module mpi_client