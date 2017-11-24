module mpi_client
      INTEGER,PARAMETER::MPI_TAG_DT_SYNC=1001
      INTEGER,PARAMETER::MPI_TAG_OUTPUT_DATA_SYNC=1002
      INTEGER,PARAMETER::MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC=2000
      INTEGER,PARAMETER::MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC=3000
      integer,save::left,right,taskid,numtasks,ierr
      
contains
      subroutine share_data_server(nfl,nfr,i)
      use parameters
      use commons
      use TRACE
      implicit none
      include "mpif.h"
      INTEGER,INTENT(IN)::nfl,nfr,i
      integer::status(MPI_STATUS_SIZE)
      INTEGER::size,sizep1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call MPI_SEND(fejer,1,MPI_LOGICAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(fmtout,1,MPI_LOGICAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(snn,nneut,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(maxstep,1,MPI_INT,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(hrmax,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(dt0,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(dthr,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(hrpr,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(grad_in,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(glat_in,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(glon_in,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(rmin,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(rmax,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(altmin,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(fbar,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(f10p7,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(ap,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(year,1,MPI_INT,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(day,1,MPI_INT,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(mmass,1,MPI_INT,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(nion1,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(nion2,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(hrinit,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(tvn0,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(tvexb0,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(ve01,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(gams,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(gamp,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(stn,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(denmin,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(alt_crit,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(cqe,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

      call MPI_SEND( dt,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

      call MPI_SEND(x0,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(y0,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(z0,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(plat,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(plon,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
      call MPI_SEND(bb0,1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
          !call MPI_SEND(xnormp,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_SEND(ynormp,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_SEND(znormp,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(arg(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(athg(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(aphig(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

          !call MPI_SEND(delsp,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vol(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_SEND(areap,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
         ! call MPI_SEND(areas,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
         ! call MPI_SEND(vnx,nzp1*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
         ! call MPI_SEND(vny,nzp1*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
         ! call MPI_SEND(vnz,nzp1*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(xdels(:,nfl:nfr,:),nz*size*2,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(xdelp(:,nfl:nfr,:),nzp1*size*2,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_SEND(vexbs(:,nfl:nfr),nzp1*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_SEND(vexbp,nz*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_SEND(vexb,nzp1*sizep1,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

            call MPI_SEND(ichem,nchem*3,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)!02 to share
            call MPI_SEND(ireact,nion*nneut*nchem,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

          call MPI_SEND(deni(:,nfl:nfr,:),nz*size*nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)!01 to share
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

          call MPI_SEND(ami(:),nion,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(amn(:),nneut,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(alpha0(:),nneut,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(aap(:),7,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)

          call MPI_SEND(cx(:,nfl:nfr),nz*size,MPI_REAL,i,0,MPI_COMM_WORLD,status,ierr)!Evaluated for each client by itself

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call MPI_RECV(fejer,1,MPI_LOGICAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(fmtout,1,MPI_LOGICAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(snn,nneut,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(maxstep,1,MPI_INT,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(hrmax,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(dt0,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(dthr,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(hrpr,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(grad_in,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(glat_in,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(glon_in,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(rmin,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(rmax,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(altmin,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(fbar,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(f10p7,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(ap,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(year,1,MPI_INT,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(day,1,MPI_INT,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(mmass,1,MPI_INT,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(nion1,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(nion2,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(hrinit,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(tvn0,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(tvexb0,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(ve01,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(gams,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(gamp,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(stn,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(denmin,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(alt_crit,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(cqe,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

      call MPI_RECV( dt,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

      call MPI_RECV(x0,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(y0,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(z0,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(plat,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(plon,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
      call MPI_RECV(bb0,1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
          !call MPI_RECV(xnormp,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_RECV(ynormp,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_RECV(znormp,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(arg,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(athg,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(aphig,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

          !call MPI_RECV(delsp,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(vol,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_RECV(areap,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_RECV(areas,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_RECV(vnx,nzp1*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_RECV(vny,nzp1*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_RECV(vnz,nzp1*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(xdels,nz*nf*2,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(xdelp,nzp1*nf*2,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_RECV(vexbs,nzp1*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_RECV(vexbp,nz*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          !call MPI_RECV(vexb,nzp1*nfp1,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

            call MPI_RECV(ichem,nchem*3,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
            call MPI_RECV(ireact,nion*nneut*nchem,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)


          call MPI_RECV(deni,nz*nf*nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)!01 to share
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
          
          call MPI_RECV(ami,nion,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(amn,nneut,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(alpha0,nneut,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)
          call MPI_RECV(aap,7,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)

          call MPI_RECV(cx,nz*nf,MPI_REAL,0,0,MPI_COMM_WORLD,status,ierr)!Evaluated for each client by itself

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

      subroutine share_output_data_client
      use parameters
      use commons
      implicit none
      include "mpif.h"
      !INTEGER,INTENT(IN)::nfl,nfr,i
      integer::status(MPI_STATUS_SIZE)
      !INTEGER::size,sizep1
            
          call MPI_SEND(deni,nz*nf*nion,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)!01 to share
          call MPI_SEND(denn,nz*nf*nneut,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(ne,nz*nf,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vsi,nz*nf*nion,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vsid,nz*nf*nion,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(sumvsi,nz*nf*nion,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(vsic,nz*nf*nion,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(te,nz*nf,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(ti,nz*nf*nion,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
          call MPI_SEND(tn,nz*nf,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
          !call MPI_SEND(u,nz*nf,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
          !call MPI_SEND(v,nz*nf,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
          !call MPI_SEND(vpi,nz*nf,MPI_REAL,0,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)


      end subroutine share_output_data_client
      subroutine share_output_data_server(nfsize,nfstindex)
      use parameters
      use commons
      implicit none
      include "mpif.h"
      INTEGER,DImEnsion(:),INTENT(IN)::nfsize,nfstindex
      integer::status(MPI_STATUS_SIZE)
      INTEGER::SRCwrkr,nfsizewrkr,nfstindexwrkr,realsizewrkr
      !INTEGER::size,sizep1
      real,save,dimension(:),allocatable::tmp_buf1
      real,save,dimension(:,:),allocatable::tmp_buf2
      real,save,dimension(:,:,:),allocatable::tmp_buf3
      real,save,dimension(:,:,:,:),allocatable::tmp_buf4

            call MPI_Probe(MPI_ANY_SOURCE, MPI_TAG_OUTPUT_DATA_SYNC, MPI_COMM_WORLD, status, ierr)
            SRCwrkr=status(MPI_SOURCE)
            
            nfsizewrkr=nfsize(SRCwrkr)
            nfstindexwrkr=nfstindex(SRCwrkr)
            realsizewrkr=nfsizewrkr+merge(1,0,SRCwrkr.ne.1)+merge(1,0,SRCwrkr.ne.numtasks-1)


            ALLOCATE(tmp_buf3(nz,realsizewrkr,nion)) 
            call MPI_RECV(tmp_buf3,nz*realsizewrkr*nion,MPI_REAL,SRCwrkr,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)!01 to share
            deni(:,nfstindexwrkr:(nfstindexwrkr+nfsizewrkr-1),:)=tmp_buf3(:,merge(1,0,SRCwrkr.ne.1)+1:merge(1,0,SRCwrkr.ne.1)+nfsizewrkr,:)
            DeALLOCATE(tmp_buf3(nz,realsizewrkr,nion)) 

            ALLOCATE(tmp_buf3(nz,realsizewrkr,nneut)) 
            call MPI_RECV(tmp_buf3,nz*realsizewrkr*nneut,MPI_REAL,SRCwrkr,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)!01 to share
            denn(:,nfstindexwrkr:(nfstindexwrkr+nfsizewrkr-1),:)=tmp_buf3(:,merge(1,0,SRCwrkr.ne.1)+1:merge(1,0,SRCwrkr.ne.1)+nfsizewrkr,:)
            DeALLOCATE(tmp_buf3(nz,realsizewrkr,nneut)) 

            ALLOCATE(tmp_buf2(nz,realsizewrkr)) 
            call MPI_RECV(tmp_buf2,nz*realsizewrkr,MPI_REAL,SRCwrkr,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
            ne(:,nfstindexwrkr:(nfstindexwrkr+nfsizewrkr-1))=tmp_buf2(:,merge(1,0,SRCwrkr.ne.1)+1:merge(1,0,SRCwrkr.ne.1)+nfsizewrkr)
            DeALLOCATE(tmp_buf2(nz,realsizewrkr)) 

            ALLOCATE(tmp_buf3(nz,realsizewrkr,nion)) 
            call MPI_RECV(tmp_buf3,nz*realsizewrkr*nion,MPI_REAL,SRCwrkr,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)!01 to share
            vsi(:,nfstindexwrkr:(nfstindexwrkr+nfsizewrkr-1),:)=tmp_buf3(:,merge(1,0,SRCwrkr.ne.1)+1:merge(1,0,SRCwrkr.ne.1)+nfsizewrkr,:)
            DeALLOCATE(tmp_buf3(nz,realsizewrkr,nion)) 

            ALLOCATE(tmp_buf3(nz,realsizewrkr,nion)) 
            call MPI_RECV(tmp_buf3,nz*realsizewrkr*nion,MPI_REAL,SRCwrkr,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)!01 to share
            vsid(:,nfstindexwrkr:(nfstindexwrkr+nfsizewrkr-1),:)=tmp_buf3(:,merge(1,0,SRCwrkr.ne.1)+1:merge(1,0,SRCwrkr.ne.1)+nfsizewrkr,:)
            DeALLOCATE(tmp_buf3(nz,realsizewrkr,nion)) 

            ALLOCATE(tmp_buf3(nz,realsizewrkr,nion)) 
            call MPI_RECV(tmp_buf3,nz*realsizewrkr*nion,MPI_REAL,SRCwrkr,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)!01 to share
            sumvsi(:,nfstindexwrkr:(nfstindexwrkr+nfsizewrkr-1),:)=tmp_buf3(:,merge(1,0,SRCwrkr.ne.1)+1:merge(1,0,SRCwrkr.ne.1)+nfsizewrkr,:)
            DeALLOCATE(tmp_buf3(nz,realsizewrkr,nion)) 

            ALLOCATE(tmp_buf3(nz,realsizewrkr,nion)) 
            call MPI_RECV(tmp_buf3,nz*realsizewrkr*nion,MPI_REAL,SRCwrkr,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)!01 to share
            vsic(:,nfstindexwrkr:(nfstindexwrkr+nfsizewrkr-1),:)=tmp_buf3(:,merge(1,0,SRCwrkr.ne.1)+1:merge(1,0,SRCwrkr.ne.1)+nfsizewrkr,:)
            DeALLOCATE(tmp_buf3(nz,realsizewrkr,nion)) 


            ALLOCATE(tmp_buf2(nz,realsizewrkr)) 
            call MPI_RECV(tmp_buf2,nz*realsizewrkr,MPI_REAL,SRCwrkr,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
            te(:,nfstindexwrkr:(nfstindexwrkr+nfsizewrkr-1))=tmp_buf2(:,merge(1,0,SRCwrkr.ne.1)+1:merge(1,0,SRCwrkr.ne.1)+nfsizewrkr)
            DeALLOCATE(tmp_buf2(nz,realsizewrkr)) 

            ALLOCATE(tmp_buf3(nz,realsizewrkr,nion)) 
            call MPI_RECV(tmp_buf3,nz*realsizewrkr*nion,MPI_REAL,SRCwrkr,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)!01 to share
            ti(:,nfstindexwrkr:(nfstindexwrkr+nfsizewrkr-1),:)=tmp_buf3(:,merge(1,0,SRCwrkr.ne.1)+1:merge(1,0,SRCwrkr.ne.1)+nfsizewrkr,:)
            DeALLOCATE(tmp_buf3(nz,realsizewrkr,nion)) 

            ALLOCATE(tmp_buf2(nz,realsizewrkr)) 
            call MPI_RECV(tmp_buf2,nz*realsizewrkr,MPI_REAL,SRCwrkr,MPI_TAG_OUTPUT_DATA_SYNC,MPI_COMM_WORLD,status,ierr)
            tn(:,nfstindexwrkr:(nfstindexwrkr+nfsizewrkr-1))=tmp_buf2(:,merge(1,0,SRCwrkr.ne.1)+1:merge(1,0,SRCwrkr.ne.1)+nfsizewrkr)
            DeALLOCATE(tmp_buf2(nz,realsizewrkr)) 




      end subroutine share_output_data_server
      subroutine share_data_btwn_clients!(nfl,nfr,i)
      use parameters
      use commons
      implicit none
      include "mpif.h"
      integer::status(MPI_STATUS_SIZE)
      INTEGER::targetleft,targetright,datanum
      INTEGER:: mpirequest

          targetleft=taskid-1
          datanum=0
          if(targetleft.ge.1) then
               call MPI_iSEND(deni(:,2,:),nz*1*nion,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD, mpirequest,ierr)
               datanum=datanum+1
               call MPI_iSEND(denn(:,2,:),nz*1*nneut,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD, mpirequest,ierr)
               datanum=datanum+1
               call MPI_iSEND(ne(:,2),nz*1,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
               datanum=datanum+1
!               call MPI_iSEND(vsi(:,2,:),nz*1*nneut,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD, mpirequest,ierr)
!               datanum=datanum+1
!               call MPI_iSEND(vsid(:,2,:),nz*1*nneut,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD, mpirequest,ierr)
!               datanum=datanum+1
!               call MPI_iSEND(sumvsi(:,2,:),nz*1*nion,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD, mpirequest,ierr)
!               datanum=datanum+1
!               call MPI_iSEND(vsic(:,2,:),nz*1*nion,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD, mpirequest,ierr)
!               datanum=datanum+1
               
               call MPI_iSEND(te(:,2),nz*1,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
               datanum=datanum+1
               call MPI_iSEND(ti(:,2,:),nz*1*nion,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD, mpirequest,ierr)
               datanum=datanum+1
               call MPI_iSEND(tn(:,2),nz*1,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
               datanum=datanum+1
          endif

          targetright=taskid+1
          datanum=0
          if(targetright.le.numtasks-1) then
               call MPI_iSEND(deni(:,nf-1,:),nz*1*nion,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
               datanum=datanum+1
               call MPI_iSEND(denn(:,nf-1,:),nz*1*nneut,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
               datanum=datanum+1
               call MPI_iSEND(ne(:,nf-1),nz*1,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
               datanum=datanum+1
!               call MPI_iSEND(vsi(:,nf-1,:),nz*1*nion,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
!               datanum=datanum+1
!               call MPI_iSEND(vsid(:,nf-1,:),nz*1*nion,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
!               datanum=datanum+1
!               call MPI_iSEND(sumvsi(:,nf-1,:),nz*1*nion,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
!               datanum=datanum+1
!               call MPI_iSEND(vsic(:,nf-1,:),nz*1*nion,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
!               datanum=datanum+1

               call MPI_iSEND(te(:,nf-1),nz*1,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
               datanum=datanum+1
               call MPI_iSEND(ti(:,nf-1,:),nz*1*nion,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
               datanum=datanum+1
               call MPI_iSEND(tn(:,nf-1),nz*1,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,mpirequest,ierr)
               datanum=datanum+1
          endif

          datanum=0
          if(targetright.le.numtasks-1) then
               call MPI_RECV(deni(:,nf,:),nz*1*nion,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
               datanum=datanum+1
               call MPI_RECV(denn(:,nf,:),nz*1*nneut,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
               datanum=datanum+1
               call MPI_RECV(ne(:,nf),nz*1,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
               datanum=datanum+1
!               call MPI_RECV(vsi(:,nf,:),nz*1*nion,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
!               datanum=datanum+1
!               call MPI_RECV(vsid(:,nf,:),nz*1*nion,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
!               datanum=datanum+1
!               call MPI_RECV(sumvsi(:,nf,:),nz*1*nion,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
!               datanum=datanum+1
!               call MPI_RECV(vsic(:,nf,:),nz*1*nion,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
!               datanum=datanum+1

               call MPI_RECV(te(:,nf),nz*1,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
               datanum=datanum+1
               call MPI_RECV(ti(:,nf,:),nz*1*nion,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
               datanum=datanum+1
               call MPI_RECV(tn(:,nf),nz*1,MPI_REAL,targetright,MPI_TAG_SHARE_CLIENT_DATA_LEFT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
               datanum=datanum+1
          endif

          datanum=0
          if(targetleft.ge.1) then
               call MPI_RECV(deni(:,1,:),nz*1*nion,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
               datanum=datanum+1
               call MPI_RECV(denn(:,1,:),nz*1*nneut,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
               datanum=datanum+1
               call MPI_RECV(ne(:,1),nz*1,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
               datanum=datanum+1
!               call MPI_RECV(vsi(:,1,:),nz*1*nion,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
!               datanum=datanum+1
!               call MPI_RECV(vsid(:,1,:),nz*1*nion,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
!               datanum=datanum+1
!               call MPI_RECV(sumvsi(:,1,:),nz*1*nion,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
!               datanum=datanum+1
!               call MPI_RECV(vsic(:,1,:),nz*1*nion,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
!               datanum=datanum+1

               call MPI_RECV(te(:,1),nz*1,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
               datanum=datanum+1
               call MPI_RECV(ti(:,1,:),nz*1*nion,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
               datanum=datanum+1
               call MPI_RECV(tn(:,1),nz*1,MPI_REAL,targetleft,MPI_TAG_SHARE_CLIENT_DATA_RIGHT_SYNC+datanum,MPI_COMM_WORLD,status,ierr)
               datanum=datanum+1
          endif
      end subroutine share_data_btwn_clients
end module mpi_client