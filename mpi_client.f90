module mpi_client
      integer,save::left,right,taskid,numtasks,ierr
      
contains
      subroutine share_data
         
         implicit none
         include "mpif.h"
         
      end subroutine share_data
end module mpi_client