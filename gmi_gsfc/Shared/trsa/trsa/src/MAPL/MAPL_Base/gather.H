! $Id: gather.H,v 1.2 2011-08-09 22:13:00 mrdamon Exp $

#ifdef NAME_
#undef NAME_
#endif

#define NAME_ MAPL_CommsGatherV_

#include "overload.macro"

  subroutine SUB_ (layout, sendbuf, sendcnt, &
                   recvbuf, recvcnts, displs, root, RC)
    type (ESMF_DELayout)                         :: layout
    TYPE_(kind=EKIND_), intent(IN   )            :: sendbuf DIMENSIONS_
    integer,            intent(IN   )            :: sendcnt
    TYPE_(kind=EKIND_), intent(  OUT)            :: recvbuf (:)
    integer,            intent(IN   )            :: recvcnts(:)
    integer,            intent(IN   )            :: displs(:)
    integer,            intent(IN   )            :: root
    integer         ,   intent(  OUT),  optional :: RC

    character(len=ESMF_MAXSTR), parameter :: IAM='MAPL_CommsGather'
    integer                               :: status
    integer                               :: comm
    type(ESMF_VM)                         :: vm

!ALT: there is no ESMF_DELayoutGet(layout, vm=vm) method
! so we get the GLOBAL vm

    call ESMF_VMGetGlobal(vm, rc=status)
    VERIFY_(STATUS)

    call ESMF_VMGet(vm, mpiCommunicator=COMM, rc=status)
    VERIFY_(STATUS)

    call mpi_gatherv( sendbuf, sendcnt,  MPITYPE_,     &
                  recvbuf, recvcnts, displs, MPITYPE_, &
                  root, comm, status )
    VERIFY_(STATUS)

#ifdef sysLinux
    call mpi_barrier(comm, status)
    VERIFY_(STATUS)
#endif

    RETURN_(ESMF_SUCCESS)
  END SUBROUTINE SUB_

#undef NAME_

#undef DIMENSIONS_
#undef RANK_
#undef VARTYPE_

