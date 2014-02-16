! $Id: allgather.H,v 1.2 2011-08-09 22:13:00 mrdamon Exp $

#ifdef NAME_
#undef NAME_
#endif

#define NAME_ MAPL_CommsAllGather_

#include "overload.macro"

  subroutine SUB_ (layout, sendbuf, sendcnt, &
                   recvbuf, recvcnt, RC)
    type (ESMF_DELayout)                         :: layout
    TYPE_,              intent(IN   )            :: sendbuf DIMENSIONS_
    integer,            intent(IN   )            :: sendcnt
    TYPE_,              intent(INOUT)            :: recvbuf DIMENSIONS_
    integer,            intent(IN   )            :: recvcnt
    integer,            intent(  OUT),  optional :: RC

    character(len=ESMF_MAXSTR), parameter :: IAM='MAPL_CommsAllGatherV'
    integer                               :: status
    integer                               :: comm
    type(ESMF_VM)                         :: vm

!ALT: there is no ESMF_DELayoutGet(layout, vm=vm) method
! so we get the GLOBAL vm

    call ESMF_VMGetGlobal(vm, rc=status)
    VERIFY_(STATUS)

    call ESMF_VMGet(vm, mpiCommunicator=COMM, rc=status)
    VERIFY_(STATUS)

    call mpi_allgather ( sendbuf, sendcnt, MPITYPE_, &
                         recvbuf, recvcnt, MPITYPE_, &
                         comm, status )
    VERIFY_(STATUS)


    RETURN_(ESMF_SUCCESS)
  END SUBROUTINE SUB_

#undef NAME_

#undef DIMENSIONS_
#undef RANK_
#undef VARTYPE_

