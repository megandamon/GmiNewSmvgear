! $Id: scatter.H,v 1.2 2011-08-09 22:12:59 mrdamon Exp $

#ifdef NAME_
#undef NAME_
#endif

#define NAME_ MAPL_CommsScatterV_

#include "overload.macro"

  subroutine SUB_ (layout, sendbuf, sendcnts, displs, &
                   recvbuf, recvcnt, root, RC)
    type (ESMF_DELayout)                         :: layout
    TYPE_(kind=EKIND_), intent(IN   )            :: sendbuf (:)
    integer,            intent(IN   )            :: sendcnts(:)
    TYPE_(kind=EKIND_), intent(  OUT)            :: recvbuf DIMENSIONS_
    integer,            intent(IN   )            :: recvcnt
    integer,            intent(IN   )            :: displs(:)
    integer,            intent(IN   )            :: root
    integer         ,   intent(  OUT),  optional :: RC

    character(len=ESMF_MAXSTR), parameter :: IAM='MAPL_CommsScatter'
    integer                               :: status
    integer                               :: comm
    type(ESMF_VM)                         :: vm

!ALT: there is no ESMF_DELayoutGet(layout, vm=vm) method
! so we get the GLOBAL vm

    call ESMF_VMGetGlobal(vm, rc=status)
    VERIFY_(STATUS)

    call ESMF_VMGet(vm, mpiCommunicator=COMM, rc=status)
    VERIFY_(STATUS)

    call mpi_scatterv( sendbuf, sendcnts, displs, MPITYPE_, &
                       recvbuf, recvcnt,  MPITYPE_,             &
                       root, comm, status )
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  END SUBROUTINE SUB_

#undef NAME_

#undef DIMENSIONS_
#undef RANK_
#undef VARTYPE_

