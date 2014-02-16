! $Id: arraygather.H,v 1.2 2011-08-09 22:12:59 mrdamon Exp $

#ifdef NAME_
#undef NAME_
#endif

#define NAME_ ArrayGather_

#include "overload.macro"

  subroutine SUB_(local_array, global_array, grid, mask, depe, hw, rc)
    TYPE_(kind=EKIND_), intent(IN   )            :: local_array DIMENSIONS_
    TYPE_(kind=EKIND_), intent(  OUT)            :: global_array DIMENSIONS_
    type (ESMF_Grid)      :: grid
    integer, optional,  intent(IN   )   :: mask(:)
    integer, optional,  intent(IN   )   :: depe
    integer, optional,  intent(IN   )   :: hw
    integer, optional, intent(OUT)     :: rc
    
! Local variables

    integer                               :: status
    character(len=ESMF_MAXSTR)            :: IAm='ArrayGather'

    type (ESMF_DELayout)  :: layout
    type(ESMF_AxisIndex), dimension(:,:), pointer :: AI
    integer, allocatable, dimension(:)            :: recvcounts, displs, kk
    integer                                       :: nDEs
    integer                                       :: sendcount
    
    integer                                       :: I, J, K, II, JJ, LX
    integer                                       :: de, deId
    integer                                       :: NX, NY
    integer                                       :: I1, IN
    integer                                       :: ibeg,iend
    integer                                       :: gridRank
#if (RANK_ > 1)
    integer                                       :: J1, JN
    integer                                       :: k1, k2
    integer                                       :: jbeg,jend
#endif
    integer                                       :: ISZ, JSZ
    integer                                       :: destPE, myhw
    TYPE_(kind=EKIND_), allocatable               :: var(:)
    
! Works only on 1D and D arrays

    ASSERT_(RANK_ <= 2)

    ISZ = size(GLOBAL_ARRAY,1)

#if (RANK_ == 2)
    JSZ = size(GLOBAL_ARRAY,2)
#else
    JSZ = 1
#endif

    if(present(depe)) then
      destPE = depe 
    else
      destPE = root
    end if

    if(present(hw)) then
      myhw = hw
    else
      myhw = 0
    end if

    call ESMF_GridGet(GRID, dimCount=gridRank, rc=STATUS)
    VERIFY_(STATUS)

    call ESMF_GridGet(grid, delayout=layout, rc=status)
    VERIFY_(STATUS)

    call ESMF_DELayoutGet(layout, deCount=nDEs, localDE = deId, rc=status)
    VERIFY_(STATUS)

    allocate (AI(nDEs,gridRank), stat=status)
    VERIFY_(STATUS)

    if (gridRank == 3) then
       call ESMF_GridGetAllAxisIndex(grid, &
                                     horzRelLoc=ESMF_CELL_CENTER, &
                                     vertRelLoc=ESMF_CELL_CELL, &
                                     globalAI=AI, rc=status)
    else
       call ESMF_GridGetAllAxisIndex(grid, &
                                     horzRelLoc=ESMF_CELL_CENTER, &
                                     globalAI=AI, rc=status)
    end if
    VERIFY_(STATUS)

    allocate (recvcounts(nDEs), displs(0:nDEs), stat=status)
    VERIFY_(STATUS)

    allocate(VAR(0:size(GLOBAL_ARRAY)-1), stat=status)
    VERIFY_(STATUS)

    displs(0) = 0
    do I = 1,nDEs
       J = I - 1
       de = J
       I1 = AI(I,1)%min
       IN = AI(I,1)%max
#if (RANK_ > 1)
       J1 = AI(I,2)%min
       JN = AI(I,2)%max
       recvcounts(I) = (IN - I1 + 1) * (JN - J1 + 1)
#else
       recvcounts(I) = (IN - I1 + 1)
#endif
       if (de == deId) then
          sendcount = recvcounts(I)      ! Count I will send
          ibeg = 1+myhw
          iend = IN-I1+1+myhw
#if (RANK_ > 1)
          jbeg = 1+myhw
          jend = JN-J1+1+myhw
#endif
       endif
       displs(I) = displs(J) + recvcounts(I)
    enddo

    if (present(mask) .or. myHW == 0) then
       call MAPL_CommsGatherV(layout, local_array, sendcount, &
                              var, recvcounts, displs, destPE, status)
    else
#if (RANK_ > 1)
       call MAPL_CommsGatherV(layout, local_array(ibeg:iend,jbeg:jend), &
                              sendcount, var, recvcounts, displs, destPE, &
                              status)
#else
       call MAPL_CommsGatherV(layout, local_array(ibeg:iend), sendcount, &
                              var, recvcounts, displs, destPE, status)
#endif
    end if
    VERIFY_(STATUS)

    if (deId == destPE) then
       if (present(mask)) then
          allocate(KK (0:nDEs-1        ), stat=status); VERIFY_(STATUS)
          KK = DISPLS(0:nDEs-1)

          do I=1,ISZ
             K = MASK(I)
             II = KK(K)
#if (RANK_ == 1)
             GLOBAL_ARRAY(I) = VAR(II)
#else
             LX = AI(K+1,1)%max - AI(K+1,1)%min + 1 
             do J=1,JSZ
                GLOBAL_ARRAY(I,J) = VAR(II+LX*(J-1))
             end do
#endif
             KK(MASK(I)) = KK(MASK(I)) + 1 
          end do

          deallocate(KK, stat=status); VERIFY_(STATUS)
       else
#if (RANK_ == 1)
          global_array = var ! ALT: I am not sure if this is correct
#else
          do I = 1,nDEs
             I1 = AI(I,1)%min
             IN = AI(I,1)%max
             J1 = AI(I,2)%min
             JN = AI(I,2)%max

             K = displs(I-1)
             do JJ=J1,JN
                do II=I1,IN
                   global_array(II,JJ) = var(K)
                   K = K+1
                end do
             end do
          end do
#endif
       end if ! if (present(mask))
    end if

    deallocate(VAR, stat=status)
    VERIFY_(STATUS)
    deallocate(recvcounts, displs, AI, stat=status)
    VERIFY_(STATUS)

    RETURN_(ESMF_SUCCESS)
  end subroutine SUB_

#undef NAME_

#undef DIMENSIONS_
#undef RANK_
#undef VARTYPE_
