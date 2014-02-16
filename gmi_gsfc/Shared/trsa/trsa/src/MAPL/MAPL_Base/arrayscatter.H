! $Id: arrayscatter.H,v 1.2 2011-08-09 22:13:00 mrdamon Exp $

#ifdef NAME_
#undef NAME_
#endif

#define NAME_ ArrayScatter_

#include "overload.macro"

  subroutine SUB_(local_array, global_array, grid, mask, depe, hw, rc)

! Mask is really a permutation on the first dimension

    TYPE_(kind=EKIND_),         intent(  OUT) :: local_array  DIMENSIONS_
    TYPE_(kind=EKIND_), target, intent(IN   ) :: global_array DIMENSIONS_
    type (ESMF_Grid)                          :: grid
    integer, optional,          intent(IN   ) :: mask(:)
    integer, optional,          intent(IN   ) :: depe
    integer, optional,          intent(IN   ) :: hw
    integer, optional,          intent(  OUT) :: rc
    
! Local variables

    integer                               :: status
    character(len=ESMF_MAXSTR)            :: IAm='ArrayScatter'

    TYPE_(kind=EKIND_),    pointer        :: myglob DIMENSIONS_ => null()
    TYPE_(kind=EKIND_),    pointer        :: VAR(:)
    type(ESMF_DELayout)                   :: LAYOUT
    type(ESMF_AxisIndex),  pointer        :: AI(:,:)
    integer, dimension(:), allocatable    :: SENDCOUNTS, DISPLS, KK
    integer                               :: nDEs
    integer                               :: recvcount
    integer                               :: I, J, K, II,  JJ, de, deId
    integer                               :: I1, IN, J1, JN
    integer                               :: gridRank
    integer                               :: LX, LY
    integer                               :: srcPE
    integer                               :: MYHW, ISZ, JSZ, ISZL

! Works only on 1D and D arrays

    ASSERT_(RANK_ <= 2)

! Optional change of source PE. Default=root

    if(present(depe)) then
       srcPE = depe 
    else
       srcPE = root
    end if

! Optional single halo width

    if(present(hw)) then
       myhw = hw
    else
       myhw = 0
    end if

! Some halo limitations

    if(myhw > 0) then
       ASSERT_(RANK_ == 2        ) ! No halo allowed on 1D
       ASSERT_(.not.present(MASK)) ! No halo allowed if 1st dim is permutted
    end if

! Get grid and layout information

    call ESMF_GridGet    (GRID,   dimCount=gridRank, rc=STATUS);VERIFY_(STATUS)
    call ESMF_GridGet    (GRID,   delayout=layout,   rc=STATUS);VERIFY_(STATUS)
    call ESMF_DELayoutGet(LAYOUT, localDE =deId,     rc=STATUS);VERIFY_(STATUS)
    call ESMF_DELayoutGet(LAYOUT, deCount =nDEs,     rc=STATUS);VERIFY_(STATUS)

    allocate (AI(1:nDEs,gridRank),  stat=status); VERIFY_(STATUS)
    allocate (sendcounts(0:nDEs-1), stat=status); VERIFY_(STATUS)

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

    ISZ = size(GLOBAL_ARRAY,1)

#if (RANK_ == 2)
       JSZ = size(GLOBAL_ARRAY,2)
#else
       JSZ = 1
#endif

! Compute count to be sent to each PE

    if(present(mask)) then
       sendcounts = 0
       do II = 1,ISZ
          sendcounts(mask(ii)) = sendcounts(mask(ii)) + 1
       enddo
       sendcounts = sendcounts*JSZ

    else
       do I = 1,nDEs
          LX = AI(I,1)%max - AI(I,1)%min + 1 + 2*MYHW
#if (RANK_ == 1)
          sendcounts(I-1) = LX
#else
          LY = AI(I,2)%max - AI(I,2)%min + 1 + 2*MYHW
          sendcounts(I-1) = LX*LY
#endif
       end do
    end if

! Count I will recieve

    recvcount = sendcounts(deID)

! Put VAR together at the srcPE 

    if (deId == srcPE) then

       allocate(DISPLS(0:nDEs          ), stat=status); VERIFY_(STATUS)

! Compute displacements into the VAR vector

       displs(0) = 0
       do I = 1,nDEs
          displs(I) = displs(I-1) + sendcounts(I-1)
       end do
!ALT       ASSERT_(displs(nDEs) == (ISZ+2*myhw)*(JSZ+2*myhw))

! If there is a halo, make a haloed copy of the global array.
!   otherwise just copy the pointer.

       myglob => global_array
       
#if (RANK_ == 2)
       if (myhw > 0) then
          allocate(myglob(1-myhw:isz+myhw,1-myhw:jsz+myhw), stat=status)
          VERIFY_(STATUS)
          myglob(1:isz,1:jsz) = GLOBAL_ARRAY

! Fill the halo (I is cyclic)
       
          do j=1,myhw
             myglob(1  -j,:) = myglob(isz-j+1,:)
             myglob(isz+j,:) = myglob(    j  ,:)
             myglob(:,1  -j) = MAPL_Undef
             myglob(:,jsz+j) = MAPL_Undef
          enddo
       endif
#endif

! Fill the VAR vector
       
       if (present(mask)) then
          allocate(VAR(0:displs(nDEs)-1), stat=status); VERIFY_(STATUS)
          allocate(KK (0:nDEs-1        ), stat=status); VERIFY_(STATUS)
          KK = DISPLS(0:nDEs-1)

          do I=1,ISZ
             K = MASK(I)
             II = KK(K)
#if (RANK_ == 1)
             VAR(II) = MYGLOB(I)
#else
             LX = AI(K+1,1)%max - AI(K+1,1)%min + 1 
             do J=1,JSZ
                VAR(II+LX*(J-1)) = MYGLOB(I,J)
             end do
#endif
             KK(MASK(I)) = KK(MASK(I)) + 1 
          end do

          deallocate(KK, stat=status); VERIFY_(STATUS)

       else

#if (RANK_ == 1)
          var => myglob
#else
          allocate(VAR(0:displs(nDEs)-1), stat=status); VERIFY_(STATUS)

          do I = 1,nDEs
             I1 = AI(I,1)%min - myhw
             IN = AI(I,1)%max + myhw
             J1 = AI(I,2)%min - myhw
             JN = AI(I,2)%max + myhw

             K = displs(I-1)
             do JJ=J1,JN
                do II=I1,IN
                   var(K) = MYglob(II,JJ)
                   K = K+1
                end do
             end do
          end do
#endif

       endif !  present(mask)

       if (myhw > 0) then
          deallocate(myglob, stat=status); VERIFY_(STATUS)
       end if
    
    end if !  I am srcPE


! Do the communications
    
    call MAPL_CommsScatterV(layout, var, sendcounts, displs, &
                            local_array, recvcount, srcPE, status)
    VERIFY_(STATUS)

! Clean-up

    if (deId == srcPE) then
       deallocate(displs, stat=status); VERIFY_(STATUS)
       if(present(mask) .or. RANK_ == 2) then
          deallocate(VAR, stat=status); VERIFY_(STATUS)
       end if
    endif

    deallocate(sendcounts, stat=status); VERIFY_(STATUS)
    deallocate(AI,         stat=status); VERIFY_(STATUS)

! All done

    RETURN_(ESMF_SUCCESS)
  end subroutine SUB_

#undef NAME_
#undef DIMENSIONS_
#undef RANK_
#undef VARTYPE_
