!------------------------------------------------------------------------------
! NASA/GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: computeMaxMinArray_mod
! 
   MODULE computeMaxMinArray_mod 
!
! !USES:
!
! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC calcMaxMinArray   ! functional

   interface calcMaxMinArray
     module procedure calcMaxMinArray2d
     module procedure calcMaxMinArray3d
   end interface
!
#     include "mpif.h"
!
!EOP
!---------------------------------------------------------------------------
CONTAINS
!---------------------------------------------------------------------------
!BOP
! 
! !IROUTINE: calcMaxMinArray3d
!
     subroutine calcMaxMinArray3d ( qname, a, pmin, pmax, im, jt, fac, &
                                    procID, commuWorld )
      implicit none
!
! !INPUT PARAMETERS:
      character*(*), intent(in) :: qname      ! name of the variable
      integer,       intent(in) :: procID     ! processor ID
      integer,       intent(in) :: commuWorld ! MPI world communicator
      integer,       intent(in) :: im         ! number of horizontal grid points (ixj)
      integer,       intent(in) :: jt         ! number of vertical levels
      real*8,        intent(in) :: a(:,:,:)
      real*8,        intent(in) :: fac        ! multiplication factor
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: pmax ! maximum value of the array
      real*8 , intent(out) :: pmin ! minimum value of the array
!EOP
!---------------------------------------------------------------------------
!BOC
      call calcMaxMinArray2d ( qname, reshape(a,(/ im, jt /)), &
                             pmin, pmax, im, jt, fac, procID, commuWorld)

      return

      end subroutine calcMaxMinArray3d
!EOC
!---------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcMaxMinArray2d
!
      subroutine calcMaxMinArray2d ( qname, a, pmin, pmax, im, jt, fac, &
                                     procID, commuWorld )

      implicit none

!
! !INPUT PARAMETERS:
      character*(*), intent(in) :: qname      ! name of the variable
      integer,       intent(in) :: procID     ! processor ID
      integer,       intent(in) :: commuWorld ! MPI world communicator
      integer,       intent(in) :: im         ! number of horizontal grid points (ixj)
      integer,       intent(in) :: jt         ! number of vertical levels
      real*8,        intent(in) :: a(im,jt)
      real*8,        intent(in) :: fac        ! multiplication factor
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: pmax ! maximum value of the array
      real*8 , intent(out) :: pmin ! minimum value of the array
!
! !LOCAL VARIABLES:
      integer :: i, j, two=2

      integer :: rootProc = 0   ! root processor ID
      real*8 :: qmin(jt), qmax(jt)
      real*8 :: pm1(2)
      real*8 :: pm_res(2)

      character(len=16) :: name
      integer :: status
!EOP
!---------------------------------------------------------------------------
!BOC

      do j=1,jt
         pmax = a(1,j)
         pmin = a(1,j)
         do i=2,im
            pmax = max(pmax, a(i,j))
            pmin = min(pmin, a(i,j))
         enddo
         qmax(j) = pmax
         qmin(j) = pmin
      enddo
!
! Now find max/min of amax/amin
!
      pmax = qmax(1)
      pmin = qmin(1)
      do j=2,jt
         pmax = max(pmax, qmax(j))
         pmin = min(pmin, qmin(j))
      enddo

      pm1(1) = pmax
      pm1(2) = -pmin

      call MPI_Reduce (pm1, pm_res, two, MPI_REAL8, MPI_MAX, rootProc, &
                       commuWorld, status)

      pmax=pm_res(1)
      pmin=-pm_res(2)

      if ( fac /= 0.0 ) then  ! trick to prevent printing
         if ( procID == rootProc ) then
            name = '            '
            name(1:len(qname)) = qname
            write(*,*) name, ' max = ', pmax*fac, ' min = ', pmin*fac
            return
         end if
      end if

      return

    end subroutine calcMaxMinArray2d
!EOC
!---------------------------------------------------------------------------
   END MODULE computeMaxMinArray_mod 
