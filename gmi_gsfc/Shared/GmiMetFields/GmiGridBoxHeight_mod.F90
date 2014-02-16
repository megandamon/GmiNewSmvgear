   module GmiGridBoxHeight_mod

      implicit none

      private
!      public  :: CalcBoxHeightCenter
!      public  :: CalcBoxHeightEdge
      public  :: CalcGridBoxHeight

   contains

!-----------------------------------------------------------------------------
!
! ROUTINE
!   CalcGridBoxHeight
!
! DESCRIPTION
!   This routine calculates the height of each grid box (m).          
!
! ARGUMENTS
!-----------------------------------------------------------------------------

      subroutine CalcGridBoxHeight  &
     &  (press3e, surfPressure, humidity, temperature, gridBoxHeight, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in ) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      real*8 , intent(in ) :: surfPressure (ilo:ihi, julo:jhi)
      real*8 , intent(in ) :: humidity   (i1:i2,   ju1:j2,   k1:k2)
      real*8 , intent(in ) :: temperature     (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in ) :: press3e    (ilo:ihi, julo:jhi, k1-1:k2)
      real*8 , intent(out) :: gridBoxHeight(i1:i2,ju1:j2,k1:k2)

      integer :: ik

!     ----------------
!     Begin execution.
!     ----------------

      do ik = k1, k2
        gridBoxHeight(:,:,ik) =  &
     &      29.3d0 * temperature(i1:i2,ju1:j2,ik) *  &
     &      (1.0d0 + (0.61d-3 * humidity(:,:,ik))) *  &
     &      (Log (surfPressure(i1:i2,ju1:j2) /  press3e(i1:i2,ju1:j2,ik  )) - &
     &       Log (surfPressure(i1:i2,ju1:j2) /  press3e(i1:i2,ju1:j2,ik-1)))
      enddo

      return

      end  subroutine CalcGridBoxHeight

  end module GmiGridBoxHeight_mod
