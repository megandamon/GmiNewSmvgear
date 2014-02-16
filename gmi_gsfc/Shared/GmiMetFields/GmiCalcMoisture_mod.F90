!-------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:
!
      module GmiCalcMoisture_mod
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: CalcMoistureChanges
!
! !AUTHOR:
! Jules Kouatchou, Jules.Kouatchou-1@nasa.gov
!
!EOP
!-------------------------------------------------------------------------
      contains
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcMoistureChanges
!
! !INTERFACE:
!
      subroutine CalcMoistureChanges &
     &  (moistq, rain, mcor, mass, i1, i2, ju1, j2, k1, k2, ivert)
!
      implicit none
!
#     include "gmi_phys_constants.h"
!
! !INPUT PARAMETERS:
      integer, intent(in ) :: i1, i2, ju1, j2, k1, k2, ivert
      real*8 , intent(in ) :: mcor   (i1:i2,ju1:j2)
      real*8 , intent(in ) :: mass   (i1:i2,ju1:j2,k1:k2)
      real*8 , intent(in ) :: rain   (i1:i2,ju1:j2,k1:k2)
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: moistq (i1:i2,ju1:j2,k1:k2)
!
! !DEFINED PARAMETERS:
      real*8, parameter :: CM2PM2    = CMPM * CMPM  ! cm^2/m^2
      real*8, parameter :: CMPMM     = 0.1d0        ! cm/mm
      real*8, parameter :: WATER_DEN = 1.0d0        ! gr/cm^3
!
! !DESCRIPTION:
!     moistq (moisture change) needs to be calculated from the difference
!     in the rain variable across the grid box.  The units need to be
!     changed from mm/day to g/kg/day.
!EOP
!-------------------------------------------------------------------------
!BOC
      moistq(:,:,k1:k2-1) =  (rain(:,:,k1+1:k2) - rain(:,:,k1:k2-1)) *  &
     &          Spread (mcor(:,:), 3, ivert-1) * CM2PM2 * CMPMM * WATER_DEN /  &
     &          mass(:,:,k1:k2-1)

      moistq(:,:,k2) = 0.0d0

      return

      end subroutine CalcMoistureChanges
!EOC
!-------------------------------------------------------------------------

end module GmiCalcMoisture_mod
