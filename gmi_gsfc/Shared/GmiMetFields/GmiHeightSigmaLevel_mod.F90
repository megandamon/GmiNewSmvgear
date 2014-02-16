  module GmiHeightSigmaLevel_mod

      implicit none
      
      private 
      public   :: CalcHeightSigmaLevel
      public   :: CalcConvection_bmass
      public   :: CalcDryDepBoxHeight

 contains

!-----------------------------------------------------------------------------
!
! ROUTINE
!   CalcHeightSigmaLevel
!
! DESCRIPTION
!   This routine computes the height of full and half-sigma levels above 
!   ground level. This is used in Llnl emission.
!
! ARGUMENTS
!   temperature         : temperature (degK)
!   press3c     : atmospheric pressure at the center of each grid box (mb)
!   humidity    : specific humidity (g/kg)
!
!-----------------------------------------------------------------------------

      subroutine CalcHeightSigmaLevel  &
     &  (temperature, humidity, za, press3e, surfPressure, &
     &   pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)

      implicit none
!     ----------------------
!     Argument declarations.
!     ----------------------

      logical, intent(in ) :: pr_diag
      integer, intent(in ) :: loc_proc
      integer, intent(in ) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      real*8 , intent(in ) :: press3e     (ilo:ihi, julo:jhi, k1-1:k2)
      real*8 , intent(in ) :: temperature (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in ) :: humidity    (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in ) :: surfPressure(ilo:ihi, julo:jhi)
      real*8 , intent(out) :: za(i1:i2, ju1:j2, k1:k2)

!     ----------------------
!     Local Variables
!     ----------------------

      real*8  :: zq(i1:i2, ju1:j2, k1-1:k2)
      integer :: ik

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'CalcHeightSigmaLevel called by ', loc_proc
      end if

!       --------------------------------------------------------------------
!       Compute the height of full and half-sigma levels above ground level.
!       --------------------------------------------------------------------

        zq(i1:i2, ju1:j2, k1-1) = 0.0d0

        do ik = k1, k2

          zq(i1:i2,ju1:j2,ik) =  &
     &      29.3d0 * temperature(i1:i2,ju1:j2,ik) *  &
     &      (1.0d0 + (0.61d-3 * humidity(i1:i2,ju1:j2,ik))) *  &
     &      (Log (surfPressure(i1:i2,ju1:j2) / press3e(i1:i2,ju1:j2,ik  )) -  &
     &       Log (surfPressure(i1:i2,ju1:j2) / press3e(i1:i2,ju1:j2,ik-1)) ) + &
     &      zq(i1:i2,ju1:j2,ik-1)

        end do

        do ik = k1, k2
          za(i1:i2,ju1:j2,ik) = 0.5d0 * (zq(i1:i2,ju1:j2,ik) +  &
     &                                   zq(i1:i2,ju1:j2,ik-1))
        end do

      return

      end subroutine CalcHeightSigmaLevel

!-----------------------------------------------------------------------------
      subroutine CalcConvection_bmass &
     &    (bmass, ps, ai, bi, pt, &
     &     i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)

      implicit none

#     include "gmi_phys_constants.h"

      integer, intent(in ) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      real*8 , intent(in ) :: pt
      real*8 , intent(in ) :: ai    (k1-1:k2)
      real*8 , intent(in ) :: bi    (k1-1:k2)
      real*8 , intent(in ) :: ps    (ilo:ihi, julo:jhi)
      real*8 , intent(out) :: bmass(i1:i2, ju1:j2, k1:k2)

      real*8  :: ap(k1:k2+1)
      real*8  :: bp(k1:k2+1)
      real*8  :: dapx, dbkx, fg
      integer :: ik, ij

      bmass(:,:,:) = 0.0d0

      do ik = k1, k2 + 1
        ap(ik) = ai(k2+1-ik)
        bp(ik) = bi(k2+1-ik)
      end do

      fg = 100.0d0 / GMI_G

      do ik = k1, k2

        dapx = fg * (ap(ik+1) - ap(ik)) * pt
        dbkx = fg * (bp(ik+1) - bp(ik))

        do ij = ju1, j2
          bmass(:,ij,ik) = dapx + (dbkx * ps(i1:i2,ij))
        end do

      end do

      return

      end subroutine CalcConvection_bmass
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
      subroutine CalcDryDepBoxHeight &
     &     (BoxHeightCenter, BoxHeightEdge, surfPressure, temperature, humidity, &
     &      press3c, press3e, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in)  :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      real*8 , intent(in)  :: surfPressure   (ilo:ihi, julo:jhi)
      real*8 , intent(in)  :: temperature    (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in ) :: press3c        (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in ) :: press3e        (ilo:ihi, julo:jhi, k1-1:k2)
      real*8 , intent(in)  :: humidity       (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(out) :: BoxHeightCenter(i1:i2,ju1:j2)
      real*8 , intent(out) :: BoxHeightEdge  (i1:i2,ju1:j2)


       BoxHeightCenter(:,:) =  &
     &      29.3d0 * temperature(i1:i2,ju1:j2,k1) *  &
     &      (1.0d0 + 0.61d-3*humidity(i1:i2,ju1:j2,k1)) *  &
     &      Log (surfPressure(i1:i2,ju1:j2) / press3c(i1:i2,ju1:j2,k1))

       BoxHeightEdge  (:,:) =  &
     &      29.3d0 * temperature(i1:i2,ju1:j2,k1) *  &
     &      (1.0d0 + 0.61d-3*humidity(i1:i2,ju1:j2,k1)) *  &
     &      Log (surfPressure(i1:i2,ju1:j2) / press3e(i1:i2,ju1:j2,k1))

      return

      end subroutine CalcDryDepBoxHeight
!-----------------------------------------------------------------------------

      end module GmiHeightSigmaLevel_mod
