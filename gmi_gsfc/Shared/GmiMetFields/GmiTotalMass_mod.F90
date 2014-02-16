!-----------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!-----------------------------------------------------------------------------
!BOP
  module GmiTotalMass_mod

      implicit none
      
      private 
      public   :: calcTotalMass
!EOP
!-----------------------------------------------------------------------------
 contains

!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcTotalMass
!
! !INTERFACE:
!
      subroutine calcTotalMass  &
     &  (pt, ai, bi, psf, mcor, mass, &
     &   pr_diag, procID, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
! !USES:
      use GmiPrintError_mod, only : GmiPrintError

      implicit none

#     include "gmi_phys_constants.h"
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      ! pressure = (ai * pt) + (bi * psf) (mb)
      real*8 , intent(in) :: pt
      ! pressure = (ai * pt) + (bi * psf), ai at zone interface
      real*8 , intent(in) :: ai  (k1-1:k2)
      ! pressure = (ai * pt) + (bi * psf), bi at zone interface
      real*8 , intent(in) :: bi  (k1-1:k2)
      ! surface pressure field at t1, known at zone centers (mb)
      real*8 , intent(in) :: psf (ilo:ihi, julo:jhi)
      ! area of grid box (m^2)
      real*8 , intent(in) :: mcor(i1:i2, ju1:j2)
!
! !OUTPUT PARAMETERS:
      ! total mass of the atmosphere within each grid box   (kg)
      real*8 , intent(out) :: mass(i1:i2, ju1:j2, k1:k2)
!
! !DESCRIPTION:
!   This routine calculates the total mass of the atmosphere within
!   each grid box.  The calculation is based on the pressure difference
!   between the top and the bottom of the grid boxes.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      integer :: il, ij, ik
      integer :: ilx, ijx, ikx
      real*8  :: dp
      real*8  :: fac
!EOP
!------------------------------------------------------------------------------
!BOC

      if (pr_diag)  Write (6,*) 'calcTotalMass called by ', procID

!     --------------------------------------------------------------
!     Required mixing ratio change is converted to term proportional
!     to required mass change by multiplying by dx * dy * dp.
!     --------------------------------------------------------------

      ilx = -999

      fac = 100.0d0 / GMI_G

      do ik = k1, k2
        do ij = ju1, j2
          do il = i1, i2

            dp = (ai(ik) - ai(ik-1)) * pt +  &
     &           (bi(ik) - bi(ik-1)) * psf(il,ij)

            mass(il,ij,ik) = Abs (dp) * mcor(il,ij) * fac

            if (mass(il,ij,ik) <= 0.0d0) then
              ilx = il
              ijx = ij
              ikx = ik
            end if

          end do
        end do
      end do

      if (ilx /= -999) then
        err_msg = 'Bad mass in calcTotalMass.'
        call GmiPrintError (err_msg, .false., 2, ilx, ijx, 0, 0.0d0, 0.0d0)
        call GmiPrintError (err_msg, .true., 1, ikx, 0, 1, mass(ilx,ijx,ikx), 0.0d0)
      end if

      return

      end subroutine calcTotalMass
!EOC
!------------------------------------------------------------------------------

      end module GmiTotalMass_mod
