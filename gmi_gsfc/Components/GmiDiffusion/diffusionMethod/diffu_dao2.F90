
!=============================================================================
!
! $Id: diffu_dao2.F90,v 1.9 2011-08-09 22:12:58 mrdamon Exp $
!
! CODE DEVELOPER
!   Dan Bergmann,   LLNL (Original code from John Walton, UMICH)
!   dbergmann@llnl.gov
!
! FILE
!   diffu_dao2.F
!
! ROUTINES
!   Do_Vert_Diffu
!
!=============================================================================


!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Vert_Diffu
!
! DESCRIPTION
!   This routine calculates vertical diffusion using an implicit technique.
!
!   Note:  Constituent mass is conserved during mixing.
!          Constant field should remain constant.
!
! ARGUMENTS
!   tdt     : model time step (s)
!   vert_diffu_coef : scalar vertical diffusion coefficient (m^2/s)
!   pbl     : planetary boundary layer depth (m)
!   kzz     : array of vertical diffusion coefficients (m^2/s)
!   press3c : atmospheric pressure at the center of each grid box (mb)
!   press3e : atmospheric pressure at the edge   of each grid box (mb)
!   psf     : surface pressure field at t1, known at zone centers (mb)
!   qq      : species concentration (mixing ratio)
!
!-----------------------------------------------------------------------------

      subroutine Do_Vert_Diffu  &
     &  (tdt, vert_diffu_coef, pbl, kzz, press3c, press3e, psf, concentration, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi, num_species)

      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      use GmiSpcConcentrationMethod_mod, only : isFixedConcentration

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in   ) :: i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi
      integer, intent(in   ) :: num_species
      real*8 , intent(in   ) :: tdt
      real*8 , intent(in   ) :: vert_diffu_coef
      real*8 , intent(in   ) :: pbl(i1:i2, ju1:j2)
      real*8 , intent(in   ) :: kzz(i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in   ) :: press3c(ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in   ) :: press3e(ilo:ihi, julo:jhi, k1-1:k2)
      real*8 , intent(in   ) :: psf(ilo:ihi, julo:jhi)
      type (t_GmiArrayBundle), intent(inout) :: concentration(num_species)


!     -----------------------
!     Parameter declarations.
!     -----------------------

      real*8, parameter :: HZERO  = 8000.0d0  ! atmospheric scale height (m)

      real*8, parameter :: HZERO2 = HZERO * HZERO  ! (m^2)

!     ----------------------
!     Variable declarations.
!     ----------------------

      integer :: ik, ic

      real*8  :: avdiff(i1:i2, ju1:j2, k1:k2-1)
      real*8  :: bvdiff(i1:i2, ju1:j2, k1:k2-1)
      real*8  :: cvdiff(i1:i2, ju1:j2, k1:k2-1)
      real*8  :: vdiff (i1:i2, ju1:j2, k1:k2-1)

      real*8  :: bb    (i1:i2, ju1:j2, k1:k2-1)
      real*8  :: cc    (i1:i2, ju1:j2, k1:k2-1)
      real*8  :: rr    (i1:i2, ju1:j2, k1:k2-1)


!     ----------------
!     Begin execution.
!     ----------------

      avdiff(:,:,:) = 0.0d0; bvdiff(:,:,:) = 0.0d0
      cvdiff(:,:,:) = 0.0d0; vdiff (:,:,:) = 0.0d0

      bb(:,:,:) = 0.0d0; cc(:,:,:) = 0.0d0
      rr(:,:,:) = 0.0d0



!     --------------------------------------------------------
!     If there was no kzz in the met file, it was set to < 0.
!     In this case calculate a vertical diffusion coefficient,
!     otherwise use kzz.
!     --------------------------------------------------------

      if (Maxval (kzz) < 0.0d0) then

!       ------------------------------------------------------
!       Set up the vertical diffusion coefficients.
!       First put height (m) into vdiff, then put simple
!       parameterization for diffusion coefficient into vdiff,
!       then put it into units of per mb
!       ------------------------------------------------------

        vdiff = 0.0d0

        do ik = 1, 10

          vdiff(:,:,ik) =  &
     &      -HZERO *  &
     &      Log (press3e(i1:i2,ju1:j2,ik) / psf(i1:i2,ju1:j2))

          where (vdiff(:,:,ik) < pbl(:,:))

            vdiff(:,:,ik) =  &
     &        vert_diffu_coef *  &
     &          (1.0d0 -  &
     &           2.0d0 *  &
     &           (Abs (vdiff(:,:,ik) - (pbl(:,:) * 0.5d0)) / pbl(:,:)))

          elsewhere

            vdiff(:,:,ik) = 0.0d0

          end where

          vdiff(:,:,ik) =  &
     &      vdiff(:,:,ik) /  &
!c   &      HZERO2 * (press3e(i1:i2,ju1:j2,ik)**2) /
     &      HZERO2 *  &
     &      (press3e(i1:i2,ju1:j2,ik) * press3e(i1:i2,ju1:j2,ik)) /  &
     &      (press3c(i1:i2,ju1:j2,ik) - press3c(i1:i2,ju1:j2,ik+1))

        end do

      else

          vdiff(:,:,k1:k2-1) =  &
     &      kzz(:,:,k1:k2-1) /  &
!c   &      HZERO2 * (press3e(i1:i2,ju1:j2,k1:k2-1)**2) /
     &      HZERO2 *  &
     &      (press3e(i1:i2,ju1:j2,k1:k2-1) *  &
     &       press3e(i1:i2,ju1:j2,k1:k2-1)) /  &
     &      (press3c(i1:i2,ju1:j2,k1:k2-1) -  &
     &       press3c(i1:i2,ju1:j2,k1+1:k2))

      end if


!     -----------------------------------------------------------
!     Construct the tri-diagonal matrix which will be used in the
!     implicit vertical diffusion solution.
!       avdiff is the lower band
!       bvdiff is the main diagonal
!       cvdiff is the upper band
!     -----------------------------------------------------------

!     -----------------------
!     First do the top layer.
!     -----------------------

      avdiff(:,:,k2-1) = 0.0d0
      bvdiff(:,:,k2-1) =  vdiff(:,:,k2-2)
      cvdiff(:,:,k2-1) = -vdiff(:,:,k2-2)


!     ----------------------
!     Now do all the layers.
!     ----------------------

      avdiff(:,:,k1+1:k2-2) = -vdiff(:,:,k1+1:k2-2)
      bvdiff(:,:,k1+1:k2-2) = (vdiff(:,:,k1+1:k2-2) +  &
     &                         vdiff(:,:,k1  :k2-3))
      cvdiff(:,:,k1+1:k2-2) = -vdiff(:,:,k1  :k2-3)


!     ------------------------------------
!     Now do the layer next to the ground.
!     ------------------------------------

      avdiff(:,:,k1) = -vdiff(:,:,k1)
      bvdiff(:,:,k1) =  vdiff(:,:,k1)
      cvdiff(:,:,k1) = 0.0d0


!     ---------------------------------------------------
!     Now convert them to the proper dimensionless units.
!     ---------------------------------------------------

      avdiff(i1:i2,ju1:j2,k1:k2-1) =  &
     &  avdiff(i1:i2,ju1:j2,k1:k2-1) * tdt /  &
     &  (press3e(i1:i2,ju1:j2,k1-1:k2-2) -  &
     &   press3e(i1:i2,ju1:j2,k1:k2-1))

      bvdiff(i1:i2,ju1:j2,k1:k2-1) =  &
     &  bvdiff(i1:i2,ju1:j2,k1:k2-1) * tdt /  &
     &  (press3e(i1:i2,ju1:j2,k1-1:k2-2) -  &
     &   press3e(i1:i2,ju1:j2,k1:k2-1))  +  &
     &  1.0d0

      cvdiff(i1:i2,ju1:j2,k1:k2-1) =  &
     &  cvdiff(i1:i2,ju1:j2,k1:k2-1) * tdt  /  &
     &  (press3e(i1:i2,ju1:j2,k1-1:k2-2) -  &
     &   press3e(i1:i2,ju1:j2,k1:k2-1))


!     ---------------------
!     Now solve the system.
!     ---------------------

      icloop: do ic = 1, num_species

!       ====================================
        if (isFixedConcentration(ic)) cycle icloop
!       ====================================

!       -----------------------------------------
!       Set up the right hand side of the system.
!       -----------------------------------------

        bb(:,:,k2-1) = bvdiff(:,:,k2-1)
        rr(:,:,k2-1) = concentration(ic)%pArray3D(:,:,k2-1)

!       ---------------------------------
!       Eliminate the lower diagonal (a).
!       ---------------------------------

        do ik = k2 - 1, k1 + 1, -1

          cc(:,:,ik) = cvdiff(:,:,ik) / bb(:,:,ik)
          rr(:,:,ik) =     rr(:,:,ik) / bb(:,:,ik)


          bb(:,:,ik-1) = bvdiff(:,:,ik-1) -  &
     &                   avdiff(:,:,ik-1) * cc(:,:,ik)
          rr(:,:,ik-1) = concentration(ic)%pArray3D(:,:,ik-1) -  &
     &                   avdiff(:,:,ik-1) * rr(:,:,ik)

        end do

!       -------------------------------
!       Solve for the new mixing ratio.
!       -------------------------------

        concentration(ic)%pArray3D(:,:,k1) = rr(:,:,k1) / bb(:,:,k1)

        do ik = k1 + 1, k2 - 1

          concentration(ic)%pArray3D(:,:,ik) =  &
     &      rr(:,:,ik) -  &
     &      concentration(ic)%pArray3D(:,:,ik-1) * cc(:,:,ik)

        end do

      end do icloop

      return

      end subroutine Do_Vert_Diffu

