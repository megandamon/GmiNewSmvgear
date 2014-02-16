!-------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:
!
      module GmiCloudVariables_mod
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:      
      private 
      public   :: CalcFractionalCloudCover
      public   :: CalcTotalCloudFraction
      public   :: CalcCloudOpticalDepth
      public   :: DiagCloudOpticalDepth
!
! !AUTHOR:
! Jules Kouatchou, Jules.Kouatchou-1@nasa.gov
!
!EOP
!-------------------------------------------------------------------------
      contains
!-----------------------------------------------------------------------------
!
! ROUTINE
!   CalcFractionalCloudCover
!
! DESCRIPTION
!   This routine calculates the cloud fraction.
!
! ARGUMENTS
!   max_cloud      : maximum overlap cloud fraction for LW
!   ran_cloud      : random  overlap cloud fraction for LW
!   fracCloudCover : fractional cloud cover
!-----------------------------------------------------------------------------

      subroutine CalcFractionalCloudCover  &
     &  (max_cloud, ran_cloud, fracCloudCover, i1, i2, ju1, j2, k1, k2)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in ) :: i1, i2, ju1, j2, k1, k2
      real*8 , intent(in ) :: max_cloud    (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in ) :: ran_cloud    (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(out) :: fracCloudCover(i1:i2, ju1:j2)

!EOP
!-------------------------------------------------------------------------
!BOC


      fracCloudCover(:,:) =  &
     &  1.0d0 -  &
     &  ((1.0d0 - Maxval (max_cloud(:,:,:), dim=3)) *  &
     &   (Product ((1.0d0 - ran_cloud(:,:,:)), dim=3)))

      return

      end subroutine CalcFractionalCloudCover
!EOC
!-------------------------------------------------------------------------
!BOP

!-----------------------------------------------------------------------------
!
! ROUTINE
!   CalcTotalCloudFraction
!
! DESCRIPTION
!   Compute total cloud fraction for combined large-scale and convective
!   clouds (GEOS-1 GCM descriptive documentation).
!
! ARGUMENTS
!   max_cloud      : maximum overlap cloud fraction for LW
!   ran_cloud      : random  overlap cloud fraction for LW
!   totalCloudFraction : fractional cloud cover
!-----------------------------------------------------------------------------

      subroutine CalcTotalCloudFraction  &
     &  (max_cloud, ran_cloud, totalCloudFraction, i1, i2, ju1, j2, k1, k2)

      implicit none

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer, intent(in ) :: i1, i2, ju1, j2, k1, k2
      real*8 , intent(in ) :: max_cloud    (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in ) :: ran_cloud    (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(out) :: totalCloudFraction(i1:i2, ju1:j2, k1:k2)

!EOP
!-------------------------------------------------------------------------
!BOC


      totalCloudFraction(:,:,:) =  &
     &  1.0d0 -  (1.0d0 - max_cloud(:,:,:)) * (1.0d0 - ran_cloud(:,:,:))

      return

      end subroutine CalcTotalCloudFraction
!EOC
!-------------------------------------------------------------------------
!BOP
!
! i!IROUTINE: CalcCloudOpticalDepth
!
! ARGUMENTS
!   max_cloud        : maximum overlap cloud fraction for LW
!   ran_cloud        : random  overlap cloud fraction for LW
!   tau_cloud        : cloud optical depth
!   temperature      : temperature (degK)
!   press3e          : atmospheric pressure at the edge of each grid box   (mb)
!   metdata_name_org : first  part of metdata_name, e.g., "NCAR"
!   surfPressure     : CTM surface pressure at t1+tdt (mb)
!   do_clear_sky     : 
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcCloudOpticalDepth
!
! !INTERFACE:
!
      subroutine CalcCloudOpticalDepth &
     &  (temperature, press3e, surfPressure, max_cloud, ran_cloud, tau_cloud, &
     &   i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      real*8 , intent(in) :: surfPressure (ilo:ihi, julo:jhi)
      real*8 , intent(in) :: temperature  (ilo:ihi, julo:jhi, k1:k2)
      real*8 , intent(in) :: press3e      (ilo:ihi, julo:jhi, k1-1:k2)
      real*8 , intent(in) :: max_cloud (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in) :: ran_cloud (i1:i2, ju1:j2, k1:k2)
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: tau_cloud (i1:i2, ju1:j2, k1:k2)
!
! !DESCRIPTION:
!    The met data did not have optical depth, so calculate it from
!    the cloud fractions, temperature and delta pressure.  This
!    algorithm comes from the GEOS-CHEM model and was given to us by
!    Bob Yantosca.
!
! !LOCAL VARIABLES:
      integer :: ij, il
      real*8  :: tauro(i1:i2, ju1:j2, k1:k2)

!EOP
!-------------------------------------------------------------------------
!BOC

      where (temperature(i1:i2,ju1:j2,:)  <= 190.66d0)  &
     &    tauro(:,:,:) = 0.0d0

      where ((temperature(i1:i2,ju1:j2,:) >  190.66d0) .and.  &
     &       (temperature(i1:i2,ju1:j2,:) <= 263.16d0))  &
     &   tauro(:,:,:) = (temperature(i1:i2,ju1:j2,:) - 190.66d0)**2 * 2.0d-6

      where ((temperature(i1:i2,ju1:j2,:) >  263.16d0) .and.  &
     &       (temperature(i1:i2,ju1:j2,:) <= 273.38d0))  &
     &   tauro(:,:,:) = (temperature(i1:i2,ju1:j2,:) * 6.95d-3) - 1.82d0

      where (temperature(i1:i2,ju1:j2,:)  >  273.38d0)  &
     &   tauro(:,:,:) = 0.08d0

      tau_cloud(:,:,:) =  (0.16d0       * max_cloud(:,:,:)) +  &
     &                    (tauro(:,:,:) * ran_cloud(:,:,:))

      do ij = ju1, j2
         do il = i1, i2
            tau_cloud(il,ij,k1:k2) =  tau_cloud(il,ij,k1:k2) *  &
                      (press3e(il,ij,k1-1:k2-1) - press3e(il,ij,k1:k2))
         end do
      end do

      return

      end subroutine CalcCloudOpticalDepth
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DiagCloudOpticalDepth
!
! !INTERFACE:
!
      subroutine DiagCloudOpticalDepth (max_cloud, ran_cloud, tau_cloud, &
                   OptDepth, i1, i2, ju1, j2, k1, k2, num_AerDust)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in)    :: i1, i2, ju1, j2, k1, k2, num_AerDust
      real*8 , intent(in)    :: max_cloud  (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in)    :: tau_cloud  (i1:i2, ju1:j2, k1:k2)
      real*8 , intent(in)    :: ran_cloud  (i1:i2, ju1:j2, k1:k2)
!
! !INPUT/OUTPUT PARAMETERS:
      REAL*8 , intent(inOut) :: OptDepth (i1:i2, ju1:j2, k1:k2, num_AerDust)

! !DESCRIPTION:
!   Computes the diagnostics for cloud optical depth.
!   It is used for aerosol/dust computations only.
!
!EOP
!-------------------------------------------------------------------------
!BOC

      !==============================================================
      ! OptDepth Diagnostic:
      !
      ! Tracer #1: Cloud optical depths
      ! Tracer #2: Max Overlap Cld Frac
      ! Tracer #3: Random Overlap Cld Frac
      !==============================================================

      !---------------------------------------------------
      ! OptDepth tracer #1: Cloud optical depths (1000 nm)
      !---------------------------------------------------
      OptDepth(:,:,:,1) = OptDepth(:,:,:,1) + tau_cloud(:,:,:)

      !---------------------------------------------------
      ! OptDepth tracer #2: Max Overlap Cld Frac
      !---------------------------------------------------
      OptDepth(:,:,:,2) = OptDepth(:,:,:,2) + max_cloud(:,:,:)

      !---------------------------------------------------
      ! OptDepth tracer #3: Random Overlap Cld Frac
      !---------------------------------------------------
      OptDepth(:,:,:,3) = OptDepth(:,:,:,3) + ran_cloud(:,:,:)

      return 
    
      end subroutine DiagCloudOpticalDepth
!EOC
!-------------------------------------------------------------------------

      end module GmiCloudVariables_mod
