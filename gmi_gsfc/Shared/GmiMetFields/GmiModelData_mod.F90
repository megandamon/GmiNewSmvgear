      module GmiModelData_mod

      use GmiPrintError_mod, only : GmiPrintError

      implicit none

      private
      public   :: Set_Vert_Dao, Calc_Layer_Dsigma

#     include "gmi_phys_constants.h"
!#     include "gmi_grid_constants.h"

!=============================================================================
      contains
!=============================================================================


!-----------------------------------------------------------------------------
!BOP
! 
! !IROUTINE:
!
! !INTERFACE:
!
      subroutine Set_Vert_Dao  &
     &  (metdata_name, ktr, ptop, pt, ai, bi, am, bm, k1, k2)

      implicit none

#     include "gmi_vert_data.h"

! !INPUT PARAMETERS:
      integer, intent(in) :: k1, k2
      ! met data netcdf file attribute "Met_Dat_Name", e.g., "NCAR_MATCH_4x5x52"
      character (len=50) :: metdata_name
!
! !OUPUT PARAMETERS:
      ! number of vertical tropospheric  levels
      integer, intent(out) :: ktr
      ! pressure at top edge of top zone (mb); (i.e., pressure at top of model grid)
      real*8 , intent(out) :: ptop
      ! pressure = (ai * pt) + (bi * psf) (mb)
      real*8 , intent(out) :: pt
      ! pressure = (ai * pt) + (bi * psf), ai at zone interface
      real*8 , intent(out) :: ai (k1-1:k2)
      ! pressure = (ai * pt) + (bi * psf), bi at zone interface
      real*8 , intent(out) :: bi (k1-1:k2)
      ! pressure = (am * pt) + (bm * psf), am at zone midpoint
      real*8 , intent(out) :: am (k1:k2)
      ! pressure = (am * pt) + (bm * psf), bm at zone midpoint
      real*8 , intent(out) :: bm (k1:k2)
!
! !DESCRIPTION:
!   This routine defines the DAO model vertical levels.
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg
      integer k

!EOP
!------------------------------------------------------------------------------
!BOC
!     =================================
      select case (Trim (metdata_name))
!     =================================

!       =========================
        case ('EC_OSLO_2x2x40')
!       =========================

          ktr   = KTR_DF40

          ptop  = PTOP_DF40
          pt    = PT_DF40

          ai(:) = ai_data_df40(:)
          bi(:) = bi_data_df40(:)

          am(:) = am_data_df40(:)
          bm(:) = bm_data_df40(:)

!       =========================
        case ('EC_OSLO_2x2x37')
!       =========================

          ktr   = KTR_DF37

          ptop  = PTOP_DF37
          pt    = PT_DF37

          ai(:) = ai_data_df37(:)
          bi(:) = bi_data_df37(:)

          am(:) = am_data_df37(:)
          bm(:) = bm_data_df37(:)

!       =========================
        case ('DAO_FVCCM_4x5x28',  &
     &        'DAO_FVDAS_4x5x28',  &
     &        'DAO_FVCCM_2.5x2x28',  &
     &        'GMAO_GCM_2x2%5x28')
!       =========================

          ktr   = KTR_DF28

          ptop  = PTOP_DF28
          pt    = PT_DF28

          ai(:) = ai_data_df28(:)
          bi(:) = bi_data_df28(:)

          am(:) = am_data_df28(:)
          bm(:) = bm_data_df28(:)

!       ===========================
        case ('GMAO_GCM_2x2%5x33')
!       ===========================

          ktr   = KTR_DF33

          ptop  = PTOP_DF33
          pt    = PT_DF33

          ai(:) = ai_data_df33(:)
          bi(:) = bi_data_df33(:)

          am(:) = am_data_df33(:)
          bm(:) = bm_data_df33(:)

!       ===========================
        case ('GMAO_GEOS4GCM_2%5x2x42',  &
     &        'GMAO_GEOS4GCM_5x4x42', 'GMAO_GEOS4DAS_5x4x42')
!       ===========================

          ktr   = KTR_DF42

          ptop  = PTOP_DF42
          pt    = PT_DF42

          ai(:) = ai_data_df42(:)
          bi(:) = bi_data_df42(:)

          am(:) = am_data_df42(:)
          bm(:) = bm_data_df42(:)

!       ===========================
        case ('GMAO_GEOS5DAS_4x5x72','GMAO_GEOS5DAS_2x2%5x72', &
     &        'GMAO_GEOS5DAS_1x1%25x72','GMAO_GEOS5DAS_0.5x0.67x72', &
     &        'GMAO_GEOS5MERRA300_2x2%5x72', 'GMAO_GEOS5CCM_2x2%5x72', &
     &        'GMAO_GEOS5MERRA_2x2%5x72', 'GMAO_GEOS5MERRA300_1x1%25x72', &
     &        'GMAO_GEOS5MERRA_1x1%25x72')
!       ===========================

          ktr   = KTR_GD72

          ptop  = PTOP_GD72
          pt    = PT_GD72

          do k=0,72
            ai(k) = ai_data_gd72(72-k)
            bi(k) = bi_data_gd72(72-k)
           enddo
          do k=1,72
            am(k) = am_data_gd72(72-k+1)
            bm(k) = bm_data_gd72(72-k+1)
           enddo

!       ===========================
!       case ('DAO_FVDAS_2x2.5x55')
!       This is for vertical degradation of GEOS5 to the GEOS4 vertical grid
        case ('DAO_FVDAS_2x2.5x55', 'GMAO_GEOS5DAS_2x2%5x55', &
    &         'GMAO_GEOS5MERRA_2x2%5x55', 'GMAO_GEOS5CCM_2x2%5x55', &
    &        'GMAO_GEOS5MERRA300_2x2%5x55')
!       ===========================

          ktr   = KTR_DF55

          ptop  = PTOP_DF55
          pt    = PT_DF55

          ai(:) = ai_data_df55(:)
          bi(:) = bi_data_df55(:)

          am(:) = am_data_df55(:)
          bm(:) = bm_data_df55(:)

!       ========================
        case ('DAO_GS_2x2.5x20',  &
     &        'DAO_GS_4x5x20')
!       ========================

          ktr   = KTR_DG20

          ptop  = PTOP_DG20
          pt    = PT_DG20

          ai(:) = ai_data_dg20(:)
          bi(:) = bi_data_dg20(:)

          am(:) = am_data_dg20(:)
          bm(:) = bm_data_dg20(:)

!       ===========================
        case ('DAO_GEOS3_1x1x48',  &
     &        'DAO_GEOS3_2x2.5x48')
!       ===========================

          ktr   = KTR_DG48

          ptop  = PTOP_DG48
          pt    = PT_DG48

          ai(:) = ai_data_dg48(:)
          bi(:) = bi_data_dg48(:)

          am(:) = am_data_dg48(:)
          bm(:) = bm_data_dg48(:)

!       ========================
        case ('DAO_GS_2x2.5x26',  &
     &        'DAO_GS_4x5x26')
!       ========================

          ktr   = KTR_DS26

          ptop  = PTOP_DS26
          pt    = PT_DS26

          ai(:) = ai_data_ds26(:)
          bi(:) = bi_data_ds26(:)

          am(:) = am_data_ds26(:)
          bm(:) = bm_data_ds26(:)

!       ========================
        case ('DAO_GS_2x2.5x29',  &
     &        'DAO_GS_4x5x29')
!       ========================

          ktr   = KTR_DS29

          ptop  = PTOP_DS29
          pt    = PT_DS29

          ai(:) = ai_data_ds29(:)
          bi(:) = bi_data_ds29(:)

          am(:) = am_data_ds29(:)
          bm(:) = bm_data_ds29(:)

!       ========================
        case ('DAO_GS_2x2.5x46',  &
     &        'DAO_GS_4x5x46')
!       ========================

          ktr   = KTR_DS46

          ptop  = PTOP_DS46
          pt    = PT_DS46

          ai(:) = ai_data_ds46(:)
          bi(:) = bi_data_ds46(:)

          am(:) = am_data_ds46(:)
          bm(:) = bm_data_ds46(:)

!       ========================
        case ('DAO_GU_2x2.5x25')
!       ========================

          ktr   = KTR_DU25

          ptop  = PTOP_DU25
          pt    = PT_DU25

          ai(:) = ai_data_du25(:)
          bi(:) = bi_data_du25(:)

          am(:) = am_data_du25(:)
          bm(:) = bm_data_du25(:)

!       ===========================
        case ('GISS_2prime_4x5x23')
!       ===========================

          ktr   = KTR_G2P23

          ptop  = PTOP_G2P23
          pt    = PT_G2P23

          ai(:) = ai_data_g2p23(:)
          bi(:) = bi_data_g2p23(:)

          am(:) = am_data_g2p23(:)
          bm(:) = bm_data_g2p23(:)

!       ===========================
        case ('GISS_2prime_4x5x46')
!       ===========================

          ktr   = KTR_G2P46

          ptop  = PTOP_G2P46
          pt    = PT_G2P46

          ai(:) = ai_data_g2p46(:)
          bi(:) = bi_data_g2p46(:)

          am(:) = am_data_g2p46(:)
          bm(:) = bm_data_g2p46(:)

!       =========================
        case ('NCAR_CCM3_4x5x18',  &
     &        'NCAR_CCM3_T42x18')
!       =========================

          ktr   = KTR_NC18

          ptop  = PTOP_NC18
          pt    = PT_NC18

          ai(:) = ai_data_nc18(:)
          bi(:) = bi_data_nc18(:)

          am(:) = am_data_nc18(:)
          bm(:) = bm_data_nc18(:)

!       =========================
        case ('NCAR_CCM2_4x5x44')
!       =========================

          ktr   = KTR_NC44

          ptop  = PTOP_NC44
          pt    = PT_NC44

          ai(:) = ai_data_nc44(:)
          bi(:) = bi_data_nc44(:)

          am(:) = am_data_nc44(:)
          bm(:) = bm_data_nc44(:)

!       ==========================
        case ('NCAR_MATCH_4x5x52')
!       ==========================

          ktr   = KTR_NM52

          ptop  = PTOP_NM52
          pt    = PT_NM52

          ai(:) = ai_data_nm52(:)
          bi(:) = bi_data_nm52(:)

          am(:) = am_data_nm52(:)
          bm(:) = bm_data_nm52(:)

!       ============
        case default
!       ============

          err_msg = 'Set_Vert_Dao problems:  ' // metdata_name
          call GmiPrintError (err_msg, .true., 1, k2, 0, 0, 0.0d0, 0.0d0)

!     ==========
      end select
!     ==========

      return

      end subroutine Set_Vert_Dao
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Layer_Dsigma
!
! !INTERFACE:
!
      subroutine Calc_Layer_Dsigma  (pt, ai, bi, dap, dbk, k1, k2)

      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: k1, k2
      ! pressure = (ai * pt) + (bi * psf) (mb)
      real*8 , intent(in) :: pt
      ! pressure = (ai * pt) + (bi * psf), ai at zone interface
      real*8 , intent(in) :: ai(k1-1:k2)
      ! pressure = (ai * pt) + (bi * psf), bi at zone interface
      real*8 , intent(in) :: bi(k1-1:k2)
!
! !OUTPUT PARAMETERS:
      ! pressure difference across layer from log pressure term
      real*8, intent(out) :: dap(k1:k2)
      ! pressure difference across layer from sigma        term
      real*8, intent(out) :: dbk(k1:k2)
!
! !DESCRIPTION:
! This routine calculates the change in hybrid pressure factors.
!
! !LOCAL VARIABLES:
      integer :: ik
!EOP
!-----------------------------------------------------------------------------
!BOC
      do ik = k1, k2
         dap(ik) = (ai(ik-1) - ai(ik)) * pt
         dbk(ik) =  bi(ik-1) - bi(ik)
      end do

      return

      end subroutine Calc_Layer_Dsigma
!EOC
!-----------------------------------------------------------------------------
      end module GmiModelData_mod
