!------------------------------------------------------------------------------
! NASA GSFC - SSSO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
      module fast_JX
!
! !USES:
      use GmiPrintError_mod, only : GmiPrintError
      use GmiASCIIoperations_mod, only : AsciiOpenRead
      use fastJX_Parameters_mod
      use fastJX_Bundle_mod, only : t_fastJXbundle

      IMPLICIT NONE

      private
      public :: Control_Fast_JX
      public :: InitializeFastJX
      public :: GetQAA_RAAinFastJX

#     include "GmiParameters.h"
!
! !DESCRIPTION:
! 
! !AUTHORS:
! Michael Prather, Oliver Wild, UCI
! Philip Cameron-Smith, LLNL
! Jules Kouatchou, Amidu Oloso, NASA GSFC
! 
! !HISTORY:
!   - Original code (adpated from Fast-J v3 distributed by Oliver Wild, UCI)
!     by Philip Cameron-Smith, LLNL
!   - Aug2004: A. Oloso
!     Code modifications based on Fast-JX v2 distributed by M. Prather
!   - 9March2005: Jules Kouatchou
!     The variable "ozone_ij" is now an optional arguments. When it is available,
!     the code uses the model climatology and the variable "do3" is set.
!     The variables "ODAER_ij" and "ODMDUST_ij" were added as arguments.
!
!
!EOP
!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
!!BOP
!
! !IROUTINE:   Control_Fast_JX
!
! DESCRIPTION
!   This routine acts as the interface between the GMI model and Fast-JX.
!
! ARGUMENTS
!
!   INPUTS
!
!-----------------------------------------------------------------------------

      subroutine Control_Fast_JX (JXbundle, k1, k2, chem_mask_khi, num_qjs,     &
                         month_gmi, jday, time_sec, fastj_offset_sec, &
                         londeg_i, latdeg_j, SZA_ij, press3c_ij,      &
                         pctm_ij, kel_ij, optdepth_ij, surf_alb_ij,   &
                         qjgmi_ij, overheadO3col_ij, ODAER_ij,        &
                         ODMDUST_ij, sflux_ij_1, sflux_ij_2,          &
                         sflux_ij_3, loc_proc, latyp_ij, ozone_ij)

      implicit none

# include "gmi_AerDust_const.h"
# include "setkin_par.h"
!# include "fast_JX_cmn_h.h"
!# include "fast_JX_cmn_t.h"
!# include "fast_JX_cmn_w.h"
# include "fast_JX_jv_cmn.h"
!
! !INPUT PARAMETERS:
      integer, intent(in) :: k1            ! first interface index [lpar = k2 - k1 + 1]
      integer, intent(in) :: k2            ! last  interface index [lpar = k2 - k1 + 1]
      integer, intent(in) :: chem_mask_khi ! number of chemistry levels     [jpnl]
      integer, intent(in) :: num_qjs       ! number of photolysis reactions [jppj]
      integer, intent(in) :: month_gmi     ! number of month (1- 12)        [month]
      integer, intent(in) :: jday          ! day    of year  (1-365)        [iday]
      integer, intent(in) :: latyp_ij      ! 
      integer, intent(in) :: loc_proc      ! local processor id
      real*8,  intent(in) :: time_sec      ! time of day in model (s, GMT)  [tau (hrs)]
      real*8,  intent(in) :: fastj_offset_sec ! offset from tau at which to do photolysis (s)
      real*8,  intent(in) :: londeg_i      ! longitude (midpoint, deg) [xgrd (rad)]
      real*8,  intent(in) :: latdeg_j      ! latitude  (midpoint, deg) [ygrd (rad), ydgrd (deg)]
      real*8,  intent(in) :: SZA_ij        ! solar zenith angle
      real*8,  intent(in) :: press3c_ij(k1:k2) ! pressure (mb)
      real*8,  intent(in) :: pctm_ij ! surface pressure (mb) [p]
      real*8,  intent(in) :: kel_ij(k1:k2) ! temperature at box centers (degK) [t]
      real*8,  intent(in) :: optdepth_ij(k1:k2) ! optical depth in box   (unitless) [od]
      real*8,  intent(in) :: surf_alb_ij        ! surface albedo     (fraction 0-1) [sa]
      real*8,  intent(in) :: ODAER_ij(k1:k2,NSADaer*nrh_b) ! optical depth for aerosol
      real*8,  intent(in) :: ODMDUST_ij(k1:k2,NSADdust)    ! optical depth for mineral dust
      real*8,  intent(in), optional ::  ozone_ij   (k1:k2) ! solar zenith angle
!
! !OUTPUT PARAMETERS:
      real*8, intent(out) :: overheadO3col_ij(k1:k2) ! overhead ozone column
      real*8, intent(out) :: qjgmi_ij(k1:chem_mask_khi, num_qjs) ! jvalues at grid centers (s^-1?) [zpj]
      real*8, intent(out) :: sflux_ij_1(k1:chem_mask_khi)     !nika
      real*8, intent(out) :: sflux_ij_2(k1:chem_mask_khi)     !nika
      real*8, intent(out) :: sflux_ij_3(k1:chem_mask_khi)     !nika
!
! !INPUT/OUTPUT PARAMETERS:
      TYPE(t_fastJXbundle), intent(inOut) :: JXbundle
!
!EOP
!-------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      character (len=75) :: err_msg

      logical, save :: first = .true.

      integer :: il, ik, j, k

      real*8  :: pi
      real*8  :: timej
      real*8  ::  amm(k1-1:k2)
      real*8  ::  bmm(k1-1:k2)

      real*8  :: zpj(chem_mask_khi-k1+1, num_qjs)
      real*8  :: flux_ij_gt_1(chem_mask_khi-k1+1)   !nika
      real*8  :: flux_ij_gt_2(chem_mask_khi-k1+1)   !nika
      real*8  :: flux_ij_gt_3(chem_mask_khi-k1+1)   !nika

      logical :: rootProc

!     -------------------------
!     Just do one column {PJC}.
!     -------------------------

      pi = 3.141592653589793d0  ! FLAG {PJC}

!     -----------------------------------------------------
!     Check that k1=1, otherwise Fast-JX may not work {PJC}.
!     -----------------------------------------------------

      if (k1 /= 1) then
         err_msg = &
           'Control_Fastj: Fast-J may not work correctly if k1 /= 1.'
         call GmiPrintError (err_msg, .true., 1, k1, 0, 0, 0.0d0, 0.0d0)
      end if

!     ----------------------------------------------------------------
!     Check dimensions against dimension maximums specified in include
!     files {PJC}.
!     ----------------------------------------------------------------

      if (JXbundle%lpar > LPAR_MAX) then
         err_msg = 'Control_Fast_JX: lpar exceeds LPAR_MAX. '//  &
                   'Increase LPAR_MAX and recompile.'
         call GmiPrintError(err_msg, .true., 2, JXbundle%lpar, LPAR_MAX, 0, 0.0d0, 0.0d0)
      end if

      if (JXbundle%jpnl > JPNL_MAX) then
         err_msg = 'Control_Fast_JX: jpnl exceeds JPNL_MAX. '//  &
                   'Increase JPNL_MAX and recompile.'
         call GmiPrintError(err_msg, .true., 2, JXbundle%jpnl, JPNL_MAX, 0, 0.0d0, 0.0d0)
      end if

      if (JXbundle%jppj > JPPJ_MAX) then
         err_msg = 'Control_Fast_JX: jppj exceeds JPPJ_MAX. '//  &
                   'Increase JPPJ_MAX and recompile.'
         call GmiPrintError(err_msg, .true., 2, JXbundle%jppj, JPPJ_MAX, 0, 0.0d0, 0.0d0)
      end if

      if (nb > NB_MAX) then
         err_msg = 'Control_Fast_JX: nb exceeds NB_MAX. '//  &
                   'Increase NB_MAX and recompile.'
         call GmiPrintError (err_msg, .true., 2, nb,   NB_MAX,   0, 0.0d0, 0.0d0)
      end if

      if (nc > NC_MAX) then
         err_msg = 'Control_Fast_JX: nc exceeds NC_MAX. '//  &
                   'Increase NC_MAX and recompile.'
         call GmiPrintError (err_msg, .true., 2, nc,   NC_MAX,   0, 0.0d0, 0.0d0)
      end if

!     -------------------------------------------------
!     Set up Fast-JX variables from GMI variables {PJC}.
!     -------------------------------------------------

      if ( present(ozone_ij) ) then
         JXbundle%do_model_clima = .true.
      else
         JXbundle%do_model_clima = .false.
      endif

      JXbundle%month = month_gmi
      JXbundle%iday  = jday
      JXbundle%tau   = time_sec         / 3600.0d0     ! convert to hours {PJC}
      timej          = fastj_offset_sec / 3600.0d0     ! convert to hours {PJC}

      JXbundle%xgrd (JXbundle%nslon) = londeg_i ! * pi / 180.0d0
      JXbundle%ygrd (JXbundle%nslat) = latdeg_j ! * pi / 180.0d0
      JXbundle%ydgrd(JXbundle%nslat) = latdeg_j

      pj(1)                    = pctm_ij
      pj(2:NB)                 = press3c_ij(k1:k2)
      pj(NB+1)                 = 0.0d0

      JXbundle%t (JXbundle%nslon, JXbundle%nslat, 1:JXbundle%lpar) = kel_ij     (k1:k2)

      ! Clouds
      JXbundle%od(JXbundle%nslon, JXbundle%nslat, 1:JXbundle%lpar) = optdepth_ij(k1:k2)

      JXbundle%lantyp(JXbundle%nslon, JXbundle%nslat) = latyp_ij

! beg bianfix Dec19, 2011
      !if (do_model_clima) then
      !   do3(1:lpar)              = ozone_ij(k1:k2)
      !endif
      if (present(ozone_ij)) then
         do3(1:JXbundle%lpar) = ozone_ij(k1:k2)
      endif
! end bianfix, Dec19,2011

      JXbundle%sa(JXbundle%nslon, JXbundle%nslat) = surf_alb_ij

      ! Aerosol OD profile [unitless]
      JXbundle%optaer(1:JXbundle%lpar,:)  = ODAER_ij(k1:k2,:)

      ! Mineral dust OD profile [unitless]
      JXbundle%optdust(1:JXbundle%lpar,:) = ODMDUST_ij(k1:k2,:)

      zpj(:,:) = -1.0d0

      sflux_ij_1(:)= 0.0d0
      sflux_ij_2(:)= 0.0d0
      sflux_ij_3(:)= 0.0d0

!     --------------------------------------------------------------------
!     call main Fast-JX routine, which actually calculates photolysis rates.
!     --------------------------------------------------------------------

      call Photo (JXbundle, zpj, timej, SZA_ij, &
                  flux_ij_gt_1, flux_ij_gt_2,flux_ij_gt_3, loc_proc, &
                  londeg_i,latdeg_j)

      qjgmi_ij(:,:) = zpj(:,:)

      sflux_ij_1(:)=flux_ij_gt_1(:)
      sflux_ij_2(:)=flux_ij_gt_2(:)
      sflux_ij_3(:)=flux_ij_gt_3(:)
      
!!!!
!!!! Overhead ozone column
!!!!
      overheadO3col_ij = 0.0d0
      overheadO3col_ij(k2) = DO3(k2)
      do il = k1, k2-1
         do ik = il+1, k2
            overheadO3col_ij(il) = overheadO3col_ij(il) + DO3(ik)
         end do
      end do


      return

      end subroutine Control_Fast_JX
!EOC
!------------------------------------------------------------------------------
!BOP
      subroutine InitializeFastJX (JXbundle, cross_section_file,rate_file, &
                           T_O3_climatology_file, num_qj_o3_to_2oh, &
                           num_qjs, chem_mask_khi, k2, k1, rootProc)

!# include "fast_JX_cmn_h.h"
!# include "fast_JX_cmn_t.h"
# include "fast_JX_jv_cmn.h"
!     cross_section_file    : X-Section quantum yield
!     rate_file             : Master rate file
!     T_O3_climatology_file : T & O3 climatology
      integer :: num_qj_o3_to_2oh, num_qjs, chem_mask_khi, k2, k1
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: cross_section_file
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: rate_file
      character (len=MAX_LENGTH_FILE_NAME), intent(in) :: T_O3_climatology_file
      logical, intent(in) :: rootProc
      TYPE(t_fastJXbundle), intent(inOut) :: JXbundle

!EOP
!------------------------------------------------------------------------------
!BOC
      JXbundle%nslat = 1
      JXbundle%nslon = 1

!     -------------------------
!     Set up Fast-JX dimensions.
!     -------------------------

      JXbundle%lpar = k2 - k1 + 1

      JXbundle%jpnl = chem_mask_khi
      JXbundle%jppj = num_qjs

      nb   = JXbundle%lpar + 1
      nc   = 2 * nb

!       ----------------------------------------
!       Initial call to Fast-JX to set things up.
!       ----------------------------------------

        call Inphot (JXbundle, cross_section_file,rate_file,  &
                     T_O3_climatology_file, num_qj_o3_to_2oh, rootProc)

      return 

      end subroutine InitializeFastJX
!------------------------------------------------------------------------------
!BOP

      subroutine GetQAA_RAAinFastJX (JXbundle, bRAA, bQAA)

!#     include "fast_JX_cmn_h.h"
#     include "fast_JX_jv_cmn.h"

      real*8 , intent(out) :: bRAA(:,:), bQAA(:,:)
      TYPE(t_fastJXbundle), intent(inOut) :: JXbundle

!EOP
!------------------------------------------------------------------------------
!BOC

      bRAA(:,:) = RAA(:,:)
      bQAA(:,:) = QAA(:,:)

      return

      end subroutine GetQAA_RAAinFastJX
!EOC
!------------------------------------------------------------------------------

# include "newcol.code"

!------------------------------------------------------------------------------
      end module fast_JX
