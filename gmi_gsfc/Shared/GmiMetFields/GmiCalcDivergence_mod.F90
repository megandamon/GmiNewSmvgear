!------------------------------------------------------------------------------
! NASA GSFC _ SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !IMODULE: GmiCalcDivergence_mod
!
! !INTERFACE:
!
      module GmiCalcDivergence_mod
!
      implicit none
!
      private
      public  :: calcDivergence, Deter_Jrange
!
! !DESCRIPTION:
!
! !AUTHOR:
! Jules Kouatchou, NASA/GSFC, Jules.Kouatchou@nasa.gov
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcDivergence
!
! !INTERFACE:
!
      subroutine calcDivergence &
     &  (do_reduction, geofac_pc, geofac, dpi, xmass, ymass, &
     &   pr_diag, loc_proc, numLonDomains, j1p, j2p, i1_gl, i2_gl, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole)
!
      implicit none
!
! !INPUT VARIABLES:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: loc_proc, numLonDomains, commu_npole, commu_spole
      integer, intent(in) :: j1p, j2p
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
                          ! set to false if called on root processor
                          ! set to true  if called on worker processors
      logical, intent(in) :: do_reduction
                          ! special geometrical factor (geofac) for Polar cap
      real*8 , intent(in) :: geofac_pc
                          ! geometrical factor for meridional advection; geofac uses
                          ! correct spherical geometry, and replaces acosp as the
                          ! meridional geometrical factor in tpcore
      real*8 , intent(in) :: geofac(ju1_gl:j2_gl)
                          ! horizontal mass flux in E-W direction (mb)
      real*8 , intent(in) :: xmass (ilo:ihi, julo:jhi, k1:k2)
                          ! horizontal mass flux in N-S direction (mb)
      real*8 , intent(in) :: ymass (ilo:ihi, julo:jhi, k1:k2)
!
! !OUTPUT VARIABLES:
                          ! divergence at a grid point; used to calculate vertical motion (mb)
      real*8 , intent(out):: dpi(i1:i2, ju1:j2, k1:k2)
!
! !DESCRIPTION:
!  Calculates the divergence.
!
! !LOCAL VARIABLES:
      integer :: il, ij
      integer :: jst, jend
      real*8 :: dxmass(i1:i2, ju1:j2, k1:k2)
      real*8 :: dymass(i1:i2, ju1:j2, k1:k2)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
      if (pr_diag) then
        Write (6,*) 'calcDivergence called by ', loc_proc
      end if
!
      dxmass(:,:,:) = 0.0d0
      dymass(:,:,:) = 0.0d0
!
!
!test      call Deter_Jrange (ju1, j2, j1p, j2p, jst, jend, ju1, j2, ju1_gl, j2_gl)
!
!     -------------------------
!     Calculate N-S divergence.
!     -------------------------
!test      do ij = jst, jend
      do ij = ju1, j2
!
        dymass(:,ij,:) = (ymass(i1:i2,ij,:) - ymass(i1:i2,ij+1,:)) * geofac(ij)
!
      end do
!
!     -------------------------
!     Calculate E-W divergence.
!     -------------------------
!
!test      do ij = jst, jend
      do ij = ju1,j2
        do il = i1, i2
!
          dxmass(il,ij,:) = xmass(il,ij,:) - xmass(il+1,ij,:)
!
        end do
      end do
!
      dpi(:,:,:) = dymass(:,:,:) + dxmass(:,:,:)
!
!     ===========================
      call Do_Divergence_Pole_Sum  &
!     ===========================
     &  (do_reduction, geofac_pc, dpi, ymass, &
     &   numLonDomains, i1_gl, i2_gl, j1p, j2p, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole)
!
!       --------------------------------------------
!       Polar cap enlarged:  copy dpi to polar ring.
!       --------------------------------------------
!
      if (j1p /= ju1_gl+1) then
!
        if (ju1 == ju1_gl) then
!
          dpi(:,ju1+1,:) = dpi(:,ju1,:)
!
        end if
!
        if (j2 == j2_gl) then
!
          dpi(:,j2-1,:)  = dpi(:,j2,:)
!
        end if
!
      end if
!
      return
!
      end subroutine calcDivergence
!EOC
!------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Deter_Jrange
!
! DESCRIPTION
!   This routine determines the proper j (i.e., latitude) range to use.
!
! ARGUMENTS
!   jst_nopole  : proper first j if not at Pole
!   jend_nopole : proper last  j if not at Pole
!   jst_pole    : proper first j if     at Pole
!   jend_pole   : proper last  j if     at Pole
!   jst         : first j index to use
!   jend        : last  j index to use
!
!-----------------------------------------------------------------------------
!
      subroutine Deter_Jrange  &
     &  (jst_nopole, jend_nopole, jst_pole, jend_pole, jst, jend, &
     &   ju1, j2, ju1_gl, j2_gl)
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in)  :: ju1, j2, ju1_gl, j2_gl
      integer, intent(in)  :: jst_nopole, jend_nopole
      integer, intent(in)  :: jst_pole,   jend_pole
!
      integer, intent(out) :: jst, jend
!
!     ----------------
!     Begin execution.
!     ----------------
!
      if (ju1 == ju1_gl) then
!
        jst = jst_pole
!
      else
!
        jst = jst_nopole
!
      end if
!
!
      if (j2 == j2_gl) then
!
        jend = jend_pole
!
      else
!
        jend = jend_nopole
!
      end if
!
!
      return
!
      end subroutine Deter_Jrange
!
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_Divergence_Pole_Sum
!
! DESCRIPTION
!   This routine sets the divergence at the Poles.
!
! ARGUMENTS
!   do_reduction : set to true  if called on Master;
!                  set to false if called by Slaves
!   geofac_pc : special geometrical factor (geofac) for Polar cap
!   dpi   : divergence at a grid point; used to calculate vertical motion (mb)
!   ymass : horizontal mass flux in N-S direction (mb)
!
!-----------------------------------------------------------------------------
!
      subroutine Do_Divergence_Pole_Sum  &
     &  (do_reduction, geofac_pc, dpi, ymass, &
     &   numLonDomains, i1_gl, i2_gl, j1p, j2p, &
     &   ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2, &
     &   commu_npole, commu_spole)
!
      use GmiReduce_mod, only : Gmi_Sum_Pole_Reduce
!
      implicit none
!
!     ----------------------
!     Argument declarations.
!     ----------------------
!
      integer, intent(in) :: numLonDomains, commu_npole, commu_spole
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, j1p, j2p
      integer, intent(in) :: ilo, ihi, julo, jhi
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      logical :: do_reduction
      real*8, intent(in)  :: geofac_pc
      real*8  :: dpi  ( i1:i2,   ju1:j2,  k1:k2)
      real*8  :: ymass(ilo:ihi, julo:jhi, k1:k2)
!
!
!     ----------------------
!     Variable declarations.
!     ----------------------
!
      integer :: il, ik
!
      real*8  :: ri2
!
      real*8  :: mean_np(k1:k2)
      real*8  :: mean_sp(k1:k2)
      real*8  :: sumnp  (k1:k2)
      real*8  :: sumsp  (k1:k2)
!
!
!     ----------------
!     Begin execution.
!     ----------------
!
      ri2 = i2_gl
!
!     ==================
      if (ju1 == ju1_gl) then
!     ==================
!
        do ik = k1, k2
!
          sumsp(ik) = 0.0d0
!
          do il = i1, i2
!
            sumsp(ik) = sumsp(ik) + ymass(il,j1p,ik)
!
          end do
!
        end do
!
        if (do_reduction .and. (numLonDomains /= 1)) then
!         ===========================
          call Gmi_Sum_Pole_Reduce  &
!         ===========================
     &      (k1, k2, sumsp, ju1, j2, ju1_gl, j2_gl, &
     &      commu_npole, commu_spole)
        end if
!
        do ik = k1, k2
!
          mean_sp(ik) = -sumsp(ik) / ri2 * geofac_pc
!
          do il = i1, i2
!
            dpi(il,ju1,ik) = mean_sp(ik)
!
          end do
!
        end do
!
!     ======
      end if
!     ======
!
!     ================
      if (j2 == j2_gl) then
!     ================
!
        do ik = k1, k2
!
          sumnp(ik) = 0.0d0
!
          do il = i1, i2
!
            sumnp(ik) = sumnp(ik) + ymass(il,j2p+1,ik)
!
          end do
!
        end do
!
        if (do_reduction .and. (numLonDomains /= 1)) then
!         ===========================
          call Gmi_Sum_Pole_Reduce  &
!         ===========================
     &      (k1, k2, sumnp, ju1, j2, ju1_gl, j2_gl, &
     &      commu_npole, commu_spole)
        end if
!
        do ik = k1, k2
!
          mean_np(ik) = sumnp(ik) / ri2 * geofac_pc
!
          do il = i1, i2
!
            dpi(il,j2,ik) = mean_np(ik)
!
          end do
!
        end do
!
!     ======
      end if
!     ======
!
!
      return
!
      end subroutine Do_Divergence_Pole_Sum
!
!
      end module GmiCalcDivergence_mod
!
