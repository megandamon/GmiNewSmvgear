!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiMassFluxes
!
! !INTERFACE:
!
      module GmiMassFluxes_mod
!
! !USES:
      use ESMF_Mod
!
     implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
     private
!.not used..     public  :: calcHorizontalMassFlux
     public  :: calcVerticalMassFlux
     public  :: convertMassFlux, calcDelta_xyt
!
#     include "gmi_phys_constants.h"
!
      real*8             , save :: dtdy      ! dt/dy      (s/m)
      real*8             , save :: dtdy5     ! 0.5 * dtdy (s/m)
      real*8             , save :: deltay    ! dy         (m)
      real*8, allocatable, save :: dtdx (:)  ! dt/dx      (s/m)
      real*8, allocatable, save :: dtdx5(:)  ! 0.5 * dtdx (s/m)
      real*8, allocatable, save :: deltax(:) ! dx         (m)
!
!
! !DESCRIPTION:
!  Routines for calculating the horizontal and vertical mass fluxes.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertMassFlux
!
! !INTERFACE:
!
      subroutine convertMassFlux (xMassFlux, yMassFlux, zMassFlux, MX, MY, MZ, &
     &                  tdt, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: ilo, ihi, julo, jhi
      real*8 , intent(in) :: tdt              ! model time step (s)
      real*8 , intent(in) :: xMassFlux(ilo:ihi,julo:jhi,k1:k2) ! horizontal mass flux in E-W
                                              ! direction (mb)
      real*8 , intent(in) :: yMassFlux(ilo:ihi,julo:jhi,k1:k2) ! horizontal mass flux in N-S
                                              ! direction (mb)
      real*8 , intent(in) :: zMassFlux(i1:i2,ju1:j2,k1:k2) ! vertical mass flux (mb)
!
! !OUTPUT PARAMETERS:
      real(ESMF_KIND_R4), intent(out) :: MX(i1:i2,ju1:j2,k1:k2) ! eastward pressure flux
                                                   ! on Cgrid (m+2 s-1 Pa)
      real(ESMF_KIND_R4), intent(out) :: MY(i1:i2,ju1:j2,k1:k2) ! northward pressure flux
                                                   ! on Cgrid (m+2 s-1 Pa)
      real(ESMF_KIND_R4), intent(out) :: MZ(i1:i2,ju1:j2,k1:k2) ! upward mass flux due to
                                                   ! dynamics (kg m-2 s-1)
!
! !DESCRIPTION:
!  This routine converts the mass fluxes calculated in the GMI code into
!  mass flux quantities acceptable by advCode.
!   \begin{eqnarray*}
!   \begin{array}{rcl}
!   MX & = & 100 \times \frac{\Delta x \Delta y}{\Delta t} xmassFlux \\
!   MY & = & 100 \times \frac{\Delta x \Delta y}{\Delta t} ymassFlux \\
!   MZ & = & 100 \times \frac{1}{g\Delta t} zMassFlux
!   \end{array}
!   \end{eqnarray*}
!
! !DEFINED PARAMETERS:
      real(ESMF_KIND_R4), parameter :: mbToPascal = 100.0
! !LOCAL VARIABLES:
      integer :: ij, ic, jc
      real*8  :: dxyt
!
!EOP
!-------------------------------------------------------------------------
!BOC
      do ij = ju1, j2
         dxyt = deltax(ij)*deltay/tdt
         MX(i1:i2,ij,k1:k2) = mbToPascal*dxyt*xMassFlux(i1:i2,ij,k2:k1:-1)
         MY(i1:i2,ij,k1:k2) = mbToPascal*dxyt*yMassFlux(i1:i2,ij,k2:k1:-1)
      end do
!
      MZ(i1:i2,ju1:j2,k1:k2) = mbToPascal*(1.0/(GMI_G*tdt))*  &
     &                       zMassFlux(i1:i2,ju1:j2,k2:k1:-1)
!
      return
!
      end subroutine convertMassFlux
!EOC
!------------------------------------------------------------------------------
!BOP
!.not used..!
!.not used..! !IROUTINE: calcHorizontalMassFlux
!.not used..!
!.not used..! !INTERFACE:
!.not used..!
!.not used..      subroutine calcHorizontalMassFlux &
!.not used..     &    (met_grid_type, tdt, cose, cosp, delpm,  &
!.not used..     &     uu, vv, xMassFlux, yMassFlux, &
!.not used..     &     pr_diag, procID, gmi_nborder, i2_gl, ju1_gl, j2_gl, &
!.not used..     &     ilo, ihi, julo, jvlo, jhi, i1, i2, ju1, j2, k1, k2)
!.not used..!
!.not used..! !USES:
!.not used..      use GmiWrapMaster_mod, only : wrapMaster_3du
!.not used..!
!.not used..      implicit none
!.not used..!
!.not used..! !INPUT PARAMETERS:
!.not used..      logical         , intent(in) :: pr_diag
!.not used..      integer         , intent(in) :: procID
!.not used..      integer         , intent(in) :: gmi_nborder
!.not used..      integer         , intent(in) :: i2_gl, ju1_gl, j2_gl
!.not used..      integer         , intent(in) :: ilo, ihi, julo, jvlo, jhi
!.not used..      integer         , intent(in) :: i1, i2, ju1, j2, k1, k2
!.not used..      character(len=1), intent(in) :: met_grid_type ! met grid type, 'A' or 'C'
!.not used..      real*8          , intent(in) :: tdt ! model time step (s)
!.not used..      real*8          , intent(in) :: cose (ju1_gl:j2_gl+1) ! cosine of grid box edges
!.not used..      real*8          , intent(in) :: cosp (ju1_gl:j2_gl) ! cosine of grid box centers
!.not used..      ! pressure thickness, the psudo-density in a hydrostatic system
!.not used..      ! at t1+tdt/2 (approximate) (mb)
!.not used..      real*8           , intent(in) :: delpm(ilo:ihi, julo:jhi, k1:k2)
!.not used..      ! wind velocity  in E-W direction at t1+tdt/2 (m/s)
!.not used..      real*8           , intent(in) :: uu   (ilo:ihi, julo:jhi, k1:k2)
!.not used..      ! wind velocity  in N-S direction at t1+tdt/2 (m/s)
!.not used..      real*8           , intent(in) :: vv   (ilo:ihi, jvlo:jhi, k1:k2)
!.not used..!
!.not used..! !OUTPUT PARAMETERS:
!.not used..      ! horizontal mass flux in E-W direction (mb)
!.not used..      real*8, intent(out) :: xMassFlux(ilo:ihi, julo:jhi, k1:k2)
!.not used..      ! horizontal mass flux in N-S direction (mb)
!.not used..      real*8, intent(out) :: yMassFlux(ilo:ihi, julo:jhi, k1:k2)
!.not used..!
!.not used..! !DESCRIPTION:
!.not used..!  This routine converts winds on A or C grid to mass flux on the C grid.
!.not used..!
!.not used..! !DEFINED PARAMETERS:
!.not used..!     ----------------------------------------------------------------------
!.not used..!     INTERP_MASS_FLUX_COMPONENTS :
!.not used..!
!.not used..!        If true, then interpolate winds and delta pressure separately from
!.not used..!        A-grid to C-grid before multiplying to get mass flux on the C-grid.
!.not used..!
!.not used..!        If false, then multiply winds and delta pressure to get mass flux
!.not used..!        on the A-grid, then interpolate mass flux on to the C-grid.
!.not used..!     ----------------------------------------------------------------------
!.not used..
!.not used..      logical, parameter :: INTERP_MASS_FLUX_COMPONENTS = .TRUE.
!.not used..!
!.not used..! !LOCAL VARIABLES:
!.not used..      logical, save :: first = .true.
!.not used..      integer :: il, ij, ik
!.not used..      real*8  :: tmp_pu1 ! layer thickness on eastern  edge of gridbox (mbar)
!.not used..      real*8  :: tmp_pu2 ! layer thickness on southern edge of gridbox (mbar)
!.not used..      real*8  :: tmp_crx ! zonal Courant number on eastern edge of gridbox
!.not used..      real*8  :: tmp_cry ! meridional Courant number on southern edge of gridbox
!.not used..!
!.not used..! !REVISION HISTORY:
!.not used..!  Initial code.
!.not used..!
!.not used..!EOP
!.not used..!-------------------------------------------------------------------------
!.not used..!BOC
!.not used..      if (pr_diag) Write (6,*) 'calcHorizontalMassFlux called by ', procID
!.not used..
!.not used..      if (first) then
!.not used..         first = .false.
!.not used..         call calcDelta_xyt (cosp, tdt, i2_gl, ju1_gl, j2_gl)
!.not used..      end if
!.not used..
!.not used..      xMassFlux(:,:,:) = 0.0d0
!.not used..      yMassFlux(:,:,:) = 0.0d0
!.not used..
!.not used..      if (INTERP_MASS_FLUX_COMPONENTS .or. (met_grid_type == 'C')) then
!.not used..
!.not used..         do ik = k1, k2
!.not used..            do ij = ju1, j2
!.not used..               do il = i1, i2
!.not used..                  if (met_grid_type == 'A') then  ! A grid
!.not used..                     tmp_crx = dtdx5(ij) * (uu(il,ij,ik) + uu(il-1,ij, ik))
!.not used..                     tmp_cry = dtdy5     * (vv(il,ij,ik) + vv(il,  ij-1,ik))
!.not used..                  else  ! C grid
!.not used..                     tmp_crx = dtdx(ij) * uu(il-1,ij,  ik)
!.not used..                     tmp_cry = dtdy     * vv(il,  ij-1,ik)
!.not used..                  endif
!.not used..
!.not used..                  ! Calculate E-W horizontal mass flux.
!.not used..
!.not used..                  tmp_pu1 = 0.5d0 * (delpm(il,ij,ik) + delpm(il-1,ij,ik))
!.not used..                  xMassFlux(il,ij,ik) = tmp_crx * tmp_pu1
!.not used..
!.not used..                  ! Calculate N-S horizontal mass flux.
!.not used..
!.not used..                  tmp_pu2 = 0.5d0 * (delpm(il,ij,ik) + delpm(il,ij-1,ik))
!.not used..                  yMassFlux(il,ij,ik) = tmp_cry * tmp_pu2 * cose(ij)
!.not used..               enddo
!.not used..            enddo
!.not used..         enddo
!.not used..
!.not used..      else
!.not used..
!.not used..         do ik = k1, k2
!.not used..            do ij = ju1, j2
!.not used..               do il = i1, i2
!.not used..                  xMassFlux(il, ij, ik) = dtdx5(ij) *  &
!.not used..     &                 ( uu(il-1,ij,ik) * delpm(il-1,ij,ik) +  &
!.not used..     &                   uu(il  ,ij,ik) * delpm(il  ,ij,ik)  )
!.not used..
!.not used..                  yMassFlux(il, ij, ik) = dtdy5 * cose(ij) *  &
!.not used..     &                 ( vv(il,ij-1,ik) * delpm(il,ij-1,ik) +  &
!.not used..     &                   vv(il,ij  ,ik) * delpm(il,ij  ,ik)  )
!.not used..
!.not used..               enddo
!.not used..            enddo
!.not used..         enddo
!.not used..
!.not used..      endif
!.not used..
!.not used..!     ------------------------------------------------
!.not used..!     Always done on Master, so only "wrap" necessary.
!.not used..!     ------------------------------------------------
!.not used..
!.not used..      call wrapMaster_3du (xMassFlux, i1, i2, ju1, j2, k1, k2,    &
!.not used..                           ilo, ihi, julo, jhi, gmi_nborder)
!.not used..      call wrapMaster_3du (yMassFlux, i1, i2, ju1, j2, k1, k2,    &
!.not used..                           ilo, ihi, julo, jhi, gmi_nborder)
!.not used..
!.not used..      return
!.not used..
!.not used..      end subroutine calcHorizontalMassFlux
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcVerticalMassFlux
!
! !INTERFACE:
!
      subroutine calcVerticalMassFlux &
     &  (dbk, dps_ctm, dpi, zMassFlux, &
     &   pr_diag, procID, i1, i2, ju1, j2, k1, k2)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: procID
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      ! difference in bi across layer - the dSigma term
      real*8 , intent(in) :: dbk(k1:k2)
      ! CTM surface pressure tendency; sum over vertical of dpi
      ! calculated from original mass fluxes (mb)
      ! dps_ctm(:,:) = Sum (dpi(:,:,:), dim=3)
      real*8 , intent(in) :: dps_ctm(i1:i2, ju1:j2)
      ! divergence at a grid point; used to calculate vertical motion (mb)
      real*8 , intent(in) :: dpi    (i1:i2, ju1:j2, k1:k2)
!
! !OUTPUT PARAMETERS:
      real*8, intent(out) :: zMassFlux(i1:i2, ju1:j2, k1:k2)
!
! !DESCRIPTION:
!  Calculates the vertical mass flux.
!
! !LOCAL VARIABLES:
      integer :: ik
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write (6,*) 'calcVerticalMassFlux called by ', procID
!
      ! Compute vertical mass flux from mass conservation.
!
      zMassFlux(:,:,k1) = dpi(:,:,k1) - (dbk(k1) * dps_ctm(i1:i2,ju1:j2))
!
      zMassFlux(:,:,k2) = 0.0d0
!
      do ik = k1 + 1, k2 - 1
         zMassFlux(:,:,ik) =  zMassFlux (:,:,ik-1) +  &
     &                dpi(:,:,ik) - (dbk(ik) * dps_ctm(i1:i2,ju1:j2))
!
      end do
!
      return
!
      end subroutine calcVerticalMassFlux
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcDelta_xyt
!
! !INTERFACE:
!
      subroutine calcDelta_xyt (cosp, tdt, i2_gl, ju1_gl, j2_gl)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer          , intent(in) :: i2_gl, ju1_gl, j2_gl
      ! model time step (s)
      real*8           , intent(in) :: tdt
      ! cosine of grid box centers
      real*8           , intent(in) :: cosp (ju1_gl:j2_gl)
!
! !DESCRIPTION:
!  Compute the variables $\Delta x$, $\Delta y$,
!  $\frac{\Delta t}{\Delta x}$, and $\frac{\Delta t}{\Delta y}$.
!
! !LOCAL VARIABLES:
      integer :: ij
      real*8  :: ri2, rj2m1
      real*8  :: dl            ! spacing in longitude (rad)
      real*8  :: dp            ! spacing in latitude  (rad)
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      Allocate (dtdx  (ju1_gl:j2_gl))
      Allocate (dtdx5 (ju1_gl:j2_gl))
      Allocate (deltax(ju1_gl:j2_gl))
      dtdx   = 0.0d0
      dtdx5  = 0.0d0
      deltax = 0.0d0
!
      ri2   = i2_gl
      rj2m1 = j2_gl - 1
!
      dl    = 2.0d0 * GMI_PI / ri2
      dp    = GMI_PI / rj2m1
!
      deltay = RADEAR * dp
      dtdy  = tdt / deltay
      dtdy5 = 0.5d0 * dtdy
!
!      dtdx  (ju1_gl) = 0.0d0
!      dtdx5 (ju1_gl) = 0.0d0
!      deltax(ju1_gl) = 0.0d0
!
!      do ij = ju1_gl + 1, j2_gl - 1
      do ij = ju1_gl, j2_gl
         deltax(ij) = dl * RADEAR * cosp(ij)
         dtdx (ij)  = tdt / deltax(ij)
         dtdx5(ij)  = 0.5d0 * dtdx(ij)
      end do
!
!      dtdx  (j2_gl) = 0.0d0
!      dtdx5 (j2_gl) = 0.0d0
!      deltax(j2_gl) = 0.0d0
!
      return
!
      end subroutine calcDelta_xyt
!EOC
!------------------------------------------------------------------------------
      end module GmiMassFluxes_mod
!
