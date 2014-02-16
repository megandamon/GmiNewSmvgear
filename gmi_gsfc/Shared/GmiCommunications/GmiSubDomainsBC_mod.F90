!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiSubDomainsBC_mod
!
! !INTERFACE:
!
      module GmiSubDomainsBC_mod
!
! !USES:
      use GmiGhostZones_mod, only : Msg3d_Ew      , Msg3d_Ns
      use GmiWrapMaster_mod, only : wrapMaster_2d, wrapMaster_3du
      use GmiWrapMaster_mod, only : wrapMaster_4d, wrapMaster_3dv
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: subDomainsBc
      public  :: Gmi_Bc2d_For_Sub
      public  :: Gmi_Bc3du_For_Sub
      public  :: Gmi_Bc3dv_For_Sub
      public  :: Gmi_Bc4d_For_Sub
!
! !DESCRIPTION:
! Worker processor communication routines for ghost zones.
!
! !AUTHOR:
!  John Tannahill, LLNL , jrt@llnl.gov
!  Jules Kouatchou, NASA GSFC, Jules.Kouatchou-1@nasa.gov
!
!EOP
!-----------------------------------------------------------------------------
      CONTAINS
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gmi_Bc2d_For_Sub
!
! !INTERFACE:
!
      subroutine Gmi_Bc2d_For_Sub  &
     &  (var_2d, msg_id, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, &
     &   ewflag, nspoleu, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain,  &
     &   msg_comm, numDomains)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: numDomains
      integer, intent(in) :: numLonDomains, numLatDomains, gmi_nborder
      integer, intent(in) :: i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi
      integer, intent(in) :: msg_comm
      integer, intent(in) :: nb, sb, eb, wb
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(julo:jhi)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
      integer, intent(in) :: msg_id
!
! !INPUT/OUTPUT PARAMETERS:
      real*8  :: var_2d(ilo:ihi, julo:jhi)
!
! !DESCRIPTION:
! This routine updates the ghost zones for the specified 2D array.
!
! !LOCAL VARIABLES:
      real*8  :: var_3d(ilo:ihi, julo:jhi, 1:1)
!
!EOP
!---------------------------------------------------------------------------
!BOC
      if (gmi_nborder /= 0) then
!
          var_3d(:,:,1) = var_2d(:,:)
!
          call subDomainsBc  &
     &      (var_3d, ewflag, nspoleu, 0,  &
     &       ilo, ihi, julo, jhi, 1, 1,  &
     &       numLonDomains, numLatDomains, msg_id, &
     &       gmi_nborder, i2_gl, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain, 0, 1, &
     &       msg_comm, numDomains)
!
          var_2d(:,:) = var_3d(:,:,1)
!
      end if
!
      return
!
      end subroutine Gmi_Bc2d_For_Sub
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gmi_Bc3du_For_Sub
!
! !INTERFACE:
!
      subroutine Gmi_Bc3du_For_Sub  &
     &  (var_3du, msg_id, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &   ewflag, nspoleu, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain,  &
     &   msg_comm, numDomains)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: numDomains
      integer, intent(in) :: numLonDomains, numLatDomains, gmi_nborder
      integer, intent(in) :: i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2
      integer, intent(in) :: msg_comm
      integer, intent(in) :: nb, sb, eb, wb
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(julo:jhi)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
      integer, intent(in) :: msg_id
!
! !INPUT/OUTPUT PARAMETERS:
      real*8  :: var_3du(ilo:ihi, julo:jhi, k1:k2) ! input 3D/julo array
!
! !DESCRIPTION:
!   This routine updates the ghost zones for the specified 3D/julo array.
!
!EOP
!---------------------------------------------------------------------------
!BOC
      if (gmi_nborder /= 0) then
!
          call subDomainsBc  &
     &      (var_3du, ewflag, nspoleu, 0,  &
     &       ilo, ihi, julo, jhi, k1, k2,  &
     &       numLonDomains, numLatDomains, msg_id, &
     &       gmi_nborder, i2_gl, nb, sb, eb, wb,  &
     &       north_domain, south_domain, east_domain, west_domain, 0, 1, &
     &       msg_comm, numDomains)
!
!
      end if
!
      return
!
      end subroutine Gmi_Bc3du_For_Sub
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gmi_Bc3dv_For_Sub
!
! !INTERFACE:
!
      subroutine Gmi_Bc3dv_For_Sub  &
     &  (var_3dv, msg_id, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, jv1, j2, i2_gl, ilo, ihi, jvlo, jhi, k1, k2, &
     &   ewflag, nspoleu, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain,  &
     &   msg_comm, numDomains)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: numLonDomains, numLatDomains, gmi_nborder
      integer, intent(in) :: numDomains
      integer, intent(in) :: i1, i2, jv1, j2, i2_gl, ilo, ihi, jvlo, jhi, k1, k2
      integer, intent(in) :: msg_comm
      integer, intent(in) :: nb, sb, eb, wb
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(jvlo:jhi)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
      integer, intent(in) :: msg_id
!
! !INPUT/OUTPUT PARAMETERS:
      real*8  :: var_3dv(ilo:ihi, jvlo:jhi, k1:k2)
!
! !DESCRIPTION:
!   This routine updates the ghost zones for the specified 3D/jvlo array.
!EOP
!---------------------------------------------------------------------------
!BOC
      if (gmi_nborder /= 0) then
!
          call subDomainsBc  &
     &      (var_3dv, ewflag, nspoleu, 0, ilo, ihi, jvlo, jhi, k1, k2,  &
     &       numLonDomains, numLatDomains, msg_id, gmi_nborder, i2_gl, &
     &       nb, sb, eb, wb, north_domain, south_domain, east_domain, &
     &       west_domain, 0, 1, msg_comm, numDomains)
!
      end if
!
      return
!
      end subroutine Gmi_Bc3dv_For_Sub
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gmi_Bc4d_For_Sub
!
! !INTERFACE:
!
      subroutine Gmi_Bc4d_For_Sub  &
     &  (var_4d, msg_ids, num_d4, &
     &   numLonDomains, numLatDomains, gmi_nborder, &
     &   i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2, &
     &   ewflag, nspoleu, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain,  &
     &   msg_comm, numDomains)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: numLonDomains, numLatDomains, gmi_nborder
      integer, intent(in) :: numDomains
      integer, intent(in) :: i1, i2, ju1, j2, i2_gl, ilo, ihi, julo, jhi, k1, k2
      integer, intent(in) :: msg_comm
      integer, intent(in) :: num_d4 ! size of 4th dimension (note always starts at 1)
      integer, intent(in) :: msg_ids(num_d4)
      integer, intent(in) :: nb, sb, eb, wb
      integer, intent(in) :: ewflag(ilo:ihi)
      integer, intent(in) :: nspoleu(julo:jhi)
      integer, intent(in) :: north_domain, south_domain, east_domain, west_domain
!
! !INPUT/OUTPUT PARAMETERS:
      real*8  :: var_4d (ilo:ihi, julo:jhi, k1:k2, num_d4)
!
! !DESCRIPTION:
!   This routine updates the ghost zones for the specified 4D array.
!
! !LOCAL VARIABLES:
      integer      :: ic
!EOP
!---------------------------------------------------------------------------
!BOC
      if (gmi_nborder /= 0) then
!
          do ic = 1, num_d4
!
            call subDomainsBc  &
     &        (var_4d(:,:,:,ic), ewflag, nspoleu, 0,  &
     &         ilo, ihi, julo, jhi, k1, k2,  &
     &         numLonDomains, numLatDomains, msg_ids(ic), &
     &         gmi_nborder, i2_gl, nb, sb, eb, wb,  &
     &         north_domain, south_domain, east_domain, west_domain, 0, 1, &
     &        msg_comm, numDomains)
!
          end do
      end if
!
      return
!
      end subroutine Gmi_Bc4d_For_Sub
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: subDomainsBc
!
! !INTERFACE:
!
      subroutine subDomainsBc  &
     &  (array, ewflag, nspole, polar_flag,  ilo, ihi, jlo, jhi, k1, k2, &
     &   numLonDomains, numLatDomains, msg_id, nborder, imax, nb, sb, eb, wb,  &
     &   north_domain, south_domain, east_domain, west_domain, nsperiod, &
     &   ewperiod, msg_comm, numDomains)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: numDomains
      integer, intent(in) :: ilo, ihi, jlo, jhi, k1, k2
               ! east-west (dateline) indicator array (precalculated and mesh-dependent)
      integer, intent(in) :: ewflag(ilo:ihi)
               ! north-south pole     indicator array (precalculated and mesh-dependent)
      integer, intent(in) :: nspole(jlo:jhi)
      integer, intent(in) :: polar_flag
               ! number of longitudinal subdomains
      integer, intent(in) :: numLonDomains
               ! number of latitudinal  subdomains
      integer, intent(in) :: numLatDomains
               ! a unique integer to indicate the message ID
      integer, intent(in) :: msg_id
               ! message passing communicator
      integer, intent(in) :: msg_comm
               ! number of ghost zones defining the border frame around the
               ! subdomain (e.g., 1 for second order, 2 for fourth order)
      integer, intent(in) :: nborder
      integer, intent(in) :: imax
               ! flag to indicate whether a subdomain is on north boundary
      integer, intent(in) :: nb
               ! flag to indicate whether a subdomain is on south boundary
      integer, intent(in) :: sb
               ! flag to indicate whether a subdomain is on east  boundary
      integer, intent(in) :: eb
               ! flag to indicate whether a subdomain is on west  boundary
      integer, intent(in) :: wb
               ! local subdomain to the north of this one
      integer, intent(in) :: north_domain
               ! local subdomain to the south of this one
      integer, intent(in) :: south_domain
               ! local subdomain to the east of this one
      integer, intent(in) :: east_domain
               ! local subdomain to the west of this one
      integer, intent(in) :: west_domain
               ! flag to indicate north-south periodicity
      integer, intent(in) :: nsperiod
               ! flag to indicate east-west   periodicity
      integer, intent(in) :: ewperiod
!
! !INPUT/OUTPUT PARAMETERS:
               ! array to apply boundary conditions to
      real*8 , intent(inOut) :: array (ilo:ihi, jlo:jhi, k1:k2)
!
! !DESCRIPTION:
!   This routine implements the horizontal boundary conditions on a real
!   array distributed across the subdomains.  These boundaries may be
!   interface boundary conditions for interior subdomains, may require a
!   polar boundary condition near the singular points, or may be a
!   reflection boundary condition in the case of subdomains with global
!   extent in one of the directions.
!
! !LOCAL VARIABLES:
      integer :: i, j
      integer :: lpe, lpw
      integer :: lpn, lps
      integer :: ebx, wbx
      integer :: nbx, sbx
      integer :: ewperiodx, weperiodx
      integer :: nsperiodx, snperiodx
      real*8  :: pass_from_east (jlo:jhi, k1:k2, nborder)
      real*8  :: pass_from_west (jlo:jhi, k1:k2, nborder)
      real*8  :: pass_from_north(ilo:ihi, k1:k2, nborder)
      real*8  :: pass_from_south(ilo:ihi, k1:k2, nborder)
      real*8  :: pass_to_east   (jlo:jhi, k1:k2, nborder)
      real*8  :: pass_to_west   (jlo:jhi, k1:k2, nborder)
      real*8  :: pass_to_north  (ilo:ihi, k1:k2, nborder)
      real*8  :: pass_to_south  (ilo:ihi, k1:k2, nborder)
!EOP
!---------------------------------------------------------------------------
!BOC
      pass_from_east (:,:,:) = 0.0d0; pass_from_west (:,:,:) = 0.0d0
      pass_from_north(:,:,:) = 0.0d0; pass_from_south(:,:,:) = 0.0d0
      pass_to_east   (:,:,:) = 0.0d0; pass_to_west   (:,:,:) = 0.0d0
      pass_to_north  (:,:,:) = 0.0d0; pass_to_south  (:,:,:) = 0.0d0
!
!     -------------------------------------------------------------
!     Zero out 2 and 3 dimensional arrays in the polar ghost zones.
!     -------------------------------------------------------------
!
      if (polar_flag == 0) then
!
        do j = jlo, jhi
!
          if (nspole(j) /= 0) then
!
            do i = ilo, ihi
              array(i,j,:) = 0.0d0
            end do
!
          end if
!
        end do
!
      else
!
        do j = jlo, jhi
!
          if (nspole(j) == -1) then
!
            do i = ilo, ihi
              array(i,j,:) = array(i,j+1,:)
            end do
!
          end if
!
          if (nspole(j) == 1) then
!
            do i = ilo, ihi
              array(i,j,:) = array(i,j-1,:)
            end do
!
          end if
!
        end do
!
      end if
!
!
!     ---------------------------------------------------------------------
!     Impose the periodic boundary condition via the ghost zones in the
!     longitudinal (meridional) direction on the 3 dimensional arrays for
!     domain decompositions involving a single longitudinal subdomain;
!     or
!     impose the interface boundary condition via a message passing call in
!     the longitudinal (meridional) direction on the 3 dimensional arrays
!     for domain decompositions involving multiple longitudinal subdomains.
!     ---------------------------------------------------------------------
!
      if ((numLonDomains == 1) .and. (ewperiod == 1)) then
!
        do i = ilo+1, ihi-1
!
          if (ewflag(i-1) < 0) then
!
            do j = jlo, jhi
              array(i-1,j,:) = array(imax-1+i,j,:)
            end do
!
          end if
!
          if (ewflag(i+1) > 0) then
!
            do j = jlo, jhi
              array(i+1,j,:) = array(i+1-imax,j,:)
            end do
!
          end if
!
        end do
!
      end if
!
!
      if (numLonDomains > 1) then
!
        lpe = east_domain
        lpw = west_domain
!
        if (lpe == -9999) then
          ebx       = 1
          ewperiodx = 0
        else
          ebx       = eb
          ewperiodx = ewperiod
        end if
!
        if (lpw == -9999) then
          wbx       = 1
          weperiodx = 0
        else
          wbx       = wb
          weperiodx = ewperiod
        end if
!
!       =============
        call Msg3d_Ew  &
!       =============
     &    (array, ilo, ihi, jlo, jhi, k1, k2,  &
     &     nborder, ebx, wbx, lpe, lpw, ewperiodx, weperiodx, msg_id,  &
     &     pass_to_east, pass_to_west, pass_from_east, pass_from_west, &
     &     msg_comm)
!
      end if
!
!
!     --------------------------------------------------------------------
!     Impose the inteface boundary condition via a message passing call in
!     the latitudinal direction on the 3 dimensional arrays for domain
!     decompositions involving multiple latitudinal subdomains.
!     --------------------------------------------------------------------
!
      if (numLatDomains > 1) then
!
        nbx = nb
        sbx = sb
!
        nsperiodx = nsperiod
        snperiodx = nsperiod
!
        lpn = north_domain
        lps = south_domain
!
        if (lpn == -9999) nbx = 1
        if (lps == -9999) sbx = 1
!
        if (lpn == -9999) nsperiodx = 0
        if (lps == -9999) snperiodx = 0
!
!       =============
        call Msg3d_Ns  &
!       =============
     &    (array, ilo, ihi, jlo, jhi, k1, k2,  &
     &     nborder, nbx, sbx, lpn, lps, nsperiodx, snperiodx, msg_id,  &
     &     pass_to_north, pass_to_south,  &
     &     pass_from_north, pass_from_south, &
     &     msg_comm)
!
      end if
!
      return
!
      end subroutine subDomainsBc
!EOC
!------------------------------------------------------------------------
!
      end module GmiSubDomainsBC_mod
!
