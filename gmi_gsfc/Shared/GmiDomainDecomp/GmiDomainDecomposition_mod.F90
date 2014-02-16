!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  GmiDomainDecomposition_mod
!
! !INTERFACE:
!
#include "MAPL_ErrLog.h"
!
      module GmiDomainDecomposition_mod
!
! !USES:
      use ESMF_Mod
      use GmiPrintError_mod, only : GmiPrintError
      use GmiGrid_mod, only : t_gmiGrid, Set_i1, Set_i2, Set_ju1, Set_jv1,     &
     &       Set_j2, Set_k1, Set_k2, Set_ilong, Set_ilat, Set_ivert,           &
     &       Set_itloop, Set_ilo, Set_ihi, Set_julo, Set_jvlo, Set_jhi
      use GmiGrid_mod, only : Get_i1, Get_i2, Get_ju1, Get_j2, Get_k1, Get_k2, &
     &       Set_i1_gl, Set_i2_gl, Set_ju1_gl, Set_j2_gl, Set_j1p, Get_j1p, &
     &       Set_jv1_gl, Set_k1_gl, Set_k2_gl, Set_ilo_gl, Set_j2p, &
     &       Set_ihi_gl, Set_julo_gl, Set_jhi_gl, Set_gmi_nborder, Set_jvlo_gl,&
     &       Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_j2_gl,          &
     &       Get_jv1_gl, Get_jv1, Get_ivert, Get_ilong, Get_ilat,   &
     &       Get_itloop, Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jvlo_gl,     &
     &       Get_jhi_gl, Get_numSpecies, Get_gmi_nborder, Set_numSpecies
      use GmiMessagePassing_mod, only : stopCode, writeMpiError
      use GmiMessagePassing_mod, only : sendInteger, receiveInteger,           &
     &       broadcastInteger
!
      implicit none
!
#     include "GmiParameters.h"
#     include "gmi_phys_constants.h"
#     include "mpif.h"
!
! !PUBLIC MEMBER FUNCTIONS:
      private
      public  :: domainDecomposition,       setSubDomainData
      public  :: Get_coscen, Get_cosedg, Get_cose, Get_cosp, Get_dlatr
      public  :: Get_geofac, Get_latdeg, Get_londeg, Get_mcor
      public  :: Get_rel_area, Get_rel_areaGlob, Get_mcorGlob, Get_geofac_pc
      public  :: Set_procID               , Get_procID
      public  :: Set_rootProc             , Get_rootProc
      public  :: Set_eastWestFlag         , Get_eastWestFlag
      public  :: Set_northSouthPole       , Get_northSouthPole
      public  :: Get_map1_u               , Get_map2_u
      public  :: Get_map1_v               , Get_map2_v
      public  :: Get_mapi_all
      public  :: Get_globEastDomain       , Get_globWestDomain
      public  :: Get_globNorthDomain      , Get_globSouthDomain
      public  :: Get_globEastBorder       , Get_globWestBorder
      public  :: Get_globNorthBorder      , Get_globSouthBorder
      public  :: Get_globMostEasternDomain, Get_globMostWesternDomain
      public  :: Get_communicatorWorld
      public  :: Get_communicatorSouthPole, Get_communicatorNorthPole
      public  :: Set_communicatorWorld
      public  :: Get_southBorder          , Get_westBorder
      public  :: Get_northBorder          , Get_eastBorder
      public  :: Get_southNeighbor        , Get_westNeighbor
      public  :: Get_northNeighbor        , Get_eastNeighbor
      public  :: Get_iAmRootProc          , Set_iAmRootProc
      public  :: Get_numDomains           , Set_numDomains
      public  :: Get_numLatDomains        , Set_numLatDomains
      public  :: Get_numLonDomains        , Set_numLonDomains
!
! !DESCRIPTION:
! This module contains routines to carry out the domain decomposition.
! It is assumed that all the processors read the namelist file to obtain
! the global domain dimensions. The root processor then performs the
! domain decomposition and send subdomain information to worker processors.
! This module populates the derived types gmiDomain (declared here) and
! gmiGrid (defined outside).
!
! It important to note that the root processor is not involved in the
! model calculations.
!
! !PUBLIC DATA MEMBERS:
!
      public  ::   t_gmiDomain
!
      type  t_gmiDomain
        private
        logical  :: iAmRootProc           ! Am I the root processor?
        integer  :: procID                ! processor idenfier
        integer  :: rootProc              ! root processor identifier
        integer  :: numDomains            ! total number of subdomains
                                          ! same as the total number of processors
                                          ! (equal to numLonDomains * numLatDomains)
        integer  :: numLonDomains         ! number of longitudinal subdomains
        integer  :: numLatDomains         ! number of latitudinal  subdomains
        integer  :: eastNeighbor          ! east  neighbor subdomain
        integer  :: westNeighbor          ! west  neighbor subdomain
        integer  :: northNeighbor         ! north neighbor subdomain
        integer  :: southNeighbor         ! south neighbor subdomain
        integer  :: eastBorder            ! east  border subdomain
        integer  :: westBorder            ! west  border subdomain
        integer  :: northBorder           ! north border subdomain
        integer  :: southBorder           ! south border subdomain
        integer  :: mostEasternDomain     ! easternmost domain
        integer  :: mostWesternDomain     ! westernmost domain
        integer  :: communicatorWorld     ! world communicator
        integer  :: communicatorSouthPole ! communicator for procs covering south pole
        integer  :: communicatorNorthPole ! communicator for procs covering north pole
        !
        !
                 ! the following two arrays define any global ghost zones,
                 ! if 0       => not a ghost zone;
                 ! if -1 or 1 => lower or upper phone zone, respectively.
                 ! Used for communications.
        integer, pointer :: eastWestFlag         (:) => null()
        integer, pointer :: northSouthPole       (:) => null()
                 ! The following five arrays define the ranges of various 1D/2D
                 ! data mappings to each processor. map1 includes border zones
                 ! while map2 does not. The interpretation of these map arrays:
                 !     (1,1,:) = southwest corner
                 !     (2,1,:) = southesat corner
                 !     (1,2,:) = northwest corner
                 !     (2,2,:) = northeast corner
        integer, pointer :: map1_u           (:,:,:) => null()
        integer, pointer :: map2_u           (:,:,:) => null()
        integer, pointer :: map1_v           (:,:,:) => null()
        integer, pointer :: map2_v           (:,:,:) => null()
        integer, pointer :: mapi_all           (:,:) => null()
                 ! global array for east  neighbor subdomains
        integer, pointer :: globEastDomain       (:) => null()
                 ! global array for west  neighbor subdomains
        integer, pointer :: globWestDomain       (:) => null()
                 ! global array for north neighbor subdomains
        integer, pointer :: globNorthDomain      (:) => null()
                 ! global array for south neighbor subdomains
        integer, pointer :: globSouthDomain      (:) => null()
                 ! global array indicating that a subdomain is at east  border
        integer, pointer :: globEastBorder       (:) => null()
                 ! global array indicating that a subdomain is at west  border
        integer, pointer :: globWestBorder       (:) => null()
                 ! global array indicating that a subdomain is at north border
        integer, pointer :: globNorthBorder      (:) => null()
                 ! global array indicating that a subdomain is at south border
        integer, pointer :: globSouthBorder      (:) => null()
                 ! easternmost domain for a given latitude
        integer, pointer :: globMostEasternDomain(:) => null()
                 ! westernmost domain for a given latitude
        integer, pointer :: globMostWesternDomain(:) => null()
                 ! processor number as a function of subdomain number
!
        real*8 , pointer :: latdeg(:) => null() ! latitude (deg)
        real*8 , pointer :: londeg(:) => null() ! longitude (deg)
        real*8 , pointer :: dlatr (:) => null() ! latitude of zone center in
                                                ! latitude direction (rad)
        real*8 , pointer :: coscen(:) => null() ! cosine of latitude of zone
                                                ! centers = cos(dlatr)
        real*8 , pointer :: cosedg(:) => null() ! cosine of latitude of zone
                                                ! edges   = cos(dlatr2)
        real*8 , pointer :: cose  (:) => null() ! cosine of grid box edges
        real*8 , pointer :: cosp  (:) => null() ! cosine of grid box centers
        real*8 , pointer :: mcor(:,:) => null() ! grid box area (m^2)
        real*8 , pointer :: mcorGlob(:,:) => null() ! grid box area (m^2)
        real*8 , pointer :: rel_area(:,:) => null() ! relative grid box area
        real*8 , pointer :: rel_areaGlob(:,:) => null() ! relative grid box area
        real*8           :: geofac_pc
        real*8 , pointer :: geofac  (:) => null() ! geometrical factor for
                                                  ! meridional advection; it uses
                                                  ! correct spherical geometry,
                                                  ! and replaces acosp as the
                                                  ! meridional geometrical factor
                                                  ! in tpcore special geometrical
                                                  ! factor (geofac) for Polar cap
      end type  t_gmiDomain
!
! !AUTHOR:
! Jules Kouatchou, NASA/GSFC, Jules.Kouatchou-1@nasa.gov
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domainDecomposition
!
! !INTERFACE:
!
      subroutine domainDecomposition(self, gmiGrid, numProcessors, config)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: numProcessors
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_gmiGrid  ), intent(inOut) :: gmiGrid
      type (Esmf_Config), intent(inOut) :: config
      type (t_gmiDomain), intent(inOut) :: self
!
! !DESCRIPTION:
! This routine performs the domain decomposition. It is assumed that all the
! processors read the namelist file and have the global domain data (contained
! in the derived type gmiGrid). Here, the root processor uses global domain
! dimensions to assign a subdomain to each processor.
!
! !LOCAL VARIABLES:
      integer :: gmi_nborder
      integer :: i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, k1_gl, k2_gl
      integer :: i1, i2, ju1, jv1, j2, k1, k2, ilong, ilat, ivert, itloop
      integer :: ilo, ihi, julo, jvlo, jhi, j1p, j2p
      integer :: ilo_gl, ihi_gl, julo_gl, jvlo_gl, jhi_gl
!
      integer, parameter :: max_tag_msg = 30
      integer            :: tag(max_tag_msg)
      integer            :: msgSize, i, domain, lproc, ierr
      integer            :: numDomains, numLonDomains, numLatDomains, numSpecies
      integer            :: eastNeighbor, westNeighbor, southNeighbor, northNeighbor
      integer            :: eastBorder, westBorder, southBorder, northBorder
      integer            :: mostEasternDomain, mostWesternDomain
      integer, allocatable :: globEastDomain       (:)
      integer, allocatable :: globWestDomain       (:)
      integer, allocatable :: globSouthDomain      (:)
      integer, allocatable :: globNorthDomain      (:)
      integer, allocatable :: globEastBorder       (:)
      integer, allocatable :: globWestBorder       (:)
      integer, allocatable :: globSouthBorder      (:)
      integer, allocatable :: globNorthBorder      (:)
      integer, allocatable :: globMostEasternDomain(:)
      integer, allocatable :: globMostWesternDomain(:)
      integer, allocatable :: eastWestFlag         (:)
      integer, allocatable :: northSouthPole       (:)
      integer :: tArray(1)
      integer :: procID, rootProc, commuWorld
      logical :: iAmRootProc
      integer :: RC, STATUS
      character(len=ESMF_MAXSTR) :: IAm , err_msg
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      IAm = "domainDecomposition"
!
      ! Get data previously set during the MPI initialization
!
      call Get_procID           (self, procID)
      call Get_iAmRootProc      (self, iAmRootProc)
      call Get_rootProc         (self, rootProc   )
      call Get_communicatorWorld(self, commuWorld)
!
      !################################
      ! Get data from the resource file
      !################################
!
      call ESMF_ConfigGetAttribute(config, self%numLonDomains, &
     &                label   = "NX:",&
     &                default = 1, rc=STATUS )
      !VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, self%numLatDomains, &
     &                label   = "NY:",&
     &                default = 1, rc=STATUS )
      !VERIFY_(STATUS)
!
      self%numDomains = self%numLonDomains * self%numLatDomains
!
!
      ! Set global domain information
!
      call ESMF_ConfigGetAttribute(config, i2_gl, &
     &                label   = "IM:",&
     &                default = 72, rc=STATUS )
      !VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, j2_gl, &
     &                label   = "JM:",&
     &                default = 26, rc=STATUS )
      !VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, k2_gl, &
     &                label   = "LM:",&
     &                default = 29, rc=STATUS )
      !VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, gmi_nborder, &
     &                label   = "gmi_nborder:",&
     &                default = 4, rc=STATUS )
      !VERIFY_(STATUS)
!
      call ESMF_ConfigGetAttribute(config, numSpecies, &
     &                label   = "numSpecies:",&
     &                default = 1, rc=STATUS )
      !VERIFY_(STATUS)
!
      i1_gl   = 1
      ju1_gl  = 1
      jv1_gl  = ju1_gl - 1
      k1_gl   = 1
      ilo_gl  = i1_gl  - gmi_nborder
      ihi_gl  = i2_gl  + gmi_nborder
      julo_gl = ju1_gl - gmi_nborder
      jvlo_gl = jv1_gl - gmi_nborder
      jhi_gl  = j2_gl  + gmi_nborder
!
      call Set_k1         (gmiGrid, k1_gl  )
      call Set_k2         (gmiGrid, k2_gl  )
      call Set_k1_gl      (gmiGrid, k1_gl  )
      call Set_k2_gl      (gmiGrid, k2_gl  )
      call Set_i1_gl      (gmiGrid, i1_gl  )
      call Set_i2_gl      (gmiGrid, i2_gl  )
      call Set_j2_gl      (gmiGrid, j2_gl  )
      call Set_ju1_gl     (gmiGrid, ju1_gl)
      call Set_jv1_gl     (gmiGrid, jv1_gl)
      call Set_ilo_gl     (gmiGrid, ilo_gl  )
      call Set_ihi_gl     (gmiGrid, ihi_gl  )
      call Set_julo_gl    (gmiGrid, julo_gl )
      call Set_jvlo_gl    (gmiGrid, jvlo_gl )
      call Set_jhi_gl     (gmiGrid, jhi_gl  )
      call Set_gmi_nborder(gmiGrid, gmi_nborder)
      call Set_numSpecies (gmiGrid, numSpecies)
!
!
!... read in polar cap size from resource file - default=3
!      j1p     = 3
      call ESMF_ConfigGetAttribute(config, j1p, &
     &                label   = "j1p:",&
     &                default = 3, rc=STATUS )
!
      j2p     = j2_gl - j1p + 1
      call Set_j1p        (gmiGrid, j1p)
      call Set_j2p        (gmiGrid, j2p)
!
      !===========================================================
      ! Does error checking to determine proper initial setting.
      !===========================================================
!
      if (self%numDomains /= numProcessors) then
         call stopCode(commuWorld, "Check the number of processors")
      end if
!
      if ((i2_gl / self%numLonDomains) < gmi_nborder) then
         err_msg = 'IM/NX problem in Check_Nlvalue.'
         call GmiPrintError(err_msg, .true., 2, i2_gl, self%numLonDomains, &
     &           0, 0.0d0, 0.0d0)
      end if
!
      if ((j2_gl / self%numLatDomains) < gmi_nborder) then
         err_msg = 'JM/NY problem in Check_Nlvalue.'
         call GmiPrintError(err_msg, .true., 2, j2_gl, self%numLatDomains, &
     &           0, 0.0d0, 0.0d0)
      end if
!
      if (gmi_nborder < 3) then
         err_msg = 'gmi_nborder problem in Check_Nlvalue.'
         call GmiPrintError (err_msg, .true., 1, j1p, 0, 0, 0.0d0, 0.0d0)
      end if
!
      if (k2_gl < 6) then
         err_msg = 'LM must be >= 6 in Check_Nlvalue.'
         call GmiPrintError (err_msg, .true., 1, k2_gl, 0, 0, 0.0d0, 0.0d0)
      end if
!
      if (numSpecies == 0) then
         err_msg = 'numSpecies problem in the resource file.'
         call GmiPrintError  &
     &       (err_msg, .true., 1, numSpecies, 0, 0, 0.0d0, 0.0d0)
      end if
!
      if (MAX_NUM_CONST_GIO /= MAX_NUM_CONST) then
         err_msg = 'MAX_NUM_CONST_GIO/MAX_NUM_CONST problem in the rc file.'
         call GmiPrintError  &
     &    (err_msg, .true., 2, MAX_NUM_CONST_GIO, MAX_NUM_CONST,  &
     &     0, 0.0d0, 0.0d0)
      end if
!
      !=======================================
      ! Allocate variables in the derived type
      !=======================================
!
      call allocateSubDomain(self)
!
      call Allocate_eastWestFlag  (self, ilo_gl, ihi_gl)
      call Allocate_northSouthPole(self, julo_gl, jhi_gl)
      allocate(eastWestFlag  (ilo_gl:ihi_gl))
      allocate(northSouthPole(julo_gl:jhi_gl))
!
       call setGhostZones (self, i1_gl, i2_gl, ju1_gl, j2_gl, gmi_nborder)
!
      !=========================================
      ! Determines subdomains by all processors
      !=========================================
!
      call setSubdomains (self, gmi_nborder, &
     &           i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl)
!
      call setSubDomainDimensions(self, gmiGrid)
!
      !==================================
      ! Initialize subgroup communicators
      !==================================
!
      call Get_ju1    (gmiGrid, ju1    )
      call Get_j2     (gmiGrid, j2     )
      call initPoleCommunicators      (self, ju1, ju1_gl, j2, j2_gl )
!
      !=====================================================================
      ! The root processor sends additional informationi to other processors
      !=====================================================================
!
      self%eastNeighbor      = self%globEastDomain(procID+1)
      self%westNeighbor      = self%globWestDomain(procID+1)
      self%northNeighbor     = self%globNorthDomain(procID+1)
      self%southNeighbor     = self%globSouthDomain(procID+1)
      self%eastBorder        = self%globEastBorder(procID+1)
      self%westBorder        = self%globWestBorder(procID+1)
      self%northBorder       = self%globNorthBorder(procID+1)
      self%southBorder       = self%globSouthBorder(procID+1)
      self%mostEasternDomain = self%globMostEasternDomain(procID+1)
      self%mostWesternDomain = self%globMostWesternDomain(procID+1)
!
      return
!
      end subroutine domainDecomposition
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setSubDomainDimensions
!
! !INTERFACE:
!
      subroutine setSubDomainDimensions(self, gmiGrid)
!
      implicit none
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_gmiGrid  ), intent(inOut) :: gmiGrid
      type (t_gmiDomain), intent(inOut) :: self
!
! !DESCRIPTION:
! This routine allows the other processors to get the domain decomposition
! data from the root processor.
!
! !LOCAL VARIABLES:
      integer :: i1, i2, ju1, jv1, j2, k1, k2, ilong, ilat, ivert, itloop
      integer :: ilo, ihi, julo, jvlo, jhi, i1_gl, i2_gl, jv1_gl, ju1_gl, j2_gl, j1p
      integer :: numDomains
      integer, allocatable :: mapi_all(:, :)
      integer, allocatable :: map1_u(:, :, :)
      integer, allocatable :: map1_v(:, :, :)
      integer, allocatable :: map2_u(:, :, :)
      integer, allocatable :: map2_v(:, :, :)
      integer :: lpx, domain, msgSize, numSpecies, gmi_nborder
      integer :: bArray(1)
      integer :: procID, rootProc, commuWorld
!EOP
!------------------------------------------------------------------------------
!BOC
      call Get_i1_gl (gmiGrid, i1_gl)
      call Get_i2_gl (gmiGrid, i2_gl)
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_jv1_gl(gmiGrid, jv1_gl)
      call Get_j2_gl (gmiGrid, j2_gl)
      call Get_k1    (gmiGrid, k1)
      call Get_k2    (gmiGrid, k2)
      call Get_j1p   (gmiGrid, j1p)
!
      call Get_procID           (self, procID    )
      call Get_rootProc         (self, rootProc  )
      call Get_numDomains       (self, numDomains)
      call Get_communicatorWorld(self, commuWorld)
!
      allocate(mapi_all(2,    numDomains))
      allocate(map1_u  (2, 2, numDomains))
      allocate(map1_v  (2, 2, numDomains))
      allocate(map2_u  (2, 2, numDomains))
      allocate(map2_v  (2, 2, numDomains))
!
      call Get_map1_u  (self, map1_u)
      call Get_map1_v  (self, map1_v)
      call Get_map2_u  (self, map2_u)
      call Get_map2_v  (self, map2_v)
      call Get_mapi_all(self, mapi_all)
!
      ilo  = map2_u(1,1,procID+1)
      ihi  = map2_u(2,1,procID+1)
      julo = map2_u(1,2,procID+1)
      jvlo = map2_v(1,2,procID+1)
      jhi  = map2_u(2,2,procID+1)
!
      i1   = map1_u(1,1,procID+1)
      i2   = map1_u(2,1,procID+1)
      ju1  = map1_u(1,2,procID+1)
      jv1  = map1_v(1,2,procID+1)
      j2   = map1_u(2,2,procID+1)
!
      ilat   = j2 - ju1 + 1
      ilong  = i2 -  i1 + 1
      ivert  = k2 -  k1 + 1
!
      itloop = ilat * ilong * ivert
!
      call Set_ilat  (gmiGrid, ilat)
      call Set_ilong (gmiGrid, ilong)
      call Set_ivert (gmiGrid, ivert)
      call Set_itloop(gmiGrid, itloop)
!
      call Set_ilo  (gmiGrid, ilo)
      call Set_ihi  (gmiGrid, ihi)
      call Set_julo (gmiGrid, julo)
      call Set_jvlo (gmiGrid, jvlo)
      call Set_jhi  (gmiGrid, jhi)
!
      call Set_i1   (gmiGrid, i1)
      call Set_i2   (gmiGrid, i2)
      call Set_ju1  (gmiGrid, ju1)
      call Set_jv1  (gmiGrid, jv1)
      call Set_j2   (gmiGrid, j2 )
!
      call Get_numSpecies (gmiGrid, numSpecies)
      call Get_gmi_nborder(gmiGrid, gmi_nborder)
!
      if (procID == rootProc) then
!
         Write (6,*) '------------------------------------'
         Write (6,*) '------Global Domain Dimensions------'
         Write (6,*) '------------------------------------'
         Write (6,*) 'i1    = ', i1_gl,  '  i2  = ', i2_gl
         Write (6,*) 'ju1   = ', ju1,    '  j2  = ', j2_gl
         Write (6,*) 'jv1   = ', jv1_gl, '  j1p = ', j1p
         Write (6,*) 'k1    = ', k1,     '  k2  = ', k2
         Write (6,*) 'numSpecies  = ', numSpecies
         Write (6,*) 'gmi_nborder = ', gmi_nborder
!
         Write (6,*) '------------------------------------'
         Write (6,*) '-------Sub Domain Dimensions--------'
         Write (6,*) '------------------------------------'
         do domain = 1, numDomains
            lpx = domain - 1
!
            i1   = map1_u(1,1,lpx+1)
            i2   = map1_u(2,1,lpx+1)
            ju1  = map1_u(1,2,lpx+1)
            jv1  = map1_v(1,2,lpx+1)
            j2   = map1_u(2,2,lpx+1)
!
            ilat   = j2 - ju1 + 1
            ilong  = i2 -  i1 + 1
            ivert  = k2 -  k1 + 1
!
            itloop = ilat * ilong * ivert
!
            Write (6,*) '------Dimensions for Processor ', lpx
            Write (6,*) 'i1    = ', i1,    '  i2     = ', i2
            Write (6,*) 'ju1   = ', ju1,   '  j2     = ', j2
            Write (6,*) 'jv1   = ', jv1,   '  j1p    = ', j1p
            Write (6,*) 'k1    = ', k1,    '  k2     = ', k2
            Write (6,*) 'ilat  = ', ilat,  '  ilong  = ', ilong
            Write (6,*) 'ivert = ', ivert, '  itloop = ', itloop
         end do
      end if
!
      return
!
      end subroutine setSubDomainDimensions
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocateSubDomain
!
! !INTERFACE:
!
    subroutine allocateSubDomain(self)
!
       implicit none
!
! !INPUT/OUTPUT PARAMETERS:
       type(t_gmiDomain), intent(inOut) :: self
!
! !DESCRIPTION:
!  Allocates subdomain related arrays.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      call Allocate_map1_u  (self)
      call Allocate_map2_u  (self)
      call Allocate_map1_v  (self)
      call Allocate_map2_v  (self)
      call Allocate_mapi_all(self)
!
      call Allocate_globEastDomain(self)
      call Allocate_globWestDomain(self)
      call Allocate_globNorthDomain(self)
      call Allocate_globSouthDomain(self)
      call Allocate_globEastBorder(self)
      call Allocate_globWestBorder(self)
      call Allocate_globNorthBorder(self)
      call Allocate_globSouthBorder(self)
      call Allocate_globMostEasternDomain(self)
      call Allocate_globMostWesternDomain(self)
!
      return
!
      end subroutine allocateSubDomain
!EOC
!-----------------------------------------------------------------------------
!BOP
! !IROUTINE: initPoleCommunicators
!
! !INTERFACE:
!
      subroutine initPoleCommunicators (self, ju1, ju1_gl, j2, j2_gl )
!
      implicit none
!
! !INPUT PARAMETERS:
      integer , intent(in   ) :: ju1, ju1_gl, j2, j2_gl
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_gmiDomain), intent(inOut) :: self
!
! !DESCRIPTION:
!   Constructs groups of processor numbers and communicators
!   used for message passing.
!
! !LOCAL VARIABLES:
      integer :: color
      integer :: ierr
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      if (ju1 == ju1_gl) then  ! South Pole.
        color = 1
      else
        color = 0
      end if
!
!     ===================
      call MPI_Comm_Split  &
!     ===================
     &  (self%communicatorWorld, color, 0, self%communicatorSouthPole, ierr)
!
      if (ierr /= MPI_SUCCESS) then
        call writeMpiError (self%communicatorWorld, .true., ierr)
      end if
!
      if (j2 == j2_gl) then  ! North Pole.
        color = 1
      else
        color = 0
      end if
!
!     ===================
      call MPI_Comm_Split  &
!     ===================
     &  (self%communicatorWorld, color, 0, self%communicatorNorthPole, ierr)
!
      if (ierr /= MPI_SUCCESS) then
        call writeMpiError (self%communicatorWorld, .true., ierr)
      end if
!
      return
!
      end subroutine initPoleCommunicators
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setSubdomains
!
! !INTERFACE:
!
      subroutine setSubdomains &
     &  (self, gmi_nborder, i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl)
!
      implicit none
!
! !INPUT PARAMETERS:
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl
      integer, intent(in) :: gmi_nborder
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_gmiDomain), intent(inOut) :: self
!
! !DESCRIPTION:
!   This routine sets up the GMI subdomains and maps.
!
! !LOCAL VARIABLES:
      integer :: il, ij, ik, lproc
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
!     ===================
      call getSubdomains  &
!     ===================
     &  (self%globWestDomain, self%globEastDomain, self%globSouthDomain,       &
     &   self%globNorthDomain, self%globWestBorder, self%globEastBorder,       &
     &   self%globSouthBorder, self%globNorthBorder, self%map1_u, self%map2_u, &
     &   gmi_nborder, self%globMostWesternDomain, self%globMostEasternDomain,  &
     &   self%numDomains, self%numLonDomains, self%numLatDomains, i1_gl,       &
     &   i2_gl, ju1_gl, j2_gl, 1)
!
!
      do ik = 1, self%numDomains
        do ij = 1, 2
          do il = 1, 2
!
            if ((il == 1) .and.  &
     &          (ij == 2) .and.  &
     &          (self%map1_u(1,2,ik) == ju1_gl)) then
!
              self%map1_v(il,ij,ik) =  &
     &          self%map1_u(il,ij,ik) + (jv1_gl - ju1_gl)
              self%map2_v(il,ij,ik) =  &
     &          self%map2_u(il,ij,ik) + (jv1_gl - ju1_gl)
!
            else
!
              self%map1_v(il,ij,ik) = self%map1_u(il,ij,ik)
              self%map2_v(il,ij,ik) = self%map2_u(il,ij,ik)
!
            end if
!
          end do
        end do
      end do
!
      do ik = 1,self%numDomains
         lproc =  ik
!
         self%mapi_all(1,lproc) = self%map1_u(1,1,ik)
         self%mapi_all(2,lproc) = self%map1_u(2,1,ik)
      end do
!
      return
!
      end subroutine setSubdomains
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setGhostZones
!
! !INTERFACE:
!
      subroutine setGhostZones (self, i1_gl, i2_gl, ju1_gl, j2_gl, gmi_nborder)
!
      implicit none
!
! !INPUT PARAMETERS:
!
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, gmi_nborder
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_gmiDomain), intent(inOut) :: self
!
! !DESCRIPTION:
!   This routine sets flags to indicate i and j ghost zones at end of arrays.
!
! !LOCAL VARIABLES:
      integer :: il, ij
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      self%eastWestFlag (:) = 0
      self%northSouthPole(:) = 0
!
      if (gmi_nborder /= 0) then
!
        do il = 1, gmi_nborder
          self%eastWestFlag (i1_gl  - il) = -1
          self%eastWestFlag (i2_gl  + il) =  1
        end do
!
        do ij = 1, gmi_nborder
          self%northSouthPole(ju1_gl - ij) = -1
          self%northSouthPole(j2_gl  + ij) =  1
        end do
!
      end if
!
      return
!
      end subroutine setGhostZones
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getSubdomains
!
! !INTERFACE:
!
      subroutine getSubdomains  &
     &  (west_domain, east_domain, south_domain, north_domain,  &
     &   wb, eb, sb, nb, map1, map2, nborder, most_west, most_east,  &
     &   numDomains, numLonDomains, numLatDomains, ilo, ihi, jlo, jhi,  &
     &   option)
!
      implicit none
!
! !INPUT PARAMETERS:
               ! total number of subdomains (numLonDomains * numLatDomains)
      integer, intent(in ) :: numDomains
               ! number of longitudinal subdomains
      integer, intent(in ) :: numLonDomains
               ! number of latitudinal  subdomains
      integer, intent(in ) :: numLatDomains
               ! the number of border zones; in finite diff. algorithms,
               ! this is determined by the order of the calculation; for
               ! second order use 1, fourth order use 2, etc.
      integer, intent(in ) :: nborder
               ! the global dimensions (without border zones)
      integer, intent(in ) :: ilo, ihi, jlo, jhi
               ! a flag to determine how to decompose the global domain
      integer, intent(in ) :: option
!
! !OUTPUT PARAMETERS:
               ! The following four variables determine each subdomain's neighboring
               ! subdomains and are used for interprocessor communication =>
               !
               ! west      neighbor subdomain
      integer, intent(out) :: west_domain (numDomains)
               ! east      neighbor subdomain
      integer, intent(out) :: east_domain (numDomains)
               ! south     neighbor subdomain
      integer, intent(out) :: south_domain(numDomains)
               ! north     neighbor subdomain
      integer, intent(out) :: north_domain(numDomains)
               ! The following four variables determine the location
               ! of a subdomain within the context of the global domain
               ! and are used for the imposition of horizontal boundary conditions:
               !
               ! flag to indicate that a subdomain is at west  border
      integer, intent(out) :: wb(numDomains)
               ! flag to indicate that a subdomain is at east  border
      integer, intent(out) :: eb(numDomains)
               ! flag to indicate that a subdomain is at south border
      integer, intent(out) :: sb(numDomains)
               ! flag to indicate that a subdomain is at north border
      integer, intent(out) :: nb(numDomains)
               ! a map defining the subdomains; map1 does not include
               ! border zones whereas map2 does.
               ! The interpretation of these map arrays is =>
               !         (1,1,numDomains) = southwest corner
               !         (2,1,numDomains) = southeast corner
               !         (1,2,numDomains) = northwest corner
               !         (2,2,numDomains) = northeast corner
      integer, intent(out) :: map1(2,2,numDomains)
      integer, intent(out) :: map2(2,2,numDomains)
               ! westernmost domain for a given latitude
      integer, intent(out) :: most_west(numDomains)
               ! easternmost domain for a given latitude
      integer, intent(out) :: most_east(numDomains)
!
! !DESCRIPTION:
!   This routine creates a subdomain decomposition determined by the user
!   input.  The default decomposition is to attempt to make the subdomains
!   as uniform as possible in each horizontal direction.  The user has the
!   option of defining the decomposition in either of these directions as
!   namelist input.
!
! !LOCAL VARIABLES:
      integer  delta_i, delta_j
      integer  extra_i, extra_j
      integer  nzone_i, nzone_j
      integer  domain                    ! local subdomain number
      integer  domain_i, domain_j
      integer  east, west, north, south
      integer  ispan, jspan
      integer  westmost
      integer :: map(2, 2, numDomains)
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
      map(:,:,:) = 0
!
!
!     ====================
      if (numDomains == 1) then
!     ====================
!
!       -----------------------------------------------------
!       Create a trivial map if only one subdomain is chosen.
!       -----------------------------------------------------
!
        wb(1) = 1
        eb(1) = 1
        nb(1) = 1
        sb(1) = 1
!
        west_domain(1)  = 1
        east_domain(1)  = 1
        north_domain(1) = 1
        south_domain(1) = 1
!
        map1(1,1,1) = ilo
        map1(2,1,1) = ihi
        map1(1,2,1) = jlo
        map1(2,2,1) = jhi
!
        map2(1,1,1) = ilo - nborder
        map2(2,1,1) = ihi + nborder
        map2(1,2,1) = jlo - nborder
        map2(2,2,1) = jhi + nborder
!
        most_west(1) = 1
        most_east(1) = 1
!
!     ====
      else
!     ====
!
!       -----------------------------------------------------------------
!       Create a map of the neighboring subdomains.  The numbering of the
!       subdomains is assumed to be monotonically increasing from west to
!       east and from south to north.
!       -----------------------------------------------------------------
!
        do domain = 1, numDomains
          wb(domain) = 0
          eb(domain) = 0
          nb(domain) = 0
          sb(domain) = 0
        end do
!
!
        domain = 0
!
        do domain_j = 1, numLatDomains
!
          south = domain_j - 2
          north = domain_j
!
          if (domain_j == 1)             south = numLatDomains - 1
          if (domain_j == numLatDomains) north = 0
!
          do domain_i = 1, numLonDomains
!
            east = domain_i + 1
            west = domain_i - 1
!
            if (domain_i == 1)             west = numLonDomains
            if (domain_i == numLonDomains) east = 1
!
            domain = domain + 1
!
            west_domain(domain) = (domain_j - 1) * numLonDomains + west
            east_domain(domain) = (domain_j - 1) * numLonDomains + east
!
            south_domain(domain) = (south * numLonDomains) + domain_i
            north_domain(domain) = (north * numLonDomains) + domain_i
!
            if (domain_i == 1)             wb(domain) = 1
            if (domain_i == numLonDomains) eb(domain) = 1
            if (domain_j == 1)             sb(domain) = 1
            if (domain_j == numLatDomains) nb(domain) = 1
!
            if (domain_i == 1)             westmost = domain
!
            most_west(domain) = westmost
            most_east(domain) = westmost + numLonDomains - 1
!
          end do
!
        end do
!
!!!!!! New additions for the n-processor decomposition
!
        west_domain (:) = west_domain (:) - 1
        east_domain (:) = east_domain (:) - 1
        south_domain(:) = south_domain(:) - 1
        north_domain(:) = north_domain(:) - 1
!
        most_west   (:) = most_west   (:) - 1
        most_east   (:) = most_east   (:) - 1
!
!!!!!! End of section
!
!       -----------------------------------------------
!       Create a map using two dimensional coordinates.
!       -----------------------------------------------
!
        nzone_i = ihi - ilo + 1
        nzone_j = jhi - jlo + 1
!
!       ------------------------------------------------------------------
!       Find the location of the of the southwest corner of the subdomains
!       for the "nearly uniform" option.
!       ------------------------------------------------------------------
!
        if (option == 1) then
!
          map(1,2,1) = jlo
          delta_j    = nzone_j / numLatDomains
          extra_j    = nzone_j - delta_j * numLatDomains + 1
          jspan      = 0
!
          do domain_j = 2, numLatDomains
!
            if (domain_j <= extra_j) jspan = delta_j + 1
            if (domain_j >  extra_j) jspan = delta_j
!
            map(1,2,domain_j)   = map(1,2,domain_j-1) + jspan
            map(2,2,domain_j-1) = map(1,2,domain_j)   - 1
!
          end do
!
          domain_j = numLatDomains + 1
!
          if (domain_j <= extra_j) jspan = delta_j + 1
          if (domain_j >  extra_j) jspan = delta_j
!
          if (numLatDomains /= 1) then
            map(2,2,numLatDomains) = map(2,2,numLatDomains-1) + jspan
          else
            map(2,2,numLatDomains) = jhi
          end if
!
          map(1,1,1) = ilo
          delta_i    = nzone_i / numLonDomains
          extra_i    = nzone_i - delta_i * numLonDomains + 1
          ispan      = 0
!
          do domain_i = 2, numLonDomains
!
            if (domain_i <= extra_i) ispan = delta_i + 1
            if (domain_i >  extra_i) ispan = delta_i
!
            map(1,1,domain_i)   = map(1,1,domain_i-1) + ispan
            map(2,1,domain_i-1) = map(1,1,domain_i)   - 1
!
          end do
!
          domain_i = numLonDomains + 1
!
          if (domain_i <= extra_i) ispan = delta_i + 1
          if (domain_i >  extra_i) ispan = delta_i
!
          if (numLonDomains /= 1) then
            map(2,1,numLonDomains) = map(2,1,numLonDomains-1) + ispan
          else
            map(2,1,numLonDomains) = ihi
          end if
!
        end if
!
!       -----------------------------------------------------------
!       Convert the two dimensional map into a one dimensional map.
!       -----------------------------------------------------------
!
        domain = 1
!
        do domain_j = 1,numLatDomains
          do domain_i = 1,numLonDomains
!
            map1(1,1,domain) = map(1,1,domain_i)
            map1(2,1,domain) = map(2,1,domain_i)
            map1(1,2,domain) = map(1,2,domain_j)
            map1(2,2,domain) = map(2,2,domain_j)
!
            map2(1,1,domain) = map1(1,1,domain) - nborder
            map2(2,1,domain) = map1(2,1,domain) + nborder
            map2(1,2,domain) = map1(1,2,domain) - nborder
            map2(2,2,domain) = map1(2,2,domain) + nborder
!
            domain = domain + 1
!
          end do
        end do
!
!     ======
      end if
!     ======
!
      return
!
      end subroutine getSubdomains
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setSubDomainData
!
! !INTERFACE:
!
      subroutine setSubDomainData(self, gmiGrid, pr_diag, metdata_name_org, &
     &               met_grid_type)
!
      implicit none
!
! !INPUT PARAMETER:
      character (len=*) , intent(in) :: metdata_name_org ! first  part of
                                                         ! metdata_name, e.g., "NCAR"
      character (len=*) , intent(in) :: met_grid_type    ! met grid type, 'A' or 'C'
      logical           , intent(in) :: pr_diag
      type (t_gmiGrid)  , intent(in) :: gmiGrid
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_gmiDomain), intent(inOut) :: self
!
! !DESCRIPTION:
!  Determines the following domain data:
!  \begin{itemize}
!  \item latitude (deg)
!  \item longitude (deg)
!  \item grid box area (m^2)
!  \item relative grid box area
!  \item
!  \end{itemize}
!  The call to this routine is made by all the processors. There is no need
!  for the root processor to send the computed data to the other ones.
!
! !LOCAL VARIABLES:
      integer :: i1, i2, ju1, j2, i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl, j1p
      real*8 , allocatable :: dlatr2(:)
!EOP
!-------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'setSubDomainData called by ', self%procID
!
      call Get_i1    (gmiGrid, i1   )
      call Get_i2    (gmiGrid, i2   )
      call Get_ju1   (gmiGrid, ju1  )
      call Get_j2    (gmiGrid, j2   )
      call Get_i1_gl (gmiGrid, i1_gl)
      call Get_i2_gl (gmiGrid, i2_gl)
      call Get_ju1_gl(gmiGrid, ju1_gl)
      call Get_jv1_gl(gmiGrid, jv1_gl)
      call Get_j2_gl (gmiGrid, j2_gl)
      call Get_j1p   (gmiGrid, j1p)
!
      allocate(self%latdeg(ju1_gl:j2_gl))
      self%latdeg = 0.0d0
!
      allocate(self%londeg( i1_gl:i2_gl))
      self%londeg = 0.0d0
!
      allocate(self%dlatr (ju1_gl:j2_gl))
      self%dlatr = 0.0d0
!
      allocate(self%coscen(ju1_gl:j2_gl))
      self%coscen = 0.0d0
!
      allocate(self%cosedg(ju1_gl:j2_gl+1))
      self%cosedg = 0.0d0
!
      allocate(self%cose  (ju1_gl:j2_gl+1))
      self%cose = 0.0d0
!
      allocate(self%cosp  (ju1_gl:j2_gl))
      self%cosp = 0.0d0
!
      allocate(self%mcor    (i1:i2,ju1:j2))
      self%mcor = 0.0d0
!
      allocate(self%mcorGlob(i1_gl:i2_gl,ju1_gl:j2_gl))
      self%mcorGlob = 0.0d0
!
      allocate(self%rel_area(i1:i2,ju1:j2))
      self%rel_area = 0.0d0
!
      allocate(self%rel_areaGlob(i1_gl:i2_gl,ju1_gl:j2_gl))
      self%rel_areaGlob = 0.0d0
!
      allocate(self%geofac  (ju1_gl:j2_gl))
      self%geofac = 0.0d0
!
      allocate (dlatr2(ju1_gl:j2_gl+1))
!
      call Set_Grid_Imp (self, pr_diag, dlatr2, metdata_name_org, met_grid_type,  &
     &       i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl)
!
      call calcGridBoxSurfArea (self, pr_diag, dlatr2, i1, i2, ju1, j2, j1p, i1_gl, &
     &               i2_gl, ju1_gl, jv1_gl, j2_gl)
!
      deallocate(dlatr2)
!
      return
!
      end subroutine setSubDomainData
!EOC
!-------------------------------------------------------------------------
!BOP
! !IROUTINE: Set_Grid_Imp
!
! !INTERFACE
!
      subroutine Set_Grid_Imp  (self, pr_diag, dlatr2, metdata_name_org,       &
     &       met_grid_type, i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl)
!
!
      implicit none
!
!#     include "gmi_grid_constants.h"
!
! !INPUT PARAMETER:
      logical          , intent(in) :: pr_diag
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl
      ! first  part of metdata_name, e.g., "NCAR"
      character (len=*) :: metdata_name_org
      ! met grid type, 'A' or 'C'
      character (len=*)  :: met_grid_type
!
! !OUTPUT PARAMETERS:
      real*8, intent(out) :: dlatr2(ju1_gl:j2_gl+1)
!
! !INPUT/OUTPUT PARAMETERS:
      type(t_gmiDomain), intent(inOut) :: self
!
! !DESCRIPTION:
!   This routine defines the grid.
!   Note that for DAO, 90.0 deg is a special fix.
!
! !LOCAL VARIABLES:
!      real*8, parameter :: TINY_DAO   = 0.01d0
      real*8 :: polecap_res
      real*8, parameter :: TINY_DLATR = 1.0d-06
      character (len=75) :: err_msg
      integer :: il, ij, ilong
      real*8  :: fac
      real*8  :: rijx
      real*8  :: ustart_lat, uend_lat,  vstart_lat, latdeg_spac
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'Set_Grid_Imp called by ', self%procID
!
      ilong = i2_gl - i1_gl + 1
      fac = 360.0d0 / ilong
      self%londeg(i1_gl) = 0.0d0
!
      do il=i1_gl+1,i2_gl
        self%londeg(il) = (il - 1) * fac
      end do
!
      latdeg_spac = 180.d0 / ((j2_gl - ju1_gl + 1.0d0) - 1.0d0)
      ustart_lat  = -90.0d0
      uend_lat    =  90.0d0
      vstart_lat  = ustart_lat + (latdeg_spac/2.0d0)
      polecap_res = (vstart_lat - ustart_lat)
!
      dlatr2(ju1_gl) = ustart_lat
      dlatr2(ju1_gl+1) = vstart_lat
      do ij=ju1_gl+2,j2_gl
        dlatr2(ij) = dlatr2(ij-1) + LATDEG_SPAC
      enddo
      dlatr2(j2_gl+1) = dlatr2(j2_gl) + polecap_res
!
      do ij=ju1_gl,j2_gl
        self%latdeg(ij) = (dlatr2(ij)+dlatr2(ij+1))/2
      enddo
!
      dlatr2(:) = dlatr2(:) * RADPDEG
!
      do ij=ju1_gl,j2_gl+1
        self%cosedg(ij) = Cos (dlatr2(ij))
      enddo
!
      do ij=ju1_gl,j2_gl
        self%dlatr(ij) = (dlatr2(ij)+dlatr2(ij+1))/2
        self%coscen(ij) = Cos(self%dlatr(ij))
      enddo
!
!      if (j2_gl == 46) then
!
!        latdeg_spac = LATDEG_SPAC_46
!        ustart_lat  = USTART_LAT_46
!        uend_lat    = UEND_LAT_46
!        vstart_lat  = VSTART_LAT_46
!
!      else if (j2_gl == 64) then
!
!        latdeg_spac = LATDEG_SPAC_64
!        ustart_lat  = USTART_LAT_64
!        uend_lat    = UEND_LAT_64
!        vstart_lat  = VSTART_LAT_64
!
!      else if (j2_gl == 91) then
!
!        latdeg_spac = LATDEG_SPAC_91
!        ustart_lat  = USTART_LAT_91
!        uend_lat    = UEND_LAT_91
!        vstart_lat  = VSTART_LAT_91
!
!      else if (j2_gl == 181) then
!
!        latdeg_spac = LATDEG_SPAC_181
!        ustart_lat  = USTART_LAT_181
!        uend_lat    = UEND_LAT_181
!        vstart_lat  = VSTART_LAT_181
!
!      else
!
!        err_msg = 'Problem with j2_gl in Set_Grid_Imp.'
!        call GmiPrintError (err_msg, .true., 1, j2_gl, 0, 0, 0.0d0, 0.0d0)
!
!      end if
!
!!     -------------------------
!!     Determine dlatr & coscen.
!!     -------------------------
!      self%dlatr (ju1_gl) = (ustart_lat + polecap_res) * RADPDEG
!!      self%dlatr (ju1_gl) = (ustart_lat + TINY_DAO) * RADPDEG
!
!      self%coscen(:)   = 0.0d0
!      self%coscen(ju1_gl) = Cos (self%dlatr(ju1_gl))
!
!      do ij = ju1_gl + 1, j2_gl - 1
!
!        rijx = ij - 1
!
!        self%dlatr(ij)  = (ustart_lat + (rijx * latdeg_spac)) * RADPDEG
!
!        if (Abs (self%dlatr(ij)) < TINY_DLATR) self%dlatr(ij) = 0.0d0
!
!        self%coscen(ij) = Cos (self%dlatr(ij))
!
!      end do
!
!      self%dlatr (j2_gl) = (uend_lat - polecap_res) * RADPDEG
!!      self%dlatr (j2_gl) = (uend_lat - TINY_DAO) * RADPDEG
!      self%coscen(j2_gl) = Cos (self%dlatr(j2_gl))
!
!!     -----------------
!!     Determine latdeg.
!!     -----------------
!
!      self%latdeg(:) = self%dlatr(:) * DEGPRAD
!
!!     --------------------------
!!     Determine dlatr2 & cosedg.
!!     --------------------------
!
!      dlatr2(:) = 0.0d0
!      self%cosedg(:) = 0.0d0
!
!      do ij = jv1_gl, j2_gl - 1
!
!        rijx = ij - 1
!
!        dlatr2(ij) = (vstart_lat + (rijx * latdeg_spac)) * RADPDEG
!
!        self%cosedg(ij) = Cos (dlatr2(ij))
!
!      end do
!
!!     ==================
!      call Set_Cose_Cosp  &
!!     ==================
!     &  (metdata_name_org, met_grid_type, self%coscen, self%cosedg, self%cose, &
!     &   self%cosp, ju1_gl, jv1_gl, j2_gl)
!
      self%cose = self%cosedg
      self%cosp = self%coscen
!
      return
!
      end subroutine Set_Grid_Imp
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcGridBoxSurfArea
!
! !INTERFACE:
!
      subroutine calcGridBoxSurfArea (self, pr_diag, dlatr2, i1, i2, ju1, j2, j1p, i1_gl, &
     &               i2_gl, ju1_gl, jv1_gl, j2_gl)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical          , intent(in) :: pr_diag
      integer, intent(in) :: i1, i2, ju1, j2, j1p
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, jv1_gl, j2_gl
      real*8 , intent(in) :: dlatr2(ju1_gl:j2_gl+1)
!
! !OUTPUT PARAMETERS:
      type (t_gmiDomain), intent(inOut) :: self
!
! !DESCRIPTION:
! This routine calculates mcor, the surface area of a grid volume;
! a mixing ratio to mass conversion term (related to horiz. area of
! grid zones -- that is, zone length in longitude * zone length in
! latitude; again, mcor(ju1) = mcor(j2) = poles and is a special
! fix = area of circular pole regions).
!
! !LOCAL VARIABLES:
!      real*8, parameter :: HALFPI = 0.5d0 * GMI_PI
      real*8, parameter :: TERM1  = 2.0d0 * GMI_PI * RADEAR * RADEAR
      integer :: ij
      real*8  :: ri2_gl, tot_area
!
!EOP
!------------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'calcGridBoxSurfArea called by ', self%procID
!
      ri2_gl = i2_gl
!
!     --------------------------
!     Determine mcorGlob(i1_gl:i2_gl,ju1_gl).
!     --------------------------
!
!      self%mcorGlob(:,ju1_gl) = TERM1 * (Cos (dlatr2(jv1_gl) + HALFPI) - 1.0d0) / ri2_gl
!
!     ---------------------------------
!     Determine globMCOR(i1_gl:i2_gl,ju1_gl+1:j2_gl-1).
!     ---------------------------------
!
!      do ij = ju1_gl + 1, j2_gl - 1
!
!        self%mcorGlob(:,ij) = TERM1 *  &
!     &    (Cos (dlatr2(ij) + HALFPI) - Cos (dlatr2(ij-1) + HALFPI)) / ri2_gl
!
!      end do
!
!     -------------------------
!     Determine mcorGlob(i1_gl:i2_gl,j2_gl).
!     -------------------------
!
!      self%mcorGlob(:,j2_gl) = TERM1 * (-1.0d0 - Cos (dlatr2(j2_gl-1) + HALFPI)) /  ri2_gl
!
      do ij=ju1_gl, j2_gl
        self%mcorGlob(:,ij) = TERM1 * (sin(dlatr2(ij+1)) - sin(dlatr2(ij))) / ri2_gl
      enddo
!
      tot_area = sum(self%mcorGlob(:,:))
      self%rel_areaGlob(:,:) = self%mcorGlob(:,:) / tot_area
!
      self%mcor    (i1:i2,ju1:j2) = self%mcorGlob(i1:i2,ju1:j2)
      self%rel_area(i1:i2,ju1:j2) = self%rel_areaGlob(i1:i2,ju1:j2)
!
      call calcAdvecFactors(self, pr_diag, i1_gl, i2_gl, ju1_gl, j2_gl, j1p)
!
      return
!
      end subroutine calcGridBoxSurfArea
!EOC
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcAdvecFactors
!
! !INTERFACE:
!
      subroutine calcAdvecFactors (self, pr_diag, i1_gl, i2_gl, ju1_gl, j2_gl, j1p)
!
      implicit none
!
! !INPUT PARAMETERS:
      logical, intent(in) :: pr_diag
      integer, intent(in) :: i1_gl, i2_gl, ju1_gl, j2_gl, j1p
      ! relative surface area of grid box (fraction)
!
! !INPUT/OUTPUT PARAMETERS:
      type (t_gmiDomain), intent(inOut) :: self
!
! !DESCRIPTION:
! Calculates the relative area of each grid box, and the
! geometrical factors used by the modified version of tpcore.  These
! geomoetrical factors DO assume that the space is regularly gridded,
! but do not assume any link between the surface area and the linear
! dimensions.
!
! !LOCAL VARIABLES:
      integer :: ij
      real*8  :: dp           ! spacing in latitude (rad)
      real*8  :: dp_pole      ! spacing in latitude at pole is 1/2 grid box
      real*8  :: ri2_gl
      real*8  :: rj2m1
      real*8  :: total_area
!EOP
!-----------------------------------------------------------------------------
!BOC
      if (pr_diag) Write(6,*) 'calcAdvecFactors called by ', self%procID
!
      ri2_gl = i2_gl
!
!     ---------------------------------------------------------
!     Calculate geometrical factor for meridional advection.
!     Note that it is assumed that all grid boxes in a latitude
!     band are the same.
!     ---------------------------------------------------------
!
      rj2m1 = j2_gl - 1
      dp    = GMI_PI / rj2m1
!
      do ij = ju1_gl, j2_gl
        self%geofac(ij) = dp / (2.0d0 * self%rel_areaGlob(1,ij) * ri2_gl)
      end do
!
!      dp_pole = dp*0.5d0
!      do ij = ju1_gl+1, j2_gl-1
!        self%geofac(ij) = dp / (2.0d0 * self%rel_areaGlob(1,ij) * ri2_gl)
!      end do
!      self%geofac(ju1_gl) = dp_pole / (2.0d0 * self%rel_areaGlob(1,ju1_gl) * ri2_gl)
!      self%geofac(j2_gl) = dp_pole / (2.0d0 * self%rel_areaGlob(1,j2_gl) * ri2_gl)
!
!
!.test that pole grid point is 1/2 size of rest of grid. 
!      self%geofac(ju1_gl) = self%geofac(ju1_gl)*0.5d0
!      self%geofac(j2_gl) = self%geofac(j2_gl)*0.5d0
!
!      if(j1p.eq.ju1_gl+1) then
!        self%geofac_pc = self%geofac(ju1_gl)
!      else
        self%geofac_pc = dp / (2.0d0*Sum(self%rel_areaGlob(1,ju1_gl:ju1_gl+1))*ri2_gl)
!      endif
!
!      print *,'geofac_pc: ', self%procID, j1p, self%geofac_pc
!      print *,'geofac: ', self%procID, j1p, self%geofac
!
      return
!
      end subroutine calcAdvecFactors
!EOC
!!-----------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Set_Cose_Cosp
!!
!! !INTERFACE:
!!
!      subroutine Set_Cose_Cosp  &
!     &  (metdata_name_org, met_grid_type, coscen, cosedg,  &
!     &   cose, cosp, ju1_gl, jv1_gl, j2_gl)
!!
!      implicit none
!!
!! !INPUT PARAMETERS:
!      ! first  part of metdata_name, e.g., "NCAR"
!      character (len=MAX_LENGTH_MET_NAME) :: metdata_name_org
!      ! met grid type, 'A' or 'C'
!      character (len=1)  :: met_grid_type
!      integer, intent(in) :: ju1_gl, jv1_gl, j2_gl
!      ! cosine of latitude of zone centers = cos(dlatr)
!      real*8 , intent(in) :: coscen(ju1_gl:j2_gl)
!      ! cosine of latitude of zone edges   = cos(dlatr2)
!      real*8 , intent(in) :: cosedg(ju1_gl:j2_gl+1)
!!
!! !OUTPUT PARAMETERS:
!      ! cosine of grid box edges
!      real*8 , intent(out) :: cose  (ju1_gl:j2_gl+1)
!      ! cosine of grid box centers
!      real*8 , intent(out) :: cosp  (ju1_gl:j2_gl)
!!
!! !DESCRIPTION:
!! This routine sets cose and cosp.
!!
!! !LOCAL VARIABLES:
!      real*8  :: dp     ! spacing in latitude (rad)
!      real*8  :: rj2m1
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!      rj2m1 = j2_gl - 1
!      dp    = GMI_PI / rj2m1
!
!      if (met_grid_type == 'A') then  ! A grid
!
!!       =========
!        call Cosa (dp, cose, cosp, ju1_gl, j2_gl)
!!       =========
!
!      else  ! C grid
!
!        if (metdata_name_org(1:3) /= 'DAO') then
!
!!         =========
!          call Cosg (cose, cosp, coscen, cosedg, ju1_gl, jv1_gl, j2_gl)
!!         =========
!
!        else
!
!!         ------------------------------------------------------------
!!         This is a special C grid case for consistency with GEOS-GCM;
!!         cos(grid edges) are average of cos(grid_centers).
!!         ------------------------------------------------------------
!
!!         =========
!          call Cosc (dp, cose, cosp, ju1_gl, j2_gl)
!!         =========
!
!        end if
!
!      end if
!
!      return
!
!      end subroutine Set_Cose_Cosp
!!EOC
!!-----------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE:
!!
!! !INTERFACE:
!!
!      subroutine Cosa  (dp, cose, cosp, ju1_gl, j2_gl)
!
!      implicit none
!
!! !INPUT PARAMETERS:
!      integer, intent(in) :: ju1_gl, j2_gl
!      ! spacing in latitude (rad)
!      real*8 , intent(in) :: dp
!!
!! !OUTPUT PARAMETERS:
!      ! cosine of grid box edges
!      real*8 , intent(out) :: cose(ju1_gl:j2_gl+1)
!      ! cosine of grid box centers
!      real*8 , intent(out) :: cosp(ju1_gl:j2_gl)
!!
!! !DESCRIPTION:
!! This routine computes the analytic cosine at cell edges (for A grid).
!!
!! !LOCAL VARIABLES:
!      integer :: ij, j2eq, j2m1
!      real*8  :: ph5, rijm1
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!      do ij = ju1_gl + 1, j2_gl
!        rijm1    = ij - 1
!        ph5      = (-0.5d0 * GMI_PI) + ((rijm1 - 0.5d0) * dp)
!        cose(ij) = Cos (ph5)
!      end do
!
!      j2m1 = j2_gl - 1
!      j2eq = (j2_gl + 1) / 2
!
!      if (j2m1 == (2 * (j2m1 / 2))) then
!
!        do ij = j2_gl, j2eq + 1, -1
!          cose(ij) = cose(j2_gl+2-ij)
!        end do
!
!      else
!
!!       ---------------------
!!       Cell edge at equator.
!!       ---------------------
!
!        cose(j2eq+1) = 1.0d0
!
!        do ij = j2_gl, j2eq + 2, -1
!          cose(ij) = cose(j2_gl+2-ij)
!        end do
!
!      end if
!
!      cosp(ju1_gl) = 0.0d0
!
!      do ij = ju1_gl + 1, j2_gl - 1
!        cosp(ij) = 0.5d0 * (cose(ij) + cose(ij+1))
!      end do
!
!      cosp(j2_gl) = 0.0d0
!
!      return
!
!      end subroutine Cosa
!!EOC
!!-----------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: Cosc
!!
!! !INTERFACE:
!!
!      subroutine Cosc (dp, cose, cosp, ju1_gl, j2_gl)
!!
!      implicit none
!!
!! !INPUT PARAMETERS:
!      integer, intent(in) :: ju1_gl, j2_gl
!      ! spacing in latitude (rad)
!      real*8 , intent(in) :: dp
!!
!! !OUTPUT PARAMETERS:
!      ! cosine of grid box edges
!      real*8 , intent(out) :: cose(ju1_gl:j2_gl+1)
!      ! cosine of grid box centers
!      real*8 , intent(out) :: cosp(ju1_gl:j2_gl)
!!
!! !DESCRIPTION:
!!   This routine defines cosines consistent with GEOS-GCM
!!   (using dycore2.0 or later; for C grid).
!!
!! !LOCAL VARIABLES:
!      integer :: ij
!      real*8  :: phi
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!      phi = -0.5d0 * GMI_PI
!
!      do ij = ju1_gl + 1, j2_gl - 1
!        phi      = phi + dp
!        cosp(ij) = Cos (phi)
!      end do
!
!      cosp(ju1_gl) = 0.0d0
!
!      do ij = ju1_gl + 1, j2_gl
!        cose(ij) = 0.5d0 * (cosp(ij) + cosp(ij-1))
!      end do
!
!      cosp(j2_gl) = 0.0d0
!
!      do ij = ju1_gl + 1, j2_gl - 1
!        cosp(ij) = 0.5d0 * (cose(ij) + cose(ij+1))
!      end do
!
!      return
!
!      end subroutine Cosc
!EOC
!------------------------------------------------------------------------------
    subroutine Set_iAmRootProc(self, iAmRootProc)
       implicit none
       logical          , intent(in   ) :: iAmRootProc
       type(t_gmiDomain), intent(inOut) :: self
       self%iAmRootProc = iAmRootProc
       return
    end subroutine Set_iAmRootProc
!-----------------------------------------------------------------------------
    subroutine Get_iAmRootProc(self, iAmRootProc)
       implicit none
       logical, intent(out) :: iAmRootProc
       type(t_gmiDomain), intent(in) :: self
       iAmRootProc = self%iAmRootProc
       return
    end subroutine Get_iAmRootProc
!-----------------------------------------------------------------------------
    subroutine Get_westBorder(self, westBorder)
       implicit none
       integer, intent(out) :: westBorder
       type(t_gmiDomain), intent(in) :: self
       westBorder = self%westBorder
       return
    end subroutine Get_westBorder
!-----------------------------------------------------------------------------
    subroutine Get_eastBorder(self, eastBorder)
       implicit none
       integer, intent(out) :: eastBorder
       type(t_gmiDomain), intent(in) :: self
       eastBorder = self%eastBorder
       return
    end subroutine Get_eastBorder
!-----------------------------------------------------------------------------
    subroutine Get_northBorder(self, northBorder)
       implicit none
       integer, intent(out) :: northBorder
       type(t_gmiDomain), intent(in) :: self
       northBorder = self%northBorder
       return
    end subroutine Get_northBorder
!-----------------------------------------------------------------------------
    subroutine Get_southBorder(self, southBorder)
       implicit none
       integer, intent(out) :: southBorder
       type(t_gmiDomain), intent(in) :: self
       southBorder = self%southBorder
       return
    end subroutine Get_southBorder
!-----------------------------------------------------------------------------
    subroutine Get_westNeighbor(self, westNeighbor)
       implicit none
       integer, intent(out) :: westNeighbor
       type(t_gmiDomain), intent(in) :: self
       westNeighbor = self%westNeighbor
       return
    end subroutine Get_westNeighbor
!-----------------------------------------------------------------------------
    subroutine Get_eastNeighbor(self, eastNeighbor)
       implicit none
       integer, intent(out) :: eastNeighbor
       type(t_gmiDomain), intent(in) :: self
       eastNeighbor = self%eastNeighbor
       return
    end subroutine Get_eastNeighbor
!-----------------------------------------------------------------------------
    subroutine Get_northNeighbor(self, northNeighbor)
       implicit none
       integer, intent(out) :: northNeighbor
       type(t_gmiDomain), intent(in) :: self
       northNeighbor = self%northNeighbor
       return
    end subroutine Get_northNeighbor
!-----------------------------------------------------------------------------
    subroutine Get_southNeighbor(self, southNeighbor)
       implicit none
       integer, intent(out) :: southNeighbor
       type(t_gmiDomain), intent(in) :: self
       southNeighbor = self%southNeighbor
       return
    end subroutine Get_southNeighbor
!-----------------------------------------------------------------------------
    subroutine Get_communicatorSouthPole(self, communicatorSouthPole)
       implicit none
       integer, intent(out) :: communicatorSouthPole
       type(t_gmiDomain), intent(in) :: self
       communicatorSouthPole = self%communicatorSouthPole
       return
    end subroutine Get_communicatorSouthPole
!-----------------------------------------------------------------------------
    subroutine Get_communicatorNorthPole(self, communicatorNorthPole)
       implicit none
       integer, intent(out) :: communicatorNorthPole
       type(t_gmiDomain), intent(in) :: self
       communicatorNorthPole = self%communicatorNorthPole
       return
    end subroutine Get_communicatorNorthPole
!-----------------------------------------------------------------------------
    subroutine Set_communicatorWorld(self, communicatorWorld)
       implicit none
       integer, intent(in) :: communicatorWorld
       type(t_gmiDomain), intent(inOut) :: self
       self%communicatorWorld = communicatorWorld
       return
    end subroutine Set_communicatorWorld
!-----------------------------------------------------------------------------
    subroutine Get_communicatorWorld(self, communicatorWorld)
       implicit none
       integer, intent(out) :: communicatorWorld
       type(t_gmiDomain), intent(in) :: self
       communicatorWorld = self%communicatorWorld
       return
    end subroutine Get_communicatorWorld
!-----------------------------------------------------------------------------
    subroutine Get_procID(self, procID)
       implicit none
       integer, intent(out) :: procID
       type(t_gmiDomain), intent(in) :: self
       procID = self%procID
       return
    end subroutine Get_procID
!-----------------------------------------------------------------------------
    subroutine Set_procID(self, procID)
       implicit none
       integer, intent(in) :: procID
       type(t_gmiDomain), intent(inOut) :: self
       self%procID = procID
       return
    end subroutine Set_procID
!-----------------------------------------------------------------------------
    subroutine Get_numLonDomains(self, numLonDomains)
       implicit none
       integer, intent(out) :: numLonDomains
       type(t_gmiDomain), intent(in) :: self
       numLonDomains = self%numLonDomains
       return
    end subroutine Get_numLonDomains
!-----------------------------------------------------------------------------
    subroutine Set_numLonDomains(self, numLonDomains)
       implicit none
       integer, intent(in) :: numLonDomains
       type(t_gmiDomain), intent(inOut) :: self
       self%numLonDomains = numLonDomains
       return
    end subroutine Set_numLonDomains
!-----------------------------------------------------------------------------
    subroutine Get_numLatDomains(self, numLatDomains)
       implicit none
       integer, intent(out) :: numLatDomains
       type(t_gmiDomain), intent(in) :: self
       numLatDomains = self%numLatDomains
       return
    end subroutine Get_numLatDomains
!-----------------------------------------------------------------------------
    subroutine Set_numLatDomains(self, numLatDomains)
       implicit none
       integer, intent(in) :: numLatDomains
       type(t_gmiDomain), intent(inOut) :: self
       self%numLatDomains = numLatDomains
       return
    end subroutine Set_numLatDomains
!-----------------------------------------------------------------------------
    subroutine Get_numDomains(self, numDomains)
       implicit none
       integer, intent(out) :: numDomains
       type(t_gmiDomain), intent(in) :: self
       numDomains = self%numDomains
       return
    end subroutine Get_numDomains
!-----------------------------------------------------------------------------
    subroutine Set_numDomains(self, numDomains)
       implicit none
       integer, intent(in) :: numDomains
       type(t_gmiDomain), intent(inOut) :: self
       self%numDomains = numDomains
       return
    end subroutine Set_numDomains
!-----------------------------------------------------------------------------
    subroutine Get_rootProc(self, rootProc)
       implicit none
       integer, intent(out) :: rootProc
       type(t_gmiDomain), intent(in) :: self
       rootProc = self%rootProc
       return
    end subroutine Get_rootProc
!-----------------------------------------------------------------------------
    subroutine Set_rootProc(self, rootProc)
       implicit none
       integer, intent(in) :: rootProc
       type(t_gmiDomain), intent(inOut) :: self
       self%rootProc = rootProc
       return
    end subroutine Set_rootProc
!-----------------------------------------------------------------------------
    subroutine Set_southBorder(self, southBorder)
       implicit none
       integer, intent(in) :: southBorder
       type(t_gmiDomain), intent(inOut) :: self
       self%southBorder = southBorder
       return
    end subroutine Set_southBorder
!-----------------------------------------------------------------------------
    subroutine Set_northBorder(self, northBorder)
       implicit none
       integer, intent(in) :: northBorder
       type(t_gmiDomain), intent(inOut) :: self
       self%northBorder = northBorder
       return
    end subroutine Set_northBorder
!-----------------------------------------------------------------------------
    subroutine Set_westBorder(self, westBorder)
       implicit none
       integer, intent(in) :: westBorder
       type(t_gmiDomain), intent(inOut) :: self
       self%westBorder = westBorder
       return
    end subroutine Set_westBorder
!-----------------------------------------------------------------------------
    subroutine Set_eastBorder(self, eastBorder)
       implicit none
       integer, intent(in) :: eastBorder
       type(t_gmiDomain), intent(inOut) :: self
       self%eastBorder = eastBorder
       return
    end subroutine Set_eastBorder
!-----------------------------------------------------------------------------
    subroutine Set_southNeighbor(self, southNeighbor)
       implicit none
       integer, intent(in) :: southNeighbor
       type(t_gmiDomain), intent(inOut) :: self
       self%southNeighbor = southNeighbor
       return
    end subroutine Set_southNeighbor
!-----------------------------------------------------------------------------
    subroutine Set_northNeighbor(self, northNeighbor)
       implicit none
       integer, intent(in) :: northNeighbor
       type(t_gmiDomain), intent(inOut) :: self
       self%northNeighbor = northNeighbor
       return
    end subroutine Set_northNeighbor
!-----------------------------------------------------------------------------
    subroutine Set_eastNeighbor(self, eastNeighbor)
       implicit none
       integer, intent(in) :: eastNeighbor
       type(t_gmiDomain), intent(inOut) :: self
       self%eastNeighbor = eastNeighbor
       return
    end subroutine Set_eastNeighbor
!-----------------------------------------------------------------------------
    subroutine Set_westNeighbor(self, westNeighbor)
       implicit none
       integer, intent(in) :: westNeighbor
       type(t_gmiDomain), intent(inOut) :: self
       self%westNeighbor = westNeighbor
       return
    end subroutine Set_westNeighbor
!-----------------------------------------------------------------------------
    subroutine Set_mostWesternDomain(self, mostWesternDomain)
       implicit none
       integer, intent(in) :: mostWesternDomain
       type(t_gmiDomain), intent(inOut) :: self
       self%mostWesternDomain = mostWesternDomain
       return
    end subroutine Set_mostWesternDomain
!-----------------------------------------------------------------------------
    subroutine Set_mostEasternDomain(self, mostEasternDomain)
       implicit none
       integer, intent(in) :: mostEasternDomain
       type(t_gmiDomain), intent(inOut) :: self
       self%mostEasternDomain = mostEasternDomain
       return
    end subroutine Set_mostEasternDomain
!-----------------------------------------------------------------------------
    subroutine Allocate_globMostWesternDomain(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%globMostWesternDomain(self%numDomains))
       self%globMostWesternDomain = 0
       return
    end subroutine Allocate_globMostWesternDomain
!-----------------------------------------------------------------------------
    subroutine Get_globMostWesternDomain(self, globMostWesternDomain)
       implicit none
       integer, intent(out) :: globMostWesternDomain(:)
       type(t_gmiDomain), intent(in) :: self
       globMostWesternDomain(:) = self%globMostWesternDomain(:)
       return
    end subroutine Get_globMostWesternDomain
!-----------------------------------------------------------------------------
    subroutine Set_globMostWesternDomain(self, globMostWesternDomain)
       implicit none
       integer, intent(in) :: globMostWesternDomain(:)
       type(t_gmiDomain), intent(inOut) :: self
       self%globMostWesternDomain(:) = globMostWesternDomain(:)
       return
    end subroutine Set_globMostWesternDomain
!-----------------------------------------------------------------------------
    subroutine Allocate_globMostEasternDomain(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%globMostEasternDomain(self%numDomains))
       self%globMostEasternDomain = 0
       return
    end subroutine Allocate_globMostEasternDomain
!-----------------------------------------------------------------------------
    subroutine Get_globMostEasternDomain(self, globMostEasternDomain)
       implicit none
       integer, intent(out) :: globMostEasternDomain(:)
       type(t_gmiDomain), intent(in) :: self
       globMostEasternDomain(:) = self%globMostEasternDomain(:)
       return
    end subroutine Get_globMostEasternDomain
!-----------------------------------------------------------------------------
    subroutine Set_globMostEasternDomain(self, globMostEasternDomain)
       implicit none
       integer, intent(in) :: globMostEasternDomain(:)
       type(t_gmiDomain), intent(inOut) :: self
       self%globMostEasternDomain(:) = globMostEasternDomain(:)
       return
    end subroutine Set_globMostEasternDomain
!-----------------------------------------------------------------------------
    subroutine Allocate_globSouthBorder(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%globSouthBorder(self%numDomains))
       self%globSouthBorder = 0
       return
    end subroutine Allocate_globSouthBorder
!-----------------------------------------------------------------------------
    subroutine Get_globSouthBorder(self, globSouthBorder)
       implicit none
       integer, intent(out) :: globSouthBorder(:)
       type(t_gmiDomain), intent(in) :: self
       globSouthBorder(:) = self%globSouthBorder(:)
       return
    end subroutine Get_globSouthBorder
!-----------------------------------------------------------------------------
    subroutine Set_globSouthBorder(self, globSouthBorder)
       implicit none
       integer, intent(in) :: globSouthBorder(:)
       type(t_gmiDomain), intent(inOut) :: self
       self%globSouthBorder(:) = globSouthBorder(:)
       return
    end subroutine Set_globSouthBorder
!-----------------------------------------------------------------------------
    subroutine Allocate_globNorthBorder(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%globNorthBorder(self%numDomains))
       self%globNorthBorder = 0
       return
    end subroutine Allocate_globNorthBorder
!-----------------------------------------------------------------------------
    subroutine Get_globNorthBorder(self, globNorthBorder)
       implicit none
       integer, intent(out) :: globNorthBorder(:)
       type(t_gmiDomain), intent(in) :: self
       globNorthBorder(:) = self%globNorthBorder(:)
       return
    end subroutine Get_globNorthBorder
!-----------------------------------------------------------------------------
    subroutine Set_globNorthBorder(self, globNorthBorder)
       implicit none
       integer, intent(in) :: globNorthBorder(:)
       type(t_gmiDomain), intent(inOut) :: self
       self%globNorthBorder(:) = globNorthBorder(:)
       return
    end subroutine Set_globNorthBorder
!-----------------------------------------------------------------------------
    subroutine Allocate_globWestBorder(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%globWestBorder(self%numDomains))
       self%globWestBorder = 0
       return
    end subroutine Allocate_globWestBorder
!-----------------------------------------------------------------------------
    subroutine Get_globWestBorder(self, globWestBorder)
       implicit none
       integer, intent(out) :: globWestBorder(:)
       type(t_gmiDomain), intent(in) :: self
       globWestBorder(:) = self%globWestBorder(:)
       return
    end subroutine Get_globWestBorder
!-----------------------------------------------------------------------------
    subroutine Set_globWestBorder(self, globWestBorder)
       implicit none
       integer, intent(in) :: globWestBorder(:)
       type(t_gmiDomain), intent(inOut) :: self
       self%globWestBorder(:) = globWestBorder(:)
       return
    end subroutine Set_globWestBorder
!-----------------------------------------------------------------------------
    subroutine Allocate_globEastBorder(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%globEastBorder(self%numDomains))
       self%globEastBorder = 0
       return
    end subroutine Allocate_globEastBorder
!-----------------------------------------------------------------------------
    subroutine Get_globEastBorder(self, globEastBorder)
       implicit none
       integer, intent(out) :: globEastBorder(:)
       type(t_gmiDomain), intent(in) :: self
       globEastBorder(:) = self%globEastBorder(:)
       return
    end subroutine Get_globEastBorder
!-----------------------------------------------------------------------------
    subroutine Set_globEastBorder(self, globEastBorder)
       implicit none
       integer, intent(in) :: globEastBorder(:)
       type(t_gmiDomain), intent(inOut) :: self
       self%globEastBorder(:) = globEastBorder(:)
       return
    end subroutine Set_globEastBorder
!-----------------------------------------------------------------------------
    subroutine Allocate_globSouthDomain(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%globSouthDomain(self%numDomains))
       self%globSouthDomain = 0
       return
    end subroutine Allocate_globSouthDomain
!-----------------------------------------------------------------------------
    subroutine Get_globSouthDomain(self, globSouthDomain)
       implicit none
       integer, intent(out) :: globSouthDomain(:)
       type(t_gmiDomain), intent(in) :: self
       globSouthDomain(:) = self%globSouthDomain(:)
       return
    end subroutine Get_globSouthDomain
!-----------------------------------------------------------------------------
    subroutine Set_globSouthDomain(self, globSouthDomain)
       implicit none
       integer, intent(in) :: globSouthDomain(:)
       type(t_gmiDomain), intent(inOut) :: self
       self%globSouthDomain(:) = globSouthDomain(:)
       return
    end subroutine Set_globSouthDomain
!-----------------------------------------------------------------------------
    subroutine Allocate_globNorthDomain(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%globNorthDomain(self%numDomains))
       self%globNorthDomain = 0
       return
    end subroutine Allocate_globNorthDomain
!-----------------------------------------------------------------------------
    subroutine Get_globNorthDomain(self, globNorthDomain)
       implicit none
       integer, intent(out) :: globNorthDomain(:)
       type(t_gmiDomain), intent(in) :: self
       globNorthDomain(:) = self%globNorthDomain(:)
       return
    end subroutine Get_globNorthDomain
!-----------------------------------------------------------------------------
    subroutine Set_globNorthDomain(self, globNorthDomain)
       implicit none
       integer, intent(in) :: globNorthDomain(:)
       type(t_gmiDomain), intent(inOut) :: self
       self%globNorthDomain(:) = globNorthDomain(:)
       return
    end subroutine Set_globNorthDomain
!-----------------------------------------------------------------------------
    subroutine Allocate_globWestDomain(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%globWestDomain(self%numDomains))
       self%globWestDomain = 0
       return
    end subroutine Allocate_globWestDomain
!-----------------------------------------------------------------------------
    subroutine Get_globWestDomain(self, globWestDomain)
       implicit none
       integer, intent(out) :: globWestDomain(:)
       type(t_gmiDomain), intent(in) :: self
       globWestDomain(:) = self%globWestDomain(:)
       return
    end subroutine Get_globWestDomain
!-----------------------------------------------------------------------------
    subroutine Set_globWestDomain(self, globWestDomain)
       implicit none
       integer, intent(in) :: globWestDomain(:)
       type(t_gmiDomain), intent(inOut) :: self
       self%globWestDomain(:) = globWestDomain(:)
       return
    end subroutine Set_globWestDomain
!-----------------------------------------------------------------------------
    subroutine Allocate_globEastDomain(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%globEastDomain(self%numDomains))
       self%globEastDomain = 0
       return
    end subroutine Allocate_globEastDomain
!-----------------------------------------------------------------------------
    subroutine Get_globEastDomain(self, globEastDomain)
       implicit none
       integer, intent(out) :: globEastDomain(:)
       type(t_gmiDomain), intent(in) :: self
       globEastDomain(:) = self%globEastDomain(:)
       return
    end subroutine Get_globEastDomain
!-----------------------------------------------------------------------------
    subroutine Set_globEastDomain(self, globEastDomain)
       implicit none
       integer, intent(in) :: globEastDomain(:)
       type(t_gmiDomain), intent(inOut) :: self
       self%globEastDomain(:) = globEastDomain(:)
       return
    end subroutine Set_globEastDomain
!-----------------------------------------------------------------------------
    subroutine Allocate_eastWestFlag(self, ilo_gl, ihi_gl)
       implicit none
       integer, intent(in) :: ilo_gl, ihi_gl
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%eastWestFlag(ilo_gl:ihi_gl))
       self%eastWestFlag = 0
       return
    end subroutine Allocate_eastWestFlag
!-----------------------------------------------------------------------------
    subroutine Get_eastWestFlag(self, eastWestFlag)
       implicit none
       integer, intent(out) :: eastWestFlag(:)
       type(t_gmiDomain), intent(in) :: self
       eastWestFlag(:) = self%eastWestFlag(:)
       return
    end subroutine Get_eastWestFlag
!-----------------------------------------------------------------------------
    subroutine Set_eastWestFlag(self, eastWestFlag)
       implicit none
       integer, intent(in) :: eastWestFlag(:)
       type(t_gmiDomain), intent(inOut) :: self
       self%eastWestFlag(:) = eastWestFlag(:)
       return
    end subroutine Set_eastWestFlag
!-----------------------------------------------------------------------------
    subroutine Allocate_northSouthPole(self, julo_gl, jhi_gl)
       implicit none
       integer, intent(in) :: julo_gl, jhi_gl
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%northSouthPole(julo_gl:jhi_gl))
       self%northSouthPole = 0
       return
    end subroutine Allocate_northSouthPole
!-----------------------------------------------------------------------------
    subroutine Get_northSouthPole(self, northSouthPole)
       implicit none
       integer, intent(out) :: northSouthPole(:)
       type(t_gmiDomain), intent(in) :: self
       northSouthPole(:) = self%northSouthPole(:)
       return
    end subroutine Get_northSouthPole
!-----------------------------------------------------------------------------
    subroutine Set_northSouthPole(self, northSouthPole)
       implicit none
       integer, intent(in) :: northSouthPole(:)
       type(t_gmiDomain), intent(inOut) :: self
       self%northSouthPole(:) = northSouthPole(:)
       return
    end subroutine Set_northSouthPole
!-----------------------------------------------------------------------------
    subroutine Allocate_mapi_all(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%mapi_all(2, self%numDomains))
       self%mapi_all = 0
       return
    end subroutine Allocate_mapi_all
!-----------------------------------------------------------------------------
    subroutine Get_mapi_all(self, mapi_all)
       implicit none
       integer, intent(out) :: mapi_all(:,:)
       type(t_gmiDomain), intent(in) :: self
       mapi_all(:,:) = self%mapi_all(:,:)
       return
    end subroutine Get_mapi_all
!-----------------------------------------------------------------------------
    subroutine Set_mapi_all(self, mapi_all)
       implicit none
       integer, intent(in) :: mapi_all(:,:)
       type(t_gmiDomain), intent(inOut) :: self
       self%mapi_all(:,:) = mapi_all(:,:)
       return
    end subroutine Set_mapi_all
!-----------------------------------------------------------------------------
    subroutine Allocate_map1_u(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%map1_u(2, 2, self%numDomains))
       self%map1_u = 0
       return
    end subroutine Allocate_map1_u
!-----------------------------------------------------------------------------
    subroutine Get_map1_u(self, map1_u)
       implicit none
       integer, intent(out) :: map1_u(:,:,:)
       type(t_gmiDomain), intent(in) :: self
       map1_u(:,:,:) = self%map1_u(:,:,:)
       return
    end subroutine Get_map1_u
!-----------------------------------------------------------------------------
    subroutine Set_map1_u(self, map1_u)
       implicit none
       integer, intent(in) :: map1_u(:,:,:)
       type(t_gmiDomain), intent(inOut) :: self
       self%map1_u(:,:,:) = map1_u(:,:,:)
       return
    end subroutine Set_map1_u
!-----------------------------------------------------------------------------
    subroutine Allocate_map2_u(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%map2_u(2, 2, self%numDomains))
       self%map2_u = 0
       return
    end subroutine Allocate_map2_u
!-----------------------------------------------------------------------------
    subroutine Get_map2_u(self, map2_u)
       implicit none
       integer, intent(out) :: map2_u(:,:,:)
       type(t_gmiDomain), intent(in) :: self
       map2_u(:,:,:) = self%map2_u(:,:,:)
       return
    end subroutine Get_map2_u
!-----------------------------------------------------------------------------
    subroutine Set_map2_u(self, map2_u)
       implicit none
       integer, intent(in) :: map2_u(:,:,:)
       type(t_gmiDomain), intent(inOut) :: self
       self%map2_u(:,:,:) = map2_u(:,:,:)
       return
    end subroutine Set_map2_u
!-----------------------------------------------------------------------------
    subroutine Allocate_map1_v(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%map1_v(2, 2, self%numDomains))
       self%map1_v = 0
       return
    end subroutine Allocate_map1_v
!-----------------------------------------------------------------------------
    subroutine Get_map1_v(self, map1_v)
       implicit none
       integer, intent(out) :: map1_v(:,:,:)
       type(t_gmiDomain), intent(in) :: self
       map1_v(:,:,:) = self%map1_v(:,:,:)
       return
    end subroutine Get_map1_v
!-----------------------------------------------------------------------------
    subroutine Set_map1_v(self, map1_v)
       implicit none
       integer, intent(in) :: map1_v(:,:,:)
       type(t_gmiDomain), intent(inOut) :: self
       self%map1_v(:,:,:) = map1_v(:,:,:)
       return
    end subroutine Set_map1_v
!-----------------------------------------------------------------------------
    subroutine Allocate_map2_v(self)
       implicit none
       type(t_gmiDomain), intent(inOut) :: self
       allocate(self%map2_v(2, 2, self%numDomains))
       self%map2_v = 0
       return
    end subroutine Allocate_map2_v
!-----------------------------------------------------------------------------
    subroutine Get_map2_v(self, map2_v)
       implicit none
       integer, intent(out) :: map2_v(:,:,:)
       type(t_gmiDomain), intent(in) :: self
       map2_v(:,:,:) = self%map2_v(:,:,:)
       return
    end subroutine Get_map2_v
!-----------------------------------------------------------------------------
    subroutine Set_map2_v(self, map2_v)
       implicit none
       integer, intent(in) :: map2_v(:,:,:)
       type(t_gmiDomain), intent(inOut) :: self
       self%map2_v(:,:,:) = map2_v(:,:,:)
       return
    end subroutine Set_map2_v
!-----------------------------------------------------------------------------
    subroutine Get_coscen(self, coscen)
       implicit none
       real*8, intent(out) :: coscen(:)
       type(t_gmiDomain), intent(in) :: self
       coscen(:) = self%coscen(:)
       return
    end subroutine Get_coscen
!-----------------------------------------------------------------------------
    subroutine Get_cosedg(self, cosedg)
       implicit none
       real*8, intent(out) :: cosedg(:)
       type(t_gmiDomain), intent(in) :: self
       cosedg(:) = self%cosedg(:)
       return
    end subroutine Get_cosedg
!-----------------------------------------------------------------------------
    subroutine Get_cose(self, cose)
       implicit none
       real*8, intent(out) :: cose(:)
       type(t_gmiDomain), intent(in) :: self
       cose(:) = self%cose(:)
       return
    end subroutine Get_cose
!-----------------------------------------------------------------------------
    subroutine Get_cosp(self, cosp)
       implicit none
       real*8, intent(out) :: cosp(:)
       type(t_gmiDomain), intent(in) :: self
       cosp(:) = self%cosp(:)
       return
    end subroutine Get_cosp
!-----------------------------------------------------------------------------
    subroutine Get_londeg(self, londeg)
       implicit none
       real*8, intent(out) :: londeg(:)
       type(t_gmiDomain), intent(in) :: self
       londeg(:) = self%londeg(:)
       return
    end subroutine Get_londeg
!-----------------------------------------------------------------------------
    subroutine Get_latdeg(self, latdeg)
       implicit none
       real*8, intent(out) :: latdeg(:)
       type(t_gmiDomain), intent(in) :: self
       latdeg(:) = self%latdeg(:)
       return
    end subroutine Get_latdeg
!-----------------------------------------------------------------------------
    subroutine Get_geofac(self, geofac)
       implicit none
       real*8, intent(out) :: geofac(:)
       type(t_gmiDomain), intent(in) :: self
       geofac(:) = self%geofac(:)
       return
    end subroutine Get_geofac
!-----------------------------------------------------------------------------
    subroutine Get_dlatr(self, dlatr)
       implicit none
       real*8, intent(out) :: dlatr(:)
       type(t_gmiDomain), intent(in) :: self
       dlatr(:) = self%dlatr(:)
       return
    end subroutine Get_dlatr
!-----------------------------------------------------------------------------
    subroutine Get_mcorGlob(self, mcorGlob)
       implicit none
       real*8, intent(out) :: mcorGlob(:,:)
       type(t_gmiDomain), intent(in) :: self
       mcorGlob(:,:) = self%mcorGlob(:,:)
       return
    end subroutine Get_mcorGlob
!-----------------------------------------------------------------------------
    subroutine Get_rel_areaGlob(self, rel_areaGlob)
       implicit none
       real*8, intent(out) :: rel_areaGlob(:,:)
       type(t_gmiDomain), intent(in) :: self
       rel_areaGlob(:,:) = self%rel_areaGlob(:,:)
       return
    end subroutine Get_rel_areaGlob
!-----------------------------------------------------------------------------
    subroutine Get_mcor(self, mcor)
       implicit none
       real*8, intent(out) :: mcor(:,:)
       type(t_gmiDomain), intent(in) :: self
       mcor(:,:) = self%mcor(:,:)
       return
    end subroutine Get_mcor
!-----------------------------------------------------------------------------
    subroutine Get_rel_area(self, rel_area)
       implicit none
       real*8, intent(out) :: rel_area(:,:)
       type(t_gmiDomain), intent(in) :: self
       rel_area(:,:) = self%rel_area(:,:)
       return
    end subroutine Get_rel_area
!-----------------------------------------------------------------------------
    subroutine Get_geofac_pc(self, geofac_pc)
       implicit none
       real*8, intent(out) :: geofac_pc
       type(t_gmiDomain), intent(in) :: self
       geofac_pc = self%geofac_pc
       return
    end subroutine Get_geofac_pc
!-----------------------------------------------------------------------------
!
      end module GmiDomainDecomposition_mod
!
