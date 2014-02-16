!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiESMFgrid_mod
!
! !INTERFACE:
!
#include "GmiESMF_ErrLog.h"
!
      module GmiESMFgrid_mod
!
! !USES:
      use GmiTimeControl_mod, only : t_GmiClock

      ! ESMF module, defines all ESMF data types and procedures
      use ESMF_Mod
      use GmiESMF_ErrorChecking_mod
      use GmiESMF_LatLonGrid_mod

      implicit none

! !PUBLIC DATA MEMBERS:
      private
      public  :: createESMFgrid
!
! !DESCRIPTION:
! Routines for creating the ESMF grid.
!
! !AUTHOR:
!  Jules Kouatchou, NASA/GSFC, Jules.Kouatchou@nasa.gov
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: createESMFgrid
!
! !INTERFACE:
!
      function createESMFgrid(config, vm, rc) result(grid)
!
      implicit none
!
#include "gmi_phys_constants.h"
!
! !INPUT PARAMETERS:
      type (ESMF_VM), intent(in) :: VM
!
! !OUTPUT PARAMETERS:
      integer, optional, intent(out) :: rc
!
! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Config), intent(inOut) :: config
!
! !RETURN VALUE:
      type(ESMF_Grid)     :: grid
!
! !DESCRIPTION:
! Create the ESMF grid.
!
! !LOCAL VARIABLES:
      integer                         :: STATUS
      integer                         :: IM_GLOBAL, JM_GLOBAL, LM_GLOBAL
      integer                         :: NX, NY, NZ
      integer                         :: LATLON, L, POLEEDGE
      character(len=ESMF_MAXSTR)      :: GRIDNAME
      real              :: deltaX, deltaY, deltaZ
      real              :: minCoord(3)
      real              :: xmin, xmax, ymin, ymax
      real                            :: LON0, LAT0
      character(len=ESMF_MAXSTR), parameter :: IAm = "createESMFgrid"
      type(ESMF_DELayout) :: layout
!EOP
!------------------------------------------------------------------------------
!BOC

! Get Decomposition from CF
!--------------------------

! !RESOURCE_ITEM: none :: Processing elements in 1st dimension
      call ESMF_ConfigGetAttribute(config, NX, label='NX:', default=1, rc=STATUS)
      VERIFY_(STATUS)

! !RESOURCE_ITEM: none :: Processing elements in 2nd dimension
      call ESMF_ConfigGetAttribute(config, NY, label='NY:', default=1, rc=STATUS)
      VERIFY_(STATUS)

! Get World problem size from CF
!-------------------------------

! !RESOURCE_ITEM: none :: Grid size in 1st dimension
      call ESMF_ConfigGetAttribute(config, IM_GLOBAL, label='IM:', rc=STATUS )
      VERIFY_(STATUS)

! !RESOURCE_ITEM: none :: Grid size in 2nd dimension
      call ESMF_ConfigGetAttribute(config, JM_GLOBAL, label='JM:', rc=STATUS )
      VERIFY_(STATUS)

! !RESOURCE_ITEM: none :: Grid size in 3rd dimension
      call ESMF_ConfigGetAttribute(config, LM_GLOBAL, label='LM:', default=1,  &
     &         rc=STATUS )
      VERIFY_(STATUS)

! !RESOURCE_ITEM: 0 or 1 :: 1->gridedge at pole; 0->gridpoint at pole
      call ESMF_ConfigGetAttribute(config, POLEEDGE, label='POLEEDGE:',        &
     &         default=0, rc=STATUS )
      VERIFY_(STATUS)

! !RESOURCE_ITEM: degrees :: Longituce of center of first gridbox
      call ESMF_ConfigGetAttribute(config, LON0, label='BEGLON:',             &
!     &         default=0., rc=STATUS )
     &         default=-90., rc=STATUS )
      VERIFY_(STATUS)

      LON0 = LON0 * (GMI_PI/180.)

! Lat-Lon Grid definition
!-------------------------

      deltaX      = 2.0*GMI_PI/IM_GLOBAL
      minCoord(1) = LON0-deltaX/2

      if (POLEEDGE==0) then
         deltaY = GMI_PI/(JM_GLOBAL-1)
         minCoord(2) = -GMI_PI/2-deltaY/2
      else
         deltaY = GMI_PI/(JM_GLOBAL  )
         minCoord(2) = -GMI_PI/2
      end if

      deltaZ = 1.0D0
      minCoord(3) = deltaZ/2

      xmin = 0.0
      xmax = 2.0*GMI_PI
      ymin = -GMI_PI/2.0
      ymax = +GMI_PI/2.0


      GRID = createESMF_LatLonGrid ("A-grid", NX , NY, &
                            IM_GLOBAL, xMin, deltaX, &
                            JM_GLOBAL, yMin, deltaY, &
                            LM_GLOBAL, rc=STATUS )
      VERIFY_(STATUS)

!      GRID = MAPL_LatLonGridCreate (      &
!                  Name = "A-grid",        &
!                  Nx     = NX ,           &
!                  Ny     = NY ,           &
!                  IM_World = IM_GLOBAL,   &
!                  BegLon   = xMin,        &
!                  DelLon   = deltaX,      &
!                  JM_World = JM_GLOBAL,   &
!                  BegLat   = yMin,        &
!                  DelLat   = deltaY,      &
!                  LM_World = LM_GLOBAL,   &
!                  rc=STATUS )
!      VERIFY_(STATUS)

      end function createESMFgrid
!EOC
!------------------------------------------------------------------------------

      end module GmiESMFgrid_mod
