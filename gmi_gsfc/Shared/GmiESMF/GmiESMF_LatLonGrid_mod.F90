#include "GmiESMF_ErrLog.h"
!------------------------------------------------------------------------------
! NASA GSFC - SSSO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiESMF_LatLonGrid_mod
!
! !INTERFACE:
!
      module GmiESMF_LatLonGrid_mod
!
! !USES:
      USE ESMF_Mod
      USE GmiESMF_ErrorChecking_mod
!
      IMPLICIT NONE
!
! !PUBLIC DATA MEMBERS:
      PRIVATE
      PUBLIC  :: getESMF_GridInterior
      PUBLIC  :: createESMF_LatLonGrid
!
! !DESCRIPTION:
! Routines/Functions to create the ESMF lat/lon grid.
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: createESMF_LatLonGrid
!
! !INTERFACE:
!
      FUNCTION createESMF_LatLonGrid (Name, Nx, Ny, &
                     IM_World, BegLon, DelLon,      &
                     JM_World, BegLat, DelLat,      &
                     LM_World, rc )                       RESULT(Grid)
!
# include "gmi_phys_constants.h"
!
! !INPUT PARAMETERS:
      character(len=*), intent(in)  :: Name ! grid type (A-grid, etc.)

      integer,          intent(in)  :: Nx, Ny          ! Layout
      integer,          intent(in)  :: IM_World        ! Zonal 
      real,             intent(in)  :: BegLon, DelLon  ! in degrees

      integer,          intent(in)  :: JM_World        ! Meridional
      real,             intent(in)  :: BegLat, DelLat  ! in degrees

      integer,          intent(in)  :: LM_World        ! Vertical
!
! !OUTPUT PARAMETERS:
      type (ESMF_Grid)                        :: Grid  ! Distributed grid
      integer,          OPTIONAL, intent(out) :: rc    ! return code
!
! !DESCRIPTION:
!  Creates a distributed ESMF grid where the horizontal
!  coordinates are regular longitudes and latitudes. The grid is
!  created on the user specified {\bf VM}, or on the current VM if the user
!  does not specify one. The layout and the coordinate information can
!  be provided with a {\tt ESMF\_Config attribute}, a resource file name
!  or specified through the argument list.
!
!EOP
!-------------------------------------------------------------------------------
!BOC

!   Internal version of the input arguments
!   ---------------------------------------
      type(ESMF_Config), pointer :: Config_
      integer           :: IM_World_
      real              :: BegLon_
      real              :: DelLon_
      integer           :: JM_World_
      real              :: BegLat_
      real              :: DelLat_
      integer           :: LM_World_
      integer           :: Nx_, Ny_, Nz_

      integer, allocatable            :: IMs(:), JMs(:), LMs(:)
      real(ESMF_KIND_R8)              :: minCoord(3)
      real(ESMF_KIND_R8)              :: deltaX, deltaY
      type (ESMF_VM), pointer         :: VM_
      integer                         :: I, J, I1, IN, J1, JN

      real(ESMF_KIND_R8), pointer     :: centerX(:,:)
      real(ESMF_KIND_R8), pointer     :: centerY(:,:)
      real(ESMF_KIND_R8), allocatable :: cornerX(:)
      real(ESMF_KIND_R8), allocatable :: cornerY(:)

      real, parameter                 :: D2R = GMI_PI / 180.

      integer                         :: STATUS
      character(len=ESMF_MAXSTR)      :: IAm='createESMF_LatLonGrid'

!                                ------

!  Initial Settings
!  --------
      Nz_     =  1      ! place holder for now

      IM_World_ = IM_World
      JM_World_ = JM_World
      LM_World_ = LM_World

      Nx_ = Nx
      Ny_ = Ny

      BegLon_ = BegLon
      DelLon_ = DelLon
      BegLat_ = BegLat
      DelLat_ = DelLat

!  Global grids
!  ------------
      if ( IM_World_ < 1 .OR. JM_World_ < 1 ) then
           STATUS = 400
           VERIFY_(STATUS)
      end if
      if ( DelLon_ < 0.0 ) then  ! convention for global grids
         if ( IM_World_ == 1 ) then
           DelLon_ = 0.0
         else
           DelLon_ = 360. / IM_World_
         end if
      end if
      if ( DelLat_ < 0.0 ) then  ! convention for global grids
         if ( JM_World_ == 1 ) then
              DelLat_ = 0.0
         else
              DelLat_ = 180. / ( JM_World_ - 1)
         end if
      end if

!  Give the IMs, JMs and LMs the default distribution
!  -------------------------------------------------------
      allocate( IMs(0:Nx_-1), JMs(0:Ny_-1), LMs(0:Nz_-1), stat=STATUS)
      VERIFY_(STATUS)
      call decomposeDim ( IM_World_, IMs, Nx_ )
      call decomposeDim ( JM_World_, JMs, Ny_ )
      call decomposeDim ( LM_World_, LMs, Nz_ )

      Grid = ESMF_GridCreateShapeTile(    &
           name=Name,                 &
           countsPerDEDim1=IMs,           &
           countsPerDEDim2=JMs,           &
           indexFlag = ESMF_INDEX_USER,   &
           gridMemLBound = (/1,1/),       &
           gridEdgeLWidth = (/0,0/),      &
           gridEdgeUWidth = (/0,0/),      &
           coordDep1 = (/1,2/),           &
           coordDep2 = (/1,2/),           &
           rc=status)
      VERIFY_(STATUS)

      call ESMF_AttributeSet(grid, name='GRID_LM', value=LM_World, rc=status)
      VERIFY_(STATUS)

!  -------------------------------------------------------------------
!  NOTE: In the remaining part of this routine it is assumed that the 
!        1st and 2nd axes correspond to lat/lon; revise this for other 
!        arrangements (say, YZ grids)
!  -------------------------------------------------------------------

!  Allocate coords at default stagger location
!  -------------------------------------------
      call ESMF_GridAddCoord(Grid, rc=status)
      VERIFY_(STATUS)

!  Compute the coordinates (the corner/center is for backward compatibility)
!  -------------------------------------------------------------------------
      deltaX      = D2R * DelLon_
      deltaY      = D2R * DelLat_
      minCoord(1) = D2R * BegLon_ - deltaX/2
      minCoord(2) = D2R * BegLat_ - deltaY/2

      allocate(cornerX(IM_World_+1),cornerY(JM_World_+1), stat=STATUS)
      VERIFY_(STATUS)

      cornerX(1) = minCoord(1)
      do i = 1,IM_World_
         cornerX(i+1) = cornerX(i) + deltaX
      enddo

      cornerY(1) = minCoord(2)
      do j = 1,JM_World_
         cornerY(j+1) = cornerY(j) + deltaY
      enddo

!  Retrieve the coordinates so we can set them
!  -------------------------------------------
      call ESMF_GridGetCoord (Grid, coordDim=1, localDE=0, &
                           staggerloc=ESMF_STAGGERLOC_CENTER, &
                           fptr=centerX, rc=status)
      VERIFY_(STATUS)

      call ESMF_GridGetCoord (Grid, coordDim=2, localDE=0, &
                           staggerloc=ESMF_STAGGERLOC_CENTER, &
                           fptr=centerY, rc=status)
      VERIFY_(STATUS)

      call getESMF_GridInterior (Grid,i1,in,j1,jn)

      do i = 1,size(centerX,1)
         centerX(i,:) = 0.5d0*(cornerX(i+i1-1)+cornerX(i+i1))
      end do

      do j = 1,size(centerY,2)
         centerY(:,j) = 0.5d0*(cornerY(j+j1-1)+cornerY(j+j1))
      enddo


!  Make sure we've got it right
!  ----------------------------
      call ESMF_GridValidate(Grid,rc=status)
      VERIFY_(STATUS)

!  Clean up
!  --------   
      deallocate(cornerY,cornerX)
      deallocate(IMs,JMs,LMs)
   
      end function createESMF_LatLonGrid
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: decomposeDim
!
! !INTERFACE:
!
      subroutine decomposeDim ( dim_world,dim,NDEs )
      implicit   none
      integer    dim_world, NDEs
      integer    dim(0:NDEs-1)
      integer    n,im,rm,nbeg,nend
!EOP
!-------------------------------------------------------------------------------
!BOC
      im = dim_world/NDEs
      rm = dim_world-NDEs*im
      do n=0,NDEs-1
                      dim(n) = im
      if( n.le.rm-1 ) dim(n) = im+1
      enddo
      RETURN
      end subroutine decomposeDim
!EOC
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getESMF_GridInterior
!
! !INTERFACE:
!
      subroutine getESMF_GridInterior(GRID, I1, IN, J1, JN)
!
! !INPUT PARAMETERS:
      type (ESMF_Grid), intent(IN) :: grid
!
! !OUTPUT PARAMETERS:
      integer, intent(OUT)         :: I1   ! left  lon point
      integer, intent(OUT)         :: IN   ! right lon point
      integer, intent(OUT)         :: J1   ! left  lat point
      integer, intent(OUT)         :: JN   ! right lat point
!
! !DESCRIPTION:
! Determines the interior grid points assigned to the processor.
!
!EOP
!-------------------------------------------------------------------------------
!BOC
! local vars
      integer                               :: status
      character(len=ESMF_MAXSTR)            :: IAm='getESMF_GridInterior'

      type (ESMF_DistGrid)                  :: distGrid
      type(ESMF_DELayout)                   :: LAYOUT
      integer,               allocatable    :: AL(:,:)
      integer,               allocatable    :: AU(:,:)
      integer                               :: nDEs
      integer                               :: deId
      integer                               :: gridRank
      integer                               :: deList(1)

      call ESMF_GridGet    (GRID, dimCount=gridRank, distGrid=distGrid, rc=STATUS)
      call ESMF_DistGridGet(distGRID, delayout=layout, rc=STATUS)
      call ESMF_DELayoutGet(layout, deCount =nDEs, localDeList=deList, rc=status)
      deId = deList(1)

      allocate (AL(gridRank,0:nDEs-1),  stat=status)
      allocate (AU(gridRank,0:nDEs-1),  stat=status)

      call ESMF_DistGridGet(distgrid, &
           minIndexPDimPDe=AL, maxIndexPDimPDe=AU, rc=status)

      I1 = AL(1, deId)
      IN = AU(1, deId)
      J1 = AL(2, deId)
      JN = AU(2, deId)
      deallocate(AU, AL)

      RETURN

      end subroutine getESMF_GridInterior
!EOC
!-------------------------------------------------------------------------


      end module GmiESMF_LatLonGrid_mod
