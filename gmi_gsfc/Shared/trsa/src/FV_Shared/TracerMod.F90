! $Id: TracerMod.F90,v 1.2 2011-08-09 22:12:59 mrdamon Exp $
!
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: TracerMod - Container for T_Tracer type and its base methods
!
! !DESCRIPTION: 
!
!------------------------------------------------------------------------------
!
!
! !INTERFACE:

  module TracerMod

! !USES:

    use shr_kind_mod,       only : r8 => shr_kind_r8, r4 => shr_kind_r4
    implicit none
    private

! !PUBLIC TYPES:

    public T_Tracers, T_Tracer_Grid

! !PUBLIC MEMBER FUNCTIONS:
    public Tracer_GridInit, Tracer_GridFinal, Tracer_GridGet, Tracer_Set, Tracer_Get

!EOP

    type T_Tracers
       logical                                   :: is_r4
       real(r8), dimension(:,:,:  ), pointer     :: content
       real(r4), dimension(:,:,:  ), pointer     :: content_r4
    end type T_Tracers

 type T_Tracer_Grid

    integer :: myid_y = 0    ! subdomain index (0-based) in latitude (y)
    integer :: myid_z = 0    ! subdomain index (0 based) in level (z)
    integer :: npryz_y  = 1  ! number of subdomains in y
    integer :: npryz_z  = 1  ! number of subdomains in z
    integer :: iam = 0       !

    integer :: mod_method = 0  ! 1 for mpi derived types with transposes
                               ! 0 for contiguous buffers

    integer :: commglobal        ! global communicator
    integer :: comm_y            ! communicator in latitude
    integer :: comm_z            ! communicator in vertical

    integer :: jfirst            ! Start latitude (inclusive)
    integer :: jlast             ! End latitude (inclusive)
    integer :: kfirst            ! Start level (inclusive)
    integer :: klast             ! End level (inclusive)

!
    integer :: ng_d              ! Halo size

!
    real(r8)                        :: dlon
    real(r8)                        :: dphi
    real(r8)                        :: acap
    real(r8)                        :: rcap
!
    real(r8), dimension(:), pointer :: cosp  ! Cos of lat angle -- vol mean
    real(r8), dimension(:), pointer :: sinp  ! Sin of lat angle -- vol mean
    real(r8), dimension(:), pointer :: cose  ! Cos at finite volume edge
    real(r8), dimension(:), pointer :: sine  ! Sin at finite volume edge
    real(r8), dimension(:), pointer :: acosp ! Reciprocal of cos of lat angle

    real(r8), dimension(:), pointer :: dlat  ! Latitude array (1..jm)
    real(r8), dimension(:), pointer :: gw    ! Gaussian weights
 end type T_Tracer_Grid

 interface Tracer_Set
   module procedure Tracer_SetR4
   module procedure Tracer_SetR8
 end interface

 interface Tracer_Get
   module procedure Tracer_GetIsR4
   module procedure Tracer_GetR4
   module procedure Tracer_GetR8
 end interface


contains

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Tracer_GridInit -- Initialize the tracer grid
    subroutine Tracer_GridInit( ae, pi, jord, im, jm, km,           &
                                jfirstyz, jlastyz, kfirstyz, klastyz,   &
                                myidyz_y, myidyz_z, npryz_y, npryz_z,   &
                                grid )
! !USES:
      use parutilitiesmodule, only : gid, commglobal, parsplit
      implicit none

! !PARAMETERS:
      real(r8), intent(in)  :: ae                  !  Radius of earth
      real(r8), intent(in)  :: pi                  !  PI
      integer, intent(in)   :: jord                !  Order of scheme in Y
      integer, intent(in)   :: im, jm, km          !  Global dims

      integer, intent(in)   :: jfirstyz            !  First Y on this PE  (YZ)
      integer, intent(in)   :: jlastyz             !  Last Y on this PE (YZ)
      integer, intent(in)   :: kfirstyz            !  First Z on this PE (YZ)
      integer, intent(in)   :: klastyz             !  Last Z on this PE (YZ)

      integer, intent(in)   :: myidyz_y            !  YZ decomp - Nr in Y
      integer, intent(in)   :: myidyz_z            !  YZ decomp - Nr in Z

      integer, intent(in)   :: npryz_y             !  YZ decomp - Nr in Y
      integer, intent(in)   :: npryz_z             !  YZ decomp - Nr in Z
      type(T_TRACER_GRID), intent(inout)  :: grid  ! Resulting grid

! !DESCRIPTION:
!
!   Initialize Lin-Rood specific variables
!
! !REVISION HISTORY:
!
!   2006.06.27   Sawyer     Creation
!   2006.07.24   Sawyer     Extensively revised
!   2007.05.21   Sawyer     Revised again at request of Max Suarez
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
      integer  :: j
      real(r8) :: dl, dp, lat, ph5

      real(r8), parameter :: TWO  = 2.0_r8
      real(r8), parameter :: ONE  = 1.0_r8
      real(r8), parameter :: HALF = 0.5_r8
      real(r8), parameter :: ZERO = 0.0_r8

      integer   :: rank_y, rank_z
      integer   :: size_y, size_z

!  Set the basic grid variables

      grid%jfirst   = jfirstyz
      grid%jlast    = jlastyz
      grid%kfirst   = kfirstyz
      grid%klast    = klastyz

      allocate(grid%cosp(jm))
      allocate(grid%sinp(jm))
      allocate(grid%cose(jm))
      allocate(grid%sine(jm))
      allocate(grid%acosp(jm))

      allocate(grid%dlat(jm))
      allocate(grid%gw(jm))

      dl  = TWO*pi / im
      dp  = pi/(jm-1)
      grid%dlon = dl
      grid%dphi = dp
      lat = -HALF*pi
      ph5 = lat - HALF*dp
      do j=1,jm
         grid%dlat(j) = lat
         grid%sine(j) = sin(ph5)                         !  WS: should be OK
!!!         ph5  = -HALF*pi + ((j-1)-HALF)*(pi/(jm-1))      !  WS: should be OK
         ph5 = ph5 + dp
         lat = lat + dp
      enddo

      grid%cosp( 1) =  ZERO
      grid%cosp(jm) =  ZERO
      grid%gw( 1)   =  ZERO
      grid%gw(jm)   =  ZERO
      do j=2,jm-1
         grid%gw(j)   = grid%sine(j+1) - grid%sine(j)
         grid%cosp(j) = grid%gw(j) / dp
      enddo


! Define cosine at edges..

      do j=2,jm
         grid%cose(j) = HALF * (grid%cosp(j-1) + grid%cosp(j))
      enddo
      grid%cose(1) = grid%cose(2)

! Define cosine at edges..

      do j=2,jm
         grid%cose(j) = HALF * (grid%cosp(j-1) + grid%cosp(j))
      enddo
      grid%cose(1) = grid%cose(2)

      grid%sinp( 1) = -ONE
      grid%sinp(jm) =  ONE

      do j=2,jm-1
         grid%sinp(j) = HALF * (grid%sine(j) + grid%sine(j+1))
      enddo

!
! Pole cap area and inverse
      grid%acap = im*(ONE+grid%sine(2)) / dp
      grid%rcap = ONE / grid%acap

      if( mod(im,2) /= 0) then
         write(6,*) 'im must be an even integer'
         stop
      endif

      do j=2,jm-1
         grid%acosp(j) = ONE / grid%cosp(j)
      enddo
      grid%acosp( 1) = grid%rcap * im
      grid%acosp(jm) = grid%rcap * im

!
! Calculate the ghost region sizes for the SPMD version (tricky stuff)
!
      grid%ng_d = min( abs(jord), 3)   ! SJL: number of max ghost latitudes
      grid%ng_d = max( grid%ng_d, 2)

      grid%npryz_y  = npryz_y
      grid%npryz_z  = npryz_z

      grid%iam      = gid
      grid%myid_z   = myidyz_z
      grid%myid_y   = myidyz_y


! Split communicators

      grid%commglobal = commglobal
      call parsplit(commglobal, grid%myid_z, gid, grid%comm_y, rank_y, size_y)
      call parsplit(commglobal, grid%myid_y, gid, grid%comm_z, rank_z, size_z)

!EOC
    end subroutine Tracer_GridInit
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Tracer_SetR4 -- Set the tracer to R4 data
    subroutine Tracer_SetR4( tracer, array )

! !USES:
      implicit none

! !PARAMETERS:
      type(T_TRACERS), intent(inout)  :: tracer        ! Tracer
      real(r4), pointer               :: array(:,:,:)  ! Pointer to values

! !DESCRIPTION:
!
!   Set the tracer to real4 values
!
! !REVISION HISTORY:
!
!   2007.06.09   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
     tracer%is_r4  = .TRUE.
     tracer%content_r4 => array
     return
!EOC
    end subroutine Tracer_SetR4
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Tracer_SetR8 -- Set the tracer to R8 data
    subroutine Tracer_SetR8( tracer, array )

! !USES:
      implicit none

! !PARAMETERS:
      type(T_TRACERS), intent(inout)  :: tracer        ! Tracer
      real(r8), pointer               :: array(:,:,:)  ! Pointer to values

! !DESCRIPTION:
!
!   Set the tracer to real8 values
!
! !REVISION HISTORY:
!
!   2007.06.09   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
     tracer%is_r4  = .FALSE.
     tracer%content => array
     return
!EOC
    end subroutine Tracer_SetR8
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Tracer_GetIsR4 -- Get IsR4 variable
    subroutine Tracer_GetIsR4( tracer, isR4 )

! !USES:
      implicit none

! !PARAMETERS:
      type(T_TRACERS), intent(in)     :: tracer        ! Tracer
      logical                         :: isR4          ! Is tracer R4

! !DESCRIPTION:
!
!   Find out whether the tracer contains REAL4 data or not
!
! !REVISION HISTORY:
!
!   2007.06.09   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
     isR4 = tracer%is_r4
     return
!EOC
    end subroutine Tracer_GetIsR4
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Tracer_GetR4 -- Get the tracer R4 data
    subroutine Tracer_GetR4( tracer, array )

! !USES:
      implicit none

! !PARAMETERS:
      type(T_TRACERS), intent(in)     :: tracer        ! Tracer
      real(r4), pointer               :: array(:,:,:)  ! Pointer to values

! !DESCRIPTION:
!
!   Get the tracer real4 values
!
! !REVISION HISTORY:
!
!   2007.06.09   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
     if ( .not. tracer%is_r4 ) then
       print *, "Warning:  Retrieving tracer R4 data, but tracer is R8"
     endif
     array => tracer%content_r4
     return
!EOC
    end subroutine Tracer_GetR4
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Tracer_GetR8 -- Get the tracer R8 data
    subroutine Tracer_GetR8( tracer, array )

! !USES:
      implicit none

! !PARAMETERS:
      type(T_TRACERS), intent(in)     :: tracer        ! Tracer
      real(r8), pointer               :: array(:,:,:)  ! Pointer to values

! !DESCRIPTION:
!
!   Get the tracer real8 values
!
! !REVISION HISTORY:
!
!   2007.06.09   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
     if ( tracer%is_r4 ) then
       print *, "Warning:  Retrieving tracer R8 data, but tracer is R4"
     endif
     array => tracer%content
     return
!EOC
    end subroutine Tracer_GetR8
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Tracer_GridGet -- Retrieve values from the tracer grid
    subroutine Tracer_GridGet( grid, ng_d, acap, rcap, dlon, dphi,      &
                               comm_y, comm_z, commglobal,              &
                               jfirstyz, jlastyz, kfirstyz, klastyz,    &
                               myidyz_y, myidyz_z, npryz_y, npryz_z,    &
                               cosp, sinp, cose, sine, acosp, gw, dlat )
! !USES:
      use parutilitiesmodule, only : gid, commglobal, parsplit
      implicit none

! !PARAMETERS:
      type(T_TRACER_GRID), intent(in)   :: grid        ! Tracer grid

      integer, optional,  intent(out)   :: ng_d        ! Halo size

      real(r8), optional, intent(out)   :: dlon        ! Delta latitude (rads)
      real(r8), optional, intent(out)   :: dphi        ! Delta latitude (rads)

      real(r8), optional, intent(out)   :: acap        ! Area of pole cap
      real(r8), optional, intent(out)   :: rcap        ! Reciprocal of cap area

      integer, optional,  intent(out)   :: commglobal  !  YZ decomp - Comm in Y
      integer, optional,  intent(out)   :: comm_y      !  YZ decomp - Comm in Y
      integer, optional,  intent(out)   :: comm_z      !  YZ decomp - Comm in Z

      integer, optional,  intent(out)   :: jfirstyz    !  First Y on this PE  (YZ)
      integer, optional,  intent(out)   :: jlastyz     !  Last Y on this PE (YZ)
      integer, optional,  intent(out)   :: kfirstyz    !  First Z on this PE (YZ)
      integer, optional,  intent(out)   :: klastyz     !  Last Z on this PE (YZ)

      integer, optional,  intent(out)   :: myidyz_y    !  YZ decomp - Nr in Y
      integer, optional,  intent(out)   :: myidyz_z    !  YZ decomp - Nr in Z
      integer, optional,  intent(out)   :: npryz_y     !  YZ decomp - Nr in Y
      integer, optional,  intent(out)   :: npryz_z     !  YZ decomp - Nr in Z

      real(r8), optional, pointer       :: cosp(:)     ! Cos of lat angle -- vol mean
      real(r8), optional, pointer       :: sinp(:)     ! Sin of lat angle -- vol mean
      real(r8), optional, pointer       :: cose(:)     ! Cos at finite volume edge
      real(r8), optional, pointer       :: sine(:)     ! Sin at finite volume edge
      real(r8), optional, pointer       :: acosp(:)    ! Reciprocal of cos of lat angle

      real(r8), optional, pointer       :: gw(:)       ! Gaussian weights
      real(r8), optional, pointer       :: dlat(:)     ! Delta latitudes (array)

! !DESCRIPTION:
!
!   Get tracer grid variables
!
! !REVISION HISTORY:
!
!   2007.05.21   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC

   if (present(ng_d))       ng_d     = grid%ng_d

   if (present(dlon))       dlon     = grid%dlon
   if (present(dphi))       dphi     = grid%dphi

   if (present(acap))       acap     = grid%acap
   if (present(rcap))       rcap     = grid%rcap

   if (present(commglobal)) commglobal = grid%commglobal
   if (present(comm_y))     comm_y   = grid%comm_y
   if (present(comm_z))     comm_z   = grid%comm_z

   if (present(jfirstyz))   jfirstyz = grid%jfirst
   if (present(jlastyz))    jlastyz  = grid%jlast
   if (present(kfirstyz))   kfirstyz = grid%kfirst
   if (present(klastyz))    klastyz  = grid%klast

   if (present(myidyz_y))   myidyz_y = grid%myid_y
   if (present(myidyz_z))   myidyz_z = grid%myid_z
   if (present(npryz_y))    npryz_y  = grid%npryz_y
   if (present(npryz_z))    npryz_z  = grid%npryz_z

   if (present(cosp))       cosp    => grid%cosp
   if (present(sinp))       sinp    => grid%sinp
   if (present(cose))       cose    => grid%cose
   if (present(sine))       sine    => grid%sine
   if (present(acosp))      acosp   => grid%acosp

   if (present(gw))         gw      => grid%gw
   if (present(dlat))       dlat    => grid%dlat

!EOC
 end subroutine Tracer_GridGet
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !ROUTINE: Tracer_GridFinal --- finalize the tracer grid
!
! !INTERFACE:
 subroutine Tracer_GridFinal(grid)
   use parutilitiesmodule, only : parfree

! !INPUT/OUTPUT PARAMETERS:  
      type(T_TRACER_GRID), intent(inout)  :: grid     ! Resulting grid

! !DESCRIPTION:
!
!   Finalize Lin-Rood specific variables
!
! !REVISION HISTORY:
!
!   2006.06.27   Sawyer     Creation
!   2007.05.22   Sawyer     Moved to this module from FVAdvMod
!
!EOP
!-----------------------------------------------------------------------
!BOC
!

 deallocate(grid%dlat)
 deallocate(grid%gw)

 deallocate(grid%cosp)
 deallocate(grid%sinp)
 deallocate(grid%cose)
 deallocate(grid%sine)
 deallocate(grid%acosp)
 call parfree( grid%comm_y )
 call parfree( grid%comm_z )
 return

!EOC
end subroutine Tracer_GridFinal
!-----------------------------------------------------------------------

end module TracerMod
