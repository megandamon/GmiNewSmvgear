!BOP
!
! !MODULE: mod_comm --- SPMD parallel decompostion/communication module
      module mod_comm
!
! !DESCRIPTION:
!
!  \paragraph{Overview}
!
!    This module contains SPMD parallelism decomposition and
!    communication routines.  This library was originally written by
!    W. Putman and S.-J. Lin for simple gridded communications in the
!    Finite-Volume General Circulation Model (FVGCM).  Most of the
!    member functions are specific to the type of gridded data, 
!    ghost communication and decompositions used in FVGCM (which 
!    are, however, very common in atmospheric models).
!
!    The module was extended for irregular communication
!    by W. Sawyer and A. Mirin. It is now
!    a more general tool and has been incorporated into the Parallel
!    Library for Grid Manipulations (PILGRIM) which is used in the
!    Community Atmospheric Model (CAM) and the Physical-space
!    Statistical Analysis System (PSAS).
!
!    **********************************************************
!    The storage associated with the irregular communications
!    is based on CAM requirements. It runs the risk of functioning
!    improperly when used in another code.
!    **********************************************************
!
!    Irregular communication is based on the {\tt blockdescriptor}
!    derived type, which defines a set of parcels which are to be
!    send to (or received from) another PE.  The irregular 
!    communication routines operate on arrays of block descriptors
!    whose length is equal to number of PEs involved in the
!    communication.  This means the irregular communication primitives
!    are merely non-blocking (potentially) all-to-all primitives.
! 
!    Primitives using MPI1 and MPI2 have been implemented, and SHMEM is
!    implemented for the Cray X1E. The module can also make use of OpenMP,
!    and certain communications (MPI2) can be multitasked.  
!
!
!  \paragraph{Use of Global Arrays}
!
!    The module uses the concept of global arrays (coined from former
!    usage of shared memory arenas in the "multi-level parallelism"
!    (MLP) paradigm).  Global arrays are merely buffers into which
!    data are packed for the transfer to other PEs and are not
!    necessarily of global extent. All such arrays are
!    1-dimensional; they are accessed as needed with offset vars.
!
!  \paragraph{Use of MPI-2 Windows}
!
!       All implementations use real*8, real*4, and integer*4 windows
!       which are used with global arrays as follows:
!
!       \begin{itemize}
!         \item   r8\_win -> ga\_r8 - for use with real*8 types
!         \item   r4\_win -> ga\_r4 - for use with real*4 types
!         \item   i4\_win -> ga\_i4 - for use with integer*4 types
!       \end{itemize}
!
!       note: MPI routines need 2 buffers per GA, ga\_<type>\_s & ga\_<type>\_r
!             ga\_<type>\_r is used for the MPI2 windows
!
!  \paragraph{Use of Shmem}
!
!       Shmem communications on the Cray X1E are supported. The collecive
!       calls assume sufficiently large global buffers.
!
!  \paragraph{Compilation}
!
!    This module contains numerous optimizations for various platforms
!    and underlying communication primitives.   To take advantage of
!    these, the various CPP tokens can be defined.
!
!    \begin{itemize}
!      \item {\tt STAND_ALONE}:  Use as stand-alone library (if
!                                defined) or as part of CAM (if 
!                                undefined)
!      \item {\tt MODCM_TIMING}: Turn on CAM timing routines (only
!                                available if compiled in CAM framework)
!      \item {\tt USE\_MPI2}:    Use MPI-2 (one-sided communication)
!                                for underlying communication.
!      \item {\tt USE\_SHMEM}:   Use SHMEM (one-sided communication)
!                                for underlying communication.
!      \item {\tt MODCM\_STATIC}:Static allocation of buffers in
!                                module (if defined), dynamic 
!                                allocation (when undefined);
!                                for MPI2 uses mpi_alloc_mem;
!                                for SHMEM uses shpalloc;
!      \item {\tt MT\_OFF}:      Explicitly turn off multithreading
!      \item {\tt PIN_CPUS}:     (SGI specific) map tasks to hardware
!                                cpus for performance (NOT RECENTLY TESTED)
!      \item {\tt LINUX}:        Compilation specific to Linux (NOT RECENTLY TESTED)
!      \item {\tt _OPENMP}:      Implicit token (controlled by
!                                compiler) to enable OpenMP
!    \end{itemize}
!
!    
!  \paragraph{Usage}
!
!    NOTE - must call PILGRIM routine parinit to initialize before
!    making any other calls.
!
!    The public members of this module are:
!
!      \begin{itemize}
!         \item {\tt mp\_init}:          Initialize module
!         \item {\tt mp\_exit}:          Exit module
!         \item {\tt mp\_send4d\_ns}:    Ghost 4D array on north/south
!         \item {\tt mp\_recv4d\_ns}:    Complete 4D N/S ghost operation
!         \item {\tt mp\_send2\_ns}:     Ghost 2 3D arrays on north/south
!         \item {\tt mp\_recv2\_ns}:     Complete 2x3D N/S ghost operation
!         \item {\tt mp\_send3d}:        Send 3D general ghost region
!         \item {\tt mp\_recv3d}:        Complete 3D general ghost operation
!         \item {\tt mp\_send3d\_2}:     Send 2x3D general ghost regions
!         \item {\tt mp\_recv3d\_2}:     Complete 2x3D general ghost operation
!         \item {\tt get\_partneroffset}:Offset for remote write
!         \item {\tt mp\_sendirr}:       Initiate all-to-all send of parcels
!         \item {\tt mp\_recvirr}:       Complete all-to-all chunk commun.
!       \end{itemize}
!
!     There are variants of some of these routines for r4 and i4 data types.
!     There are other public routines, but these are only used internally
!     in PILGRIM, and they should not be called by user applications.
!
! !REVISION HISTORY:
!    2001.09.01   Lin
!    2002.04.16   Putman  Modified for Global Array code
!    2002.04.16   Putman  Added ProTeX documentation
!    2002.05.28   Putman  Added use of precision module
!    2003.06.24   Sawyer  Minor additions for use with mod_irreg
!    2004.01.08   Sawyer  Removed older functionality, no longer needed
!    2004.02.10   Mirin   Major restructuring and simplification. Documentation
!    2004.03.06   Sawyer  Additional documentation; cosmetics
!    2005.03.20   Sawyer  Added extensive support for real*4
!    2005.10.12   Worley  Improved vectorization of buffer copies and general clean-up
!    2006.05.15   Mirin   Make dynamic allocation the default; general clean-up.
! !USES:
#if defined( STAND_ALONE )
# define iulog 6
#else
      use cam_logfile, only: iulog
#endif
!
#include "pilgrim.h"
!
! Mod_comm has option for stand-alone use as well as within CAM
!

#if defined ( SPMD )

#if defined( STAND_ALONE )
# define r8 selected_real_kind(12)
# define r4 selected_real_kind( 6)
# define i8 selected_int_kind(13)
# define i4 selected_int_kind( 6)
# define PLON        144
# define PLAT         91
# define PLEV         26
# define PCNST         1
#else
      use shr_kind_mod, only : r8 => shr_kind_r8, r4 => shr_kind_r4,  &
                               i8 => shr_kind_i8, i4 => shr_kind_i4
#endif
#if defined( MODCM_TIMING )
      use perf_mod
#endif

#if defined(UNICOSMP)
# if defined (USE_MPI2) || defined (USE_SHMEM)
# define MT_OFF
# endif
#endif

      implicit none

!
! Shmem option implemented on top of MPI for Cray X1E
!
#include "mpif.h"
#if defined ( USE_SHMEM )
# include "mpp/shmem.fh"
#endif

! !PUBLIC MEMBER FUNCTIONS:
      public mp_init, mp_exit,                                             &
             mp_send4d_ns, mp_recv4d_ns, mp_send4d_ns_r4, mp_recv4d_ns_r4, &
             mp_send2_ns, mp_recv2_ns, mp_send3d_2, mp_recv3d_2,           &
             mp_send3d, mp_recv3d, mp_sendirr, mp_recvirr,                 &
             mp_sendirr_r4, mp_recvirr_r4, mp_sendirr_i4, mp_recvirr_i4,   &
             mp_barrier, get_partneroffset, mp_r8, mp_r4, mp_i4
      public modcam_method, modcam_geopk, modcam_gatscat, modcam_npryz
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
! Used for mpi 2-sided communications
      public sqest, rqest, nsend, nread
#endif

!------------------------------------------------------------------------------
!  type declaration for describing an arbitrary number of contiguous parcels
!  this is for irregular communications
!------------------------------------------------------------------------------
      type blockdescriptor
         integer              :: method             ! transpose method
         integer              :: type               ! Ptr to MPI derived type
         integer, pointer     :: displacements(:)   ! Offsets in local segment
         integer, pointer     :: blocksizes(:)      ! Block sizes to transfer
         integer              :: partneroffset      ! Aggregated partner offset
         integer              :: partnertype        ! Ptr to partner's MPI derived type
         integer              :: Nparcels           ! size( displacements )
         integer              :: Tot_Size           ! sum ( blocksizes )
      end type blockdescriptor

! transpose methods (method)
!      0 for contiguous temporary buffer
!      (MPI-1) 1 for direct communication (derived types)
!      (MPI-2) 1 or 2 for direct put into contiguous temporary target window
!                 1 for OpenMP over contiguous segments destined for target task
!                 2 for OpenMP over target tasks, using derived types for each target

!     MAX_TRf is defined in pilgrim.h and is the maximum number of simultaneous
!      messages (indexed by BegTrf and EndTrf).
      INTEGER, ALLOCATABLE, SAVE :: InHandle(:, :)
      INTEGER, ALLOCATABLE, SAVE :: OutHandle(:, :)
      INTEGER, SAVE :: BegTrf = 0  ! Ongoing overlapped begintransfer #
      INTEGER, SAVE :: EndTrf = 0  ! Ongoing overlapped endtransfer #

! !PUBLIC DATA MEMBERS:
      integer, SAVE:: gid                         ! PE id
      integer(i4), SAVE:: masterpro = 0           ! Master process id 
      integer(i4), SAVE:: numpro                  ! Permanent No. of PEs
      integer(i4), SAVE:: numcomm                 ! Local No. of PEs
      integer(i4), SAVE:: numcpu                  ! No. of threads
      integer, SAVE:: commglobal                  ! Global Communicator
      integer, SAVE:: Max_Nparcels = 0            ! Maximum number of parcels in
                                                  !  single blockdescriptor

!------------------------------------------------------------------------------
!  Local parameters for use with MPI-2 and MPI-1
!------------------------------------------------------------------------------
      integer, parameter:: nbuf = 2               ! Max No. of sends per call
      integer, parameter:: nghost = 3             ! No. of ghost indices
      integer, parameter:: max_nq = 1             ! No. of tracers simultaneously
                                                  !  border communicated; can be
                                                  !  overridden with dynamic storage
      integer, parameter:: max_trac = PCNST       ! No. of tracers
      integer, parameter:: max_call = 2           ! Max No. of back-to-back...
                                                  ! ...mp_send calls
      integer, parameter:: idimsize = PLON*nghost*(PLEV+1)*max_nq
                                                  ! Size of MPI buffer region
                                                  ! in mp_send/mp_recv calls, used
                                                  ! to determine offset in GA
      integer, parameter:: platg = PLAT + 2*nghost

#if defined(USE_SHMEM)
      integer, parameter :: mp_r4 = r4
      integer, parameter :: mp_r8 = r8
      integer, parameter :: mp_i4 = i4
      integer, SAVE::  reduce_sync(SHMEM_REDUCE_SYNC_SIZE)
      integer, SAVE::  collect_sync(SHMEM_COLLECT_SYNC_SIZE)
      integer, SAVE::  bcast_sync(SHMEM_BCAST_SYNC_SIZE)
      integer, SAVE::  barrier_sync(SHMEM_BARRIER_SYNC_SIZE)
#else
      integer, parameter :: mp_r4 = MPI_REAL
      integer, parameter :: mp_r8 = MPI_DOUBLE_PRECISION
      integer, parameter :: mp_i4 = MPI_INTEGER
#endif

      integer ierror

!------------------------------------------------------------------------------
!  Local variables for use with MPI-2, MPI-1 and SHMEM
!------------------------------------------------------------------------------
      integer, SAVE:: ncall_r, ncall_s

      integer, SAVE:: sizet1, sizer8, sizer4, sizei4

! CAM-specific variables
      integer, SAVE:: tracmax, tracbmax, dpvarmax, totvar
      integer, SAVE:: phys_transpose_mod
      integer, SAVE:: idimsizz
      integer, SAVE:: modcam_method, modcam_geopk, modcam_gatscat
      integer, SAVE:: modcam_npryz(4)
      integer, parameter :: phys_transpose_modmin = 11
      integer, parameter :: phys_transpose_vars = 7
      data phys_transpose_mod / -1 /
      data modcam_method / -1 /
      data modcam_geopk / -1 /
      data modcam_gatscat / -1 /
      data modcam_npryz / -1, -1, -1, -1 /
!
! tracmax is the maximum number of tracers simultaneously transposed within dynamics (set to 1)
!    (except in dynamics-physics transposes)
! tracbmax is the maximum number of tracers simultaneously border communicated
! dpvarmax is the number of variables communicated in dynamics-physics transposes 
! totvar is the maximum number of variables simultaneously transposed
! phys_transpose_mod is the communication method for dynamics/physics transposes; admissable values
!     are >= phys_transpose_modmin; it is communicated from CAM when such transposes
!     are requested.
! phys_transpose_vars is the number of non-tracer variables transposed between dynamics and
!     physics instantiations in CAM.
! modcam_method, modcam_geopk and modcam_gatscat correspond to mod_method, mod_geopk and
!     mod_gatscat in CAM.
! modcam_npryz corresponds to npr_yz in CAM.

!------------------------------------------------------------------------------
!  Variables to control global array locations and window synchronization
!------------------------------------------------------------------------------
      integer win_count                 ! Counts No. of windows in use
      integer lastwin                   ! ID of last synch'd window
      integer pkgs_per_pro              ! No. of MPI-2 packages per PE with OpenMP
      integer igosouth, igonorth        ! Index of latitudinal send direction
      integer ifromsouth, ifromnorth    ! Index of latitudinal recv direction

!------------------------------------------------------------------------------
!  Local type declaration for mp_windows
!------------------------------------------------------------------------------
      type window
         integer :: id            ! Window id
         integer :: size          ! Size of global window (point based)
         integer :: ncall_s       ! Count send calls on window
         integer :: ncall_r       ! Count recv calls on window
#if defined(USE_MPI2)

!  The lines immediately below assume that MPI-2 is not implemented
!  on the specified architectures. This faulty assumption is historical
!  and must be modified for machines that do support MPI-2.
!  We are supporting MPI-2 at the compile stage but not at the
!  execution stage.
# if defined(LINUX)
#  define MPI_OFFSET_KIND 8
#  define MPI_ADDRESS_KIND 8
# endif

         integer(kind=MPI_ADDRESS_KIND) :: offset_s
         integer(kind=MPI_ADDRESS_KIND) :: offset_r
#else
         integer :: offset_s      ! Starting position in GA send
         integer :: offset_r      ! Starting position in GA recv
#endif
         integer :: dest          ! For use with send calls
         integer :: src           ! For use with recv calls
         integer :: size_r        ! Size of incoming message
     end type window

!------------------------------------------------------------------------------
! Beginning Global Array variable declaration:
!------------------------------------------------------------------------------

      type (window) :: r8_win
      type (window) :: r4_win
      type (window) :: i4_win

!
!  Note: t1 (below) and r8 windows can be safely joined provided there is no
!    overlap between transpose and border communications
!

! Upper bound on ratio of local to average storage over subdomains.
!
! This takes into account different sized domain decompositions. Also,
! ghost points are not figured into the computation of the average. In an
! extreme situation, the average number of latitudes in a subdomain could
! equal 3, but with one subdomain having 4. Considering 3 ghost points on
! each side, that would give us 10 latitudes versus the expected 3, requiring
! alloc_slack_factor of 10/3. Throw in some unevenness in the vertical, and
! it could be larger. Ghost points can also occur in the vertical due to
! edge quantities, but that does not occur simultaneously with ghost points
! in latitude.

      real*8, parameter :: alloc_slack_factor = 8.d0

!
!  SHMEM variable declarations
!

#if defined(USE_SHMEM)
      integer ga_ptr
!
! Join t1 and r8 windows
!
# define ga_t1_r ga_r8_r
# define ga_t1_s ga_r8_s
# define t1_win  r8_win

# if defined(MODCM_STATIC)

      real(r8), SAVE, TARGET:: ga_r8a_r(MAX( PLON*platg*(PLEV+1)*max_nq, idimsize*nbuf*max_call ))
      real(r8), SAVE, TARGET:: ga_r8a_s(MAX( PLON*platg*(PLEV+1)*max_nq, idimsize*nbuf*max_call ))

      real(r8), SAVE, TARGET:: ga_r8b_r(MAX( PLON*platg*(PLEV+1)*max_nq, idimsize*nbuf*max_call ))
      real(r8), SAVE, TARGET:: ga_r8b_s(MAX( PLON*platg*(PLEV+1)*max_nq, idimsize*nbuf*max_call ))

      real(r8), DIMENSION(:), SAVE, POINTER:: ga_r8_r
      real(r8), DIMENSION(:), SAVE, POINTER:: ga_r8_s

      real(r4), SAVE, TARGET:: ga_r4a_r(MAX( PLON*platg*(PLEV+1)*max_nq, idimsize*nbuf*max_call ))
      real(r4), SAVE, TARGET:: ga_r4a_s(MAX( PLON*platg*(PLEV+1)*max_nq, idimsize*nbuf*max_call ))

      real(r4), SAVE, TARGET:: ga_r4b_r(MAX( PLON*platg*(PLEV+1)*max_nq, idimsize*nbuf*max_call ))
      real(r4), SAVE, TARGET:: ga_r4b_s(MAX( PLON*platg*(PLEV+1)*max_nq, idimsize*nbuf*max_call ))

      real(r4), DIMENSION(:), SAVE, POINTER:: ga_r4_r
      real(r4), DIMENSION(:), SAVE, POINTER:: ga_r4_s

      integer(i4), SAVE, TARGET:: ga_i4a_r(PLON*PLAT*PLEV)
      integer(i4), SAVE, TARGET:: ga_i4a_s(PLON*PLAT*PLEV)

      integer(i4), SAVE, TARGET:: ga_i4b_r(PLON*PLAT*PLEV)
      integer(i4), SAVE, TARGET:: ga_i4b_s(PLON*PLAT*PLEV)

      integer(i4), DIMENSION(:), SAVE, POINTER:: ga_i4_r
      integer(i4), DIMENSION(:), SAVE, POINTER:: ga_i4_s
# else
      real(r8) ga_r8a_r
      real(r8) ga_r8a_s
      real(r8) ga_r8b_r
      real(r8) ga_r8b_s
      real(r8) ga_r8_r
      real(r8) ga_r8_s
      real(r4) ga_r4a_r
      real(r4) ga_r4a_s
      real(r4) ga_r4b_r
      real(r4) ga_r4b_s
      real(r4) ga_r4_r
      real(r4) ga_r4_s
      integer(i4) ga_i4a_r
      integer(i4) ga_i4a_s
      integer(i4) ga_i4b_r
      integer(i4) ga_i4b_s
      integer(i4) ga_i4_r
      integer(i4) ga_i4_s
!
! Cray pointers required for shpalloc
!
      pointer (pa_r8a_r, ga_r8a_r(1))
      pointer (pa_r8a_s, ga_r8a_s(1))
      pointer (pa_r8b_r, ga_r8b_r(1))
      pointer (pa_r8b_s, ga_r8b_s(1))
      pointer (pa_r8_r,  ga_r8_r(1))
      pointer (pa_r8_s,  ga_r8_s(1))
      pointer (pa_r4a_r, ga_r4a_r(1))
      pointer (pa_r4a_s, ga_r4a_s(1))
      pointer (pa_r4b_r, ga_r4b_r(1))
      pointer (pa_r4b_s, ga_r4b_s(1))
      pointer (pa_r4_r,  ga_r4_r(1))
      pointer (pa_r4_s,  ga_r4_s(1))
      pointer (pa_i4a_r, ga_i4a_r(1))
      pointer (pa_i4a_s, ga_i4a_s(1))
      pointer (pa_i4b_r, ga_i4b_r(1))
      pointer (pa_i4b_s, ga_i4b_s(1))
      pointer (pa_i4_r,  ga_i4_r(1))
      pointer (pa_i4_s,  ga_i4_s(1))

      save pa_r8a_r
      save pa_r8a_s
      save pa_r8b_r
      save pa_r8b_s
      save pa_r8_r
      save pa_r8_s
      save pa_r4a_r
      save pa_r4a_s
      save pa_r4b_r
      save pa_r4b_s
      save pa_r4_r
      save pa_r4_s
      save pa_i4a_r
      save pa_i4a_s
      save pa_i4b_r
      save pa_i4b_s
      save pa_i4_r
      save pa_i4_s
# endif
#endif

!
!  MPI-1 and MPI-2 window variable declarations
!

#if !defined(USE_SHMEM)
      type (window) :: t1_win
# if defined(MODCM_STATIC)
      real(r8),    SAVE:: ga_t1_r(idimsize*nbuf*max_call)
      real(r8),    SAVE:: ga_t1_s(idimsize*nbuf*max_call)
      real(r8),    SAVE:: ga_r8_r(PLON*platg*(PLEV+1)*max_nq)
      real(r8),    SAVE:: ga_r8_s(PLON*platg*(PLEV+1)*max_nq)
      real(r4),    SAVE:: ga_r4_r(PLON*platg*(PLEV+1)*max_nq)
      real(r4),    SAVE:: ga_r4_s(PLON*platg*(PLEV+1)*max_nq)
      integer(i4), SAVE:: ga_i4_r(PLON*PLAT*PLEV)
      integer(i4), SAVE:: ga_i4_s(PLON*PLAT*PLEV)
!
!  t1 and r8 windows could be joined through equivalence statements,
!    but that might degrade performance
!
# else
#  if defined(USE_MPI2)
      real(r8) ga_t1_r
      real(r8) ga_t1_s
      real(r8) ga_r8_r
      real(r8) ga_r8_s
      real(r4) ga_r4_r
      real(r4) ga_r4_s
      integer(i4) ga_i4_r
      integer(i4) ga_i4_s
!
! Cray pointers required for mpi_alloc_mem
!
      pointer (pa_t1_r, ga_t1_r(1))
      pointer (pa_t1_s, ga_t1_s(1))
      pointer (pa_r8_r, ga_r8_r(1))
      pointer (pa_r8_s, ga_r8_s(1))
      pointer (pa_r4_r, ga_r4_r(1))
      pointer (pa_r4_s, ga_r4_s(1))
      pointer (pa_i4_r, ga_i4_r(1))
      pointer (pa_i4_s, ga_i4_s(1))
      save pa_t1_r
      save pa_t1_s
      save pa_r8_r
      save pa_r8_s
      save pa_r4_r
      save pa_r4_s
      save pa_i4_r
      save pa_i4_s
#  else
      real(r8), allocatable,    SAVE:: ga_t1_r(:)
      real(r8), allocatable,    SAVE:: ga_t1_s(:)
      real(r8), allocatable,    SAVE:: ga_r8_r(:)
      real(r8), allocatable,    SAVE:: ga_r8_s(:)
      real(r4), allocatable,    SAVE:: ga_r4_r(:)
      real(r4), allocatable,    SAVE:: ga_r4_s(:)
      integer(i4), allocatable, SAVE:: ga_i4_r(:)
      integer(i4), allocatable, SAVE:: ga_i4_s(:)
#  endif
# endif
#endif

!
!  MPI-2 auxiliary variable declarations
!

#if defined(USE_MPI2)

!------------------------------------------------------------------------------
!  The lines immediately below assume that MPI-2 is not implemented
!  on the specified architectures. This faulty assumption is historical
!  and must be modified for machines that do support MPI-2.
!  We are supporting MPI-2 at the compile stage but not at the
!  execution stage.
!
# if defined(LINUX)
# define MPI_ADDRESS_KIND 8
      integer, parameter:: MPI_MODE_NOCHECK    = 0 
      integer, parameter:: MPI_MODE_NOSTORE    = 0    
      integer, parameter:: MPI_MODE_NOPUT      = 0 
      integer, parameter:: MPI_MODE_NOPRECEDE  = 0 
      integer, parameter:: MPI_MODE_NOSUCCEED  = 0 
# endif
!------------------------------------------------------------------------------

      integer(kind=MPI_ADDRESS_KIND) bsize
      integer, SAVE:: Status(MPI_STATUS_SIZE)
#endif

!
!  MPI-1 auxiliary variable declarations
!

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
      integer, SAVE:: nsend                   ! Number of messages out-going
      integer, SAVE:: nrecv                   ! Number of messages in-coming
      integer, SAVE:: nread                   ! Number of messages read
      integer, allocatable, SAVE:: sqest(:)   ! Request handler for sends
      integer, allocatable, SAVE:: rqest(:)   ! Request handler for recvs
      integer, SAVE:: bsize             
      integer, SAVE:: Status(MPI_STATUS_SIZE)
      integer, allocatable, SAVE:: Stats(:)
#endif

!
!EOP
!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_init --- Initialize SPMD parallel communication
!
! !INTERFACE:
      subroutine mp_init( comm, npryzxy, mod_method, mod_geopk, mod_gatscat )
!
! !INPUT PARAMETERS:
      integer, optional :: comm                        ! communicator
      integer, optional, intent(in) :: npryzxy(4)      ! 2D decomposition
      integer, optional, intent(in) :: mod_method      ! CAM optimization
      integer, optional, intent(in) :: mod_geopk       ! CAM optimization
      integer, optional, intent(in) :: mod_gatscat     ! CAM optimization
! !DESCRIPTION:
!
!     Initialize SPMD parallel communication.  It is recommended that
!     COMM (main communicator) and NPRYZXY (2D decomposition) be set.
!
!     Set the mod* variables only if you are acquainted with their 
!     meaning (default is 0).
!
! !REVISION HISTORY: 
!    2001.09.01   Lin
!    2002.02.15   Putman        Modified for Global Array code
!    2002.04.09   Putman        Added ProTeX documentation
!    2002.08.06   Sawyer        Added optional communicator input argument
!    2006.06.15   Sawyer        Added CAM-dependent optional arguments
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer mysize
#if defined(USE_MPI2)
      integer (KIND=MPI_ADDRESS_KIND) sizet18
      integer (KIND=MPI_ADDRESS_KIND) sizer88
      integer (KIND=MPI_ADDRESS_KIND) sizer48
      integer (KIND=MPI_ADDRESS_KIND) sizei48
#endif
#if defined(USE_SHMEM)
      integer sizetrm
#endif
      integer using_window, vertical_lines
      integer local_dynamic_storage
      real*8 geopkrat, one

! Initialize MPI; allow for general communicator
      if ( present(comm) ) then
        call mpi_start( comm )
      else
        call mpi_start( MPI_COMM_WORLD )
      endif
! Initialize OpenMP
      call omp_start

!
! Adopt 2D decomposition if provided.
!
      modcam_npryz = (/ 1,1,1,1 /)    ! Default value (sequential)
      if ( present( npryzxy ) ) modcam_npryz(1:4) = npryzxy(1:4)
      if (gid .eq. 0) then
        write(iulog,*) 'MOD_COMM - modcam_npryz = ', &
               modcam_npryz(1), modcam_npryz(2),     &
               modcam_npryz(3), modcam_npryz(4)
      endif

!
! Set CAM optimization variables
!
! modcam_method refers to irregular communications for transposes
! modcam_geopk refers to irregular communications for the geopotential
! modcam_gatscat refers to irregular communications for gather/scatters
! For any of these, a value of 0 means source data will be gathered into a contiguous
!  buffer (window), communicated to a contiguous buffer (window) in the target, and
!  then scattered to its final destination; a value of 1 means MPI derived types will
!  be used (hence not requiring window storage).
      modcam_method  = 0   ! Default value
      modcam_geopk   = 0   ! Default value
      modcam_gatscat = 0   ! Default value
      if ( present( mod_method ) )  modcam_method  = mod_method
      if ( present( mod_geopk ) )   modcam_geopk   = mod_geopk
      if ( present( mod_gatscat ) ) modcam_gatscat = mod_gatscat

      if (gid .eq. 0) then
        write(iulog,*) 'MOD_COMM - modcam_method modcam_geopk modcam_gatscat = ',    &
        modcam_method, modcam_geopk, modcam_gatscat
      endif

!
! End CAM optimizations
!

      win_count = 0
!
!*************************************************************************
! local_dynamic_storage is set to 1 when window storage is based on locally dimensioned
!  arrays, 0 otherwise; this occurs when modcam_gatscat equals 1, as it is only the
!  gather/scatters that require global storage.
!*************************************************************************
!
      local_dynamic_storage = 0
# if !defined(USE_SHMEM) && !defined(USE_MPI2)
      if (modcam_gatscat .eq. 1) local_dynamic_storage = 1
# endif

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
! Used for MPI 2-sided communications
      allocate( sqest(MAX(nbuf,numpro)*max_call) )
      allocate( rqest(MAX(nbuf,numpro)*max_call) )
      allocate( Stats(MAX(nbuf,numpro)*max_call*MPI_STATUS_SIZE) )
      allocate( InHandle(MAX(nbuf,numpro),MAX_TRF) )
      allocate( OutHandle(MAX(nbuf,numpro),MAX_TRF) )
#endif

      idimsizz = idimsize
#if defined(MODCM_STATIC)
      if (gid .eq. 0) write(iulog,*) 'Using static allocation in mod_comm'
#else
      if (gid .eq. 0) write(iulog,*) 'Using dynamic allocation in mod_comm'
      if (local_dynamic_storage .eq. 1) then
         if (gid .eq. 0) write(iulog,*) 'Using local dynamic storage for mod_comm window'
      else
         if (gid .eq. 0) write(iulog,*) 'Using global dynamic storage for mod_comm window'
      endif
!
! Dynamically allocate target global arrays
!
!*************************************************************************
! Compute extent to which required window storage for geopotential computation
!   exceeds that of transpose - relevant only for local dynamic storage,
!   since with global storage there will be enough space anyway; also,
!   this applies only when using window; further, this applies only when
!   the CAM variable geopktrans equals 1, though we do not test for that here.
! The geopotential calculation sends a latitude line to every other process
!   either vertically above or below the given process; there can be
!   at most modcam_npryz(2)-1 such target processes; compared to transposes
!   (which send all vertical lines), the amount of data sent is expressed
!   as the ratio geopkrat; our concern is making the window (whose size
!   is computed based on transposes) large enough, so we must multiply its
!   size by geopkrat; we never shrink the window, so geopkrat >= 1.
!*************************************************************************
      using_window = 1   !  This is a local variable
# if !defined(USE_SHMEM) && !defined(USE_MPI2)
      if (modcam_geopk .eq. 1) using_window = 0
# endif
      one = real(1,r8)
      geopkrat = one
      if (using_window .eq. 1 .and. local_dynamic_storage .eq. 1) then
         vertical_lines = ceiling(real(PLEV,r8)/real(modcam_npryz(2),r8))
         geopkrat = real(modcam_npryz(2)-1,r8)/real(vertical_lines,r8)
         geopkrat = min(geopkrat,one)
      endif
      if (gid .eq. 0) write(iulog,*) 'Mod_comm - geopkrat = ', geopkrat

!*************************************************************************
! beginning of CAM totvar computation
!*************************************************************************

! CAM contains two kinds of transposes. The most commonly referred to transposes
!  refer to those which connect the xy and yz decompositions. Depending on
!  the physics decomposition, CAM might additionally compute transposes between
!  the dynamics and physics; this depends on the variable phys_loadbalance.
!  Furthermore, these transposes might or might not be computed using mod_comm.
!  The former transposes are generally performed one variable at a time; the
!  latter transposes combine all variables to be transposed, including the
!  full complement of tracers. The maximum number of variables to be 
!  simultaneously subject to irregular communications is dependent on
!  whether or not mod_comm is used to compute dynamics-physics transposes
!  and could depend on the number of tracers.

! Compute maximum number of variables to be simultaneously subject
!  to irregular communications (e.g., transposed variables based on CAM)
!  and store in the variable 'totvar'.

! Tracmax is the number of tracers simultaneously transposed within dynamics;
! Tracbmax is the number of tracers simultaneously border comunicated within trac2d;
!  both of these are currently hardwired to 1.
      tracmax = 1
      tracbmax = 1
      totvar = tracmax

! Now consider dynamics-physics transposes in CAM dp_coupling (dpvarmax)
!  If phys_transpose_mod is still -1, that means it has not been updated
!  by CAM and hence mod_comm will not be used for dynamics-physics transposes.
! (NOTE: phys_transpose_mod is computed in phys_grid_setopts in phys_grid.F90.)

! Also note that the logic involving phys_transpose_mod and phys_transpose_modmin
!  must remain consistent with the coding in phys_grid.F90. Additionally,
!  phys_transpose_vars must remain consistent with the coding in dp_coupling.F90.
!  (See above declaration and initialization for CAM-specific variables.)

! (begin dpvarmax calculation)

      if (phys_transpose_mod .eq. -1) then
         if (gid .eq. 0) write(iulog,*)       &
           '(MOD_COMM) - mod_comm not being used for dynamcis-physics transposes'
         dpvarmax = 0
!
! If phys_transpose_mod is >= phys_transpose_modmin, that is a signal that mod_comm is to be used
!  for dynamics/physics transposes in CAM. In that case, one must allocate enough window
!  storage for those transposes. Presently, the number of such simultaneously transposed
!  variables equals phys_transpose_vars plus the number of constituents.
!
      elseif (phys_transpose_mod .ge. phys_transpose_modmin) then
         dpvarmax = phys_transpose_vars + max_trac
      else
         dpvarmax = 0
      endif

! (end dpvarmax calculation)

! totvar is the maximum of (1) the number of tracers to be simultaneously transposed
!  within the dynamics, and (2) the number of variables to be transposed between
!  dynamics and physics instantiations in CAM

      totvar = max(totvar, dpvarmax)

!*************************************************************************
! end of CAM totvar computation
!*************************************************************************

! totvar is used to determine window size, so when window is not needed, it is
!  reset to 1

      using_window = 1   !  This is a local variable
# if !defined(USE_SHMEM) && !defined(USE_MPI2)
      if (modcam_method .eq. 1) using_window = 0
# endif
      if (using_window .eq. 0) totvar = 1

      if (gid .eq. 0) write(iulog,*) 'Mod_comm - tracmax dpvarmax totvar tracbmax = ',     &
          tracmax, dpvarmax, totvar, tracbmax 

      idimsizz = (idimsize/max_nq)*tracbmax
      sizet1 = idimsizz*nbuf*max_call
      sizer8 = PLON*platg*(PLEV+1)*totvar*geopkrat
      sizer4 = PLON*platg*(PLEV+1)*totvar*geopkrat
      sizei4 = PLON*PLAT*PLEV

! Compute local storage requirement for irregular communications by dividing
!    global requirement by the number of tasks. Allow slack factor to account
!    for nonuniformity of decomposition and ghost zones. Not valid for global
!    operations such as gathers and scatters when local windows are used.
      if (local_dynamic_storage .eq. 1) then
         sizer8 = ceiling( alloc_slack_factor*real(sizer8,r8)/real(numpro,r8) )
         sizer4 = ceiling( alloc_slack_factor*real(sizer4,r8)/real(numpro,r8) )
         sizei4 = ceiling( alloc_slack_factor*real(sizei4,r8)/real(numpro,r8) )
      endif
# if defined ( NOR4 )
      sizer4 = 1
      if (gid .eq. 0) write(iulog,*) 'Mod_comm - r4 windows disabled'
# endif

! Allocate global storage

# if defined(USE_MPI2)
      sizet18 = 8*sizet1
      sizer88 = 8*sizer8
      sizer48 = 8*sizer4
      sizei48 = 8*sizei4
      call mpi_alloc_mem(sizet18, mpi_info_null, pa_t1_r, ierror)
      if (ierror .ne. 0) write(iulog,*) 'MPI_ALLOC_MEM error'
      call mpi_alloc_mem(sizet18, mpi_info_null, pa_t1_s, ierror)
      if (ierror .ne. 0) write(iulog,*) 'MPI_ALLOC_MEM error'
      call mpi_alloc_mem(sizer88, mpi_info_null, pa_r8_r, ierror)
      if (ierror .ne. 0) write(iulog,*) 'MPI_ALLOC_MEM error'
      call mpi_alloc_mem(sizer88, mpi_info_null, pa_r8_s, ierror)
      if (ierror .ne. 0) write(iulog,*) 'MPI_ALLOC_MEM error'
      call mpi_alloc_mem(sizer48, mpi_info_null, pa_r4_r, ierror)
      if (ierror .ne. 0) write(iulog,*) 'MPI_ALLOC_MEM error'
      call mpi_alloc_mem(sizer48, mpi_info_null, pa_r4_s, ierror)
      if (ierror .ne. 0) write(iulog,*) 'MPI_ALLOC_MEM error'
      call mpi_alloc_mem(sizei48, mpi_info_null, pa_i4_r, ierror)
      if (ierror .ne. 0) write(iulog,*) 'MPI_ALLOC_MEM error'
      call mpi_alloc_mem(sizei48, mpi_info_null, pa_i4_s, ierror)
      if (ierror .ne. 0) write(iulog,*) 'MPI_ALLOC_MEM error'
# elif defined(USE_SHMEM)
      sizetrm=max(sizet1,sizer8)
! Allocation units are in terms of integer words (hence factor of 2)
      call shpalloc(pa_r8a_r, 2*sizetrm, ierror, 1)
      if (ierror .ne. 0) write(iulog,*) 'shpalloc error'
      call shpalloc(pa_r8a_s, 2*sizetrm, ierror, 1)
      if (ierror .ne. 0) write(iulog,*) 'shpalloc error'
      call shpalloc(pa_r8b_r, 2*sizetrm, ierror, 1)
      if (ierror .ne. 0) write(iulog,*) 'shpalloc error'
      call shpalloc(pa_r8b_s, 2*sizetrm, ierror, 1)
      if (ierror .ne. 0) write(iulog,*) 'shpalloc error'
      call shpalloc(pa_r4a_r, sizer4,  ierror, 1)
      if (ierror .ne. 0) write(iulog,*) 'shpalloc error'
      call shpalloc(pa_r4a_s, sizer4,  ierror, 1)
      if (ierror .ne. 0) write(iulog,*) 'shpalloc error'
      call shpalloc(pa_r4b_r, sizer4,  ierror, 1)
      if (ierror .ne. 0) write(iulog,*) 'shpalloc error'
      call shpalloc(pa_r4b_s, sizer4,  ierror, 1)
      if (ierror .ne. 0) write(iulog,*) 'shpalloc error'
      call shpalloc(pa_i4a_r, sizei4,  ierror, 1)
      if (ierror .ne. 0) write(iulog,*) 'shpalloc error'
      call shpalloc(pa_i4a_s, sizei4,  ierror, 1)
      if (ierror .ne. 0) write(iulog,*) 'shpalloc error'
      call shpalloc(pa_i4b_r, sizei4,  ierror, 1)
      if (ierror .ne. 0) write(iulog,*) 'shpalloc error'
      call shpalloc(pa_i4b_s, sizei4,  ierror, 1)
      if (ierror .ne. 0) write(iulog,*) 'shpalloc error'
      write(iulog,*) 'Shpalloc completed'
# else
      allocate( ga_t1_r(sizet1) )
      allocate( ga_t1_s(sizet1) )
      allocate( ga_r8_r(sizer8) )
      allocate( ga_r8_s(sizer8) )
      allocate( ga_r4_r(sizer4) )
      allocate( ga_r4_s(sizer4) )
      allocate( ga_i4_r(sizei4) )
      allocate( ga_i4_s(sizei4) )
# endif
#endif

#if defined(USE_SHMEM)
        ga_ptr = 1
# if defined(MODCM_STATIC)
        ga_r8_r => ga_r8a_r
        ga_r8_s => ga_r8a_s
        ga_r4_r => ga_r4a_r
        ga_r4_s => ga_r4a_s
        ga_i4_r => ga_i4a_r
        ga_i4_s => ga_i4a_s
# else
        pa_r8_r = pa_r8a_r
        pa_r8_s = pa_r8a_s
        pa_r4_r = pa_r4a_r
        pa_r4_s = pa_r4a_s
        pa_i4_r = pa_i4a_r
        pa_i4_s = pa_i4a_s
# endif
        reduce_sync = SHMEM_SYNC_VALUE
        collect_sync = SHMEM_SYNC_VALUE
        bcast_sync = SHMEM_SYNC_VALUE
        barrier_sync = SHMEM_SYNC_VALUE
#endif

! Initialize windows

#if !defined(USE_SHMEM)
# if defined(MODCM_STATIC)
        mysize = idimsize*nbuf*max_call
# else
        mysize = sizet1
# endif
        call win_init_r8(comm, t1_win, ga_t1_r, mysize)
        if (gid .eq. 0) write(iulog,*) 'Mod_comm t1_win window size = ', mysize
#endif

#if defined(MODCM_STATIC)
        mysize = PLON*platg*(PLEV+1)*max_nq
#else
        mysize = sizer8
#endif
        call win_init_r8(comm, r8_win, ga_r8_r, mysize)
        if (gid .eq. 0) write(iulog,*) 'Mod_comm r8_win window size = ', mysize

#if defined(MODCM_STATIC)
        mysize = PLON*platg*(PLEV+1)*max_nq
#else
        mysize = sizer4
#endif
        call win_init_r4(comm, r4_win, ga_r4_r, mysize)
        if (gid .eq. 0) write(iulog,*) 'Mod_comm r4_win window size = ', mysize

#if defined(MODCM_STATIC)
        mysize = PLON*PLAT*PLEV
#else
        mysize = sizei4
#endif
        call win_init_i4(comm, i4_win, ga_i4_r, mysize)
        if (gid .eq. 0) write(iulog,*) 'Mod_comm i4_win window size = ', mysize

        igosouth   = 0
        igonorth   = 1
        ifromsouth = 1
        ifromnorth = 0

        ncall_s = 0
        ncall_r = 0
        lastwin = r8_win%id

#if defined(USE_SHMEM)
        if (gid .eq. 0) write(iulog,*) 'Using Shmem communications in mod_comm'
#elif defined(USE_MPI2)
        if (gid .eq. 0) write(iulog,*) 'Using Mpi2 communications in mod_comm'
#else
        if (gid .eq. 0) write(iulog,*) 'Using Mpi1 communications in mod_comm'
#endif
!EOC
      end subroutine mp_init
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_exit --- End SPMD parallel communication
!
! !INTERFACE:
      subroutine mp_exit( comm )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
! !DESCRIPTION:
!
!     End SPMD parallel communication
!
! !REVISION HISTORY: 
!    2001.09.01   Lin
!    2002.02.15   Putman        Modified for Global Array code
!    2002.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if !defined(USE_SHMEM)
# if defined(USE_MPI2)
        call MPI_WIN_FREE( t1_win%id, ierror )
        call MPI_WIN_FREE( r8_win%id, ierror )
        call MPI_WIN_FREE( i4_win%id, ierror )
# endif
        call MPI_FINALIZE (ierror)
#elif defined(USE_SHMEM)
        call mp_barrier( comm )
#endif
        return
!EOC
      end subroutine mp_exit
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: omp_start --- Start openMP parallelism
!
! !INTERFACE:
      subroutine omp_start
! !DESCRIPTION:
!
!     Start openMP parallelism
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer ios, n, nowpro, nowcpu

! Compute number of OpenMP threads

#if defined(_OPENMP)

        integer omp_get_num_threads
!$omp parallel
        numcpu = omp_get_num_threads()
!$omp end parallel

#else
        numcpu = 1
#endif

#if defined(USE_MPI2) || defined(USE_SHMEM)
# if defined(MT_OFF)
        pkgs_per_pro = 1
# else
        pkgs_per_pro = numcpu
# endif
#endif

! Old SGI coding below - may not work (AAM - 5/2006)
#if defined(PIN_CPUS)
!$omp parallel do private(n,nowcpu)
        nowpro = gid
        do n=1,numcpu
          nowcpu = n + (nowpro) * numcpu-1
          call mp_assign_to_cpu(nowcpu)
        enddo
#endif
!EOC
      end subroutine omp_start
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mpi_start --- Start MPI parallelism
!
! !INTERFACE:
      subroutine mpi_start( comm )
! !INPUT PARAMETERS:
      integer :: comm      !  communicator
! !DESCRIPTION:
!
!     Start MPI parallelism
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!    02.08.06   Sawyer  Added communicator input arguments
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        logical flag
        integer npthreads

        call MPI_INITIALIZED( flag, ierror )
        if ( .not. flag ) then
          call MPI_INIT( ierror )
          comm = MPI_COMM_WORLD
        endif

! Check for multiple threading capability with MPI-2
#if defined(USE_MPI2)
        call MPI_QUERY_THREAD(npthreads, ierror)
# if !defined(MT_OFF)
        if (npthreads /= MPI_THREAD_MULTIPLE) then
          write(iulog,*) gid, 'did not provide MPI_THREAD_MULTIPLE. ', &
                'Change to MPI_THREAD_MULTIPLE with MPI_INIT_THREAD ', &
                'for multi-threading MPI2, or define MT_OFF - exiting'
          call MPI_FINALIZE(ierror)
          call exit(1)
        endif
# endif
! Don't allow MPI-2 under Linux
# if defined(LINUX)
        write(iulog,*) 'USE_MPI2 not supported under LINUX - exiting'
        call MPI_FINALIZE(ierror)
        call exit(1)
# endif
#endif

! Shmem is meant to go with UNICOSMP and is not multi-threaded
#if defined (USE_SHMEM)
# if defined(UNICOSMP)
#  if !defined(MT_OFF)
        write(iulog,*) 'USE_SHMEM requires MT_OFF - exiting'
        call MPI_FINALIZE(ierror)
        call exit(1)
#  endif
# else
        write(iulog,*) 'USE_SHMEM requires UNICOSMP - exiting'
        call MPI_FINALIZE(ierror)
        call exit(1)
# endif
#endif

        call MPI_COMM_RANK (comm, gid, ierror)
        call MPI_COMM_SIZE (comm, numpro, ierror)
        call MPI_COMM_DUP  (comm, commglobal, ierror)
!EOC
      end subroutine mpi_start
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: win_init_r8 --- Initialize real*8 communication window
!
! !INTERFACE:
      subroutine win_init_r8(comm, win, ga, isize)
! !INPUT PARAMETERS:
        integer, intent(in) :: comm      !  communicator
        integer, intent(in) :: isize
        real(r8), intent(in) :: ga(isize)
! !OUTPUT PARAMETERS:
        type (window), intent(inout) :: win
! !DESCRIPTION:
!
!     Initialize real*8 communication window
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer mp_size, info

        win_count = win_count + 1

#if defined(USE_MPI2)
        call MPI_TYPE_SIZE(mp_r8, mp_size, ierror)
        bsize = isize*mp_size
        info = MPI_INFO_NULL
        call MPI_WIN_CREATE(ga, bsize, mp_size, info, comm, &
                            win%id, ierror)
# if defined(UNICOSMP)
! Cray X1E insufficiency
        call MPI_WIN_FENCE(0, win%id, ierror)
# else
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, win%id, ierror)
# endif
#else
        win%id = win_count
#endif
        win%size = isize
        win%ncall_s = 0
        win%ncall_r = 0
!EOC
      end subroutine win_init_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: win_init_r4 --- Initialize real*4 communication window
!
! !INTERFACE:
      subroutine win_init_r4(comm, win, ga, isize)
! !INPUT PARAMETERS:
        integer, intent(in) :: comm      !  communicator
        integer, intent(in) :: isize
        real(r4), intent(in) :: ga(isize)
! !OUTPUT PARAMETERS:
        type (window), intent(inout) :: win
! !DESCRIPTION:
!
!     Initialize real*4 communication window
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer mp_size, info

        win_count = win_count + 1

#if defined(USE_MPI2)
        call MPI_TYPE_SIZE(mp_r4, mp_size, ierror)
        bsize = isize*mp_size
        info = MPI_INFO_NULL
        call MPI_WIN_CREATE(ga, bsize, mp_size, info, comm, &
                            win%id, ierror)
# if defined(UNICOSMP)
! Cray X1E insufficiency
        call MPI_WIN_FENCE(0, win%id, ierror)
# else
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, win%id, ierror)
# endif
#else
        win%id = win_count
#endif
        win%size = isize
        win%ncall_s = 0
        win%ncall_r = 0
!EOC
      end subroutine win_init_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: win_init_i4 --- Initialize integer*4 communication window
!
! !INTERFACE:
      subroutine win_init_i4(comm, win, ga, isize)
! !INPUT PARAMETERS:
        integer, intent(in) :: comm      !  communicator
        integer, intent(in) :: isize
        integer(i4), intent(in) :: ga(isize)
! !OUTPUT PARAMETERS:
        type (window), intent(inout) :: win
! !DESCRIPTION:
!
!     Initialize integer*4 communication window
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer mp_size, info

        win_count = win_count + 1

#if defined(USE_MPI2)
        call MPI_TYPE_SIZE(mp_i4, mp_size, ierror)
        bsize = isize*mp_size
        info = MPI_INFO_NULL
        call MPI_WIN_CREATE(ga, bsize, mp_size, info, comm, &
                            win%id, ierror)
# if defined(UNICOSMP)
! Cray X1E insufficiency
        call MPI_WIN_FENCE(0, win%id, ierror)
# else
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, win%id, ierror)
# endif
#else
        win%id = win_count
#endif
        win%size = isize
        win%ncall_s = 0
        win%ncall_r = 0
!EOC
      end subroutine win_init_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send4d_ns --- Send 4d north/south ghost latitudes (real*8)
!
! !INTERFACE:
      subroutine mp_send4d_ns(comm, im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                              ng_s, ng_n, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_s      ! southern zones to ghost 
      integer, intent(in):: ng_n      ! northern zones to ghost 
      real(r8), intent(in):: q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Send 4d north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman        Modified for Global Arrays code   
!    2002.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: gidu

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_RANK (comm, gidu, ierror)

      call Win_Open(comm, t1_win)

! Send to south
      if ( jfirst > 1 ) then
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
        t1_win%src = gidu - 1
        t1_win%offset_r = ifromsouth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        t1_win%size_r = im*ng_s*(klast-kfirst+1)*nq
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
#endif
        t1_win%dest = gidu - 1
        t1_win%offset_s = igosouth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8(comm, q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jfirst, jfirst+ng_n-1, kfirst, klast, 1, nq,   &
                         ga_t1_s, ga_t1_r )
      endif
! Send to north
      if ( jlast < jm ) then
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
        t1_win%src = gidu + 1
        t1_win%offset_r = ifromnorth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        t1_win%size_r = im*ng_n*(klast-kfirst+1)*nq
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
#endif
        t1_win%dest = gidu + 1
        t1_win%offset_s = igonorth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8(comm, q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jlast-ng_s+1, jlast, kfirst, klast, 1, nq,     &
                         ga_t1_s, ga_t1_r )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send4d_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv4d_ns --- Receive 4d north/south ghost latitudes (real*8)
!
! !INTERFACE:
      subroutine mp_recv4d_ns(comm, im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                              ng_s, ng_n, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_s      ! southern zones to ghost 
      integer, intent(in):: ng_n      ! northern zones to ghost 
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Receive 4d north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman        Modified for Global Arrays code   
!    2002.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: gidu

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_RANK (comm, gidu, ierror)

      call Win_Close(comm, t1_win)


! Recv from south
      if ( jfirst > 1 ) then
        t1_win%src  = gidu-1
        t1_win%offset_r = ifromsouth*idimsizz + (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8(comm, q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jfirst-ng_s, jfirst-1,   kfirst, klast, 1, nq, &
                         ga_t1_r  )
      endif
! Recv from north
      if ( jlast < jm ) then
        t1_win%src  = gidu+1
        t1_win%offset_r = ifromnorth*idimsizz + (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8(comm, q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jlast+1,     jlast+ng_n, kfirst, klast, 1, nq, &
                         ga_t1_r  )
      endif

      call Win_Finalize(comm, t1_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv4d_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send4d_ns_r4 --- Send 4d north/south ghost latitudes (real*4)
!
! !INTERFACE:
      subroutine mp_send4d_ns_r4(comm, im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                                 ng_s, ng_n, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_s      ! southern zones to ghost 
      integer, intent(in):: ng_n      ! northern zones to ghost 
      real(r4), intent(in):: q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Send 4d north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2005.03.20   Sawyer        Creation from mp_send4d_ns
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: gidu

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: mp_send4d_ns_r4 - r4 windows disabled - exiting'
        call exit(1)
#endif

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_RANK (comm, gidu, ierror)

      call Win_Open(comm, r4_win)

! Send to south
      if ( jfirst > 1 ) then
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
        r4_win%src = gidu - 1
        r4_win%offset_r = ifromsouth*idimsizz + (r4_win%ncall_s-1)*idimsizz*nbuf
        r4_win%size_r = im*ng_s*(klast-kfirst+1)*nq
        call Ga_RecvInit_r4(comm, r4_win, ga_r4_r)
#endif
        r4_win%dest = gidu - 1
        r4_win%offset_s = igosouth*idimsizz + (r4_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r4(comm, q, r4_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jfirst, jfirst+ng_n-1, kfirst, klast, 1, nq,   &
                         ga_r4_s, ga_r4_r )
      endif
! Send to north
      if ( jlast < jm ) then
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
        r4_win%src = gidu + 1
        r4_win%offset_r = ifromnorth*idimsizz + (r4_win%ncall_s-1)*idimsizz*nbuf
        r4_win%size_r = im*ng_n*(klast-kfirst+1)*nq
        call Ga_RecvInit_r4(comm, r4_win, ga_r4_r)
#endif
        r4_win%dest = gidu + 1
        r4_win%offset_s = igonorth*idimsizz + (r4_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r4(comm, q, r4_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jlast-ng_s+1, jlast, kfirst, klast, 1, nq,     &
                         ga_r4_s, ga_r4_r )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send4d_ns_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv4d_ns_r4 --- Receive 4d north/south ghost latitudes (real*4)
!
! !INTERFACE:
      subroutine mp_recv4d_ns_r4(comm, im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                              ng_s, ng_n, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_s      ! southern zones to ghost 
      integer, intent(in):: ng_n      ! northern zones to ghost 
! !OUTPUT PARAMETERS:
      real(r4), intent(inout):: q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Receive 4d north/south ghost latitudes (real*4)
!
! !REVISION HISTORY: 
!    2005.03.20   Sawyer        Creation from mp_recv4d_ns
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: gidu

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: mp_recv4d_ns_r4 - r4 windows disabled - exiting'
        call exit(1)
#endif

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_RANK (comm, gidu, ierror)

      call Win_Close(comm, r4_win)


! Recv from south
      if ( jfirst > 1 ) then
        r4_win%src  = gidu-1
        r4_win%offset_r = ifromsouth*idimsizz + (r4_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r4(comm, q, r4_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jfirst-ng_s, jfirst-1,   kfirst, klast, 1, nq, &
                         ga_r4_r  )
      endif
! Recv from north
      if ( jlast < jm ) then
        r4_win%src  = gidu+1
        r4_win%offset_r = ifromnorth*idimsizz + (r4_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r4(comm, q, r4_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jlast+1,     jlast+ng_n, kfirst, klast, 1, nq, &
                         ga_r4_r  )
      endif

      call Win_Finalize(comm, r4_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv4d_ns_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send2_ns --- Send 2 variables north/south ghost latitudes
!
! !INTERFACE:
      subroutine mp_send2_ns(comm, im, jm, km, jfirst, jlast, kfirst, klast, &
                             nd, q1, q2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast !careful: klast might be klast+1
      integer, intent(in):: nd
      real(r8), intent(in):: q1(im,jfirst-nd:jlast+nd,kfirst:klast) 
      real(r8), intent(in):: q2(im,jfirst-nd:jlast+nd,kfirst:klast) 
!
! !DESCRIPTION:
!
!     Send 2 variables north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman        Modified for Global Arrays code   
!    2002.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer :: gidu

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_RANK (comm, gidu, ierror)

      call Win_Open(comm, t1_win)

! Send to south
      if ( jfirst > 1 ) then
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
        t1_win%src  = gidu - 1
        t1_win%size_r = im*(klast-kfirst+1)
        t1_win%offset_r = ifromsouth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
#endif
        t1_win%dest = gidu - 1
        t1_win%offset_s = igosouth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8( comm, q1, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, jfirst,    jfirst,   kfirst, klast, 1, 1, &
                          ga_t1_s, ga_t1_r  )
        t1_win%offset_s = t1_win%offset_s + im*(klast-kfirst+1)
        call Ga_Put4d_r8( comm, q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, jfirst,    jfirst,   kfirst, klast, 2, 2, &
                          ga_t1_s, ga_t1_r  )
      endif
! Send to north
      if ( jlast < jm ) then
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
        t1_win%src  = gidu + 1
        t1_win%size_r = im*(klast-kfirst+1)
        t1_win%offset_r = ifromnorth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
#endif
        t1_win%dest = gidu + 1
        t1_win%offset_s = igonorth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8( comm, q1, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, jlast,     jlast,    kfirst, klast, 1, 1, &
                          ga_t1_s, ga_t1_r  )
        t1_win%offset_s = t1_win%offset_s + im*(klast-kfirst+1)
        call Ga_Put4d_r8( comm, q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, jlast,     jlast,    kfirst, klast, 2, 2, &
                          ga_t1_s, ga_t1_r  )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send2_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv2_ns --- Receive 2 variables north/south ghost latitudes
!
! !INTERFACE:
      subroutine mp_recv2_ns(comm, im, jm, km, jfirst, jlast, kfirst, klast, &
                             nd, q1, q2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast !careful: klast might be klast+1
      integer, intent(in):: nd
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: q1(im,jfirst-nd:jlast+nd,kfirst:klast) 
      real(r8), intent(inout):: q2(im,jfirst-nd:jlast+nd,kfirst:klast)
!
! !DESCRIPTION:
!
!     Receive 2 variables north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman        Modified for Global Arrays code   
!    2002.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer j
      integer :: gidu

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_RANK (comm, gidu, ierror)

      call Win_Close(comm, t1_win)

! Recv from south
      if ( jfirst > 1 ) then
        j = jfirst - 1
        t1_win%src  = gidu - 1
        t1_win%offset_r = ifromsouth*idimsizz + (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8( comm, q1, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, j,         j,        kfirst, klast, 1, 1, &
                          ga_t1_r  )
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_Get4d_r8( comm, q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, j,         j,        kfirst, klast, 2, 2, &
                          ga_t1_r  )
      endif
! Recv from north
      if ( jlast < jm ) then
        j = jlast + 1
        t1_win%src  = gidu + 1
        t1_win%offset_r = ifromnorth*idimsizz + (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8( comm, q1, t1_win, im, jm, km, 2, & 
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, j,         j,        kfirst, klast, 1, 1, &
                          ga_t1_r  )
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_Get4d_r8( comm, q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, j,         j,        kfirst, klast, 2, 2, &
                          ga_t1_r  )
      endif

      call Win_Finalize(comm, t1_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv2_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send3d --- Send ghost region
!
! !INTERFACE:
      subroutine mp_send3d(comm, dest, src, im, jm, km, if, il, jf, jl, kf, kl, &
                                                  i1, i2, j1, j2, k1, k2, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: dest, src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl 
      integer, intent(in):: i1, i2, j1, j2, k1, k2
      real(r8), intent(in):: q(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Send a general 3d real*8 ghost region
!
! !REVISION HISTORY:
!    02.04.15   Putman  
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      call Win_Open(comm, t1_win)

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
! Init Recv src
      if ( src >= 0 .and. src < numcomm ) then     ! is PE in valid range?
        t1_win%src = src
        t1_win%size_r = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)  ! chunk size
        t1_win%offset_r = (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
      endif
#endif
! Send ghost region
      if ( dest >= 0 .and. dest < numcomm ) then
        t1_win%dest = dest
        t1_win%offset_s = (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8( comm, q, t1_win, im, jm, km, 1, &
                          if, il, jf, jl, kf, kl, 1, 1,  &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_s, ga_t1_r  )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv3d --- Recv ghost region
!
! !INTERFACE:
      subroutine mp_recv3d(comm, src, im, jm, km, if, il, jf, jl, kf, kl, &
                                                i1, i2, j1, j2, k1, k2, qout)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl 
      integer, intent(in):: i1, i2, j1, j2, k1, k2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: qout(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Recv a general 3d real*8 ghost region
!
! !REVISION HISTORY:
!    02.04.15   Putman  
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      call Win_Close(comm, t1_win)

! Recv from src
      if ( src >= 0 .and. src < numcomm ) then     ! is PE in valid range?
        t1_win%src  = src
        t1_win%offset_r = (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8( comm, qout, t1_win, im, jm, km, 1, &
                          if, il, jf, jl, kf, kl, 1, 1,   &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_r  )
      endif

      call Win_Finalize(comm, t1_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send3d_2 --- Send 2 ghost regions
!
! !INTERFACE:
      subroutine mp_send3d_2(comm, dest, src, im, jm, km, if, il, jf, jl, kf, kl, &
                                        i1, i2, j1, j2, k1, k2, q1, q2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: dest, src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl
      integer, intent(in):: i1, i2, j1, j2, k1, k2
      real(r8), intent(in):: q1(if:il, jf:jl, kf:kl)
      real(r8), intent(in):: q2(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Send two general 3d real*8 ghost region
!
! !REVISION HISTORY:
!    02.04.15   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      call Win_Open(comm, t1_win)

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
! Init Recv src
      if ( src >= 0 .and. src < numcomm ) then     ! is PE in valid range?
        t1_win%src = src
        t1_win%size_r = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)  ! chunk size
        t1_win%offset_r = (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
        t1_win%offset_r = t1_win%offset_r + t1_win%size_r 
        call Ga_RecvInit_r8(comm, t1_win, ga_t1_r)
      endif
#endif
! Send ghost region
      if ( dest >= 0 .and. dest < numcomm ) then
        t1_win%dest = dest
        t1_win%offset_s = (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8( comm, q1, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 1, 1,  &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_s, ga_t1_r  )
        t1_win%offset_s = t1_win%offset_s + (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
        call Ga_Put4d_r8( comm, q2, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 2, 2,  &
                          i1, i2, j1, j2, k1, k2, 2, 2, ga_t1_s, ga_t1_r  )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send3d_2
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv3d_2 --- Recv 2 ghost regions
!
! !INTERFACE:
      subroutine mp_recv3d_2(comm, src, im, jm, km, if, il, jf, jl, kf, kl, &
                                  i1, i2, j1, j2, k1, k2, qout1, qout2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: comm      !  communicator
      integer, intent(in):: src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl
      integer, intent(in):: i1, i2, j1, j2, k1, k2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: qout1(if:il, jf:jl, kf:kl)
      real(r8), intent(inout):: qout2(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Recv two general 3d real*8 ghost regions
!
! !REVISION HISTORY:
!    02.04.15   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      call Win_Close(comm, t1_win)

! Recv from src
      if ( src >= 0 .and. src < numcomm ) then     ! is PE in valid range?
        t1_win%src  = src
        t1_win%offset_r = (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8( comm, qout1, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 1, 1,   &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_r  )
        t1_win%offset_r = t1_win%offset_r + (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
        call Ga_Get4d_r8( comm, qout2, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 2, 2,   &
                          i1, i2, j1, j2, k1, k2, 2, 2, ga_t1_r  )
      endif

      call Win_Finalize(comm, t1_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv3d_2
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_barrier --- Synchronize all SPMD processes
!
! !INTERFACE:
      subroutine mp_barrier (comm)
!
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
! !DESCRIPTION:
!
!     Synchronize all SPMD processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined (USE_SHMEM)
        call MPI_COMM_SIZE (comm, numcomm, ierror)
        call SHMEM_BARRIER(masterpro, 0, numcomm, barrier_sync)
#else
        call MPI_BARRIER(comm, ierror)
#endif

!EOC
      end subroutine mp_barrier
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Win_Open --- Open a communication window
!
! !INTERFACE:
      subroutine Win_Open(comm, win)
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
! !OUTPUT PARAMETERS:
      type(window), intent(inout):: win
!
! !DESCRIPTION:
!
!     Begin a communication epoch, by opening a comm window.
!     Update number of send calls on the window (win%ncall_s).
!     Barrier synchronzize if necessary.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      win%ncall_s = win%ncall_s + 1
      ncall_s = ncall_s + 1

#if defined(USE_SHMEM)
      if (ncall_s == 1) then
         if (ga_ptr == 1) then 
            ga_ptr = 2
# if defined(MODCM_STATIC)
            ga_r8_r => ga_r8b_r
            ga_r8_s => ga_r8b_s
            ga_r4_r => ga_r4b_r
            ga_r4_s => ga_r4b_s
            ga_i4_r => ga_i4b_r
            ga_i4_s => ga_i4b_s
# else
            pa_r8_r = pa_r8b_r
            pa_r8_s = pa_r8b_s
            pa_r4_r = pa_r4b_r
            pa_r4_s = pa_r4b_s
            pa_i4_r = pa_i4b_r
            pa_i4_s = pa_i4b_s
# endif
         else
            ga_ptr = 1
# if defined(MODCM_STATIC)
            ga_r8_r => ga_r8a_r
            ga_r8_s => ga_r8a_s
            ga_r4_r => ga_r4a_r
            ga_r4_s => ga_r4a_s
            ga_i4_r => ga_i4a_r
            ga_i4_s => ga_i4a_s
# else
            pa_r8_r = pa_r8a_r
            pa_r8_s = pa_r8a_s
            pa_r4_r = pa_r4a_r
            pa_r4_s = pa_r4a_s
            pa_i4_r = pa_i4a_r
            pa_i4_s = pa_i4a_s
# endif
         endif
      endif
#endif

!EOC
      end subroutine Win_Open
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Win_Close --- Close a communication window
!
! !INTERFACE:
      subroutine Win_Close(comm, win)
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
! !OUTPUT PARAMETERS:
      type(window), intent(inout):: win
!
! !DESCRIPTION:
!
!     End a communication epoch, by closing a comm window.
!     Update number of receive calls on the window (win%ncall_r).
!     Barrier synchronzize if necessary.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      win%ncall_r = win%ncall_r + 1
      ncall_r = ncall_r + 1
#if defined(USE_SHMEM)
      if (ncall_r == 1) then
          call mp_barrier(comm)
      endif
#endif
#if defined(USE_MPI2)
      if (win%ncall_r == 1) then
#if defined(UNICOSMP)
          call MPI_WIN_FENCE(0, win%id, ierror)
#else
          call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                             win%id, ierror)
#endif
      endif
#endif

!EOC
      end subroutine Win_Close
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Win_Finalize --- Reset a communication window after a comm epoch.
!
! !INTERFACE:
      subroutine Win_Finalize(comm, win)
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
! !OUTPUT PARAMETERS:
      type(window), intent(inout):: win
!
! !DESCRIPTION:
!
!     Complete a communication epoch and reset a comm window.
!     Update global lastwin with win%id.
!     Barrier synchronzize if necessary.
!
! !REVISION HISTORY:
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if (ncall_s == ncall_r) then
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
        call MPI_WAITALL(nsend, sqest, Stats, ierror)
        nsend = 0
        nrecv = 0
        nread = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif

      if (win%ncall_s == win%ncall_r) then
#if defined(USE_MPI2)
# if defined(UNICOSMP)
        call MPI_WIN_FENCE(0, win%id, ierror)
# else
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, win%id, ierror)
# endif
#endif
        lastwin = win%id
        win%ncall_s = 0
        win%ncall_r = 0
      endif

!EOC
      end subroutine Win_Finalize
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Put4d_r8 --- Write to real*8 4d global array
!
! !INTERFACE:
      subroutine Ga_Put4d_r8 ( comm, q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga_s, ga_r )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(in)  :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r8), intent(in)  :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: ga_s(win%size)
      real(r8), intent(inout):: ga_r(win%size)
!
! !DESCRIPTION:
!
!     Write to real*8 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

      integer send_tag, qsize
#if defined(USE_MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#if defined(USE_SHMEM)
      integer p, tmpsize, mysize, mydisp
#endif
      integer :: gidu

      call MPI_COMM_RANK (comm, gidu, ierror)

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length 
      ijk_length = i_length*j_length*k_length

! Begin Non-Blocking Sends ( MPI-2,MPI-1 )
        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_s) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                ga_s(inc+i) = q(i,j,k,iq)
              enddo
            enddo
          enddo
        enddo

      qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*(nq2-nq1+1)
#if defined(USE_MPI2) || defined(USE_SHMEM)
      tmpsize = ceiling(real(qsize,r8)/real(pkgs_per_pro,r8))
!$omp parallel do private(p,mysize,mydisp,ierror)
      do p=1,MIN(pkgs_per_pro,qsize)
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        mydisp = win%offset_s + (p-1)*tmpsize
# if defined(USE_MPI2)
        call MPI_PUT(ga_s(mydisp+1), mysize, mp_r8, &
                     win%dest, mydisp, mysize, mp_r8, &
                     win%id, ierror)
# else
           call SHMEM_PUT64(ga_r(mydisp+1), ga_s(mydisp+1), &
                                 mysize, win%dest)
# endif
      enddo
#else
      send_tag = gidu
      nsend = nsend + 1
      call MPI_ISEND(ga_s(win%offset_s+1), qsize, mp_r8, win%dest, &
                     send_tag, comm, sqest(nsend), ierror)
#endif

!EOC
      end subroutine Ga_Put4d_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_RecvInit_r8 --- Initiate real*8 Non-Blocking receive
!
! !INTERFACE:
      subroutine Ga_RecvInit_r8( comm, win, ga )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(in)    :: win           ! Global Array Window
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: ga(win%size)
!
! !DESCRIPTION:
!
!     Initiate real*8 Non-Blocking receive ( MPI-1 ).
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!    03.06.06   Sawyer        Added else clause
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer qsize, recv_tag

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
      if (win%size >= win%offset_r + win%size_r) then
        recv_tag = win%src
        qsize    = win%size_r
        nrecv    = nrecv + 1
        call MPI_IRECV(ga(win%offset_r+1), qsize, mp_r8, win%src, &
                       recv_tag, comm, rqest(nrecv), ierror)
      else
        write(iulog,*) "Fatal ga_recvinit_r8: receive window out of space - exiting"
        write(iulog,*) 'gid win%size win%offset_r win%size_r = ', gid,  &
                  win%size, win%offset_r, win%size_r
        call exit(1)
      endif
#endif

!EOC
      end subroutine Ga_RecvInit_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Get4d_r8 --- Read from real*8 4d global array
!
! !INTERFACE:
      subroutine Ga_Get4d_r8 ( comm, q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(in)    :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r8), intent(in)  :: ga(win%size)
! !OUTPUT PARAMETERS:
      real(r8), intent(inout) :: q(ifrom:ito, jfrom:jto, kfrom:kto, nqfrom:nqto)
!
! !DESCRIPTION:
!
!     Read from real*8 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
      nread = nread + 1
      call MPI_WAIT(rqest(nread), Status, ierror)
#endif

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_r) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                q(i,j,k,iq) = ga(inc+i)
              enddo
            enddo
          enddo
        enddo

!EOC
      end subroutine Ga_Get4d_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Put4d_r4 --- Write to real*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Put4d_r4 ( comm, q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga_s, ga_r )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(in)  :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r4), intent(in)  :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
! !OUTPUT PARAMETERS:
      real(r4), intent(inout):: ga_s(win%size)
      real(r4), intent(inout):: ga_r(win%size)
!
! !DESCRIPTION:
!
!     Write to real*4 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

      integer send_tag, qsize
#if defined(USE_MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#if defined(USE_SHMEM)
      integer p, tmpsize, mysize, mydisp
#endif
      integer :: gidu

      call MPI_COMM_RANK (comm, gidu, ierror)

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: Ga_Put4d_r4 - r4 windows disabled - exiting'
        call exit(1)
#endif

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length 
      ijk_length = i_length*j_length*k_length

! Begin Non-Blocking Sends ( MPI-2,MPI-1 )
        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_s) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                ga_s(inc+i) = q(i,j,k,iq)
              enddo
            enddo
          enddo
        enddo

      qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*(nq2-nq1+1)
#if defined(USE_MPI2) || defined(USE_SHMEM)
      tmpsize = ceiling(real(qsize,r8)/real(pkgs_per_pro,r8))
!$omp parallel do private(p,mysize,mydisp,ierror)
      do p=1,MIN(pkgs_per_pro,qsize)
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        mydisp = win%offset_s + (p-1)*tmpsize
# if defined(USE_MPI2)
        call MPI_PUT(ga_s(mydisp+1), mysize, mp_r4, &
                     win%dest, mydisp, mysize, mp_r4, &
                     win%id, ierror)
# else
           call SHMEM_PUT32(ga_r(mydisp+1), ga_s(mydisp+1), &
                                 mysize, win%dest)
# endif
      enddo
#else
      send_tag = gidu
      nsend = nsend + 1
      call MPI_ISEND(ga_s(win%offset_s+1), qsize, mp_r4, win%dest, &
                     send_tag, comm, sqest(nsend), ierror)
#endif

!EOC
      end subroutine Ga_Put4d_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_RecvInit_r4 --- Initiate real*4 Non-Blocking receive
!
! !INTERFACE:
      subroutine Ga_RecvInit_r4( comm, win, ga )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(in)    :: win           ! Global Array Window
! !OUTPUT PARAMETERS:
      real(r4), intent(inout):: ga(win%size)
!
! !DESCRIPTION:
!
!     Initiate real*8 Non-Blocking receive ( MPI-1 ).
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!    03.06.06   Sawyer        Added else clause
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer qsize, recv_tag

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: Ga_RecvInit_r4 - r4 windows disabled - exiting'
        call exit(1)
#endif

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
      if (win%size >= win%offset_r + win%size_r) then
        recv_tag = win%src
        qsize    = win%size_r
        nrecv    = nrecv + 1
        call MPI_IRECV(ga(win%offset_r+1), qsize, mp_r4, win%src, &
                       recv_tag, comm, rqest(nrecv), ierror)
      else
        write(iulog,*) "Fatal ga_recvinit_r4: receive window out of space - exiting"
        write(iulog,*) 'gid win%size win%offset_r win%size_r = ', gid,  &
                  win%size, win%offset_r, win%size_r
        call exit(1)
      endif
#endif

!EOC
      end subroutine Ga_RecvInit_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Get4d_r4 --- Read from real*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Get4d_r4 ( comm, q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(in)    :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r4), intent(in)  :: ga(win%size)
! !OUTPUT PARAMETERS:
      real(r4), intent(inout) :: q(ifrom:ito, jfrom:jto, kfrom:kto, nqfrom:nqto)
!
! !DESCRIPTION:
!
!     Read from real*8 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: Ga_Get4d_r4 - r4 windows disabled - exiting'
        call exit(1)
#endif

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
      nread = nread + 1
      call MPI_WAIT(rqest(nread), Status, ierror)
#endif

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_r) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                q(i,j,k,iq) = ga(inc+i)
              enddo
            enddo
          enddo
        enddo

!EOC
      end subroutine Ga_Get4d_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Put4d_i4 --- Write to integer*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Put4d_i4 ( comm, q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga_s, ga_r )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(in)  :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      integer(i4), intent(in)  :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout):: ga_s(win%size)
      integer(i4), intent(inout):: ga_r(win%size)
!
! !DESCRIPTION:
!
!     Write to integer*4 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

      integer send_tag, qsize
#if defined(USE_MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#if defined(USE_SHMEM)
      integer p, tmpsize, mysize, mydisp
#endif
      integer :: gidu

      call MPI_COMM_RANK (comm, gidu, ierror)

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

! Begin Non-Blocking Sends ( MPI-2,MPI-1 )
        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_s) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                ga_s(inc+i) = q(i,j,k,iq)
              enddo
            enddo
          enddo
        enddo

      qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*(nq2-nq1+1)
#if defined(USE_MPI2) || defined(USE_SHMEM)
      tmpsize = ceiling(real(qsize,r8)/real(pkgs_per_pro,r8))
!$omp parallel do private(p,mysize,mydisp,ierror)
      do p=1,MIN(pkgs_per_pro,qsize)
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        mydisp = win%offset_s + (p-1)*tmpsize
# if defined(USE_MPI2)
        call MPI_PUT(ga_s(mydisp+1), mysize, mp_i4, &
                     win%dest, mydisp, mysize, mp_i4, &
                     win%id, ierror)
# else
           call SHMEM_PUT32(ga_r(mydisp+1), ga_s(mydisp+1), &
                                 mysize, win%dest)
# endif
      enddo
#else
      send_tag = gidu
      nsend = nsend + 1
      call MPI_ISEND(ga_s(win%offset_s+1), qsize, mp_i4, win%dest, &
                     send_tag, comm, sqest(nsend), ierror)
#endif

!EOC
      end subroutine Ga_Put4d_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_RecvInit_i4 --- Initiate integer*4 Non-Blocking receive
!
! !INTERFACE:
      subroutine Ga_RecvInit_i4( comm, win, ga )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(in)    :: win           ! Global Array Window
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout):: ga(win%size)
!
! !DESCRIPTION:
!
!     Initiate integer*4 Non-Blocking receive ( MPI-1 ).
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!    06.05.21   Mirin         Added else clause
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer qsize, recv_tag

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
      if (win%size >= win%offset_r + win%size_r) then
        recv_tag = win%src
        qsize    = win%size_r
        nrecv    = nrecv + 1
        call MPI_IRECV(ga(win%offset_r+1), qsize, mp_i4, win%src, &
                       recv_tag, comm, rqest(nrecv), ierror)
      else
        write(iulog,*) "Fatal ga_recvinit_i4: receive window out of space - exiting"
        write(iulog,*) 'gid win%size win%offset_r win%size_r = ', gid,  &
                  win%size, win%offset_r, win%size_r
        call exit(1)
      endif
#endif
!EOC
      end subroutine Ga_RecvInit_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Get4d_i4 --- Read from integer*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Get4d_i4 ( comm, q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      type(window), intent(in)    :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      integer(i4), intent(in)  :: ga(win%size)
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: q(ifrom:ito, jfrom:jto, kfrom:kto, nqfrom:nqto)
!
! !DESCRIPTION:
!
!     Read from integer*4 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman        Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
      nread = nread + 1
      call MPI_WAIT(rqest(nread), Status, ierror)
#endif

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1
      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_r) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                q(i,j,k,iq) = ga(inc+i)
              enddo
            enddo
          enddo
        enddo

!EOC
      end subroutine Ga_Get4d_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Broadcast_r8 --- Broadcast an real*8 1d global array
!
! !INTERFACE:
      subroutine Ga_Broadcast_r8 ( comm, q, isize )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      integer, intent(in)  :: isize
! !OUTPUT PARAMETERS:
      real(r8), intent(inout) :: q(isize)
!
! !DESCRIPTION:
!
!     Broadcast an real*8 1d global array.
!
! !REVISION HISTORY:
!    03.04.02        Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer i
#if defined(USE_SHMEM)
      integer p, tmpsize, mysize, mydisp
#endif

#if defined(USE_SHMEM)
      call MPI_COMM_SIZE (comm, numcomm, ierror)
      mysize = isize
      call mp_barrier(comm)
      do i=1,mysize
         ga_r8_s(i) = q(i)
      enddo
      call SHMEM_BROADCAST64(ga_r8_r, ga_r8_s, mysize, masterpro, masterpro,   &
                             0, numcomm, bcast_sync)
      if (gid /= 0) then
         do i=1,isize
            q(i) = ga_r8_r(i)
         enddo
      endif
#endif

#if !defined(USE_SHMEM)
      call MPI_BCAST(q, isize, mp_r8, 0, comm, ierror)
#endif

!EOC
      end subroutine Ga_Broadcast_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Broadcast_r4 --- Broadcast an real*4 1d global array
!
! !INTERFACE:
      subroutine Ga_Broadcast_r4 ( comm, q, isize )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      integer, intent(in)  :: isize
! !OUTPUT PARAMETERS:
      real(r4), intent(inout) :: q(isize)
!
! !DESCRIPTION:
!
!     Broadcast an real*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02        Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer i
#if defined(USE_SHMEM)
      integer p, tmpsize, mysize, mydisp
#endif

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: Ga_Broadcast_r4 - r4 windows disabled - exiting'
        call exit(1)
#endif

#if defined(USE_SHMEM)
      call MPI_COMM_SIZE (comm, numcomm, ierror)
      mysize = isize
      call mp_barrier(comm)
      do i=1,mysize
         ga_r4_s(i) = q(i)
      enddo
      call SHMEM_BROADCAST32(ga_r4_r, ga_r4_s, mysize, masterpro, masterpro,   &
                             0, numcomm, bcast_sync)
      if (gid /= 0) then
         do i=1,isize
            q(i) = ga_r4_r(i)
         enddo
      endif
#endif

#if !defined(USE_SHMEM)
      call MPI_BCAST(q, isize, mp_r4, 0, comm, ierror)
#endif

!EOC
      end subroutine Ga_Broadcast_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Broadcast_i4 --- Broadcast an integer*4 1d global array
!
! !INTERFACE:
      subroutine Ga_Broadcast_i4 ( comm, q, isize )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      integer, intent(in)  :: isize
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: q(isize)
!
! !DESCRIPTION:
!
!     Broadcast an integer*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02        Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer i
#if defined(USE_SHMEM)
      integer mysize
#endif

#if defined(USE_SHMEM)
      call MPI_COMM_SIZE (comm, numcomm, ierror)
      mysize = isize
      call mp_barrier(comm)
      do i=1,mysize
         ga_i4_s(i) = q(i)
      enddo
      call SHMEM_BROADCAST32(ga_i4_r, ga_i4_s, mysize, masterpro, masterpro,   &
                             0, numcomm, bcast_sync)
      if (gid /= 0) then
         do i=1,isize
            q(i) = ga_i4_r(i)
         enddo
      endif
#endif

#if !defined(USE_SHMEM)
      call MPI_BCAST(q, isize, mp_i4, 0, comm, ierror)
#endif

!EOC
      end subroutine Ga_Broadcast_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_AllToAll_r8 --- All to All of an real*8 1d global array
!
! !INTERFACE:
      subroutine Ga_AllToAll_r8 ( comm, q, Gsize, Lsize, istart )
! !INPUT PARAMETERS:
      integer, intent(in) :: comm      !  communicator
      integer, intent(in)  :: Gsize    ! Global size of array
      integer, intent(in)  :: Lsize    ! size of Local portion
      integer, intent(in)  :: istart   ! starting point
! !OUTPUT PARAMETERS:
      real(r8), intent(inout) :: q(Gsize)
!
! !DESCRIPTION:
!
!     All to All of a real*8 1d global array.
!
! !REVISION HISTORY:
!    03.04.02        Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer n, i
#if defined(USE_SHMEM)
      integer my_Lsize
#endif

#if defined(USE_SHMEM)
      call MPI_COMM_SIZE (comm, numcomm, ierror)
      my_Lsize = Lsize
      call mp_barrier(comm)
      do i=1,Lsize
         ga_r8_s(i) = q(i+istart-1)
      enddo
      call SHMEM_COLLECT8(ga_r8_r, ga_r8_s, my_Lsize, masterpro, 0, numcomm, collect_sync)
      do i=1,Gsize
         q(i) = ga_r8_r(i)
      enddo
#endif
#if !defined(USE_SHMEM)
      call MPI_ALLGATHER(q(istart), Lsize, mp_r8, q, Lsize, mp_r8, comm, ierror)
#endif

!EOC
      end subroutine Ga_AllToAll_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_AllToAll_r4 --- All to All of an real*4 1d global array
!
! !INTERFACE:
      subroutine Ga_AllToAll_r4 ( comm, q, Gsize, Lsize, istart )
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      integer, intent(in)  :: Gsize   ! Global size of array
      integer, intent(in)  :: Lsize   ! size of Local portion
      integer, intent(in)  :: istart  ! starting point
! !OUTPUT PARAMETERS:
      real(r4), intent(inout) :: q(Gsize)
!
! !DESCRIPTION:
!
!     All to All of an real*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02        Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer n, i
#if defined(USE_SHMEM)
      integer my_Lsize
#endif

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: Ga_AllToAll_r4 - r4 windows disabled - exiting'
        call exit(1)
#endif

#if defined(USE_SHMEM)
      call MPI_COMM_SIZE (comm, numcomm, ierror)
      my_Lsize = Lsize
      call mp_barrier(comm)
      do i=1,Lsize
         ga_r4_s(i) = q(i+istart-1)
      enddo
      call SHMEM_COLLECT4(ga_r4_r, ga_r4_s, my_Lsize, masterpro, 0, numcomm, collect_sync)
      do i=1,Gsize
         q(i) = ga_r4_r(i)
      enddo
#endif
#if !defined(USE_SHMEM)
      call MPI_ALLGATHER(q(istart), Lsize, mp_r4, q, Lsize, mp_r4, comm, ierror)
#endif

!EOC
      end subroutine Ga_AllToAll_r4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_AllToAll_i4 --- All to All of an integer*4 1d global array
!
! !INTERFACE:
      subroutine Ga_AllToAll_i4 ( comm, q, Gsize, Lsize, istart )
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      integer, intent(in)  :: Gsize   ! Global size of array
      integer, intent(in)  :: Lsize   ! size of Local portion
      integer, intent(in)  :: istart  ! starting point
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: q(Gsize)
!
! !DESCRIPTION:
!
!     All to All of an integer*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02        Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer n, i
#if defined(USE_SHMEM)
      integer my_Lsize
#endif

#if defined(USE_SHMEM)
      call MPI_COMM_SIZE (comm, numcomm, ierror)
      my_Lsize = Lsize
      call mp_barrier(comm)
      do i=1,Lsize
         ga_i4_s(i) = q(i+istart-1)
      enddo
      call SHMEM_COLLECT4(ga_i4_r, ga_i4_s, my_Lsize, masterpro, 0, numcomm, collect_sync)
      do i=1,Gsize
         q(i) = ga_i4_r(i)
      enddo
#endif
#if !defined(USE_SHMEM)
      call MPI_ALLGATHER(q(istart), Lsize, mp_i4, q, Lsize, mp_i4, comm, ierror)
#endif

!EOC
      end subroutine Ga_AllToAll_i4
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: get_partneroffset --- Computes partneroffset/type from descriptor
!
! !INTERFACE:
      subroutine get_partneroffset ( comm, send_bl, recv_bl )

! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
! !INPUT/OUTPUT PARAMETERS:
      type(blockdescriptor), intent(inout)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(inout)  :: recv_bl(:) ! receive blocks

!
! !DESCRIPTION:
!     Compute partneroffsets/types from other blockdescriptor
!     information.  Used exclusively for irregular communication 
!     in PILGRIM.
!
! !REVISION HISTORY: 
!    03.10.31   Mirin       Creation
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      integer :: i, j, k, ns, pos, por, numpsq, ierror
      integer :: ami(numpro,numpro), am(numpro,numpro)
      integer mod_method, mpsh_flag, num_s, num_r

      mpsh_flag = 0
#if defined (USE_MPI2) || defined (USE_SHMEM)
      mpsh_flag = 1
#endif

      num_s = size(send_bl)
      num_r = size(recv_bl)

      do j = 1, num_s
         send_bl(j)%partneroffset = 0
         send_bl(j)%partnertype = MPI_DATATYPE_NULL
      enddo
      do j = 1, num_r
         recv_bl(j)%partneroffset = 0
         recv_bl(j)%partnertype = MPI_DATATYPE_NULL
      enddo

    if (mpsh_flag .eq. 1) then

      mod_method = recv_bl(1)%method

! partner offset (dedicated window)
! am(i,j) is number of words sent from i to j (1-based indices)
! compute global table by using reduction on local parts
      ami(:,:) = 0
      am(:,:) = 0
      i = gid + 1
      numpsq = numcomm*numcomm
     if (i .le. num_r) then
      do j = 1, num_s
         ns = send_bl(j)%Nparcels
         do k = 1, ns
            ami(i,j) = ami(i,j) + send_bl(j)%blocksizes(k)
         enddo
      enddo
     endif

      call mpi_allreduce(ami, am, numpsq, mpi_integer, MPI_MAX, comm, ierror)

     if (i .le. num_r) then
      do j = 1, num_s
         send_bl(j)%partneroffset = 0
         pos = 0
         do k = 1, i-1
            pos = pos + am(k,j)
         enddo
         if (ami(i,j) .ne. 0) send_bl(j)%partneroffset = pos
      enddo
     endif
     if (i .le. num_s) then
      do j = 1, num_r
         recv_bl(j)%partneroffset = 0
         por = 0
         do k = 1, i-1
            por = por + am(j,k)
         enddo
         if (ami(j,i) .ne. 0) recv_bl(j)%partneroffset = por
      enddo
     endif

    endif  ! mpsh_flag .eq. 1

      end subroutine get_partneroffset
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_sendirr --- Write r8 contiguous parcels to global array
!
! !INTERFACE:
      subroutine mp_sendirr ( comm, q, send_bl, recv_bl, qout )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      real(r8), intent(in) :: q(*)                     ! local data segment

! !OUTPUT PARAMETERS:
      real(r8), intent(out) :: qout(*)                 ! local output segment
!
! !DESCRIPTION:
!     Send a number of contiguous parcels destined to an arbitrary set 
!     of PEs.  This is basically the non-blocking start of an
!     all-to-all communcation primitive.  This fundamental
!     routine forms a basis for higher level primitives.
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Use partneroffset
!    03.06.24   Sawyer      Integrated Use_Mpi_Types; added qout
!    04.02.24   Mirin       Various mpi2 options
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer ipe, qsize, offset, blocksize, nparcels, offset_s, offset_r, ierr, mod_method, mpi2_flag
      integer p, mysize, nthpc, minsize, nthrd, pn, pt, tmpsize
#if defined (USE_MPI2)
      integer (kind=MPI_ADDRESS_KIND) target_disp
#endif
      integer i, j, send_tag, recv_tag, num_s, num_r
      integer :: offset_v (Max_Nparcels)

!
!     initialize window

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

!     num_s = 0 if this processes is not part of the sending decomposition
      num_s = size(send_bl)
      if (send_bl(1)%Nparcels == -1) then
         num_s = 0
      endif

!     num_r = 0 if this processes is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      mpi2_flag = 0
#if defined (USE_MPI2)
      mpi2_flag = 1
#endif

      mod_method = recv_bl(1)%method
#if defined (USE_SHMEM)
      if (mod_method .ne. 0) then
         write(iulog,*) 'mp_sendirr with shmem - invalid mod_method = ', mod_method, ' - exiting'
         call exit(1)
      endif
#endif

      call Win_Open(comm, r8_win)
    if (mod_method .gt. 0) then

#if defined (USE_MPI2)
     if (mod_method .eq. 1) then
!  directly put contiguous segments into global target window
      do ipe = 1, num_s
         qsize = send_bl(ipe)%Tot_Size
         if (qsize .ne. 0) then
!  thread over parcels, but take into consideration number of
!    threads versus number of parcels and minimum size of parcel
            nparcels = send_bl(ipe)%Nparcels
            minsize = send_bl(ipe)%blocksizes(1)
            offset_v(1) = 0
            do p = 2, nparcels
                minsize = min(minsize, send_bl(ipe)%blocksizes(p))
                offset_v(p) = offset_v(p-1) + send_bl(ipe)%blocksizes(p-1)
            enddo
            r8_win%dest = ipe-1
            nthpc = ceiling(real(pkgs_per_pro,r8)/real(nparcels,r8))
            nthrd = min(nthpc,minsize)      !  number of threads per block
            r8_win%offset_r = send_bl(ipe)%partneroffset
!$omp parallel do private(pn,p,pt,tmpsize,mysize,offset_s,ierr)
            do pn = 1, nparcels*nthrd
               p = (pn-1)/nthrd + 1        ! block number
               pt = mod(pn-1,nthrd) + 1    ! thread number
               tmpsize = ceiling(real(send_bl(ipe)%blocksizes(p),r8)/real(nthrd,r8))
               mysize = min(tmpsize, max(send_bl(ipe)%blocksizes(p)-tmpsize*(pt-1),0))
               offset_s = send_bl(ipe)%displacements(p)
               call MPI_PUT(q(offset_s+(pt-1)*tmpsize+1), mysize,       &
                            mp_r8, r8_win%dest,                         &
                            r8_win%offset_r+offset_v(p)+(pt-1)*tmpsize,   &
                            mysize, mp_r8, r8_win%id, ierr)
            enddo
         endif
      enddo
     elseif (mod_method .eq. 2) then
! directly put derived types into global target window
!   thread over targets
!$omp parallel do private(ipe,mysize,ierr,target_disp)
        do ipe = 1, num_s
           if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
               mysize = send_bl(ipe)%Tot_Size
               target_disp = send_bl(ipe)%partneroffset
               call MPI_PUT(q, 1, send_bl(ipe)%type, ipe-1, target_disp,   &
                            mysize, mp_r8, r8_win%id, ierr)
           endif
        enddo
     endif
#elif defined(USE_SHMEM)
#else
!
! mpi-1 with derived types
! Increment the ongoing transfer number
      BegTrf = MOD(BegTrf,MAX_TRF) + 1
!
! MPI: Irecv over all processes
!
      OutHandle(:,BegTrf) = MPI_REQUEST_NULL
      do ipe=1, num_r

!
! Receive the buffers with MPI_Irecv. Non-blocking
!
        if ( recv_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_irecv( qout, 1, recv_bl(ipe)%type, ipe-1, ipe-1,     &
                          comm, OutHandle(ipe,BegTrf), ierr )
        endif
      enddo

!
! MPI: Isend over all processes
!
      InHandle(:,BegTrf) = MPI_REQUEST_NULL
      do ipe=1, num_s

!
! Send the individual buffers with non-blocking sends
!
        if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_isend( q, 1, send_bl(ipe)%type, ipe-1, gid,        &
                          comm, InHandle(ipe,BegTrf), ierr )
        endif
      enddo
#endif
    else

! temporary contiguous buffers
! for MPI-1, issue call to receive data in global receive buffer
      offset_s = 0
      offset_r = 0

      do ipe=1, num_r
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
         r8_win%size_r = recv_bl(ipe)%Tot_Size
         if (r8_win%size_r .ne. 0) then
            r8_win%offset_r = offset_r
            offset_r = offset_r + r8_win%size_r
            r8_win%src = ipe-1
            if (r8_win%size >= r8_win%offset_r + r8_win%size_r) then
              recv_tag = r8_win%src
              qsize    = r8_win%size_r
              nrecv    = nrecv + 1
              call MPI_IRECV(ga_r8_r(r8_win%offset_r+1), qsize, mp_r8, r8_win%src, &
                             recv_tag, comm, rqest(nrecv), ierror)
            else
              write(iulog,*) "Fatal mp_sendirr: receive window out of space - exiting"
              write(iulog,*) 'gid ipe r8_win%size r8_win%offset_r r8_win%size_r = ', gid,  &
                        ipe, r8_win%size, r8_win%offset_r, r8_win%size_r
              call exit(1)
            endif
         endif
#endif
      enddo
! gather data into global send buffer
      do ipe=1, num_s
         qsize = send_bl(ipe)%Tot_Size
         if (qsize .ne. 0) then
            r8_win%dest = ipe-1
            r8_win%offset_s = offset_s
            offset_s = offset_s + qsize
            if (offset_s .gt. r8_win%size) then
              write(iulog,*) "Fatal mp_sendirr: send window out of space - exiting"
              write(iulog,*) 'gid ipe offset_s r8_win%size = ', gid,  &
                        ipe, offset_s, r8_win%size
              call exit(1)
            endif

            offset_v(1) = r8_win%offset_s
            do j = 2, send_bl(ipe)%nparcels
               offset_v(j) = offset_v(j-1) + send_bl(ipe)%blocksizes(j-1)
            enddo

!dir$ concurrent
!dir$ preferstream
            do j = 1, send_bl(ipe)%nparcels
!dir$ concurrent
               do i = 1, send_bl(ipe)%blocksizes(j)
                  ga_r8_s(offset_v(j)+i) = q(send_bl(ipe)%displacements(j)+i)
               enddo
            enddo

#if defined (USE_MPI2) || defined (USE_SHMEM)
! one-sided puts for MPI-2 or SHMEM
            r8_win%offset_r = send_bl(ipe)%partneroffset
            if ( r8_win%offset_r+qsize > r8_win%size ) then
               write(iulog,*) "Fatal mp_sendirr: target window out of space - exiting"
               write(iulog,*) 'gid ipe r8_win%offset_r qsize r8_win%size = ', gid,      &
                         ipe, r8_win%offset_r, qsize, r8_win%size
               call exit(1)
            endif
            tmpsize = ceiling(real(qsize,r8)/real(pkgs_per_pro,r8))
!$omp parallel do private(p,mysize,ierr)
            do p=1,MIN(pkgs_per_pro,qsize)
              mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
# if defined (USE_MPI2)
              call MPI_PUT(ga_r8_s(r8_win%offset_s+(p-1)*tmpsize+1), mysize, mp_r8,     &
                           r8_win%dest, r8_win%offset_r+(p-1)*tmpsize, mysize, mp_r8, &
                           r8_win%id, ierr)
# else
              call SHMEM_PUT64(ga_r8_r(r8_win%offset_r+(p-1)*tmpsize+1),    &
                           ga_r8_s(r8_win%offset_s+(p-1)*tmpsize+1),        &
                           mysize, r8_win%dest)
# endif
            enddo
#else
! nonblocking send for MPI-1
            send_tag = gid
            nsend = nsend + 1
            call MPI_ISEND(ga_r8_s(r8_win%offset_s+1), qsize, mp_r8, r8_win%dest, &
                           send_tag, comm, sqest(nsend), ierr)
#endif
         endif
      enddo

    endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

      end subroutine mp_sendirr
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recvirr --- Read r8 contiguous parcels to global array
!
! !INTERFACE:
      subroutine mp_recvirr ( comm, qout, recv_bl )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
! !INPUT PARAMETERS:
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! Global Array Window
! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: qout(*)               ! local data segment
!
! !DESCRIPTION:
!
!     Complete transfer of a generalized region initiated by {\tt mp\_sendirr}.
!
! !REVISION HISTORY:
!    02.08.15   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Now using packed arrays for MPI2
!    04.02.24   Mirin       Various mpi2 options
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ipe, blocksize, offset_r, mod_method, mpi2_flag
      integer Ierr
      integer InStats(numpro*MPI_STATUS_SIZE)
      integer OutStats(numpro*MPI_STATUS_SIZE)
      integer i, j, num_r
      integer :: offset_v (Max_Nparcels)

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

!     num_r = 0 if this processes is not part of the receiving decomposition
      num_r = size(recv_bl)
      if (recv_bl(1)%Nparcels == -1) then
         num_r = 0
      endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      mpi2_flag = 0
#if defined (USE_MPI2)
      mpi2_flag = 1
#endif
      mod_method = recv_bl(1)%method
#if defined (USE_SHMEM)
      if (mod_method .ne. 0) then
         write(iulog,*) 'mp_recvirr with shmem - invalid mod_method = ', mod_method, ' - exiting'
         call exit(1)
      endif
#endif

    call Win_Close(comm, r8_win)
    if (mod_method .gt. 0 .and. mpi2_flag .eq. 0) then

! mpi-1 derived types
      EndTrf = MOD(EndTrf,MAX_TRF) + 1

! I cannot tell if this is a sending process. Either way, unused handles
!   are set to MPI_REQUEST_NULL, so this call should work (AAM).

      CALL MPI_WAITALL( numcomm, InHandle(:,EndTrf), InStats, Ierr )
      if (num_r .gt. 0) then
         CALL MPI_WAITALL( numcomm, OutHandle(:,EndTrf), OutStats, Ierr )
      endif

    else

! temporary contiguous buffer / global window
! scatter data from global receive buffer to final destination
      offset_r = 0
      do ipe=1, num_r
         r8_win%size_r = recv_bl(ipe)%Tot_Size
         if (r8_win%size_r .ne. 0) then
            r8_win%offset_r = offset_r
            offset_r = offset_r + r8_win%size_r
            if (offset_r .gt. r8_win%size) then
              write(iulog,*) "Fatal mp_recvirr: receive window out of space - exiting"
              write(iulog,*) 'gid ipe offset_r r8_win%size = ', gid,  &
                        ipe, offset_r, r8_win%size
              call exit(1)
            endif

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
            nread = nread + 1
            call MPI_WAIT(rqest(nread), Status, ierr)
#endif
            offset_v(1) = r8_win%offset_r
            do j = 2, recv_bl(ipe)%Nparcels
               offset_v(j) = offset_v(j-1) + recv_bl(ipe)%blocksizes(j-1)
            enddo

!dir$ concurrent
!dir$ preferstream
            do j = 1, recv_bl(ipe)%Nparcels
!dir$ concurrent
               do i = 1, recv_bl(ipe)%blocksizes(j)
                  qout(recv_bl(ipe)%displacements(j)+i) = ga_r8_r(offset_v(j)+i)
               enddo
            enddo

         endif
      enddo
    endif

! Contains MPI_WAITALL(nsend,sqest,...), which should really only be there for
!   contiguous buffer (mod_method=0) case; however, for non-contiguous
!   buffer case (mod_method=1), nsend is 0, so it shouldn't matter.

    call Win_Finalize(comm, r8_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recvirr
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_sendirr_r4 --- Write r4 contiguous parcels to global array
!
! !INTERFACE:
      subroutine mp_sendirr_r4 ( comm, q, send_bl, recv_bl, qout )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      real(r4), intent(in) :: q(*)                     ! local data segment

! !OUTPUT PARAMETERS:
      real(r4), intent(out) :: qout(*)                 ! local output segment
!
! !DESCRIPTION:
!     Send a number of contiguous parcels destined to an arbitrary set 
!     of PEs.  This is basically the non-blocking start of an
!     all-to-all communcation primitive.  This fundamental
!     routine forms a basis for higher level primitives.
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Use partneroffset
!    03.06.24   Sawyer      Integrated Use_Mpi_Types; added qout
!    04.02.24   Mirin       Various mpi2 options
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer ipe, qsize, offset, blocksize, nparcels, offset_s, offset_r, ierr, mod_method, mpi2_flag
      integer p, mysize, nthpc, minsize, nthrd, pn, pt, tmpsize
#if defined (USE_MPI2)
      integer (kind=MPI_ADDRESS_KIND) target_disp
#endif
      integer i, j, send_tag, recv_tag
      integer :: offset_v (Max_Nparcels)

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: mp_sendirr_r4 - r4 windows disabled - exiting'
        call exit(1)
#endif

!
!     initialize window

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      mpi2_flag = 0
#if defined (USE_MPI2)
      mpi2_flag = 1
#endif

      mod_method = recv_bl(1)%method
#if defined (USE_SHMEM)
      if (mod_method .ne. 0) then
         write(iulog,*) 'mp_sendirr_r4 with shmem - invalid mod_method = ', mod_method, ' - exiting'
         call exit(1)
      endif
#endif

      call Win_Open(comm, r4_win)
    if (mod_method .gt. 0) then

#if defined (USE_MPI2)
     if (mod_method .eq. 1) then
!  directly put contiguous segments into global target window
      do ipe = 1, numcomm
         qsize = send_bl(ipe)%Tot_Size
         if (qsize .ne. 0) then
!  thread over parcels, but take into consideration number of
!    threads versus number of parcels and minimum size of parcel
            nparcels = send_bl(ipe)%Nparcels
            minsize = send_bl(ipe)%blocksizes(1)
            offset_v(1) = 0
            do p = 2, nparcels
                minsize = min(minsize, send_bl(ipe)%blocksizes(p))
                offset_v(p) = offset_v(p-1) + send_bl(ipe)%blocksizes(p-1)
            enddo
            r4_win%dest = ipe-1
            nthpc = ceiling(real(pkgs_per_pro,r8)/real(nparcels,r8))
            nthrd = min(nthpc,minsize)      !  number of threads per block
            r4_win%offset_r = send_bl(ipe)%partneroffset
!$omp parallel do private(pn,p,pt,tmpsize,mysize,offset_s,ierr)
            do pn = 1, nparcels*nthrd
               p = (pn-1)/nthrd + 1        ! block number
               pt = mod(pn-1,nthrd) + 1    ! thread number
               tmpsize = ceiling(real(send_bl(ipe)%blocksizes(p),r8)/real(nthrd,r8))
               mysize = min(tmpsize, max(send_bl(ipe)%blocksizes(p)-tmpsize*(pt-1),0))
               offset_s = send_bl(ipe)%displacements(p)
               call MPI_PUT(q(offset_s+(pt-1)*tmpsize+1), mysize,       &
                            mp_r4, r4_win%dest,                         &
                            r4_win%offset_r+offset_v(p)+(pt-1)*tmpsize,   &
                            mysize, mp_r4, r4_win%id, ierr)
            enddo
         endif
      enddo
     elseif (mod_method .eq. 2) then
! directly put derived types into global target window
!   thread over targets
!$omp parallel do private(ipe,mysize,ierr,target_disp)
        do ipe = 1, numcomm
           if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
               mysize = send_bl(ipe)%Tot_Size
               target_disp = send_bl(ipe)%partneroffset
               call MPI_PUT(q, 1, send_bl(ipe)%type, ipe-1, target_disp,   &
                            mysize, mp_r4, r4_win%id, ierr)
           endif
        enddo
     endif
#elif defined(USE_SHMEM)
#else
!
! mpi-1 with derived types
! Increment the ongoing transfer number
      BegTrf = MOD(BegTrf,MAX_TRF) + 1
!
! MPI: Irecv over all processes
!
      do ipe=1, numcomm

!
! Receive the buffers with MPI_Irecv. Non-blocking
!
        OutHandle(ipe,BegTrf) = MPI_REQUEST_NULL
        if ( recv_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_irecv( qout, 1, recv_bl(ipe)%type, ipe-1, ipe-1,     &
                          comm, OutHandle(ipe,BegTrf), ierr )
        endif
      enddo

!
! MPI: Isend over all processes
!
      do ipe=1, numcomm

!
! Send the individual buffers with non-blocking sends
!
        InHandle(ipe,BegTrf) = MPI_REQUEST_NULL
        if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_isend( q, 1, send_bl(ipe)%type, ipe-1, gid,        &
                          comm, InHandle(ipe,BegTrf), ierr )
        endif
      enddo
#endif
    else

! temporary contiguous buffers
! for MPI-1, issue call to receive data in global receive buffer
      offset_s = 0
      offset_r = 0

      do ipe=1, numcomm
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
         r4_win%size_r = recv_bl(ipe)%Tot_Size
         if (r4_win%size_r .ne. 0) then
            r4_win%offset_r = offset_r
            offset_r = offset_r + r4_win%size_r
            r4_win%src = ipe-1
            if (r4_win%size >= r4_win%offset_r + r4_win%size_r) then
              recv_tag = r4_win%src
              qsize    = r4_win%size_r
              nrecv    = nrecv + 1
              call MPI_IRECV(ga_r4_r(r4_win%offset_r+1), qsize, mp_r4, r4_win%src, &
                             recv_tag, comm, rqest(nrecv), ierror)
            else
              write(iulog,*) "Fatal mp_sendirr_r4: receive window out of space - exiting"
              write(iulog,*) 'gid ipe r4_win%size r4_win%offset_r r4_win%size_r = ', gid,  &
                        ipe, r4_win%size, r4_win%offset_r, r4_win%size_r
              call exit(1)
            endif
         endif
#endif
! gather data into global send buffer
         qsize = send_bl(ipe)%Tot_Size
         if (qsize .ne. 0) then
            r4_win%dest = ipe-1
            r4_win%offset_s = offset_s
            offset_s = offset_s + qsize
            if (offset_s .gt. r4_win%size) then
              write(iulog,*) "Fatal mp_sendirr_r4: send window out of space - exiting"
              write(iulog,*) 'gid ipe offset_s r4_win%size = ', gid,  &
                        ipe, offset_s, r4_win%size
              call exit(1)
            endif

            offset_v(1) = r4_win%offset_s
            do j = 2, send_bl(ipe)%nparcels
               offset_v(j) = offset_v(j-1) + send_bl(ipe)%blocksizes(j-1)
            enddo

!dir$ concurrent
!dir$ preferstream
            do j = 1, send_bl(ipe)%nparcels
!dir$ concurrent
               do i = 1, send_bl(ipe)%blocksizes(j)
                  ga_r4_s(offset_v(j)+i) = q(send_bl(ipe)%displacements(j)+i)
               enddo
            enddo

#if defined (USE_MPI2) || defined (USE_SHMEM)
! one-sided puts for MPI-2 or SHMEM
            r4_win%offset_r = send_bl(ipe)%partneroffset
            if ( r4_win%offset_r+qsize > r4_win%size ) then
               write(iulog,*) "Fatal mp_sendirr_r4: target window out of space - exiting"
               write(iulog,*) 'gid ipe r4_win%offset_r qsize r4_win%size = ', gid,      &
                         ipe, r4_win%offset_r, qsize, r4_win%size
               call exit(1)
            endif
            tmpsize = ceiling(real(qsize,r8)/real(pkgs_per_pro,r8))
!$omp parallel do private(p,mysize,ierr)
            do p=1,MIN(pkgs_per_pro,qsize)
              mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
# if defined (USE_MPI2)
              call MPI_PUT(ga_r4_s(r4_win%offset_s+(p-1)*tmpsize+1), mysize, mp_r4,     &
                           r4_win%dest, r4_win%offset_r+(p-1)*tmpsize, mysize, mp_r4, &
                           r4_win%id, ierr)
# else
              call SHMEM_PUT32(ga_r4_r(r4_win%offset_r+(p-1)*tmpsize+1),    &
                           ga_r4_s(r4_win%offset_s+(p-1)*tmpsize+1),        &
                           mysize, r4_win%dest)
# endif
            enddo
#else
! nonblocking send for MPI-1
            send_tag = gid
            nsend = nsend + 1
            call MPI_ISEND(ga_r4_s(r4_win%offset_s+1), qsize, mp_r4, r4_win%dest, &
                           send_tag, comm, sqest(nsend), ierr)
#endif
         endif
      enddo

    endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

      end subroutine mp_sendirr_r4
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recvirr_r4 --- Read r4 contiguous parcels to global array
!
! !INTERFACE:
      subroutine mp_recvirr_r4 ( comm, qout, recv_bl )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! Global Array Window
! !INPUT/OUTPUT PARAMETERS:
      real(r4), intent(inout) :: qout(*)               ! local data segment
!
! !DESCRIPTION:
!
!     Complete transfer of a generalized region initiated by {\tt mp\_sendirr\_r4}.
!
! !REVISION HISTORY:
!    02.08.15   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Now using packed arrays for MPI2
!    04.02.24   Mirin       Various mpi2 options
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ipe, blocksize, offset_r, mod_method, mpi2_flag
      integer Ierr
      integer InStats(numpro*MPI_STATUS_SIZE)
      integer OutStats(numpro*MPI_STATUS_SIZE)
      integer i, j
      integer :: offset_v (Max_Nparcels)

#if defined ( NOR4 )
        write(iulog,*) 'Mod_comm: mp_recvirr_r4 - r4 windows disabled - exiting'
        call exit(1)
#endif

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      mpi2_flag = 0
#if defined (USE_MPI2)
      mpi2_flag = 1
#endif
      mod_method = recv_bl(1)%method
#if defined (USE_SHMEM)
      if (mod_method .ne. 0) then
         write(iulog,*) 'mp_recvirr_r4 with shmem - invalid mod_method = ', mod_method, ' - exiting'
         call exit(1)
      endif
#endif

    call Win_Close(comm, r4_win)
    if (mod_method .gt. 0 .and. mpi2_flag .eq. 0) then

! mpi-1 derived types
      EndTrf = MOD(EndTrf,MAX_TRF) + 1
      CALL MPI_WAITALL( numcomm, InHandle(:,EndTrf), InStats, Ierr )
      CALL MPI_WAITALL( numcomm, OutHandle(:,EndTrf), OutStats, Ierr )

    else

! temporary contiguous buffer / global window
! scatter data from global receive buffer to final destination
      offset_r = 0
      do ipe=1, numcomm
         r4_win%size_r = recv_bl(ipe)%Tot_Size
         if (r4_win%size_r .ne. 0) then
            r4_win%offset_r = offset_r
            offset_r = offset_r + r4_win%size_r
            if (offset_r .gt. r4_win%size) then
              write(iulog,*) "Fatal mp_recvirr_r4: receive window out of space - exiting"
              write(iulog,*) 'gid ipe offset_r r4_win%size = ', gid,  &
                        ipe, offset_r, r4_win%size
              call exit(1)
            endif

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
            nread = nread + 1
            call MPI_WAIT(rqest(nread), Status, ierr)
#endif
            offset_v(1) = r4_win%offset_r
            do j = 2, recv_bl(ipe)%Nparcels
               offset_v(j) = offset_v(j-1) + recv_bl(ipe)%blocksizes(j-1)
            enddo

!dir$ concurrent
!dir$ preferstream
            do j = 1, recv_bl(ipe)%Nparcels
!dir$ concurrent
               do i = 1, recv_bl(ipe)%blocksizes(j)
                  qout(recv_bl(ipe)%displacements(j)+i) = ga_r4_r(offset_v(j)+i)
               enddo
            enddo

         endif
      enddo
    endif
    call Win_Finalize(comm, r4_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recvirr_r4
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_sendirr_i4 --- Write i4 contiguous parcels to global array
!
! !INTERFACE:
      subroutine mp_sendirr_i4 ( comm, q, send_bl, recv_bl, qout )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      integer(i4), intent(in) :: q(*)                     ! local data segment

! !OUTPUT PARAMETERS:
      integer(i4), intent(out) :: qout(*)                 ! local output segment
!
! !DESCRIPTION:
!     Send a number of contiguous parcels destined to an arbitrary set 
!     of PEs.  This is basically the non-blocking start of an
!     all-to-all communcation primitive.  This fundamental
!     routine forms a basis for higher level primitives.
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Use partneroffset
!    03.06.24   Sawyer      Integrated Use_Mpi_Types; added qout
!    04.02.24   Mirin       Various mpi2 options
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer ipe, qsize, offset, blocksize, nparcels, offset_s, offset_r, ierr, mod_method, mpi2_flag
      integer p, mysize, nthpc, minsize, nthrd, pn, pt, tmpsize
#if defined (USE_MPI2)
      integer (kind=MPI_ADDRESS_KIND) target_disp
#endif
      integer i, j, send_tag, recv_tag
      integer :: offset_v (Max_Nparcels)

!
!     initialize window

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      mpi2_flag = 0
#if defined (USE_MPI2)
      mpi2_flag = 1
#endif

      mod_method = recv_bl(1)%method
#if defined (USE_SHMEM)
      if (mod_method .ne. 0) then
         write(iulog,*) 'mp_sendirr_i4 with shmem - invalid mod_method = ', mod_method, ' - exiting'
         call exit(1)
      endif
#endif

      call Win_Open(comm, i4_win)
    if (mod_method .gt. 0) then

#if defined (USE_MPI2)
     if (mod_method .eq. 1) then
!  directly put contiguous segments into global target window
      do ipe = 1, numcomm
         qsize = send_bl(ipe)%Tot_Size
         if (qsize .ne. 0) then
!  thread over parcels, but take into consideration number of
!    threads versus number of parcels and minimum size of parcel
            nparcels = send_bl(ipe)%Nparcels
            minsize = send_bl(ipe)%blocksizes(1)
            offset_v(1) = 0
            do p = 2, nparcels
                minsize = min(minsize, send_bl(ipe)%blocksizes(p))
                offset_v(p) = offset_v(p-1) + send_bl(ipe)%blocksizes(p-1)
            enddo
            i4_win%dest = ipe-1
            nthpc = ceiling(real(pkgs_per_pro,r8)/real(nparcels,r8))
            nthrd = min(nthpc,minsize)      !  number of threads per block
            i4_win%offset_r = send_bl(ipe)%partneroffset
!$omp parallel do private(pn,p,pt,tmpsize,mysize,offset_s,ierr)
            do pn = 1, nparcels*nthrd
               p = (pn-1)/nthrd + 1        ! block number
               pt = mod(pn-1,nthrd) + 1    ! thread number
               tmpsize = ceiling(real(send_bl(ipe)%blocksizes(p),r8)/real(nthrd,r8))
               mysize = min(tmpsize, max(send_bl(ipe)%blocksizes(p)-tmpsize*(pt-1),0))
               offset_s = send_bl(ipe)%displacements(p)
               call MPI_PUT(q(offset_s+(pt-1)*tmpsize+1), mysize,       &
                            mp_i4, i4_win%dest,                         &
                            i4_win%offset_r+offset_v(p)+(pt-1)*tmpsize,   &
                            mysize, mp_i4, i4_win%id, ierr)
            enddo
         endif
      enddo
     elseif (mod_method .eq. 2) then
! directly put derived types into global target window
!   thread over targets
!$omp parallel do private(ipe,mysize,ierr,target_disp)
        do ipe = 1, numcomm
           if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
               mysize = send_bl(ipe)%Tot_Size
               target_disp = send_bl(ipe)%partneroffset
               call MPI_PUT(q, 1, send_bl(ipe)%type, ipe-1, target_disp,   &
                            mysize, mp_i4, i4_win%id, ierr)
           endif
        enddo
     endif
#elif defined(USE_SHMEM)
#else
!
! mpi-1 with derived types
! Increment the ongoing transfer number
      BegTrf = MOD(BegTrf,MAX_TRF) + 1
!
! MPI: Irecv over all processes
!
      do ipe=1, numcomm

!
! Receive the buffers with MPI_Irecv. Non-blocking
!
        OutHandle(ipe,BegTrf) = MPI_REQUEST_NULL
        if ( recv_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_irecv( qout, 1, recv_bl(ipe)%type, ipe-1, ipe-1,     &
                          comm, OutHandle(ipe,BegTrf), ierr )
        endif
      enddo

!
! MPI: Isend over all processes
!
      do ipe=1, numcomm

!
! Send the individual buffers with non-blocking sends
!
        InHandle(ipe,BegTrf) = MPI_REQUEST_NULL
        if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_isend( q, 1, send_bl(ipe)%type, ipe-1, gid,        &
                          comm, InHandle(ipe,BegTrf), ierr )
        endif
      enddo
#endif
    else

! temporary contiguous buffers
! for MPI-1, issue call to receive data in global receive buffer
      offset_s = 0
      offset_r = 0

      do ipe=1, numcomm
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
         i4_win%size_r = recv_bl(ipe)%Tot_Size
         if (i4_win%size_r .ne. 0) then
            i4_win%offset_r = offset_r
            offset_r = offset_r + i4_win%size_r
            i4_win%src = ipe-1
            if (i4_win%size >= i4_win%offset_r + i4_win%size_r) then
              recv_tag = i4_win%src
              qsize    = i4_win%size_r
              nrecv    = nrecv + 1
              call MPI_IRECV(ga_i4_r(i4_win%offset_r+1), qsize, mp_i4, i4_win%src, &
                             recv_tag, comm, rqest(nrecv), ierror)
            else
              write(iulog,*) "Fatal mp_sendirr_i4: receive window out of space - exiting"
              write(iulog,*) 'gid ipe i4_win%size i4_win%offset_r i4_win%size_r = ', gid,  &
                        ipe, i4_win%size, i4_win%offset_r, i4_win%size_r
              call exit(1)
            endif
         endif
#endif
! gather data into global send buffer
         qsize = send_bl(ipe)%Tot_Size
         if (qsize .ne. 0) then
            i4_win%dest = ipe-1
            i4_win%offset_s = offset_s
            offset_s = offset_s + qsize
            if (offset_s .gt. i4_win%size) then
              write(iulog,*) "Fatal mp_sendirr_i4: send window out of space - exiting"
              write(iulog,*) 'gid ipe offset_s i4_win%size = ', gid,  &
                        ipe, offset_s, i4_win%size
              call exit(1)
            endif

            offset_v(1) = i4_win%offset_s
            do j = 2, send_bl(ipe)%nparcels
               offset_v(j) = offset_v(j-1) + send_bl(ipe)%blocksizes(j-1)
            enddo

!dir$ concurrent
!dir$ preferstream
            do j = 1, send_bl(ipe)%nparcels
!dir$ concurrent
               do i = 1, send_bl(ipe)%blocksizes(j)
                  ga_i4_s(offset_v(j)+i) = q(send_bl(ipe)%displacements(j)+i)
               enddo
            enddo

#if defined (USE_MPI2) || defined (USE_SHMEM)
! one-sided puts for MPI-2 or SHMEM
            i4_win%offset_r = send_bl(ipe)%partneroffset
            if ( i4_win%offset_r+qsize > i4_win%size ) then
               write(iulog,*) "Fatal mp_sendirr_i4: target window out of space - exiting"
               write(iulog,*) 'gid ipe i4_win%offset_r qsize i4_win%size = ', gid,      &
                         ipe, i4_win%offset_r, qsize, i4_win%size
               call exit(1)
            endif
            tmpsize = ceiling(real(qsize,r8)/real(pkgs_per_pro,r8))
!$omp parallel do private(p,mysize,ierr)
            do p=1,MIN(pkgs_per_pro,qsize)
              mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
# if defined (USE_MPI2)
              call MPI_PUT(ga_i4_s(i4_win%offset_s+(p-1)*tmpsize+1), mysize, mp_i4,     &
                           i4_win%dest, i4_win%offset_r+(p-1)*tmpsize, mysize, mp_i4, &
                           i4_win%id, ierr)
# else
              call SHMEM_PUT32(ga_i4_r(i4_win%offset_r+(p-1)*tmpsize+1),    &
                           ga_i4_s(i4_win%offset_s+(p-1)*tmpsize+1),        &
                           mysize, i4_win%dest)
# endif
            enddo
#else
! nonblocking send for MPI-1
            send_tag = gid
            nsend = nsend + 1
            call MPI_ISEND(ga_i4_s(i4_win%offset_s+1), qsize, mp_i4, i4_win%dest, &
                           send_tag, comm, sqest(nsend), ierr)
#endif
         endif
      enddo

    endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

      end subroutine mp_sendirr_i4
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recvirr_i4 --- Read i4 contiguous parcels to global array
!
! !INTERFACE:
      subroutine mp_recvirr_i4 ( comm, qout, recv_bl )
 
! !INPUT PARAMETERS:
      integer, intent(in)  :: comm      !  communicator
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! Global Array Window
! !INPUT/OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: qout(*)               ! local data segment
!
! !DESCRIPTION:
!
!     Complete transfer of a generalized region initiated by {\tt mp\_sendirr\_i4}
!
! !REVISION HISTORY:
!    02.08.15   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Now using packed arrays for MPI2
!    04.02.24   Mirin       Various mpi2 options
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ipe, blocksize, offset_r, mod_method, mpi2_flag
      integer Ierr
      integer InStats(numpro*MPI_STATUS_SIZE)
      integer OutStats(numpro*MPI_STATUS_SIZE)
      integer i, j
      integer :: offset_v (Max_Nparcels)

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call MPI_COMM_SIZE (comm, numcomm, ierror)

      mpi2_flag = 0
#if defined (USE_MPI2)
      mpi2_flag = 1
#endif
      mod_method = recv_bl(1)%method
#if defined (USE_SHMEM)
      if (mod_method .ne. 0) then
         write(iulog,*) 'mp_recvirr_i4 with shmem - invalid mod_method = ', mod_method, ' - exiting'
         call exit(1)
      endif
#endif

    call Win_Close(comm, i4_win)
    if (mod_method .gt. 0 .and. mpi2_flag .eq. 0) then

! mpi-1 derived types
      EndTrf = MOD(EndTrf,MAX_TRF) + 1
      CALL MPI_WAITALL( numcomm, InHandle(:,EndTrf), InStats, Ierr )
      CALL MPI_WAITALL( numcomm, OutHandle(:,EndTrf), OutStats, Ierr )

    else

! temporary contiguous buffer / global window
! scatter data from global receive buffer to final destination
      offset_r = 0
      do ipe=1, numcomm
         i4_win%size_r = recv_bl(ipe)%Tot_Size
         if (i4_win%size_r .ne. 0) then
            i4_win%offset_r = offset_r
            offset_r = offset_r + i4_win%size_r
            if (offset_r .gt. i4_win%size) then
              write(iulog,*) "Fatal mp_recvirr_i4: receive window out of space - exiting"
              write(iulog,*) 'gid ipe offset_r i4_win%size = ', gid,  &
                        ipe, offset_r, i4_win%size
              call exit(1)
            endif

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
            nread = nread + 1
            call MPI_WAIT(rqest(nread), Status, ierr)
#endif
            offset_v(1) = i4_win%offset_r
            do j = 2, recv_bl(ipe)%Nparcels
               offset_v(j) = offset_v(j-1) + recv_bl(ipe)%blocksizes(j-1)
            enddo

!dir$ concurrent
!dir$ preferstream
            do j = 1, recv_bl(ipe)%Nparcels
!dir$ concurrent
               do i = 1, recv_bl(ipe)%blocksizes(j)
                  qout(recv_bl(ipe)%displacements(j)+i) = ga_i4_r(offset_v(j)+i)
               enddo
            enddo

         endif
      enddo
    endif
    call Win_Finalize(comm, i4_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recvirr_i4
!------------------------------------------------------------------------------

#endif
      end module mod_comm

