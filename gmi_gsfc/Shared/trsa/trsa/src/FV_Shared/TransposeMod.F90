module TransposeMod
!BOP
!
! !MODULE: TransposeMod --- manages routes between XY and YZ decompositions
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4

!
   use decompmodule, only: decomptype
   use ghostmodule, only:  ghosttype
   use parutilitiesmodule, only: parpatterntype, REAL4, INT4

! !PUBLIC MEMBER FUNCTIONS:
   public TransposeInit, TransposeFinal

! !PUBLIC DATA MEMBERS:

   public T_TRANSPOSE

!
! T_TRANSPOSE contains information about the transpose from XY
! to/from YZ decomposition
!
 
  type T_TRANSPOSE

!
! PILGRIM communication information (was in spmd_dyn)
!
    type (parpatterntype) :: g2l_xy, l2g_xy, yz2xy, xy2yz
    type (parpatterntype) :: ng_xy2yz, ng_yz2xy, jp1_xy2yz, jp1_yz2xy

!
! END PILGRIM communication information
!

  end type T_TRANSPOSE

!
! !DESCRIPTION:
!
!      This module provides variables which are specific to the Lin-Rood
!      dynamical core.  Most of them were previously SAVE variables in 
!      different routines and were set with an "if (first)" statement.
!
!      \begin{tabular}{|l|l|} \hline \hline
!        TransposeInit  &  Initialize the transpose variables \\ \hline
!        TransposeFinal &  Deallocate the transpose variables \\ \hline 
!                                \hline
!      \end{tabular}
!
! !REVISION HISTORY:
!   06.07.26   Sawyer     Creation (from dynamics_vars.F90)
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: TransposeInit --- initialize the transpose variables
!
! !INTERFACE: 
   subroutine TransposeInit( im, jm, km, ng, method,                &
                             ifirstxy, ilastxy, jfirstxy, jlastxy,  &
                             jfirst, jlast, kfirst, klast,          &
                             nprxy_x, nprxy_y, npryz_y, npryz_z,    &
                             imxy, jmxy, jmyz, kmyz,                &
                             transposes )
! !USES:
      use decompmodule, only: decompcreate, decompfree
      use ghostmodule, only:  ghostcreate, ghostfree
      use parutilitiesmodule, only: commglobal, gid, parpatterncreate
      implicit none

! !INPUT PARAMETERS:
      integer, intent(in)   :: im, jm, km         ! Global dimensions
      integer, intent(in)   :: ng                 ! Ghost width
      integer, intent(in)   :: method             ! Communication method
      integer, intent(in)   :: ifirstxy, ilastxy  !  Interval
      integer, intent(in)   :: jfirstxy, jlastxy  !  Interval
      integer, intent(in)   :: jfirst, jlast      !  Interval
      integer, intent(in)   :: kfirst, klast      !  Interval
      integer, intent(in)   :: nprxy_x, nprxy_y, npryz_y, npryz_z
      integer, dimension(:), intent(in)   :: imxy
      integer, dimension(:), intent(in)   :: jmxy
      integer, dimension(:), intent(in)   :: jmyz
      integer, dimension(:), intent(in)   :: kmyz

! !INPUT/OUTPUT PARAMETERS:
      type(T_TRANSPOSE), intent(inout)   :: transposes ! Result

! !DESCRIPTION:
!
!   Initialize transpose variables
!
! !REVISION HISTORY:
!
!   06.07.26   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!

! !LOCAL VARIABLES:
      type(decomptype) :: decompxy, decompyz, global, global2d, decompxy2d
      type(ghosttype)  :: ghost_ng, ghost_jp1
      integer, allocatable :: xdist(:), ydist(:)
      integer  :: rc

!
! Create PILGRIM decompositions (see decompmodule)
!

      allocate( xdist(nprxy_x), ydist(nprxy_y) )
      xdist = 0
      ydist = 0
      xdist(1) = im
      ydist(1) = jm

      call decompcreate( nprxy_x, nprxy_y, xdist, ydist, global2d )
      call decompcreate( nprxy_x, nprxy_y, 1, xdist, ydist, (/km/), global )
      deallocate( xdist, ydist )

      call decompcreate( nprxy_x, nprxy_y, imxy, jmxy, decompxy2d )
      call decompcreate( nprxy_x, nprxy_y, 1, imxy, jmxy, (/km/), decompxy )
      call decompcreate( 1, npryz_y, npryz_z, (/im/), jmyz, kmyz, decompyz )

      call ghostcreate( decompyz, gid, im, 1, im, .true.,                   &
                        jm, jfirst-ng, jlast+ng, .false.,                   &
                        km, kfirst, klast, .false., ghost_ng )
 
      call ghostcreate( decompyz, gid, im, 1, im, .true.,                   &
                        jm, jfirst, jlast+1, .false.,                       &
                        km, kfirst, klast, .false., ghost_jp1 )

! Initialize transposes
!
      call parpatterncreate(commglobal, global2d, decompxy2d,               &
                            transposes%g2l_xy, mod_method=method)

      call parpatterncreate(commglobal, decompxy2d, global2d,               &
                            transposes%l2g_xy, mod_method=method)

      call parpatterncreate(commglobal, decompyz, decompxy,                 &
                            transposes%yz2xy, mod_method=method)

      call parpatterncreate(commglobal, decompxy, decompyz,                 &
                            transposes%xy2yz, mod_method=method)

      call parpatterncreate(commglobal, ghost_ng, decompxy,                 &
                            transposes%ng_yz2xy, mod_method=method)

      call parpatterncreate(commglobal, decompxy, ghost_ng,                 &
                            transposes%ng_xy2yz, mod_method=method)

      call parpatterncreate(commglobal, ghost_jp1, decompxy,                &
                            transposes%jp1_yz2xy, mod_method=method)

      call parpatterncreate(commglobal, decompxy, ghost_jp1,                &
                            transposes%jp1_xy2yz, mod_method=method)

      call decompfree( global2d )
      call decompfree( decompxy2d )
      call decompfree( global )
      call decompfree( decompxy )
      call decompfree( decompyz )
      call ghostfree(  ghost_ng )
      call ghostfree(  ghost_jp1)

      return

!EOC
   end subroutine TransposeInit
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!BOP
! !IROUTINE: TransposeFinal -- Finalize transpose variables
!
! !INTERFACE: 
   subroutine TransposeFinal(transposes)

! !USES:
      use parutilitiesmodule, only: commglobal, parpatternfree
      implicit none

! !INPUT/OUTPUT PARAMETERS:
      type(T_TRANSPOSE), intent(inout)  :: transposes 


! !DESCRIPTION:
!
! Cleans up (deallocate) transpose variables
!
! !REVISION HISTORY:
!
!   06.07.26   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! Temporary data structures

! Clean transposes
!
    call parpatternfree(commglobal, transposes%g2l_xy)
    call parpatternfree(commglobal, transposes%l2g_xy)
    call parpatternfree(commglobal, transposes%yz2xy)
    call parpatternfree(commglobal, transposes%xy2yz)
    call parpatternfree(commglobal, transposes%ng_yz2xy)
    call parpatternfree(commglobal, transposes%ng_xy2yz)
    call parpatternfree(commglobal, transposes%jp1_yz2xy)
    call parpatternfree(commglobal, transposes%jp1_xy2yz)

   return
!EOC
   end subroutine TransposeFinal
!-----------------------------------------------------------------------

end module TransposeMod

