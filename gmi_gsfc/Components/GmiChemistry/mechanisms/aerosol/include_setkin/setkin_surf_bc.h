!=======================================================================
!
! $Id: setkin_surf_bc.h,v 1.2 2011-08-08 17:46:52 mrdamon Exp $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_surf_bc.h
!
! DESCRIPTION
!   This include file contains information about surface boundary
!   conditions.
!
!  Chemistry input file:    2:20 PM 1/14/2003
!  Reaction dictionary:     NMHC reactions.db
!  Setkin files generated:  Tue Feb 11 19:35:44 2003
!
!=======================================================================

!     -----------------------
!     Parameter declarations.
!     -----------------------

!     -------------------------------------------------------
!     NUM_SBC : number of surface bounday conditions to reset
!     K_SBC   : max k to which species will be adjusted
!     -------------------------------------------------------

      integer, parameter :: NUM_SBC = 1

      integer, parameter :: K_SBC   = 1

!     ------------------
!     Integer variables.
!     ------------------

      integer, save :: sbc_map(NUM_SBC) =  &
     &  (/ 0 /)

!                                  --^--

