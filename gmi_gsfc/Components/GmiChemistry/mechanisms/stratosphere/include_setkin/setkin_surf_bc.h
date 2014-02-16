!=======================================================================
!
! $Id: setkin_surf_bc.h,v 1.2 2011-08-09 22:12:58 mrdamon Exp $
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
!  Chemistry input file:    GMIS2 4:19 PM 1/23/2003
!  Reaction dictionary:     JPL00_reactions.db
!  Setkin files generated:  Fri Jan 24 00:31:24 2003
!
!=======================================================================

!     -----------------------
!     Parameter declarations.
!     -----------------------

!     -------------------------------------------------------
!     NUM_SBC : number of surface bounday conditions to reset
!     K_SBC   : max k to which species will be adjusted
!     -------------------------------------------------------

      integer, parameter :: NUM_SBC = 2

      integer, parameter :: K_SBC   = 2

!     ------------------
!     Integer variables.
!     ------------------

      integer, save :: sbc_map(NUM_SBC) =  &
     &  (/ 17 ,22 /)

!                                  --^--

