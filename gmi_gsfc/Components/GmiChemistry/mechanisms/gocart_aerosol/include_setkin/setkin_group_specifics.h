!=======================================================================
!
! $Id: setkin_group_specifics.h,v 1.1 2006-07-03 02:54:25 kouatch Exp $
!
! CODE DEVELOPER
!   Peter Connell, LLNL
!   connell2@llnl.gov
!
! FILE
!   setkin_group_specifics.h
!
! DESCRIPTION
!   This include file contains information about grouping species for
!   transport purposes, that is specific to the mechanism. In general,
!   this file can not be generated automatically, but must be manually
!   produced.
!
!  Chemistry input file:    2:20 PM 1/14/2003
!  Reaction dictionary:     NMHC reactions.db
!  Setkin files generated:  Tue Feb 11 19:35:44 2003
!
!=======================================================================

      qq2(:,:,:)               = 0.0d0

      max_bry_adjust           = 0.0d0
      max_cly_adjust           = 0.0d0
      max_nox_adjust           = 0.0d0

      select case (ig)

!.... YOy

      end select

