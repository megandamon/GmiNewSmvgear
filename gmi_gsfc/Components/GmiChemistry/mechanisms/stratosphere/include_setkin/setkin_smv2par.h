!=======================================================================
!
! $Id: setkin_smv2par.h,v 1.2 2011-08-09 22:12:58 mrdamon Exp $
!
! FILE
!   setkin_smv2par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Chemistry input file:    GMIS2 4:19 PM 1/23/2003
!  Reaction dictionary:     JPL00_reactions.db
!  Setkin files generated:  Fri Jan 24 00:31:24 2003
!
!========1=========2=========3=========4=========5=========6=========7==

      integer  &
     &  SK_IGAS  &
     & ,SK_IPHOT  &
     & ,SK_ITHERM  &
     & ,SK_NACT

      parameter (SK_IGAS   =  55)
      parameter (SK_IPHOT  =  44)
      parameter (SK_ITHERM = 122)
      parameter (SK_NACT   =  52)

!                                  --^--

