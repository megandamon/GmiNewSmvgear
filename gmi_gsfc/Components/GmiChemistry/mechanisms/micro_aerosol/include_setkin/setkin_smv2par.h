!=======================================================================
!
! $Id: setkin_smv2par.h,v 1.2 2011-08-09 22:02:42 mrdamon Exp $
!
! FILE
!   setkin_smv2par.h
!   12 JUN 02 - PSC
!
! DESCRIPTION
!   This include file sets symbolic constants for the photochemical
!   mechanism.
!
!  Chemistry input file:    2:20 PM 1/14/2003
!  Reaction dictionary:     NMHC reactions.db
!  Setkin files generated:  Tue Feb 11 19:35:44 2003
!
!========1=========2=========3=========4=========5=========6=========7==

      integer  &
     &  SK_IGAS  &
     & ,SK_IPHOT  &
     & ,SK_ITHERM  &
     & ,SK_NACT

      parameter (SK_IGAS   =  91)
      parameter (SK_IPHOT  =  61)
      parameter (SK_ITHERM = 200)
      parameter (SK_NACT   =  88)

!                                  --^--

