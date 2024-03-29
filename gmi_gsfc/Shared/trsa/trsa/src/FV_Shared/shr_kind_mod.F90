!===============================================================================
! CVS: $Id: shr_kind_mod.F90,v 1.2 2011-08-09 22:13:00 mrdamon Exp $
! CVS: $Source: /cvsroot/gmi/gmi_gsfc/Shared/trsa/trsa/src/FV_Shared/shr_kind_mod.F90,v $
! CVS: $Name:  $
!===============================================================================

MODULE shr_kind_mod

   !----------------------------------------------------------------------------
   ! precision/kind constants add data public
   !----------------------------------------------------------------------------
   public
#ifdef NO_R16
   integer,parameter :: SHR_KIND_R16= selected_real_kind(12) ! 8 byte real
#else
   integer,parameter :: SHR_KIND_R16= selected_real_kind(24) ! 16 byte real
#endif
   integer,parameter :: SHR_KIND_R8 = selected_real_kind(12) ! 8 byte real
   integer,parameter :: SHR_KIND_R4 = selected_real_kind( 6) ! 4 byte real
   integer,parameter :: SHR_KIND_RN = kind(1.0)              ! native real
   integer,parameter :: SHR_KIND_I8 = selected_int_kind (13) ! 8 byte integer
   integer,parameter :: SHR_KIND_I4 = selected_int_kind ( 6) ! 4 byte integer
   integer,parameter :: SHR_KIND_IN = kind(1)                ! native integer
   integer,parameter :: SHR_KIND_CL = 256                    ! long char
   integer,parameter :: SHR_KIND_CS = 80                     ! short char

END MODULE shr_kind_mod
