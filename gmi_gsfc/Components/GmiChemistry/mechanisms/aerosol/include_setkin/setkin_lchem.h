!=======================================================================
!
! $Id: setkin_lchem.h,v 1.2 2011-08-08 17:46:52 mrdamon Exp $
!
! FILE
!   setkin_lchem.h - character labels for species and reactions
!             (setkin_lchem.h)
!   10 JUN 98 - PSC
!
! DESCRIPTION
!   Include file that provides ascii strings identifying reactions
!   and species
!
!  Chemistry input file generated: 2/23/00 10:12AM
!  Reaction dictionary:            3/29/99 2:55PM
!  Setkin files generated:         Wed Feb 23 10:30:41 2000
!
!   "setkin_lchem.h" was modified for sulfchem by Xiaohong Liu
!
!=======================================================================

      integer i

      character*16, save :: lchemvar(NSP),  ldynvar(NDYN)
      character*80, save :: lqkchem(NUM_K), lqjchem(NUM_J)

!
!.... All species labels
!
      data lchemvar(1) /"H2O2"/
      data lchemvar(2) /"fSO2"/
      data lchemvar(3) /"nSO2"/
      data lchemvar(4) /"nDMS"/
      data lchemvar(5) /"fSO4a"/
      data lchemvar(6) /"fSO4n1"/
      data lchemvar(7) /"fSO4n2"/
      data lchemvar(8) /"fSO4n3"/
      data lchemvar(9) /"nSO4a"/
      data lchemvar(10) /"nSO4n1"/
      data lchemvar(11) /"nSO4n2"/
      data lchemvar(12) /"nSO4n3"/
      data lchemvar(13) /"nOC"/
      data lchemvar(14) /"fOC"/
      data lchemvar(15) /"fBC"/
      data lchemvar(16) /"bOC"/
      data lchemvar(17) /"bBC"/
      data lchemvar(18) /"dust1"/
      data lchemvar(19) /"dust2"/
      data lchemvar(20) /"dust3"/
      data lchemvar(21) /"dust4"/
      data lchemvar(22) /"sslt1"/
      data lchemvar(23) /"sslt2"/
      data lchemvar(24) /"sslt3"/
      data lchemvar(25) /"sslt4"/
      data lchemvar(26) /"O3"/
      data lchemvar(27) /"OH"/
      data lchemvar(28) /"HO2"/
      data lchemvar(29) /"NO3"/
      data lchemvar(30) /"ad"/
!
!.... Dynamic (transported) species labels
!
      data ldynvar(1) /"H2O2"/
      data ldynvar(2) /"fSO2"/
      data ldynvar(3) /"nSO2"/
      data ldynvar(4) /"nDMS"/
      data ldynvar(5) /"fSO4a"/
      data ldynvar(6) /"fSO4n1"/
      data ldynvar(7) /"fSO4n2"/
      data ldynvar(8) /"fSO4n3"/
      data ldynvar(9) /"nSO4a"/
      data ldynvar(10) /"nSO4n1"/
      data ldynvar(11) /"nSO4n2"/
      data ldynvar(12) /"nSO4n3"/
      data ldynvar(13) /"nOC"/
      data ldynvar(14) /"fOC"/
      data ldynvar(15) /"fBC"/
      data ldynvar(16) /"bOC"/
      data ldynvar(17) /"bBC"/
      data ldynvar(18) /"dust1"/
      data ldynvar(19) /"dust2"/
      data ldynvar(20) /"dust3"/
      data ldynvar(21) /"dust4"/
      data ldynvar(22) /"sslt1"/
      data ldynvar(23) /"sslt2"/
      data ldynvar(24) /"sslt3"/
      data ldynvar(25) /"sslt4"/
      data ldynvar(26) /"O3"/
      data ldynvar(27) /"OH"/
      data ldynvar(28) /"HO2"/
      data ldynvar(29) /"NO3"/

!
!.... Thermal reaction labels
!
      data (lqkchem(i), i=1,8) /  &
     & '???',  &
     & '???',  &
     & '???',  &
     & '???',  &
     & '???',  &
     & '???',  &
     & '???',  &
     & '???' /

      data (lqjchem(i), i=1,1) /  &
     & '???' /
