!=======================================================================
!
! $Id: setkin_lchem.h,v 1.2 2011-08-09 22:02:42 mrdamon Exp $
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
      data lchemvar(5) /"SO4g"/
      data lchemvar(6) /"SO4m1"/
      data lchemvar(7) /"SO4n1"/
      data lchemvar(8) /"SO4m2"/
      data lchemvar(9) /"SO4n2"/
      data lchemvar(10) /"SO4nOC"/
      data lchemvar(11) /"SO4fOC"/
      data lchemvar(12) /"SO4fBC"/
      data lchemvar(13) /"SO4bOC"/
      data lchemvar(14) /"SO4bBC"/
      data lchemvar(15) /"SO4d1"/
      data lchemvar(16) /"SO4d2"/
      data lchemvar(17) /"SO4d3"/
      data lchemvar(18) /"SO4d4"/
      data lchemvar(19) /"SO4s1"/
      data lchemvar(20) /"SO4s2"/
      data lchemvar(21) /"SO4s3"/
      data lchemvar(22) /"SO4s4"/
      data lchemvar(23) /"nOC"/
      data lchemvar(24) /"fOC"/
      data lchemvar(25) /"fBC"/
      data lchemvar(26) /"bOC"/
      data lchemvar(27) /"bBC"/
      data lchemvar(28) /"dust1"/
      data lchemvar(29) /"dust2"/
      data lchemvar(30) /"dust3"/
      data lchemvar(31) /"dust4"/
      data lchemvar(32) /"sslt1"/
      data lchemvar(33) /"sslt2"/
      data lchemvar(34) /"sslt3"/
      data lchemvar(35) /"sslt4"/
      data lchemvar(36) /"O3"/
      data lchemvar(37) /"OH"/
      data lchemvar(38) /"HO2"/
      data lchemvar(39) /"NO3"/
      data lchemvar(40) /"ad"/
!
!.... Dynamic (transported) species labels
!
      data ldynvar(1) /"H2O2"/
      data ldynvar(2) /"fSO2"/
      data ldynvar(3) /"nSO2"/
      data ldynvar(4) /"nDMS"/
      data ldynvar(5) /"SO4g"/
      data ldynvar(6) /"SO4m1"/
      data ldynvar(7) /"SO4n1"/
      data ldynvar(8) /"SO4m2"/
      data ldynvar(9) /"SO4n2"/
      data ldynvar(10) /"SO4nOC"/
      data ldynvar(11) /"SO4fOC"/
      data ldynvar(12) /"SO4fBC"/
      data ldynvar(13) /"SO4bOC"/
      data ldynvar(14) /"SO4bBC"/
      data ldynvar(15) /"SO4d1"/
      data ldynvar(16) /"SO4d2"/
      data ldynvar(17) /"SO4d3"/
      data ldynvar(18) /"SO4d4"/
      data ldynvar(19) /"SO4s1"/
      data ldynvar(20) /"SO4s2"/
      data ldynvar(21) /"SO4s3"/
      data ldynvar(22) /"SO4s4"/
      data ldynvar(23) /"nOC"/
      data ldynvar(24) /"fOC"/
      data ldynvar(25) /"fBC"/
      data ldynvar(26) /"bOC"/
      data ldynvar(27) /"bBC"/
      data ldynvar(28) /"dust1"/
      data ldynvar(29) /"dust2"/
      data ldynvar(30) /"dust3"/
      data ldynvar(31) /"dust4"/
      data ldynvar(32) /"sslt1"/
      data ldynvar(33) /"sslt2"/
      data ldynvar(34) /"sslt3"/
      data ldynvar(35) /"sslt4"/
      data ldynvar(36) /"O3"/
      data ldynvar(37) /"OH"/
      data ldynvar(38) /"HO2"/
      data ldynvar(39) /"NO3"/
!
!.... Thermal reaction labels
!
      data (lqkchem(i), i=1,8) /  &
     & 'DMS + OH = SO2',  &
     & 'SO2 + OH = H2SO4',  &
     & 'HO2 + HO2 = H2O2 + O2',  &
     & 'H2O2 + OH = H2O + HO2',  &
     & 'SO2(aq) + H2O2(aq) = H2SO4(aq)',  &
     & 'SO2(aq) + O3(aq) = H2SO4(aq) + O2',  &
     & 'HO2 + HO2 + H2O = H2O2 + O2 + H2O',  &
     & 'DMS + NO3 = SO2 +' /

!
!.... Photolytic reaction labels
!
      data (lqjchem(i), i=1,1) /  &
     & 'H2O2 + hv = 2 OH' /

