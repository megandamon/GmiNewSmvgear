!=======================================================================
!
! DESCRIPTION
!   Include file that provides ascii strings identifying reactions
!   and species
!
!
!=======================================================================


      integer i

      character*16 lchemvar(NSP), ldynvar(NDYN)
      character*80 lqkchem(NUM_K),  lqjchem(NUM_J)
!
!.... All species labels
!
      data lchemvar( 1) /"Age"/
      data lchemvar( 2) /"e90"/
      data lchemvar( 3) /"tm25"/
      data lchemvar( 4) /"Rn-222"/
      data lchemvar( 5) /"Pb-210"/
      data lchemvar( 6) /"Pb-210S"/
      data lchemvar( 7) /"Be-7"/
      data lchemvar( 8) /"Be-10"/
      data lchemvar( 9) /"Be-7S"/
      data lchemvar(10) /"Be-10S"/
      data lchemvar(11) /"CH3I"/
      data lchemvar(12) /"fCO2"/
      data lchemvar(13) /"LINOZ"/
      data lchemvar(14) /"SYNOZ"/
      data lchemvar(15) /"SF6"/
      data lchemvar(16) /"CLOCK"/
      data lchemvar(17) /"Uniform"/
      data lchemvar(18) /"Strat_O3"/
!
!.... Dynamic (transported) species labels
!
      data ldynvar( 1) /"Age"/
      data ldynvar( 2) /"e90"/
      data ldynvar( 3) /"tm25"/
      data ldynvar( 4) /"Rn-222"/
      data ldynvar( 5) /"Pb-210"/
      data ldynvar( 6) /"Pb-210S"/
      data ldynvar( 7) /"Be-7"/
      data ldynvar( 8) /"Be-10"/
      data ldynvar( 9) /"Be-7S"/
      data ldynvar(10) /"Be-10S"/
      data ldynvar(11) /"CH3I"/
      data ldynvar(12) /"fCO2"/
      data ldynvar(13) /"LINOZ"/
      data ldynvar(14) /"SYNOZ"/
      data ldynvar(15) /"SF6"/
      data ldynvar(16) /"CLOCK"/
      data ldynvar(17) /"Uniform"/
      data ldynvar(18) /"Strat_O3"/

!                                  --^--

