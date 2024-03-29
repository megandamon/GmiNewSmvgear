!
!  Wind data header file required by CTM - provide just enough stuff here
!  to allow model subroutines to be run as a stand-alone code.
!
!                                                 Oliver (9 July 99)
!-----------------------------------------------------------------------
!
      real*8  P(ipar,jpar)         !  Surface pressure
!      real*8  T(ipar,jpar,lpar+1)  !  Temperature profile
      real*8  T(ipar,jpar,lpar_max+1)  !  Temperature profile
!
!      real*8  OD(ipar,jpar,lpar+1) !  Optical Depth profile
      real*8  OD(ipar,jpar,lpar_max) !  Optical Depth profile
      real*8  SA(ipar,jpar)        !  Surface Albedo
!
!!      real*8  optaer(ipar,jpar,lpar_max,NSADaer*nrh_b) ! Opt. Depth profile for aerosols
!!      real*8  optdust(ipar,jpar,lpar_max,NSADdust)     ! Opt. Depth profile for dust
      real*8  optaer(lpar_max,NSADaer*nrh_b) ! Opt. Depth profile for aerosols
      real*8  optdust(lpar_max,NSADdust)     ! Opt. Depth profile for dust
      integer lantyp(ipar,jpar)    !  Land type
!
      common/metdat/P,T,OD,SA, optaer, optdust, lantyp
!
!-----------------------------------------------------------------------
