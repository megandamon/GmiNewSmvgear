!=======================================================================
!  Code Developer
!    Nicholas Meskhidze, GA TECH
!    nmeskhidze@eas.gatech.edu
!
!     DESCRIPTION
!
! *** include file "parametr.inc"
! *** this file contains the declarations of the global constants
!     and variables.
!
! *** written by athanasios nenes
!
!=======================================================================
!
!      implicit double precision (a-h,o-z)
!
      integer, parameter :: nsmx_par = 2000
      integer, parameter :: npgauss  = 10
      integer, parameter :: maxit_par= 100      ! max iterations for solution

      real*8 , parameter :: eps_par  = 1.0d-6   ! convergence criterion
      real*8 , parameter :: amw_par  = 18d-3    ! water molecular weight
      real*8 , parameter :: ama_par  = 29d-3    ! air molecular weight
      real*8 , parameter :: grav_par = 9.81d0   ! g constant
      real*8 , parameter :: rgas_par = 8.31d0   ! universal gas constant
      real*8 , parameter :: accom_par= 0.06     ! default accommodation coef
      real*8 , parameter :: pi_par   = 3.1415927d0             ! some constants
      real*8 , parameter :: zero_par = 0d0
      real*8 , parameter :: great_par= 1d30
      real*8 , parameter :: sq2pi_par= 2.5066282746d0

      real*8 , parameter :: denw_par = 1d3             ! water density
      real*8 , parameter :: dhv_par  = 2.25d6          ! water enthalpy of vaporization
      real*8 , parameter :: cpair_par= 1.0061d3        ! air cp
 !
      real*8  ::   xgs_par, wgs_par
      logical ::  crit2, ccnspst
      real*8  ::  wparcel, temp_par, pres_par, temp_AG, pres_AG

      common /inputs/ wparcel, temp_par, pres_par, temp_AG, pres_AG
!
      integer ::     nmd_par
      real*8  ::     sg_par,  tp_par, dpg_par,sig_par,  &
     &               vhf_par, ams_par,dens_par, deni_par,  &
     &               amfs_par, dbig_par,  &
     &               dlow_par

      real*8  ::     sg_AG,  tp_AG,dpg_AG, sig_AG, vhf_AG, ams_AG, &
     &               dens_AG,deni_AG, amfs_AG 

!
      common/ccnspc/ sg_par(nsmx_par), tp_par(nsmx_par),  &
     &               dpg_par(nsmx_par),sig_par(nsmx_par),  &
     &               vhf_par(nsmx_par), ams_par(nsmx_par),  &
     &               dens_par(nsmx_par), deni_par(nsmx_par),  &
     &               amfs_par(nsmx_par), dbig_par,  &
     &               dlow_par,nmd_par
!
      common/ccnscAG/  sg_AG(nsmx_par),  tp_AG(nsmx_par),  &
     &                dpg_AG(nsmx_par), sig_AG(nsmx_par),  &
     &                vhf_AG(nsmx_par), ams_AG(nsmx_par),  &
     &               dens_AG(nsmx_par),deni_AG(nsmx_par),  &
     &               amfs_AG(nsmx_par)
!
      real*8  ::      akoh_par, ssplt_par, alfa_par, bet1_par,  &
     &                bet2_par, crit2_par, ccnspst_par
!
      common/actvpr/  akoh_par, ssplt_par, alfa_par, bet1_par,  &
     &                bet2_par, crit2_par, ccnspst_par
!
      real*8  ::       akoh_AG, alfa_AG, bet1_AG, bet2_AG  
!
      common/actvprAG/ akoh_AG, alfa_AG, bet1_AG, bet2_AG  

!
      real*8  ::       &
     &                dhv_AG, aka_par,  &
     &                dv_par,  psat_par, dair_par, surt_par, dv_AG
!
      common/thermo/  &
     &                aka_par,  &
     &                dv_par,  psat_par, dair_par, surt_par, dhv_AG,dv_AG
!
      real*8  ::      wpdbg, pddbg, nadbg,smdbg
!
      common/pardbg/ wpdbg(npgauss), pddbg(npgauss),  &
     &                nadbg(npgauss), smdbg(npgauss)
!
!      integer ::      niter_par
!      common/slnpar/ niter_par
!
      common/gaussl/ xgs_par(npgauss), wgs_par(npgauss)
!


