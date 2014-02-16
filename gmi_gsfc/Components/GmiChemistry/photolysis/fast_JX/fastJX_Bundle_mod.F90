!------------------------------------------------------------------------------
! NASA GSFC - SSSO Code 610.3
!------------------------------------------------------------------------------

      module fastJX_Bundle_mod

      USE fastJX_Parameters_mod

      private
      PUBLIC  :: t_fastJXbundle


# include "gmi_AerDust_const.h"
# include "setkin_par.h"

      TYPE t_fastJXbundle
           !---------------------------------------------------------------
           ! fast_JX_jv_cmn.h
           INTEGER :: NB ! Number of levels in CTM plus one for above model top
           INTEGER :: NC ! Number of levels in the fundamental Fast-J grid
           integer :: ncat_acetone_trop,ncat_acetone_stra ! Added for GMI to {AOO, 8/04}
           integer :: ncat_met_vinyl_ketone               ! Added for GMI to {AOO, 8/04}
           integer :: ncat_met_ethyl_ketone               ! Added for GMI to {AOO, 8/04}
           integer :: ncat_methyl_glyoxal                 ! Added for GMI to {AOO, 8/04}
           INTEGER :: jind(JPPJ_MAX), jadsub(NC_MAX), nhz, hzind(NH)
           INTEGER :: NJVAL, NW1, NW2, MIEDX(MX), NAA, NLBATM, npdep, jpdep(NS)
           REAL*8  :: dtaumax, szamax, zj(JPNL_MAX,JPPJ_MAX), jfacta(JPPJ_MAX)
           REAL*8  :: dtausub, dsubdiv
           REAL*8  :: RAD, RFLECT, SZA, U0, TANHT, ZZHT
           REAL*8  :: TREF(51,18,12), OREF(51,18,12), BREF(51)
           REAL*8  :: TJ(NB_MAX), PJ(NB_MAX+1), DM(NB_MAX), DO3(NB_MAX),  &
                      DBC(NB_MAX), Z(NB_MAX), AER(MX,NB_MAX),  &
                      AMF(NB_MAX,NB_MAX)
           REAL*8  :: zpdep(NW,3)
           REAL*8  :: WBIN(NW+1), WL(NW), FL(NW), QO2(NW,3), QO3(NW,3),  &
                      Q1D(NW,3), QQQ(NW,2,NS-3), QRAYL(NW+1), TQQ(3,NS),  &
                      FFF(NW,jpnl_max), VALJ(NS), WAA(4,NP), QAA(4,NP),  &
                      PAA(8,4,NP), RAA(4,NP), SSA(4,NP), QBC(NW)
           REAL*8  :: FF1(NW,jpnl_max) ! Flux@(TOA)*transmittance
           REAL*8  :: FF2(NW,jpnl_max) ! Flux*(Transmittance^2)
           REAL*8  :: FF3(NW,jpnl_max)
           CHARACTER*20 :: TITLEA(NP)
           CHARACTER*78 :: TITLE0
           CHARACTER*7  :: TITLEJ(3,NS)
           CHARACTER*7  :: jlabel(JPPJ_MAX)
           CHARACTER*7  :: hzlab(NH)

           !----------------------------------------------------------------------------
           !fast_JX_cmn_h.h

           integer :: lpar
           integer :: jpnl, jppj
           real*8  :: xgrd(ipar)         !  Longitude (midpoint, radians)
           real*8  :: ygrd(jpar)         ! Latitude  (midpoint, radians)
           real*8  :: ydgrd(jpar)        ! Latitude  (midpoint, degrees)
           real*8  :: etaa(lpar_max+1)   ! Eta(a) value for level boundaries
           real*8  :: etab(lpar_max+1)   ! Eta(b) value for level boundaries
           real*8  :: tau                ! Time of Day (hours, GMT)
           integer :: month              ! Number of month (1-12)
           integer :: iday               ! Day of year
           logical :: do_model_clima     ! determines if climatology data come
                                         ! from the model.

           !----------------------------------------------------------------------------
           ! fast_JX_cmn_w.h

           real*8  :: P(ipar,jpar)                   ! Surface pressure
           real*8  :: T(ipar,jpar,lpar_max+1)        ! Temperature profile
           real*8  :: OD(ipar,jpar,lpar_max)         ! Optical Depth profile
           real*8  :: SA(ipar,jpar)                  ! Surface Albedo
           real*8  :: optaer(lpar_max,NSADaer*nrh_b) ! Opt. Depth profile for aerosols
           real*8  :: optdust(lpar_max,NSADdust)     ! Opt. Depth profile for dust
           integer :: lantyp(ipar,jpar)              !  Land type

           !----------------------------------------------------------------------------
           ! fast_JX_jv_mie.h

           REAL*8  :: ZREFL,ZFLUX,RADIUS,ZU0
           INTEGER :: ND,N,M,MFIT

           REAL*8  :: A(M__)
           REAL*8  :: C1(M__)
           REAL*8  :: H(M__)
           REAL*8  :: EMU(M__)
           REAL*8  :: V1(M__)
           REAL*8  :: WT(M__)
           REAL*8  :: B(M__,M__)
           REAL*8  :: AA(M__,M__)
           REAL*8  :: CC(M__,M__)
           REAL*8  :: S(M__,M__)
           REAL*8  :: W(M__,M__)
           REAL*8  :: U1(M__,M__)
           REAL*8  :: PM(M__,2*M__)
           REAL*8  :: PM0(2*M__)
           REAL*8  :: POMEGA(2*M__,N__)
           REAL*8  :: ZTAU(N__)
           REAL*8  :: FZ(N__)
           REAL*8  :: FJ(N__)
           REAL*8  :: DD(M__,M__,N__)
           REAL*8  :: RR(M__,N__)

           !----------------------------------------------------------------------------
           !fast_JX_cmn_t.h

           integer :: nslat               !  Latitude of current profile point
           integer :: nslon               !  Longitude of current profile point 

           !----------------------------------------------------------------------------

      END TYPE t_fastJXbundle

      end module fastJX_Bundle_mod
