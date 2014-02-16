!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: m_FastjCTM_FastjLinking
!
! !INTERFACE:
!
      module m_FastjCTM_FastjLinking
!
      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
      private
      public  :: PHOTOJ
!
! !DESCRIPTION:
!
! !AUTHOR:
!  Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
      CONTAINS

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: PHOTOJ
!
! !INTERFACE:
!
      subroutine PHOTOJ(UTIME,IDAY,ILNG,JLAT, SZA,ZPJ)
!
! !USES:
!
      use m_FastjMetVars      , only : ETAA, ETAB, DM, TJ, DO3, SA
      use m_FastjMetVars      , only : P, ZH, YGRD, XGRD
      use m_FastjMetVars      , only : NAER1, NAER2, NAER3, NCLDX
      use m_FastjMetVars      , only : DAER1, DAER2, DAER3, ODCLD
      use m_FastjJvaluesVars  , only : RAD, NAA, JTAUMX, ZZHT, ATAU, ATAU0
      use m_FastjJvaluesVars  , only : JIND, JFACTA, NJVAL, NW1, NW2, NRATJ
      use m_FastjJvaluesVars  , only : WL, FL, QRAYL, TQQ, QO2, QO3
      use m_FastjCTMparameters, only : L_, L1_, L2_, W_, X_
      use m_FastjCTMparameters, only : JVL_, JVN_, SZAMAX
      use m_FastjMIEparameters, only : N_
      use m_FastjCoreFastj    , only : EXTRAL
      use m_FastjCoreFastj    , only : SOLARZ
      use m_FastjCoreFastj    , only : OPMIE
      use m_FastjCoreFastj    , only : SPHERE
      use m_FastjCoreFastj    , only : FLINT
      use m_FastjCoreFastj    , only : JRATET
      use m_FastjCoreFastj    , only : JP_ATM
!
      implicit none
!
! !INPUT PARAMETERS:
      real*8, intent(in)  ::  UTIME
      integer,intent(in)  ::  IDAY,ILNG,JLAT
!
! !OUTPUT PARAMETERS:
      real*8 , intent(out) :: SZA
      real*8 , intent(out) :: ZPJ(JVL_,JVN_) !2-D array of J's indexed to CTM chemistry!
!
! !DESCRIPTION:
!  PHOTOJ is the gateway to fast-JX calculations:
!        only access to CTM 3-D GLOBAL arrays
!        sets up the 1-D column arrays for calculating J's
!   \begin{verbatim}
!     AVGF   Attenuation of beam at each level for each wavelength
!     FFF    Actinic flux at each desired level
!     XQO2   Absorption cross-section of O2
!     XQO3   Absorption cross-section of O3
!   \end{verbatim}
!
! !LOCAL VARIABLES:
!--------key amtospheric data needed to solve plane-parallel J---------
      real*8, dimension(5,L1_) :: AER
      integer,dimension(5,L1_) :: ADX
      real*8, dimension(L1_)   :: ABX, TTJ,DDJ,ZZJ,ZHL
      real*8, dimension(L1_+1) :: PPJ
      integer,dimension(L2_+1) :: JXTRA

      real*8                   :: RFLECT,U0,SOLF
      real*8                   :: AMF(L1_,L1_)
!------------key arrays AFTER solving for J's---------------------------
      real*8  FFF(W_,JVL_),VALJ(X_)
      real*8  VALJL(JVL_,NJVAL) !2-D array of J_s returned by JRATET

      integer  I,J,K,L,M,KM
      real*8   AVGF(L_),XQO3,XQO2    ,WAVE, TTT
!
! !AUTHOR:
! Michael Prather
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-----------------------------------------------------------------------
!BOC
      ZPJ(:,:) = 0.d0
      FFF(:,:) = 0.d0

!-----------------------------------------------------------------------
      call SOLARZ(UTIME,IDAY,YGRD(JLAT),XGRD(ILNG), SZA,U0,SOLF)
!-----------------------------------------------------------------------
      SOLF = 1.d0

!---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
!                        or         99.                  80 km
      if (SZA .gt. SZAMAX) goto 99

!---load the amtospheric column data
      do L = 1,L1_
        PPJ(L) = ETAA(L) + ETAB(L)*P(ILNG,JLAT)
        TTJ(L) = TJ(ILNG,JLAT,L)
        DDJ(L) = DM(ILNG,JLAT,L)
        ZZJ(L) = DO3(ILNG,JLAT,L)
      enddo
        PPJ(L1_+1) = 0.d0

!---calculate spherical weighting functions (AMF: Air Mass Factor)
      do L = 1,L1_
        ZHL(L) = ZH(ILNG,JLAT,L)
      enddo

!-----------------------------------------------------------------------
      call SPHERE(U0,RAD,ZHL,ZZHT,AMF,L1_)
!-----------------------------------------------------------------------

!---load the profiles of aerosols & clouds: treated as the same from here
      do L = 1,L1_
        AER(1,L)  = DAER1(ILNG,JLAT,L)  ! Opt. Depth aerosol 1 in layer L
        AER(2,L)  = DAER2(ILNG,JLAT,L)  ! Opt. Depth aerosol 2 in layer L
        AER(3,L)  = DAER3(ILNG,JLAT,L)  ! Opt. Depth aerosol 3 in layer L
        AER(4,L)  = ODCLD(ILNG,JLAT,L)  ! cloud Opt. Depth in L
        AER(5,L)  = 0.d0                ! save space for Rayleigh
        ADX(1,L)  = NAER1(ILNG,JLAT,L)  ! index for aerosol 1 at layer L
        ADX(2,L)  = NAER2(ILNG,JLAT,L)  ! index for aerosol 2 at layer L
        ADX(3,L)  = NAER3(ILNG,JLAT,L)  ! index for aerosol 3 at layer L
        ADX(4,L)  = NCLDX(ILNG,JLAT,L)  ! index for cloud at layer L
        ADX(5,L)  = 1                   ! index for Rayleigh phase in L
      enddo

      do L = 1,L1_
       do I=1,4
        ADX(I,L) = min(NAA, max(0, ADX(I,L)))
       enddo
      enddo

!---Now given the aerosol+cloud OD/layer in visible (600 nm) can calculate
!        how to add additonal levels at top of clouds (now uses log spacing)
!-----------------------------------------------------------------------
      call EXTRAL(AER,ADX,L1_,L2_,N_,JTAUMX,ATAU,ATAU0, JXTRA)
!-----------------------------------------------------------------------

!---set surface reflectance
        RFLECT = SA(ILNG,JLAT)
        RFLECT = max(0.d0,min(1.d0,RFLECT))

!---Loop over all wavelength bins to calc mean actinic flux AVGF(L)
      do K = NW1,NW2
        WAVE = WL(K)
!---Pick nearest Mie wavelength, no interpolation--------------
                               KM=1  ! use 300 nm aerosol properties for <355 nm
        if( WAVE .gt. 355.d0 ) KM=2  ! use 400 nm prop for 355-500 nm
        if( WAVE .gt. 500.d0 ) KM=3
        if( WAVE .gt. 800.d0 ) KM=4 

!---Loop over CTM layers L=1:L1_ = 1:L_+1,
!     values at L1_=L_+1 are a pseudo layer above the top CTM layer (L_)
        do L = 1,L1_
         TTT     = TTJ(L)
         XQO3 = FLINT(TTT,TQQ(1,2),TQQ(2,2),TQQ(3,2) &
     &                      ,QO3(K,1),QO3(K,2),QO3(K,3))
         XQO2 = FLINT(TTT,TQQ(1,1),TQQ(2,1),TQQ(3,1) &
     &                      ,QO2(K,1),QO2(K,2),QO2(K,3))

         ABX(L) = XQO3*ZZJ(L) + XQO2*DDJ(L)*0.20948d0

         AER(5,L) = DDJ(L)*QRAYL(K)
        enddo

!-----------------------------------------------------------------------
      call OPMIE(K,KM, WAVE, ABX,AER,ADX,U0,RFLECT,AMF,JXTRA,AVGF)
!-----------------------------------------------------------------------

        do L = 1,JVL_
          FFF(K,L) = FFF(K,L) + SOLF*FL(K)*AVGF(L)
        enddo

      enddo       ! loop over wavelngth K

!-----------------------------------------------------------------------
      call JRATET(PPJ,TTJ,FFF, VALJL)
!-----------------------------------------------------------------------

!---map the J-values from fast-JX onto ASAD ones (use JIND & JFACTA)
      do L = 1,JVL_
        do J = 1,NRATJ
          if (JIND(J).gt.0) then 
            ZPJ(L,J) = VALJL(L,JIND(J))*JFACTA(J)
          else
            ZPJ(L,J) = 0.d0
          endif
        enddo
      enddo

!---diagnostics that are NOT returned to the CTM code

!-----------------------------------------------------------------------
      write(6,*) 'fast-JX(5.3b)--------------internal print------------'

      call JP_ATM(PPJ,TTJ,DDJ,ZZJ,ZHL,ZZHT,AER,ABX,ADX,JXTRA)

!---Print solar flux terms
      write(6,'(a,4f10.4)') ' fast-JX:  SZA/u0/albedo/SOL*fact', &
     &  SZA,U0,RFLECT,SOLF
      write(6,'(a5,18i9)')   ' bin:',(K, K=NW1,NW2)
      write(6,'(a5,18f9.2)') ' wvl:',(WL(K), K=NW1,NW2)
      do L = JVL_,1,-1
        write(6,'(i3,2x,18f9.6)') L,(FFF(K,L)/FL(K),K=NW1,NW2)
      enddo
      write(6,*) 'fast-JX(5.3b)--------------internal print------------'
!-----------------------------------------------------------------------
      
   99 continue

      return

      end subroutine PHOTOJ
!EOC
!-----------------------------------------------------------------------

      end module m_FastjCTM_FastjLinking
