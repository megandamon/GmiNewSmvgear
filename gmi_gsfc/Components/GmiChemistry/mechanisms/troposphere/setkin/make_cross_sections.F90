!=======================================================================
!
! $Id: make_cross_sections.F90,v 1.2 2011-08-09 22:12:58 mrdamon Exp $
!
! ROUTINE
!   csq_table - GMI version (make_cross_sections.F)
!   19 JUN 02 - PSC
!
! DESCRIPTION
!
!
!
!
!  Chemistry input file:    4:00 PM 10/28/2006
!  Reaction dictionary:     GMI_Trop_rxns_85species_JPL06.db
!  Setkin files genelookupd:  Mon Oct 30 20:16:31 2006
!
!=======================================================================
      program csq_table

      real*8  &
     &  temperature (200)  &
     & ,dummy_pres  (200)

      real*8  &
     &  sflx    (126,200)  &
     & ,sqr     (126,200)

      real*8  &
     &  dummy

      integer  &
     &  i, j  &
     & ,ichnl  &
     & ,nbin0  &
     & ,npres0

      do i             = 1, 200
        temperature(i) = i + 149.0d0
      end do

      sflx(:,:)        = 1.0d0

      ichnl  = 1
      nbin0  = 126
      npres0 = 200

      pitest = 3.14159
      dummy  = 0.0d0

      open(unit=10,file="cross_sections.txt",status="new")

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 223: O3 + hv = 2 OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqo3(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 224: NO2 + hv = NO + O3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqno2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 225: H2O2 + hv = 2 OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqh2o2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 226: MP + hv = CH2O + HO2 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 227: CH2O + hv = CO + 2 HO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch2o(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 228: CH2O + hv = CO + H2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch2o(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 229: HNO3 + hv = NO2 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqhno3(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 230: HNO2 + hv = NO + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqhono(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 231: HNO4 + hv = NO3 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqho2no2(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 232: NO3 + hv = NO2 + O3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqno3(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 233: NO3 + hv = NO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqno3(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 234: N2O5 + hv = NO2 + NO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqn2o5(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 235: N2O5 + hv = NO + NO3 + O3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqn2o5(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 236: HNO4 + hv = HO2 + NO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqho2no2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 237: ALD2 + hv = CO + HO2 + MO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3cho(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 238: ALD2 + hv = CH4 + CO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3cho(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 239: PAN + hv = MCO3 + NO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqpan(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 240: RCHO + hv = CO + ETO2 + HO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqc2h5cho(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 241: ACET + hv = MCO3 + MO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqacetone(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 242: MEK + hv = ETO2 + MCO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqmek(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 243: GLYC + hv = CH2O + CO + 2 HO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3cho(nbin0,npres0,temperature,3,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 244: GLYX + hv = 2 CO + H2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqglyoxal(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 245: GLYX + hv = 2 CO + 2 HO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqglyoxal(nbin0,npres0,temperature,3,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 246: GLYX + hv = CH2O + CO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqglyoxal(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 247: MGLY + hv = CO + HO2 + MCO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqmgly(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 248: MGLY + hv = ALD2 + CO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqmgly(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 249: MVK + hv = CO + PRPE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqmek(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 250: MVK + hv = CH2O + CO + HO2 + MCO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqmek(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 251: MVK + hv = MAO3 + MO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqmek(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 252: MACR + hv = HO2 + MAO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqc2h5cho(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 253: MACR + hv =  0.20 CH2O + CO +  1.80 HO2 +  0.20 MCO3 +  0.80 MGLY
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqc2h5cho(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 254: HAC + hv = CH2O + HO2 + MCO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqacetone(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 255: INPN + hv = HO2 + NO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3o2no2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 256: PRPN + hv = HO2 + NO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3o2no2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 257: ETP + hv = ALD2 + HO2 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 258: RA3P + hv = HO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 259: RB3P + hv = HO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 260: R4P + hv = HO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 261: PP + hv = HO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 262: RP + hv = ALD2 + HO2 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 263: GP + hv = CH2O + HO2 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 264: RIP + hv =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR +  0.40 MVK + OH +  0.14 RIO1
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 265: IAP + hv =  0.67 CO +  0.26 GLYC +  0.19 H2 +  0.36 HAC + HO2 +  0.58 MGLY + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 266: ISNP + hv = HO2 + NO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3o2no2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 267: VRP + hv =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3 +  0.28 MGLY + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 268: MRP + hv =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 MGLY + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 269: MAOP + hv = HO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 270: R4N2 + hv =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O2 +  0.32 ETO2 +  0.27 HO2 +  0.19 MEK +  0.18 MO2 + NO2 +  0.13 RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3o2no2(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 271: MAP + hv = MO2 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

      write(10,10) pitest

  10  format (5e12.5)
      stop

      end
