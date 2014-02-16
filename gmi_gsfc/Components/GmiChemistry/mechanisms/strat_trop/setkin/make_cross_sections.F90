!=======================================================================
!
! $Id: make_cross_sections.F90,v 1.2 2011-08-09 22:52:11 mrdamon Exp $
!
! ROUTINE
!   csq_table - GMI version (make_cross_sections.F90)
!   19 JUN 02 - PSC
!
! DESCRIPTION
!   
!   
!
!
!  Chemistry input file:    10/2006
!  Reaction dictionary:     GMI_Combo_rxns_124species_SO2_JPL06.db
!  Setkin files genelookupd:  Thu Nov 19 21:21:30 2009
!
!=======================================================================
      program csq_table

      real*8 &
     &  temperature (200) &
     & ,dummy_pres  (200)

      real*8 &
     &  sflx    (126,200) &
     & ,sqr     (126,200)

      real*8 &
     &  dummy

      integer &
     &  i, j &
     & ,ichnl &
     & ,nbin0 &
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
!cc 322: O2 + hv = 2 O
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cc      call sqo2(1)

      write(10,10) ((dummy,i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 323: O3 + hv = O + O2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqo3phi(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 324: O3 + hv = O1D + O2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqo3phi(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 325: HO2 + hv = O + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqho2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 326: H2O + hv = H + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqh2o(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 327: NO + hv = N + O
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!cc      call sqno(1)

      write(10,10) ((dummy,i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 328: N2O + hv = N2 + O1D
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqn2o(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 329: NO2 + hv = NO + O
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqno2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 330: H2O2 + hv = 2 OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqh2o2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 331: MP + hv = CH2O + HO2 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 332: CH2O + hv = CO + H + HO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch2o(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 333: CH2O + hv = CO + H2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch2o(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 334: HNO3 + hv = NO2 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqhno3(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 335: HNO2 + hv = NO + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqhono(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 336: HNO4 + hv = NO3 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqho2no2(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 337: NO3 + hv = NO2 + O3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqno3(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 338: NO3 + hv = NO + O2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqno3(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 339: N2O5 + hv = NO2 + NO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqn2o5(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 340: N2O5 + hv = NO + NO3 + O3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqn2o5(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 341: HNO4 + hv = HO2 + NO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqho2no2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 342: Cl2 + hv = 2 Cl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqcl2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 343: OClO + hv = ClO + O
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqoclo(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 344: Cl2O2 + hv = 2 Cl + O2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqcl2o2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 345: HOCl + hv = Cl + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqhocl(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 346: ClONO2 + hv = Cl + NO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqclono2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 347: ClONO2 + hv = ClO + NO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqclono2(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 348: BrCl + hv = Br + Cl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqbrcl(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 349: BrO + hv = Br + O
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqbro(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 350: HOBr + hv = Br + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqhobr(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 351: BrONO2 + hv = Br + NO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqbrono2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 352: BrONO2 + hv = BrO + NO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqbrono2(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 353: CH3Cl + hv = Cl + MO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3cl(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 354: CCl4 + hv = 4 Cl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqccl4(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 355: CH3CCl3 + hv = 3 Cl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqf140(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 356: CFCl3 + hv = 3 Cl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqf11(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 357: CF2Cl2 + hv = 2 Cl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqf12(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 358: CFC113 + hv = 3 Cl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqf113(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 359: CFC114 + hv = 2 Cl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqf114(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 360: CFC115 + hv = Cl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqf115(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 361: HCFC141b + hv = 2 Cl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqhf141b(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 362: HCFC142b + hv = Cl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqhf142b(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 363: CH3Br + hv = Br + MO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3br(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 364: CF3Br + hv = Br
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqh1301(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 365: CF2Br2 + hv = 2 Br
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqh1202(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 366: H2402 + hv = 2 Br
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqh2402(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 367: CF2ClBr + hv = Br + Cl
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqh1211(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 368: ALD2 + hv = CO + HO2 + MO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3cho(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 369: ALD2 + hv = CH4 + CO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3cho(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 370: PAN + hv = MCO3 + NO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqpan(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 371: RCHO + hv = CO + ETO2 + HO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqc2h5cho(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 372: ACET + hv = MCO3 + MO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqacetone(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 373: MEK + hv = ETO2 + MCO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqmek(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 374: GLYC + hv = CH2O + CO + 2 HO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3cho(nbin0,npres0,temperature,3,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 375: GLYX + hv = 2 CO + H2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqglyoxal(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 376: GLYX + hv = 2 CO + 2 HO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqglyoxal(nbin0,npres0,temperature,3,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 377: GLYX + hv = CH2O + CO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqglyoxal(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 378: MGLY + hv = CO + HO2 + MCO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqmgly(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 379: MGLY + hv = ALD2 + CO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqmgly(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 380: MVK + hv = CO + PRPE
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqmek(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 381: MVK + hv = CH2O + CO + HO2 + MCO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqmek(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 382: MVK + hv = MAO3 + MO2
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqmek(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 383: MACR + hv = HO2 + MAO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqc2h5cho(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 384: MACR + hv =  0.20 CH2O + CO +  1.80 HO2 +  0.20 MCO3 +  0.80 MGLY
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqc2h5cho(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 385: HAC + hv = CH2O + HO2 + MCO3
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqacetone(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 386: INPN + hv = HO2 + NO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3o2no2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 387: PRPN + hv = HO2 + NO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3o2no2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 388: ETP + hv = ALD2 + HO2 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 389: RA3P + hv = HO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 390: RB3P + hv = HO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 391: R4P + hv = HO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 392: PP + hv = HO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 393: RP + hv = ALD2 + HO2 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 394: GP + hv = CH2O + HO2 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 395: RIP + hv =  0.69 CH2O +  0.86 HO2 +  0.13 IALD +  0.29 MACR +  0.40 MVK + OH +  0.14 RIO1
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 396: IAP + hv =  0.67 CO +  0.26 GLYC +  0.19 H2 +  0.36 HAC + HO2 +  0.58 MGLY + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 397: ISNP + hv = HO2 + NO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3o2no2(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 398: VRP + hv =  0.28 CH2O +  0.72 GLYC +  0.28 HO2 +  0.72 MCO3 +  0.28 MGLY + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 399: MRP + hv =  0.17 CH2O +  0.83 CO +  0.83 HAC + HO2 +  0.17 MGLY + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 400: MAOP + hv = HO2 + OH + RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 401: R4N2 + hv =  0.05 A3O2 +  0.32 ACET +  0.32 ALD2 +  0.18 B3O2 +  0.32 ETO2 +  0.27 HO2 +  0.19 MEK +  0.18 MO2 + NO2 +  0.13 RCHO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3o2no2(nbin0,npres0,temperature,2,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc 402: MAP + hv = MO2 + OH
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call sqch3ooh(nbin0,npres0,temperature,1,sqr)

      write(10,10) (((sqr(i,j)),i=1,126),j=1,200)

      write(10,10) pitest

  10  format (5e12.5)
      stop

      end
