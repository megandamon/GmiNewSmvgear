!=============================================================================
!
! CODE DEVELOPER
!   originated from ???
!   modified by Stephen Steenrod
!
! FILE
!   nh3chem_solver.F90
!
! ROUTINES
!   Do_NH3_Chem 
!
! HISTORY
!   - 14 Apr 2010 Steve Steenrod
!          initial development
!=============================================================================

!-----------------------------------------------------------------------------
!
! ROUTINE
!   Do_NH3_Solver 
!
! DESCRIPTION
!
!     *****************************************************************
!     * carries out the ammonium chemistry calculation
!     *****************************************************************
!
!     *----------------------------------------------------------------
!
! ARGUMENTS
!   i1, i2, ju1, j2, k1, k2: array dimensions
!   kel           : temperature (K)
!   concentration : species concentration, known at zone centers
!                   (molecules/cm^3 for gas-species and g/g for aerosols)
!   relhum        : relative humidity [0 - 1]
!
!-----------------------------------------------------------------------------

      subroutine Do_NH3_Solver  &
     &   (loc_proc, pr_diag, i1, i2, ju1, j2, k1, k2, ilo, ihi, julo, jhi,   &
     &    kel, relhum, concentration)

! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

      implicit none

#     include "setkin_par.h"
#     include "setkin_mw.h"
#     include "gmi_phys_constants.h"

      integer, intent(in) :: loc_proc
      logical, intent(in) :: pr_diag
      integer, intent(in) :: i1,  i2,  ju1,  j2, k1, k2
      integer, intent(in) :: ilo, ihi, julo, jhi
      real*8, intent(in)  :: kel (ilo:ihi, julo:jhi, k1:k2)
      real*8, intent(in)  :: relhum   (i1:i2, ju1:j2, k1:k2)
      type (t_GmiArrayBundle), intent(inOut) :: concentration(NSP)

!... local variables
      integer :: i, j, k
      real*8 :: SO4, GHNO3, GNH3, RH, TEMP, ASO4, AHSO4, ANO3, ANH4, AH2O
      real*8 :: SO43d(i1:i2, ju1:j2, k1:k2)
!      real*8, parameter :: mw_SO4 = 9.8074d+01

#ifdef nonZeroInd_tracers
      integer, parameter :: INH3_1     = 1
      integer, parameter :: INH4A_1    = 1
      integer, parameter :: IHNO3_1    = 1
      integer, parameter :: INO3A_1    = 1
      integer, parameter :: IFSO4a_1   = 1
      integer, parameter :: INSO4a_1   = 1
      integer, parameter :: IFSO4n1_1  = 1
      integer, parameter :: IFSO4n2_1  = 1
      integer, parameter :: IFSO4n3_1  = 1
      integer, parameter :: INSO4n1_1  = 1
      integer, parameter :: INSO4n2_1  = 1
      integer, parameter :: INSO4n3_1  = 1
      integer, parameter :: INO3An1_1  = 1
      integer, parameter :: IMGAS_1    = 1
#else
      integer, parameter :: INH3_1     = 1 ! INH3
      integer, parameter :: INH4A_1    = 1 ! INH4A
      integer, parameter :: IHNO3_1    = 1 ! IHNO3
      integer, parameter :: INO3A_1    = 1 ! INO3A
      integer, parameter :: IFSO4a_1   = 1 ! IFSO4a
      integer, parameter :: INSO4a_1   = 1 ! INSO4a
      integer, parameter :: IFSO4n1_1  = 1 ! IFSO4n1
      integer, parameter :: IFSO4n2_1  = 1 ! IFSO4n2
      integer, parameter :: IFSO4n3_1  = 1 ! IFSO4n3
      integer, parameter :: INSO4n1_1  = 1 ! INSO4n1
      integer, parameter :: INSO4n2_1  = 1 ! INSO4n2
      integer, parameter :: INSO4n3_1  = 1 ! INSO4n3
      integer, parameter :: INO3An1_1  = 1 ! INO3An1
      integer, parameter :: IMGAS_1    = IMGAS
#endif

!     ----------------
!     Begin execution.
!     ----------------

      if (pr_diag) then
        Write (6,*) 'Do_NH3_Solver called by ', loc_proc
      end if

    
!... get total sulfate into micrograms/m^3
!... these are in volume mixing ratio
      SO43d(:,:,:) = concentration(IFSO4n1_1)%pArray3D(:,:,:) +  &
     &               concentration(IFSO4n2_1)%pArray3D(:,:,:) +  &
     &               concentration(IFSO4n3_1)%pArray3D(:,:,:) +  &
     &               concentration(INSO4n1_1)%pArray3D(:,:,:) +  &
     &               concentration(INSO4n2_1)%pArray3D(:,:,:) +  &
     &               concentration(INSO4n3_1)%pArray3D(:,:,:)

!... presently need to convert above SO4 aerosols to mass mixing ratio
      SO43d(:,:,:) = SO43d(:,:,:) * mw_data(IFSO4n1_1) / MWTAIR

!... add in aqueous sulfate (already in mass mixing ratio)
      SO43d(:,:,:) = concentration(IFSO4a_1)%pArray3D(:,:,:) +  &
     &               concentration(INSO4a_1)%pArray3D(:,:,:) + SO43d(:,:,:)

!... convert to micrograms/m^3
      SO43d(:,:,:) = SO43d(:,:,:) *(concentration(IMGAS_1)%pArray3D(:,:,:)* MWTAIR/AVOGAD) * 1.0d12
      
      do k=k1,k2
        do j=ju1,j2
          do i=i1,i2
!.. convertion to 10-6 g/m^3 from 

!... aerosols

!     * Calculate water vapor concentration in molecules/cm^3
! humidity(ijk) (g/kg) * * MWTAIR / MWTH2O / GPKG * xr(ijk, IMGAS)

!  GNO3  : Nitric Acid vapor in MICROGRAMS/M**3 as nitric acid 
!            GHNO3 = concentration(IHNO3)%pArray3D(i,j,k)  *mw_data(IHNO3)/mw_data(IMGAS)
            GHNO3 = concentration(IHNO3_1)%pArray3D(i,j,k)  &
     &             * (concentration(IMGAS_1)%pArray3D(i,j,k)* mw_data(IHNO3_1)/AVOGAD) * 1.0d12
     
!  GNH3  : Gas phase ammonia in MICROGRAMS/M**3 
!            GNH3 = concentration(INH3)%pArray3D(i,j,k)  *mw_data(INH3)/mw_data(IMGAS)
            GNH3 = concentration(INH3_1)%pArray3D(i,j,k)  &
     &             * (concentration(IMGAS_1)%pArray3D(i,j,k)* mw_data(INH3_1)/AVOGAD) * 1.0d12
     
!  RH    : Fractional relative humidity 
            RH = relhum(i,j,k)
            
!  TEMP  : Temperature in Kelvin 
            TEMP = kel(i,j,k)
            
!  SO4  : Total Aerosol phase sulfate in MICROGRAMS/M**3
            SO4 = SO43d(i,j,k)

!  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 
            if(INO3A .ne. 0) then
               ANO3 = concentration(INO3A_1)%pArray3D(i,j,k)  &
     &              *(concentration(IMGAS_1)%pArray3D(i,j,k)* MWTAIR/AVOGAD) * 1.0d12
            endif
            if(INO3An1 .ne. 0) then
               ANO3 = concentration(INO3An1_1)%pArray3D(i,j,k)  &
     &              *(concentration(IMGAS_1)%pArray3D(i,j,k)* mw_data(INO3An1_1)/AVOGAD) * 1.0d12
            endif

!  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3 
            ANH4 = concentration(INH4A_1)%pArray3D(i,j,k)   &
     &             *(concentration(IMGAS_1)%pArray3D(i,j,k)* MWTAIR/AVOGAD) * 1.0d12

!  AHSO4 : Aerosol phase in bisulfate in MICROGRAMS/M**3 [rjp, 12/12/01]
            AHSO4 = 1e-30

!  AH2O  : Aerosol phase water in MICROGRAMS/M**3 
            AH2O = 1e-30

!... call ammonia solver           
            call RPMARES( SO4,  GHNO3,  GNH3, RH,   TEMP,  &
     &                    ASO4, AHSO4, ANO3, AH2O, ANH4 )

!  GNH3  : Gas phase ammonia in MICROGRAMS/M**3  convert to volume mixing ratio
!            concentration(INH3)%pArray3D(i,j,k) = GNH3 *mw_data(IMGAS)/mw_data(INH3)
            concentration(INH3_1)%pArray3D(i,j,k) = GNH3  &
     &              / ( concentration(IMGAS_1)%pArray3D(i,j,k)* mw_data(INH3_1)/AVOGAD * 1.0d12 )

!  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 convert to mass mixing ratio
            if(INO3A .ne. 0) then
               concentration(INO3A_1)%pArray3D(i,j,k) = ANO3  &
     &              / ( concentration(IMGAS_1)%pArray3D(i,j,k)* MWTAIR/AVOGAD * 1.0d12 )
            endif
            if(INO3An1 .ne. 0) then
               concentration(INO3An1_1)%pArray3D(i,j,k) = ANO3  &
     &              / ( concentration(IMGAS_1)%pArray3D(i,j,k)* mw_data(INO3An1_1)/AVOGAD * 1.0d12 )
            endif

!  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3  convert to mass mixing ratio
            concentration(INH4A_1)%pArray3D(i,j,k) = ANH4  &
     &              / ( concentration(IMGAS_1)%pArray3D(i,j,k)* MWTAIR/AVOGAD * 1.0d12 )

!!  HNO3  : Gas phase nitric acid in MICROGRAMS/M**3  convert to vol mixing ratio
!            concentration(IHNO3)%pArray3D(i,j,k) = GHNO3  &
!     &              / ( concentration(IMGAS)%pArray3D(i,j,k)* mw_data(IHNO3)/AVOGAD * 1.0d12 )


          enddo
        enddo
      enddo

      return

      end
