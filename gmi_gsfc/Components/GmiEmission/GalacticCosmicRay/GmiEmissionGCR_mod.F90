  module GmiEmissionGCR_mod

  use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle

  implicit none

!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: Add_Emiss_GCR
  
  contains

!=============================================================================
!-----------------------------------------------------------------------------
!
! ROUTINE
!   Add_Emiss_Gsfc
!
! DESCRIPTION
!   This routine adds emissions to const.
!
!     1) Take emissions of moecules/s and multiply by the time step to get total
!        molecules of emission over time step.
!
!     2) Divide by mass of the zone to obtain emissions in terms of mixing
!        ratio; also multiply by the ratio of molecular weight of air to
!        molecular weight of the chemical emission to get volume mixing
!        ratio from mass mixing ratio.
!
!     3) Add emitted mixing ratio amount to existing mixing ratio of const.
!
! ARGUMENTS
!   mass           : total mass of the atmosphere within each grid box (kg)
!   const          : species concentration, known at zone centers
!                    (mixing ratio)
!   tdt            : model time step
!
! CODE DEVELOPER
!   Stephen Steenrod ; GSFC
!   stephen.d.steenrod@nasa.gov
!-----------------------------------------------------------------------------

  subroutine Add_Emiss_GCR  &
   (gcr_sunspot, gcr_slope, gcr_aintcp, concentration, IMGAS, tdt,  &
    pnox, GCR_NOx, mass, pr_diag, loc_proc, i1, i2, ju1, j2, k1, k2, num_species)

!     ----------------------
!     Argument declarations.
!     ----------------------

  real*8 , intent(in ) :: gcr_sunspot                 ! Galactic Cosmic Ray parameter
  real*8 , intent(in ) :: gcr_slope (ju1:j2, k1:k2)   ! Galactic Cosmic Ray parameter
  real*8 , intent(in ) :: gcr_aintcp (ju1:j2, k1:k2)  ! Galactic Cosmic Ray parameter
  type (t_GmiArrayBundle), intent(in) :: concentration(num_species)
  integer, intent(in ) :: IMGAS
  real*8 , intent(in ) :: tdt                          ! time step
  real*8 , intent(out) :: pnox ( i1:i2, ju1:j2, k1:k2) ! Galactic Cosmic Ray N/NO production (ppv/s)
  real*8 , intent(out) :: GCR_NOx (i1:i2, ju1:j2, k1:k2)
  real*8 , intent(in)  :: mass (i1:i2, ju1:j2, k1:k2)
  logical, intent(in ) :: pr_diag
  integer, intent(in ) :: loc_proc, i1, i2, ju1, j2, k1, k2, num_species

! ----------------------
! Variable declarations.
! ----------------------

  integer :: il, ij, ik

  real*8  :: work, ppmtoppv


! ----------------
! Begin execution.
! ----------------

  if (pr_diag) Write (6,*) 'Add_Emiss_GCR called by ', loc_proc

!... ppmtoppv --> converts from mass of air to mass of N.
  ppmtoppv = 28.97 / 14.

! -------------------------------------------
! Galactic Cosmic Ray emission of N and NO
! -------------------------------------------
  do ik = k1, k2
    do ij = ju1, j2
!... Calc total emission of NOx (input molec/cm^3/s - convert to molec/cm^3)
      work = gcr_sunspot * gcr_slope(ij,ik) + gcr_aintcp(ij,ik)

      do il = i1, i2
!.... production of NOx to mixing ratio
        pnox(il,ij,ik) = work  / concentration(IMGAS)%pArray3D(il,ij,ik)

      enddo
    enddo
  enddo
!... diagnostic for Galactic Cosmic Ray NOx emissions
!...  Nitrogen production in  kg/timestep
  GCR_NOx(:,:,:) =  (pnox(:,:,:)*mass(:,:,:)*tdt) / ppmtoppv

  return
  end subroutine Add_Emiss_GCR
!-------------------------------------------------------------------------
!EOC
  end module GmiEmissionGCR_mod
