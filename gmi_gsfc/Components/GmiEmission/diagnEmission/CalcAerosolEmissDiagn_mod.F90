!-------------------------------------------------------------------------
!  NASA/GFSC, SIVO, Code 610.3
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  CalcAerosolEmissDiagn_mod
!
! !INTERFACE:
!
module CalcAerosolEmissDiagn_mod
!
    implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
    private
    public  :: calcAerosolEmissDiagn
    public  :: setAerosolSurfEmissMap
!
! !DESCRIPTION:
!  Routines for producing emission diagnostics for aerosols.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  8February2007 Initial code.
!
!EOP
!-------------------------------------------------------------------------
    contains

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setAerosolSurfEmissMap
!
! !INTERFACE:
!
  subroutine setAerosolSurfEmissMap (aerosolSurfEmissMap, emiss_map_aero, &
                        emiss_map_dust, naero, ndust)

  implicit none

# include "setkin_par.h"
!
! !INPUT PARAMETERS:
  integer, intent(in) :: ndust
                         ! number of dust species
  integer, intent(in) :: naero
                         ! number of aerosol species
  integer, intent(in) :: emiss_map_aero(naero)
                         ! mapping of aerosol emiss number to const species #
                         ! OC, BC and sea salt.
  integer, intent(in) :: emiss_map_dust(ndust)
                         ! mapping of dust    emiss number to const species #
                         ! DUST1, DUST2, DUST3, DUST4, and/or DUST5.
!
! !OUTPUT PARAMETERS:
  integer, intent(out) :: aerosolSurfEmissMap(ndust+naero+5)
                         ! mapping of SurfEmiss number for aerosols to const species #
!
! !DESCRIPTION:
!  Map the aerosol sufrace emission species indices to the const array.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
  aerosolSurfEmissMap(1)                     = IFSO2
  aerosolSurfEmissMap(2)                     = INSO2
  aerosolSurfEmissMap(3)                     = INDMS
  aerosolSurfEmissMap(4)                     = INSO4A
  aerosolSurfEmissMap(5)                     = IFSO4A

  ! Mapping for aerosol species (OC, BC and sea salt)
  aerosolSurfEmissMap(6:naero+5)             = emiss_map_aero(1:naero)

  ! Mapping for dust    species (DUST1, DUST2, DUST3, DUST4, and/or DUST5)
  aerosolSurfEmissMap(naero+6:naero+ndust+5) = emiss_map_dust(1:ndust)

  return

  end subroutine setAerosolSurfEmissMap
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calcAerosolEmissDiagn
!
! !INTERFACE:
!
  subroutine  calcAerosolEmissDiagn(aerosolEmiss3D, aerosolSurfEmiss, aerosolSurfEmissMap, & 
                    emissionArray, emiss_dust, emiss_aero, tdt, mcor, &
                    pr_surf_emiss, pr_emiss_3d, &
                    emiss_aero_opt, emiss_dust_opt, emiss_map, num_species, &
                    ndust, naero, i1, i2, ju1, j2, k1, k2, num_emiss)
!
! !USES:
!
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
!
      implicit none
! 
! !INPUT PARAMETERS:
      integer, intent(in) :: i1, i2, ju1, j2, k1, k2
      integer, intent(in) :: num_emiss, naero, ndust, num_species
      integer, intent(in) :: emiss_aero_opt, emiss_dust_opt
      integer, intent(in) :: emiss_map(num_emiss)
      integer, intent(in) :: aerosolSurfEmissMap(ndust+naero+5)
      real*8 , intent(in) :: mcor(i1:i2, ju1:j2)
      real*8 , intent(in) :: tdt
      real*8 , intent(in) :: emiss_dust(i1:i2, ju1:j2, ndust)
      real*8 , intent(in) :: emiss_aero(i1:i2, ju1:j2, naero)
      logical, intent(in) :: pr_surf_emiss, pr_emiss_3d
      type (t_GmiArrayBundle), intent(in) :: emissionArray(num_emiss)
!
! !INPUT/OUTPUT PARAMETERS:
      real*8 , intent(inout) :: aerosolSurfEmiss(i1:i2, ju1:j2, ndust+naero+5)
                   ! Array storing the surface emission diagnostics for 
                   ! aerosols. The number 5 is for species
                   ! FSO2, NSO2, DMS, NSO4A, and FSO4A.
      real*8 , intent(inout) :: aerosolEmiss3D(i1:i2, ju1:j2, k1:k2, 5)
                   ! Array storing the 3D emission diagnostics for 
                   ! aerosols. The number 5 is for species
                   ! FSO2, NSO2, DMS, NSO4A, and FSO4A.
!
! !DESCRIPTION:
! This routine calculates the 3D and surface emission for aerosols.
!
! !LOCAL VARIABLES:
      integer       :: ic, ix, ik, inum
      integer, save :: in_emiss_map(5) = -1
      logical, save :: first = .true.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
      if (first) then
         ! Identify the aerosol surface emission species
         ! (FSO2, NSO2, DMS, NSO4A, FSO4A) that are also emitted.
         inum = 0
         do ix= 1, num_emiss
            if (emiss_map(ix) > 0) then
               inum = inum + 1
               do ic = 1, 5
                  if (emiss_map(ix) == aerosolSurfEmissMap(ic)) &
                     in_emiss_map(ic) = inum
               end do
            end if
         end do
         first = .false.
      end if

      if (pr_surf_emiss) then
         ! Surface emission for FSO2, NSO2, DMS, NSO4A, FSO4A
         ! Only perform the calculation if the species is emitted.
         do ic = 1, 5
            ix = in_emiss_map(ic)
            if (ix > 0) then
               aerosolSurfEmiss(:,:,ic) = aerosolSurfEmiss(:,:,ic) +  &
     &             (emissionArray(ix)%pArray3D(:,:,1) * tdt / mcor(:,:))
            end if
         end do

         ! Surface emission for BC, OC and sea salt
         if (emiss_aero_opt > 0) then
            do ic = 6, naero+5
               ix = ic - 5
               aerosolSurfEmiss(:,:,ic) = aerosolSurfEmiss(:,:,ic) + &
     &                (emiss_aero(:,:,ix) * tdt / mcor(:,:))
            end do
         end if

         ! Surface emission for dust
         if (emiss_dust_opt > 0) then
            do ic = naero+6, ndust+naero+5
               ix = ic - naero -5
               aerosolSurfEmiss(:,:,ic) = aerosolSurfEmiss(:,:,ic) + &
     &                (emiss_dust(:,:,ix) * tdt / mcor(:,:))
            end do
         end if
      end if

      if (pr_emiss_3d) then
         ! 3D emission for FSO2, NSO2, DMS, NSO4A, FSO4A
         ! Only perform the calculation if the species is emitted.
         do ic = 1, 5
            ix = in_emiss_map(ic)
            if (ix > 0) then
               do ik = k1, k2
                  aerosolEmiss3D(:,:,ik,ic) = aerosolEmiss3D(:,:,ik,ic) +  &
     &               (emissionArray(ix)%pArray3D(:,:,ik) * tdt / mcor(:,:))
               end do
            end if
         end do
      end if

      return

      end subroutine calcAerosolEmissDiagn
!EOC
!-----------------------------------------------------------------------------------

end module CalcAerosolEmissDiagn_mod
