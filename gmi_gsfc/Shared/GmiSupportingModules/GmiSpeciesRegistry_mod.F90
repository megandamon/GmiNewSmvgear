!------------------------------------------------------------------------------
! NASA GSFC - SIVO Code 610.3
!------------------------------------------------------------------------------
!BOP
!
#include "GmiESMF_ErrLog.h"
!
 module GmiSpeciesRegistry_mod

 use ESMF_Mod
 use GmiESMF_ErrorChecking_mod
 use GmiPrintError_mod        , only : GmiPrintError
 use GmiStringManipulation_mod, only : stringLowerCase, constructListNames
 use GmiESMFrcFileReading_mod , only : rcEsmfReadTable, reconstructPhrase

 implicit none

 private
 public  :: getSpeciesIndex, getSpeciesPosition
 public  :: setNumberSpecies, set_labelsSpecies
 public  :: set_molWeightSpecies, get_molWeightSpecies
 public  :: getSpeciesMolWeight
 public  :: UNKNOWN_SPECIES

# include "GmiParameters.h"
# include "setkin_par.h"
# include "setkin_lchem.h"
# include "setkin_mw.h"

 integer, parameter :: UNKNOWN_SPECIES         = -2

 integer, save :: numSpecies
 real*8                     , save , pointer :: molWeightSpecies(:) => null()
 character (len=ESMF_MAXSTR), save , pointer :: labelsSpecies   (:) => null()
!EOP
!------------------------------------------------------------------------------
 contains
!------------------------------------------------------------------------------
!BOP
!
 subroutine setNumberSpecies(num_species)

 integer, intent(in) :: num_species
!EOP
!------------------------------------------------------------------------------
!BOC

 numSpecies = num_species

 return

 end subroutine setNumberSpecies
!EOC
!------------------------------------------------------------------------------
!BOP
!
 subroutine set_labelsSpecies(config)

 type (Esmf_Config), intent(inOut) :: config
 integer :: rc, STATUS
 character(len=MAX_STRING_LENGTH) :: tempListNames
 character(len=ESMF_MAXSTR), parameter :: IAm = "set_labelsSpecies"
!EOP
!------------------------------------------------------------------------------
!BOC

      call ESMF_ConfigGetAttribute(config, numSpecies, &
     &                label   = "numSpecies:",&
     &                default = 1, rc=STATUS )
      VERIFY_(STATUS)

      allocate(labelsSpecies(numSpecies))

      call ESMF_ConfigFindLabel(config, label="const_labels::", rc=STATUS)

      if (STATUS == ESMF_SUCCESS) then
         call rcEsmfReadTable (config, tempListNames, &
     &              "const_labels::", rc=STATUS)
         VERIFY_(STATUS)
         labelsSpecies(1:numSpecies) = ''
         call constructListNames(labelsSpecies, tempListNames)
      else 
         labelsSpecies(1:numSpecies) = lchemvar(1:numSpecies)
      end if

 return

 end subroutine set_labelsSpecies
!EOC
!------------------------------------------------------------------------------
!BOP
!
 subroutine set_molWeightSpecies(config)

 type (Esmf_Config), intent(inOut) :: config
 integer :: rc, STATUS
 character(len=ESMF_MAXSTR), parameter :: IAm = "set_molWeightSpecies"
!EOP
!------------------------------------------------------------------------------
!BOC

      allocate(molWeightSpecies(numSpecies))

      call ESMF_ConfigFindLabel(config, label="mw::", rc=STATUS)

      if (STATUS == ESMF_SUCCESS) then
         call rcEsmfReadTable (config, molWeightSpecies, &
     &              "mw::", rc=STATUS)
         VERIFY_(STATUS)
      else
         molWeightSpecies(1:numSpecies) = mw_data(1:numSpecies)
      end if

 return

 end subroutine set_molWeightSpecies
!EOC
!------------------------------------------------------------------------------
!BOP
!
 subroutine get_molWeightSpecies(mw)

 real*8, intent(out) :: mw(numSpecies)
!EOP
!------------------------------------------------------------------------------
!BOC

 allocate(molWeightSpecies(numSpecies))
 
 mw(1:numSpecies) = molWeightSpecies(1:numSpecies)

 return

 end subroutine get_molWeightSpecies
!EOC
!------------------------------------------------------------------------------
!BOC
!
! !IROUTINE: getSpeciesMolWeight
!
! !INTERFACE:
!
      function getSpeciesMolWeight(name) result(molWeight)
!
      implicit none
!
! !INPUT PARAMETERS:
      character (len=*), intent(in) :: name
!
! !RETURN VALUE:
      real*8 molWeight
!
! !DESCRIPTION:
! Given a species name, this function returns its molecular weight.
!
! !LOCAL VARIABLES:
      integer :: index
!EOP
!------------------------------------------------------------------------------
!BOC      
      index     = getSpeciesIndex (name )    ! get the species index
      molWeight = molWeightSpecies(index)
      
      return

      end function getSpeciesMolWeight
!EOC      
!------------------------------------------------------------------------------
!BOP
!
 function getSpeciesIndex(name) result(index)
 character (len=*), intent(in) :: name
 integer                      :: index
 character (len=MAX_LENGTH_ERROR_MSG)  :: err_msg

 integer :: ii
!EOP
!------------------------------------------------------------------------------
!BOC
 index = UNKNOWN_SPECIES

 if (trim(stringLowerCase(name)) == 'xxx') then
    index = -1  ! this is mainly used for emiss_map
 else
    do ii = 1, numSpecies
       if (trim(stringLowerCase(labelsSpecies(ii))) == trim(stringLowerCase(name))) then
          index = ii
          exit
       end if
    end do
 end if

 if (index == UNKNOWN_SPECIES) then
    err_msg = 'The species does not exist to index: '// name
    call GmiPrintError(err_msg, .true., 1, index, 0, 0, 0.0d0, 0.0d0)
 end if

 end function getSpeciesIndex
!EOC
!------------------------------------------------------------------------------
!BOP
!
 function getSpeciesPosition(speciesName, speciesList, numSpc) result(index)
 
 implicit none

! !INPUT PARAMETERS:
  integer          , intent(in) :: numSpc
  character (len=*), intent(in) :: speciesName
  character (len=*), intent(in) :: speciesList(numSpc)
!
! !RETURN VALUE:
 integer                  :: index
!
! !DESCRIPTION:
! Giving a name and a list of species names, this function returns the
! index prosition of the name in the list.
!
! !LOCAL VARIABLES:
  character (len=MAX_LENGTH_ERROR_MSG)  :: err_msg
  integer :: ii
!EOP
!------------------------------------------------------------------------------
!BOC
  index = UNKNOWN_SPECIES

  do ii = 1, numSpc
     if (trim(stringLowerCase(speciesList(ii))) == trim(stringLowerCase(speciesName))) then
        index = ii
        exit
     end if
  end do

 if (index == UNKNOWN_SPECIES) then
    err_msg = 'The species does not exist: '// speciesName
    call GmiPrintError(err_msg, .true., 1, index, 0, 0, 0.0d0, 0.0d0)
 end if

 end function getSpeciesPosition
!EOC
!------------------------------------------------------------------------------

 end module GmiSpeciesRegistry_mod
