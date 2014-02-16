      module GmiShipEmission_mod
!
      implicit none
!
      private
      public   :: calcShipEmission
!
      CONTAINS
!
!=============================================================================
!
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE: calcShipEmission
!
! !INTERFACE:
!
      subroutine calcShipEmission (emiss_o3, emiss_ozone, prevRecord, &
                        curRecord, latdeg, jno2val, emissionArray, &
                        o3_index, i1, i2, ju1, j2, k1, k2, ju1_gl, &
                        j2_gl, num_emiss)
!
! !USES:
      use GmiArrayBundlePointer_mod, only : t_GmiArrayBundle
      use GmiSpeciesRegistry_mod   , only : getSpeciesMolWeight
!
      implicit none
!
! !INPUT PARAMETERS:
      integer                , intent(in   ) :: i1, i2, ju1, j2, k1, k2
      integer                , intent(in   ) :: ju1_gl, j2_gl
      integer                , intent(in   ) :: o3_index, curRecord
      integer                , intent(in   ) :: num_emiss
      real*8                 , intent(in   ) :: latdeg (ju1_gl:j2_gl)
      real*8                 , intent(in   ) :: jno2val(i1:i2, ju1:j2)
      real*8                 , intent(inOut) :: emiss_o3(i1:i2, ju1:j2)
      type (t_GmiArrayBundle), intent(inOut) :: emissionArray(num_emiss)
      integer                , intent(inOut) :: prevRecord
      real*8                 , intent(inOut) :: emiss_ozone(i1:i2, ju1:j2)
!
! !DESCRIPTION:
! This routine calculates the ship emissions.
! Only able to preprocess ship emissions if doing surface emissions inside
! the chemistry solver (smvgear). Otherwise, NO is directly emitted to code,
! which causes excessive O3 production.
! Preprocess ship NO emissions as:
!     \begin{itemize}
!     \item E(hno3) = E(no)
!     \item E(o3) = E(no)*10*(jno2/0.0095)^2
!     \end{itemize}
! Following loosely
!     Chen et al., An investigation of the chemistry of ship emission plumes
!     during ITCT 2002, J. Geophys. Res., Vol. 110, No. D10,
!     D10S90 10.1029/2004JD005236, 20 May 2005.
!
! !LOCAL VARIABLES:
      integer :: i,j
      real*8  :: tempVar
      real*8  :: mwNO      ! molecular weight of NO
      real*8  :: mwO3      ! molecular weight of O3
!
! !REVISION HISTORY:
!  Initial code. Bryan Duncan 2-28-06
!
!EOP
!-------------------------------------------------------------------------
!BOC
      mwNO = getSpeciesMolWeight('NO')
      mwO3 = getSpeciesMolWeight('O3')

      ! Check if new record (month or day).
      if (prevRecord /= curRecord) then
         emiss_ozone(:,:) = emissionArray(o3_index)%pArray3D(:,:,1)
         prevRecord       = curRecord
      end if

      emiss_o3(:,:) = 0.0d0

      do j=ju1,j2
         if ((latdeg(j).lt.60.0d0) .and. (latdeg(j).gt.-60.0d0)) then
            do i=i1,i2
               if (jno2val(i,j).le.0.0095d0) then
                  tempVar       = jno2val(i,j)/0.0095d0
                  emiss_o3(i,j) = 10.0d0 * emiss_ozone(i,j) * tempVar*tempVar
               else
                  emiss_o3(i,j) = emiss_ozone(i,j) * 10.0d0
               endif
            enddo
         endif
      enddo

      ! convert kg NO/s to molec NO/s by dividing by mw of NO (0.030kg/mole)
      ! convert molec o3/s to kg o3/s by multiplying by mw of o3 (0.048kg/mole)

      emiss_o3(:,:) = (mwO3/mwNO)* emiss_o3(:,:)

      emissionArray(o3_index)%pArray3D(:,:,1) = emiss_o3(:,:)

      return

      end subroutine calcShipEmission
!EOC
!-------------------------------------------------------------------------
      end module GmiShipEmission_mod
!
!
