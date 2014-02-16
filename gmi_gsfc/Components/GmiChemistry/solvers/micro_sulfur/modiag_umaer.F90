      module modiag_umaer
!-----------------------------------------------------------------------
!---->define accuracy
!-----------------------------------------------------------------------
      use precision
      implicit real(kdef) (a-h,o-z), integer (i-n)
!-----------------------------------------------------------------------
!---->diagnostic quantities
!
!     so4nucn   - number of newly nucleated aerosol
!                 particles since the beginning of
!                 the simulation                       in [#particles/m3]
!     so4nucm   - mass of newly nucleated aerosol
!                 particles since the beginning of
!                 the simulation                       in [#so4 molecules/m3]
!     so4con    - condensation since the beginning
!                 of the simulation                    in [#so4 molecules/m3]
!     so4coag   - change in number concentration
!                 since the beginning of the simulation
!                 due to coagulation                   in [#particles/m3]
!     so4grvn   - change in number concentration
!                 since the beginning of the simulation
!                 due to gravitational settling        in [#particles/m3]
!     so4grvm   - change in mass concentration
!                 since the beginning of the simulation
!                 due to gravitational settling        in [#so4 molecules/m3]
!-----------------------------------------------------------------------
      real(kdef), allocatable, target, save ::  &
     &                      so4nucn(:),so4nucm(:),  &
     &                      so4con (:),so4conn(:),  &
     &                      so4coag(:),so4coagn(:),  &
     &                      so4grvn(:),so4grvm(:),  &
     &                      begnuc(:,:)

      end module modiag_umaer
