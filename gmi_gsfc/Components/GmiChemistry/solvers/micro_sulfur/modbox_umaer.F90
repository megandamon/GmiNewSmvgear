      module modbox_umaer
!-----------------------------------------------------------------------
!---->define accuracy
!-----------------------------------------------------------------------
      use precision
      implicit real(kdef) (a-h,o-z), integer (i-n)
!-----------------------------------------------------------------------
!---->parcel quantities unchanged during time integration
!
!     nvec          vector/loop length
!
!     h2osat        h2o saturation pressure               [#molec/m3]
!     h2ogas        h2o gas concentration                 [#molec/m3]
!     so4sat        h2so4 saturation pressure             [#molec/m3]
!     so4mfrac      so4 mass fraction in aerosol          [1]
!     so4frac       so4 mole fraction in aerosol          [1]
!     so4dens       aerosol density                       [kg/m3]
!     so4kelvin     factor for Kelvin effect              [m]
!     pathair       mean free path of air                 [m]
!     pathso4       mean free path of h2so4 in air        [m]
!     diffso4       diffusion coeff. for h2so4 in air     [m2/s]
!     viscos        dynamic viscosity of air              [Ns/m2]
!     so4crit       h2so4 conc with nuc' rate of 1/cm3    [#molec/m3]
!     so4sol        h2so4 sat pressure over solution      [#molec/m3]
!
!     pmsmerg       limit particle mass for merging       [#molec]
!     condnon       condensation coef's for non-so4 aer   [m3/#molec]
!                   incl. correction for collision geometry
!                   following Fuchs and Sutugin (1971)
!     xnonum        non-so4 aerosol number concentrations [#part/m3]
!     rwetnon       humid volume mean radius              [m]
!     xnonkelvin    Kelvin effect for non-sulfate aerosol [1]
!     diffnon       coefficient for coagulation
!                   for non-sulfate aerosol
!     so4sig        assumed geometric standard deviation  [1]
!-----------------------------------------------------------------------
      integer, parameter :: nvec=511

      real(kdef), save :: h2osat(nvec),h2ogas(nvec),so4sat(nvec),  &
     &                    so4mfrac(nvec),so4frac(nvec),  &
     &                    so4dens(nvec),so4kelvin(nvec),  &
     &                    pathair(nvec),pathso4(nvec),  &
     &                    diffso4(nvec),viscos(nvec),  &
     &                    so4crit(nvec),so4sol(nvec)

      real(kdef), allocatable, save :: pmsmerg(:,:),  &
     &                    condnon(:,:),xnonum(:,:),rwetnon(:,:),  &
     &                    xnonkelvin(:,:),diffnon(:,:)

      real(kdef), allocatable, save :: so4sig(:,:),xnonsig(:,:),  &
     &                    fracsrc(:,:),fracsrcnon(:,:),  &
     &                    so4alpha(:,:)

#ifdef DIAG
!-----------------------------------------------------------------------
!---->arrays for diagnostic quantities needed during time integration
!-----------------------------------------------------------------------
      real(kdef), pointer, save :: dxnucn(:),dxnucm(:),  &
     &                             dxcon (:),dxconn(:),  &
     &                             dxcoag(:),dxcoagn(:),  &
     &                             dxgrvn(:),dxgrvm(:)
      real(kdef), pointer, save :: fxnucn(:),fxnucm(:),  &
     &                             fxcon (:),fxconn(:),  &
     &                             fxcoag(:),fxcoagn(:),  &
     &                             fxgrvn(:),fxgrvm(:)

      real(kdef), target, save ::  &
     &                    dynucn(nvec,2),dynucm(nvec,2),  &
     &                    dycon (nvec,2),dyconn(nvec,2),  &
     &                    dycoag(nvec,2),dycoagn(nvec,2),  &
     &                    dygrvn(nvec,2),dygrvm(nvec,2)
      real(kdef), target, save ::  &
     &                    fynucn(nvec,3),fynucm(nvec,3),  &
     &                    fycon (nvec,3),fyconn(nvec,3),  &
     &                    fycoag(nvec,3),fycoagn(nvec,3),  &
     &                    fygrvn(nvec,3),fygrvm(nvec,3)

      real(kdef), save :: fsavenucn(nvec),fsavenucm(nvec),  &
     &                    fsavecon (nvec),fsaveconn(nvec),  &
     &                    fsavecoag(nvec),fsavecoagn(nvec),  &
     &                    fsavegrvn(nvec),fsavegrvm(nvec)

      real(kdef), save :: xnucn(nvec),xnucm(nvec),  &
     &                    xcon (nvec),xconn(nvec),  &
     &                    xcoag(nvec),xcoagn(nvec),  &
     &                    xgrvn(nvec),xgrvm(nvec),  &
     &                    xbeg(2,nvec)

#endif
      end module modbox_umaer
