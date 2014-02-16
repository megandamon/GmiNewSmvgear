      module modpar_umaer
!-----------------------------------------------------------------------
!---->kdef             default accuracy, set in module precision
!-----------------------------------------------------------------------
      use precision
      implicit real(kdef) (a-h,o-z), integer (i-n)
!-----------------------------------------------------------------------
!---->problem size
!
!     nso4             number of sulfate aerosol modes
!     nmomso4          number of sulfate aerosol moments
!     nnon             number of non-sulfate aerosol modes
!     nmomnon          number of non-sulfate aerosol moments
!     ng               number of lognormal distr. in background aerosol
!     naer             number of modes and moments for sulfate aerosol
!     ntot             total number of prognostic variables
!-----------------------------------------------------------------------
      integer, save :: nso4,nmomso4,nnon,nmomnon,ng
      integer, save :: naer,ntot
!-----------------------------------------------------------------------
!---->r<n>             real numbers
!-----------------------------------------------------------------------
      real(kdef), parameter  :: r0=0.d0, r1=1.d0, r2=2.d0
!-----------------------------------------------------------------------
!---->time stepping
!
!     dtmin            minimum time step                [s]
!     dtmax            maximum time step                [s]
!     itermax          maximum number of iterations per period
!     xscalemin(ntot)  minimum allowed absolute error   [#/m3]
!-----------------------------------------------------------------------
      real(kdef), parameter :: dtmin=1.d-2,  &
     &                         dtmax=3600.d0

      integer,    parameter :: itermax=300000

      real(kdef), allocatable, save ::xscalemin(:)
!-----------------------------------------------------------------------
!---->physical constants and numbers
!
!     pi               pi                           [1]
!     epsilon          epsilon                      [1]
!     avno             Avogadro number              [mol-1]
!     rgas             gas constant                 [J/K/mol]
!     boltz            Boltzmann constant           [J/K]
!     <spec>mol        molecular weights            [kg]
!     so4min           mass of so4 molecule         [kg/molecule]
!     grav             gravitational constant       [m/s2]
!-----------------------------------------------------------------------
      real(kdef), parameter :: pi     =3.141592653589793d0,  &
     &                         epsilon=1.d-30,  &
     &                         avno   =6.02252d23,  &
     &                         rgas   =8.314d0,  &
     &                         boltz  =1.380485d-23

      real(kdef), parameter :: drymol =28.966d-3,  &
     &                         h2omol =18.017d-3,  &
     &                         so4mol =98.080d-3

      real(kdef), parameter :: grav=9.81d0

      real(kdef), parameter :: so4min=so4mol/avno
!-----------------------------------------------------------------------
!---->merging
!
!     nmerg            number of bins for table lookup
!                      of partial integral over lognormal distribution
!
!     fracmerg(nmerg+1)      fraction merged particles            [1]
!     radgmerg(nmerg+1,nso4) geometric radius that corresponds to fracmerg [m]
!     radmerg(nso4)          limit radius of merged particles     [m]
!     pmsfac(nso4)           factor for mass of a merged particle [#molec m3/kg]
!-----------------------------------------------------------------------
      integer, parameter            :: nmerg=20
      real(kdef), allocatable, save :: fracmerg(:),radgmerg(:,:),  &
     &                                 radmerg(:),pmsfac(:)
!-----------------------------------------------------------------------
!---->accumulation mode
!
!     naccu            number of bins for table lookup
!                      of partial integral over lognormal distribution
!
!     fracaccu(naccu+1)      fraction particles in accumulation mode  [1]
!     radgaccu(naccu+1,nso4) geometric radius that corresponds to fracaccu [m]
!     radaccu(2)             lower & upper radius of accu. mode       [m]
!     racculim               limit radius for which size distribution [m]
!                            is in the middle of accumulation mode
!-----------------------------------------------------------------------
      integer, parameter            :: naccu=20
      real(kdef), allocatable, save :: fracaccu(:,:),radgaccu(:,:,:),  &
     &                                 fraccnon(:),radgaccnon(:,:)
      real(kdef), save              :: radaccu(2),racculim
!-----------------------------------------------------------------------
!---->aerosol properties
!
!     radatkin           Atkin radius                   [m]
!     radmin             minimun radius of a nucleus    [m]
!
!     sulfate aerosol in near SO2 source regions:
!     treated as direct particle emissions
!     radgsrc(ng)        geometric radius               [m]
!     siggsrc(ng)        geometric standard deviation   [1]
!     pmssrc(ng)         sulfate particle mass          [#molec/particle]
!
!     background sulfate aerosol:
!     radgso4(ng)        geometric radius               [m]
!     siggso4(ng)        geometric standard deviation   [1]
!     fracso4(ng)        relative fraction              [1]
!     radvso4            dry volume mean radius         [m]
!
!     background non-sulfate aerosol:
!     radgnon(ng,nnon)   geometric radius               [m]
!     siggnon(ng,nnon)   geometric standard deviation   [1]
!     fracnon(ng,nnon)   relative fraction              [1]
!     rhonon(nnon)       aerosol density                [kg/m3]
!     xmolnon(nnon)      molecular weight               [kg/mol]
!     pmsnon(nnon)       particle mass                  [kg/particle]
!     radvnon(nnon)      volume mean radius             [m]
!
!---->accommodation coefficient
!
!     etaso4(nso4)       for sulfate aerosol            [1]
!     etanon(nnon)       for non-sulfate aerosol        [1]
!
!---->correction factor for condensation
!     for using the volume mean radius instead of the mean radius
!    (alphaso4           for sulfate aerosol            [1])
!     alphanon(nnon)     for non-sulfate aerosol        [1]
!
!---->crh[1-4]           coefficient for humidity growth
!                        of sea salt aerosol
!-----------------------------------------------------------------------
      real(kdef),parameter :: radatkin=0.05d-6,  &
     &                        radmin=3.d-10

      real(kdef),save :: radvso4

      real(kdef), allocatable, save :: radgsrc(:),siggsrc(:),  &
     &                    radgso4(:),siggso4(:),fracso4(:),  &
     &                    radgnon(:,:),siggnon(:,:),fracnon(:,:),  &
     &                    rhonon(:),xmolnon(:)
      integer, allocatable, save :: nrh_flg(:)

      real(kdef), allocatable, save :: pmssrc(:)
      real(kdef), allocatable, save :: pmsnon(:),radvnon(:)
      real(kdef), allocatable, save :: etaso4(:),etanon(:)
      real(kdef), allocatable, save :: alphanon(:)

!     real(kdef), parameter :: crh1= 1.10413d0,
!    $                         crh2= 3.07900d0,
!    $                         crh3= 3.65124d-14,
!    $                         crh4=-1.42400d0

      real(kdef), allocatable, save :: crh1(:),  &
     &                                 crh2(:),  &
     &                                 crh3(:),  &
     &                                 crh4(:)

!$$$      real(kdef), parameter :: crh1= 2.07014d0,
!$$$     $                         crh2= 3.07900d0,
!$$$     $                         crh3= 5.18131d-17,
!$$$     $                         crh4=-1.42400d0

      end module modpar_umaer
