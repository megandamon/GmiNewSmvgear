!-------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiGrid_mod
!
! !INTERFACE:
!
  module GmiGrid_mod
!
! !USES:
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  :: Get_i1, Get_i2, Get_ju1, Get_jv1, Get_j2, Get_k1, Get_k2
  public  :: Get_i1_gl, Get_i2_gl, Get_ju1_gl, Get_jv1_gl, Get_j2_gl, Get_k1_gl, Get_k2_gl
  public  :: Get_ilo, Get_ihi, Get_julo, Get_jvlo, Get_jhi
  public  :: Get_ilo_gl, Get_ihi_gl, Get_julo_gl, Get_jvlo_gl, Get_jhi_gl
  public  :: Get_ilong, Get_ilat, Get_ivert, Get_itloop
  public  :: Get_j1p, Get_j2p, Get_gmi_nborder, Get_numSpecies

  public  :: Set_i1, Set_i2, Set_ju1, Set_jv1, Set_j2, Set_k1, Set_k2
  public  :: Set_i1_gl, Set_i2_gl, Set_ju1_gl, Set_jv1_gl, Set_j2_gl, Set_k1_gl, Set_k2_gl
  public  :: Set_ilo, Set_ihi, Set_julo, Set_jvlo, Set_jhi
  public  :: Set_ilo_gl, Set_ihi_gl, Set_julo_gl, Set_jvlo_gl, Set_jhi_gl
  public  :: Set_ilong, Set_ilat, Set_ivert, Set_itloop
  public  :: Set_j1p, Set_j2p, Set_gmi_nborder, Set_numSpecies

! !PUBLIC DATA MEMBERS:

  public  :: t_gmiGrid

  type t_gmiGrid
    private
    integer :: numSpecies  ! Total number of species

    integer :: i1_gl       ! index of first global longitude      (no ghost zones)
    integer :: i2_gl       ! index of last  global longitude      (no ghost zones)
    integer :: ju1_gl      ! index of first global "u"   latitude (no ghost zones)
    integer :: jv1_gl      ! index of first global "v"   latitude (no ghost zones)
    integer :: j2_gl       ! index of last  global "u&v" latitude (no ghost zones)
    integer :: k1_gl       ! index of first global altitude       (no ghost zones)
    integer :: k2_gl       ! index of last  global altitude       (no ghost zones)

    integer :: ilo_gl      ! i1_gl  - gmi_nborder (has ghost zones)
    integer :: ihi_gl      ! i2_gl  + gmi_nborder (has ghost zones)
    integer :: julo_gl     ! ju1_gl - gmi_nborder (has ghost zones)
    integer :: jvlo_gl     ! jv1_gl - gmi_nborder (has ghost zones)
    integer :: jhi_gl      ! j2_gl  + gmi_nborder (has ghost zones)

    integer :: i1          ! index of first local longitude      (no ghost zones)
    integer :: i2          ! index of last  local longitude      (no ghost zones)
    integer :: ju1         ! index of first local "u"   latitude (no ghost zones)
    integer :: jv1         ! index of first local "v"   latitude (no ghost zones)
    integer :: j2          ! index of last  local "u&v" latitude (no ghost zones)
    integer :: k1          ! index of first local altitude       (no ghost zones)
    integer :: k2          ! index of last  local altitude       (no ghost zones)

    integer :: ilo         ! i1  - gmi_nborder (has ghost zones)
    integer :: ihi         ! i2  + gmi_nborder (has ghost zones)
    integer :: julo        ! ju1 - gmi_nborder (has ghost zones)
    integer :: jvlo        ! jv1 - gmi_nborder (has ghost zones)
    integer :: jhi         ! j2  + gmi_nborder (has ghost zones)

    integer :: ilat        ! number of latitudes
    integer :: ilong       ! number of longitudes
    integer :: ivert       ! number of vertical layers
    integer :: itloop      ! number of zones (ilat * ilong * ivert)

    integer :: gmi_nborder ! number of longitude and latitude ghost zones

    integer :: j1p         ! determines size of the Polar cap
    integer :: j2p         ! j2_gl - j1p + 1
  end type t_gmiGrid

! !DESCRIPTION:
!  This module defined a derived type with member variables contain grid 
!  related information.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
! 5April2007 Initial implementation.
!
!EOP
!-------------------------------------------------------------------------

  CONTAINS

!-------------------------------------------------------------------------
  subroutine Get_numSpecies (self, numSpecies)
    implicit none
    integer         , intent(out)  :: numSpecies
    type (t_gmiGrid), intent(in)   :: self
    numSpecies = self%numSpecies
    return
  end subroutine Get_numSpecies
!------------------------------------------------------------------------
  subroutine Set_numSpecies (self, numSpecies)
    implicit none
    integer         , intent(in)  :: numSpecies
    type (t_gmiGrid), intent(inOut)   :: self
    self%numSpecies = numSpecies
    return
  end subroutine Set_numSpecies
!------------------------------------------------------------------------
  subroutine Get_i1 (self, i1)
    implicit none
    integer         , intent(out)  :: i1
    type (t_gmiGrid), intent(in)   :: self
    i1 = self%i1
    return
  end subroutine Get_i1
!------------------------------------------------------------------------
  subroutine Set_i1 (self, i1)
    implicit none
    integer         , intent(in)  :: i1
    type (t_gmiGrid), intent(inOut)   :: self
    self%i1 = i1
    return
  end subroutine Set_i1
!------------------------------------------------------------------------
  subroutine Get_i2 (self, i2)
    implicit none
    integer         , intent(out)  :: i2
    type (t_gmiGrid), intent(in)   :: self
    i2 = self%i2
    return
  end subroutine Get_i2
!------------------------------------------------------------------------
  subroutine Set_i2 (self, i2)
    implicit none
    integer         , intent(in)  :: i2
    type (t_gmiGrid), intent(inOut)   :: self
    self%i2 = i2
    return
  end subroutine Set_i2
!------------------------------------------------------------------------
  subroutine Get_ju1 (self, ju1)
    implicit none
    integer         , intent(out)  :: ju1
    type (t_gmiGrid), intent(in)   :: self
    ju1 = self%ju1
    return
  end subroutine Get_ju1
!------------------------------------------------------------------------
  subroutine Set_ju1 (self, ju1)
    implicit none
    integer         , intent(in)  :: ju1
    type (t_gmiGrid), intent(inOut)   :: self
    self%ju1 = ju1
    return
  end subroutine Set_ju1
!------------------------------------------------------------------------
  subroutine Get_jv1 (self, jv1)
    implicit none
    integer         , intent(out)  :: jv1
    type (t_gmiGrid), intent(in)   :: self
    jv1 = self%jv1
    return
  end subroutine Get_jv1
!------------------------------------------------------------------------
  subroutine Set_jv1 (self, jv1)
    implicit none
    integer         , intent(in)  :: jv1
    type (t_gmiGrid), intent(inOut)   :: self
    self%jv1 = jv1
    return
  end subroutine Set_jv1
!------------------------------------------------------------------------
  subroutine Get_j2 (self, j2)
    implicit none
    integer         , intent(out)  :: j2
    type (t_gmiGrid), intent(in)   :: self
    j2 = self%j2
    return
  end subroutine Get_j2
!------------------------------------------------------------------------
  subroutine Set_j2 (self, j2)
    implicit none
    integer         , intent(in)  :: j2
    type (t_gmiGrid), intent(inOut)   :: self
    self%j2 = j2
    return
  end subroutine Set_j2
!------------------------------------------------------------------------
  subroutine Get_k1 (self, k1)
    implicit none
    integer         , intent(out)  :: k1
    type (t_gmiGrid), intent(in)   :: self
    k1 = self%k1
    return
  end subroutine Get_k1
!------------------------------------------------------------------------
  subroutine Set_k1 (self, k1)
    implicit none
    integer         , intent(in)  :: k1
    type (t_gmiGrid), intent(inOut)   :: self
    self%k1 = k1
    return
  end subroutine Set_k1
!------------------------------------------------------------------------
  subroutine Get_k2 (self, k2)
    implicit none
    integer         , intent(out)  :: k2
    type (t_gmiGrid), intent(in)   :: self
    k2 = self%k2
    return
  end subroutine Get_k2
!------------------------------------------------------------------------
  subroutine Set_k2 (self, k2)
    implicit none
    integer         , intent(in)  :: k2
    type (t_gmiGrid), intent(inOut)   :: self
    self%k2 = k2
    return
  end subroutine Set_k2
!------------------------------------------------------------------------
  subroutine Get_i1_gl (self, i1_gl)
    implicit none
    integer         , intent(out)  :: i1_gl
    type (t_gmiGrid), intent(in)   :: self
    i1_gl = self%i1_gl
    return
  end subroutine Get_i1_gl
!------------------------------------------------------------------------
  subroutine Set_i1_gl (self, i1_gl)
    implicit none
    integer         , intent(in)  :: i1_gl
    type (t_gmiGrid), intent(inOut)   :: self
    self%i1_gl = i1_gl
    return
  end subroutine Set_i1_gl
!------------------------------------------------------------------------
  subroutine Get_i2_gl (self, i2_gl)
    implicit none
    integer         , intent(out)  :: i2_gl
    type (t_gmiGrid), intent(in)   :: self
    i2_gl = self%i2_gl
    return
  end subroutine Get_i2_gl
!------------------------------------------------------------------------
  subroutine Set_i2_gl (self, i2_gl)
    implicit none
    integer         , intent(in)  :: i2_gl
    type (t_gmiGrid), intent(inOut)   :: self
    self%i2_gl = i2_gl
    return
  end subroutine Set_i2_gl
!------------------------------------------------------------------------
  subroutine Get_ju1_gl (self, ju1_gl)
    implicit none
    integer         , intent(out)  :: ju1_gl
    type (t_gmiGrid), intent(in)   :: self
    ju1_gl = self%ju1_gl
    return
  end subroutine Get_ju1_gl
!------------------------------------------------------------------------
  subroutine Set_ju1_gl (self, ju1_gl)
    implicit none
    integer         , intent(in)  :: ju1_gl
    type (t_gmiGrid), intent(inOut)   :: self
    self%ju1_gl = ju1_gl
    return
  end subroutine Set_ju1_gl
!------------------------------------------------------------------------
  subroutine Get_jv1_gl (self, jv1_gl)
    implicit none
    integer         , intent(out)  :: jv1_gl
    type (t_gmiGrid), intent(in)   :: self
    jv1_gl = self%jv1_gl
    return
  end subroutine Get_jv1_gl
!------------------------------------------------------------------------
  subroutine Set_jv1_gl (self, jv1_gl)
    implicit none
    integer         , intent(in)  :: jv1_gl
    type (t_gmiGrid), intent(inOut)   :: self
    self%jv1_gl = jv1_gl
    return
  end subroutine Set_jv1_gl
!------------------------------------------------------------------------
  subroutine Get_j2_gl (self, j2_gl)
    implicit none
    integer         , intent(out)  :: j2_gl
    type (t_gmiGrid), intent(in)   :: self
    j2_gl = self%j2_gl
    return
  end subroutine Get_j2_gl
!------------------------------------------------------------------------
  subroutine Set_j2_gl (self, j2_gl)
    implicit none
    integer         , intent(in)  :: j2_gl
    type (t_gmiGrid), intent(inOut)   :: self
    self%j2_gl = j2_gl
    return
  end subroutine Set_j2_gl
!------------------------------------------------------------------------
  subroutine Get_k1_gl (self, k1_gl)
    implicit none
    integer         , intent(out)  :: k1_gl
    type (t_gmiGrid), intent(in)   :: self
    k1_gl = self%k1_gl
    return
  end subroutine Get_k1_gl
!------------------------------------------------------------------------
  subroutine Set_k1_gl (self, k1_gl)
    implicit none
    integer         , intent(in)  :: k1_gl
    type (t_gmiGrid), intent(inOut)   :: self
    self%k1_gl = k1_gl
    return
  end subroutine Set_k1_gl
!------------------------------------------------------------------------
  subroutine Get_k2_gl (self, k2_gl)
    implicit none
    integer         , intent(out)  :: k2_gl
    type (t_gmiGrid), intent(in)   :: self
    k2_gl = self%k2_gl
    return
  end subroutine Get_k2_gl
!------------------------------------------------------------------------
  subroutine Set_k2_gl (self, k2_gl)
    implicit none
    integer         , intent(in)  :: k2_gl
    type (t_gmiGrid), intent(inOut)   :: self
    self%k2_gl = k2_gl
    return
  end subroutine Set_k2_gl
!------------------------------------------------------------------------
  subroutine Get_ilo (self, ilo)
    implicit none
    integer         , intent(out)  :: ilo
    type (t_gmiGrid), intent(in)   :: self
    ilo = self%ilo
    return
  end subroutine Get_ilo
!------------------------------------------------------------------------
  subroutine Set_ilo (self, ilo)
    implicit none
    integer         , intent(in)  :: ilo
    type (t_gmiGrid), intent(inOut)   :: self
    self%ilo = ilo
    return
  end subroutine Set_ilo
!------------------------------------------------------------------------
  subroutine Get_ihi (self, ihi)
    implicit none
    integer         , intent(out)  :: ihi
    type (t_gmiGrid), intent(in)   :: self
    ihi = self%ihi
    return
  end subroutine Get_ihi
!------------------------------------------------------------------------
  subroutine Set_ihi (self, ihi)
    implicit none
    integer         , intent(in)  :: ihi
    type (t_gmiGrid), intent(inOut)   :: self
    self%ihi = ihi
    return
  end subroutine Set_ihi
!------------------------------------------------------------------------
  subroutine Get_julo (self, julo)
    implicit none
    integer         , intent(out)  :: julo
    type (t_gmiGrid), intent(in)   :: self
    julo = self%julo
    return
  end subroutine Get_julo
!------------------------------------------------------------------------
  subroutine Set_julo (self, julo)
    implicit none
    integer         , intent(in)  :: julo
    type (t_gmiGrid), intent(inOut)   :: self
    self%julo = julo
    return
  end subroutine Set_julo
!------------------------------------------------------------------------
  subroutine Get_jvlo (self, jvlo)
    implicit none
    integer         , intent(out)  :: jvlo
    type (t_gmiGrid), intent(in)   :: self
    jvlo = self%jvlo
    return
  end subroutine Get_jvlo
!------------------------------------------------------------------------
  subroutine Set_jvlo (self, jvlo)
    implicit none
    integer         , intent(in)  :: jvlo
    type (t_gmiGrid), intent(inOut)   :: self
    self%jvlo = jvlo
    return
  end subroutine Set_jvlo
!------------------------------------------------------------------------
  subroutine Get_jhi (self, jhi)
    implicit none
    integer         , intent(out)  :: jhi
    type (t_gmiGrid), intent(in)   :: self
    jhi = self%jhi
    return
  end subroutine Get_jhi
!------------------------------------------------------------------------
  subroutine Set_jhi (self, jhi)
    implicit none
    integer         , intent(in)  :: jhi
    type (t_gmiGrid), intent(inOut)   :: self
    self%jhi = jhi
    return
  end subroutine Set_jhi
!------------------------------------------------------------------------
  subroutine Get_ilo_gl (self, ilo_gl)
    implicit none
    integer         , intent(out)  :: ilo_gl
    type (t_gmiGrid), intent(in)   :: self
    ilo_gl = self%ilo_gl
    return
  end subroutine Get_ilo_gl
!------------------------------------------------------------------------
  subroutine Set_ilo_gl (self, ilo_gl)
    implicit none
    integer         , intent(in)  :: ilo_gl
    type (t_gmiGrid), intent(inOut)   :: self
    self%ilo_gl = ilo_gl
    return
  end subroutine Set_ilo_gl
!------------------------------------------------------------------------
  subroutine Get_ihi_gl (self, ihi_gl)
    implicit none
    integer         , intent(out)  :: ihi_gl
    type (t_gmiGrid), intent(in)   :: self
    ihi_gl = self%ihi_gl
    return
  end subroutine Get_ihi_gl
!------------------------------------------------------------------------
  subroutine Set_ihi_gl (self, ihi_gl)
    implicit none
    integer         , intent(in)  :: ihi_gl
    type (t_gmiGrid), intent(inOut)   :: self
    self%ihi_gl = ihi_gl
    return
  end subroutine Set_ihi_gl
!------------------------------------------------------------------------
  subroutine Get_julo_gl (self, julo_gl)
    implicit none
    integer         , intent(out)  :: julo_gl
    type (t_gmiGrid), intent(in)   :: self
    julo_gl = self%julo_gl
    return
  end subroutine Get_julo_gl
!------------------------------------------------------------------------
  subroutine Set_julo_gl (self, julo_gl)
    implicit none
    integer         , intent(in)  :: julo_gl
    type (t_gmiGrid), intent(inOut)   :: self
    self%julo_gl = julo_gl
    return
  end subroutine Set_julo_gl
!------------------------------------------------------------------------
  subroutine Get_jvlo_gl (self, jvlo_gl)
    implicit none
    integer         , intent(out)  :: jvlo_gl
    type (t_gmiGrid), intent(in)   :: self
    jvlo_gl = self%jvlo_gl
    return
  end subroutine Get_jvlo_gl
!------------------------------------------------------------------------
  subroutine Set_jvlo_gl (self, jvlo_gl)
    implicit none
    integer         , intent(in)  :: jvlo_gl
    type (t_gmiGrid), intent(inOut)   :: self
    self%jvlo_gl = jvlo_gl
    return
  end subroutine Set_jvlo_gl
!------------------------------------------------------------------------
  subroutine Get_jhi_gl (self, jhi_gl)
    implicit none
    integer         , intent(out)  :: jhi_gl
    type (t_gmiGrid), intent(in)   :: self
    jhi_gl = self%jhi_gl
    return
  end subroutine Get_jhi_gl
!------------------------------------------------------------------------
  subroutine Set_jhi_gl (self, jhi_gl)
    implicit none
    integer         , intent(in)  :: jhi_gl
    type (t_gmiGrid), intent(inOut)   :: self
    self%jhi_gl = jhi_gl
    return
  end subroutine Set_jhi_gl
!------------------------------------------------------------------------
  subroutine Get_ilong (self, ilong)
    implicit none
    integer         , intent(out)  :: ilong
    type (t_gmiGrid), intent(in)   :: self
    ilong = self%ilong
    return
  end subroutine Get_ilong
!------------------------------------------------------------------------
  subroutine Set_ilong (self, ilong)
    implicit none
    integer         , intent(in)  :: ilong
    type (t_gmiGrid), intent(inOut)   :: self
    self%ilong = ilong
    return
  end subroutine Set_ilong
!------------------------------------------------------------------------
  subroutine Get_ilat (self, ilat)
    implicit none
    integer         , intent(out)  :: ilat
    type (t_gmiGrid), intent(in)   :: self
    ilat = self%ilat
    return
  end subroutine Get_ilat
!------------------------------------------------------------------------
  subroutine Set_ilat (self, ilat)
    implicit none
    integer         , intent(in)  :: ilat
    type (t_gmiGrid), intent(inOut)   :: self
    self%ilat = ilat
    return
  end subroutine Set_ilat
!------------------------------------------------------------------------
  subroutine Get_ivert (self, ivert)
    implicit none
    integer         , intent(out)  :: ivert
    type (t_gmiGrid), intent(in)   :: self
    ivert = self%ivert
    return
  end subroutine Get_ivert
!------------------------------------------------------------------------
  subroutine Set_ivert (self, ivert)
    implicit none
    integer         , intent(in)  :: ivert
    type (t_gmiGrid), intent(inOut)   :: self
    self%ivert = ivert
    return
  end subroutine Set_ivert
!------------------------------------------------------------------------
  subroutine Get_itloop (self, itloop)
    implicit none
    integer         , intent(out)  :: itloop
    type (t_gmiGrid), intent(in)   :: self
    itloop = self%itloop
    return
  end subroutine Get_itloop
!------------------------------------------------------------------------
  subroutine Set_itloop (self, itloop)
    implicit none
    integer         , intent(in)  :: itloop
    type (t_gmiGrid), intent(inOut)   :: self
    self%itloop = itloop
    return
  end subroutine Set_itloop
!------------------------------------------------------------------------
  subroutine Get_j1p (self, j1p)
    implicit none
    integer         , intent(out)  :: j1p
    type (t_gmiGrid), intent(in)   :: self
    j1p = self%j1p
    return
  end subroutine Get_j1p
!------------------------------------------------------------------------
  subroutine Set_j1p (self, j1p)
    implicit none
    integer         , intent(in)  :: j1p
    type (t_gmiGrid), intent(inOut)   :: self
    self%j1p = j1p
    return
  end subroutine Set_j1p
!------------------------------------------------------------------------
  subroutine Get_j2p (self, j2p)
    implicit none
    integer         , intent(out)  :: j2p
    type (t_gmiGrid), intent(in)   :: self
    j2p = self%j2p
    return
  end subroutine Get_j2p
!------------------------------------------------------------------------
  subroutine Set_j2p (self, j2p)
    implicit none
    integer         , intent(in)  :: j2p
    type (t_gmiGrid), intent(inOut)   :: self
    self%j2p = j2p
    return
  end subroutine Set_j2p
!------------------------------------------------------------------------
  subroutine Get_gmi_nborder (self, gmi_nborder)
    implicit none
    integer         , intent(out)  :: gmi_nborder
    type (t_gmiGrid), intent(in)   :: self
    gmi_nborder = self%gmi_nborder
    return
  end subroutine Get_gmi_nborder
!------------------------------------------------------------------------
  subroutine Set_gmi_nborder (self, gmi_nborder)
    implicit none
    integer         , intent(in)  :: gmi_nborder
    type (t_gmiGrid), intent(inOut)   :: self
    self%gmi_nborder = gmi_nborder
    return
  end subroutine Set_gmi_nborder
!------------------------------------------------------------------------


  end module GmiGrid_mod
