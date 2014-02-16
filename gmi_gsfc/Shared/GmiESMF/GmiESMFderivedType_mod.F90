

      module GmiESMFderivedType_mod

      ! ESMF module, defines all ESMF data types and procedures
      use ESMF_Mod

      implicit none

! !PUBLIC DATA MEMBERS:
      private
      public  :: t_gmiESMF
      
      ! ESMF GMI Component derived type
      TYPE t_gmiESMF
           type(ESMF_State   ) :: stateImp
           type(ESMF_State   ) :: stateExp
           type(ESMF_GridComp) :: compGridded
      end TYPE t_gmiESMF

      end module GmiESMFderivedType_mod
