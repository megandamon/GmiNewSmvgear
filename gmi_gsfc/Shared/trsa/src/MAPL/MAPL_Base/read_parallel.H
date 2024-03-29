! $Id: read_parallel.H,v 1.2 2011-08-09 22:12:59 mrdamon Exp $

#ifdef NAME_
#undef NAME_
#endif

#define NAME_ READ_PARALLEL_

#include "overload.macro"

subroutine SUB_ ( layout, DATA, UNIT, FORMAT, RC)

  type (ESMF_DELayout)                         :: layout
#if (VARTYPE_ > 0)
  TYPE_(kind=EKIND_), intent(  OUT)            :: data DIMENSIONS_
#else
  character(LEN=*),   intent(  OUT)            :: data DIMENSIONS_
#endif
  integer,            intent(in   ),  optional :: UNIT
  character(len=*),   intent(in   ),  optional :: FORMAT
  integer         ,   intent(  out),  optional :: RC

  character(len=ESMF_MAXSTR) :: FORMATTED
  character(LEN=ESMF_MAXSTR) :: FILENAME
  logical                :: IS_NAMED
  integer                :: USABLE_UNIT
  integer                :: IOSTAT
  integer                :: status
  character(len=ESMF_MAXSTR) :: IAM='READ_PARALLEL'
#if (RANK_ == 1 && VARTYPE_ == 4)
  integer                               :: nretries
#endif

  if (present(unit)) then
     USABLE_UNIT = unit
  else
     USABLE_UNIT = 5 ! fortran stdin
  end if

  if(USABLE_UNIT < 0 ) then

#if (RANK_ > 2 || VARTYPE_ <= 0)

     ASSERT_(.FALSE.)

#else

     ASSERT_(-USABLE_UNIT<=LAST_UNIT)
     munit => MEM_units(-USABLE_UNIT)
     munit%prevrec = munit%prevrec + 1
#if (RANK_ > 0)
     ASSERT_(associated(munit%Records(munit%prevrec)%TKR_))
#endif
     data = munit%Records(munit%prevrec)%TKR_

#endif

  else

     if (MAPL_AM_I_ROOT(layout)) then
        if (present(UNIT)) then
           inquire(unit=USABLE_UNIT, formatted=FORMATTED)
#if (RANK_ == 1 && VARTYPE_ == 4)
           nretries = 0
10         continue
#endif
           if   (FORMATTED == "YES") then
              if (present(FORMAT)) then; read(USABLE_UNIT, FORMAT, IOSTAT=IOSTAT) DATA
              else;    read(USABLE_UNIT, *, IOSTAT=IOSTAT) DATA
              end if
           elseif(FORMATTED == "NO") then
	      read(USABLE_UNIT, IOSTAT=IOSTAT) DATA
           end if
#if (RANK_ == 1 && VARTYPE_ == 4)
           if (iostat /= 0) then
              nretries = nretries + 1
              print *,trim(IAM), ' read fails on attempt ',nretries
              if (nretries < 3) goto 10
           endif
#endif

           if (IOSTAT < 0) then
              inquire(unit=USABLE_UNIT, NAMED=IS_NAMED, NAME=FILENAME)
              if (.not. IS_NAMED) then
                 FILENAME = 'UNNAMED'
              end if
              print *, "Premature end of file ",FILENAME
              RETURN_(ESMF_FAILURE)
           end if
        else
           if (present(FORMAT)) then; read(*, FORMAT ) DATA
           else; read          *, DATA
           end if
        end if
     end if

#if (RANK_ == 0)
#if (VARTYPE_ == 0)
     call MAPL_CommsBcast(layout, data, len(data), root, status)
#else
     call MAPL_CommsBcast(layout, data, 1, root, status)
#endif
#else
     call MAPL_CommsBcast(layout, data, size(data), root, status)
#endif
     VERIFY_(status)

  end if

  RETURN_(ESMF_SUCCESS)
END SUBROUTINE SUB_

!---------------------------
#undef NAME_

#undef DIMENSIONS_
#undef RANK_
#undef VARTYPE_
