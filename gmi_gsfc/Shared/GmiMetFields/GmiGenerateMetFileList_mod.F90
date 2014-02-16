!-----------------------------------------------------------------------------
! NASA GSFC - SSSO Code 610.3
!-----------------------------------------------------------------------------
!BOP
!
! !IMODULE: GmiGenerateMetFileList_mod
!
! !INTERFACE:
!
module GmiGenerateMetFileList_mod
!
! !USES:
      USE GmiTimeControl_mod, only : incrementDate

      implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
      PRIVATE
      PUBLIC  :: createMetFileList
!
#     include "GmiParameters.h"
#     include "gmi_time_constants.h"

      INTEGER, PARAMETER :: stderr = 6
!
! !DESCRIPTION:
! Provides supporting routines to be able to automatically generate the metFields
! file list from a file name template.
!
!EOP
!------------------------------------------------------------------------------
      CONTAINS
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: createMetFileList
!
! !INTERFACE:
!
      SUBROUTINE createMetFileList(metFileList, metFileTemplate, begDate, &
                       timeStamp, lenRun, doMEGANemiss, numListItems)
!
! !INPUT PARAMETERS:
      INTEGER, intent(in) :: lenRun ! length of the experiment (days)
      INTEGER, intent(in) :: begDate ! reference date (YYYYMMDD) for the first metFields file
      INTEGER, intent(in) :: timeStamp ! reference time when each metField file was produced
      LOGICAL, intent(in) :: doMEGANemiss ! do MEGAN emissions?
      CHARACTER(len=*), intent(in) :: metFileTemplate
!
! !INPUT/OUTPUT PARAMETERS:
      INTEGER, INTENT(out) :: numListItems ! number of files in the list
      CHARACTER(len=*), intent(out) :: metFileList(MAX_INFILE_NAMES)
!
! !DESCRIPTION:
! Given 
! \bi
! \item a file name template
! \item the metFields beginning date
! \item a reference time
! \item the length of the experiments
! \item whether MEGAN emissions is done or not
! \ei
! this subroutine automatically generates
!
! !LOCAL VARIABLES:
      INTEGER :: i, stat
      INTEGER :: currentDay, nextDay, previous15Days(15)
!EOP
!------------------------------------------------------------------------------
!BOC

      IF (doMEGANemiss) THEN
         numListItems = 15 + lenRun + 1

         ! Get the previous 15 days
         previous15Days = getPrevious15Days(begDate)

         DO i = 15, 1, -1
            call getFileNameFromTemplate(metFileList(15-i+1), metFileTemplate, &
                      nymd=previous15Days(i), &
                      nhms=timeStamp, &
                      stat = stat)
         ENDDO

         ! Files for the model run
         currentDay = begDate
         DO i = 1, lenRun+1
            call getFileNameFromTemplate(metFileList(i+15), metFileTemplate, &
                      nymd=currentDay, &
                      nhms=timeStamp, &
                      stat = stat)
            nextDay = getNextDay(currentDay)
            currentDay = nextDay
         END DO
         
      ELSE
         numListItems = lenRun + 1
         ! Files for the model run
         currentDay = begDate
         DO i = 1, lenRun+1
            call getFileNameFromTemplate(metFileList(i), metFileTemplate, &
                      nymd=currentDay, &
                      nhms=timeStamp, &
                      stat = stat)
            nextDay = getNextDay(currentDay)
            currentDay = nextDay
         END DO
         
      ENDIF

      RETURN

      END SUBROUTINE createMetFileList
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getNextDay
!
! !INTERFACE:
!
      FUNCTION getNextDay(curDay) result(nextDay)
!
! !INPUT PARAMETERS:
      INTEGER, intent(in) :: curDay ! current day (YYYYMMDD)
!
! !RETURNED VALUE:
      INTEGER :: nextDay ! next day (YYYYMMDD)
!
! !DESCRIPTION:
! Given the current day, returns the next day.
!EOC
!------------------------------------------------------------------------------
!BOC
      nextDay = incrementDate(curDay, 1)

      END FUNCTION getNextDay
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: update15DayList
!
! !INTERFACE:
!
      SUBROUTINE update15DayList(previous15Days) 
!
! !INPUT/OUTPUT PARAMETERS:
      INTEGER, intent(inOut) :: previous15Days(15) ! (YYYYMMDD)
!
! !DESCRIPTION:
! Updates the list of the previous 15 days.
!
      INTEGER :: id, curDay
!EOC
!------------------------------------------------------------------------------
!BOC
        curDay = previous15Days(1) 
        DO id = 15, 2, -1
           previous15Days(id) = previous15Days(id-1)
        END DO
        previous15Days(1) = incrementDate(curDay, 1)

        END SUBROUTINE update15DayList
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: getPrevious15Days
!
! !INTERFACE:
!
      FUNCTION getPrevious15Days(curDay) result(previous15Days)
!
! !INPUT PARAMETERS:
      INTEGER, intent(in) :: curDay ! current day (YYYYMMDD)
!
! !RETURNED VALUE:
      INTEGER :: previous15Days(15) !  (YYYYMMDD)
!
! !DESCRIPTION:
! Given the current day, returns the previous 15 days.
!
      INTEGER :: id
!EOC
!------------------------------------------------------------------------------
!BOC
        previous15Days(1) = incrementDate(curDay, -1)
        DO id = 2, 15
           previous15Days(id) = incrementDate(previous15Days(id-1), -1)
        END DO

        END FUNCTION getPrevious15Days
!EOC
!------------------------------------------------------------------------------
!BOP 
!
! !IROUTINE: getFileNameFromTemplate - evaluate a GrADS style string template
!
! !DESCRIPTION:
!
! !INTERFACE:
!
      subroutine getFileNameFromTemplate(str,tmpl,xid,nymd,nhms,stat)
!
      character(len=*),intent(out) :: str
      character(len=*),intent(in ) :: tmpl
      character(len=*),optional,intent(in) :: xid
      integer,optional,intent(in)  :: nymd
      integer,optional,intent(in)  :: nhms
      integer,optional,intent(out) :: stat

!EOP
!-------------------------------------------------------------------------

      character(len=*),parameter :: myname_='::getFileNameFromTemplate'

      integer :: iy4,iy2,imo,idy
      integer :: ihr,imn
      integer :: i,i1,i2,m,k
      integer :: ln_tmpl,ln_str
      integer :: istp,kstp
      integer :: ier

      character(len=1) :: c0,c1,c2
      character(len=4) :: sbuf
!________________________________________
        ! Determine iyr, imo, and idy
  iy4=-1
  iy2=-1
  imo=-1
  idy=-1
  if(present(nymd)) then
        if(nymd <= 0) then
          call perr(myname_,'nymd <= 0',nymd)
          if(.not.present(stat)) call die(myname_)
          stat=1
          return
        endif

    i=nymd
    iy4=i/10000
    iy2=mod(iy4,100)
      i=mod(i,10000)
    imo=i/100
      i=mod(i,100)
    idy=i
  endif
!________________________________________
        ! Determine ihr and imn
  ihr=-1
  imn=-1
  if(present(nhms)) then
        if(nhms < 0) then
          call perr(myname_,'nhms < 0',nhms)
          if(.not.present(stat)) call die(myname_)
          stat=1
          return
        endif

    i=nhms
    ihr=i/10000
      i=mod(i,10000)
    imn=i/100
  endif
!________________________________________

  ln_tmpl=len_trim(tmpl)        ! size of the format template
  ln_str =len(str)              ! size of the output string
!________________________________________

  if(present(stat)) stat=0

str=""

i=0; istp=1
k=1; kstp=1

do while( i+istp <= ln_tmpl )   ! A loop over all tokens in (tmpl)

  if(k>ln_Str) exit     ! truncate the output here.

  i=i+istp
  c0=tmpl(i:i)

  select case(c0)
  case ("$")
    call genv_(tmpl,ln_tmpl,i,istp,str,ln_str,k,ier)
    if(ier/=0) then
      call perr(myname_,'genv_("'//tmpl(i:ln_tmpl)//'"',ier)
      if(.not.present(stat)) call die(myname_)
      stat=1
      return
    endif

  case ("%")
        !________________________________________

    c1=""
    i1=i+1
    if(i1 <= ln_Tmpl) c1=tmpl(i1:i1)
        !________________________________________

    select case(c1)

    case("s")
      if(.not.present(xid)) then
        write(stderr,'(2a)') myname_,   &
                ': optional argument expected, "xid="'
        if(.not.present(stat)) call die(myname_)
        stat=1
        return
      endif

      istp=2
      m=min(k+len_trim(xid)-1,ln_str)
      str(k:m)=xid
      k=m+1
      cycle

    case("%","$")

      istp=2
      str(k:k)=c1
      k=k+1     ! kstp=1
      cycle

    case default

      c2=""
      i2=i+2
      if(i2 <= ln_Tmpl) c2=tmpl(i2:i2)
        !________________________________________

      select case(c1//c2)

      case("y4","y2","m1","m2","mc","Mc","MC","d1","d2")
        if(.not.present(nymd)) then
          write(stderr,'(2a)') myname_, &
                ': optional argument expected, "nymd="'
          if(.not.present(stat)) call die(myname_)
          stat=1
          return
        endif
        istp=3

      case("h1","h2","h3","n2")
        if(.not.present(nhms)) then
          write(stderr,'(2a)') myname_, &
                ': optional argument expected, "nhms="'
          if(.not.present(stat)) call die(myname_)
          stat=1
          return
        endif
        istp=3

      case default

        write(stderr,'(4a)') myname_,   &
          ': invalid template entry, "',trim(tmpl(i:)),'"'
        if(.not.present(stat)) call die(myname_)
        stat=2
        return

      end select          ! case(c1//c2)
    end select          ! case(c1)
        !________________________________________

    select case(c1)

    case("y")
      select case(c2)
      case("2")
        write(sbuf,'(i2.2)') iy2
        kstp=2
      case("4")
        write(sbuf,'(i4.4)') iy4
        kstp=4
      case default
        write(stderr,'(4a)') myname_,   &
          ': invalid template entry, "',trim(tmpl(i:)),'"'
        if(.not.present(stat)) call die(myname_)
        stat=2
        return
      end select

    case("m")
      select case(c2)
      case("1")
        if(imo < 10) then
          write(sbuf,'(i1)') imo
          kstp=1
        else
          write(sbuf,'(i2)') imo
          kstp=2
        endif
      case("2")
        write(sbuf,'(i2.2)') imo
        kstp=2
      case("c")
        sbuf=mon_lc(imo)
        kstp=3
      case default
        write(stderr,'(4a)') myname_,   &
          ': invalid template entry, "',trim(tmpl(i:)),'"'
        if(.not.present(stat)) call die(myname_)
        stat=2
        return
      end select

    case("M")
      select case(c2)
      case("c")
        sbuf=mon_wd(imo)
        kstp=3
      case("C")
        sbuf=mon_uc(imo)
        kstp=3
      case default
        write(stderr,'(4a)') myname_,   &
          ': invalid template entry, "',trim(tmpl(i:)),'"'
        if(.not.present(stat)) call die(myname_)
        stat=2
        return
      end select

    case("d")
      select case(c2)
      case("1")
        if(idy < 10) then
          write(sbuf,'(i1)') idy
          kstp=1
        else
          write(sbuf,'(i2)') idy
          kstp=2
        endif
      case("2")
        write(sbuf,'(i2.2)') idy
        kstp=2
      case default
        write(stderr,'(4a)') myname_,   &
          ': invalid template entry, "',trim(tmpl(i:)),'"'
        if(.not.present(stat)) call die(myname_)
        stat=2
        return
      end select

    case("h")
      select case(c2)
      case("1")
        if(ihr < 10) then
          write(sbuf,'(i1)') ihr
          kstp=1
        else
          write(sbuf,'(i2)') ihr
          kstp=2
        endif
      case("2")
        write(sbuf,'(i2.2)') ihr
        kstp=2
      case("3")
        write(sbuf,'(i3.3)') ihr
        kstp=3
      case default
        write(stderr,'(4a)') myname_,   &
          ': invalid template entry, "',trim(tmpl(i:)),'"'
        if(.not.present(stat)) call die(myname_)
        stat=2
        return
      end select

    case("n")
      select case(c2)
      case("2")
        write(sbuf,'(i2.2)') imn
        kstp=2
      case default
        write(stderr,'(4a)') myname_,   &
          ': invalid template entry, "',trim(tmpl(i:)),'"'
        if(.not.present(stat)) call die(myname_)
        stat=2
        return
      end select

    case default
        write(stderr,'(4a)') myname_,   &
          ': invalid template entry, "',trim(tmpl(i:)),'"'
        if(.not.present(stat)) call die(myname_)
        stat=2
      return
    end select  ! case(c1)

    m=min(k+kstp-1,ln_Str)
    str(k:m)=sbuf
    k=m+1

  case default

    istp=1
    str(k:k)=tmpl(i:i)
    k=k+1

  end select    ! case(c0)
end do

contains
      subroutine genv_(tmpl,lnt,i,istp,str,lns,k,ier)
      character(len=*),intent(in) :: tmpl
      integer,intent(in)  :: lnt
      integer,intent(in)  :: i
      integer,intent(out) :: istp
      character(len=*),intent(inout) :: str
      integer         ,intent(in)    :: lns
      integer         ,intent(inout) :: k
      integer,intent(out) :: ier

      integer :: j,jb,je
      integer :: l,m
      logical :: bracket,more
      character(len=256) :: env

      j=i+1         ! skip "$"
      ier=0

      if (j>lnt) then
         ier=1
         return
      endif

      bracket = tmpl(j:j)=='{'
      if (bracket) j=j+1

        ! There is at least one a letter (including "_") to start a
        ! variable name

      select case(tmpl(j:j))
      case ("A":"Z","a":"z","_")
      case default
        ier=2
        return
      end select

      jb=j
      je=j

      if (bracket) then

         more=.true.
         do while(more)
            select case(tmpl(j:j))
            case ("A":"Z","a":"z","_","0":"9")
              je=j
              j=j+1
            case ("}")        ! End if "}" or eos
              j=j+1
              exit
            case default
              ier=3
              return
            end select
            more=j<=lnt
         enddo

      else

         more=.true.
         do while(more)
            select case(tmpl(j:j))
            case ("A":"Z","a":"z","_","0":"9")
              je=j
              j=j+1
            case default
              exit
            end select
            more=j<=lnt
         enddo
      endif
    
      istp=j-i

      call getenv(tmpl(jb:je),env)
      l=len_trim(env)
      m=min(k+l-1,lns)
      str(k:m)=env
      k=m+1

      end subroutine genv_
!EOC
!------------------------------------------------------------------------------
      end subroutine getFileNameFromTemplate

      subroutine perr ( name, msg, i )
      character(len=*) :: name, msg
      integer :: i
      print *, trim(name)//':'//trim(msg), i
      return
      end subroutine perr

      subroutine die ( name )
      character(len=*) :: name
      print *, trim(name)//': fatal error in CFIO, but not aborting...'
      return
      end subroutine die

end module GmiGenerateMetFileList_mod
