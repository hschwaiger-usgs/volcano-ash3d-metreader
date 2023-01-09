!##############################################################################
!##############################################################################
!
! MetSonde
!
! This stand-alone program calculates vertical temperature profile from the GFS
! and NCEP 50-year reanalysis products.  This program assumes the GFS and NCEP
! data are stored in the default locations using the scripts in autorun_scripts,
! i.e. in /data/WindFiles/
! There are 6 required command-line arguments:
!       lon        : real    : longitude of sonde point
!       lat        : real    : latitude of sonde point
!       YYYY       : integer : year of sonde
!       MM         : integer : month of sonde
!       DD         : integer : day of sonde
!       H.H        : real    : hour of sonde
!       [WINDROOT] : char    : path of windroot, if different from default
!
! MetSonde -169.9468 52.8217 2022 8 29 5.5
!
! Note: This program only works interpolating on to a lat/lon point and
! interpolating to a particular time for GFS and NCEP data.  If you want to
! probe other types of wind files or individual files (without time
! interpolation), you can use the program probe_Met, also in this directory.
!
!##############################################################################

      program MetSonde

      use MetReader

      implicit none

      ! These are the variables that must be set in the input file
      real(kind=4)        :: inlon
      real(kind=4)        :: inlat
      integer             :: inyear
      integer             :: inmonth
      integer             :: inday
      real(kind=8)        :: inhour
      character(len=100)  :: WINDROOT

      integer             :: FC_freq = 12
      integer             :: GFS_Archive_Days = 14
      integer             :: GFS_FC_TotHours = 198  ! This is the max anticipated
      integer             :: nxmax,nymax,nzmax
      real(kind=4),dimension(:)    ,allocatable :: lon_grid
      real(kind=4),dimension(:)    ,allocatable :: lat_grid
      real(kind=4),dimension(:)    ,allocatable :: z_cc
      logical             :: IsPeriodic

      integer :: BaseYear = 1900
      logical :: useLeap  = .true.

      INTERFACE
        subroutine Read_ComdLine(inlon,inlat, &
                      inyear,inmonth,inday,inhour,&
                      WINDROOT)
          integer,parameter   :: sp        = 4 ! single precision
          integer,parameter   :: dp        = 8 ! double precision
          real(kind=sp)      ,intent(out) :: inlon
          real(kind=sp)      ,intent(out) :: inlat
          integer            ,intent(out) :: inyear
          integer            ,intent(out) :: inmonth
          integer            ,intent(out) :: inday
          real(kind=dp)      ,intent(out) :: inhour
          character(len=100) ,intent(out) :: WINDROOT
        end subroutine
        subroutine GetWindFile(inyear,inmonth,inday,inhour,WINDROOT,FC_freq,GFS_Archive_Days,GFS_FC_TotHours)
          integer,parameter   :: dp        = 8 ! double precision
          integer            ,intent(in) :: inyear
          integer            ,intent(in) :: inmonth
          integer            ,intent(in) :: inday
          real(kind=dp)      ,intent(in) :: inhour
          character(len=100) ,intent(in) :: WINDROOT
          integer            ,intent(in) :: FC_freq
          integer            ,intent(in) :: GFS_Archive_Days
          integer            ,intent(in) :: GFS_FC_TotHours
        end subroutine
        subroutine GetMetProfile(inlon,inlat,inyear,inmonth,inday,inhour)
          integer,parameter   :: sp        = 4 ! single precision
          integer,parameter   :: dp        = 8 ! double precision
          real(kind=sp), intent(in)  :: inlon
          real(kind=sp), intent(in)  :: inlat
          integer      , intent(in)  :: inyear
          integer      , intent(in)  :: inmonth
          integer      , intent(in)  :: inday
          real(kind=dp), intent(in)  :: inhour
        end subroutine
      END INTERFACE

      call Read_ComdLine(inlon,inlat, &
                      inyear,inmonth,inday,inhour,&
                      WINDROOT)

      ! Now set up the computational grid
      nxmax = 3 ! 
      nymax = 3 ! 
      nzmax = 2 ! This is not really used in this utility
      allocate(lon_grid(nxmax)); lon_grid(1:3) = (/inlon-0.5_4,inlon,inlon+0.5_4/)
      allocate(lat_grid(nymax)); lat_grid(1:3) = (/inlat-0.5_4,inlat,inlat+0.5_4/)
      allocate(z_cc(nzmax))    ; z_cc(1:2) = (/0.0_4, 10.0_4/)
      IsPeriodic = .false.

      ! Make sure user MetReader is using the same calendar
      MR_BaseYear = BaseYear
      MR_useLeap  = useLeap

      write(MR_global_production,*)"Set up winfile data structure"
      call GetWindFile(inyear,inmonth,inday,inhour,WINDROOT,FC_freq,GFS_Archive_Days,GFS_FC_TotHours)

      write(MR_global_production,*)"Setting up wind grids"
      call MR_Initialize_Met_Grids(nxmax,nymax,nzmax,&
                              lon_grid(1:nxmax), &
                              lat_grid(1:nymax), &
                              z_cc(1:nzmax)    , &
                              IsPeriodic)

      write(MR_global_production,*)"Interpolating profile onto ",inlon,inlat
      call GetMetProfile(inlon,inlat,inyear,inmonth,inday,inhour)

      write(MR_global_production,*)"Program ended normally."

      end program MetSonde

!##############################################################################
!
!     Subroutines
!
!##############################################################################

!##############################################################################
! Read_ComdLine
!
!  This subroutine will parse the command-line options and set up the run
!  specifications. 
!
!##############################################################################

      subroutine Read_ComdLine(inlon,inlat, &
                    inyear,inmonth,inday,inhour,&
                    WINDROOT)

      use MetReader
      use projection

      implicit none

      ! These are the variables that must be set in the command line

      real(kind=4)       ,intent(out) :: inlon
      real(kind=4)       ,intent(out) :: inlat
      integer            ,intent(out) :: inyear
      integer            ,intent(out) :: inmonth
      integer            ,intent(out) :: inday
      real(kind=8)       ,intent(out) :: inhour
      character(len=100) ,intent(out) :: WINDROOT

      integer             :: nargs
      integer             :: status
      character (len=100) :: arg

      integer :: iprojflag
      real(kind=8) :: lambda0,phi0,phi1,phi2,k0,radius_earth
      logical :: IsLatLon

!     TEST READ COMMAND LINE ARGUMENTS
      nargs = command_argument_count()
      if (nargs.lt.6) then
        ! We need at least 6 parameter to define the run.
        write(MR_global_info,*)"Too few command-line arguments:"
        write(MR_global_info,*)"  Usage: MetSonde lon lat YYYY MM DD HH.H [WIND_ROOT]"
        write(MR_global_info,*)"           lon       = longitude of start point"
        write(MR_global_info,*)"           lat       = latitude of start point"
        write(MR_global_info,*)"           YYYY      = start year"
        write(MR_global_info,*)"           MM        = start month"
        write(MR_global_info,*)"           DD        = start day"
        write(MR_global_info,*)"           HH.H      = start hour"
        write(MR_global_info,*)"           WIND_ROOT = path to windfiles (optional)"
        stop 1
      else
        call get_command_argument(1, arg, status)
        read(arg,*)inlon
        if(inlon.lt.-360.0)then
          write(MR_global_error,*)"ERROR: Longitude must be gt -360"
          stop 1
        endif
        if(inlon.lt.0.0_4.or.inlon.gt.360.0_4)inlon=mod(inlon+360.0_4,360.0_4)
        call get_command_argument(2, arg, status)
        read(arg,*)inlat
        call get_command_argument(3, arg, status)
        read(arg,*)inyear
        call get_command_argument(4, arg, status)
        read(arg,*)inmonth
        call get_command_argument(5, arg, status)
        read(arg,*)inday
        call get_command_argument(6, arg, status)
        read(arg,*)inhour
        if(nargs.ge.7)then
          call get_command_argument(7, arg, status)
          WINDROOT = TRIM(arg)
        else
          WINDROOT='/data/WindFiles'
        endif
      endif

      IsLatLon = .true.
        ! These are all dummy values since the comp grid is lon/lat
      iprojflag = 1
      lambda0      = -105.0
      phi0         = 90.0
      phi1         = 90.0
      phi2         = 90.0
      k0           = 0.933
      radius_earth = 6371.229
      call MR_Set_CompProjection(IsLatLon,iprojflag,lambda0,phi0,phi1,phi2,&
                                 k0,radius_earth)

      ! write out values of parameters defining the run
      write(MR_global_info,*)"inlon              = ",real(inlon,kind=4)
      write(MR_global_info,*)"inlat              = ",real(inlat,kind=4)
      write(MR_global_info,*)"inyear             = ",inyear
      write(MR_global_info,*)"inmonth            = ",inmonth
      write(MR_global_info,*)"inday              = ",inday
      write(MR_global_info,*)"inhour             = ",real(inhour,kind=4)
      !write(MR_global_info,*)"IsLatLon           = ",IsLatLon
      !write(MR_global_info,*)"iw                 = ",iw
      !write(MR_global_info,*)"iwf                = ",iwf
      !write(MR_global_info,*)"igrid              = ",igrid
      !write(MR_global_info,*)"idf                = ",idf
      !write(MR_global_info,*)"FC_freq            = ",FC_freq
      !write(MR_global_info,*)"GFS_Archive_Days   = ",GFS_Archive_Days
      !write(MR_global_info,*)"iwfiles            = ",iwfiles
      write(MR_global_info,*)"--------------------------------------------------------------"

      end subroutine

!##############################################################################
!##############################################################################
!
!  GetWindFile
!
!  This subroutine sets the list of windfiles to be used in the calculation.
!  These will be through an assessment of the current GFS and NCEP files on
!  the system.
!
!##############################################################################

      subroutine GetWindFile(inyear,inmonth,inday,inhour,&
                             WINDROOT, &
                             FC_freq,GFS_Archive_Days,GFS_FC_TotHours)

      use MetReader

      implicit none

      integer            ,intent(in) :: inyear
      integer            ,intent(in) :: inmonth
      integer            ,intent(in) :: inday
      real(kind=8)       ,intent(in) :: inhour
      character(len=100) ,intent(in) :: WINDROOT
      integer            ,intent(in) :: FC_freq
      integer            ,intent(in) :: GFS_Archive_Days
      integer            ,intent(in) :: GFS_FC_TotHours

      !integer            :: HS_YearOfEvent
      !integer            :: HS_MonthOfEvent
      !integer            :: HS_DayOfEvent
      !real(kind=8)       :: HS_HourOfDay

      character(len=8)   :: date
      character(LEN=10)  :: time2         ! time argument used to get current
                                          !  date and time.
      character(len=5)   :: zone          ! variables used by the date_and_time subroutine
      integer            :: values(8)     ! return values from date_and_time
      integer            :: timezone      ! timezone of grid relative to UTC

      real(kind=8)      :: StartHour
      real(kind=8)      :: RunStartHour    ! Current UTC time, in hours since MR_BaseYear
      character(len=17) :: RunStartHour_ch
      real(kind=8)      :: Probe_StartHour

      integer      :: RunStartYear
      integer      :: RunStartMonth
      integer      :: RunStartDay
      integer      :: RunStartHr

      integer      :: FC_Package_year
      integer      :: FC_Package_month
      integer      :: FC_Package_day
      integer      :: FC_Package_hour
      integer      :: FC_hour_int
      real(kind=8) :: FC_Package_StartHour
      real(kind=8) :: FC_Archive_StartHour

      integer      :: iw,iwf,igrid,iwfiles
      real(kind=8)      :: Simtime_in_hours = 0.0
      character(len=47) :: string1,string2

      integer :: i,ii
      integer :: FC_year,FC_mon,FC_day
      real(kind=8) :: FC_hour,FC_Package_hour_dp,FC_intvl
      integer :: NumFCpackages

      character (len=130):: testfile
      real(kind=8)       :: FCStartHour
      real(kind=8)       :: FCEndHour
      logical            :: IsThere

      logical,dimension(:),allocatable :: GFS_candidate
      integer,dimension(:),allocatable :: GFS_FC_step_avail
      integer :: OptimalPackageNum

      INTERFACE
        character (len=13) function HS_yyyymmddhhmm_since(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_yyyymmddhhmm_since
        integer function HS_YearOfEvent(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_YearOfEvent
        integer function HS_MonthOfEvent(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_MonthOfEvent
        integer function HS_DayOfEvent(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_DayOfEvent
        real(kind=8) function HS_HourOfDay(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_HourOfDay
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer     ,intent(in) :: iyear
          integer     ,intent(in) :: imonth
          integer     ,intent(in) :: iday
          real(kind=8),intent(in) :: hours
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_hours_since_baseyear
      END INTERFACE

       ! Get the UTC time that the program is called
       !   This will be used to determine if gfs or NCEP winds are to be used
      call date_and_time(date,time2,zone,values)
      read(zone,'(i3)') timezone
        ! FIND TIME IN UTC
      StartHour = real(values(5)-timezone,kind=8) + &
                  real(values(6)/60.0,kind=8)
        ! find time in hours since BaseYear
      RunStartHour = HS_hours_since_baseyear(values(1),values(2),values(3),&
                                             StartHour,MR_BaseYear,MR_useLeap)
        ! get character string
      RunStartHour_ch = HS_yyyymmddhhmm_since(RunStartHour,MR_BaseYear,MR_useLeap)
      read(RunStartHour_ch,'(i4)') RunStartYear
      read(RunStartHour_ch,'(4x,i2)') RunStartMonth
      read(RunStartHour_ch,'(6x,i2)') RunStartDay
      read(RunStartHour_ch,'(8x,i2)') RunStartHr

      ! Find the start time given on the command line
      Probe_StartHour = HS_hours_since_baseyear(inyear,inmonth,inday,inhour,&
                                                MR_BaseYear,MR_useLeap)
      MR_Comp_StartHour     = Probe_StartHour
      MR_Comp_Time_in_hours = Simtime_in_hours

      write(*,*)"Checking"
      if(RunStartHour-Probe_StartHour.gt.24.0_8*GFS_Archive_Days)then
        ! NCEP case
        write(MR_global_info,*)"Requested start time is too old for GFS archive."
        write(MR_global_info,*)"Start time is older than the hardwired threshold of ",&
                  GFS_Archive_Days," days"
        write(MR_global_info,*)"Using NCEP 50-year Reanalysis"
        iw  = 5
        iwf = 25
        igrid   = 0
        iwfiles = 1

        call MR_Allocate_FullMetFileList(iw,iwf,igrid,2,iwfiles)

        do i=1,MR_iwindfiles
          write(MR_windfiles(i),*)trim(ADJUSTL(WINDROOT)), '/NCEP'
        enddo

      elseif(inyear.lt.1948)then
        ! Too old for NCEP runs, must use control file
        write(MR_global_info,*)"Requested start time is too old for NCEP Reanalysis."
        stop 1
      else
        ! GFS case
        write(MR_global_info,*)"Requested start time is within the GFS archive."
        write(MR_global_info,*)"Using archived global forecast data (GFS 0.5-degree)"
        if(RunStartHour-Probe_StartHour.lt.0.0_8)then
          ! GFS case for future run
          write(MR_global_info,*)"Requested start time is later than current time,"
          write(MR_global_info,*)"but it might fit in the current forecast package."
          if (Probe_StartHour-RunStartHour.ge.real(GFS_FC_TotHours,kind=8))then
            write(MR_global_info,*)" Run cannot complete with the current FC package."
            stop 1
          endif
        endif
        FC_intvl = 3.0_8

        ! calculate the number of forecast packages stored on system
        NumFCpackages = GFS_Archive_Days * (24/FC_freq)
        allocate(GFS_candidate(NumFCpackages))
        allocate(GFS_FC_step_avail(NumFCpackages))
        GFS_candidate = .true.
        write(MR_global_info,*)"Approximate Number of FC packages on system:",NumFCpackages
        write(MR_global_info,*)"Checking to see which packages might work for the"
        write(MR_global_info,*)"requested start-time"

        ! First, get the start hour of the forecast package immediately before
        ! the execution time
        FC_Package_year    = HS_YearOfEvent(RunStartHour,MR_BaseYear,MR_useLeap)
        FC_Package_month   = HS_MonthOfEvent(RunStartHour,MR_BaseYear,MR_useLeap)
        FC_Package_day     = HS_DayOfEvent(RunStartHour,MR_BaseYear,MR_useLeap)
        FC_Package_hour_dp = HS_HourOfDay(RunStartHour,MR_BaseYear,MR_useLeap)
        ! now round down the hour to the forecast package increment
        FC_Package_hour    = floor(FC_Package_hour_dp/FC_freq) * FC_freq
        ! and get the hours since base time for that FC package
        FC_Package_StartHour = HS_hours_since_baseyear(FC_Package_year,&
                                 FC_Package_month,FC_Package_day,real(FC_Package_hour,kind=8),&
                                 MR_BaseYear,MR_useLeap)
        ! estimate the start time of the oldest forecast package on the system
        FC_Archive_StartHour = FC_Package_StartHour - GFS_Archive_Days*24.0_8

        ! Loop through all the packages and check which ones might span the needed
        ! time range
        write(*,*)'Looping backward through packages'
        do i = NumFCpackages,1,-1
          write(*,*)"Package # ",i
          ! Get the start hour of this forecast package
          FCStartHour = FC_Archive_StartHour + real((i-1)*FC_freq,kind=8)
          if (FCStartHour.gt.Probe_StartHour)then
            ! This package starts after the requested time so dismiss it
            GFS_candidate(i)     = .false.
            GFS_FC_step_avail(i) = 0
            write(*,*)"   Package starts too late"
            cycle
          endif

          FCEndHour = FCStartHour + real(GFS_FC_TotHours,kind=8)
          if (FCEndHour.lt.Probe_StartHour)then
            ! This package ends before the needed time so dismiss it
            GFS_candidate(i)     = .false.
            GFS_FC_step_avail(i) = 0
            write(*,*)"   Package ends too early"
            cycle
          endif

          FC_year = HS_YearOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
          FC_mon  = HS_MonthOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
          FC_day  = HS_DayOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
          FC_hour = HS_HourOfDay(FCStartHour,MR_BaseYear,MR_useLeap)
          FC_Package_hour = floor(FC_hour/FC_freq) * FC_freq
          write(*,*)'   This one could work. Testing each file. ', &
                     FC_year,FC_mon,FC_day,real(FC_hour,kind=4)

          FC_hour_int = 0

          ! Now we loop through all the files in this package assuming 198 hours
          ! of forecast data and 3 hour steps (67 files overall)
            do ii = 1,67
              if(.not.GFS_candidate(i)) cycle  ! This will be true until a
                                               ! needed file is missing
              FC_hour_int = nint((ii-1)*FC_intvl)
              if (FCStartHour+FC_hour_int.lt.Probe_StartHour-FC_intvl) cycle ! to early
              if (FCStartHour+FC_hour_int.gt.Probe_StartHour+FC_intvl)   cycle ! to late

              ! if we are here, then we are inspecting a file that would be needed for
              ! the requested time span, starting with the most recent forecast package
              ! See if the file actually exists
              write(string1,'(a9,I4.4,I2.2,I2.2,I2.2,a1)')'/gfs/gfs.', &
                            FC_year,FC_mon,FC_day,FC_Package_hour,'/'
              write(string2,'(I4.4,I2.2,I2.2,I2.2,a2,I3.3,a3)')&
                            FC_year,FC_mon,FC_day,FC_Package_hour, &
                            '.f',FC_hour_int,'.nc'

              write(testfile,*)trim(ADJUSTL(WINDROOT)), &
                                   trim(ADJUSTL(string1)), &
                                   trim(ADJUSTL(string2))
              inquire( file=trim(adjustl(testfile)), exist=IsThere )
              if (IsThere)then
                GFS_FC_step_avail(i) = ii
                FCEndHour = FCStartHour + real(FC_hour_int,kind=8)
                write(MR_global_info,*)"     Found: ",trim(adjustl(testfile))
              else
                ! If a needed file is not available, mark the whole package as not a
                ! candidate
                GFS_candidate(i) = .false.
                cycle
              endif
            enddo

          ! If after all that testing, if we still have a valid candidate forecast package,
          ! set this package as the one we will use and exit the do loop
          if (GFS_candidate(i)) then
            OptimalPackageNum = i
            exit
          endif

        enddo ! 1,NumFCpackages

        if (OptimalPackageNum.eq.0)then
          write(MR_global_info,*)"No GFS package available on this system will span the"
          write(MR_global_info,*)"requested simulation time.  Exiting"
          write(MR_global_info,*)"  Probe start time     = ",Probe_StartHour
          write(MR_global_info,*)"  FC_Archive_StartHour =",FC_Archive_StartHour
          write(MR_global_info,*)"  FCStartHour          = ",FCStartHour
          write(MR_global_info,*)"  FCEndHour            = ",FCEndHour
          stop 1
        endif

        iw      = 4
        iwf     = 20
        igrid   = 0
        iwfiles = 34
        call MR_Allocate_FullMetFileList(iw,iwf,igrid,2,iwfiles)

          !  Now list the windfiles we will use and copy to MR_windfiles
          FCStartHour = FC_Archive_StartHour + real((OptimalPackageNum-1)*FC_freq,kind=8)
          FC_year = HS_YearOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
          FC_mon  = HS_MonthOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
          FC_day  = HS_DayOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
          FC_hour = HS_HourOfDay(FCStartHour,MR_BaseYear,MR_useLeap)
          FC_Package_hour = floor(FC_hour/FC_freq) * FC_freq

          do i=1,MR_iwindfiles

            FCStartHour = FC_Archive_StartHour + real((OptimalPackageNum-1)*FC_freq,kind=8)

            FC_year = HS_YearOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
            FC_mon  = HS_MonthOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
            FC_day  = HS_DayOfEvent(FCStartHour,MR_BaseYear,MR_useLeap)
            FC_hour = HS_HourOfDay(FCStartHour,MR_BaseYear,MR_useLeap)
            FC_Package_hour = floor(FC_hour/FC_freq) * FC_freq
            FC_hour_int = nint((i-1)*FC_intvl)
            write(string1,'(a9,I4.4,I2.2,I2.2,I2.2,a1)')'/gfs/gfs.', &
                          FC_year,FC_mon,FC_day,FC_Package_hour,'/'
            write(string2,'(I4.4,I2.2,I2.2,I2.2,a2,I3.3,a3)')&
                          FC_year,FC_mon,FC_day,FC_Package_hour, &
                          '.f',FC_hour_int,'.nc'

            write(MR_windfiles(i),*)trim(ADJUSTL(WINDROOT)), &
                                 trim(ADJUSTL(string1)), &
                                 trim(ADJUSTL(string2))
          enddo

      endif

        ! Check for existance and compatibility with simulation time requirements
      call MR_Read_Met_DimVars(FC_year)

      call MR_Set_Met_Times(Probe_StartHour, Simtime_in_hours)

      write(MR_global_production,*)"Sonde time: ",inyear,inmonth,inday,inhour
      write(MR_global_production,*)"Now       : ",RunStartYear,RunStartMonth,RunStartDay,RunStartHr
      write(MR_global_production,*)"FC   time : ",inyear,inmonth,inday,FC_Package_hour

      end subroutine GetWindFile

!##############################################################################
!##############################################################################
!
!  GetMetProfile
!
!  This subroutine interpolates the data onto the requested point/time and
!  writes out the profile to GFS_prof.dat.
!
!##############################################################################

      subroutine GetMetProfile(inlon,inlat,inyear,inmonth,inday,inhour)

      use MetReader

      implicit none

      real(kind=4), intent(in)  :: inlon
      real(kind=4), intent(in)  :: inlat
      integer     , intent(in)  :: inyear
      integer     , intent(in)  :: inmonth
      integer     , intent(in)  :: inday
      real(kind=8), intent(in)  :: inhour

      real(kind=8)       :: Probe_StartHour
      integer            :: ivar,i
      real(kind=4),dimension(:,:,:),allocatable :: AirTemp_meso_last_step_MetP_sp
      real(kind=4),dimension(:,:,:),allocatable :: AirTemp_meso_next_step_MetP_sp
      real(kind=4) :: tfrac,tc,xfrac,xc,yfrac,yc
      real(kind=4) :: a1,a2,a3,a4

      real(kind=4),dimension(:),allocatable :: GPHprof1,GPHprof2,GPHprof
      real(kind=4),dimension(:),allocatable :: tempprof1,tempprof2,tempprof
      real(kind=4) :: TropoH,TropoP,TropoT
      real(kind=4) :: lapse_1,lapse_2,lapse_3

      INTERFACE
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer     ,intent(in) :: iyear
          integer     ,intent(in) :: imonth
          integer     ,intent(in) :: iday
          real(kind=8),intent(in) :: hours
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_hours_since_baseyear
      END INTERFACE
      allocate(AirTemp_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(AirTemp_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(GPHprof(np_fullmet))
      allocate(GPHprof1(np_fullmet))
      allocate(GPHprof2(np_fullmet))
      allocate(tempprof(np_fullmet))
      allocate(tempprof1(np_fullmet))
      allocate(tempprof2(np_fullmet))

      ! First load the Met grids for Geopotential
      MR_iMetStep_Now = 1 ! This is initialized to 0
      call MR_Read_HGT_arrays(MR_iMetStep_Now)
      ivar = 5 ! Temperature
      call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now)
      AirTemp_meso_last_step_MetP_sp = MR_dum3d_MetP

      call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now+1)
      AirTemp_meso_next_step_MetP_sp = MR_dum3d_MetP

      
      ! Get the fractional time between forecast steps
      Probe_StartHour = HS_hours_since_baseyear(inyear,inmonth,inday,inhour,&
                                                MR_BaseYear,MR_useLeap)
      tfrac = real((Probe_StartHour-MR_MetStep_Hour_since_baseyear(1)) / &
               MR_MetStep_Interval(MR_iMetStep_Now),kind=4)
      tc    = 1.0_4-tfrac
      ! Get the fractional position in cell and corner weights
      xfrac=(inlon - x_submet_sp(1))/dx_met_const
      yfrac=(inlat - y_submet_sp(1))/dy_met_const
      xc = 1.0_4-xfrac
      yc = 1.0_4-yfrac
      a1=xc*yc
      a2=xfrac*yc
      a3=xfrac*yfrac
      a4=yfrac*xc
      GPHprof1  = a1*MR_geoH_metP_last(1,1,:) + &
                  a2*MR_geoH_metP_last(2,1,:) + &
                  a3*MR_geoH_metP_last(2,2,:) + &
                  a4*MR_geoH_metP_last(1,2,:)

      tempprof1 = a1*AirTemp_meso_last_step_MetP_sp(1,1,:) + &
                  a2*AirTemp_meso_last_step_MetP_sp(2,1,:) + &
                  a3*AirTemp_meso_last_step_MetP_sp(2,2,:) + &
                  a4*AirTemp_meso_last_step_MetP_sp(1,2,:)
      GPHprof2  = a1*MR_geoH_metP_next(1,1,:) + &
                  a2*MR_geoH_metP_next(2,1,:) + &
                  a3*MR_geoH_metP_next(2,2,:) + &
                  a4*MR_geoH_metP_next(1,2,:)

      tempprof2 = a1*AirTemp_meso_next_step_MetP_sp(1,1,:) + &
                  a2*AirTemp_meso_next_step_MetP_sp(2,1,:) + &
                  a3*AirTemp_meso_next_step_MetP_sp(2,2,:) + &
                  a4*AirTemp_meso_next_step_MetP_sp(1,2,:)

      open(unit=20,file='GFS_prof.dat')
      do i = 1,np_fullmet
        GPHprof(i) = tc*GPHprof1(i) + tfrac*GPHprof2(i)
        tempprof(i) = tc*tempprof1(i) + tfrac*tempprof2(i)
        write(20,*)GPHprof(i),p_fullmet_sp(i),tempprof(i)-273.0
      enddo
      close(20)

      ! Get Height of tropopause by calculating lapse rate
      do i = 2,np_fullmet-2
        lapse_1 = (tempprof(i-1)-tempprof(i  ))/(GPHprof(i  )-GPHprof(i-1))
        lapse_2 = (tempprof(i  )-tempprof(i+1))/(GPHprof(i+1)-GPHprof(i  ))
        lapse_3 = (tempprof(i+1)-tempprof(i+2))/(GPHprof(i+2)-GPHprof(i+1))
        if(lapse_1.gt.0.002.and.&
           lapse_2.lt.0.002.and.&
           lapse_3.lt.0.002)then
          TropoH = GPHprof(i)
          TropoT = tempprof(i)
          TropoP = p_fullmet_sp(i)
          exit
        endif
      enddo
      write(MR_global_production,*)"Tropopause Height, Temp, Pressure"
      write(MR_global_production,*)TropoH,TropoT,TropoP


      end subroutine GetMetProfile

!##############################################################################
!##############################################################################
