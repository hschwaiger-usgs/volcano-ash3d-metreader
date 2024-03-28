!##############################################################################
!##############################################################################
!
! MetCheck
!
! This stand-alone program tests windfiles for corrupt values.  The program
! loops through all points and checks the five variables (GPH, U, V, W, T)
! and verifies that there are no values outside of expected limits.  If any
! are outside the limits, the windfile is assumed to be corrupt and the program
! exits.  If all points are successfully read and 'valid', then a one-line
! file is written (MetCheck_log.txt) containing:
!   filename : step_in_file : year : days_in_year
! e.g.
! 2022092800.f099.nc  :            1  :         2022  :    275.125000
!
! There are 3 required and one optional command-line arguments:
!       iwf        : integer : product identifier
!       idf        : integer : data format (ASCII,netcdf,grib)
!       filename   : char    : name of file to test
!       [year]     : integer : year, needed to specify past NCEP data
!
!##############################################################################
!##############################################################################

      program MetCheck

      use MetReader,       only : &
         MR_nio,VB,verbosity_error,verbosity_info,verbosity_production,errlog,outlog,&
         MR_BaseYear,MR_useLeap,MR_nio,MR_Comp_StartHour,MR_Comp_Time_in_hours,&
         IsLatLon_MetGrid,Met_iprojflag,Met_k0,Met_lam0,Met_phi0,Met_phi1,Met_phi2,&
         Met_Re,MR_useCompGrid,MR_useCompTime,nt_fullmet,nx_fullmet,ny_fullmet,&
         MR_windfiles,MR_windfile_stephour,MR_dum3d_metP,nlevs_fullmet,&
         Met_dim_IsAvailable,Met_var_zdim_idx,MR_windfiles,&
           MR_Allocate_FullMetFileList,MR_windfile_starthour,&
           MR_Read_Met_DimVars,&
           MR_Set_Met_Times,&
           MR_Set_CompProjection,&
           MR_Initialize_Met_Grids,&
           MR_Read_HGT_arrays,&
           MR_Read_3d_MetP_Variable,&
           MR_Reset_Memory,&
           MR_FileIO_Error_Handler

      implicit none

      real(kind=4), parameter  :: H_MIN = -1000.0_4
      real(kind=4), parameter  :: H_MAX = 80000.0_4
      real(kind=4), parameter  :: U_MIN = -200.0_4
      real(kind=4), parameter  :: U_MAX =  200.0_4
      real(kind=4), parameter  :: V_MIN = -200.0_4
      real(kind=4), parameter  :: V_MAX =  200.0_4
      real(kind=4), parameter  :: W_MIN = -40.0_4
      real(kind=4), parameter  :: W_MAX =  40.0_4
      real(kind=4), parameter  :: T_MIN = 130.0_4
      real(kind=4), parameter  :: T_MAX = 350.0_4

      integer             :: nargs
      integer             :: iostatus
      integer             :: inlen
      character(len=120)  :: iomessage
      character(len=100)  :: arg

      real(kind=4)        :: inlon,inlat

      character(len=100)  :: infile1
      integer             :: nxmax,nymax,nzmax
      real(kind=4),dimension(:)    ,allocatable :: lon_grid
      real(kind=4),dimension(:)    ,allocatable :: lat_grid
      real(kind=4),dimension(:)    ,allocatable :: z_cc
      logical             :: IsPeriodic

      integer :: idf, igrid, iw, iwf, iwfiles
      integer :: ivar
      integer :: iy
      integer :: i,j,p,np,imetstep
      real(kind=4) :: v1,v2,tmp

      integer :: iprojflag
      real(kind=8) :: lambda0,phi0,phi1,phi2,k0,radius_earth
      logical :: IsLatLon

      integer :: BaseYear = 1900
      logical :: useLeap  = .true.
        ! min and max possible for this var (for error-checking)
      real(kind=4)    ,dimension(50,2) :: Met_var_MinMax
      integer,dimension(8) :: values
      integer :: Current_Year
      real(kind=8) :: hsince
      integer      :: idx

      integer :: io                           ! Index for output streams
      character(len=80)  :: linebuffer080

      INTERFACE
        subroutine Print_Usage
        end subroutine Print_Usage
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer     ,intent(in) :: iyear
          integer     ,intent(in) :: imonth
          integer     ,intent(in) :: iday
          real(kind=8),intent(in) :: hours
          integer     ,intent(in) :: byear
          logical     ,intent(in) :: useLeaps
        end function HS_hours_since_baseyear
        integer function HS_YearOfEvent(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_YearOfEvent
        integer function HS_DayOfEvent(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_DayOfEvent
        integer function HS_DayOfYear(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_DayOfYear
        real(kind=8) function HS_HourOfDay(HoursSince,byear,useLeaps)
          real(kind=8)   ,intent(in) ::  HoursSince
          integer        ,intent(in) ::  byear
          logical        ,intent(in) ::  useLeaps
        end function HS_HourOfDay
      END INTERFACE

      Met_var_MinMax(1,1:2) = (/ H_MIN , H_MAX  /)  ! GPH
      Met_var_MinMax(2,1:2) = (/ U_MIN , U_MAX  /)  ! U
      Met_var_MinMax(3,1:2) = (/ V_MIN , V_MAX  /)  ! V
      Met_var_MinMax(4,1:2) = (/ W_MIN , W_MAX  /)  ! W
      Met_var_MinMax(5,1:2) = (/ T_MIN , T_MAX  /)  ! T

      ! Make user MetReader is using the same calendar
      MR_BaseYear = BaseYear
      MR_useLeap  = useLeap

!     TEST READ COMMAND LINE ARGUMENTS
      nargs = command_argument_count()
      if (nargs.lt.3) then
        ! Not enough arguments; print usage and exit
        call Print_Usage
      else
        ! Get wind product
        call get_command_argument(1, arg, length=inlen, status=iostatus)
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read first command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)iwf
        linebuffer080(1:inlen) = arg
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check iwf
        if(iwf.ne.3.and.&
           iwf.ne.4.and.&
           iwf.ne.5.and.&
           iwf.ne.6.and.&
           iwf.ne.7.and.&
           iwf.ne.8.and.&
           iwf.ne.10.and.&
           iwf.ne.11.and.&
           iwf.ne.12.and.&
           iwf.ne.13.and.&
           iwf.ne.20.and.&
           iwf.ne.21.and.&
           iwf.ne.22.and.&
           iwf.ne.23.and.&
           iwf.ne.24.and.&
           iwf.ne.25.and.&
           iwf.ne.26.and.&
           iwf.ne.27.and.&
           iwf.ne.28.and.&
           iwf.ne.29.and.&
           iwf.ne.30.and.&
           iwf.ne.32.and.&
           iwf.ne.33.and.&
           iwf.ne.41.and.&
           iwf.ne.42)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: windformat not recognized"
            write(errlog(io),*)" iwf = ",iwf
          endif;enddo
          stop 1
        endif

        ! Get data format ID
        call get_command_argument(2, arg, length=inlen, status=iostatus)
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read second command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)idf
        linebuffer080(1:inlen) = arg
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check idf
        if(idf.ne.2.and.&
           idf.ne.3)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: Only netcdf (2) and grib (3) supported"
            write(errlog(io),*)" idf = ",idf
          endif;enddo
          stop 1
        endif

        call get_command_argument(3, arg, length=inlen, status=iostatus)
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read third command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif

        infile1 = trim(arg)

        ! Now check for optional  argument for year (use for checking NCEP files
        call date_and_time(VALUES=values)
        Current_Year = values(1)
        if(nargs.gt.3)then
          call get_command_argument(4, arg, length=inlen, status=iostatus)
          if(iostatus.ne.0)then
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)"MR ERROR : Could not read fourth command-line argument"
              write(errlog(io),*)" arg = ",arg
            endif;enddo
            call Print_Usage
          endif
          read(arg,*,iostat=iostatus,iomsg=iomessage)iy
          linebuffer080(1:inlen) = arg
          if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        else
          iy = Current_Year
        endif
        !MR_Comp_StartHour     = HS_hours_since_baseyear(iy,1,1,0.0_8,1900,.True.)
        !MR_Comp_Time_in_hours = 1.0
      endif
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Set up windfile data structure"
      endif;enddo
      if(iwf.eq.25)then
        iw      = 5
        MR_Comp_StartHour     = HS_hours_since_baseyear(iy,1,1,0.0_8,1900,.True.)
        MR_Comp_Time_in_hours = 1.0
      else
        iw      = 4
      endif
      igrid   = 0  ! this will get reset in MR_Allocate_FullMetFileList
      iwfiles = 1
      call MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)
      write(MR_windfiles(1),*)trim(adjustl(infile1))

        ! Check for windfile existance and read their sizes
      call MR_Read_Met_DimVars(iy)

      ! This allows us to pass zero's for starttime and duration and get
      ! the whole metstep list
      MR_useCompTime = .false.
      call MR_Set_Met_Times(0.0_8,0.0_8)

      ! These are dummy values.
      nxmax = 2 !
      nymax = 2 !

      nzmax = 2 ! This is not really used in this utility
      inlon = 0.0_4
      inlat = 0.0_4
      allocate(lon_grid(nxmax)); lon_grid(1:nxmax) = (/inlon-0.5_4,inlon+0.5_4/)
      allocate(lat_grid(nymax)); lat_grid(1:nymax) = (/inlat-0.5_4,inlat+0.5_4/)

      allocate(z_cc(nzmax))    ; z_cc(1:2) = (/0.0_4, 10.0_4/)
      IsPeriodic = .false.
      IsPeriodic = .true.

      ! We just want to access the met grid, so set our 'comp' grid to the
      ! same projection
      IsLatLon     = IsLatLon_MetGrid
      iprojflag    = Met_iprojflag
      lambda0      = Met_lam0
      phi0         = Met_phi0
      phi1         = Met_phi1
      phi2         = Met_phi2
      k0           = Met_k0
      radius_earth = Met_Re
      call MR_Set_CompProjection(IsLatLon,iprojflag,lambda0,phi0,phi1,phi2,&
                                 k0,radius_earth)

      ! This program needs no interpolation so we need to explicitly
      ! declare that we do not need a computational grid defined
      MR_useCompGrid = .false.

      call MR_Initialize_Met_Grids(nxmax,nymax,nzmax,&
                              lon_grid(1:nxmax), &
                              lat_grid(1:nymax), &
                              z_cc(1:nzmax)    , &
                              IsPeriodic)
      ! We need to populate the HGT arrays or the vertical velocity read will fail
      ! since the conversion relies on dp/dz
      call MR_Read_HGT_arrays(1)  ! Just fill with the geopotential height at step 1

      open(unit=19,iostat=iostatus,file='MetCheck_log.txt',status='old')
      if (iostatus == 0) close(19, status='delete')

      open(unit=19,file='MetCheck_log.txt',status='replace')

      do imetstep = 1,nt_fullmet
        do ivar=1,5
          idx = Met_var_zdim_idx(ivar)
          if(Met_dim_IsAvailable(ivar).eqv..true.)then
            np = nlevs_fullmet(idx)
            v1 = Met_var_MinMax(ivar,1)
            v2 = Met_var_MinMax(ivar,2)
            call MR_Read_3d_MetP_Variable(ivar,imetstep)
            do i=1,nx_fullmet
              do j=1,ny_fullmet
                do p=1,np
                  tmp=MR_dum3d_metP(i,j,p)
                  if(tmp.lt.v1.or.tmp.gt.v2)then
                    do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                      write(errlog(io),*)"MR ERROR reading value for ivar=",ivar
                      write(errlog(io),*)"         at i,j,p = ",i,j,p
                      write(errlog(io),*)"         Value read = ",tmp
                      write(errlog(io),*)"         Min / Max = ",v1,v2
                    endif;enddo
                    close(19)
                    stop 1
                  endif
                enddo
              enddo
            enddo
          endif
        enddo

        if(iwf.eq.25)then
          write(19,*)adjustl(trim(MR_windfiles(1)))," : ", &
                     imetstep                      ," : ", &
                     iy                            ," : ", &
                     real((imetstep-1),kind=4)*6.0_4/24_4
        else
          hsince = MR_windfile_starthour(1)+MR_windfile_stephour(1,imetstep)
          write(19,*)adjustl(trim(MR_windfiles(1)))," : ", &
                     imetstep                      ," : ", &
                     HS_YearOfEvent(hsince,MR_BaseYear,MR_useLeap)," : ",&
                     real(HS_DayOfYear(hsince,MR_BaseYear,MR_useLeap)+ &
                       HS_HourOfDay(hsince,MR_BaseYear,MR_useLeap)/24.0_8,kind=4)
        endif

      enddo

      close(19)

      call MR_Reset_Memory

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"Program ended normally."
      endif;enddo

      end program MetCheck

!##############################################################################
!
!  Print_Usage
!
!  This subroutine is called if there is an error reading the command-line.
!  Expected usage is written to stdout and the program exits.
!
!##############################################################################

      subroutine Print_Usage

      use MetReader,       only : &
         MR_nio,VB,verbosity_error,errlog,VB

      implicit none

      integer :: io

        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: not enough command-line arguments."
          write(errlog(io),*)"  Usage: MetCheck iwf idf filename [year]"
          write(errlog(io),*)"   where "
          write(errlog(io),*)"     iwf =  3 NARR3D NAM221 32 km North America files"
          write(errlog(io),*)"     iwf =  4 NARR3D NAM221 32 km North America files"
          write(errlog(io),*)"     iwf =  5 NAM216 AK 45km"
          write(errlog(io),*)"     iwf =  6 NAM Regional 90 km grid 104"
          write(errlog(io),*)"     iwf =  7 CONUS 212 40km"
          write(errlog(io),*)"     iwf =  8 CONUS 218 (12km)"
          write(errlog(io),*)"     iwf = 10 NAM 242 11.25 km AK"
          write(errlog(io),*)"     iwf = 11 NAM 196 2.5 km HI"
          write(errlog(io),*)"     iwf = 12 NAM 198 5.953 km AK"
          write(errlog(io),*)"     iwf = 13 NAM 91 2.976 km AK"
          write(errlog(io),*)"     iwf = 20 GFS 0.5"
          write(errlog(io),*)"     iwf = 21 GFS 1.0"
          write(errlog(io),*)"     iwf = 22 GFS 0.25"
          write(errlog(io),*)"     iwf = 23 NCEP / DOE reanalysis 2.5 degree files"
          write(errlog(io),*)"     iwf = 24 NASA-MERRA-2 reanalysis 0.625x0.5 degree files"
          write(errlog(io),*)"     iwf = 25 NCEP/NCAR reanalysis 2.5 degree files"
          write(errlog(io),*)"     iwf = 27 NOAA-CIRES reanalysis 2.5 degree files"
          write(errlog(io),*)"     iwf = 28 ECMWF Interim Reanalysis (ERA-Interim)"
          write(errlog(io),*)"     iwf = 29 ECMWF ERA5"
          write(errlog(io),*)"     iwf = 30 ECMWF 20-Century (ERA-20C)"
          write(errlog(io),*)"     iwf = 32 Air Force Weather Agency subcenter = 0"
          write(errlog(io),*)"     iwf = 33 CCSM3.0 Community Atmosphere Model (CAM)"
          write(errlog(io),*)"     iwf = 40 NASA-GEOS Cp"
          write(errlog(io),*)"     iwf = 41 NASA-GEOS Np"
          write(errlog(io),*)"     iwf = 50 WRF"
          write(errlog(io),*)"     "
          write(errlog(io),*)"     idf = 2 NetCDF"
          write(errlog(io),*)"     idf = 3 grib"
          write(errlog(io),*)"     "
          write(errlog(io),*)"     filename = name of file of root directory for NCEP"
          write(errlog(io),*)"     "
          write(errlog(io),*)"     [year] = year for NCEP tests"
          write(errlog(io),*)"               This is optional and defaults to current year"
          write(errlog(io),*)"               if not provided."
        endif;enddo
        stop 1

      end subroutine Print_Usage

!##############################################################################
