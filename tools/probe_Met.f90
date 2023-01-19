!##############################################################################
!##############################################################################
!
!  probe_Met
!
!    This is a stand-alone program that exports a vertical profile of 
!    pressure variables from windfiles
!    This program takes 15+ command-line arguments:
!       file         : string  : name of input file
!       timestep     : integer : 
!       llflag       : integer : 0 for using windfile grid; 1 for forcing Lat/Lon
!       lon/x        : real    : longitude (or x) of sonde point
!       lat/y        : real    : latitude (or y) of sonde point
!       trunc flag   : char    : truncation flag (T or F)
!       nvars        : integer : number of pressure variables to export
!       varID(nvars) : integers: variable ID's to read and export
!       iw           : integer : windfile format code (3,4,5)
!       iwf          : integer : windfile product code
!       idf          : integer : igrid (2 for nc, 3 for grib)
!       year         : integer : only needed for iw = 5 files
!       month        : integer : 
!       day          : integer : 
!       hour         : real    : 

!               if true, then truncate x and y to gridpoint and return
!               the native values
!
!    To probe a 0.5-degree GFS file in nc format for GPH only, use:
!      probe_Met 2020061000_5.f006.nc 1 0 190.055 52.8222 T 1 1 4 20 2
!     same, but adding u,v,t
!      probe_Met 2020061000_5.f006.nc 1 0 190.055 52.8222 T 4 1 2 3 5 4 20 2
!    To probe a NAM 91 grid over AK in nc format with a LL coordinate, use:
!      probe_Met nam.tm06.grib2 1 1 190.055 52.8222 F 4 1 2 3 5 4 13 3
!     with a projected coordinate:
!      probe_Met nam.tm06.grib2 1 0 -1363.94 -3758.61 F 4 1 2 3 5 4 13 3
!    To  probe the NCEP 2.5-degree data (or other iw=5), we need the
!    full date
!      probe_Met 2020061000_5.f006.nc 1 1 190.055 52.8222 T 1 1 4 20 2 2018 1 1 0.0
!    Output is written to the file sonde_ztuvrp.dat
!
! Note: to test against the 0.5 GFS data using ncks, use the following commands
!         52.5N is index 76; 190.0W is index 381
! ncks -C -v isobaric3 -H -Q -s '%f ' 2020061000_5.f006.nc > GFS50_P.dat
! ncks -C -d lon,190.0 -d lat,52.5 -H -Q -s '%f ' -v Geopotential_height_isobaric 2020061000_5.f006.nc  > GFS50_H.dat
! ncks -C -d lon,190.0 -d lat,52.5 -H -Q -s '%f ' -v Temperature_isobaric 2020061000_5.f006.nc > GFS50_T.dat
! ncks -C -d lon,190.0 -d lat,52.5 -H -Q -s '%f ' -v u-component_of_wind_isobaric 2020061000_5.f006.nc > GFS50_U.dat
! ncks -C -d lon,190.0 -d lat,52.5 -H -Q -s '%f ' -v v-component_of_wind_isobaric 2020061000_5.f006.nc > GFS50_V.dat
!
!         -1363.94 at index 423; -3758.61 at index 354
! ncks -C -v isobaric -H -Q -s '%f ' nam.tm06.grib2.nc > NAM_P.dat
! ncks -C -d x,-1363.525269 -d y,-3759.574707 -H -Q -s '%f ' -v Geopotential_height_isobaric nam.tm06.grib2.nc  > NAM_H.dat
! ncks -C -d x,-1363.525269 -d y,-3759.574707 -H -Q -s '%f ' -v Temperature_isobaric nam.tm06.grib2.nc > NAM_T.dat
! ncks -C -d x,-1363.525269 -d y,-3759.574707 -H -Q -s '%f ' -v u-component_of_wind_isobaric nam.tm06.grib2.nc > NAM_U.dat
! ncks -C -d x,-1363.525269 -d y,-3759.574707 -H -Q -s '%f ' -v v-component_of_wind_isobaric nam.tm06.grib2.nc > NAM_V.dat

!##############################################################################


      program probe_Met

      use MetReader

      implicit none

      integer             :: nargs
      integer             :: status
      character (len=130) :: arg

      character(len=100)  :: infile
      integer             :: intstep
      integer             :: inLLflag
      real(kind=4)        :: inlon,inlat
      character           :: intrunc
      integer             :: invars
      integer,dimension(:)    ,allocatable :: invarlist
      integer             :: iw,iwf,idf,igrid,iwfiles
      integer             :: inyear,inmonth,inday
      real(kind=8)        :: inhour

      integer             :: nxmax,nymax,nzmax
      real(kind=4),dimension(:)    ,allocatable :: lon_grid
      real(kind=4),dimension(:)    ,allocatable :: lat_grid
      real(kind=4),dimension(:)    ,allocatable :: z_cc
      logical             :: IsPeriodic

      integer :: BaseYear = 1900
      logical :: useLeap  = .true.
      logical :: IsThere
      logical :: Truncate = .false.
      integer :: i
      real(kind=8) :: steptime

      INTERFACE
        real(kind=8) function HS_HourOfDay(HoursSince,byear,useLeaps)
          real(kind=8)          :: HoursSince
          integer               :: byear
          logical               :: useLeaps
        end function HS_HourOfDay
        integer function HS_YearOfEvent(HoursSince,byear,useLeaps)
          real(kind=8)          :: HoursSince
          integer               :: byear
          logical               :: useLeaps
        end function HS_YearOfEvent
        integer function HS_MonthOfEvent(HoursSince,byear,useLeaps)
          real(kind=8)          :: HoursSince
          integer               :: byear
          logical               :: useLeaps
        end function HS_MonthOfEvent

        integer function HS_DayOfEvent(HoursSince,byear,useLeaps)
          real(kind=8)          :: HoursSince
          integer               :: byear
          logical               :: useLeaps
        end function HS_DayOfEvent
        integer function HS_DayOfYear(HoursSince,byear,useLeaps)
          real(kind=8)          :: HoursSince
          integer               :: byear
          logical               :: useLeaps
        end function HS_DayOfYear

        subroutine GetMetProfile(invars,invarlist)
          integer,parameter   :: sp        = 4 ! single precision
          integer,parameter   :: dp        = 8 ! double precision
          integer             :: invars
          integer,dimension(invars) :: invarlist
        end subroutine
      END INTERFACE

      ! Make sure user MetReader is using the same calendar
      MR_BaseYear = BaseYear
      MR_useLeap  = useLeap
      MR_useCompH = .false.

!     TEST READ COMMAND LINE ARGUMENTS
      nargs = command_argument_count()
      if (nargs.lt.9) then
        write(MR_global_error,*)"ERROR: insufficient command-line arguments"
        write(MR_global_error,*)"Usage: probe_Met filename tstep llflag lon lat trunc nvars ",&
                                "vars(nvars) iw iwf (inyear inmonth inday inhour)"
        write(MR_global_error,*)"       file         : string  : name of input file"
        write(MR_global_error,*)"       timestep     : integer : step in file"
        write(MR_global_error,*)"       llflag       : integer : 0 for using windfile grid; 1 for"
        write(MR_global_error,*)"                                  forcing Lat/Lon"
        write(MR_global_error,*)"       lon/x        : real    : longitude (or x) of sonde point"
        write(MR_global_error,*)"       lat/y        : real    : latitude (or y) of sonde point"
        write(MR_global_error,*)"       trunc flag   : char    : truncation flag (T or F)"
        write(MR_global_error,*)"                                  T = coordinates truncated to nearest met node"
        write(MR_global_error,*)"                                  F = values interpolated onto given coordinate"
        write(MR_global_error,*)"       nvars        : integer : number of pressure variables to export"
        write(MR_global_error,*)"       varID(nvars) : integers: variable ID's to read and export"
        write(MR_global_error,*)"                        1 = GPH (km)"
        write(MR_global_error,*)"                        2 = U (m/s)"
        write(MR_global_error,*)"                        3 = V (m/s)"
        write(MR_global_error,*)"                        4 = W (m/s)"
        write(MR_global_error,*)"                        5 = T (K)"
        write(MR_global_error,*)"       iw           : integer : windfile format code (3,4,5)"
        write(MR_global_error,*)"       iwf          : integer : windfile product code"
        write(MR_global_error,*)"       idf          : integer : igrid (2 for nc, 3 for grib)"
        write(MR_global_error,*)"       year         : integer : only needed for iw = 5 files"
        write(MR_global_error,*)"       month        : integer : "
        write(MR_global_error,*)"       day          : integer : "
        write(MR_global_error,*)"       hour         : real    : "
        write(MR_global_error,*)"  "
        write(MR_global_error,*)"     To probe a 0.5-degree GFS file in nc format for GPH only, use:"
        write(MR_global_error,*)"       probe_Met 2020061000_5.f006.nc 1 0 190.055 52.8222 T 1 1 4 20 2"
        write(MR_global_error,*)"      same, but adding u,v,t"
        write(MR_global_error,*)"       probe_Met 2020061000_5.f006.nc 1 0 190.055 52.8222 T 4 1 2 3 5 4 20 2"
        write(MR_global_error,*)"     To probe a NAM 91 grid over AK in nc format with a LL coordinate, use:"
        write(MR_global_error,*)"       probe_Met nam.tm06.grib2 1 1 190.055 52.8222 F 4 1 2 3 5 4 13 3"
        write(MR_global_error,*)"      with a projected coordinate:"
        write(MR_global_error,*)"       probe_Met nam.tm06.grib2 1 0 -1363.94 -3758.61 F 4 1 2 3 5 4 13 3"
        write(MR_global_error,*)"     To  probe the NCEP 2.5-degree data (or other iw=5), we need the"
        write(MR_global_error,*)"     full date"
        write(MR_global_error,*)"       probe_Met 2020061000_5.f006.nc 1 1 190.055 52.8222 T 1 1 4 20 2 2018 1 1 0.0"
        write(MR_global_error,*)"     Output is written to the file NWP_prof.dat"

        stop 1
      else
        ! Get file name or windroot for iw=5
        call get_command_argument(1, arg, status)
        read(arg,*)infile
        inquire( file=infile, exist=IsThere )
        if(.not.IsThere)then
          write(MR_global_error,*)"ERROR: Cannot find input file"
          stop 1
        endif
        write(MR_global_info,*)"infile = ",infile
        ! Get time step to use
        call get_command_argument(2, arg, status)
        read(arg,*)intstep
        write(MR_global_info,*)"intstep = ",intstep

        ! Get Lat/Lon flag
        call get_command_argument(3, arg, status)
        read(arg,*)inLLflag
        if(inLLflag.eq.0)then
          write(MR_global_info,*)"Lat/Lon flag not set."
          write(MR_global_info,*)"Using native grid of windfile (which may be LL)."
        elseif(inLLflag.eq.1)then
          write(MR_global_info,*)"Lat/Lon flag is set."
          write(MR_global_info,*)"Using a Lat/Lon grid regardless of the windfile."
        else
          write(MR_global_error,*)"ERROR: inLLflag must be either 0 or 1"
          stop 1
        endif
        ! Get lon and lat
        call get_command_argument(4, arg, status)
        read(arg,*)inlon
        ! round to the nearest third decimel
        inlon = nint(inlon * 1000.0_4) *0.001_4
        write(MR_global_info,*)"inlon = ",inlon
        call get_command_argument(5, arg, status)
        read(arg,*)inlat
        inlat = nint(inlat * 1000.0_4) *0.001_4
        write(MR_global_info,*)"inlat = ",inlat

        ! Get truncation flag
        call get_command_argument(6, arg, status)
        read(arg,*)intrunc
        if(intrunc.eq.'t'.or.intrunc.eq.'T')then
          Truncate = .true.
        else
          Truncate = .false.
        endif
        if(Truncate)then
          write(MR_global_info,*)"Truncate = ",Truncate, "output will be on met node"
        else
          write(MR_global_info,*)"Truncate = ",Truncate, &
                    "output will be on interpolated point"
        endif

        ! Get number of variables to read/export
        call get_command_argument(7, arg, status)
        read(arg,*)invars
        write(MR_global_info,*)"invars = ",invars
        write(MR_global_info,*)"Allocating var list of length ",invars
        allocate(invarlist(invars))

        ! Get variable IDs
        do i=1,invars
          call get_command_argument(7+i, arg, status)
          read(arg,*)invarlist(i)
          write(MR_global_info,*)" Requested variable ID ",invarlist(i)
        enddo

        ! Get wind format, wind product, and data format ID's
        call get_command_argument(8+invars, arg, status)
        read(arg,*)iw
        call get_command_argument(9+invars, arg, status)
        read(arg,*)iwf
        call get_command_argument(10+invars, arg, status)
        read(arg,*)idf
        write(MR_global_info,*)"Expecting wind format, NWP product ID and data type:",iw,iwf,idf

        ! Get year, month, day, hour
        if(nargs.ge.11+invars)then
          call get_command_argument(11+invars, arg, status)
          read(arg,*)inyear
        else
          inyear=0
        endif
        if(nargs.ge.12+invars)then
          call get_command_argument(12+invars, arg, status)
          read(arg,*)inmonth
        else
          inmonth=1
        endif
        if(nargs.ge.13+invars)then
          call get_command_argument(13+invars, arg, status)
          read(arg,*)inday
        else
          inday=1
        endif
        if(nargs.ge.14+invars)then
          call get_command_argument(14+invars, arg, status)
          read(arg,*)inhour
        else
          inhour=0.0
        endif
        if(inyear.gt.0)then
          write(MR_global_info,*)"Read date/hour as ",inyear,inmonth,inday,inhour
          write(MR_global_info,*)"This will overwrite the requested timestep above"
          write(MR_global_info,*)"only when using iw=5"
        else
          write(MR_global_info,*)"Optional arguments for date/hour not provided"
        endif

      endif

      if(Truncate)then
        write(MR_global_info,*)"Truncating profile onto ",inlon,inlat
      else
        write(MR_global_info,*)"Interpolating profile onto ",inlon,inlat
      endif

      igrid   = 0 ! 
      iwfiles = 1 ! single command-line argument or folder path for iw=5
      call MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)
      MR_windfiles(1) = trim(adjustl(infile))
      if(inyear.gt.0)then
        call MR_Read_Met_DimVars(inyear)
      else
        call MR_Read_Met_DimVars()
      endif
      ! The x and y Met grids are now read and loaded.

      ! Determine where inlon and inlat fall on this grid
      if(inLLflag.eq.0)then
        ! Do not force a Lon/Lat grid, just use the wind grid
        ! Map case is either
        !(1) Both Comp Grid and Met grids are Lat/Lon
        !(2) Both Comp Grid and Met grids are the same projection
        IsLatLon_CompGrid = IsLatLon_MetGrid
        Comp_iprojflag    = Met_iprojflag
      else
        ! Using a Lon/Lat comp grid
        ! Map case is either
        !(1) Both Comp Grid and Met grids are Lat/Lon
        !(4) Met Grid is projected and Comp grid is Lat/Lon
        IsLatLon_CompGrid = .true.
        Comp_iprojflag    = 0
      endif
      if(IsLatLon_CompGrid)then
        if(inlon.lt.-360.0)then
          write(MR_global_error,*)"ERROR: Longitude must be gt -360"
          stop 1
        endif
        if(inlon.lt.0.0_4.or.inlon.gt.360.0_4)inlon=mod(inlon+360.0_4,360.0_4)
      endif

      if(Truncate)then
        write(MR_global_info,*)&
          "Truncating input coordinate to grid point on Met grid."
        do i=1,nx_fullmet-1
          if(x_fullmet_sp(i).le.inlon.and.x_fullmet_sp(i+1).gt.inlon)then
            inlon = x_fullmet_sp(i)
          endif
        enddo
        if(y_inverted)then
          do i=1,ny_fullmet-1
            if(y_fullmet_sp(i).gt.inlat.and.y_fullmet_sp(i+1).le.inlat)then
              inlat = y_fullmet_sp(i+1)
            endif
          enddo
        else
          do i=1,ny_fullmet-1
            if(y_fullmet_sp(i).le.inlat.and.y_fullmet_sp(i+1).gt.inlat)then
              inlat = y_fullmet_sp(i)
            endif
          enddo
        endif
        write(MR_global_info,*)"Coordinate reset to ",inlon,inlat
      endif

      nxmax = 3 ! 
      nymax = 3 ! 
      nzmax = 2 ! This is not really used in this utility
      allocate(lon_grid(nxmax))
      allocate(lat_grid(nymax))

      if(IsLatLon_CompGrid)then
        lon_grid(1:3) = (/inlon-0.5_4,inlon,inlon+0.5_4/)
        lat_grid(1:3) = (/inlat-0.5_4,inlat,inlat+0.5_4/)
      else
        lon_grid(1:3) = (/inlon-1.0_4*abs(dx_met_const),inlon,inlon+1.0_4*abs(dx_met_const)/)
        lat_grid(1:3) = (/inlat-1.0_4*abs(dy_met_const),inlat,inlat+1.0_4*abs(dy_met_const)/)
      endif
      write(MR_global_info,*)"lon grid = ",lon_grid
      write(MR_global_info,*)"lat_grid = ",lat_grid

      allocate(z_cc(nzmax))    ; z_cc(1:2) = (/0.0_4, 10.0_4/)

      IsPeriodic = .false.
      call MR_Set_CompProjection(isLatLon_CompGrid,Comp_iprojflag, & ! if LL, only these matter
                                 Met_lam0,&
                                 Met_phi0,Met_phi1,Met_phi2,&
                                 Met_k0,Met_Re)

      write(MR_global_info,*)"Setting up wind grids"
      call MR_Initialize_Met_Grids(nxmax,nymax,nzmax,&
                              lon_grid(1:nxmax), &
                              lat_grid(1:nymax), &
                              z_cc(1:nzmax)    , &
                              IsPeriodic)

      ! Sort out what to do about the time
      steptime = MR_windfile_starthour(1) + MR_windfile_stephour(1,intstep)

      call MR_Set_Met_Times(steptime,0.0_8)

      do i=1,50
        write(MR_global_info,*)i,Met_var_NC_names(i),Met_var_IsAvailable(i)
      enddo
      
      call GetMetProfile(invars,invarlist)

      write(MR_global_info,*)"Program ended normally."

      end program probe_Met

!##############################################################################
!
!     Subroutines
!
!##############################################################################

!##############################################################################
!##############################################################################
!
!  GetMetProfile
!
!  This subroutine interpolates the data onto the requested point/time and
!  writes out the profile to NWP_prof.dat.  This interpolation is accomplished
!  by using a 3x3 computational grid centered on the point in question, then
!  writing out the middle point values.
!
!##############################################################################

      subroutine GetMetProfile(invars,invarlist)

      use MetReader

      implicit none

      integer                   :: invars
      integer,dimension(invars) :: invarlist

      integer                                 :: ivar,i,iv

      real(kind=4),dimension(:,:),allocatable :: outvars
      real(kind=4),dimension(:)  ,allocatable :: u,v

      allocate(outvars(invars,np_fullmet))
      ! We only need u and v for Map_Case = 5, but allocate anyway since
      ! these are small
      allocate(u(np_fullmet))
      allocate(v(np_fullmet))
      write(MR_global_info,*)" Inside GetMetProfile"

      ! First load the Met grids for Geopotential
      MR_iMetStep_Now = 1 ! This is initialized to 0
      call MR_Read_HGT_arrays(MR_iMetStep_Now)

      do iv = 1,invars
        ivar = invarlist(iv)
        write(MR_global_info,*)"Map_Case  and varID = ",Map_Case,iv
        if (Map_Case.eq.4.and. & ! Only need to do this if Met=proj and
                                 ! comp is LL
            (iv.eq.2.or.iv.eq.3))then
          write(MR_global_info,*)"MR Calling MR_Rotate_UV_GR2ER_Met to rotate vec"
          ! if a vector call, store u and v to rotate
          call MR_Rotate_UV_GR2ER_Met(MR_iMetStep_Now,.true.)
          u(1:np_fullmet) = MR_dum3d_compP(2,2,1:np_fullmet)
          v(1:np_fullmet) = MR_dum3d_compP_2(2,2,1:np_fullmet)
          if(iv.eq.2)outvars(iv,1:np_fullmet) =u(1:np_fullmet)
          if(iv.eq.3)outvars(iv,1:np_fullmet) =v(1:np_fullmet)
        else
          write(MR_global_info,*)"Calling MR_Read_3d_Met_Variable_to_CompP"
          call MR_Read_3d_Met_Variable_to_CompP(ivar,MR_iMetStep_Now)
          outvars(iv,1:np_fullmet) = MR_dum3d_compP(2,2,1:np_fullmet)
        endif

      enddo

      open(unit=20,file='NWP_prof.dat')
      do i = 1,np_fullmet
        write(20,*)p_fullmet_sp(i),outvars(:,i)
      enddo
      close(20)

      end subroutine GetMetProfile

!##############################################################################
!##############################################################################

