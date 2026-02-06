!##############################################################################
!##############################################################################
!
!  MetProbe
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
!      MetProbe 2020061000_5.f006.nc 1 0 190.055 52.8222 T 1 1 4 20 2
!     same, but adding u,v,t
!      MetProbe 2020061000_5.f006.nc 1 0 190.055 52.8222 T 4 1 2 3 5 4 20 2
!    To probe a NAM 91 grid over AK in nc format with a LL coordinate, use:
!      MetProbe nam.tm06.grib2 1 1 190.055 52.8222 F 4 1 2 3 5 4 13 3
!     with a projected coordinate:
!      MetProbe nam.tm06.grib2 1 0 -1363.94 -3758.61 F 4 1 2 3 5 4 13 3
!    To  probe the NCEP 2.5-degree data (or other iw=5), we need the
!    full date
!      MetProbe 2020061000_5.f006.nc 1 1 190.055 52.8222 T 1 1 4 20 2 2018 1 1 0.0
!    Output is written to the file NWP_prof.dat
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

      program MetProbe

      use MetReader,       only : &
         MR_nio,MR_VB,outlog,errlog,verbosity_error,verbosity_info,&
         Comp_iprojflag,MR_MAXVARS,&
         Met_iprojflag,Met_k0,Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_Re,&
         dx_met_const,dy_met_const,IsLatLon_CompGrid,MR_Comp_StartHour,&
         IsLatLon_MetGrid,MR_BaseYear,MR_useLeap,MR_useCompH,nx_fullmet,ny_fullmet,&
         y_inverted,MR_windfiles,x_fullmet_sp,y_fullmet_sp,&
         MR_windfile_starthour,MR_windfile_stephour,Met_var_NC_names,Met_var_IsAvailable,&
           MR_Read_Met_DimVars,&
           MR_Allocate_FullMetFileList,&
           MR_Set_CompProjection,&
           MR_Initialize_Met_Grids,&
           MR_Set_Met_Times,&
           MR_Reset_Memory,&
           MR_FileIO_Error_Handler

      use, intrinsic :: iso_fortran_env, only : &
         real32,real64,error_unit

      implicit none
      !implicit none (type, external)

        ! These single and double precision parameters should be 4 and 8
      integer, parameter :: sp = real32  ! selected_real_kind( 6,   37) ! single precision
      integer, parameter :: dp = real64  ! selected_real_kind(15,  307) ! double precision

      integer             :: nargs
      integer             :: iostatus
      integer             :: inlen
      character(len=120)  :: iomessage
      character(len=130)  :: arg

      character(len=100)  :: infile
      integer             :: intstep
      integer             :: inLLflag
      real(kind=sp)        :: inlon,inlat
      character           :: intrunc
      integer             :: invars
      integer,dimension(:)    ,allocatable :: invarlist
      integer             :: iw,iwf,idf,igrid,iwfiles
      integer             :: inyear,inmonth,inday
      real(kind=dp)        :: inhour

      integer             :: nxmax,nymax,nzmax
      real(kind=sp),dimension(:)    ,allocatable :: lon_grid
      real(kind=sp),dimension(:)    ,allocatable :: lat_grid
      real(kind=sp),dimension(:)    ,allocatable :: z_cc
      logical             :: IsPeriodic

      integer :: BaseYear = 1900
      logical :: useLeap  = .true.
      logical :: IsThere
      logical :: Truncate = .false.
      integer :: i
      integer :: il
      integer :: dum_i
      real(kind=dp) :: steptime

      integer :: io                           ! Index for output streams
      character(len=80)  :: linebuffer080

      INTERFACE
        subroutine Print_Usage
          implicit none
          !implicit none (type, external)
        end subroutine Print_Usage
        real(kind=8) function HS_HourOfDay(HoursSince,byear,useLeaps)
          implicit none
          !implicit none (type, external)
          integer        ,parameter   :: dp        = 8 ! double precision
          real(kind=dp)  ,intent(in)  :: HoursSince
          integer        ,intent(in)  :: byear
          logical        ,intent(in)  :: useLeaps
        end function HS_HourOfDay
        integer function HS_YearOfEvent(HoursSince,byear,useLeaps)
          implicit none
          !implicit none (type, external)
          integer        ,parameter   :: dp        = 8 ! double precision
          real(kind=dp)  ,intent(in)  :: HoursSince
          integer        ,intent(in)  :: byear
          logical        ,intent(in)  :: useLeaps
        end function HS_YearOfEvent
        integer function HS_MonthOfEvent(HoursSince,byear,useLeaps)
          implicit none
          !implicit none (type, external)
          integer        ,parameter   :: dp        = 8 ! double precision
          real(kind=dp)  ,intent(in)  :: HoursSince
          integer        ,intent(in)  :: byear
          logical        ,intent(in)  :: useLeaps
        end function HS_MonthOfEvent
        integer function HS_DayOfEvent(HoursSince,byear,useLeaps)
          implicit none
          !implicit none (type, external)
          integer        ,parameter   :: dp        = 8 ! double precision
          real(kind=dp)  ,intent(in)  :: HoursSince
          integer        ,intent(in)  :: byear
          logical        ,intent(in)  :: useLeaps
        end function HS_DayOfEvent
        integer function HS_DayOfYear(HoursSince,byear,useLeaps)
          implicit none
          !implicit none (type, external)
          integer        ,parameter   :: dp        = 8 ! double precision
          real(kind=dp)  ,intent(in)  :: HoursSince
          integer        ,intent(in)  :: byear
          logical        ,intent(in)  :: useLeaps
        end function HS_DayOfYear

        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          implicit none
          !implicit none (type, external)
          integer        ,parameter   :: dp        = 8 ! double precision
          integer        ,intent(in)  :: iyear
          integer        ,intent(in)  :: imonth
          integer        ,intent(in)  :: iday
          real(kind=dp)  ,intent(in)  :: hours
          integer        ,intent(in)  :: byear
          logical        ,intent(in)  :: useLeaps
        end function HS_hours_since_baseyear

        subroutine GetMetProfile(invars,invarlist)
          implicit none
          !implicit none (type, external)
          integer                  ,intent(in) :: invars
          integer,dimension(invars),intent(in) :: invarlist
        end subroutine GetMetProfile
      END INTERFACE

      ! Make sure user MetReader is using the same calendar
      MR_BaseYear = BaseYear
      MR_useLeap  = useLeap
      MR_useCompH = .false.

!     TEST READ COMMAND LINE ARGUMENTS
      nargs = command_argument_count()
      if (nargs < 9) then
        ! Not enough arguments; print usage and exit
        call Print_Usage
      else
        ! Get file name or windroot for iw=5
        call get_command_argument(1, arg, length=inlen, status=iostatus)
        if(iostatus /= 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read first command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        infile = trim(adjustl(arg))
        ! Error-check infile
        inquire( file=trim(adjustl(infile)), exist=IsThere )
        if(.not.IsThere)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR: Cannot find input file"
            write(errlog(io),'(a)')infile
          endif;enddo
          stop 1
        endif
        do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
          write(outlog(io),*)"infile = ",infile
        endif;enddo
        ! Get time step to use
        call get_command_argument(2, arg, length=inlen, status=iostatus)
        if(iostatus /= 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read second command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)intstep
        linebuffer080 = arg(1:80)
        if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check intstep
        if(intstep < 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR: Time-step should be a positive integer."
            write(errlog(io),*)" intstep = ",intstep
          endif;enddo
          stop 1
        else
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)"Time step of windfile = ",intstep
          endif;enddo
        endif

        ! Get Lat/Lon flag
        call get_command_argument(3, arg, length=inlen, status=iostatus)
        if(iostatus /= 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read third command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)inLLflag
        linebuffer080 = arg(1:80)
        if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check inLLflag
        if(inLLflag == 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)"Lat/Lon flag not set."
            write(outlog(io),*)"Using native grid of windfile (which may be LL)."
          endif;enddo
        elseif(inLLflag == 1)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)"Lat/Lon flag is set."
            write(outlog(io),*)"Using a Lat/Lon grid regardless of the windfile."
          endif;enddo
        else
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR: inLLflag must be either 0 or 1"
          endif;enddo
          stop 1
        endif
        ! Get lon
        call get_command_argument(4, arg, length=inlen, status=iostatus)
        if(iostatus /= 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read fourth command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)inlon
        linebuffer080 = arg(1:80)
        if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check inlon
        if(inlon < -360.0_sp.or.&
           inlon > 360.0_sp)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)"MR ERROR: longitude not in range -180->360"
            write(errlog(io),*)" inlon = ",inlon
          endif;enddo
          stop 1
        endif
        ! round to the nearest third decimel
        inlon = nint(inlon * 1000.0_sp)*0.001_sp
        do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
          write(outlog(io),*)"inlon = ",inlon
        endif;enddo
        ! Get lat
        call get_command_argument(5, arg, length=inlen, status=iostatus)
        if(iostatus /= 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read fifth command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)inlat
        linebuffer080 = arg(1:80)
        if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check inlat
        if(inlat < -90.0_sp.or.&
           inlat > 90.0_sp)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)"MR ERROR: latitude not in range -90->90"
            write(errlog(io),*)" inlat = ",inlat
          endif;enddo
          stop 1
        endif
        ! round to the nearest third decimel
        inlat = nint(inlat * 1000.0_sp) *0.001_sp
        do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
          write(outlog(io),*)"inlat = ",inlat
        endif;enddo

        ! Get truncation flag
        call get_command_argument(6, arg, length=inlen, status=iostatus)
        if(iostatus /= 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read sixth command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)intrunc
        linebuffer080 = arg(1:80)
        if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! No error-checking on intrunc; check for T, false otherwise
        if(intrunc == 't'.or.intrunc == 'T')then
          Truncate = .true.
        else
          Truncate = .false.
        endif
        if(Truncate)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)"Truncation argument read as: ",intrunc
            write(outlog(io),*)"Truncate = ",Truncate, "output will be on met node"
          endif;enddo
        else
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)"Truncation argument read as: ",intrunc
            write(outlog(io),*)"Truncate = ",Truncate, &
                      "output will be on interpolated point"
          endif;enddo
        endif

        ! Get number of variables to read/export
        call get_command_argument(7, arg, length=inlen, status=iostatus)
        if(iostatus /= 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read seventh command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)invars
        linebuffer080 = arg(1:80)
        if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check invars
        if(invars < 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR: number of variables must be positive"
            write(errlog(io),*)" invars = ",invars
          endif;enddo
          stop 1
        endif
        if(invars > 10)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR: number of variables must be no more than 10"
            write(errlog(io),*)" invars = ",invars
          endif;enddo
          stop 1
        endif
        do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
          write(outlog(io),*)"invars = ",invars
          write(outlog(io),*)"Allocating var list of length ",invars
        endif;enddo
        allocate(invarlist(invars))

        ! Get variable IDs
        do i=1,invars
          dum_i = 7+i
          call get_command_argument(dum_i, arg, length=inlen, status=iostatus)
          if(iostatus /= 0)then
            do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
              write(errlog(io),*)"MR ERROR : Could not read 7+i command-line argument"
              write(errlog(io),*)" arg = ",arg
            endif;enddo
            call Print_Usage
          endif
          read(arg,*,iostat=iostatus,iomsg=iomessage)invarlist(i)
          linebuffer080 = arg(1:80)
          if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
          ! Error-check invarlist(i)
          if(invarlist(i) <= 0.or.invarlist(i) > MR_MAXVARS)then
            do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
              write(errlog(io),*)"MR ERROR: Variable ID must be in range 1-50"
              write(errlog(io),*)" i invarlist(i) = ",i,invarlist(i)
            endif;enddo
            stop 1
          endif
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)" Requested variable ID ",invarlist(i)
          endif;enddo
        enddo

        ! Get wind format
        call get_command_argument(8+invars, arg, length=inlen, status=iostatus)
        if(iostatus /= 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read 8+nvars command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)iw
        linebuffer080 = arg(1:80)
        if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check iw
        if(iw /= 1.and.&
           iw /= 2.and.&
           iw /= 3.and.&
           iw /= 4.and.&
           iw /= 5)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR: wind format must be one of 1,2,3,4, or 5"
            write(errlog(io),*)" iw = ",iw
          endif;enddo
          stop 1
        else
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)" iw = ",iw
          endif;enddo
        endif

        ! Get wind product
        call get_command_argument(9+invars, arg, length=inlen, status=iostatus)
        if(iostatus /= 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read 9+nvars command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)iwf
        linebuffer080 = arg(1:80)
        if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check iwf
        if(iwf /= 0.and.&
           iwf /= 1.and.&
           iwf /= 2.and.&
           iwf /= 3.and.&
           iwf /= 4.and.&
           iwf /= 5.and.&
           iwf /= 6.and.&
           iwf /= 7.and.&
           iwf /= 8.and.&
           iwf /= 10.and.&
           iwf /= 11.and.&
           iwf /= 12.and.&
           iwf /= 13.and.&
           iwf /= 20.and.&
           iwf /= 21.and.&
           iwf /= 22.and.&
           iwf /= 23.and.&
           iwf /= 24.and.&
           iwf /= 25.and.&
           iwf /= 26.and.&
           iwf /= 27.and.&
           iwf /= 28.and.&
           iwf /= 29.and.&
           iwf /= 30.and.&
           iwf /= 32.and.&
           iwf /= 33.and.&
           iwf /= 41.and.&
           iwf /= 42)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR: windformat not recognized"
            write(errlog(io),*)" iwf = ",iwf
          endif;enddo
          stop 1
        else
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)" iwf = ",iwf
          endif;enddo
        endif

        ! Get data format ID
        call get_command_argument(10+invars, arg, length=inlen, status=iostatus)
        if(iostatus /= 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR : Could not read 10+nvars command-line argument"
            write(errlog(io),*)" arg = ",arg
          endif;enddo
          call Print_Usage
        endif
        read(arg,*,iostat=iostatus,iomsg=iomessage)idf
        linebuffer080 = arg(1:80)
        if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        ! Error-check idf
        if(idf /= 1.and.&
           idf /= 2.and.&
           idf /= 3)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR: wind data format must be one of 1,2, or 3"
            write(errlog(io),*)"          where 1 = ASCII"
            write(errlog(io),*)"                2 = NetCDF"
            write(errlog(io),*)"                3 = grib"
            write(errlog(io),*)" idf = ",idf
          endif;enddo
          stop 1
        else
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            if(idf == 1)then
              write(outlog(io),*)" idf = ",idf, "ASCII"
            elseif(idf == 2)then
              write(outlog(io),*)" idf = ",idf, "NetCDF"
            elseif(idf == 3)then
              write(outlog(io),*)" idf = ",idf, "grib"
            endif
          endif;enddo
        endif

        do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
          write(outlog(io),*)"Expecting wind format, NWP product ID and data type:",iw,iwf,idf
        endif;enddo

        ! Get year, month, day, hour if provided (otherwise use intstep provided above)
        if(nargs >= 11+invars)then
          call get_command_argument(11+invars, arg, length=inlen, status=iostatus)
          if(iostatus /= 0)then
            do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
              write(errlog(io),*)"MR ERROR : Could not read 11+nvars command-line argument"
              write(errlog(io),*)" arg = ",arg
            endif;enddo
            call Print_Usage
          endif
          read(arg,*,iostat=iostatus,iomsg=iomessage)inyear
          linebuffer080 = arg(1:80)
          if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        else
          inyear=0
        endif
        if(nargs >= 12+invars)then
          call get_command_argument(12+invars, arg, length=inlen, status=iostatus)
          if(iostatus /= 0)then
            do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
              write(errlog(io),*)"MR ERROR : Could not read 12+nvars command-line argument"
              write(errlog(io),*)" arg = ",arg
            endif;enddo
            call Print_Usage
          endif
          read(arg,*,iostat=iostatus,iomsg=iomessage)inmonth
          linebuffer080 = arg(1:80)
          if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
          ! Error-check inmonth
          if(inmonth < 1.or.&
             inmonth > 12)then
            do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
              write(errlog(io),*)"MR ERROR: month must be in range 1-12"
              write(errlog(io),*)" inmonth = ",inmonth
            endif;enddo
            stop 1
          endif
        else
          inmonth=1
        endif
        if(nargs >= 13+invars)then
          call get_command_argument(13+invars, arg, length=inlen, status=iostatus)
          if(iostatus /= 0)then
            do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
              write(errlog(io),*)"MR ERROR : Could not read 13+nvars command-line argument"
              write(errlog(io),*)" arg = ",arg
            endif;enddo
            call Print_Usage
          endif
          read(arg,*,iostat=iostatus,iomsg=iomessage)inday
          linebuffer080 = arg(1:80)
          if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
          ! Error-check inday
          if(inday < 1.or.&
             inday > 31)then
            do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
              write(errlog(io),*)"MR ERROR: day must be in range 1-31"
              write(errlog(io),*)" inday = ",inday
            endif;enddo
            stop 1
          endif
        else
          inday=1
        endif
        if(nargs >= 14+invars)then
          call get_command_argument(14+invars, arg, length=inlen, status=iostatus)
          if(iostatus /= 0)then
            do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
              write(errlog(io),*)"MR ERROR : Could not read 14+nvars command-line argument"
              write(errlog(io),*)" arg = ",arg
            endif;enddo
            call Print_Usage
          endif
          read(arg,*,iostat=iostatus,iomsg=iomessage)inhour
          linebuffer080 = arg(1:80)
          if(iostatus /= 0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
          ! Error-check inhour
          if(inhour < 0.0_dp.or.&
             inhour > 24.0_dp)then
            do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
              write(errlog(io),*)"MR ERROR: hour must be in range 0.0-24.0"
              write(errlog(io),*)" inhour = ",inhour
            endif;enddo
            stop 1
          endif
        else
          inhour=0.0_dp
        endif
        if(inyear > 0)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)"Read date/hour as ",inyear,inmonth,inday,inhour
            write(outlog(io),*)"This will overwrite the requested timestep above"
            write(outlog(io),*)"only when using iw=5"
          endif;enddo
          MR_Comp_StartHour = HS_hours_since_baseyear(inyear,inmonth,inday,inhour,&
                                                MR_BaseYear,MR_useLeap)
        else
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)"Optional arguments for date/hour not provided"
            write(outlog(io),*)"The time step provided above will be used.",intstep
          endif;enddo
        endif

      endif

      if(Truncate)then
        do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
          write(outlog(io),*)"Truncating profile onto ",inlon,inlat
        endif;enddo
      else
        do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
          write(outlog(io),*)"Interpolating profile onto ",inlon,inlat
        endif;enddo
      endif

      igrid   = 0 !
      iwfiles = 1 ! single command-line argument or folder path for iw=5
      call MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)
      MR_windfiles(1) = trim(adjustl(infile))
      if(inyear > 0)then
        call MR_Read_Met_DimVars(inyear)
      else
        call MR_Read_Met_DimVars()
      endif
      ! The x and y Met grids are now read and loaded.

      ! Determine where inlon and inlat fall on this grid
      if(inLLflag == 0)then
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
        if(inlon < -360.0_sp)then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
            write(errlog(io),*)"MR ERROR: Longitude must be gt -360"
          endif;enddo
          stop 1
        endif
        if(inlon < 0.0_sp.or.inlon > 360.0_sp)inlon=mod(inlon+360.0_sp,360.0_sp)
      endif

      if(Truncate)then
        do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
          write(outlog(io),*)&
            "Truncating input coordinate to grid point on Met grid."
        endif;enddo
        do i=1,nx_fullmet-1
          if(x_fullmet_sp(i) <= inlon.and.x_fullmet_sp(i+1) > inlon)then
            inlon = x_fullmet_sp(i)
          endif
        enddo
        if(y_inverted)then
          do i=1,ny_fullmet-1
            if(y_fullmet_sp(i) > inlat.and.y_fullmet_sp(i+1) <= inlat)then
              inlat = y_fullmet_sp(i+1)
            endif
          enddo
        else
          do i=1,ny_fullmet-1
            if(y_fullmet_sp(i) <= inlat.and.y_fullmet_sp(i+1) > inlat)then
              inlat = y_fullmet_sp(i)
            endif
          enddo
        endif
        do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
          write(outlog(io),*)"Coordinate reset to ",inlon,inlat
        endif;enddo
      endif

      nxmax = 3 !
      nymax = 3 !
      nzmax = 2 ! This is not really used in this utility
      allocate(lon_grid(nxmax))
      allocate(lat_grid(nymax))

      if(IsLatLon_CompGrid)then
        lon_grid(1:3) = [inlon-0.5_sp,inlon,inlon+0.5_sp]
        lat_grid(1:3) = [inlat-0.5_sp,inlat,inlat+0.5_sp]
      else
        lon_grid(1:3) = [inlon-1.0_sp*abs(dx_met_const),inlon,inlon+1.0_sp*abs(dx_met_const)]
        lat_grid(1:3) = [inlat-1.0_sp*abs(dy_met_const),inlat,inlat+1.0_sp*abs(dy_met_const)]
      endif
      do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
        write(outlog(io),*)"lon grid = ",lon_grid
        write(outlog(io),*)"lat_grid = ",lat_grid
      endif;enddo

      allocate(z_cc(nzmax))    ; z_cc(1:2) = [0.0_sp, 10.0_sp]

      IsPeriodic = .false.
      call MR_Set_CompProjection(isLatLon_CompGrid,Comp_iprojflag, & ! if LL, only these matter
                                 Met_lam0,&
                                 Met_phi0,Met_phi1,Met_phi2,&
                                 Met_k0,Met_Re)

      do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
        write(outlog(io),*)"Setting up wind grids"
      endif;enddo
      call MR_Initialize_Met_Grids(nxmax,nymax,nzmax,&
                              lon_grid(1:nxmax), &
                              lat_grid(1:nymax), &
                              z_cc(1:nzmax)    , &
                              IsPeriodic)

      ! Sort out what to do about the time
      steptime = MR_windfile_starthour(1) + MR_windfile_stephour(1,intstep)

      call MR_Set_Met_Times(steptime,0.0_dp)

      do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
        write(outlog(io),*)"Checking full variable list for availability in this file."
      endif;enddo
      do i=1,50
        do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
          il = len(trim(adjustl(Met_var_NC_names(i))))
          if(il > 0)then
            write(outlog(io),*)i,Met_var_NC_names(i),Met_var_IsAvailable(i)
          else
            write(outlog(io),*)i,&
         "     --Variable unspecified--                                                   ",&
         Met_var_IsAvailable(i)
          endif
        endif;enddo
      enddo

      call GetMetProfile(invars,invarlist)

      call MR_Reset_Memory
      if(allocated(invarlist)) deallocate(invarlist)
      if(allocated(     z_cc)) deallocate(z_cc)
      if(allocated( lon_grid)) deallocate(lon_grid)
      if(allocated( lat_grid)) deallocate(lat_grid)

      do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
        write(outlog(io),*)"Program ended normally."
      endif;enddo

      end program MetProbe

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

      use MetReader,       only : &
         MR_nio,MR_VB,outlog,verbosity_info,&
         MR_dum3d_compP,MR_dum3d_compP_2,np_fullmet,Map_Case,MR_iMetStep_Now,&
         p_fullmet_sp,&
           MR_Read_HGT_arrays,&
           MR_Rotate_UV_GR2ER_Met,&
           MR_Read_3d_Met_Variable_to_CompP

      use, intrinsic :: iso_fortran_env, only : &
         real32,real64,error_unit

      implicit none
      !implicit none (type, external)

        ! These single and double precision parameters should be 4 and 8
      integer, parameter :: sp = real32  ! selected_real_kind( 6,   37) ! single precision
      !integer, parameter :: dp = real64  ! selected_real_kind(15,  307) ! double precision

      integer                  ,intent(in) :: invars
      integer,dimension(invars),intent(in) :: invarlist

      integer                              :: ivar,i,iv

      real(kind=sp),dimension(:,:),allocatable :: outvars
      real(kind=sp),dimension(:)  ,allocatable :: u,v
      !character        :: invarchar
      !character(len=8) :: frmtstr
      integer :: io                           ! Index for output streams

      allocate(outvars(invars,np_fullmet))
      ! We only need u and v for Map_Case = 5, but allocate anyway since
      ! these are small
      allocate(u(np_fullmet))
      allocate(v(np_fullmet))
      do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
        write(outlog(io),*)" Inside GetMetProfile"
      endif;enddo

      ! First load the Met grids for Geopotential
      MR_iMetStep_Now = 1 ! This is initialized to 0
      call MR_Read_HGT_arrays(MR_iMetStep_Now)

      do iv = 1,invars
        ivar = invarlist(iv)
        do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
          write(outlog(io),*)"Map_Case  and varID = ",Map_Case,iv
        endif;enddo
        if (Map_Case == 4.and. & ! Only need to do this if Met=proj and
                                 ! comp is LL
            (iv == 2.or.iv == 3))then
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)"MR Calling MR_Rotate_UV_GR2ER_Met to rotate vec"
          endif;enddo
          ! if a vector call, store u and v to rotate
          call MR_Rotate_UV_GR2ER_Met(MR_iMetStep_Now,.true.)
          u(1:np_fullmet) = MR_dum3d_compP(2,2,1:np_fullmet)
          v(1:np_fullmet) = MR_dum3d_compP_2(2,2,1:np_fullmet)
          if(iv == 2)outvars(iv,1:np_fullmet) =u(1:np_fullmet)
          if(iv == 3)outvars(iv,1:np_fullmet) =v(1:np_fullmet)
        else
          do io=1,MR_nio;if(MR_VB(io) <= verbosity_info)then
            write(outlog(io),*)"Calling MR_Read_3d_Met_Variable_to_CompP"
          endif;enddo
          call MR_Read_3d_Met_Variable_to_CompP(ivar,MR_iMetStep_Now)
          outvars(iv,1:np_fullmet) = MR_dum3d_compP(2,2,1:np_fullmet)
        endif

      enddo

      !write(invarchar,'(i1)')invars+1
      !frmtstr = '(' // invarchar // 'F16.5)'
      open(unit=20,file='NWP_prof.dat')
      do i = 1,np_fullmet
        !write(20,frmtstr)p_fullmet_sp(i),outvars(1:invars,i)
        write(20,107)p_fullmet_sp(i),outvars(1:invars,i)
      enddo
107   format(11F16.5)
      close(20)

      end subroutine GetMetProfile

!##############################################################################
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
         MR_nio,MR_VB,verbosity_error,errlog,MR_VB

      implicit none

      integer :: io

        do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
          write(errlog(io),*)"MR ERROR: insufficient command-line arguments"
          write(errlog(io),*)"Usage: MetProbe filename tstep llflag lon lat trunc nvars ",&
                                  "vars(nvars) iw iwf (inyear inmonth inday inhour)"
          write(errlog(io),*)"       file         : string  : name of input file"
          write(errlog(io),*)"       timestep     : integer : step in file"
          write(errlog(io),*)"       llflag       : integer : 0 for using windfile grid; 1 for"
          write(errlog(io),*)"                                  forcing Lat/Lon"
          write(errlog(io),*)"       lon/x        : real    : longitude (or x) of sonde point"
          write(errlog(io),*)"       lat/y        : real    : latitude (or y) of sonde point"
          write(errlog(io),*)"       trunc flag   : char    : truncation flag (T or F)"
          write(errlog(io),*)"                                  T = coordinates truncated to nearest met node"
          write(errlog(io),*)"                                  F = values interpolated onto given coordinate"
          write(errlog(io),*)"       nvars        : integer : number of pressure variables to export"
          write(errlog(io),*)"       varID(nvars) : integers: variable ID's to read and export"
          write(errlog(io),*)"                        1 = GPH (km)"
          write(errlog(io),*)"                        2 = U (m/s)"
          write(errlog(io),*)"                        3 = V (m/s)"
          write(errlog(io),*)"                        4 = W (m/s)"
          write(errlog(io),*)"                        5 = T (K)"
          write(errlog(io),*)"       iw           : integer : windfile format code (3,4,5)"
          write(errlog(io),*)"       iwf          : integer : windfile product code"
          write(errlog(io),*)"       idf          : integer : igrid (2 for nc, 3 for grib)"
          write(errlog(io),*)"       year         : integer : only needed for iw = 5 files"
          write(errlog(io),*)"       month        : integer : "
          write(errlog(io),*)"       day          : integer : "
          write(errlog(io),*)"       hour         : real    : "
          write(errlog(io),*)"  "
          write(errlog(io),*)"     To probe a 0.5-degree GFS file in nc format for GPH only, use:"
          write(errlog(io),*)"       MetProbe 2020061000_5.f006.nc 1 0 190.055 52.8222 T 1 1 4 20 2"
          write(errlog(io),*)"      same, but adding u,v,t"
          write(errlog(io),*)"       MetProbe 2020061000_5.f006.nc 1 0 190.055 52.8222 T 4 1 2 3 5 4 20 2"
          write(errlog(io),*)"     To probe a NAM 91 grid over AK in nc format with a LL coordinate, use:"
          write(errlog(io),*)"       MetProbe nam.tm06.grib2 1 1 190.055 52.8222 F 4 1 2 3 5 4 13 3"
          write(errlog(io),*)"      with a projected coordinate:"
          write(errlog(io),*)"       MetProbe nam.tm06.grib2 1 0 -1363.94 -3758.61 F 4 1 2 3 5 4 13 3"
          write(errlog(io),*)"     To  probe the NCEP 2.5-degree data (or other iw=5), we need the"
          write(errlog(io),*)"     full date"
          write(errlog(io),*)"       MetProbe 2020061000_5.f006.nc 1 1 190.055 52.8222 T 1 1 4 20 2 2018 1 1 0.0"
          write(errlog(io),*)"     Output is written to the file NWP_prof.dat"
        endif;enddo
        stop 1

      end subroutine Print_Usage

!##############################################################################

