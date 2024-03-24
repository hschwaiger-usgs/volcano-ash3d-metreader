!  Here is a script for downloading the Radiosonde data in Alaska
!#!/bin/bash
!
!StationNum=("70414" "70316" "70326" "70350" "70273")
!StationCode=("PASY" "PACD" "PAKN" "PADQ" "PANC")
!StationName=("Shemya Afb" "Cold Bay" "King Salmon" "Kodiak" "Anchorage")
!#yearmonthday=`date -u +%Y%m%d`                        #current year, month &
!day (e.g. 20110119)
!Y=`date -u +%Y`
!M=`date -u +%m`
!DD=`date -u +%d`
!Y=2015
!M=02
!DD=24
!HH=(00 12)
!
!cd /data/windfiles/MetProfiles
!
!for (( si=0;si<=4;si++))
!do
!  #TYPE=TEXT%3AUNMERGED
! echo
!"http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&YEAR=${Y}&MONTH=${M}&FROM=${DD}${HH[$1]}&TO=${DD}${HH[$1]}&STNM=${StationNum[si]}"
!> tmp.lnk
! wget -i tmp.lnk -O out.html
! cp out.html RAW/${StationCode[si]}_${Y}${M}${DD}${HH[$1]}_raw.dat
!
! sed '1,10d' out.html > headless.dat
! sed '/PRE/,$d' headless.dat > ${StationCode[si]}_${Y}${M}${DD}${HH[$1]}.dat
! rm tmp.lnk headless.dat out.html
!
!done
!
!PRES:   Atmospheric Pressure    [hPa]
!HGHT:   Geopotential Height     [meter]
!TEMP:   Temperature     [Celsius]
!DWPT:   Dewpoint Temperature    [Celsius]
!RELH:   Relative Humidity       [%]
!MIXR:   Mixing Ratio    [gram/kilogram]
!DRCT:   Wind Direction  [degrees true]
!SKNT:   Wind Speed      [knot]
!THTA:   Potential Temperature   [kelvin]
!THTE:   Equivalent Potential Temperature        [kelvin]
!THTV:   Virtual Potential Temperature   [kelvin]


!##############################################################################
!
!     MR_Read_Met_DimVars_ASCII_1d
!
!     Called once from MR_Read_Met_DimVars
!
!     This subroutine needs to determine the full Met grid (i.e. determine the
!     pressure grid and the spatial grid).  It also needs to determine which
!     variables are available.  Generally, the actual region used (e.g. the submet)
!     is determined through MR_Initialize_Met_Grids after the comp grid is
!     defined.  Values are normally only read on this submet grid.  For 1d
!     ASCII cases, all sonde data are loaded into memory in this subroutine.
!     MR_Initialize_Met_Grids will be used to define the mapping of sonde data
!     to the comp grid.
!
!     This subroutine also allocates and populates the time arrays for the sonde
!     data, something usually done by MR_Read_Met_Times_[netcdf,grib,etc]

!     Note: Every sonde profile will have its own pressure levels, however,
!           the 1d array, p_fullmet_sp, is allocated as the intersection of all
!           the individual pressure profiles so that programs can have data on the
!           metP grid
!     
!     Allocated the dummy arrays for storing met data on met and computational
!     grids:  MR_dum2d_met(nx_submet,ny_submet)
!             MR_dum3d_metP(nx_submet,ny_submet,np_fullmet)
!             MR_dum3d_metH(nx_submet,ny_submet,nz_comp)
!             MR_dum2d_comp(nx_comp,ny_comp)
!             MR_dum3d_compH(nx_comp,ny_comp,nz_comp)
!             MR_windfile_starthour(MR_Snd_nt_fullmet)
!             MR_windfile_stephour(MR_Snd_nt_fullmet,1)
!
!##############################################################################

      subroutine MR_Read_Met_DimVars_ASCII_1d

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_error,verbosity_info,verbosity_production,&
         x_fullmet_sp,y_fullmet_sp,p_fullmet_sp,Snd_idx,MR_nSnd_Locs,MR_iwindfiles,MR_windfile_stephour,&
         MR_windfiles_nt_fullmet,MR_windfile_starthour,Met_dim_IsAvailable,Met_var_IsAvailable,&
         MR_SndVars_metP,MR_SndVarsID,MR_Snd_np_fullmet,np_fullmet,z_approx,MR_windfiles,&
         Snd_Have_PT,nt_fullmet,MR_BaseYear,MR_useLeap,MR_Snd_nt_fullmet,MR_Snd_nvars,Snd_Have_Coord,&
         MR_Max_geoH_metP_predicted,MR_iwind,MR_iwindformat,IsGlobal_MetGrid,IsLatLon_MetGrid,IsRegular_MetGrid,&
         Met_iprojflag,Met_k0,Met_lam0,Met_lam1,Met_lam2,Met_phi0,Met_phi1,Met_phi2,Met_Re,MR_EPS_SMALL,&
         MR_Comp_StartYear,MR_Comp_StartMonth,MR_nstat,&
           MR_Z_US_StdAtm,&
           MR_Temp_US_StdAtm,&
           MR_Pres_US_StdAtm, &
           MR_FileIO_Error_Handler

      use projection,      only : &
         PJ_iprojflag,PJ_k0,PJ_lam0,PJ_lam1,PJ_lam2,PJ_phi0,PJ_phi1,PJ_phi2,PJ_Re,&
           PJ_Set_Proj_Params

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision
      real(kind=dp), parameter :: PI        = 3.141592653589793_dp
      real(kind=dp), parameter :: DEG2RAD   = 1.7453292519943295e-2_dp
      integer, parameter :: MAX_ROWS  = 300 ! maximum number of row of data
      integer, parameter :: fid       = 120

      integer :: iostatus
      integer :: ioerr
      character(len=120) :: iomessage
      integer :: istr1,istr2,istr3
      integer :: ic,il,iil,iv

      integer :: nlev,ulev
      integer :: iw_idx
      integer :: iloc, itime
      integer :: p_lidx, p_tidx
      real(kind=sp) :: p_maxtop, p_top
      integer :: k

      real(kind=sp),dimension(:),allocatable :: WindVelocity
      real(kind=sp),dimension(:),allocatable :: WindDirection
      real(kind=sp),dimension(16) :: pres_Snd_tmp
      real(kind=sp) :: WindTime
      integer :: iw,iws
      real(kind=sp) :: rvalue1,rvalue2,rvalue3 
      real(kind=sp),dimension(10) :: rvalues
      integer       :: ivalue1,ivalue2,ivalue3,ivalue4,ivalue5
      real(kind=sp) :: rvalue1_o,rvalue3_o
      integer       :: ivalue2_o,ivalue4_o,ivalue5_o
      integer       :: ilatlonflag

      character(len=80)  :: linebuffer080
      character(len=80)  :: linebuffer080_2
      character(len=80)  :: linebuffer080_3
      character(len=80)  :: linebuffer080_4
      character(len=7)   :: field_str
      character(len=6),dimension(53) :: GTSstr
      character(len=6)   :: dumstr1,dumstr2,dumstr3,dumstr4
      integer :: dum_int
      integer,dimension(:),allocatable :: H_tmp
      integer,dimension(:),allocatable :: T_tmp
      integer :: il_bot,il_top,ill
      real(kind=sp) :: vbot,vtop,zbot,ztop,zhere
      integer :: scl_idx
      real(kind=sp),dimension(7) :: scl_m,scl_a
      real(kind=sp) :: SurfPres,SurfTemp,SurfDewPoint,SurfWindDir,SurfWindSpeed
      integer :: SurfTemp_int
      logical :: In_hPa = .true.
      integer :: indx1,indx2,indx3
      integer :: ncols
      logical :: IsWindDirectSpeed
      logical :: IsCustVarOrder
      integer,dimension(:),allocatable :: SndColReadOrder
      character(len=3),dimension(0:50) :: MR_SndVarsName

      !real(kind=dp)      :: HS_hours_since_baseyear  ! function that calculates hours
                                                     !  since base year
      integer           :: Stat_ID
      integer           :: Stat_idx
      real(kind=sp)     :: Stat_elev
      integer           :: idx,idx2
      real(kind=sp)     :: fac
      integer           :: substr_pos1
      character(len=8)  :: date
      character(len=10) :: time2
      character(len=5)  :: zone
      integer           :: values(8)
      logical           :: IsGTS     = .false.
      logical           :: IsRUCNOAA = .false.
      logical           :: HasTTCC   = .false.
      integer           :: DayOfMonth, SndHour

      integer :: io                           ! Index for output streams

      INTERFACE
        subroutine MR_Load_Radiosonde_Station_Data
        end subroutine MR_Load_Radiosonde_Station_Data
        subroutine MR_Get_Radiosonde_Station_Coord(StatID,StatIdx,Stat_lon,Stat_lat,Stat_elv)
          integer,      parameter   :: sp        = 4 ! single precision
          integer,      intent(in)  :: StatID
          integer,      intent(out) :: StatIdx
          real(kind=sp),intent(out) :: Stat_lon,Stat_lat,Stat_elv
        end subroutine MR_Get_Radiosonde_Station_Coord
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer            :: iyear
          integer            :: imonth
          integer            :: iday
          real(kind=8)       :: hours
          integer            :: byear
          logical            :: useLeaps
        end function HS_hours_since_baseyear
      END INTERFACE

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------                MR_Read_Met_DimVars_ASCII_1d       ----------"
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

!------------------------------------------------------------------------------
!    MR_iwind.eq.1
!--------------------------------------------------------------------
!      MR_iwindformat.eq.1 (custom 1-d profile)
! L1 string header line
! L2 time(hr) nlev [ncol] [ivar(ncol)]
! L3 x/lon y/lat [LLflag] []
!  data
!    if no optional params, then
!      3 col -> alt(m) windspeed(m/s) winddir(E of N wind FROM)
!      5 col -> alt(m) windspeed(m/s) winddir(E of N wind FROM) pres(hPa) T(degC)
!      ncol  -> variables can be in any order but must include:
!                  1) either pres or alt or both (if one is given, the other is filled from stdAtm)
!                  2) either vx/vy (vector componets) or speed/dir (with 'dir' as wind from)
!                If T is not provided, it is auto-filled from StdAtm
!-------------
!Input wind file for Ashfall tests
!2    41                             #wind time, number of levels
!-122.18   46.20                     #Longitude, latitude of wind sounding (fictitious)
!0       10.000   90.00    101300.    15.0          z, windspeed, direction, p, T
!1000    10.000   90.00     89846.     8.5
!2000    10.000   90.00     79464.     2.0
!--------------------------------------------------------------------
!      MR_iwindformat.eq.2 (radiosonde format as downloaded from "http://weather.uwyo.edu)
!        This can either be in the Raw or WMO/GTS format, or the 'list' format.
!  Example WMO/GTS
! TTAA 56001 72694 99009 22868 29006 00137 20667 28507 92801 13662 25507 85505
!06857 21003 70077 01471 23023 50573 12971 24538 40739 25780 25040 30939 43564
!23045 25060 52563 23558 20200 61761 22548 15380 565// 24035 10637 567// 22531
!88185 62963 23051 77221 23563 40311 31313 58208 82326 51515 10164 00012 10194
!26506 22517
! TTBB 56008 72694 00009 22868 11999 20467 22979 18466 33806 02640 44802 02444
!55797 03459 66793 03658 77784 03256 88776 03456 99758 02057 11751 03070 22739
!02890 33717 01481 44713 01672 55702 01273 66698 01468 77665 00760 88660 00563
!99651 01161 11646 01165 22636 01170 33629 01365 44610 02363 55599 02965 66550
!07169 77523 10367 88483 15170 99456 18185 11427 21583 22372 29980 33304 42964
!44253 51964 55219 59158 66207 60960 77185 62963 88182 60965 99174 58373 11169
!59179 22167 59181 33166 589// 44164 579// 55160 591// 66155 567// 77150 565//
!88148 557// 99139 569// 11135 559// 22116 575// 33113 567// 44106 563// 55103
!577// 66102 567// 31313 58208 82326 41414 21701
! PPBB 56008 72694 90012 29006 27508 26008 90346 26007 28006 24509 90789 21517
!20519 22522 91246 25525 25027 25030 92057 24537 24540 25046 93017 23543 23045
!23563 939// 22548 9436/ 23554 24034 950// 23040
! TTCC 56001 72694 70864 571// 20017 50078 545// 34004 30407 497// 02010 20675
!471// 05510 10145 339// 07017 88999 77999 31313 58208 82326
!
!  If the string 'TTAA' is present, then this coded format is assumed, otherwise
!  the 'list' format is assumed.
!  The header is read line by line until a valid row of numbers is read with
!  columns expected in the order below.
!
!<HTML>
!<TITLE>University of Wyoming - Radiosonde Data</TITLE>
!<LINK REL="StyleSheet" HREF="/resources/select.css" TYPE="text/css">
!<BODY BGCOLOR="white">
!<H2>70273 PANC Anchorage Observations at 00Z 30 May 2017</H2>
!<PRE>
!-----------------------------------------------------------------------------
!   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
!    hPa     m      C      C      %    g/kg    deg   knot     K      K      K
!-----------------------------------------------------------------------------
! 1008.0     40    7.6    4.9     83   5.41    240      3  280.1  295.2  281.0
! 1003.0     91    5.8    3.8     87   5.03    243      4  278.7  292.7  279.6
! 1000.0    121    5.8    3.7     86   5.01    245      4  278.9  292.9  279.8
!:
!:
!:
!   13.3  29566  -38.5  -66.3      4   0.38     65     16  805.6  809.9  805.8
!   12.2  30175  -37.1  -66.2      3   0.42    110     23  831.2  836.1  831.4
!   11.7  30480  -36.4  -66.2      3   0.45    115     12  844.4  849.5  844.6
!   11.2  30785  -35.7  -66.1      3   0.47     90     20  857.7  863.2  857.9
!   10.7  31072  -35.1  -66.1      3   0.49                870.4  876.3  870.6
!</PRE><H3>Station information and sounding indices</H3><PRE>
!                         Station identifier: PANC
!                             Station number: 70273
!                           Observation time: 170530/0000
!                           Station latitude: 61.16
!                          Station longitude: -150.01
!                          Station elevation: 40.0
!                            Showalter index: 1.68
!                               Lifted index: 2.33
!    LIFT computed using virtual temperature: 2.32
!                                SWEAT index: 228.00
!                                    K index: 21.60
!                         Cross totals index: 28.10
!                      Vertical totals index: 29.80
!                        Totals totals index: 57.90
!      Convective Available Potential Energy: 0.00
!             CAPE using virtual temperature: 0.00
!                      Convective Inhibition: 0.00
!             CINS using virtual temperature: 0.00
!                     Bulk Richardson Number: 0.00
!          Bulk Richardson Number using CAPV: 0.00
!  Temp [K] of the Lifted Condensation Level: 274.14
!Pres [hPa] of the Lifted Condensation Level: 934.80
!     Mean mixed layer potential temperature: 279.49
!              Mean mixed layer mixing ratio: 4.44
!              1000 hPa to 500 hPa thickness: 5289.00
!Precipitable water [mm] for entire sounding: 11.68
!</PRE>

      write(linebuffer080  ,'(a80)') repeat(' ',80)
      write(linebuffer080_2,'(a80)') repeat(' ',80)
      write(linebuffer080_3,'(a80)') repeat(' ',80)
      write(linebuffer080_4,'(a80)') repeat(' ',80)

      MR_SndVarsName(:) = "   "
      MR_SndVarsName( 0) = "P  "
      MR_SndVarsName( 1) = "H  "
      MR_SndVarsName( 2) = "U  "
      MR_SndVarsName( 3) = "V  "
      MR_SndVarsName( 4) = "W  "
      MR_SndVarsName( 5) = "T  "
       ! Note that there is a deviation here from the master list: Met_var_NC_names
       !  Met_var:: 6 is for 3d pressure
       !  Met_var:: 7 is for pressure vert. velocity
      MR_SndVarsName( 6) = "Wsp"
      MR_SndVarsName( 7) = "Wdr"

      MR_SndVarsName(30) = "RH "
      MR_SndVarsName(31) = "SH "
      MR_SndVarsName(32) = "QL "
      MR_SndVarsName(33) = "QI "
      MR_SndVarsName(44) = "P0 "
      IsCustVarOrder    = .false.
      IsWindDirectSpeed = .true.

      Met_dim_IsAvailable = .false.
      Met_var_IsAvailable = .false.
      MR_nstat = min(MR_nSnd_Locs,MR_nstat)
      ! This is the start of a huge if statement that determines the type of ASCII input
      ! and loads the data.  All data for all time steps are loaded here into the
      ! variable MR_SndVars_metP(MR_nSnd_Locs,MR_Snd_nt_fullmet,MR_Snd_nvars,nrows)
      ! The if statement sorts into the following cases:
      !   1 : ascii sonde file with either 3, 5, or a custom number of columns of data
      !   2 : radio sonde data either the clear text variety or the WMO/GTS encoded
      if(MR_iwind.eq.1.and.MR_iwindformat.eq.1)then
        ! We are reading just one windfile with the following format
        !  If three values per line:
        !    Elevation(m) , Velocity(m/s) , Direction(degree E of N)
        !  If five values per line:
        !    Elevation(m) , Velocity(m/s) , Direction(degree E of N) , Pressure(hPa) ,
        !     Temperature(C)
        !  Or a custom list of variables if specified on line 2

        ! x and y fullmet array will just be the coordinates in the order listed
        allocate(x_fullmet_sp(MR_nSnd_Locs))
        allocate(y_fullmet_sp(MR_nSnd_Locs))
        allocate(Snd_idx(MR_nSnd_Locs))
        ! The number of time steps/file is fixed to 1
        nt_fullmet = 1
        ! There may be more windfiles than time steps if there are multiple sonde locations,
        ! but we will need a time stamp on each windfile
        allocate(MR_windfile_starthour(MR_iwindfiles))
        allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))

        MR_windfile_starthour = 0.0_dp ! This is the initialization; will be set below
        MR_windfile_stephour  = 0.0_dp ! This will be the final value since all files
                                       ! have one step and so no offset.
        MR_windfiles_nt_fullmet(:) = nt_fullmet

        do itime = 1,MR_Snd_nt_fullmet
          do iloc = 1,MR_nSnd_Locs
            iw_idx = (itime-1)*MR_nSnd_Locs + iloc
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Opening sonde file ",iw_idx,&
                                   trim(adjustl(MR_windfiles(iw_idx)))
            endif;enddo

            open(unit=fid,file=trim(adjustl(MR_windfiles(iw_idx))), status='old',action='read',err=1971)
            ! skip over first line (Comment line, maybe the location name)
            read(fid,*,iostat=iostatus,iomsg=iomessage)linebuffer080
            if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Reading L2: time(hr) nlev [ncol] [ivar(ncol)]
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
            if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
            ! Assume we can read at least two values (a real and an integer)
            read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) rvalue1, ivalue1
            if(iostatus.ne.0)then
              do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                if(iostatus.lt.0)then
                  write(errlog(io),*)'MR ERROR:  EOR encountered'
                else
                  write(errlog(io),*)'MR ERROR:  Error reading line 2 of ASCII windfile'
                  write(errlog(io),*)'           Expecting to read: rvalue1, ivalue1'
                  write(errlog(io),*)'           From the following line from the file: '
                  write(errlog(io),*)linebuffer080
                  write(errlog(io),*)'MR System Message: ',trim(adjustl(iomessage))
                endif
              endif;enddo
              stop 1
            endif
            WindTime = rvalue1
            MR_windfile_starthour(iw_idx) = WindTime
            nlev     = ivalue1
            ! Try for three values [ncol]
            read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) rvalue1, ivalue1, ivalue2
            if(iostatus.eq.0)then
              ! Success: third value is the number of variables
              !    Note: We need at least 5 variables for P,H,U,V,T
              ! First check if this is the first file read, otherwise do not allocate
              if(iw_idx.eq.1)then
                IsCustVarOrder = .true.
                MR_Snd_nvars = max(5,ivalue2)
                allocate(MR_SndVarsID(MR_Snd_nvars))    ! This is the storage oder
                MR_SndVarsID(1) = 0 ! P
                MR_SndVarsID(2) = 1 ! H
                MR_SndVarsID(3) = 2 ! U
                MR_SndVarsID(4) = 3 ! V
                MR_SndVarsID(5) = 5 ! T
                                    ! There might be more columns to map
                allocate(SndColReadOrder(MR_Snd_nvars)) ! This is the read oder
              endif
              read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) &
                    rvalue1,ivalue1, ivalue2, SndColReadOrder(1:MR_Snd_nvars)
              do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)" Parsing sonde info line as:"
                write(outlog(io),*)"  wind time: ",rvalue1
                write(outlog(io),*)"  nlevels  : ",ivalue1
                write(outlog(io),*)"  nvars    : ",ivalue2
                write(outlog(io),*)"  ivars(1:nvar) : ",SndColReadOrder(1:MR_Snd_nvars)
                write(outlog(io),*)"1-d ASCII file should contain ",ivalue2," columns"
              endif;enddo
              do iv = 1,MR_Snd_nvars
                do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                  write(outlog(io),*)" Column ",iv,MR_SndVarsName(SndColReadOrder(iv))
                endif;enddo
                if (SndColReadOrder(iv).eq.0)then
                  Met_dim_IsAvailable(2) = .true.  ! P
                elseif (SndColReadOrder(iv).eq.1)then
                  Met_var_IsAvailable(1) = .true.  ! GPH
                elseif (SndColReadOrder(iv).eq.2)then
                  Met_var_IsAvailable(2) = .true.  ! U
                elseif (SndColReadOrder(iv).eq.3)then
                  Met_var_IsAvailable(3) = .true.  ! V
                elseif (SndColReadOrder(iv).eq.4)then  
                  Met_var_IsAvailable(4) = .true.  ! W
                  Met_var_IsAvailable(7) = .true.  ! VVP
                elseif (SndColReadOrder(iv).eq.5)then
                  Met_var_IsAvailable(5) = .true.  ! T
                elseif (SndColReadOrder(iv).eq.6)then
                  Met_var_IsAvailable(8) = .true.  ! Wsp
                elseif (SndColReadOrder(iv).eq.7)then
                  Met_var_IsAvailable(9) = .true.  ! Wdr
                endif
              enddo
              ! Now checking how velocities are provided.  If coordinates are given, they take
              ! precedence and are used, with direction and speed ignored.  Otherwise, use
              ! direction and speed.
              if(Met_var_IsAvailable(2).and.Met_var_IsAvailable(3))then
                IsWindDirectSpeed = .false.
                if(Met_var_IsAvailable(2).and.Met_var_IsAvailable(8))then
                  do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                    write(outlog(io),*)"MR WARNING: Both U/V and Speed/Direction provided."
                    write(outlog(io),*)"            Ignoring Speed/Direction"
                  endif;enddo
                endif
              else
                IsWindDirectSpeed = .true.
              endif
              if (.not.(Met_dim_IsAvailable(2).or.Met_var_IsAvailable(1)))then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*)"MR ERROR:  No height variable given"
                endif;enddo
                stop 1
              elseif(.not.(Met_var_IsAvailable(2).and.Met_var_IsAvailable(3)).and. &
                     .not.(Met_var_IsAvailable(8).and.Met_var_IsAvailable(9)))then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*)"MR ERROR:  No U/V or Wsp/Wdr variable given"
                endif;enddo
                stop 1
              endif
            else ! end of branch where 3 values are read in
              ! If no list of variables is provided, we will still need a list up to 5
              if(iw_idx.eq.1)then
                MR_Snd_nvars = 5
                allocate(MR_SndVarsID(MR_Snd_nvars))
                ! This maps the column index of MR_SndVars_metP to ivar
                MR_SndVarsID(1) = 0 ! P
                MR_SndVarsID(2) = 1 ! H
                MR_SndVarsID(3) = 2 ! U
                MR_SndVarsID(4) = 3 ! V
                MR_SndVarsID(5) = 5 ! T
              endif
            endif
            ! Now we know what the columns will mean (3-col, 5-col, custom)
              ! This variable will hold all the sonde data, up to MAX_ROWS rows
              ! The order of the columns will be P, H, U, V, T + extras if given
            if(iw_idx.eq.1)then
              allocate(MR_SndVars_metP(MR_nSnd_Locs,MR_Snd_nt_fullmet,MR_Snd_nvars,MAX_ROWS))
              allocate(MR_Snd_np_fullmet(MR_nSnd_Locs,MR_Snd_nt_fullmet))
              MR_SndVars_metP   = 0.0_sp
              MR_Snd_np_fullmet = 0
            endif
            MR_Snd_np_fullmet(iloc,itime) = nlev
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Reading L3: x/lon y/lat [LLflag] [proj flags]
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
            if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
            read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)&
                   rvalue1, rvalue2  ! These are the coordinates of the sonde point
            if(iostatus.ne.0)then
              do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                if(iostatus.lt.0)then
                  write(errlog(io),*)'MR ERROR:  EOF encountered'
                else
                  write(errlog(io),*)'MR ERROR:  Error reading line 3 of ASCII windfile'
                  write(errlog(io),*)'           Expecting to read: rvalue1, rvalue2'
                  write(errlog(io),*)'           From the follwing line from the file: '
                  write(errlog(io),*)linebuffer080
                endif
              endif;enddo
              stop 1
            endif
            x_fullmet_sp(iloc) = rvalue1
            y_fullmet_sp(iloc) = rvalue2
            ! Try for three values
            read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) rvalue1,rvalue2, ivalue1
            if(iostatus.eq.0)then
              ! A third parameter is present: first value of projection line
              Snd_Have_Coord = .true.
              ilatlonflag = ivalue1
            else
              ! A third parameter is not present: assume met coordinates are Lon/Lat
              Snd_Have_Coord = .false.
              ilatlonflag = 1
              IsLatLon_MetGrid = .true.
            endif
            if(Snd_Have_Coord)then
              if(ilatlonflag.eq.1)then
                ! We are using a Lon/Lat grid.  No need to read anything
                ! else on this line.
                IsLatLon_MetGrid  = .true.
                IsGlobal_MetGrid  = .false.
                IsRegular_MetGrid = .false.
                ! The rest of this is ignored
                Met_iprojflag     = 4
                Met_lam0          =  265.0_8
                Met_phi0          =  25.0_8
                Met_phi1          =  25.0_8
                Met_phi2          =  25.0_8
                Met_k0            =  0.933_8
                Met_Re            =  6371.229_8
              else
                ! Grid in not Lon/Lat
                ! First try to read the projection code
                read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) rvalue1,rvalue2, &
                                                ivalue1, ivalue2
                if(iostatus.eq.0)then
                  Met_iprojflag = ivalue2
                  if(Met_iprojflag.ne.0)then
                    ! we have a geographic projection
                    ! Try to read the full projection line
                    ! We know the projection flag will be 0 since this is the branch for
                    ! that so we know ' 0 ' will be in the string, but the first or
                    ! second coordinate could have this string too.  Need a better way.

           ! HFS: compress repeated white spaces then look for the second space

                    indx1 = index(linebuffer080,' 0 ')
                    indx2 = index(linebuffer080,'#')   ! Check for the comment marker
                    indx3 = len_trim(linebuffer080)    ! get length of string
                    if(indx2.gt.0)then
                      linebuffer080_2(1:indx1+indx2) = linebuffer080(indx1+1:indx2-1)
                    else
                      linebuffer080_2(1:indx1+1+indx3) = linebuffer080(indx1+1:indx3)
                    endif
                    call PJ_Set_Proj_Params(linebuffer080_2)
                    IsLatLon_MetGrid  = .false.
                    IsGlobal_MetGrid  = .false.
                    IsRegular_MetGrid = .false.
                    Met_iprojflag = PJ_iprojflag
                    Met_k0        = PJ_k0
                    Met_Re        = PJ_Re
                    Met_lam0      = PJ_lam0
                    Met_lam1      = PJ_lam1
                    Met_lam2      = PJ_lam2
                    Met_phi0      = PJ_phi0
                    Met_phi1      = PJ_phi1
                    Met_phi2      = PJ_phi2
                  endif
                else
                  ! If we have a cartesian grid, but no projection code,
                  ! assume no geographic projection
                  Met_iprojflag=0
                endif
              endif
            else
              ! No projection info given, assume Lon/Lat
              IsLatLon_MetGrid = .true.
            endif
            ! Finished projection parameters,

            ! Allocating temporary space for wind data
            allocate( WindVelocity(nlev));  WindVelocity = 0.0_sp
            allocate(WindDirection(nlev)); WindDirection = 0.0_sp
             !Read elevation (m), speed, direction at each level
             ! Speed is given in m/s (multiply by 0.514444444 to convert
             ! from knots to m/s)
             ! Direction is given as degrees east of north and specifies the
             ! direction FROM which the wind is blowing.
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Reading L4-EOF: ncols of data in nlev rows
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do il=1,nlev
              read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
              if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
              ! For the first line, we need to determine the number of columns.
              if (il.eq.1)then
                ! Verify we can at least read 3
                ncols = 0
                read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)rvalues(1:3)
                if(iostatus.eq.0)then
                  ncols=3
                  ! Now try 4
                  read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)rvalues(1:4)
                  if(iostatus.eq.0)then
                    ncols=4
                    ! Now try 5
                    read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)rvalues(1:5)
                    if(iostatus.eq.0)then
                      ncols=5
                      ! Now try 6
                      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)rvalues(1:6)
                      if(iostatus.eq.0)then
                        ncols=6
                        ! Now try 7
                        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)rvalues(1:7)
                        if(iostatus.eq.0)then
                          ncols=7
                          ! Now try 8
                          read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)rvalues(1:8)
                          if(iostatus.eq.0)then
                            ncols=8
                            ! Now try 9
                            read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)rvalues(1:9)
                            if(iostatus.eq.0)then
                              ncols=9
                              ! Now try 10
                              read(linebuffer080,*,iostat=iostatus,iomsg=iomessage)rvalues(1:10)
                              if(iostatus.eq.0)then
                                ncols=10
                              endif
                            endif
                          endif
                        endif
                      endif
                    endif
                  endif
                else
                  do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                    write(errlog(io),*)&
                        "MR ERROR:  ASCII wind files must have at least 3 column of data."
                  endif;enddo
                  stop 1
                endif
              endif ! il.eq.1
              rvalues(:) = -1.99_sp
              ! read ncols of data on this row
              read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) rvalues(1:ncols)

              if(.not.IsCustVarOrder.and.ncols.eq.3)then
                ! If this is a 3-column file without the custom ordering, read as
                !   Altitude (m), wind speed (m/s), wind direction (deg E of N)
                IsWindDirectSpeed = .true.
                MR_SndVars_metP(iloc,itime,2,il) = rvalues(1)*1.0e-3_sp  ! convert to km
                WindVelocity(il)  = rvalues(2)
                WindDirection(il) = rvalues(3)
                ! Calculate pressure in Pa and T in K
                MR_SndVars_metP(iloc,itime,1,il) = MR_Pres_US_StdAtm(MR_SndVars_metP(iloc,itime,2,il))&
                                                   * 100.0_sp
                MR_SndVars_metP(iloc,itime,5,il) = MR_Temp_US_StdAtm(MR_SndVars_metP(iloc,itime,2,il))
              elseif(.not.IsCustVarOrder.and.ncols.eq.5)then
                ! If this is a 5-column file without the custom ordering, read as
                  ! Recall the first five columns are P,H,U,V,T
                IsWindDirectSpeed = .true.
                Snd_Have_PT = .true.
                ! Five values were successfully read, interpret as:
                !   Altitude (m), wind speed (m/s), wind direction (deg E of N)
                !   Pressure (hPa), Temperature (C)
                MR_SndVars_metP(iloc,itime,2,il) = rvalues(1)*1.0e-3_sp  ! convert to km
                WindVelocity(il)  = rvalues(2)
                WindDirection(il) = rvalues(3)
                if(iw_idx.eq.1.and.&         ! For the first file
                   il.eq.1.and.    &         ! and the first level
                   rvalues(4).gt.1500.0_sp)then ! test pressure value
                  ! Check if pressure is greater than the expected 1013 hPa.
                  ! If so, assume pressure is in Pa
                  In_hPa = .false.
                endif
                if(In_hPa)then
                  MR_SndVars_metP(iloc,itime,1,il) = rvalues(4) * 100.0_sp ! convert to Pa
                else
                  MR_SndVars_metP(iloc,itime,1,il) = rvalues(4)            ! already in Pa
                endif
                MR_SndVars_metP(iloc,itime,5,il) = rvalues(5) + 273.0_sp   ! convert to K

              elseif(IsCustVarOrder)then
                ! This is the customized list of data columns
                do ic=1,ncols
                  ! Unfortunately, we currently need to map the variable ID to the index
                  if(SndColReadOrder(ic).eq.0)then       ! pressure
                    if(iw_idx.eq.1.and.&         ! For the first file
                       il.eq.1.and.    &         ! and the first level
                       rvalues(ic).gt.1500.0_sp)then ! test pressure value
                      In_hPa = .false.
                    else
                      In_hPa = .true.
                    endif
                    if(In_hPa)then
                      MR_SndVars_metP(iloc,itime,1,il) = rvalues(ic) * 100.0_sp
                    else
                      MR_SndVars_metP(iloc,itime,1,il) = rvalues(ic)
                    endif
                  elseif(SndColReadOrder(ic).eq.1)then   ! Altitude (m)
                    MR_SndVars_metP(iloc,itime,2,il) = rvalues(ic)*1.0e-3_sp  ! convert to km
                  elseif(SndColReadOrder(ic).eq.2)then   ! U (m/s)
                    MR_SndVars_metP(iloc,itime,3,il)  = rvalues(ic)
                  elseif(SndColReadOrder(ic).eq.3)then   ! V (m/s)
                    MR_SndVars_metP(iloc,itime,4,il)  = rvalues(ic)
                  !elseif(SndColReadOrder(ic).eq.4)then   ! Vert velocity (m/s)
                  !  MR_SndVars_metP(iloc,itime,?,il)  = rvalues(i)
                  elseif(SndColReadOrder(ic).eq.5)then   ! Temperature (C)
                    MR_SndVars_metP(iloc,itime,5,il) = rvalues(ic) + 273.0_sp   ! convert to K
                  elseif(SndColReadOrder(ic).eq.6)then   ! Wind speed (m/s)
                    WindVelocity(il)  = rvalues(ic)
                  elseif(SndColReadOrder(ic).eq.7)then   ! Wind direction (from deg E of N)
                    WindDirection(il)  = rvalues(ic)
                  else
                    do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                      write(errlog(io),*)&
                        "MR WARNING: Ignoring data for variable ",MR_SndVarsName(SndColReadOrder(ic))
                    endif;enddo
                  endif
                enddo
              endif
              ! Finished interpreting the row of data

              if(IsWindDirectSpeed)then
                  ! If the wind data is provided as speed and direction, break it into U,V components
                MR_SndVars_metP(iloc,itime,3,il) = &
                  real(WindVelocity(il)*sin(pi + DEG2RAD*WindDirection(il)),kind=sp)
                MR_SndVars_metP(iloc,itime,4,il) = &
                  real(WindVelocity(il)*cos(pi + DEG2RAD*WindDirection(il)),kind=sp)
              else
                  ! Otherwise, calculate direction and speed from the components
                WindVelocity(il) = &
                  MR_SndVars_metP(iloc,itime,3,il)*MR_SndVars_metP(iloc,itime,3,il) + &
                  MR_SndVars_metP(iloc,itime,4,il)*MR_SndVars_metP(iloc,itime,4,il)
                WindVelocity(il) = sqrt(WindVelocity(il))
                WindDirection(il) = &
                  -real(atan2(MR_SndVars_metP(iloc,itime,4,il),MR_SndVars_metP(iloc,itime,3,il))/DEG2RAD,kind=sp) + &
                  90.0_sp
              endif
              !Met_dim_IsAvailable(1) = .true.  ! Time
              Met_dim_IsAvailable(2) = .true.  ! P
              !Met_dim_IsAvailable(3) = .true.  ! y
              !Met_dim_IsAvailable(4) = .true.  ! x
    
              Met_var_IsAvailable(1) = .true.  ! GPH
              Met_var_IsAvailable(2) = .true.  ! U
              Met_var_IsAvailable(3) = .true.  ! V
              Met_var_IsAvailable(5) = .true.  ! T
              if(il.eq.1)then
                do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                  write(outlog(io),*)"Sonde values loaded."
                  write(outlog(io),203)MR_SndVarsName(1:7)
                endif;enddo
              endif
              do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                write(outlog(io),204)il,real(MR_SndVars_metP(iloc,itime,1:5,il),kind=4),&
                                       WindVelocity(il),WindDirection(il)
              endif;enddo
 203          format(3x,'Level',6x,a3,'(Pa)', &
                                7x,a3,'(km)', &
                                3x,a3,'(m/s)', &
                                2x,a3,'(m/s)', &
                                1x,a3,'(K)', &
                                3x,a3,'(m/s)', &
                                2x,a3,'(degree E of N)')
 204          format(6x,i2,5x,7f10.3)
            enddo ! il=1,nlev
            close(fid)
            ! Finished reading all the data for this file

            deallocate( WindVelocity)
            deallocate(WindDirection)

          enddo ! iloc = 1,MR_nSnd_Locs
        enddo ! itime = 1,MR_Snd_nt_fullmet
        ! Finished reading all the data of all the files

        ! Here we look for the highest pressure value (lowest altitude) of all the pressure
        ! values to be used for setting the master pressure array (p_fullmet_sp)
        p_maxtop = 0.0_sp
        p_tidx = 0
        p_lidx = 0
        do itime = 1,MR_Snd_nt_fullmet
          do iloc = 1,MR_nSnd_Locs
            p_top =  MR_SndVars_metP(iloc,itime,1,MR_Snd_np_fullmet(iloc,itime))
            if(p_top.gt.p_maxtop)then
              p_maxtop = MR_Snd_np_fullmet(iloc,itime)
              p_tidx = itime
              p_lidx = iloc
            endif
          enddo
        enddo
        np_fullmet = MR_Snd_np_fullmet(p_lidx,p_tidx)
        !np_fullmet_Vz = np_fullmet
        !np_fullmet_RH = np_fullmet
        allocate(p_fullmet_sp(np_fullmet))
        !allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        !allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet)=MR_SndVars_metP(p_lidx,p_tidx,1,1:np_fullmet)
        !p_fullmet_Vz_sp = p_fullmet_sp
        !p_fullmet_RH_sp = p_fullmet_sp
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)/100.0_sp)
        allocate(z_approx(np_fullmet))
        do k=1,np_fullmet
          ! Calculate heights for US Std Atmos while pressures are still in mbars
          ! or hPa
          z_approx(k) = MR_Z_US_StdAtm(p_fullmet_sp(k))
        enddo
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elseif(MR_iwind.eq.1.and.MR_iwindformat.eq.2)then
        ! We are reading radiosonde data from http://weather.uwyo.edu/ or https://ruc.noaa.gov/raobs/
        ! First, these are always given in lon/lat, earth-relative coordinates:
        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .false.
        ! The rest of this is ignored
        Met_iprojflag     = 4
        Met_lam0          =  265.0_8
        Met_phi0          =  25.0_8
        Met_phi1          =  25.0_8
        Met_phi2          =  25.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

        ! Load the internal database of radiosonde stations
        call MR_Load_Radiosonde_Station_Data

        ! Allocate arrays for iGridCode points at iwindfiles times, and MAX_ROWS levels

        ! x and y fullmet array will just be the coordinates in the order listed
        allocate(x_fullmet_sp(MR_nSnd_Locs))
        allocate(y_fullmet_sp(MR_nSnd_Locs))
        allocate(Snd_idx(MR_nSnd_Locs))
        ! The number of time steps/file is fixed to 1
        nt_fullmet = 1
        ! There may be more windfiles than time steps if there are multiple sonde locations,
        ! but we will need a time stamp on each windfile
        allocate(MR_windfile_starthour(MR_iwindfiles))
        allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))

        MR_windfile_starthour = 0.0_dp ! This is the initialization; will be set below
        MR_windfile_stephour  = 0.0_dp ! This will be the final value since all files
                                       ! have one step and so no offset.
        MR_windfiles_nt_fullmet(:) = nt_fullmet

        do itime = 1,MR_Snd_nt_fullmet
          do iloc = 1,MR_nSnd_Locs
            iw_idx = (itime-1)*MR_nSnd_Locs + iloc
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Opening sonde file ",iw_idx,&
                                     trim(adjustl(MR_windfiles(iw_idx)))
            endif;enddo

            open(unit=fid,file=trim(adjustl(MR_windfiles(iw_idx))),status='old',action='read',err=1971)
            if(iw_idx.eq.1)then
              MR_Snd_nvars = 5
              allocate(MR_SndVarsID(MR_Snd_nvars))    ! This is the storage order
              MR_SndVarsID(1) = 0 ! P
              MR_SndVarsID(2) = 1 ! H
              MR_SndVarsID(3) = 2 ! U
              MR_SndVarsID(4) = 3 ! V
              MR_SndVarsID(5) = 5 ! T
              ! We only allocate 16 rows here because we will only read the TTAA and TTCC mandatory levels
              ! If you want to use the other information in the radiosonde, shape the data into iwindformat=1
              ! Note that this might be reduced from 16 if any of the data files do not extend up to 10 mb.
              nlev = 16
              allocate(MR_SndVars_metP(MR_nSnd_Locs,MR_Snd_nt_fullmet,MR_Snd_nvars,nlev))
              allocate(MR_Snd_np_fullmet(MR_nSnd_Locs,MR_Snd_nt_fullmet))
              allocate(H_tmp(nlev))
              allocate(T_tmp(nlev))
              MR_SndVars_metP = 0.0_sp
              MR_Snd_np_fullmet(1,1) = 0  ! We initialize this to 0 since we don't know how high
                                          ! the sonde went
              pres_Snd_tmp(1:nlev) = &
               (/1000.0_sp, 925.0_sp, 850.0_sp, 700.0_sp, 500.0_sp, &
                  400.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, &
                  100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp,  20.0_sp, &  ! Note: the TTCC part starts at 70 hPa
                   10.0_sp/)                                           !  Not all sondes go this high!
            endif
            ! Allocating temporary space for wind data
            allocate( WindVelocity(nlev));  WindVelocity = 0.0_sp
            allocate(WindDirection(nlev)); WindDirection = 0.0_sp

            ! First search the whole radiosonde file for the string 'TTAA'.  If this
            ! is present, then assume the file is in FAA604 (WMO/GTS or RAW) format
            idx   = 0
            iostatus = 0
            do while (iostatus.ge.0)
              read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
              ! No aborting on read-error here since we will read this file until the EOF
              idx =index(linebuffer080,"TTAA") ! Look for string identifying WMO/GTS or RAW (both have this)
              idx2=index(linebuffer080,"TTCC") ! Look for string indicating extended data
              if (idx.gt.0) then
                IsGTS = .true.
              endif
              if (idx2.gt.0) then
                HasTTCC = .true.
              endif
            enddo

            if(IsGTS)then
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !  Reading data GTS format
              !  This is either "Text: Raw" from http://weather.uwyo.edu/upperair/sounding.html
              !1. TTAA 56001 72694 99009 22868 29006 00137 20667 28507 92801 13662 25507 85505
              !  or "FAA604 format (WMO/GTS)" from https://ruc.noaa.gov/raobs/
              !1. SLEMANSLE
              !2. 72694 TTAA  56001 72694 99009 22868 29006 00137 20667 28507
              !  Note that the format is slightly different so we need to find out
              ! which it is
              GTSstr(1:53)="//////"  ! Note: each 'element' of GTSstr(:) is really a 6-char string
              rewind(fid)
              read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
              if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
              idx =index(linebuffer080,"TTAA")
              if(idx.gt.0)then
                ! 'TTAA' is in the first line, this is a file from Uni. Wyoming
                IsRUCNOAA = .false.
                ! the format for the TTAA block is 3 lines of 13 6-char strings
                read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080_2
                if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080_2,iomessage)
                read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080_3
                if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080_3,iomessage)
                 ! the dumstr accommodates the 'TTAA'
                read(linebuffer080  ,155,iostat=iostatus,iomsg=iomessage)dumstr1,GTSstr(1:12)
                if(iostatus.ne.0)then
                  do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                    if(iostatus.lt.0)then
                      write(errlog(io),*)'MR ERROR:  EOR encountered'
                    else
                      write(errlog(io),*)'MR ERROR:  Error reading GTS format sonde file'
                      write(errlog(io),*)'           Expecting to read: dumstr1,GTSstr(1:12)'
                      write(errlog(io),*)'           with format: 13a6'
                      write(errlog(io),*)'           From the following line from the file: '
                      write(errlog(io),*)linebuffer080
                      write(errlog(io),*)'MR System Message: '
                      write(errlog(io),*)iomessage
                    endif
                  endif;enddo
                  stop 1
                endif
                read(linebuffer080_2,155,iostat=iostatus,iomsg=iomessage)GTSstr(13:25)
                if(iostatus.ne.0)then
                  do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                    if(iostatus.lt.0)then
                      write(errlog(io),*)'MR ERROR:  EOR encountered'
                    else
                      write(errlog(io),*)'MR ERROR:  Error reading GTS format sonde file'
                      write(errlog(io),*)'           Expecting to read: GTSstr(13:25)'
                      write(errlog(io),*)'           with format: 13a6'
                      write(errlog(io),*)'           From the following line from the file: '
                      write(errlog(io),*)linebuffer080
                      write(errlog(io),*)'MR System Message: '
                      write(errlog(io),*)iomessage
                    endif
                  endif;enddo
                  stop 1
                endif
                read(linebuffer080_3,155,iostat=iostatus,iomsg=iomessage)GTSstr(26:38)
                if(iostatus.ne.0)then
                  do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                    if(iostatus.lt.0)then
                      write(errlog(io),*)'MR ERROR:  EOR encountered'
                    else
                      write(errlog(io),*)'MR ERROR:  Error reading GTS format sonde file'
                      write(errlog(io),*)'           Expecting to read: GTSstr(26:38)'
                      write(errlog(io),*)'           with format: 13a6'
                      write(errlog(io),*)'           From the following line from the file: '
                      write(errlog(io),*)linebuffer080
                      write(errlog(io),*)'MR System Message: '
                      write(errlog(io),*)iomessage
                    endif
                  endif;enddo
                  stop 1
                endif

                ! Now we need to load the TTCC string blocks if the sonde went high enough
                if(HasTTCC)then
                  idx2=0
                  do while (idx2.eq.0.and.iostatus.ge.0)
                    read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
                    ! No aborting on read-error here since we will read this file until there is an error
                    idx2=index(linebuffer080,"TTCC")
                    if (idx2.gt.0) then
                      exit
                    endif
                  enddo
                  read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080_2
                  if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080_2,iomessage)
                  read(linebuffer080  ,155,iostat=iostatus,iomsg=iomessage)&
                        dumstr1,dumstr2,dumstr3,GTSstr(39:48) ! the dumstr accommodate the first
                  if(iostatus.ne.0)then
                    do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                      if(iostatus.lt.0)then
                        write(errlog(io),*)'MR ERROR:  EOR encountered'
                      else
                        write(errlog(io),*)'MR ERROR:  Error reading GTS format sonde file'
                        write(errlog(io),*)'           Expecting to read: dumstr1,dumstr2,dumstr3,GTSstr(39:48)'
                        write(errlog(io),*)'           with format: 13a6'
                        write(errlog(io),*)'           From the following line from the file: '
                        write(errlog(io),*)linebuffer080
                        write(errlog(io),*)'MR System Message: '
                        write(errlog(io),*)iomessage
                      endif
                    endif;enddo
                    stop 1
                  endif

                  read(linebuffer080_2,155,iostat=iostatus,iomsg=iomessage)&
                        GTSstr(49:53)                         !   three blocks
                  if(iostatus.ne.0)then
                    do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                      if(iostatus.lt.0)then
                        write(errlog(io),*)'MR ERROR:  EOR encountered'
                      else
                        write(errlog(io),*)'MR ERROR:  Error reading GTS format sonde file'
                        write(errlog(io),*)'           Expecting to read: GTSstr(49:53)'
                        write(errlog(io),*)'           with format: 13a6'
                        write(errlog(io),*)'           From the following line from the file: '
                        write(errlog(io),*)linebuffer080_2
                        write(errlog(io),*)'MR System Message: '
                        write(errlog(io),*)iomessage
                      endif
                    endif;enddo
                    stop 1
                  endif
                endif

              else
                ! TTAA not in the first line, try the second
                read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
                idx =index(linebuffer080,"TTAA")
                if(idx.gt.0)then
                  ! this is a file from NOAA Rapid Update Cycle ruc.noaa.gov
                  IsRUCNOAA = .true.
                  ! the format for the TTAA block is 4 lines of 10 6-char strings
                  read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080_2
                  if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080_2,iomessage)
                  read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080_3
                  if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080_3,iomessage)
                  read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080_4
                  if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080_4,iomessage)

                                      ! here we need the space for the repeat StatID and TTAA
                  read(linebuffer080  ,154,iostat=iostatus,iomsg=iomessage)dumstr1,dumstr2,GTSstr(1:8)
                  if(iostatus.ne.0)then
                    do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                      if(iostatus.lt.0)then
                        write(errlog(io),*)'MR ERROR:  EOR encountered'
                      else
                        write(errlog(io),*)'MR ERROR:  Error reading NOAA Rapid Update Cycle format sonde file'
                        write(errlog(io),*)'           Expecting to read: dumstr1,dumstr2,GTSstr(1:8)'
                        write(errlog(io),*)'           with format: 10a6'
                        write(errlog(io),*)'           From the following line from the file: '
                        write(errlog(io),*)linebuffer080
                        write(errlog(io),*)'MR System Message: '
                        write(errlog(io),*)iomessage
                      endif
                    endif;enddo
                    stop 1
                  endif
                  read(linebuffer080_2,154,iostat=iostatus,iomsg=iomessage)GTSstr( 9:18)
                  if(iostatus.ne.0)then
                    do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                      if(iostatus.lt.0)then
                        write(errlog(io),*)'MR ERROR:  EOR encountered'
                      else
                        write(errlog(io),*)'MR ERROR:  Error reading NOAA Rapid Update Cycle format sonde file'
                        write(errlog(io),*)'           Expecting to read: GTSstr( 9:18)'
                        write(errlog(io),*)'           with format: 10a6'
                        write(errlog(io),*)'           From the following line from the file: '
                        write(errlog(io),*)linebuffer080_2
                        write(errlog(io),*)'MR System Message: '
                        write(errlog(io),*)iomessage
                      endif
                    endif;enddo
                    stop 1
                  endif
                  read(linebuffer080_3,154,iostat=iostatus,iomsg=iomessage)GTSstr(19:28)
                  if(iostatus.ne.0)then
                    do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                      if(iostatus.lt.0)then
                        write(errlog(io),*)'MR ERROR:  EOR encountered'
                      else
                        write(errlog(io),*)'MR ERROR:  Error reading NOAA Rapid Update Cycle format sonde file'
                        write(errlog(io),*)'           Expecting to read: GTSstr(19:28)'
                        write(errlog(io),*)'           with format: 10a6'
                        write(errlog(io),*)'           From the following line from the file: '
                        write(errlog(io),*)linebuffer080_3
                        write(errlog(io),*)'MR System Message: '
                        write(errlog(io),*)iomessage
                      endif
                    endif;enddo
                    stop 1
                  endif
                  read(linebuffer080_4,154,iostat=iostatus,iomsg=iomessage)GTSstr(29:38)
                  if(iostatus.ne.0)then
                    do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                      if(iostatus.lt.0)then
                        write(errlog(io),*)'MR ERROR:  EOR encountered'
                      else
                        write(errlog(io),*)'MR ERROR:  Error reading NOAA Rapid Update Cycle format sonde file'
                        write(errlog(io),*)'           Expecting to read: dumstr1,dumstr2,GTSstr(1:8)'
                        write(errlog(io),*)'           with format: 10a6'
                        write(errlog(io),*)'           From the following line from the file: '
                        write(errlog(io),*)linebuffer080_4
                        write(errlog(io),*)'MR System Message: '
                        write(errlog(io),*)iomessage
                      endif
                    endif;enddo
                    stop 1
                  endif

                  ! Now we need to load the TTCC string blocks if the sonde went high enough
                  if(HasTTCC)then
                    idx2=0
                    do while (idx2.eq.0.and.iostatus.ge.0)
                      read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
                      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
                      idx2=index(linebuffer080,"TTCC")
                      if (idx2.gt.0) then
                        exit
                      endif
                    enddo
                    read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080_2
                    if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080_2,iomessage)
                    read(linebuffer080  ,154,iostat=iostatus,iomsg=iomessage) &
                         dumstr1,dumstr2,dumstr3,dumstr4,GTSstr(39:44) 
                    if(iostatus.ne.0)then
                      do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                        if(iostatus.lt.0)then
                          write(errlog(io),*)'MR ERROR:  EOR encountered'
                        else
                          write(errlog(io),*)'MR ERROR:  Error reading NOAA Rapid Update Cycle format sonde file'
                          write(errlog(io),*)'           Expecting to read: dumstr1,dumstr2,dumstr3,dumstr4,GTSstr(39:44)'
                          write(errlog(io),*)'           with format: 10a6'
                          write(errlog(io),*)'           From the following line from the file: '
                          write(errlog(io),*)linebuffer080
                          write(errlog(io),*)'MR System Message: '
                          write(errlog(io),*)iomessage
                        endif
                      endif;enddo
                      stop 1
                    endif
                    read(linebuffer080_2,154,iostat=iostatus,iomsg=iomessage)GTSstr(45:53)                         
                    if(iostatus.ne.0)then
                      do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                        if(iostatus.lt.0)then
                          write(errlog(io),*)'MR ERROR:  EOR encountered'
                        else
                          write(errlog(io),*)'MR ERROR:  Error reading NOAA Rapid Update Cycle format sonde file'
                          write(errlog(io),*)'           Expecting to read: GTSstr(45:53)'
                          write(errlog(io),*)'           with format: 10a6'
                          write(errlog(io),*)'           From the following line from the file: '
                          write(errlog(io),*)linebuffer080_2
                          write(errlog(io),*)'MR System Message: '
                          write(errlog(io),*)iomessage
                        endif
                      endif;enddo
                      stop 1
                    endif
                  endif
                else
                  ! This file has a 'TTAA' in it, but is not one of the expected
                  ! formats (first or second line).  Abort.
                  do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                    write(errlog(io),*)"MR ERROR: Cannot determine the format of the"
                    write(errlog(io),*)"          radiosonde data."
                  endif;enddo
                  stop 1
                endif
              endif
              ! Now interpret the strings
              ! Block B: day/time
              read(GTSstr(1)(1:2),*,iostat=iostatus,iomsg=iomessage)DayOfMonth
              DayOfMonth = DayOfMonth-50
              read(GTSstr(1)(3:4),*,iostat=iostatus,iomsg=iomessage)SndHour
              ! There is no year in this file, so assume year for requested
              ! event.  Needs to be set by calling program prior to this point.
              if(MR_Comp_StartYear.lt.1900)then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  write(errlog(io),*)"MR ERROR: Trying to use sonde data in raw format"
                  write(errlog(io),*)"          with an unknown year."
                endif;enddo
                stop 1
              endif

              MR_windfile_starthour(iw_idx) = &
                 HS_hours_since_baseyear(MR_Comp_StartYear,MR_Comp_StartMonth,DayOfMonth,&
                                         real(SndHour,kind=dp),MR_BaseYear,MR_useLeap)

              ! Block A: station identifier
              read(GTSstr(2),*,iostat=iostatus,iomsg=iomessage)Stat_ID
              call MR_Get_Radiosonde_Station_Coord(Stat_ID, Stat_idx,&
                                       y_fullmet_sp(iloc),x_fullmet_sp(iloc),&
                                       Stat_elev)
              Snd_idx(iloc) = Stat_idx
              ! Block C: surface pressure
              if(GTSstr(3)(1:1).ne.'/')then
                read(GTSstr(3)(1:2),*,iostat=iostatus,iomsg=iomessage)dum_int ! should be 99 indicating surface
              else
                dum_int = -9999
              endif
              if(GTSstr(3)(3:3).ne.'/')then
                read(GTSstr(3)(3:5),*,iostat=iostatus,iomsg=iomessage)dum_int ! surface pressure in mb
              else
                dum_int = -9999
              endif
              SurfPres = real(dum_int,kind=sp)

              ! Block D: temperature/dew-point for last pressure
              if(GTSstr(4)(1:1).ne.'/')then
                read(GTSstr(4)(1:3),*,iostat=iostatus,iomsg=iomessage)SurfTemp_int  ! temperature is characters 1-3, dew-point 4-5
              else
                SurfTemp_int = -9999
              endif
              if (mod(SurfTemp_int,2).eq.0)then
                SurfTemp =  0.1_sp*real(SurfTemp_int,kind=sp)
              else
                SurfTemp = -0.1_sp*real(SurfTemp_int,kind=sp)
              endif
              SurfTemp = SurfTemp 
              if(GTSstr(4)(4:4).ne.'/')then
                read(GTSstr(4)(4:5),*,iostat=iostatus,iomsg=iomessage)dum_int
              else
                dum_int = -9999
              endif
              SurfDewPoint = real(dum_int,kind=sp)

              ! Block E: Wind direction and speed for last pressure
              if(GTSstr(5)(1:1).ne.'/')then
                read(GTSstr(5)(1:3),*,iostat=iostatus,iomsg=iomessage)dum_int  ! Wind direction is char 1-3
              else
                dum_int = -9999
              endif
              SurfWindDir = real(dum_int,kind=sp)
              if(GTSstr(5)(4:4).ne.'/')then
                read(GTSstr(5)(4:5),*,iostat=iostatus,iomsg=iomessage)dum_int  ! Wind speed (knts) is char 4-5
              else
                dum_int = -9999
              endif
              SurfWindSpeed = real(dum_int,kind=sp)

              if(HasTTCC)then
                ulev = nlev
              else
                ulev = 11
              endif
              MR_Snd_np_fullmet(iloc,itime) = ulev
              do il = 1,ulev
                istr1 = 6 + (il-1)*3    ! block for pressure/height
                istr2 = 6 + (il-1)*3 +1 ! block for temperature/dew-point
                istr3 = 6 + (il-1)*3 +2 ! block for wind direction/speed

                !  Block 1 of the 3-block level data (pressure, height)
                if(GTSstr(istr1)(1:1).ne.'/')then
                  read(GTSstr(istr1)(1:2),*,iostat=iostatus,iomsg=iomessage)dum_int ! this is a pressure indicator that we could double-check
                else
                  dum_int = -9999
                endif
                MR_SndVars_metP(iloc,itime,1,il) = pres_Snd_tmp(il)
                if(GTSstr(istr1)(3:3).ne.'/')then
                  read(GTSstr(istr1)(3:5),*,iostat=iostatus,iomsg=iomessage)H_tmp(il) ! height a.s.l
                else
                  H_tmp(il) = -9999
                endif
                !  Block 2 of the 3-block level data (temperature, dew point (not read))
                if(GTSstr(istr2)(1:1).ne.'/')then
                  read(GTSstr(istr2)(1:3),*,iostat=iostatus,iomsg=iomessage)T_tmp(il)  ! temperature in C is characters 1-3, dew-point 4-5
                else
                  T_tmp(il) = -9999
                endif
                !  Block 3 of the 3-block level data (wind direction, wind speed)
                if(GTSstr(istr3)(1:1).ne.'/')then
                  read(GTSstr(istr3)(1:3),*,iostat=iostatus,iomsg=iomessage)dum_int  ! Wind direction is char 1-3
                else
                  dum_int = -9999
                endif
                WindDirection(il) = real(dum_int,kind=sp)
                if(GTSstr(istr3)(4:4).ne.'/')then
                  read(GTSstr(istr3)(4:5),*,iostat=iostatus,iomsg=iomessage)dum_int  ! Wind speed (knts) is char 4-5
                else
                  dum_int = -9999
                endif
                WindVelocity(il) = real(dum_int,kind=sp)
              enddo

              !--------------------------------
              ! Now convert to the proper units.
              !  Pressure: No error-checking
              MR_SndVars_metP(iloc,itime,1,:) = MR_SndVars_metP(iloc,itime,1,:) * 100.0_sp ! convert to Pa
              ! Height:  also assume valid
              !  do the block-specific adjustments to get the correct m blocks for 1000 mb 
              if(H_tmp(1).ge.500)then   ! kludge for negative heights
                MR_SndVars_metP(iloc,itime,2,1)  = -1.0_sp * real(H_tmp(1)-500,kind=sp)
              else
                MR_SndVars_metP(iloc,itime,2,1)  = real(H_tmp(1),kind=sp)
              endif
              ! The heights are adjusted with a multiplicative and additive factor.  If the
              ! three-digit value is less then the previous, then we are in the next scaling
              ! bracket
              scl_idx = 1
              scl_m(1:7) = (/1.0_sp,     1.0_sp,   1.0_sp,   10.0_sp, &
                            10.0_sp,    10.0_sp,   10.0_sp/)
              scl_a(1:7) = (/0.0_sp,  1000.0_sp,  3000.0_sp,  0.0_sp, &
                         10000.0_sp, 20000.0_sp, 30000.0_sp/)
              do il = 2,ulev
                if (H_tmp(il).lt.H_tmp(il-1))then
                  ! check if we need to increment the scaling bracket
                  scl_idx = scl_idx + 1
                endif
                if (il.eq.5)scl_idx=4 ! for the 500 mb level, make sure we are using decimeters
                MR_SndVars_metP(iloc,itime,2,il) = scl_a(scl_idx) + scl_m(scl_idx)* &
                                                    real(H_tmp(il),kind=sp)
              enddo
              MR_SndVars_metP(iloc,itime,2,:) = MR_SndVars_metP(iloc,itime,2,:)*1.0e-3_sp  ! convert m to km

              ! Temperature: do error checking
              do il = 1,ulev
                if(T_tmp(il).ne.-9999)then
                  if (mod(dum_int,2).eq.0)then
                    MR_SndVars_metP(iloc,itime,5,il) =  0.1_sp*real(T_tmp(il),kind=sp)
                  else
                    MR_SndVars_metP(iloc,itime,5,il) = -0.1_sp*real(T_tmp(il),kind=sp)
                  endif
                else
                  MR_SndVars_metP(iloc,itime,5,il) = real(T_tmp(il),kind=sp)
                endif

!   HFS             write(*,*)MR_SndVars_metP(iloc,itime,:,il)
              enddo
              do il = 1,ulev
                ! Set NaN's below the surface to the surface value
                if(pres_Snd_tmp(il).gt.SurfPres.or.il.eq.1)then
                  if(T_tmp(il).eq.-9999) &
                     MR_SndVars_metP(iloc,itime,5,il) = SurfTemp
                else
                  ! For NaN's above the surface, interpolate to the next valid value
                  if(T_tmp(il).eq.-9999)then
                    il_bot=il-1
                    do ill=il,ulev
                      if(T_tmp(ill).ne.-9999)then
                        il_top = ill
                        exit
                      endif
                    enddo
                    zhere=MR_SndVars_metP(iloc,itime,2,il)
                    zbot=MR_SndVars_metP(iloc,itime,2,il_bot)
                    ztop=MR_SndVars_metP(iloc,itime,2,il_top)
                    vbot=MR_SndVars_metP(iloc,itime,5,il_bot)
                    vtop=MR_SndVars_metP(iloc,itime,5,il_top)
                    MR_SndVars_metP(iloc,itime,5,il) = vbot + (vtop-vbot)*(zhere-zbot)/(ztop-zbot)
                  endif
                endif
              enddo
              MR_SndVars_metP(iloc,itime,5,:) = MR_SndVars_metP(iloc,itime,5,:) + 273.0_sp   ! convert C to K

              ! Wind velocity: do error checking
              do il = 1,ulev
                ! Set NaN's below the surface to zero
                if(pres_Snd_tmp(il).gt.SurfPres.or.il.eq.1)then
                  if(WindVelocity(il).lt.-9998.0_sp) &
                     WindVelocity(il) = 0.0_sp
                else
                  ! For NaN's above the surface, interpolate to the next valid value
                  if(WindVelocity(il).lt.-9998.0_sp)then
                    il_bot=il-1
                    do ill=il,ulev
                      if(WindVelocity(ill).gt.-9998.0_sp)then
                        il_top = ill
                        exit
                      endif
                    enddo
                    zhere=MR_SndVars_metP(iloc,itime,2,il)
                    zbot=MR_SndVars_metP(iloc,itime,2,il_bot)
                    ztop=MR_SndVars_metP(iloc,itime,2,il_top)
                    vbot=MR_SndVars_metP(iloc,itime,5,il_bot)
                    vtop=MR_SndVars_metP(iloc,itime,5,il_top)
                    WindVelocity(il) = vbot + (vtop-vbot)*(zhere-zbot)/(ztop-zbot)
                  endif
                endif
              enddo
              WindVelocity(:)   = WindVelocity(:)*0.514444444_sp                             ! convert Knot to m/s

              ! Wind direction: do error checking
              do il = 1,ulev
                ! Set NaN's below the surface to zero
                if(pres_Snd_tmp(il).gt.SurfPres.or.il.eq.1)then
                  if(WindDirection(il).lt.-9998.0_sp) &
                     WindDirection(il) = 0.0_sp
                else
                  ! For NaN's above the surface, interpolate to the next valid value
                  if(WindDirection(il).lt.-9998.0_sp)then
                    il_bot=il-1
                    do ill=il,ulev
                      if(WindDirection(ill).gt.-9998.0_sp)then
                        il_top = ill
                        exit
                      endif
                    enddo
                    zhere=MR_SndVars_metP(iloc,itime,2,il)
                    zbot=MR_SndVars_metP(iloc,itime,2,il_bot)
                    ztop=MR_SndVars_metP(iloc,itime,2,il_top)
                    vbot=MR_SndVars_metP(iloc,itime,5,il_bot)
                    vtop=MR_SndVars_metP(iloc,itime,5,il_top)
                    WindDirection(il) = vbot + (vtop-vbot)*(zhere-zbot)/(ztop-zbot)
                  endif
                endif
              enddo

              MR_SndVars_metP(iloc,itime,3,:) = &
                real(WindVelocity(:)*sin(pi + DEG2RAD*WindDirection(:)),kind=sp)
              MR_SndVars_metP(iloc,itime,4,:) = &
                real(WindVelocity(:)*cos(pi + DEG2RAD*WindDirection(:)),kind=sp)

              do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)"==================================================================="
                write(outlog(io),140)
   140          format('Pressure (Pa)','    Height (km)', &
                       '        Vx (m/s)','        Vy (m/s)', &
                       '       Temp. (C)','  Wind Vel (m/s)','       Wind Dir.')
                do iil=1,16
                  write(outlog(io),141)MR_SndVars_metP(iloc,itime,1,iil),&
                            MR_SndVars_metP(iloc,itime,2,iil),&
                            MR_SndVars_metP(iloc,itime,3,iil),&
                            MR_SndVars_metP(iloc,itime,4,iil),&
                            MR_SndVars_metP(iloc,itime,5,iil),&
                            WindVelocity(iil),&
                            WindDirection(iil)
                enddo
   141          format(3x,7f10.3)
                write(outlog(io),*)"==================================================================="
              endif;enddo
            else  ! end of GTS format section
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !  Reading data Textlist format from http://weather.uwyo.edu/
              ! Now start reading the sonde file
              ! First we read the lines until we get a valid number
              rewind(fid)
              read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
              if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
              read(linebuffer080,150,iostat=iostatus,iomsg=iomessage) &
                   rvalue1, ivalue2, rvalue3, ivalue4, ivalue5
              do while (iostatus.ne.0)
                read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
                if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
                ! Sometimes the bottom pressure levels have only pressure and HGHT if p is below
                ! the ground surface. Formatted read of blank columns will fill the variable
                ! with a zero, unfortunately.  Here we initialize Fill_Value then only read
                ! those columns that have something other than 7-spaces.
                !   PRES   HGHT   TEMP   DRCT   SKNT
                !    hPa     m      C     deg   knot
                rvalue1 = -9999.0_sp
                ivalue2 = -9999
                rvalue3 = -9999.0_sp
                ivalue4 = -9999
                ivalue5 = -9999
                field_str(1:7) = linebuffer080(1:7)
                if(len(trim(adjustl(field_str))).gt.0)&
                  read(field_str,*,iostat=iostatus,iomsg=iomessage)rvalue1
                if(iostatus.eq.0)then
                  field_str(1:7) = linebuffer080(8:14)
                  if(len(trim(adjustl(field_str))).gt.0)&
                    read(field_str,*,iostat=ioerr,iomsg=iomessage)ivalue2
                  field_str(1:7) = linebuffer080(15:21)
                  if(len(trim(adjustl(field_str))).gt.0)&
                    read(field_str,*,iostat=ioerr,iomsg=iomessage)rvalue3
                  field_str(1:7) = linebuffer080(43:49)
                  if(len(trim(adjustl(field_str))).gt.0)&
                    read(field_str,*,iostat=ioerr,iomsg=iomessage)ivalue4
                  field_str(1:7) = linebuffer080(50:56)
                  if(len(trim(adjustl(field_str))).gt.0)&
                    read(field_str,*,iostat=ioerr,iomsg=iomessage)ivalue5
                endif
              enddo
              ! We should be at the first pressure level here.
              ! Initialize old values to the values we just read
              rvalue1_o = rvalue1
              ivalue2_o = ivalue2
              rvalue3_o = rvalue3 ! Might be 0.0/Fill_Value if pres level is lower surface
              ivalue4_o = ivalue4 !  ditto
              ivalue5_o = ivalue5 !  ditto

              ! The line buffer contains numeric values.  Assume the first value is pressure and
              ! compare with the mandatory levels
              ! Turns out some 'mandatory' levels might be absent, so we track the last level read
              ! and interpolater if needed. 
              il = 0  ! counter for the number of data lines read
              iil = 1 ! index for which mandatory pressure level we are currently at
              ! Plan to read up to MAX_ROWS
              do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)&
                   "==================================================================="
              endif;enddo
              do while (il.ne.MAX_ROWS.and. &  ! Assume there are no more than MAX_ROWS of data
                        iil.le.nlev)           ! Do not bother reading past 10 hPa (nlev=16)

                if(iil.ge.2)then
                  ! If the last manditory level has NaN's, copy this level to the last
                  if(MR_SndVars_metP(iloc,itime,1,iil-1).lt.-9990.0_sp)&
                     MR_SndVars_metP(iloc,itime,1,iil-1) = rvalue1
                  if(MR_SndVars_metP(iloc,itime,2,iil-1).lt.-9990.0_sp)&
                     MR_SndVars_metP(iloc,itime,2,iil-1) = real(ivalue2,kind=4)
                  if(MR_SndVars_metP(iloc,itime,3,iil-1).lt.-9990.0_sp)&
                     MR_SndVars_metP(iloc,itime,3,iil-1) = real(ivalue4,kind=4)
                  if(MR_SndVars_metP(iloc,itime,4,iil-1).lt.-9990.0_sp)&
                     MR_SndVars_metP(iloc,itime,4,iil-1) = real(ivalue5,kind=4)
                  if(MR_SndVars_metP(iloc,itime,5,iil-1).lt.-9990.0_sp)&
                     MR_SndVars_metP(iloc,itime,5,iil-1) = rvalue3
                endif

                ! Check if we've skipped the mandatory pressure level.  If so, interpolate
                ! Check if the pressure level we are looking for is greater than what we just read
                if(il.gt.1.and.pres_Snd_tmp(iil)-rvalue1.gt.MR_EPS_SMALL)then
                  ! The level we need is missing. Find the pressure step and interpolate
                  ! Note: we are doing a linear interpolation in pressure where it should be logarithmic
                  !       HFS Fix this
                  ! Should be log interpolation in pressure to get the best altitude, then linear
                  ! interpolation between altitude levels
                  fac = (rvalue1_o-pres_Snd_tmp(iil))/(rvalue1_o-rvalue1)
                  rvalue1 = pres_Snd_tmp(iil)                          ! PRES
                  ivalue2 = -9999
                  rvalue3 = -9999.0_sp
                  ivalue4 = -9999
                  ivalue5 = -9999
!                  ivalue2 = ivalue2_o - int(fac*(ivalue2_o - ivalue2)) ! HGHT
!                  rvalue3 = rvalue3_o -     fac*(rvalue3_o - rvalue3)  ! TEMP
!                  ivalue4 = ivalue4_o - int(fac*(ivalue4_o - ivalue4)) ! DRCT
!                  ivalue5 = ivalue5_o - int(fac*(ivalue5_o - ivalue5)) ! SKNT
                endif
                if (abs(pres_Snd_tmp(iil)-rvalue1).lt.1.0_sp) then ! Sometimes the last level is 10.5
                  ! found the next mandatory level
                  MR_SndVars_metP(iloc,itime,1,iil) = rvalue1
                  MR_SndVars_metP(iloc,itime,2,iil) = real(ivalue2,kind=4)
                  if(rvalue3.lt.-9990.0_sp)then
                     MR_SndVars_metP(iloc,itime,5,iil) = rvalue3_o
                  else
                    MR_SndVars_metP(iloc,itime,5,iil) = rvalue3
                  endif
                  WindVelocity(iil)                 = real(ivalue5,kind=4)
                  WindDirection(iil)                = real(ivalue4,kind=4)
                  ! Just put velocity and direction in the vx,vy slots for now until we clean the list
                  if(ivalue4.lt.-9990)then
                    MR_SndVars_metP(iloc,itime,3,iil) = real(ivalue4_o,kind=4)
                  else
                    MR_SndVars_metP(iloc,itime,3,iil) = real(ivalue4,kind=4)
                  endif

                  if(ivalue5.lt.-9990)then
                    MR_SndVars_metP(iloc,itime,4,iil) = real(ivalue5_o,kind=4)
                  else
                    MR_SndVars_metP(iloc,itime,4,iil) = real(ivalue5,kind=4)
                  endif
   
                  !do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                  !  write(outlog(io),'(7f13.5)')&
                  !          MR_SndVars_metP(iloc,itime,1,iil),& ! PRES
                  !          MR_SndVars_metP(iloc,itime,2,iil),& ! HGHT
                  !          MR_SndVars_metP(iloc,itime,3,iil),& ! TEMP
                  !          MR_SndVars_metP(iloc,itime,4,iil),& ! DRCT
                  !          MR_SndVars_metP(iloc,itime,5,iil),& ! SKNT
                  !          WindVelocity(iil),&
                  !          WindDirection(iil)
                  !endif;enddo
                  MR_Snd_np_fullmet(iloc,itime) = iil  ! This keeps getting reassigned with each
                                                       ! successful read
                  iil = iil + 1
                else
                  ! this is not a mandatory level; read the next line
                  rvalue1_o = rvalue1
                  ivalue2_o = ivalue2
                  rvalue3_o = rvalue3
                  ivalue4_o = ivalue4
                  ivalue5_o = ivalue5
                  read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
                  if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
                  il = il + 1
                  ! Check if we've read past the data and are in the Station footer
                  substr_pos1 = index(linebuffer080,'Station')
                  if(substr_pos1.ne.0)then
                    do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                      write(outlog(io),*)&
                         "MR WARNING: Read past the end of data without finding the last"
                      write(outlog(io),*)"             mandatory pressure level"
                    endif;enddo
                    nlev = iil-1
                    do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                      write(outlog(io),*)"             Resetting nlev to : ",nlev
                    endif;enddo
                    cycle
                  endif

                  rvalue1 = -9999.0_sp
                  ivalue2 = -9999
                  rvalue3 = -9999.0_sp
                  ivalue4 = -9999
                  ivalue5 = -9999
                  field_str(1:7) = linebuffer080(1:7)
                  if(len(trim(adjustl(field_str))).gt.0)&
                    read(field_str,*,iostat=iostatus,iomsg=iomessage)rvalue1
                  field_str(1:7) = linebuffer080(8:14)
                  if(len(trim(adjustl(field_str))).gt.0)&
                    read(field_str,*,iostat=iostatus,iomsg=iomessage)ivalue2
                  field_str(1:7) = linebuffer080(15:21)
                  if(len(trim(adjustl(field_str))).gt.0)&
                    read(field_str,*,iostat=iostatus,iomsg=iomessage)rvalue3
                  field_str(1:7) = linebuffer080(43:49)
                  if(len(trim(adjustl(field_str))).gt.0)&
                    read(field_str,*,iostat=iostatus,iomsg=iomessage)ivalue4
                  field_str(1:7) = linebuffer080(50:56)
                  if(len(trim(adjustl(field_str))).gt.0)&
                    read(field_str,*,iostat=iostatus,iomsg=iomessage)ivalue5
                  if(iostatus.ne.0)then
                    do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                      if(iostatus.lt.0)then
                        write(errlog(io),*)'MR ERROR:  EOR encountered'
                      else
                        write(errlog(io),*)'MR ERROR:  Error reading sonde in textlist format UWYO'
                        write(errlog(io),*)'           Expecting to read: rvalue1, ivalue2, rvalue3, ivalue4, ivalue5'
                        write(errlog(io),*)'           with format: 1x,f7.1,i7,f7.1,21x,i7,i7'
                        write(errlog(io),*)'           From the following line from the file: '
                        write(errlog(io),*)linebuffer080
                        write(errlog(io),*)'MR System Message: '
                        write(errlog(io),*)iomessage
                      endif
                    endif;enddo
                    stop 1
                  endif

                endif ! manditory level or not
              enddo  ! while (il.ne.MAX_ROWS.and.iil.le.nlev)
              ! Now we have read all the data from this file and hopefully filled
              ! any empty values

              ! Now convert to the expected units
              do iil=1,nlev
                  MR_SndVars_metP(iloc,itime,1,iil) = MR_SndVars_metP(iloc,itime,1,iil) * 100.0_sp  ! mb->Ppa
                  MR_SndVars_metP(iloc,itime,2,iil) = MR_SndVars_metP(iloc,itime,2,iil) *1.0e-3_sp  !  m->km
                  MR_SndVars_metP(iloc,itime,5,iil) = MR_SndVars_metP(iloc,itime,5,iil) + 273.0_sp  !  C->K
                  WindVelocity(iil)   = MR_SndVars_metP(iloc,itime,3,iil)*0.514444444_sp            ! knt->m/s
                  WindDirection(iil)  = MR_SndVars_metP(iloc,itime,4,iil)                           ! deg
                  MR_SndVars_metP(iloc,itime,3,iil) = &
                    real(WindVelocity(iil)*sin(pi + DEG2RAD*WindDirection(iil)),kind=sp)            ! Vx in m/s
                  MR_SndVars_metP(iloc,itime,4,iil) = &
                    real(WindVelocity(iil)*cos(pi + DEG2RAD*WindDirection(iil)),kind=sp)            ! Vy in m/s

                  do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                    write(outlog(io),'(7f13.5)')&
                            MR_SndVars_metP(iloc,itime,1,iil),& ! PRES Pa
                            MR_SndVars_metP(iloc,itime,2,iil),& ! HGHT km
                            MR_SndVars_metP(iloc,itime,3,iil),& ! Vx   m/s
                            MR_SndVars_metP(iloc,itime,4,iil),& ! Vy   m/s
                            MR_SndVars_metP(iloc,itime,5,iil),& ! TEMP C
                            WindVelocity(iil),&
                            WindDirection(iil)
                  endif;enddo
              enddo
              do io=1,MR_nio;if(VB(io).le.verbosity_info)then
                write(outlog(io),*)&
                   "==================================================================="
              endif;enddo
   150  format(1x,f7.1,i7,f7.1,21x,i7,i7)
   151  format(45x,i2,i2,i2,1x,i2)
   153  format(45x,i5)
   154  format(10a6)
   155  format(13a6)

              ! Now look for the station ID from the beginning of the file
              rewind(fid)
              read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
              if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
              do while (index(linebuffer080,"Station number").eq.0)
                read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
                if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
              enddo
              read(linebuffer080,153,iostat=iostatus,iomsg=iomessage)Stat_ID
              if(iostatus.ne.0)then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  if(iostatus.lt.0)then
                    write(errlog(io),*)'MR ERROR:  EOR encountered'
                  else
                    write(errlog(io),*)'MR ERROR:  Error reading sonde in textlist format UWYO'
                    write(errlog(io),*)'           Expecting to read: Stat_ID'
                    write(errlog(io),*)'           with format: 45x,i5'
                    write(errlog(io),*)'           From the following line from the file: '
                    write(errlog(io),*)linebuffer080
                    write(errlog(io),*)'MR System Message: '
                    write(errlog(io),*)iomessage
                  endif
                endif;enddo
                stop 1
              endif
              call MR_Get_Radiosonde_Station_Coord(Stat_ID, Stat_idx, &
                                       y_fullmet_sp(iloc),x_fullmet_sp(iloc),&
                                       Stat_elev)
              Snd_idx(iloc) = Stat_idx
              ! Now look for the observation time from the beginning of the file
              rewind(fid)
              read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
              if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
              do while (index(linebuffer080,"Observation time").eq.0)
                read(fid,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
                if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
              enddo
              read(linebuffer080,151,iostat=iostatus,iomsg=iomessage)ivalue1,ivalue2,ivalue3,ivalue4
              if(iostatus.ne.0)then
                do io=1,MR_nio;if(VB(io).le.verbosity_error)then
                  if(iostatus.lt.0)then
                    write(errlog(io),*)'MR ERROR:  EOR encountered'
                  else
                    write(errlog(io),*)'MR ERROR:  Error reading sonde in textlist format UWYO'
                    write(errlog(io),*)'           Expecting to read Observation time: ivalue1,ivalue2,ivalue3,ivalue4'
                    write(errlog(io),*)'           with format: 45x,i2,i2,i2,1x,i2'
                    write(errlog(io),*)'           From the following line from the file: '
                    write(errlog(io),*)linebuffer080
                    write(errlog(io),*)'MR System Message: '
                    write(errlog(io),*)iomessage
                  endif
                endif;enddo
                stop 1
              endif
              ! Need to find out if the 2-digit year is for the current century
              call date_and_time(date,time2,zone,values)
              if(2000+ivalue1.gt.values(1))then
                ivalue1=ivalue1+1900
              else
                ivalue1=ivalue1+2000
              endif
              MR_windfile_starthour(iw_idx) = HS_hours_since_baseyear(ivalue1,ivalue2,ivalue3,&
                                              real(ivalue4,kind=dp),MR_BaseYear,MR_useLeap)

            endif ! IsGTS
            close(fid)
            ! Finished reading all the data for this file

            deallocate( WindVelocity)
            deallocate(WindDirection)

          enddo
        enddo
        Met_dim_IsAvailable(2) = .true.  ! P
        Met_var_IsAvailable(1) = .true.  ! GPH
        Met_var_IsAvailable(2) = .true.  ! U
        Met_var_IsAvailable(3) = .true.  ! V
        Met_var_IsAvailable(5) = .true.  ! T
        Met_var_IsAvailable(8) = .true.  ! Wsp
        Met_var_IsAvailable(9) = .true.  ! Wdr

        ! Here we look for the highest pressure value (lowest altitude) of all the pressure
        ! values to be used for setting the master pressure array (p_fullmet_sp)
        p_maxtop = 0.0_sp
        p_tidx = 0
        p_lidx = 0
        do itime = 1,MR_Snd_nt_fullmet
          do iloc = 1,MR_nSnd_Locs
            p_top =  MR_SndVars_metP(iloc,itime,1,MR_Snd_np_fullmet(iloc,itime))
            if(p_top.gt.p_maxtop)then
              p_maxtop = MR_Snd_np_fullmet(iloc,itime)
              p_tidx = itime
              p_lidx = iloc
            endif
          enddo
        enddo
        np_fullmet = MR_Snd_np_fullmet(p_lidx,p_tidx)
        allocate(p_fullmet_sp(np_fullmet))
        p_fullmet_sp(1:np_fullmet)=MR_SndVars_metP(p_lidx,p_tidx,1,1:np_fullmet)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)/100.0_sp)
        allocate(z_approx(np_fullmet))
        do k=1,np_fullmet
          ! Calculate heights for US Std Atmos while pressures are still in mbars
          ! or hPa
          z_approx(k) = MR_Z_US_StdAtm(p_fullmet_sp(k))
        enddo
      else
        ! Neither MR_iwind.eq.1.and.MR_iwindformat.eq.1 nor 
        !         MR_iwind.eq.1.and.MR_iwindformat.eq.2
        ! Not sure what to do
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR:  MR_iwind=1, but MR_iwindformat is not 1 or 2."
        endif;enddo
        stop 1
      endif

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Sonde data read:"
        write(outlog(io),*)"  number of points     = ",MR_nSnd_Locs
        write(outlog(io),*)"  number of time steps = ",MR_Snd_nt_fullmet
      endif;enddo

      ! Finished setting up the start time of each wind file in HoursSince : MR_windfile_starthour(iw)
      !  and the forecast (offset from start of file) for each step        : MR_windfile_stephour(iw,iwstep)
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"File, step, Ref, Offset, HoursSince"
        do iw = 1,MR_iwindfiles
          do iws = 1,nt_fullmet
            write(outlog(io),'(2x,i3,5x,i3,3f13.2)')iw,iws,&
                             MR_windfile_starthour(iw),&
                             MR_windfile_stephour(iw,iws),&
                             MR_windfile_starthour(iw)+MR_windfile_stephour(iw,iws)
          enddo
        enddo
      endif;enddo

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)&
           "-----------------------------------------------------------------------"
      endif;enddo

      return

1971  do io=1,MR_nio;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error: cannot find file input wind file.',&
                    '  Program stopped.'
      endif;enddo
      stop 1

      end subroutine MR_Read_Met_DimVars_ASCII_1d

!##############################################################################
!
!     MR_Set_MetComp_Grids_ASCII_1d
!
!     Called once from MR_Initialize_Met_Grids
!
!     This subroutine sets the mapping between the network of sonde points and
!     the computational grid
!
!     Allocates the dummy arrays for storing met data on met and computational grids
!
!     Sets: 
!           n[t,x,y,p]_met     :: sets the size of the dimensions of the sub-met grid
!           [x,y,p]_met_sp     :: arrays holding dimension values of the sub-met grid
!           MR_dum2d_met(nx_submet,ny_submet)
!           MR_dum3d_metP(nx_submet,ny_submet,np_fullmet)
!           MR_dum3d_metH(nx_submet,ny_submet,nz_comp)
!           MR_dum2d_comp(nx_comp,ny_comp)
!           MR_dum3d_compH(nx_comp,ny_comp,nz_comp)
!           CompPoint_on_subMet_idx
!           bilin_map_wgt
!           CompPoint_X_on_Met_sp
!           CompPoint_Y_on_Met_sp
!
!##############################################################################

      subroutine MR_Set_MetComp_Grids_ASCII_1d

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_error,verbosity_production,&
         x_submet_sp,y_submet_sp,MR_dx_submet,MR_dy_submet,MR_dum3d_metP,MR_dum3d_metH,&
         MR_dum2d_comp_int,MR_dum2d_comp,MR_dum3d_compH,MR_dum3d_compP,&
         MR_geoH_metP_last,MR_geoH_metP_next,MR_Snd2Comp_map_wgt,MR_Snd2Comp_map_idx,&
         MR_dum3d_metP,nx_submet,ny_submet,nx_comp,ny_comp,nz_comp,MR_iwind,MR_iwindformat,&
         MR_pexp,MR_nstat,MR_nSnd_Locs,np_fullmet

      use projection,      only : &
           PJ_Set_Proj_Params

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer :: io                           ! Index for output streams

      INTERFACE
        subroutine MR_Get_Radiosonde_Stations_InDomain
        end subroutine MR_Get_Radiosonde_Stations_InDomain
        subroutine MR_Calc_Snd_Weights_NearNeigh
        end subroutine MR_Calc_Snd_Weights_NearNeigh
        subroutine MR_Calc_Snd_Weights_InvDistWeight(pexp,nstat)
          integer,      parameter   :: sp        = 4 ! single precision
          real(kind=sp), intent(in) :: pexp
          integer,       intent(in), optional :: nstat
        end subroutine MR_Calc_Snd_Weights_InvDistWeight
      END INTERFACE

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)"-----------------------------------------------------------------------"
        write(outlog(io),*)"----------                          MR_Set_MetComp_Grids_ASCII_1d -----"
        write(outlog(io),*)"-----------------------------------------------------------------------"
      endif;enddo

      if(MR_iwind.eq.1.and.MR_iwindformat.eq.1.or.&
         MR_iwind.eq.1.and.MR_iwindformat.eq.2)then
          ! 1d ascii file or Radiosonde case
        nx_submet     = MR_nSnd_Locs ! put all sonde data in x-dimension
        ny_submet     = 1
        allocate(x_submet_sp(1:nx_submet));  x_submet_sp(:)  = 0.0_sp
        allocate(y_submet_sp(1:ny_submet));  y_submet_sp(:)  = 0.0_sp
        allocate(MR_dx_submet(1:nx_submet)); MR_dx_submet(:) = 0.0_sp
        allocate(MR_dy_submet(1:ny_submet)); MR_dy_submet(:) = 0.0_sp
        allocate(MR_dum3d_metP(nx_submet,ny_submet,np_fullmet))
        allocate(MR_dum3d_metH(nx_submet,ny_submet,nz_comp))
        allocate(MR_dum2d_comp_int(nx_comp,ny_comp))
        allocate(MR_dum2d_comp(nx_comp,ny_comp))
        allocate(MR_dum3d_compH(nx_comp,ny_comp,nz_comp))
        allocate(MR_dum3d_compP(nx_comp,ny_comp,np_fullmet))
        allocate(MR_geoH_metP_last(nx_submet,ny_submet,np_fullmet)) ! This will hold elevation
        allocate(MR_geoH_metP_next(nx_submet,ny_submet,np_fullmet)) ! Copy of geoH_metP_last
        allocate(MR_Snd2Comp_map_wgt(nx_comp,ny_comp,MR_nSnd_Locs))
        allocate(MR_Snd2Comp_map_idx(nx_comp,ny_comp,MR_nSnd_Locs)) ! contains the indeces of the
                                                                    ! MR_nstat stations used in weighting
        ! For the radiosonde case, compare computational grid with internal station list and 
        ! list available stations in the domain.
        if(MR_iwindformat.eq.2)then
          call MR_Get_Radiosonde_Stations_InDomain
        endif

        if(MR_nSnd_Locs.eq.1)then
          ! Easy case: all weights are on the one sonde point
          MR_Snd2Comp_map_wgt = 0.0_sp
          MR_Snd2Comp_map_wgt(1:nx_comp,1:ny_comp,1:MR_nSnd_Locs) = 1.0_sp
          MR_Snd2Comp_map_idx(1:nx_comp,1:ny_comp,1:MR_nSnd_Locs) = 1
        else
          ! Lots of algorithms to consider: 
          !  Triangulation would be a good one if we had heaps of stations; assume we don't
          !  Nearest Neighbor is easiest.
          !call MR_Calc_Snd_Weights_NearNeigh
          !  Inverse Distance is more flexible and can duplicate NN if MR_nstat=1
          call MR_Calc_Snd_Weights_InvDistWeight(MR_pexp,MR_nstat)
        endif

      else
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR : Unknown ASCII wind file format."
          write(errlog(io),*)"        MR_iwind = ",MR_iwind
          write(errlog(io),*)"  MR_iwindformat = ",MR_iwindformat
        endif;enddo
        stop 1
      endif

      do io=1,MR_nio;if(VB(io).le.verbosity_production)then
        write(outlog(io),*)&
           "-----------------------------------------------------------------------"
      endif;enddo

      end subroutine MR_Set_MetComp_Grids_ASCII_1d


!##############################################################################
!
!     MR_Set_Met_Times_ASCII
!
!     Called once from 
!
!     This subroutine 
!
!     Allocated the dummy arrays for storing met data on met and computational
!
!##############################################################################
!
!      subroutine MR_Set_Met_Times_ASCII
!
!      use MetReader
!      use projection
!
!      implicit none
!
!      integer, parameter :: sp        = 4 ! single precision
!      integer, parameter :: dp        = 8 ! double precision
!
!      If(MR_iwind.eq.1.and.MR_iwindformat.eq.1)then
!          ! 1d ascii file
!        allocate(MR_MetStep_File(2))
!        allocate(MR_MetStep_tindex(2))
!        allocate(MR_MetStep_Hour_since_baseyear(2))
!        allocate(MR_MetStep_Interval(1))
!        MR_MetSteps_Total    = 2
!        MR_ForecastInterval = 10000.0_dp       ! give it a long forecast interval
!        MR_MetStep_Interval(1) = MR_ForecastInterval
!        MR_MetStep_File(1)   = MR_windfiles(1) ! Currently, only one sounding
!        MR_MetStep_tindex(1) = 1            ! Only one time with sounding file
!        MR_MetStep_Hour_since_baseyear(1) = MR_Comp_StartHour - 1.0_sp
!
!        MR_MetStep_File(2)   = MR_windfiles(1) ! Currently, only one sounding
!        MR_MetStep_tindex(2) = 2            ! Only one time with sounding file
!        MR_MetStep_Hour_since_baseyear(2) = MR_Comp_StartHour + MR_Comp_Time_in_hours + 1.0_sp
!      elseif(MR_iwind.eq.1.and.MR_iwindformat.eq.2)then
!        ! Radiosonde case
!      elseif(MR_iwind.eq.2.and.MR_iwindformat.eq.1)then
!        ! 3d ascii windfiles
!      else
!        write(outlog(io),*)"Unknown ASCII wind file format."
!        write(outlog(io),*)"  MR_iwind = ",MR_iwind
!        write(outlog(io),*)"MR_iwindformat = ",MR_iwindformat
!        stop 1
!      endif
!
!      end subroutine MR_Set_Met_Times_ASCII


!##############################################################################
!
!     MR_Read_MetP_Variable_ASCII
!
!     Called from MR_Read_HGT_arrays and once from Read_3d_MetP_Variable.
!
!     Sets MR_dum3d_metP, MR_dum2d_met, or MR_dum2d_met_int as appropriate
!
!##############################################################################

      subroutine MR_Read_MetP_Variable_ASCII_1d(ivar,istep)

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_info,verbosity_error,&
         Met_var_IsAvailable,MR_Snd_nvars,MR_SndVarsID,MR_dum3d_metP,nx_submet,ny_submet,&
         MR_SndVars_metP,np_fullmet,MR_nSnd_Locs,MR_MetStep_findex

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer,intent(in) :: ivar
      integer,intent(in) :: istep

      integer :: icol, i, j
      integer :: itime

      integer :: io                           ! Index for output streams

      ! itime is not just the file index since we might have multiple stations per tstep
      itime = (MR_MetStep_findex(istep)-1)/MR_nSnd_Locs + 1
            ! these variables are set in MR_Read_Met_DimVars_ASCII_1d
      !  P,H,U,V,T
      if(Met_var_IsAvailable(ivar))then
        icol = 0
        do i=1,MR_Snd_nvars
          if(ivar.eq.MR_SndVarsID(i))then
            icol = i ! Found the column with the variable we want
            exit
          endif
        enddo
        if(icol.eq.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: Could not find variable"
          endif;enddo
          stop 1
        endif

        ! Now loop over all the submet points and build the value needed by multiplying by the
        ! weights determined in MR_Set_MetComp_Grids_ASCII_1d
        ! Recall that in contrast to the other data types, ASCII files are fully loaded into memory
        ! in var (MR_SndVars_metP) in the initial reading of the files above in
        ! MR_Read_Met_DimVars_ASCII_1d.
        MR_dum3d_metP(:,:,:) = 0.0_sp
        do i = 1,nx_submet
          do j = 1,ny_submet
            ! This just copies all the current met data to MR_dum3d_metP.
            ! Interpolation will be elsewhere
            ! First index of MR_SndVars_metP is iloc, but for ASCII, we have ny_submet=1
            ! and 1:nx_submet = 1:MR_nSnd_Locs
            MR_dum3d_metP(i,j,1:np_fullmet) = MR_SndVars_metP(i,itime,icol,1:np_fullmet)
          enddo
        enddo
      else
        ! W is typically not provided
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Attempting to read unavailable variable: ",ivar
        endif;enddo
        MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_fullmet) = 0.0_sp
      endif

      return

      end subroutine MR_Read_MetP_Variable_ASCII_1d

!##############################################################################
!
!     MR_Regrid_MetSonde2Comp
!
!     This subroutine interpolates the sonde data onto the computational grid
!     using the weights calculated
!
!##############################################################################

      subroutine MR_Regrid_MetSonde2Comp(nx1,ny1,wrk_met,nx2,ny2,wrk_comp)

      use MetReader,       only : &
         MR_Snd2Comp_map_wgt,MR_nSnd_Locs

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer                         ,intent(in)  :: nx1,ny1 ! for sondes, nx1=nloc, ny1=1
      real(kind=sp),dimension(nx1,ny1),intent(in)  :: wrk_met
      integer                         ,intent(in)  :: nx2,ny2
      real(kind=sp),dimension(nx2,ny2),intent(out) :: wrk_comp

      integer :: i,j,iloc

      wrk_comp = 0.0_sp
      do i = 1,nx2
        do j = 1,ny2
          do iloc = 1,MR_nSnd_Locs
            wrk_comp(i,j) = wrk_comp(i,j) + wrk_met(iloc,1) * MR_Snd2Comp_map_wgt(i,j,iloc)
          enddo
        enddo
      enddo

      end subroutine MR_Regrid_MetSonde2Comp

!##############################################################################
!
!     MR_Load_Radiosonde_Station_Data
!
!     This subroutine simply loads the radiosonde data
!
!##############################################################################

      subroutine MR_Load_Radiosonde_Station_Data()

      use MetReader,       only : &
         MR_Snd_cd,MR_Snd_id,MR_Snd_lt,MR_Snd_ln,MR_Snd_el,MR_Snd_lnm,&
         MR_Snd_st,MR_Snd_ct,num_RadSnd_Stat

      implicit none

      integer, parameter :: sp             = 4

      integer, parameter :: MAX_STAT_NUM   = 1200

      character(len=4 ),dimension(MAX_STAT_NUM) :: cd   ! Station code
      integer          ,dimension(MAX_STAT_NUM) :: id   ! WMO station ID
      real(kind=sp)    ,dimension(MAX_STAT_NUM) :: lt   ! Station latitude
      real(kind=sp)    ,dimension(MAX_STAT_NUM) :: ln   ! Station longitude
      real(kind=sp)    ,dimension(MAX_STAT_NUM) :: el   ! Station elevation
      character(len=25),dimension(MAX_STAT_NUM) :: lnm  ! Station long name
      character(len=2 ),dimension(MAX_STAT_NUM) :: st   ! Station state
      character(len=2 ),dimension(MAX_STAT_NUM) :: ct   ! Station country
      integer :: i

      ! Load the station data.  This is from two files available on 
      ! https://ruc.noaa.gov/raobs/General_Information.html
      ! stat2000.txt (for US sites) and international_sites.txt for the rest of the world
      ! These files are modified to only contain stations currently active
      ! Here we just load the local variables with short names (so we can fit each station on one line)
      ! and will copy to the MR global variables at the end of the subroutine.
      i = 0
      ! https://ruc.noaa.gov/raobs/stat2000.txt
i=i+1;cd(i)=" NPC";id(i)=78730;lt(i)= 14.05;ln(i)=  -83.57;el(i)=  20;lnm(i)="PUERTO CABEZAS(MIL)      ";st(i)="99";ct(i)="NI"
i=i+1;cd(i)=" EPZ";id(i)=72364;lt(i)= 31.90;ln(i)= -106.70;el(i)=1257;lnm(i)="SANTA TERESA             ";st(i)="NM";ct(i)="US"
i=i+1;cd(i)=" VEF";id(i)=72388;lt(i)= 36.05;ln(i)= -115.18;el(i)= 693;lnm(i)="LAS VEGAS                ";st(i)="NV";ct(i)="US"
i=i+1;cd(i)=" YUM";id(i)=72280;lt(i)= 32.87;ln(i)= -114.33;el(i)= 131;lnm(i)="YUMA/US ARMY MET TEAM    ";st(i)="AZ";ct(i)="US"
i=i+1;cd(i)=" DRA";id(i)=72387;lt(i)= 36.62;ln(i)= -116.02;el(i)=1007;lnm(i)="DESERT ROCK/MERCURY      ";st(i)="NV";ct(i)="US"
i=i+1;cd(i)=" NKX";id(i)=72293;lt(i)= 32.87;ln(i)= -117.15;el(i)= 134;lnm(i)="MIRAMAR NAS              ";st(i)="CA";ct(i)="US"
i=i+1;cd(i)=" EDW";id(i)=72381;lt(i)= 34.90;ln(i)= -117.92;el(i)= 724;lnm(i)="EDWARDS/AFB - UPPER AIR  ";st(i)="CA";ct(i)="US"
i=i+1;cd(i)=" REV";id(i)=72489;lt(i)= 39.57;ln(i)= -119.80;el(i)=1516;lnm(i)="RENO                     ";st(i)="NV";ct(i)="US"
i=i+1;cd(i)=" LCH";id(i)=72240;lt(i)= 30.12;ln(i)=  -93.22;el(i)=   5;lnm(i)="LAKE CHARLES             ";st(i)="LA";ct(i)="US"
i=i+1;cd(i)=" JAN";id(i)=72235;lt(i)= 32.32;ln(i)=  -90.07;el(i)=  91;lnm(i)="JACKSON/THOMPSON FLD     ";st(i)="MS";ct(i)="US"
i=i+1;cd(i)=" OUN";id(i)=72357;lt(i)= 35.23;ln(i)=  -97.47;el(i)= 362;lnm(i)="NORMAN                   ";st(i)="OK";ct(i)="US"
i=i+1;cd(i)=" LZK";id(i)=72340;lt(i)= 34.83;ln(i)=  -92.27;el(i)= 172;lnm(i)="N LITTLE ROCK            ";st(i)="AR";ct(i)="US"
i=i+1;cd(i)=" FWD";id(i)=72249;lt(i)= 32.80;ln(i)=  -97.30;el(i)= 196;lnm(i)="FT WORTH                 ";st(i)="TX";ct(i)="US"
i=i+1;cd(i)=" TFX";id(i)=72776;lt(i)= 47.45;ln(i)= -111.38;el(i)=1130;lnm(i)="GREAT FALLS              ";st(i)="MT";ct(i)="US"
i=i+1;cd(i)=" LKN";id(i)=72582;lt(i)= 40.87;ln(i)= -115.73;el(i)=1608;lnm(i)="ELKO                     ";st(i)="NV";ct(i)="US"
i=i+1;cd(i)=" OTX";id(i)=72786;lt(i)= 47.68;ln(i)= -117.63;el(i)= 728;lnm(i)="SPOKANE INTNL APT        ";st(i)="WA";ct(i)="US"
i=i+1;cd(i)=" YMW";id(i)=71722;lt(i)= 46.38;ln(i)=  -75.97;el(i)= 170;lnm(i)="MANIWAKI                 ";st(i)="PQ";ct(i)="CA"
i=i+1;cd(i)=" DTX";id(i)=72632;lt(i)= 42.70;ln(i)=  -83.47;el(i)= 329;lnm(i)="DETROIT/PONTIAC          ";st(i)="MI";ct(i)="US"
i=i+1;cd(i)=" ILX";id(i)=74560;lt(i)= 40.15;ln(i)=  -89.33;el(i)= 178;lnm(i)="LINCOLN-LOGAN COUNTY AP  ";st(i)="IL";ct(i)="US"
i=i+1;cd(i)=" APX";id(i)=72634;lt(i)= 44.91;ln(i)=  -84.72;el(i)= 448;lnm(i)="GAYLORD / ALPENA         ";st(i)="MI";ct(i)="US"
i=i+1;cd(i)=" ROL";id(i)=78762;lt(i)=  9.98;ln(i)=  -84.22;el(i)= 920;lnm(i)="SAN JOSE/JUAN SANTA MARIA";st(i)="99";ct(i)="CR"
i=i+1;cd(i)=" BDI";id(i)=78954;lt(i)= 13.07;ln(i)=  -59.50;el(i)=  47;lnm(i)="SEAWELL APT              ";st(i)="99";ct(i)="BB"
i=i+1;cd(i)=" SDQ";id(i)=78486;lt(i)= 18.47;ln(i)=  -69.88;el(i)=  14;lnm(i)="SANTO DOMINGO            ";st(i)="99";ct(i)="DO"
i=i+1;cd(i)=" KPP";id(i)=78970;lt(i)= 10.58;ln(i)=  -61.35;el(i)=  12;lnm(i)="TRINIDAD/PIARCO IAP      ";st(i)="99";ct(i)="TT"
i=i+1;cd(i)=" JSJ";id(i)=78526;lt(i)= 18.43;ln(i)=  -66.00;el(i)=   3;lnm(i)="SAN JUAN/ISLA VERDE      ";st(i)="PR";ct(i)="US"
i=i+1;cd(i)=" FFR";id(i)=78897;lt(i)= 16.27;ln(i)=  -61.52;el(i)=   8;lnm(i)="POINT A PITRE/RAIZET     ";st(i)="99";ct(i)="GP"
i=i+1;cd(i)=" ACC";id(i)=78988;lt(i)= 12.20;ln(i)=  -68.97;el(i)=  54;lnm(i)="CURACAO/WILLEMSTAD       ";st(i)="99";ct(i)="AN"
i=i+1;cd(i)=" ACM";id(i)=78866;lt(i)= 18.05;ln(i)=  -63.12;el(i)=   3;lnm(i)="SINT MARTIN/JULIANA      ";st(i)="99";ct(i)="AN"
i=i+1;cd(i)=" KJP";id(i)=78397;lt(i)= 17.93;ln(i)=  -76.78;el(i)=   1;lnm(i)="KINGSTON/PALISADOES      ";st(i)="99";ct(i)="JM"
i=i+1;cd(i)=" KCR";id(i)=78384;lt(i)= 19.30;ln(i)=  -81.37;el(i)=   3;lnm(i)="GRAND CAYMAN             ";st(i)="99";ct(i)="KY"
i=i+1;cd(i)=" ZBZ";id(i)=78583;lt(i)= 17.53;ln(i)=  -88.30;el(i)=   5;lnm(i)="BELIZE                   ";st(i)="99";ct(i)="BZ"
i=i+1;cd(i)=" MEX";id(i)=76679;lt(i)= 19.43;ln(i)=  -99.07;el(i)=2309;lnm(i)="MEXICO CITY/INT APT      ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" VER";id(i)=76692;lt(i)= 19.17;ln(i)=  -96.12;el(i)=  13;lnm(i)="VERACRUZ                 ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" 999";id(i)=76805;lt(i)= 16.76;ln(i)=  -99.93;el(i)=   3;lnm(i)="ACCAPULCO                ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" YNN";id(i)=78073;lt(i)= 25.05;ln(i)=  -77.47;el(i)=   2;lnm(i)="NASSAU APT               ";st(i)="99";ct(i)="BS"
i=i+1;cd(i)=" EYW";id(i)=72201;lt(i)= 24.55;ln(i)=  -81.75;el(i)=   1;lnm(i)="KEY WEST INT AP          ";st(i)="FL";ct(i)="US"
i=i+1;cd(i)=" TBW";id(i)=72210;lt(i)= 27.70;ln(i)=  -82.40;el(i)=  13;lnm(i)="TAMPA BAY/RUSKIN         ";st(i)="FL";ct(i)="US"
i=i+1;cd(i)=" XMR";id(i)=74794;lt(i)= 28.48;ln(i)=  -80.55;el(i)=   5;lnm(i)="CAPE KENNEDY             ";st(i)="FL";ct(i)="US"
i=i+1;cd(i)=" MID";id(i)=76644;lt(i)= 20.95;ln(i)=  -89.65;el(i)=  11;lnm(i)="MERIDA IAP               ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" BRO";id(i)=72250;lt(i)= 25.90;ln(i)=  -97.43;el(i)=   7;lnm(i)="BROWNSVILLE              ";st(i)="TX";ct(i)="US"
i=i+1;cd(i)=" CRP";id(i)=72251;lt(i)= 27.77;ln(i)=  -97.50;el(i)=  14;lnm(i)="CORPUS CHRISTI           ";st(i)="TX";ct(i)="US"
i=i+1;cd(i)=" XKF";id(i)=78016;lt(i)= 32.37;ln(i)=  -64.68;el(i)=  37;lnm(i)="BERMUDA/(MCKINDLY AFB)   ";st(i)="99";ct(i)="BS"
i=i+1;cd(i)=" APG";id(i)=74002;lt(i)= 39.47;ln(i)=  -76.07;el(i)=   5;lnm(i)="PHILLIPS AFB, ABERDEEN   ";st(i)="MD";ct(i)="US"
i=i+1;cd(i)=" GSO";id(i)=72317;lt(i)= 36.08;ln(i)=  -79.95;el(i)= 277;lnm(i)="GREENSBORO               ";st(i)="NC";ct(i)="US"
i=i+1;cd(i)=" ILN";id(i)=72426;lt(i)= 39.42;ln(i)=  -83.82;el(i)= 317;lnm(i)="WILMINGTON               ";st(i)="OH";ct(i)="US"
i=i+1;cd(i)=" VPS";id(i)=72221;lt(i)= 30.52;ln(i)=  -86.58;el(i)=  20;lnm(i)="VALPARAISO/ELGIN AFB     ";st(i)="FL";ct(i)="US"
i=i+1;cd(i)=" CHS";id(i)=72208;lt(i)= 32.90;ln(i)=  -80.03;el(i)=  15;lnm(i)="CHARLESTON               ";st(i)="SC";ct(i)="US"
i=i+1;cd(i)=" JAX";id(i)=72206;lt(i)= 30.43;ln(i)=  -81.70;el(i)=  10;lnm(i)="JACKSONVILLE             ";st(i)="FL";ct(i)="US"
i=i+1;cd(i)=" BNA";id(i)=72327;lt(i)= 36.25;ln(i)=  -86.57;el(i)= 180;lnm(i)="NASHVILLE                ";st(i)="TN";ct(i)="US"
i=i+1;cd(i)=" FSI";id(i)=72355;lt(i)= 34.65;ln(i)=  -98.40;el(i)= 362;lnm(i)="FORT SILL                ";st(i)="OK";ct(i)="US"
i=i+1;cd(i)=" SHV";id(i)=72248;lt(i)= 32.45;ln(i)=  -93.83;el(i)=  84;lnm(i)="SHREVEPORT REGIONAL AP   ";st(i)="LA";ct(i)="US"
i=i+1;cd(i)=" DDC";id(i)=72451;lt(i)= 37.77;ln(i)=  -99.97;el(i)= 791;lnm(i)="DODGE CITY               ";st(i)="KS";ct(i)="US"
i=i+1;cd(i)=" SGF";id(i)=72440;lt(i)= 37.23;ln(i)=  -93.40;el(i)= 394;lnm(i)="SPRINGFIELD REGIONAL AP  ";st(i)="MO";ct(i)="US"
i=i+1;cd(i)=" TOP";id(i)=72456;lt(i)= 39.07;ln(i)=  -95.62;el(i)= 268;lnm(i)="TOPEKA                   ";st(i)="KS";ct(i)="US"
i=i+1;cd(i)=" YJT";id(i)=71815;lt(i)= 48.53;ln(i)=  -58.55;el(i)=  60;lnm(i)="STEPHENVILLE/HARMON AFB  ";st(i)="NF";ct(i)="CA"
i=i+1;cd(i)=" YYT";id(i)=71801;lt(i)= 47.67;ln(i)=  -52.75;el(i)= 140;lnm(i)="TORBAY/ST JOHNS          ";st(i)="NF";ct(i)="CA"
i=i+1;cd(i)=" CAR";id(i)=72712;lt(i)= 46.87;ln(i)=  -68.02;el(i)= 191;lnm(i)="CARIBOU                  ";st(i)="ME";ct(i)="US"
i=i+1;cd(i)=" YSA";id(i)=71600;lt(i)= 43.93;ln(i)=  -60.02;el(i)=   4;lnm(i)="SABLE ISLAND             ";st(i)="NS";ct(i)="CA"
i=i+1;cd(i)=" CHH";id(i)=74494;lt(i)= 41.67;ln(i)=  -69.97;el(i)=  16;lnm(i)="CHATHAM                  ";st(i)="MA";ct(i)="US"
i=i+1;cd(i)=" YCX";id(i)=71701;lt(i)= 45.83;ln(i)=  -66.43;el(i)=  52;lnm(i)="GAGETOWN                 ";st(i)="NB";ct(i)="CA"
i=i+1;cd(i)=" BUF";id(i)=72528;lt(i)= 42.93;ln(i)=  -78.73;el(i)= 218;lnm(i)="BUFFALO/GRTR ARPT        ";st(i)="NY";ct(i)="US"
i=i+1;cd(i)=" GRB";id(i)=72645;lt(i)= 44.48;ln(i)=  -88.13;el(i)= 210;lnm(i)="GREEN BAY                ";st(i)="WI";ct(i)="US"
i=i+1;cd(i)=" INL";id(i)=72747;lt(i)= 48.57;ln(i)=  -93.38;el(i)= 359;lnm(i)="INTERNATIONAL FALLS      ";st(i)="MN";ct(i)="US"
i=i+1;cd(i)=" ABR";id(i)=72659;lt(i)= 45.45;ln(i)=  -98.42;el(i)= 397;lnm(i)="ABERDEEN                 ";st(i)="SD";ct(i)="US"
i=i+1;cd(i)=" YYR";id(i)=71816;lt(i)= 53.30;ln(i)=  -60.37;el(i)=  36;lnm(i)="GOOSE/GOOSE BAY          ";st(i)="NF";ct(i)="CA"
i=i+1;cd(i)=" YZV";id(i)=71811;lt(i)= 50.22;ln(i)=  -66.27;el(i)=  52;lnm(i)="SEPT ILES (UA)           ";st(i)="PQ";ct(i)="CA"
i=i+1;cd(i)=" YVP";id(i)=71906;lt(i)= 58.10;ln(i)=  -68.42;el(i)=  60;lnm(i)="KUUJJUAQ (UA)            ";st(i)="PQ";ct(i)="CA"
i=i+1;cd(i)=" YPH";id(i)=71907;lt(i)= 58.45;ln(i)=  -78.12;el(i)=  26;lnm(i)="INUKJUAK                 ";st(i)="PQ";ct(i)="CA"
i=i+1;cd(i)=" YAH";id(i)=71823;lt(i)= 53.75;ln(i)=  -73.67;el(i)= 307;lnm(i)="LA GRANDE IV             ";st(i)="PQ";ct(i)="CA"
i=i+1;cd(i)=" YMO";id(i)=71836;lt(i)= 51.27;ln(i)=  -80.65;el(i)=  10;lnm(i)="MOOSONEE                 ";st(i)="PQ";ct(i)="CA"
i=i+1;cd(i)=" YYQ";id(i)=71913;lt(i)= 58.75;ln(i)=  -94.07;el(i)=  29;lnm(i)="CHURCHILL                ";st(i)="MB";ct(i)="CA"
i=i+1;cd(i)=" WPL";id(i)=71845;lt(i)= 51.47;ln(i)=  -90.20;el(i)= 386;lnm(i)="PICKLE LAKE              ";st(i)="ON";ct(i)="CA"
i=i+1;cd(i)=" YVN";id(i)=71909;lt(i)= 63.75;ln(i)=  -68.55;el(i)=  35;lnm(i)="IQALUIT (UA)             ";st(i)="NW";ct(i)="CA"
i=i+1;cd(i)=" YZS";id(i)=71915;lt(i)= 64.20;ln(i)=  -83.37;el(i)=  57;lnm(i)="CORAL HARBOUR            ";st(i)="NW";ct(i)="CA"
i=i+1;cd(i)=" YUX";id(i)=71081;lt(i)= 68.78;ln(i)=  -81.25;el(i)=   7;lnm(i)="HALL BEACH/HALL LK       ";st(i)="NW";ct(i)="CA"
i=i+1;cd(i)=" YBK";id(i)=71926;lt(i)= 64.30;ln(i)=  -96.00;el(i)=  49;lnm(i)="BAKER LAKE (UA)          ";st(i)="NW";ct(i)="CA"
i=i+1;cd(i)=" YRB";id(i)=71924;lt(i)= 74.72;ln(i)=  -94.98;el(i)=  40;lnm(i)="RESOLUTE                 ";st(i)="NW";ct(i)="CA"
i=i+1;cd(i)=" YLT";id(i)=71082;lt(i)= 82.50;ln(i)=  -62.33;el(i)=  66;lnm(i)="ALERT                    ";st(i)="NW";ct(i)="CA"
i=i+1;cd(i)=" YEU";id(i)=71917;lt(i)= 79.98;ln(i)=  -85.93;el(i)=  10;lnm(i)="EUREKA                   ";st(i)="NW";ct(i)="CA"
i=i+1;cd(i)=" 999";id(i)=76654;lt(i)= 19.07;ln(i)= -104.33;el(i)=   3;lnm(i)="MANZANILLO               ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" ITO";id(i)=91285;lt(i)= 19.72;ln(i)= -155.07;el(i)=  10;lnm(i)="HILO                     ";st(i)="HI";ct(i)="US"
i=i+1;cd(i)=" MCV";id(i)=76225;lt(i)= 28.70;ln(i)= -106.07;el(i)=1428;lnm(i)="CHIHUAHUA                ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" MZT";id(i)=76458;lt(i)= 23.18;ln(i)= -106.42;el(i)=   4;lnm(i)="MAZATLAN SINALOA         ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" DRT";id(i)=72261;lt(i)= 29.37;ln(i)= -100.92;el(i)= 313;lnm(i)="DEL RIO                  ";st(i)="TX";ct(i)="US"
i=i+1;cd(i)=" MTY";id(i)=76394;lt(i)= 25.87;ln(i)= -100.20;el(i)= 450;lnm(i)="MONTERREY                ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" 999";id(i)=76612;lt(i)= 20.68;ln(i)= -103.33;el(i)=1551;lnm(i)="GUADALAJARA              ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" LAP";id(i)=76405;lt(i)= 24.07;ln(i)= -110.33;el(i)=  14;lnm(i)="LA PAZ/DE LEON           ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" GYM";id(i)=76256;lt(i)= 27.95;ln(i)= -110.80;el(i)=  12;lnm(i)="EMPALME SONORA           ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" YEV";id(i)=71957;lt(i)= 68.32;ln(i)= -133.53;el(i)= 103;lnm(i)="INUVIK (UA)              ";st(i)="NW";ct(i)="CA"
i=i+1;cd(i)=" LIH";id(i)=91165;lt(i)= 21.98;ln(i)= -159.35;el(i)=  36;lnm(i)="LIHUE/KAUAI              ";st(i)="HI";ct(i)="US"
i=i+1;cd(i)=" MAF";id(i)=72265;lt(i)= 31.93;ln(i)= -102.20;el(i)= 873;lnm(i)="MIDLAND                  ";st(i)="TX";ct(i)="US"
i=i+1;cd(i)=" AMA";id(i)=72363;lt(i)= 35.23;ln(i)= -101.70;el(i)=1095;lnm(i)="AMARILLO                 ";st(i)="TX";ct(i)="US"
i=i+1;cd(i)=" ABQ";id(i)=72365;lt(i)= 35.05;ln(i)= -106.62;el(i)=1619;lnm(i)="ALBUQUERQUE              ";st(i)="NM";ct(i)="US"
i=i+1;cd(i)=" DNR";id(i)=72469;lt(i)= 39.77;ln(i)= -104.88;el(i)=1611;lnm(i)="DENVER/STAPLETON ARPT    ";st(i)="CO";ct(i)="US"
i=i+1;cd(i)=" GJT";id(i)=72476;lt(i)= 39.12;ln(i)= -108.53;el(i)=1472;lnm(i)="GRAND JUNCTION           ";st(i)="CO";ct(i)="US"
i=i+1;cd(i)=" TWC";id(i)=72274;lt(i)= 32.23;ln(i)= -110.96;el(i)= 753;lnm(i)="TUCSON                   ";st(i)="AZ";ct(i)="US"
i=i+1;cd(i)=" OAK";id(i)=72493;lt(i)= 37.75;ln(i)= -122.22;el(i)=   6;lnm(i)="OAKLAND  INT AP          ";st(i)="CA";ct(i)="US"
i=i+1;cd(i)=" BIS";id(i)=72764;lt(i)= 46.77;ln(i)= -100.75;el(i)= 503;lnm(i)="BISMARCK                 ";st(i)="ND";ct(i)="US"
i=i+1;cd(i)=" LBF";id(i)=72562;lt(i)= 41.13;ln(i)= -100.68;el(i)= 847;lnm(i)="NORTH PLATTE             ";st(i)="NE";ct(i)="US"
i=i+1;cd(i)=" RIW";id(i)=72672;lt(i)= 43.06;ln(i)= -108.47;el(i)=1688;lnm(i)="RIVERTON                 ";st(i)="WY";ct(i)="US"
i=i+1;cd(i)=" SLC";id(i)=72572;lt(i)= 40.77;ln(i)= -111.97;el(i)=1288;lnm(i)="SALT LAKE CITY           ";st(i)="UT";ct(i)="US"
i=i+1;cd(i)=" BOI";id(i)=72681;lt(i)= 43.57;ln(i)= -116.22;el(i)= 871;lnm(i)="BOISE                    ";st(i)="ID";ct(i)="US"
i=i+1;cd(i)=" MFR";id(i)=72597;lt(i)= 42.37;ln(i)= -122.87;el(i)= 397;lnm(i)="MEDFORD                  ";st(i)="OR";ct(i)="US"
i=i+1;cd(i)=" SLE";id(i)=72694;lt(i)= 44.92;ln(i)= -123.02;el(i)=  61;lnm(i)="SALEM                    ";st(i)="OR";ct(i)="US"
i=i+1;cd(i)=" YQD";id(i)=71867;lt(i)= 53.97;ln(i)= -101.10;el(i)= 273;lnm(i)="THE PAS                  ";st(i)="MB";ct(i)="CA"
i=i+1;cd(i)=" WSE";id(i)=71119;lt(i)= 53.55;ln(i)= -114.10;el(i)= 766;lnm(i)="EDMONTON/STONY PLAIN     ";st(i)="AB";ct(i)="CA"
i=i+1;cd(i)=" YBP";id(i)=71126;lt(i)= 50.63;ln(i)= -111.90;el(i)= 759;lnm(i)="BROOKS                   ";st(i)="AB";ct(i)="CA"
i=i+1;cd(i)=" WIQ";id(i)=71124;lt(i)= 54.80;ln(i)= -110.08;el(i)= 703;lnm(i)="PRIMROSE LAKE            ";st(i)="AB";ct(i)="CA"
i=i+1;cd(i)=" ZXS";id(i)=71908;lt(i)= 53.90;ln(i)= -122.80;el(i)= 601;lnm(i)="PRINCE GEORGE            ";st(i)="BC";ct(i)="CA"
i=i+1;cd(i)=" YZT";id(i)=71109;lt(i)= 50.68;ln(i)= -127.37;el(i)=  17;lnm(i)="PORT HARDY               ";st(i)="BC";ct(i)="CA"
i=i+1;cd(i)=" YYE";id(i)=71945;lt(i)= 58.83;ln(i)= -122.60;el(i)= 377;lnm(i)="FORT NELSON UA           ";st(i)="BC";ct(i)="CA"
i=i+1;cd(i)=" ANN";id(i)=70398;lt(i)= 55.03;ln(i)= -131.57;el(i)=  37;lnm(i)="ANNETTE ISLAND           ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" YAK";id(i)=70361;lt(i)= 59.52;ln(i)= -139.67;el(i)=  10;lnm(i)="YAKUTAT                  ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" ADQ";id(i)=70350;lt(i)= 57.75;ln(i)= -152.48;el(i)=   4;lnm(i)="KODIAK                   ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" AKN";id(i)=70326;lt(i)= 58.68;ln(i)= -156.65;el(i)=  15;lnm(i)="KING SALMON              ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" CDB";id(i)=70316;lt(i)= 55.20;ln(i)= -162.72;el(i)=  30;lnm(i)="COLD BAY                 ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" SNP";id(i)=70308;lt(i)= 57.15;ln(i)= -170.22;el(i)=  10;lnm(i)="ST PAUL ISLAND           ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" YCB";id(i)=71925;lt(i)= 69.10;ln(i)= -105.12;el(i)=  25;lnm(i)="CAMBRIDGE BAY            ";st(i)="NW";ct(i)="CA"
i=i+1;cd(i)=" YSM";id(i)=71934;lt(i)= 60.03;ln(i)= -111.95;el(i)= 203;lnm(i)="FT SMITH (UA)            ";st(i)="NW";ct(i)="CA"
i=i+1;cd(i)=" YVQ";id(i)=71043;lt(i)= 65.28;ln(i)= -126.75;el(i)=  95;lnm(i)="NORMAN WELLS (UA)        ";st(i)="NW";ct(i)="CA"
i=i+1;cd(i)=" YXY";id(i)=71964;lt(i)= 60.72;ln(i)= -135.07;el(i)= 704;lnm(i)="WHITEHORSE               ";st(i)="YK";ct(i)="CA"
i=i+1;cd(i)=" ANC";id(i)=70273;lt(i)= 61.17;ln(i)= -150.02;el(i)=  45;lnm(i)="ANCHORAGE IAP/PT. CAMPBE ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" FAI";id(i)=70261;lt(i)= 64.82;ln(i)= -147.87;el(i)= 135;lnm(i)="FAIRBANKS                ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" MCG";id(i)=70231;lt(i)= 62.97;ln(i)= -155.62;el(i)= 103;lnm(i)="MCGRATH                  ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" BET";id(i)=70219;lt(i)= 60.78;ln(i)= -161.80;el(i)=  36;lnm(i)="BETHEL                   ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" OTZ";id(i)=70133;lt(i)= 66.87;ln(i)= -162.63;el(i)=   5;lnm(i)="KOTZEBUE                 ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" OME";id(i)=70200;lt(i)= 64.50;ln(i)= -165.43;el(i)=   5;lnm(i)="NOME AP                  ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" BRW";id(i)=70026;lt(i)= 71.30;ln(i)= -156.78;el(i)=  12;lnm(i)="POINT BARROW             ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" BRW";id(i)=70027;lt(i)= 71.32;ln(i)= -156.65;el(i)=   3;lnm(i)="BARROW/POINT BARROW      ";st(i)="AK";ct(i)="US"
i=i+1;cd(i)=" FGZ";id(i)=72376;lt(i)= 35.23;ln(i)= -111.82;el(i)=2179;lnm(i)="FLAGSTAFF/BELLEMT (ARMY) ";st(i)="AZ";ct(i)="US"
i=i+1;cd(i)=" SIL";id(i)=72233;lt(i)= 30.33;ln(i)=  -89.82;el(i)=   8;lnm(i)="SLIDELL                  ";st(i)="LA";ct(i)="US"
i=i+1;cd(i)=" FFC";id(i)=72215;lt(i)= 33.35;ln(i)=  -84.56;el(i)= 246;lnm(i)="PEACHTREE CITY           ";st(i)="GA";ct(i)="US"
i=i+1;cd(i)=" BMX";id(i)=72230;lt(i)= 33.10;ln(i)=  -86.70;el(i)= 178;lnm(i)="BIRMINGHAM (SHELBY APT)  ";st(i)="AL";ct(i)="US"
i=i+1;cd(i)=" RNK";id(i)=72318;lt(i)= 37.20;ln(i)=  -80.41;el(i)= 648;lnm(i)="ROANOKE/BLACKSBURG       ";st(i)="VA";ct(i)="US"
i=i+1;cd(i)=" GYX";id(i)=74389;lt(i)= 43.89;ln(i)=  -70.25;el(i)= 125;lnm(i)="GRAY                     ";st(i)="ME";ct(i)="US"
i=i+1;cd(i)=" ALY";id(i)=72518;lt(i)= 42.69;ln(i)=  -73.83;el(i)=  94;lnm(i)="ALBANY                   ";st(i)="NY";ct(i)="US"
i=i+1;cd(i)=" MFL";id(i)=72202;lt(i)= 25.75;ln(i)=  -80.38;el(i)=   4;lnm(i)="MIAMI/FL INTL UNIV       ";st(i)="FL";ct(i)="US"
i=i+1;cd(i)=" NTD";id(i)=72391;lt(i)= 34.10;ln(i)= -119.12;el(i)=   2;lnm(i)="POINT MUGU               ";st(i)="CA";ct(i)="US"
i=i+1;cd(i)=" NSI";id(i)=72291;lt(i)= 33.25;ln(i)= -119.45;el(i)=  14;lnm(i)="SAN NICOLAS ISLAND/SITE1 ";st(i)="CA";ct(i)="US"
i=i+1;cd(i)=" VBG";id(i)=72393;lt(i)= 34.75;ln(i)= -120.57;el(i)= 100;lnm(i)="VANDENBERG               ";st(i)="CA";ct(i)="US"
i=i+1;cd(i)=" IAD";id(i)=72403;lt(i)= 38.98;ln(i)=  -77.47;el(i)=  85;lnm(i)="STERLING(WASH DULLES)    ";st(i)="VA";ct(i)="US"
i=i+1;cd(i)=" WAL";id(i)=72402;lt(i)= 37.93;ln(i)=  -75.48;el(i)=  13;lnm(i)="WALLOPS ISLAND           ";st(i)="VA";ct(i)="US"
i=i+1;cd(i)=" MHX";id(i)=72305;lt(i)= 34.70;ln(i)=  -76.80;el(i)=  11;lnm(i)="MOREHEAD CITY/NEWPORT    ";st(i)="NC";ct(i)="US"
i=i+1;cd(i)=" TLH";id(i)=72214;lt(i)= 30.45;ln(i)=  -84.30;el(i)=  52;lnm(i)="TALLAHASEE               ";st(i)="FL";ct(i)="US"
i=i+1;cd(i)=" GGW";id(i)=72768;lt(i)= 48.20;ln(i)= -106.62;el(i)= 693;lnm(i)="GLASGOW                  ";st(i)="MT";ct(i)="US"
i=i+1;cd(i)=" UNR";id(i)=72662;lt(i)= 44.07;ln(i)= -103.21;el(i)=1037;lnm(i)="RAPID CITY               ";st(i)="SD";ct(i)="US"
i=i+1;cd(i)=" YLW";id(i)=71203;lt(i)= 49.97;ln(i)= -119.38;el(i)= 454;lnm(i)="KELOWNA APT              ";st(i)="BC";ct(i)="CA"
i=i+1;cd(i)=" UIL";id(i)=72797;lt(i)= 47.95;ln(i)= -124.55;el(i)=  56;lnm(i)="QUILLAYUTE               ";st(i)="WA";ct(i)="US"
i=i+1;cd(i)=" WQI";id(i)=71603;lt(i)= 43.87;ln(i)=  -66.05;el(i)=   9;lnm(i)="YARMOUTH                 ";st(i)="NS";ct(i)="CA"
i=i+1;cd(i)=" OKX";id(i)=72501;lt(i)= 40.87;ln(i)=  -72.87;el(i)=  20;lnm(i)="BROOKHAVEN               ";st(i)="NY";ct(i)="US"
i=i+1;cd(i)=" PIT";id(i)=72520;lt(i)= 40.53;ln(i)=  -80.23;el(i)= 360;lnm(i)="PITTSBURGH/MOON TOWNSHIP ";st(i)="PA";ct(i)="US"
i=i+1;cd(i)=" OAX";id(i)=72558;lt(i)= 41.32;ln(i)=  -96.37;el(i)= 350;lnm(i)="OMAHA/VALLEY             ";st(i)="NE";ct(i)="US"
i=i+1;cd(i)=" DVN";id(i)=74455;lt(i)= 41.60;ln(i)=  -90.57;el(i)= 229;lnm(i)="DAVENPORT MUNICIPAL AP   ";st(i)="IA";ct(i)="US"
i=i+1;cd(i)=" MPX";id(i)=72649;lt(i)= 44.83;ln(i)=  -93.55;el(i)= 287;lnm(i)="MINNEAPOLIS              ";st(i)="MN";ct(i)="US"
i=i+1;cd(i)=" WEP";id(i)=71230;lt(i)= 55.13;ln(i)= -122.96;el(i)= 659;lnm(i)="WHISTLER OLYMPIC SITE    ";st(i)="BC";ct(i)="CA"
i=i+1;cd(i)=" XBK";id(i)=71569;lt(i)= 50.20;ln(i)= -104.70;el(i)= 580;lnm(i)="BRATTS LAKE              ";st(i)="SA";ct(i)="CA"
i=i+1;cd(i)=" 999";id(i)=71802;lt(i)= 47.52;ln(i)=  -52.78;el(i)=  49;lnm(i)="SAINT LAWRENCE           ";st(i)="NF";ct(i)="CA"
i=i+1;cd(i)=" CXW";id(i)=71843;lt(i)= 49.88;ln(i)=  -97.13;el(i)= 263;lnm(i)="WINNIPEG                 ";st(i)="99";ct(i)="CA"
i=i+1;cd(i)=" YQW";id(i)=71876;lt(i)= 52.77;ln(i)= -108.25;el(i)= 548;lnm(i)="NORTH BATTLEFORD         ";st(i)="99";ct(i)="CA"
i=i+1;cd(i)=" YQQ";id(i)=71893;lt(i)= 49.72;ln(i)= -124.90;el(i)=  24;lnm(i)="COMOX                    ";st(i)="BC";ct(i)="CA"
i=i+1;cd(i)=" 999";id(i)=74001;lt(i)= 34.60;ln(i)=  -86.62;el(i)= 174;lnm(i)="REDSTONE ARSENAL         ";st(i)="AL";ct(i)="US"
i=i+1;cd(i)=" 999";id(i)=74005;lt(i)= 33.53;ln(i)= -114.03;el(i)= 331;lnm(i)="LUKE AIR FORCE BASE      ";st(i)="AZ";ct(i)="US"
i=i+1;cd(i)=" 999";id(i)=74006;lt(i)= 32.86;ln(i)= -114.03;el(i)= 331;lnm(i)="YUMA PROVING GROUND      ";st(i)="AZ";ct(i)="US"
i=i+1;cd(i)=" 999";id(i)=74626;lt(i)= 33.45;ln(i)= -111.95;el(i)= 384;lnm(i)="PHOENIX                  ";st(i)="AZ";ct(i)="US"
i=i+1;cd(i)=" LMN";id(i)=74646;lt(i)= 36.68;ln(i)=  -97.47;el(i)= 306;lnm(i)="LAMONT                   ";st(i)="OK";ct(i)="US"
i=i+1;cd(i)=" 999";id(i)=76526;lt(i)= 22.75;ln(i)= -102.51;el(i)=2265;lnm(i)="GUADALUPE, ZACHATECAS    ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" CUN";id(i)=76595;lt(i)= 21.03;ln(i)=  -86.85;el(i)=   5;lnm(i)="CANCUN                   ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" 999";id(i)=76743;lt(i)= 17.98;ln(i)=  -92.92;el(i)=   6;lnm(i)="VILLAHERMOSA             ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" 999";id(i)=76856;lt(i)= 15.65;ln(i)=  -96.50;el(i)=  43;lnm(i)="PUERTO ANGEL             ";st(i)="99";ct(i)="MX"
i=i+1;cd(i)=" 999";id(i)=78083;lt(i)= 24.70;ln(i)=  -77.77;el(i)=   3;lnm(i)="ANDROS TOWN              ";st(i)="99";ct(i)="BA"
i=i+1;cd(i)=" 999";id(i)=78387;lt(i)= 18.25;ln(i)=  -78.37;el(i)=  10;lnm(i)="NEGRIL POINT (LH)        ";st(i)="99";ct(i)="JM"
i=i+1;cd(i)=" 999";id(i)=78485;lt(i)= 18.42;ln(i)=  -69.66;el(i)=  18;lnm(i)="CAUCEDO/DE LAS           ";st(i)="99";ct(i)="DO"
i=i+1;cd(i)=" 999";id(i)=78588;lt(i)= 17.18;ln(i)=  -87.50;el(i)=   1;lnm(i)="HALF MOON CAYE           ";st(i)="99";ct(i)="BZ"
i=i+1;cd(i)=" 999";id(i)=78808;lt(i)=  8.98;ln(i)=  -79.55;el(i)=  66;lnm(i)="ALBROOK                  ";st(i)="99";ct(i)="PM"
      ! https://ruc.noaa.gov/raobs/international_sites.txt
i=i+1;cd(i)="ENJA";id(i)=01001;lt(i)= 70.93;ln(i)=  -8.67;el(i)=   9;lnm(i)="JAN MAYEN(NOR-NAVY)      ";st(i)="  ";ct(i)="NO"
i=i+1;cd(i)="ENAS";id(i)=01004;lt(i)= 78.92;ln(i)=  11.93;el(i)=   8;lnm(i)="NY-ALESUND II            ";st(i)="  ";ct(i)="NO"
i=i+1;cd(i)="ENAN";id(i)=01010;lt(i)= 69.30;ln(i)=  16.15;el(i)=  14;lnm(i)="ANDOYA/ANDENES(AFB)      ";st(i)="  ";ct(i)="NO"
i=i+1;cd(i)="ENBJ";id(i)=01028;lt(i)= 74.52;ln(i)=  19.02;el(i)=  18;lnm(i)="BJORNOYA ISLAND          ";st(i)="  ";ct(i)="NO"
i=i+1;cd(i)="ENBO";id(i)=01152;lt(i)= 67.25;ln(i)=  14.40;el(i)=   8;lnm(i)="BODO VI (CIV/MIL)        ";st(i)="  ";ct(i)="NO"
i=i+1;cd(i)="ENOL";id(i)=01241;lt(i)= 63.70;ln(i)=   9.60;el(i)=   7;lnm(i)="ORLAND III(NOR-AFB)      ";st(i)="  ";ct(i)="NO"
i=i+1;cd(i)="ENEK";id(i)=01400;lt(i)= 56.53;ln(i)=   3.22;el(i)=  46;lnm(i)="EKOFISK                  ";st(i)="  ";ct(i)="NO"
i=i+1;cd(i)="ENZV";id(i)=01415;lt(i)= 58.87;ln(i)=   5.67;el(i)=  33;lnm(i)="STAVANGER/SOLA(AFB)      ";st(i)="  ";ct(i)="NO"
i=i+1;cd(i)="9999";id(i)=01492;lt(i)= 59.95;ln(i)=  10.72;el(i)=  96;lnm(i)="OSLO-BLINDERN            ";st(i)="  ";ct(i)="NO"
i=i+1;cd(i)="ESPA";id(i)=02185;lt(i)= 65.55;ln(i)=  22.13;el(i)=  16;lnm(i)="LULEA/KALLAX (AFB)       ";st(i)="  ";ct(i)="SE"
i=i+1;cd(i)="9999";id(i)=02225;lt(i)= 63.18;ln(i)=  14.50;el(i)= 366;lnm(i)="OSTERSUND                ";st(i)="  ";ct(i)="SE"
i=i+1;cd(i)="ESNN";id(i)=02365;lt(i)= 62.53;ln(i)=  17.47;el(i)=   6;lnm(i)="SUNDSVALL/HARNOSAND      ";st(i)="  ";ct(i)="SE"
i=i+1;cd(i)="ESGG";id(i)=02527;lt(i)= 57.67;ln(i)=  12.30;el(i)= 155;lnm(i)="GOTEBORG/LANDVETTER      ";st(i)="  ";ct(i)="SE"
i=i+1;cd(i)="ESQV";id(i)=02591;lt(i)= 57.65;ln(i)=  18.35;el(i)=  45;lnm(i)="VISBY AERO STATION       ";st(i)="  ";ct(i)="SE"
i=i+1;cd(i)="EFSO";id(i)=02836;lt(i)= 67.37;ln(i)=  26.65;el(i)= 178;lnm(i)="SODANKYLA                ";st(i)="  ";ct(i)="FI"
i=i+1;cd(i)="EFJY";id(i)=02935;lt(i)= 62.40;ln(i)=  25.67;el(i)= 139;lnm(i)="JYVASKYLA (MIL/CIV)      ";st(i)="  ";ct(i)="FI"
i=i+1;cd(i)="9999";id(i)=02963;lt(i)= 60.82;ln(i)=  23.50;el(i)= 103;lnm(i)="JOKIOINEN                ";st(i)="  ";ct(i)="FI"
i=i+1;cd(i)="9999";id(i)=03005;lt(i)= 60.13;ln(i)=  -1.18;el(i)=  82;lnm(i)="LERWICK/SHETLAND IS      ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=03023;lt(i)= 57.33;ln(i)=  -7.37;el(i)=  10;lnm(i)="SOUTH UIST RANGE         ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="EGPO";id(i)=03026;lt(i)= 58.22;ln(i)=  -6.32;el(i)=   9;lnm(i)="STORNOWAY                ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="EGQK";id(i)=03066;lt(i)= 57.65;ln(i)=  -3.57;el(i)=   7;lnm(i)="KINLOSS RAF              ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="EGQL";id(i)=03171;lt(i)= 56.38;ln(i)=  -2.87;el(i)=  12;lnm(i)="LEUCHARS RAF             ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=03238;lt(i)= 55.01;ln(i)=  -1.87;el(i)= 141;lnm(i)="ALBEMARLE                ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="EGQM";id(i)=03240;lt(i)= 55.42;ln(i)=  -1.60;el(i)=  23;lnm(i)="BOULMER                  ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="EGXE";id(i)=03257;lt(i)= 54.30;ln(i)=  -1.53;el(i)=  40;lnm(i)="LEEMING RAF              ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=03317;lt(i)= 53.57;ln(i)=  -3.05;el(i)=  11;lnm(i)="WOODVALE                 ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=03354;lt(i)= 53.00;ln(i)=  -1.25;el(i)= 117;lnm(i)="NOTTINGHAM WX CTR        ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="EGXW";id(i)=03377;lt(i)= 53.17;ln(i)=  -0.52;el(i)=  70;lnm(i)="WADDINGTON RAF           ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="EGOS";id(i)=03414;lt(i)= 52.80;ln(i)=  -2.67;el(i)=  76;lnm(i)="SHAWBURY RAF             ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=03495;lt(i)= 52.77;ln(i)=   1.35;el(i)=  17;lnm(i)="COLTISHALL               ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=03496;lt(i)= 52.68;ln(i)=   1.68;el(i)=  13;lnm(i)="HEMSBY                   ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="EGUC";id(i)=03502;lt(i)= 52.13;ln(i)=  -4.57;el(i)= 132;lnm(i)="ABERPORTH                ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="EGUW";id(i)=03590;lt(i)= 52.12;ln(i)=   0.97;el(i)=  87;lnm(i)="WATTISHAM RAF            ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="EGVN";id(i)=03649;lt(i)= 51.75;ln(i)=  -1.58;el(i)=  88;lnm(i)="BRIZE NORTON RAF         ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=03693;lt(i)= 51.55;ln(i)=   0.83;el(i)=   3;lnm(i)="SHOEBURYNESS             ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=03743;lt(i)= 51.20;ln(i)=  -1.80;el(i)= 132;lnm(i)="LARKHILL                 ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=03808;lt(i)= 50.22;ln(i)=  -5.32;el(i)=  87;lnm(i)="CAMBORNE                 ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=03882;lt(i)= 50.90;ln(i)=   0.32;el(i)=  54;lnm(i)="HERSTOMONCEUX WEST       ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=03918;lt(i)= 54.30;ln(i)=  -6.19;el(i)=  15;lnm(i)="CASTOR BAY               ";st(i)="  ";ct(i)="IE"
i=i+1;cd(i)="9999";id(i)=03920;lt(i)= 54.48;ln(i)=  -6.10;el(i)=  38;lnm(i)="HILLSBOROUGH             ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=03953;lt(i)= 51.93;ln(i)= -10.25;el(i)=   9;lnm(i)="VALENTIA OBSERVATORY     ";st(i)="  ";ct(i)="IE"
i=i+1;cd(i)="BIKF";id(i)=04018;lt(i)= 63.97;ln(i)= -22.60;el(i)=  54;lnm(i)="KEFLAVIK (CIV\NAS)       ";st(i)="  ";ct(i)="IS"
i=i+1;cd(i)="BIEG";id(i)=04089;lt(i)= 65.28;ln(i)= -14.37;el(i)=  38;lnm(i)="EGILSSTADIR              ";st(i)="  ";ct(i)="IS"
i=i+1;cd(i)="BGTL";id(i)=04202;lt(i)= 76.53;ln(i)= -68.75;el(i)=  77;lnm(i)="THULE AB                 ";st(i)="  ";ct(i)="GL"
i=i+1;cd(i)="BGEM";id(i)=04220;lt(i)= 68.70;ln(i)= -52.75;el(i)=  40;lnm(i)="EGEDESMINDE/AUSIAT       ";st(i)="  ";ct(i)="GL"
i=i+1;cd(i)="BGBW";id(i)=04270;lt(i)= 61.18;ln(i)= -45.43;el(i)=   4;lnm(i)="NARSSARSSUAQ             ";st(i)="  ";ct(i)="GL"
i=i+1;cd(i)="BGDH";id(i)=04320;lt(i)= 76.77;ln(i)= -18.67;el(i)=  12;lnm(i)="DANMARKSHAVN (PORT)      ";st(i)="  ";ct(i)="GL"
i=i+1;cd(i)="BGSC";id(i)=04339;lt(i)= 70.48;ln(i)= -21.97;el(i)=  66;lnm(i)="SCORESBYSUND             ";st(i)="  ";ct(i)="GL"
i=i+1;cd(i)="BGAM";id(i)=04360;lt(i)= 65.60;ln(i)= -37.63;el(i)=  50;lnm(i)="ANGMAGSSALIK             ";st(i)="  ";ct(i)="GL"
i=i+1;cd(i)="9999";id(i)=04417;lt(i)= 72.58;ln(i)= -38.46;el(i)=3255;lnm(i)="GEOSUMMIT                ";st(i)="  ";ct(i)="GL"
i=i+1;cd(i)="9999";id(i)=06011;lt(i)= 62.02;ln(i)=  -6.77;el(i)=  55;lnm(i)="THORSHAVN (PORT)         ";st(i)="  ";ct(i)="DK"
i=i+1;cd(i)="9999";id(i)=06030;lt(i)= 57.10;ln(i)=   9.87;el(i)=   3;lnm(i)="ALBORG                   ";st(i)="  ";ct(i)="DN"
i=i+1;cd(i)="EKKA";id(i)=06060;lt(i)= 56.30;ln(i)=   9.12;el(i)=  53;lnm(i)="KARUP                    ";st(i)="  ";ct(i)="DN"
i=i+1;cd(i)="9999";id(i)=06181;lt(i)= 55.77;ln(i)=  12.53;el(i)=  40;lnm(i)="KOEBENHAVN\JAEGERSBORG   ";st(i)="  ";ct(i)="DK"
i=i+1;cd(i)="9999";id(i)=06210;lt(i)= 52.17;ln(i)=   4.42;el(i)=   1;lnm(i)="VALKENBURG               ";st(i)="  ";ct(i)="NL"
i=i+1;cd(i)="EHKD";id(i)=06235;lt(i)= 52.92;ln(i)=   4.78;el(i)=  14;lnm(i)="DE KOOY (NAVY)           ";st(i)="  ";ct(i)="NL"
i=i+1;cd(i)="EHDB";id(i)=06260;lt(i)= 52.10;ln(i)=   5.18;el(i)=   2;lnm(i)="DE BILT                  ";st(i)="  ";ct(i)="NL"
i=i+1;cd(i)="EBFN";id(i)=06400;lt(i)= 51.08;ln(i)=   2.65;el(i)=   9;lnm(i)="KOKSIJDE (BEL-AFB)       ";st(i)="  ";ct(i)="BE"
i=i+1;cd(i)="9999";id(i)=06447;lt(i)= 50.80;ln(i)=   4.35;el(i)= 104;lnm(i)="UCCLE/UKKLE              ";st(i)="  ";ct(i)="BE"
i=i+1;cd(i)="9999";id(i)=06458;lt(i)= 50.75;ln(i)=   4.77;el(i)= 127;lnm(i)="BEAUVECHAIN              ";st(i)="  ";ct(i)="BE"
i=i+1;cd(i)="ESBU";id(i)=06476;lt(i)= 50.03;ln(i)=   5.40;el(i)= 558;lnm(i)="SAINT HUBERT (MIL)       ";st(i)="  ";ct(i)="BE"
i=i+1;cd(i)="EBLB";id(i)=06496;lt(i)= 50.47;ln(i)=   6.18;el(i)= 570;lnm(i)="ELSENBORN (MIL)          ";st(i)="  ";ct(i)="BE"
i=i+1;cd(i)="LSMP";id(i)=06610;lt(i)= 46.82;ln(i)=   6.95;el(i)= 490;lnm(i)="PAYERNE (MIL/AUT)        ";st(i)="  ";ct(i)="CH"
i=i+1;cd(i)="9999";id(i)=06842;lt(i)= 47.38;ln(i)=  -9.65;el(i)= 408;lnm(i)="DIEPOLDSAV               ";st(i)="  ";ct(i)="CH"
i=i+1;cd(i)="9999";id(i)=06843;lt(i)= 46.70;ln(i)=   9.42;el(i)= 865;lnm(i)="MASEIN                   ";st(i)="  ";ct(i)="CH"
i=i+1;cd(i)="LFRB";id(i)=07110;lt(i)= 48.45;ln(i)=  -4.42;el(i)=  99;lnm(i)="BREST/GUIPAVAS           ";st(i)="  ";ct(i)="FR"
i=i+1;cd(i)="9999";id(i)=07145;lt(i)= 48.77;ln(i)=   2.02;el(i)= 168;lnm(i)="TRAPPES (AUT)            ";st(i)="  ";ct(i)="FR"
i=i+1;cd(i)="LFSN";id(i)=07180;lt(i)= 48.68;ln(i)=   6.22;el(i)= 225;lnm(i)="NANCY/ESSEY              ";st(i)="  ";ct(i)="FR"
i=i+1;cd(i)="LFLL";id(i)=07481;lt(i)= 45.73;ln(i)=   5.08;el(i)= 248;lnm(i)="LYON/SATOLAS             ";st(i)="  ";ct(i)="FR"
i=i+1;cd(i)="LFBD";id(i)=07510;lt(i)= 44.83;ln(i)=  -0.70;el(i)=  45;lnm(i)="BORDEAUX/MERIGNAC        ";st(i)="  ";ct(i)="FR"
i=i+1;cd(i)="LFME";id(i)=07645;lt(i)= 43.87;ln(i)=   4.40;el(i)=  60;lnm(i)="NIMES/COURBESSAC         ";st(i)="  ";ct(i)="FR"
i=i+1;cd(i)="LFMN";id(i)=07690;lt(i)= 43.65;ln(i)=   7.20;el(i)=  10;lnm(i)="NICE/COTE D'AZUR         ";st(i)="  ";ct(i)="FR"
i=i+1;cd(i)="LFKJ";id(i)=07761;lt(i)= 41.92;ln(i)=   8.80;el(i)=   5;lnm(i)="AJACCIO/CAMPO ORO        ";st(i)="  ";ct(i)="FR"
i=i+1;cd(i)="9999";id(i)=08001;lt(i)= 43.37;ln(i)=  -8.42;el(i)=  67;lnm(i)="LA CORUNA CITY           ";st(i)="  ";ct(i)="ES"
i=i+1;cd(i)="9999";id(i)=08023;lt(i)= 43.47;ln(i)=  -3.82;el(i)=  65;lnm(i)="SANTANDER CITY           ";st(i)="  ";ct(i)="ES"
i=i+1;cd(i)="LEZG";id(i)=08160;lt(i)= 41.67;ln(i)=  -1.02;el(i)= 258;lnm(i)="ZARAGOZA (MIL/CIV)&      ";st(i)="  ";ct(i)="ES"
i=i+1;cd(i)="9999";id(i)=08190;lt(i)= 41.38;ln(i)=   2.12;el(i)=  91;lnm(i)="BARCELONA SERVEI         ";st(i)="  ";ct(i)="ES"
i=i+1;cd(i)="LEMD";id(i)=08221;lt(i)= 40.47;ln(i)=  -3.58;el(i)= 638;lnm(i)="MADRID/BARAJAS           ";st(i)="  ";ct(i)="ES"
i=i+1;cd(i)="9999";id(i)=08301;lt(i)= 39.55;ln(i)=   2.62;el(i)=   6;lnm(i)="PALMA DE MALLORCA        ";st(i)="  ";ct(i)="ES"
i=i+1;cd(i)="9999";id(i)=08302;lt(i)= 39.60;ln(i)=   2.70;el(i)=  45;lnm(i)="PALMA/SON BENET          ";st(i)="  ";ct(i)="ES"
i=i+1;cd(i)="9999";id(i)=08430;lt(i)= 38.00;ln(i)=  -1.17;el(i)=  62;lnm(i)="MURCIA CITY              ";st(i)="  ";ct(i)="ES"
i=i+1;cd(i)="LXGB";id(i)=08495;lt(i)= 36.15;ln(i)=  -5.33;el(i)=   3;lnm(i)="GIBRALTAR (CIV/MIL)      ";st(i)="  ";ct(i)="GI"
i=i+1;cd(i)="9999";id(i)=08507;lt(i)= 39.10;ln(i)= -28.03;el(i)=  28;lnm(i)="GRACIOSA AERODROME       ";st(i)="  ";ct(i)="PT"
i=i+1;cd(i)="9999";id(i)=08508;lt(i)= 38.73;ln(i)= -27.07;el(i)= 113;lnm(i)="LAJES/SANTA RITA         ";st(i)="  ";ct(i)="PT"
i=i+1;cd(i)="9999";id(i)=08522;lt(i)= 32.63;ln(i)= -16.90;el(i)=  56;lnm(i)="FUNCHAL                  ";st(i)="  ";ct(i)="PT"
i=i+1;cd(i)="9999";id(i)=08579;lt(i)= 38.77;ln(i)=  -9.13;el(i)= 105;lnm(i)="LISBON/GAGO COUTINHO     ";st(i)="  ";ct(i)="PT"
i=i+1;cd(i)="9999";id(i)=08589;lt(i)= 14.90;ln(i)= -23.52;el(i)=  27;lnm(i)="PRAIA                    ";st(i)="  ";ct(i)="CV"
i=i+1;cd(i)="GVAC";id(i)=08594;lt(i)= 16.73;ln(i)= -22.95;el(i)=  55;lnm(i)="SAL ISL/AMILCAR CAB      ";st(i)="  ";ct(i)="CV"
i=i+1;cd(i)="9999";id(i)=10035;lt(i)= 54.53;ln(i)=   9.55;el(i)=  48;lnm(i)="SCHLESWIG                ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="9999";id(i)=10113;lt(i)= 53.72;ln(i)=   7.15;el(i)=  16;lnm(i)="NORDERNEY                ";st(i)="  ";ct(i)="DL"
i=i+1;cd(i)="9999";id(i)=10184;lt(i)= 54.10;ln(i)=  13.40;el(i)=   6;lnm(i)="GREIFSWALD               ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="EDWE";id(i)=10200;lt(i)= 53.38;ln(i)=   7.23;el(i)=   1;lnm(i)="EMDEN - FLUGPLATZ        ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="ETGB";id(i)=10238;lt(i)= 52.82;ln(i)=   9.93;el(i)=  69;lnm(i)="BERGEN (MIL)             ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="9999";id(i)=10304;lt(i)= 52.73;ln(i)=   7.33;el(i)=  41;lnm(i)="MEPPEN                   ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="9999";id(i)=10393;lt(i)= 52.22;ln(i)=  14.12;el(i)= 115;lnm(i)="LIN                      ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="EDZE";id(i)=10410;lt(i)= 51.40;ln(i)=   6.97;el(i)= 153;lnm(i)="ESSEN/MULHEIM            ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="9999";id(i)=10437;lt(i)= 51.13;ln(i)=   9.28;el(i)= 223;lnm(i)="FRITZLAR - KASSELER WARTE";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="9999";id(i)=10468;lt(i)= 51.55;ln(i)=  12.07;el(i)= 106;lnm(i)="OPPIN                    ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="9999";id(i)=10548;lt(i)= 50.57;ln(i)=  10.37;el(i)= 450;lnm(i)="MEININGEN                ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="ETGI";id(i)=10618;lt(i)= 49.70;ln(i)=   7.33;el(i)= 376;lnm(i)="IDAR-OBERSTEIN(MIL)      ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="9999";id(i)=10739;lt(i)= 48.83;ln(i)=   9.20;el(i)= 315;lnm(i)="STUTTGART/SCHNARRENBERG  ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="ETGK";id(i)=10771;lt(i)= 49.43;ln(i)=  11.90;el(i)= 419;lnm(i)="KUEMMERSBRUCK            ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="9999";id(i)=10828;lt(i)= 48.10;ln(i)=   9.25;el(i)= 646;lnm(i)="SIGMARINGEN              ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="9999";id(i)=10868;lt(i)= 48.25;ln(i)=  11.55;el(i)= 489;lnm(i)="MUENCHEN/OBERSCHLEISSHEIM";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="9999";id(i)=10954;lt(i)= 47.83;ln(i)=  10.87;el(i)= 757;lnm(i)="ALTENSTADT/SHONGUA       ";st(i)="  ";ct(i)="DL"
i=i+1;cd(i)="9999";id(i)=10962;lt(i)= 47.80;ln(i)=  11.02;el(i)= 986;lnm(i)="HOHENPEISSENBERG         ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="LOWL";id(i)=11010;lt(i)= 48.23;ln(i)=  14.20;el(i)= 313;lnm(i)="LINZ (CIV/MIL)           ";st(i)="  ";ct(i)="AT"
i=i+1;cd(i)="9999";id(i)=11035;lt(i)= 48.25;ln(i)=  16.37;el(i)= 200;lnm(i)="WIEN/HOHE WARTE          ";st(i)="  ";ct(i)="AT"
i=i+1;cd(i)="LOWI";id(i)=11120;lt(i)= 47.27;ln(i)=  11.35;el(i)= 581;lnm(i)="INNSBRUCK AIRPORT        ";st(i)="  ";ct(i)="AT"
i=i+1;cd(i)="LOWG";id(i)=11240;lt(i)= 47.00;ln(i)=  15.43;el(i)= 347;lnm(i)="GRAZ (MIL/CIV)           ";st(i)="  ";ct(i)="AT"
i=i+1;cd(i)="9999";id(i)=11520;lt(i)= 50.00;ln(i)=  14.45;el(i)= 303;lnm(i)="PRAGUE/LIBUS             ";st(i)="  ";ct(i)="CZ"
i=i+1;cd(i)="9999";id(i)=11722;lt(i)= 49.08;ln(i)=  16.62;el(i)= 195;lnm(i)="BRNO REBESOVICE          ";st(i)="  ";ct(i)="CZ"
i=i+1;cd(i)="9999";id(i)=11747;lt(i)= 49.45;ln(i)=  17.13;el(i)= 216;lnm(i)="PROSTEJOV                ";st(i)="  ";ct(i)="CZ"
i=i+1;cd(i)="9999";id(i)=11952;lt(i)= 49.03;ln(i)=  20.32;el(i)= 701;lnm(i)="POPRAD/GANOVCE           ";st(i)="  ";ct(i)="SK"
i=i+1;cd(i)="9999";id(i)=12120;lt(i)= 54.75;ln(i)=  17.53;el(i)=   2;lnm(i)="LEBA                     ";st(i)="  ";ct(i)="PL"
i=i+1;cd(i)="9999";id(i)=12374;lt(i)= 52.40;ln(i)=  20.97;el(i)=  96;lnm(i)="LEGIONOWO                ";st(i)="  ";ct(i)="PL"
i=i+1;cd(i)="9999";id(i)=12425;lt(i)= 51.13;ln(i)=  16.98;el(i)= 116;lnm(i)="WROCLAW/MALY GADOW       ";st(i)="  ";ct(i)="PL"
i=i+1;cd(i)="9999";id(i)=12843;lt(i)= 47.43;ln(i)=  19.18;el(i)= 139;lnm(i)="BUDAPEST/LORINC          ";st(i)="  ";ct(i)="HU"
i=i+1;cd(i)="LHUD";id(i)=12982;lt(i)= 46.25;ln(i)=  20.10;el(i)=  84;lnm(i)="SZEGED (AUT)             ";st(i)="  ";ct(i)="HU"
i=i+1;cd(i)="9999";id(i)=13275;lt(i)= 44.77;ln(i)=  20.42;el(i)= 203;lnm(i)="BEOGRAD KOSUTNJAK        ";st(i)="  ";ct(i)="YU"
i=i+1;cd(i)="LWSK";id(i)=13586;lt(i)= 41.97;ln(i)=  21.65;el(i)= 239;lnm(i)="SKOPJE/PETROVAC          ";st(i)="  ";ct(i)="YU"
i=i+1;cd(i)="9999";id(i)=14015;lt(i)= 46.07;ln(i)=  14.52;el(i)= 316;lnm(i)="LJUBLJANA\BEZIGRAD       ";st(i)="  ";ct(i)="SI"
i=i+1;cd(i)="LDDD";id(i)=14240;lt(i)= 45.82;ln(i)=  16.03;el(i)= 123;lnm(i)="ZAGREB\MAKSIMIR          ";st(i)="  ";ct(i)="HR"
i=i+1;cd(i)="9999";id(i)=14430;lt(i)= 44.09;ln(i)=  15.35;el(i)=  84;lnm(i)="ZADAR                    ";st(i)="  ";ct(i)="CR"
i=i+1;cd(i)="LRCL";id(i)=15120;lt(i)= 46.78;ln(i)=  23.57;el(i)= 413;lnm(i)="CLUJ-NAPOCA/SOMESEN      ";st(i)="  ";ct(i)="RO"
i=i+1;cd(i)="LRBS";id(i)=15420;lt(i)= 44.50;ln(i)=  26.13;el(i)=  91;lnm(i)="BUCHAREST/BANEASA        ";st(i)="  ";ct(i)="RO"
i=i+1;cd(i)="9999";id(i)=15480;lt(i)= 44.18;ln(i)=  28.67;el(i)=  14;lnm(i)="CONSTANTA                ";st(i)="  ";ct(i)="RO"
i=i+1;cd(i)="LBSF";id(i)=15614;lt(i)= 42.65;ln(i)=  23.38;el(i)= 595;lnm(i)="SOFIA (OBSERVATORY)      ";st(i)="  ";ct(i)="BG"
i=i+1;cd(i)="LIPD";id(i)=16044;lt(i)= 46.04;ln(i)=  13.18;el(i)=  94;lnm(i)="UDINE/CAMPOFORMIDO       ";st(i)="  ";ct(i)="IT"
i=i+1;cd(i)="LIPI";id(i)=16045;lt(i)= 45.97;ln(i)=  13.05;el(i)=  52;lnm(i)="UDINE/RIVOLTO            ";st(i)="  ";ct(i)="IT"
i=i+1;cd(i)="LIML";id(i)=16080;lt(i)= 45.43;ln(i)=   9.28;el(i)= 103;lnm(i)="MILANO/LINATE            ";st(i)="  ";ct(i)="IT"
i=i+1;cd(i)="9999";id(i)=16087;lt(i)= 45.45;ln(i)=  11.00;el(i)=  90;lnm(i)="VERONA                   ";st(i)="  ";ct(i)="IT"
i=i+1;cd(i)="9999";id(i)=16113;lt(i)= 44.53;ln(i)=   7.62;el(i)= 386;lnm(i)="CUNEO - LEVALDIGI        ";st(i)="  ";ct(i)="IT"
i=i+1;cd(i)="9999";id(i)=16121;lt(i)= 44.40;ln(i)=   8.90;el(i)=  45;lnm(i)="GENOVA                   ";st(i)="  ";ct(i)="IT"
i=i+1;cd(i)="9999";id(i)=16144;lt(i)= 44.65;ln(i)=  11.62;el(i)=  38;lnm(i)="S. PIETRO CAPOFIUME      ";st(i)="  ";ct(i)="IT"
i=i+1;cd(i)="LIRE";id(i)=16245;lt(i)= 41.65;ln(i)=  12.43;el(i)=  12;lnm(i)="PRATICA DI MARE(AB)      ";st(i)="  ";ct(i)="IT"
i=i+1;cd(i)="LIBR";id(i)=16320;lt(i)= 40.65;ln(i)=  17.95;el(i)=  10;lnm(i)="BRINDISI/CASALE AFB      ";st(i)="  ";ct(i)="IT"
i=i+1;cd(i)="LICT";id(i)=16429;lt(i)= 37.92;ln(i)=  12.50;el(i)=  14;lnm(i)="TRAPANI/BIRGI (AFB)      ";st(i)="  ";ct(i)="IT"
i=i+1;cd(i)="LIED";id(i)=16546;lt(i)= 39.35;ln(i)=   8.97;el(i)=  28;lnm(i)="DECIMOMANNU              ";st(i)="  ";ct(i)="IY"
i=i+1;cd(i)="LIEE";id(i)=16560;lt(i)= 39.25;ln(i)=   9.05;el(i)=   7;lnm(i)="CAGLIARI/ELMAS(AFB)      ";st(i)="  ";ct(i)="IT"
i=i+1;cd(i)="LGTS";id(i)=16622;lt(i)= 40.52;ln(i)=  22.97;el(i)=   4;lnm(i)="THESSALONIKI/MIKRA       ";st(i)="  ";ct(i)="GR"
i=i+1;cd(i)="LGAT";id(i)=16716;lt(i)= 37.90;ln(i)=  23.73;el(i)=  14;lnm(i)="ATHENS/HELLENKION        ";st(i)="  ";ct(i)="GR"
i=i+1;cd(i)="LGIR";id(i)=16754;lt(i)= 35.33;ln(i)=  25.18;el(i)=  39;lnm(i)="IRAKLION (CIV/AFB)       ";st(i)="  ";ct(i)="GR"
i=i+1;cd(i)="9999";id(i)=17030;lt(i)= 41.28;ln(i)=  36.33;el(i)=   4;lnm(i)="SAMSUN CITY              ";st(i)="  ";ct(i)="TR"
i=i+1;cd(i)="9999";id(i)=17062;lt(i)= 40.97;ln(i)=  29.08;el(i)=  33;lnm(i)="ISTANBUL/GOZTEPE         ";st(i)="  ";ct(i)="TR"
i=i+1;cd(i)="9999";id(i)=17064;lt(i)= 40.90;ln(i)=  29.15;el(i)=  18;lnm(i)="ISTANBUL BOLGE (KARTAL)  ";st(i)="  ";ct(i)="TR"
i=i+1;cd(i)="9999";id(i)=17095;lt(i)= 39.90;ln(i)=  41.28;el(i)=1952;lnm(i)="ERZURUM BOLGE            ";st(i)="  ";ct(i)="TR"
i=i+1;cd(i)="9999";id(i)=17130;lt(i)= 39.95;ln(i)=  32.88;el(i)= 894;lnm(i)="ANKARA/CENTRAL           ";st(i)="  ";ct(i)="TR"
i=i+1;cd(i)="9999";id(i)=17220;lt(i)= 38.43;ln(i)=  27.17;el(i)=  25;lnm(i)="IZMIR                    ";st(i)="  ";ct(i)="TR"
i=i+1;cd(i)="LTBM";id(i)=17240;lt(i)= 37.75;ln(i)=  30.55;el(i)= 997;lnm(i)="ISPARTA                  ";st(i)="  ";ct(i)="TR"
i=i+1;cd(i)="LTCC";id(i)=17281;lt(i)= 37.88;ln(i)=  40.18;el(i)= 677;lnm(i)="DIYARBAKIR(CIV/AFB)      ";st(i)="  ";ct(i)="TR"
i=i+1;cd(i)="9999";id(i)=17351;lt(i)= 37.05;ln(i)=  35.35;el(i)=  28;lnm(i)="ADANA                    ";st(i)="  ";ct(i)="TR"
i=i+1;cd(i)="LTAF";id(i)=17352;lt(i)= 36.98;ln(i)=  35.30;el(i)=  20;lnm(i)="ADANA/SAKIRPASA          ";st(i)="  ";ct(i)="TR"
i=i+1;cd(i)="9999";id(i)=17600;lt(i)= 34.72;ln(i)=  32.48;el(i)=  11;lnm(i)="PAPHOS                   ";st(i)="  ";ct(i)="CY"
i=i+1;cd(i)="LCRA";id(i)=17601;lt(i)= 34.58;ln(i)=  32.98;el(i)=  23;lnm(i)="AKROTIRI (RAF)           ";st(i)="  ";ct(i)="CY"
i=i+1;cd(i)="LCNC";id(i)=17607;lt(i)= 35.15;ln(i)=  33.40;el(i)= 161;lnm(i)="NICOSIA/ATHALASSA        ";st(i)="  ";ct(i)="CY"
i=i+1;cd(i)="9999";id(i)=17609;lt(i)= 34.88;ln(i)=  33.63;el(i)=   2;lnm(i)="LARNACA                  ";st(i)="  ";ct(i)="CY"
i=i+1;cd(i)="9999";id(i)=20046;lt(i)= 80.62;ln(i)=  58.05;el(i)=  22;lnm(i)="POLARGMO                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=20292;lt(i)= 77.72;ln(i)= 104.30;el(i)=  15;lnm(i)="FEDOROVA                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=20674;lt(i)= 73.53;ln(i)=  80.40;el(i)=  47;lnm(i)="DIKSON ISLAND            ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=20744;lt(i)= 72.38;ln(i)=  52.73;el(i)=  16;lnm(i)="MALYE KARMAKULY          ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=20891;lt(i)= 72.00;ln(i)= 102.45;el(i)=  33;lnm(i)="HATANGA                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=20892;lt(i)= 72.00;ln(i)= 102.57;el(i)=  26;lnm(i)="HATANGA                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=21432;lt(i)= 76.00;ln(i)= 137.90;el(i)=  10;lnm(i)="KOTEL'NYJ ISLAND         ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=21824;lt(i)= 71.58;ln(i)= 128.92;el(i)=   8;lnm(i)="TIKSI                    ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=21946;lt(i)= 70.62;ln(i)= 147.90;el(i)=  61;lnm(i)="COKURDAH                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=21982;lt(i)= 70.98;ln(i)= 178.48;el(i)=   3;lnm(i)="VRANGELJA ISLAND         ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="ULMM";id(i)=22113;lt(i)= 68.97;ln(i)=  33.05;el(i)=  51;lnm(i)="MURMANSK                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=22217;lt(i)= 67.15;ln(i)=  32.35;el(i)=  26;lnm(i)="KANDALAKSA               ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=22271;lt(i)= 67.88;ln(i)=  44.13;el(i)=  16;lnm(i)="SOJNA                    ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=22522;lt(i)= 64.98;ln(i)=  34.80;el(i)=  10;lnm(i)="KEM (PORT)               ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=22543;lt(i)= 64.62;ln(i)=  40.51;el(i)=   6;lnm(i)="ARHANGEL'SK              ";st(i)="  ";ct(i)="RS"
i=i+1;cd(i)="ULAA";id(i)=22550;lt(i)= 64.53;ln(i)=  40.58;el(i)=   5;lnm(i)="ARHANGEL'SK              ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=22820;lt(i)= 61.82;ln(i)=  34.27;el(i)= 109;lnm(i)="PETROZAVODSK             ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=22845;lt(i)= 61.50;ln(i)=  38.93;el(i)= 121;lnm(i)="KARGOPOL'                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=23078;lt(i)= 69.33;ln(i)=  88.10;el(i)=  62;lnm(i)="NORILSK                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=23205;lt(i)= 67.65;ln(i)=  53.02;el(i)=   7;lnm(i)="NAR'JAN-MAR              ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=23330;lt(i)= 66.53;ln(i)=  66.53;el(i)=  16;lnm(i)="SALEHARD                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=23415;lt(i)= 65.11;ln(i)=  57.10;el(i)=  58;lnm(i)="PECORA                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=23418;lt(i)= 65.12;ln(i)=  57.10;el(i)=  56;lnm(i)="PECORA                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=23472;lt(i)= 65.78;ln(i)=  87.95;el(i)=  32;lnm(i)="TURUHANSK                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UUYY";id(i)=23802;lt(i)= 61.67;ln(i)=  50.78;el(i)= 117;lnm(i)="SYKTYVKAR                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UUYY";id(i)=23804;lt(i)= 61.72;ln(i)=  50.83;el(i)= 119;lnm(i)="SYKTYVKAR                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=23884;lt(i)= 61.60;ln(i)=  90.00;el(i)=  63;lnm(i)="BOR                      ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=23921;lt(i)= 60.68;ln(i)=  60.43;el(i)= 101;lnm(i)="IVDEL'                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="USHH";id(i)=23933;lt(i)= 60.97;ln(i)=  69.07;el(i)=  40;lnm(i)="HANTY-MANSIJSK           ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=23955;lt(i)= 60.43;ln(i)=  77.87;el(i)=  47;lnm(i)="ALEKSANDROVSKOE          ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=24125;lt(i)= 68.50;ln(i)= 112.43;el(i)= 220;lnm(i)="OLENEK                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=24266;lt(i)= 67.55;ln(i)= 133.38;el(i)= 137;lnm(i)="VERHOJANSK               ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=24343;lt(i)= 66.77;ln(i)= 123.40;el(i)=  93;lnm(i)="ZIGANSK                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=24507;lt(i)= 64.27;ln(i)= 100.23;el(i)= 186;lnm(i)="TURA                     ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=24641;lt(i)= 63.77;ln(i)= 121.62;el(i)= 107;lnm(i)="VILJUJSK                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=24688;lt(i)= 63.27;ln(i)= 143.15;el(i)= 745;lnm(i)="OJMJAKON                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=24726;lt(i)= 62.55;ln(i)= 113.88;el(i)= 347;lnm(i)="MIRNY                    ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=24908;lt(i)= 60.33;ln(i)= 102.27;el(i)= 260;lnm(i)="VANAVARA                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=24944;lt(i)= 60.40;ln(i)= 120.42;el(i)= 226;lnm(i)="OLEKMINSK                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UEEE";id(i)=24959;lt(i)= 62.08;ln(i)= 129.75;el(i)= 103;lnm(i)="YAKUTSK AIRPORT          ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=25123;lt(i)= 68.80;ln(i)= 161.28;el(i)=  32;lnm(i)="CERSKIJ                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=25400;lt(i)= 65.73;ln(i)= 150.90;el(i)=  43;lnm(i)="ZYRJANKA                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=25428;lt(i)= 65.23;ln(i)= 160.50;el(i)= 265;lnm(i)="SCERBAKOVO               ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UHMA";id(i)=25563;lt(i)= 64.78;ln(i)= 177.57;el(i)=  62;lnm(i)="ANADYR'                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=25703;lt(i)= 62.92;ln(i)= 152.42;el(i)= 207;lnm(i)="SEJMCAN                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UHMM";id(i)=25913;lt(i)= 59.58;ln(i)= 150.78;el(i)= 118;lnm(i)="MAGADAN                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=25954;lt(i)= 60.35;ln(i)= 166.00;el(i)=   3;lnm(i)="KORF                     ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="ULTT";id(i)=26038;lt(i)= 59.45;ln(i)=  24.80;el(i)=  44;lnm(i)="TALLIN                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="ULLI";id(i)=26063;lt(i)= 59.97;ln(i)=  30.30;el(i)=   4;lnm(i)="ST. PETERBURG            ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=26298;lt(i)= 57.90;ln(i)=  34.05;el(i)= 178;lnm(i)="BOLOGOE                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UMRR";id(i)=26422;lt(i)= 56.97;ln(i)=  24.07;el(i)=   7;lnm(i)="RIGA                     ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="ULOL";id(i)=26477;lt(i)= 56.38;ln(i)=  30.60;el(i)=  98;lnm(i)="VELIKIE LUKI             ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=26629;lt(i)= 54.88;ln(i)=  23.88;el(i)=  75;lnm(i)="KAUNAS                   ";st(i)="  ";ct(i)="BY"
i=i+1;cd(i)="9999";id(i)=26702;lt(i)= 54.70;ln(i)=  20.62;el(i)=  21;lnm(i)="KALININGRAD              ";st(i)="  ";ct(i)="BY"
i=i+1;cd(i)="9999";id(i)=26781;lt(i)= 54.75;ln(i)=  32.07;el(i)= 241;lnm(i)="SMOLENSK                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="ULWW";id(i)=27037;lt(i)= 59.23;ln(i)=  39.87;el(i)= 131;lnm(i)="VOLOGDA                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="ULWW";id(i)=27038;lt(i)= 59.32;ln(i)=  39.93;el(i)= 130;lnm(i)="VOLDOGA                  ";st(i)="  ";ct(i)="RS"
i=i+1;cd(i)="9999";id(i)=27199;lt(i)= 58.60;ln(i)=  49.63;el(i)= 158;lnm(i)="KIROV                    ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=27459;lt(i)= 56.27;ln(i)=  44.00;el(i)= 157;lnm(i)="NIZNIJ NOVGOROD          ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=27595;lt(i)= 55.78;ln(i)=  49.18;el(i)= 116;lnm(i)="KAZAN'                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=27612;lt(i)= 55.75;ln(i)=  37.57;el(i)= 156;lnm(i)="MOSCOW                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=27707;lt(i)= 54.12;ln(i)=  35.33;el(i)= 239;lnm(i)="SUHINICI                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=27730;lt(i)= 54.63;ln(i)=  39.70;el(i)= 158;lnm(i)="RJAZAN'                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UWPP";id(i)=27962;lt(i)= 53.13;ln(i)=  45.02;el(i)= 174;lnm(i)="PENZA                    ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=27995;lt(i)= 52.98;ln(i)=  49.43;el(i)=  45;lnm(i)="BEZENCUKSKAJA            ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=28225;lt(i)= 57.95;ln(i)=  56.20;el(i)= 170;lnm(i)="PERM                     ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=28275;lt(i)= 58.15;ln(i)=  68.18;el(i)=  44;lnm(i)="TOBOLSK                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=28445;lt(i)= 56.73;ln(i)=  61.07;el(i)= 290;lnm(i)="VERHNEE DUBROVO          ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=28661;lt(i)= 55.47;ln(i)=  65.40;el(i)=  79;lnm(i)="KURGAN                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=28698;lt(i)= 54.93;ln(i)=  73.40;el(i)= 123;lnm(i)="OMSK                     ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=28722;lt(i)= 54.75;ln(i)=  56.00;el(i)= 105;lnm(i)="UFA                      ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=28951;lt(i)= 59.32;ln(i)=  63.62;el(i)= 171;lnm(i)="KOSTANAI                 ";st(i)="  ";ct(i)="RA"
i=i+1;cd(i)="9999";id(i)=28952;lt(i)= 53.22;ln(i)=  63.62;el(i)= 156;lnm(i)="KUSTANAJ                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=29231;lt(i)= 58.30;ln(i)=  82.90;el(i)=  76;lnm(i)="KOLPASEV                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=29263;lt(i)= 58.45;ln(i)=  92.15;el(i)=  78;lnm(i)="YENISEYSK                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=29282;lt(i)= 58.42;ln(i)=  97.40;el(i)= 134;lnm(i)="BOGUCANY                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=29572;lt(i)= 56.18;ln(i)=  92.62;el(i)= 296;lnm(i)="EMEL'JANOVO              ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=29612;lt(i)= 55.37;ln(i)=  78.37;el(i)= 120;lnm(i)="BARABINSK                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UNNN";id(i)=29634;lt(i)= 55.03;ln(i)=  82.90;el(i)= 176;lnm(i)="NOVOSIBIRSK              ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UINN";id(i)=29698;lt(i)= 54.88;ln(i)=  99.03;el(i)= 410;lnm(i)="NIZNE-UDINSK             ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=29839;lt(i)= 53.35;ln(i)=  83.82;el(i)= 159;lnm(i)="BARNAUL                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=29862;lt(i)= 53.77;ln(i)=  91.32;el(i)= 256;lnm(i)="HAKASSKATA               ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=30054;lt(i)= 59.45;ln(i)= 112.58;el(i)= 193;lnm(i)="VITIM                    ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UIKK";id(i)=30230;lt(i)= 57.77;ln(i)= 108.12;el(i)= 258;lnm(i)="KIRENSK                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=30309;lt(i)= 56.07;ln(i)= 101.83;el(i)= 489;lnm(i)="BRATSK                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=30372;lt(i)= 56.92;ln(i)= 118.37;el(i)= 711;lnm(i)="CARA                     ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UIUB";id(i)=30554;lt(i)= 54.47;ln(i)= 113.13;el(i)= 995;lnm(i)="BOGDARIN                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=30635;lt(i)= 53.42;ln(i)= 109.02;el(i)= 457;lnm(i)="UST-BARGUZIN             ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=30673;lt(i)= 53.73;ln(i)= 119.78;el(i)= 619;lnm(i)="MOGOCA                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=30715;lt(i)= 52.48;ln(i)= 103.85;el(i)= 450;lnm(i)="ANGARSK                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=30758;lt(i)= 52.02;ln(i)= 113.33;el(i)= 685;lnm(i)="CHITA                    ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=30935;lt(i)= 50.37;ln(i)= 108.75;el(i)= 770;lnm(i)="KRASNYJ CIKOJ            ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=30965;lt(i)= 50.38;ln(i)= 116.52;el(i)= 684;lnm(i)="BORZJA                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=31004;lt(i)= 58.62;ln(i)= 125.37;el(i)= 682;lnm(i)="ALDAN                    ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=31088;lt(i)= 59.37;ln(i)= 143.20;el(i)=   6;lnm(i)="OHOTSK                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=31168;lt(i)= 56.45;ln(i)= 138.15;el(i)=   9;lnm(i)="AJAN                     ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=31300;lt(i)= 53.75;ln(i)= 127.23;el(i)= 232;lnm(i)="ZEJA                     ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=31369;lt(i)= 53.15;ln(i)= 140.70;el(i)=  68;lnm(i)="NIKOLAEVSK-ON-AM         ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=31510;lt(i)= 50.27;ln(i)= 127.50;el(i)= 137;lnm(i)="BLAGOVESCENSK            ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=31538;lt(i)= 50.07;ln(i)= 132.13;el(i)= 349;lnm(i)="SUTUR                    ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=31736;lt(i)= 48.53;ln(i)= 135.23;el(i)=  72;lnm(i)="HABAROVSK                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=31770;lt(i)= 48.97;ln(i)= 140.30;el(i)=  22;lnm(i)="SOVETSKAJA GAVAN         ";st(i)="  ";ct(i)="RA"
i=i+1;cd(i)="9999";id(i)=31873;lt(i)= 45.87;ln(i)= 133.73;el(i)= 107;lnm(i)="DAL'NERECENSK            ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=31977;lt(i)= 43.27;ln(i)= 132.05;el(i)=  82;lnm(i)="SAD-GOROD                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=32061;lt(i)= 50.90;ln(i)= 142.17;el(i)=  31;lnm(i)="ALEKSANDROVSK-SAHALINSKIJ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=32098;lt(i)= 49.22;ln(i)= 143.10;el(i)=   4;lnm(i)="PORONAJSK                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UHSS";id(i)=32150;lt(i)= 46.92;ln(i)= 142.73;el(i)=  31;lnm(i)="JUZNO-SAHALINSK          ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=32215;lt(i)= 50.68;ln(i)= 156.13;el(i)=  23;lnm(i)="SEVERO-KURILSK           ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=32389;lt(i)= 56.32;ln(i)= 160.83;el(i)=  28;lnm(i)="KLJUCI                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=32477;lt(i)= 54.30;ln(i)= 155.97;el(i)=  25;lnm(i)="SOBOLEVO                 ";st(i)="  ";ct(i)="RA"
i=i+1;cd(i)="UHPP";id(i)=32540;lt(i)= 52.97;ln(i)= 158.75;el(i)=  24;lnm(i)="PETROPAVLOVSK-KAMCHATSKIJ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=32618;lt(i)= 55.20;ln(i)= 165.98;el(i)=  14;lnm(i)="OSTROV BERINGA           ";st(i)="  ";ct(i)="RA"
i=i+1;cd(i)="UMBB";id(i)=33008;lt(i)= 52.12;ln(i)=  23.68;el(i)= 144;lnm(i)="BREST                    ";st(i)="  ";ct(i)="BY"
i=i+1;cd(i)="9999";id(i)=33041;lt(i)= 52.45;ln(i)=  31.00;el(i)= 127;lnm(i)="GOMEL'                   ";st(i)="  ";ct(i)="BY"
i=i+1;cd(i)="9999";id(i)=33317;lt(i)= 50.17;ln(i)=  27.05;el(i)= 278;lnm(i)="SEPETOVKA                ";st(i)="  ";ct(i)="UA"
i=i+1;cd(i)="UKKK";id(i)=33345;lt(i)= 50.40;ln(i)=  30.45;el(i)= 167;lnm(i)="KIEV/ZHULYANY            ";st(i)="  ";ct(i)="UA"
i=i+1;cd(i)="UKLL";id(i)=33393;lt(i)= 49.82;ln(i)=  23.95;el(i)= 325;lnm(i)="L'VIV                    ";st(i)="  ";ct(i)="UA"
i=i+1;cd(i)="9999";id(i)=33631;lt(i)= 48.63;ln(i)=  22.27;el(i)= 118;lnm(i)="UZGOROD                  ";st(i)="  ";ct(i)="UA"
i=i+1;cd(i)="9999";id(i)=33658;lt(i)= 48.37;ln(i)=  25.90;el(i)= 246;lnm(i)="CHERNIVTSI               ";st(i)="  ";ct(i)="UA"
i=i+1;cd(i)="UKDR";id(i)=33791;lt(i)= 47.92;ln(i)=  33.22;el(i)= 100;lnm(i)="LOZUVATKA                ";st(i)="  ";ct(i)="UA"
i=i+1;cd(i)="UKII";id(i)=33815;lt(i)= 47.02;ln(i)=  28.87;el(i)= 173;lnm(i)="KISHINEV                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UKOO";id(i)=33837;lt(i)= 46.43;ln(i)=  30.77;el(i)=  43;lnm(i)="ODESSA/TSENTRALNY        ";st(i)="  ";ct(i)="UA"
i=i+1;cd(i)="UKFF";id(i)=33946;lt(i)= 45.02;ln(i)=  33.98;el(i)= 187;lnm(i)="SIMFEROPOL'              ";st(i)="  ";ct(i)="UA"
i=i+1;cd(i)="9999";id(i)=33966;lt(i)= 45.05;ln(i)=  34.60;el(i)= 183;lnm(i)="KRYMSKAYA                ";st(i)="  ";ct(i)="UA"
i=i+1;cd(i)="9999";id(i)=34009;lt(i)= 51.77;ln(i)=  36.17;el(i)= 247;lnm(i)="KURSK                    ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UUOO";id(i)=34122;lt(i)= 51.67;ln(i)=  39.27;el(i)= 104;lnm(i)="VORONEZ                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UWSS";id(i)=34172;lt(i)= 51.57;ln(i)=  46.05;el(i)= 156;lnm(i)="SARATOV                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=34247;lt(i)= 50.42;ln(i)=  41.05;el(i)=  93;lnm(i)="KALAC                    ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=34300;lt(i)= 49.93;ln(i)=  36.28;el(i)= 152;lnm(i)="KHARKOV                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="URWW";id(i)=34560;lt(i)= 48.68;ln(i)=  44.35;el(i)= 145;lnm(i)="VOLGOGRAD                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="URWW";id(i)=34467;lt(i)= 48.58;ln(i)=  44.35;el(i)= 141;lnm(i)="VOLGOGRAD                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="URRR";id(i)=34731;lt(i)= 47.25;ln(i)=  39.82;el(i)=  77;lnm(i)="ROSTOV-NA-DONU           ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=34858;lt(i)= 45.92;ln(i)=  43.35;el(i)=  87;lnm(i)="DIVNOE                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=34880;lt(i)= 46.27;ln(i)=  48.03;el(i)= -18;lnm(i)="ASTRAHAN'                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=34882;lt(i)= 46.27;ln(i)=  47.98;el(i)= -17;lnm(i)="ASTRAHAN'                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=35121;lt(i)= 51.78;ln(i)=  55.22;el(i)= 109;lnm(i)="ORENBURG                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UATT";id(i)=35229;lt(i)= 50.28;ln(i)=  57.15;el(i)= 224;lnm(i)="AKTJUBINSK               ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=35394;lt(i)= 49.80;ln(i)=  73.13;el(i)= 553;lnm(i)="KARAGANDA                ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=35671;lt(i)= 47.80;ln(i)=  67.72;el(i)= 345;lnm(i)="DZEZKAZGAN               ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=35700;lt(i)= 47.02;ln(i)=  51.85;el(i)= -15;lnm(i)="GUR'EV                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=36003;lt(i)= 52.28;ln(i)=  76.95;el(i)= 123;lnm(i)="PAVLODAR                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=36096;lt(i)= 51.67;ln(i)=  94.38;el(i)= 629;lnm(i)="KYZYL                    ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UAAA";id(i)=36870;lt(i)= 43.23;ln(i)=  76.93;el(i)= 847;lnm(i)="ALMA-ATA                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=36872;lt(i)= 43.36;ln(i)=  77.00;el(i)= 663;lnm(i)="ALMATY                   ";st(i)="  ";ct(i)="RA"
i=i+1;cd(i)="9999";id(i)=37011;lt(i)= 44.09;ln(i)=  39.03;el(i)=  95;lnm(i)="TUAPSE (STREAM)          ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=37018;lt(i)= 44.10;ln(i)=  39.07;el(i)=  41;lnm(i)="TUAPSE (STREAM)          ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="URMM";id(i)=37054;lt(i)= 44.22;ln(i)=  43.10;el(i)= 314;lnm(i)="MINERALYNE VODY          ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="URMM";id(i)=37055;lt(i)= 44.22;ln(i)=  43.10;el(i)= 313;lnm(i)="MINERALYNE VODY          ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=37099;lt(i)= 43.58;ln(i)=  39.77;el(i)=  57;lnm(i)="RAZDOL/SOCI              ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=37259;lt(i)= 43.01;ln(i)=  47.48;el(i)=   2;lnm(i)="MAHACHKALA               ";st(i)="  ";ct(i)="RS"
i=i+1;cd(i)="UGGG";id(i)=37549;lt(i)= 41.68;ln(i)=  44.95;el(i)= 490;lnm(i)="TBILISI/NOVOALEXEYE      ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UGEE";id(i)=37789;lt(i)= 40.13;ln(i)=  44.47;el(i)= 890;lnm(i)="YEREVAN/ZVARTNOTS        ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=37860;lt(i)= 40.53;ln(i)=  50.00;el(i)=  28;lnm(i)="MASHTAGA                 ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=38062;lt(i)= 44.85;ln(i)=  65.50;el(i)= 130;lnm(i)="KZYLORDA                 ";st(i)="  ";ct(i)="RA"
i=i+1;cd(i)="9999";id(i)=38064;lt(i)= 44.77;ln(i)=  65.52;el(i)= 133;lnm(i)="KYZYLORDA                ";st(i)="  ";ct(i)="RA"
i=i+1;cd(i)="9999";id(i)=38341;lt(i)= 42.85;ln(i)=  71.38;el(i)= 653;lnm(i)="DZAMBUL                  ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UAFF";id(i)=38353;lt(i)= 42.85;ln(i)=  74.53;el(i)= 760;lnm(i)="FRUNZE                   ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="UTTT";id(i)=38457;lt(i)= 41.27;ln(i)=  69.27;el(i)= 489;lnm(i)="TASKENT/YUZNI            ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=38507;lt(i)= 40.03;ln(i)=  52.98;el(i)=  89;lnm(i)="KRASNOVODSK              ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="OSDI";id(i)=40080;lt(i)= 33.42;ln(i)=  36.52;el(i)= 609;lnm(i)="DAMASCUS (CIV/MIL)       ";st(i)="  ";ct(i)="SY"
i=i+1;cd(i)="OLBA";id(i)=40100;lt(i)= 33.82;ln(i)=  35.48;el(i)=  19;lnm(i)="BEIRUT (CIV/MIL)         ";st(i)="  ";ct(i)="LB"
i=i+1;cd(i)="9999";id(i)=40179;lt(i)= 32.00;ln(i)=  34.82;el(i)=  30;lnm(i)="BET DAGAN                ";st(i)="  ";ct(i)="IL"
i=i+1;cd(i)="9999";id(i)=40186;lt(i)= 31.87;ln(i)=  34.68;el(i)=  20;lnm(i)="ASHDOD NORTH             ";st(i)="  ";ct(i)="IL"
i=i+1;cd(i)="OJMF";id(i)=40265;lt(i)= 32.37;ln(i)=  36.27;el(i)= 687;lnm(i)="MAFRAQ (JOR-AFB)         ";st(i)="  ";ct(i)="JO"
i=i+1;cd(i)="OEPA";id(i)=40373;lt(i)= 28.33;ln(i)=  46.12;el(i)= 355;lnm(i)="HAFR AL-BATIN ARPT       ";st(i)="  ";ct(i)="SA"
i=i+1;cd(i)="OETB";id(i)=40375;lt(i)= 28.37;ln(i)=  36.58;el(i)= 770;lnm(i)="TABUK (SAUD-AFB)         ";st(i)="  ";ct(i)="SA"
i=i+1;cd(i)="OEHL";id(i)=40394;lt(i)= 27.52;ln(i)=  41.73;el(i)=1013;lnm(i)="HAIL                     ";st(i)="  ";ct(i)="SA"
i=i+1;cd(i)="OEMF";id(i)=40417;lt(i)= 26.45;ln(i)=  49.82;el(i)=  12;lnm(i)="KING FAHAD INTERNATIONAL ";st(i)="AP";ct(i)="SA"
i=i+1;cd(i)="OEMA";id(i)=40430;lt(i)= 24.55;ln(i)=  39.72;el(i)= 636;lnm(i)="MADINAH INTL ARPT        ";st(i)="  ";ct(i)="SA"
i=i+1;cd(i)="OERK";id(i)=40437;lt(i)= 24.93;ln(i)=  46.72;el(i)= 612;lnm(i)="RIYADH/KING KHALID       ";st(i)="  ";ct(i)="SA"
i=i+1;cd(i)="OKBK";id(i)=40582;lt(i)= 29.22;ln(i)=  47.98;el(i)=  55;lnm(i)="KUWAIT INTERNATIONAL AIRP";st(i)="OR";ct(i)="KW"
i=i+1;cd(i)="ORBB";id(i)=40650;lt(i)= 33.23;ln(i)=  44.23;el(i)=  34;lnm(i)="BAGHDAD/SIRSENK/BAM      ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="OITT";id(i)=40706;lt(i)= 38.08;ln(i)=  46.28;el(i)=1361;lnm(i)="TABRIZ(IRAN-AB/CIV)      ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="OIMM";id(i)=40745;lt(i)= 36.27;ln(i)=  59.63;el(i)= 980;lnm(i)="MASHHAD (AFB/CIV)        ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="OIII";id(i)=40754;lt(i)= 35.68;ln(i)=  51.35;el(i)=1191;lnm(i)="TEHRAN/MEHRABAD AFB      ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="OICC";id(i)=40766;lt(i)= 34.27;ln(i)=  47.12;el(i)=1322;lnm(i)="KERMANSHAH/BAKTARAN      ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="OIFM";id(i)=40800;lt(i)= 32.62;ln(i)=  51.67;el(i)=1590;lnm(i)="ESFAHAN (CIV/AFB)        ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="OIMB";id(i)=40809;lt(i)= 32.87;ln(i)=  59.20;el(i)=1491;lnm(i)="BIRJAND                  ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="OIMB";id(i)=40811;lt(i)= 31.33;ln(i)=  48.66;el(i)=  22;lnm(i)="BIRJAND                  ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="9999";id(i)=40821;lt(i)= 31.90;ln(i)=  54.28;el(i)=1237;lnm(i)="YAZD                     ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="OIKK";id(i)=40841;lt(i)= 30.25;ln(i)=  56.97;el(i)=1754;lnm(i)="KERMAN                   ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="OISS";id(i)=40848;lt(i)= 29.53;ln(i)=  52.58;el(i)=1491;lnm(i)="SHIRAZ (CIV/AFB)         ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="9999";id(i)=40856;lt(i)= 29.47;ln(i)=  60.88;el(i)=1379;lnm(i)="ZAHEDAN                  ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="OIKB";id(i)=40875;lt(i)= 27.22;ln(i)=  56.37;el(i)=  10;lnm(i)="BANDAR ABBAS INTL        ";st(i)="  ";ct(i)="IR"
i=i+1;cd(i)="OAMS";id(i)=40911;lt(i)= 36.70;ln(i)=  67.20;el(i)= 378;lnm(i)="MAZAR-I-SHARIF           ";st(i)="  ";ct(i)="AH"
i=i+1;cd(i)="OAUZ";id(i)=40913;lt(i)= 36.67;ln(i)=  68.92;el(i)= 433;lnm(i)="KUNDUZ                   ";st(i)="  ";ct(i)="AH"
i=i+1;cd(i)="OAKB";id(i)=40948;lt(i)= 34.55;ln(i)=  69.22;el(i)=1791;lnm(i)="KABUL INTL (MIL)         ";st(i)="  ";ct(i)="AH"
i=i+1;cd(i)="9999";id(i)=40980;lt(i)= 32.83;ln(i)=  67.78;el(i)=2000;lnm(i)="MOKUR                    ";st(i)="  ";ct(i)="AH"
i=i+1;cd(i)="OAKN";id(i)=40990;lt(i)= 31.50;ln(i)=  65.85;el(i)=1010;lnm(i)="KANDAHAR                 ";st(i)="  ";ct(i)="AH"
i=i+1;cd(i)="OEJN";id(i)=41024;lt(i)= 21.67;ln(i)=  39.15;el(i)=  17;lnm(i)="JEDDAH/ABDUL AZIZ        ";st(i)="  ";ct(i)="SA"
i=i+1;cd(i)="OEAB";id(i)=41112;lt(i)= 18.23;ln(i)=  42.65;el(i)=2084;lnm(i)="ABHA                     ";st(i)="  ";ct(i)="SA"
i=i+1;cd(i)="9999";id(i)=41170;lt(i)= 25.27;ln(i)=  51.55;el(i)=  10;lnm(i)="DOHA                     ";st(i)="  ";ct(i)="QT"
i=i+1;cd(i)="OMAA";id(i)=41217;lt(i)= 24.43;ln(i)=  54.65;el(i)=  27;lnm(i)="ABU DHABI INTL           ";st(i)="  ";ct(i)="AE"
i=i+1;cd(i)="OOMS";id(i)=41256;lt(i)= 23.58;ln(i)=  58.28;el(i)=  17;lnm(i)="SEEB INTL/MUSCAT         ";st(i)="  ";ct(i)="OM"
i=i+1;cd(i)="OOSA";id(i)=41316;lt(i)= 17.03;ln(i)=  54.08;el(i)=  17;lnm(i)="SALALAH                  ";st(i)="  ";ct(i)="OM"
i=i+1;cd(i)="9999";id(i)=41480;lt(i)= 12.83;ln(i)=  45.03;el(i)=   4;lnm(i)="ADEN/KHORMAKSAR          ";st(i)="  ";ct(i)="DY"
i=i+1;cd(i)="9999";id(i)=41506;lt(i)= 35.85;ln(i)=  71.83;el(i)=1494;lnm(i)="CHITRAL                  ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41515;lt(i)= 35.40;ln(i)=  71.78;el(i)=1480;lnm(i)="DROSH                    ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="OPPS";id(i)=41530;lt(i)= 34.02;ln(i)=  71.58;el(i)= 360;lnm(i)="PESHAWAR (CIV/MIL)       ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="OPSR";id(i)=41594;lt(i)= 32.05;ln(i)=  72.67;el(i)= 188;lnm(i)="SARGODHA (PAK-AFB)       ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41598;lt(i)= 32.93;ln(i)=  73.72;el(i)= 233;lnm(i)="JHELUM                   ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41624;lt(i)= 31.82;ln(i)=  70.92;el(i)= 174;lnm(i)="DERA ISMAIL KHAN         ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="OPLH";id(i)=41640;lt(i)= 31.55;ln(i)=  74.33;el(i)= 215;lnm(i)="LAHORE/WALTON            ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41661;lt(i)= 30.27;ln(i)=  66.92;el(i)=1621;lnm(i)="QUETTA                   ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41675;lt(i)= 30.20;ln(i)=  71.43;el(i)= 123;lnm(i)="MULTAN                   ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41712;lt(i)= 28.88;ln(i)=  64.40;el(i)= 849;lnm(i)="DALBANDIN                ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41715;lt(i)= 28.30;ln(i)=  68.47;el(i)=  56;lnm(i)="JACOBABAD                ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41718;lt(i)= 28.65;ln(i)=  70.68;el(i)=  87;lnm(i)="KHANPUR                  ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41739;lt(i)= 26.97;ln(i)=  64.10;el(i)= 981;lnm(i)="PANJGUR                  ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41744;lt(i)= 27.83;ln(i)=  66.63;el(i)=1225;lnm(i)="KHUZDAR                  ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41749;lt(i)= 26.25;ln(i)=  68.37;el(i)=  37;lnm(i)="NAWABSHAH                ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41756;lt(i)= 25.07;ln(i)=  61.80;el(i)=  56;lnm(i)="JIWANI                   ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41768;lt(i)= 25.52;ln(i)=  69.78;el(i)=   6;lnm(i)="CHHOR                    ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="OPKC";id(i)=41780;lt(i)= 24.90;ln(i)=  67.13;el(i)=  22;lnm(i)="KARACHI INTL ARPT        ";st(i)="  ";ct(i)="PK"
i=i+1;cd(i)="9999";id(i)=41859;lt(i)= 25.73;ln(i)=  89.23;el(i)=  34;lnm(i)="RANGPUR                  ";st(i)="  ";ct(i)="BW"
i=i+1;cd(i)="VGBG";id(i)=41883;lt(i)= 24.85;ln(i)=  89.37;el(i)=  20;lnm(i)="BOGRA                    ";st(i)="  ";ct(i)="BD"
i=i+1;cd(i)="9999";id(i)=41891;lt(i)= 24.90;ln(i)=  91.88;el(i)=  35;lnm(i)="SYLHET                   ";st(i)="  ";ct(i)="BW"
i=i+1;cd(i)="9999";id(i)=41907;lt(i)= 24.13;ln(i)=  89.05;el(i)=  14;lnm(i)="ISHURDI                  ";st(i)="  ";ct(i)="BW"
i=i+1;cd(i)="9999";id(i)=41923;lt(i)= 23.78;ln(i)=  91.18;el(i)=   9;lnm(i)="COMILLA                  ";st(i)="  ";ct(i)="BD"
i=i+1;cd(i)="9999";id(i)=41936;lt(i)= 23.18;ln(i)=  89.17;el(i)=   7;lnm(i)="JESSORE                  ";st(i)="  ";ct(i)="BW"
i=i+1;cd(i)="9999";id(i)=41943;lt(i)= 23.03;ln(i)=  91.42;el(i)=   8;lnm(i)="FENI                     ";st(i)="  ";ct(i)="BW"
i=i+1;cd(i)="9999";id(i)=41950;lt(i)= 22.75;ln(i)=  90.37;el(i)=   4;lnm(i)="BARISAL                  ";st(i)="  ";ct(i)="BW"
i=i+1;cd(i)="9999";id(i)=41977;lt(i)= 22.35;ln(i)=  91.82;el(i)=  34;lnm(i)="CHITTAGONG/AMBAGAN       ";st(i)="  ";ct(i)="BD"
i=i+1;cd(i)="9999";id(i)=41992;lt(i)= 21.43;ln(i)=  91.93;el(i)=   4;lnm(i)="COX'S BAZAR              ";st(i)="  ";ct(i)="BW"
i=i+1;cd(i)="9999";id(i)=42027;lt(i)= 34.08;ln(i)=  74.83;el(i)=1587;lnm(i)="SRINAGAR                 ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=42101;lt(i)= 30.33;ln(i)=  76.47;el(i)= 251;lnm(i)="PATIALA                  ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VIDD";id(i)=42182;lt(i)= 28.58;ln(i)=  77.20;el(i)= 216;lnm(i)="DELHI/SAFDARJUNG         ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VEMN";id(i)=42314;lt(i)= 27.48;ln(i)=  95.02;el(i)= 111;lnm(i)="DIBRUGARH/MOHANBARI      ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VIJO";id(i)=42339;lt(i)= 26.30;ln(i)=  73.02;el(i)= 224;lnm(i)="JODHPUR (IN-AFB)         ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VIGR";id(i)=42361;lt(i)= 26.23;ln(i)=  78.25;el(i)= 207;lnm(i)="GWALIOR (IN-AFB)         ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VILK";id(i)=42369;lt(i)= 26.75;ln(i)=  80.88;el(i)= 128;lnm(i)="LUCKNOW/AMAUSI           ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VEGK";id(i)=42379;lt(i)= 26.75;ln(i)=  83.37;el(i)=  77;lnm(i)="GORAKHPUR (IN-AFB)       ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=42397;lt(i)= 26.67;ln(i)=  88.37;el(i)= 123;lnm(i)="SILIGURI                 ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VEGT";id(i)=42410;lt(i)= 26.10;ln(i)=  91.58;el(i)=  54;lnm(i)="GAUHATI (IN-AFB)         ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VEPT";id(i)=42492;lt(i)= 25.60;ln(i)=  85.10;el(i)=  60;lnm(i)="PATNA                    ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=42591;lt(i)= 24.75;ln(i)=  84.95;el(i)= 116;lnm(i)="GAVA                     ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=42623;lt(i)= 24.67;ln(i)=  93.90;el(i)= 774;lnm(i)="IMPHAL                   ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VAAH";id(i)=42647;lt(i)= 23.07;ln(i)=  72.63;el(i)=  55;lnm(i)="AHMADABAD                ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VABP";id(i)=42667;lt(i)= 23.28;ln(i)=  77.35;el(i)= 523;lnm(i)="BHOPAL/BAIRAGARH         ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VERC";id(i)=42701;lt(i)= 23.32;ln(i)=  85.32;el(i)= 652;lnm(i)="RANCHI                   ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VEAT";id(i)=42724;lt(i)= 23.88;ln(i)=  91.25;el(i)=  16;lnm(i)="AGARTALA                 ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=42798;lt(i)= 22.82;ln(i)=  86.18;el(i)= 142;lnm(i)="JAMSHEDPUR               ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VECC";id(i)=42809;lt(i)= 22.65;ln(i)=  88.45;el(i)=   6;lnm(i)="CALCUTTA/DUM DUM         ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VANP";id(i)=42867;lt(i)= 21.10;ln(i)=  79.05;el(i)= 310;lnm(i)="NAGPUR SONEGAON AFB      ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=42874;lt(i)= 21.22;ln(i)=  81.67;el(i)= 298;lnm(i)="PBO RAIPUR               ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=42895;lt(i)= 21.50;ln(i)=  86.93;el(i)=  20;lnm(i)="BALASORE                 ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VEBS";id(i)=42971;lt(i)= 20.25;ln(i)=  85.83;el(i)=  46;lnm(i)="BHUBANESWAR              ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VABB";id(i)=43003;lt(i)= 19.12;ln(i)=  72.85;el(i)=  14;lnm(i)="BOMBAY/SANTACRUZ         ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VAAU";id(i)=43014;lt(i)= 19.85;ln(i)=  75.40;el(i)= 579;lnm(i)="AURANGABAD AIRPORT       ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=43041;lt(i)= 19.08;ln(i)=  82.03;el(i)= 553;lnm(i)="JAGDALPUR                ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VOHY";id(i)=43128;lt(i)= 17.45;ln(i)=  78.47;el(i)= 545;lnm(i)="HYDERABAD (CIV/MIL)      ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=43150;lt(i)= 17.70;ln(i)=  83.30;el(i)=  66;lnm(i)="VISHAKHAPATNAM/WALT      ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=43185;lt(i)= 16.20;ln(i)=  81.15;el(i)=   3;lnm(i)="MACHILIPATNAM            ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=43192;lt(i)= 15.48;ln(i)=  73.82;el(i)=  60;lnm(i)="GOA/PANJIM               ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VOMM";id(i)=43279;lt(i)= 13.00;ln(i)=  80.18;el(i)=  16;lnm(i)="MADRAS/MINAMBAKKAM       ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=43285;lt(i)= 12.95;ln(i)=  74.83;el(i)=  31;lnm(i)="PANAMBUR                 ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=43295;lt(i)= 12.97;ln(i)=  77.58;el(i)= 921;lnm(i)="BANGALORE                ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=43311;lt(i)= 11.12;ln(i)=  72.73;el(i)=   4;lnm(i)="AMINI DIVI               ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VEPB";id(i)=43333;lt(i)= 11.67;ln(i)=  92.72;el(i)=  79;lnm(i)="PORT BLAIR               ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=43346;lt(i)= 10.92;ln(i)=  79.83;el(i)=   7;lnm(i)="KARAIKAL                 ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="VOCC";id(i)=43353;lt(i)=  9.95;ln(i)=  76.27;el(i)=   3;lnm(i)="COCHIN (IN-NAVY)         ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=43369;lt(i)=  8.30;ln(i)=  73.15;el(i)=   2;lnm(i)="MINICOY ISLAND           ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=43371;lt(i)=  8.48;ln(i)=  76.95;el(i)=  64;lnm(i)="THIRUVANANTHAPURAM       ";st(i)="  ";ct(i)="IN"
i=i+1;cd(i)="9999";id(i)=43418;lt(i)=  8.58;ln(i)=  81.25;el(i)=   7;lnm(i)="TRINCOMALEE              ";st(i)="  ";ct(i)="SB"
i=i+1;cd(i)="9999";id(i)=43424;lt(i)=  8.03;ln(i)=  79.83;el(i)=   2;lnm(i)="PUTTALAM                 ";st(i)="  ";ct(i)="SB"
i=i+1;cd(i)="9999";id(i)=43466;lt(i)=  6.90;ln(i)=  79.87;el(i)=   7;lnm(i)="COLOMBO                  ";st(i)="  ";ct(i)="LK"
i=i+1;cd(i)="9999";id(i)=43497;lt(i)=  6.12;ln(i)=  81.13;el(i)=  20;lnm(i)="HANBANTOTA               ";st(i)="  ";ct(i)="SB"
i=i+1;cd(i)="9999";id(i)=43555;lt(i)=  4.20;ln(i)=  73.53;el(i)=   2;lnm(i)="MALE                     ";st(i)="  ";ct(i)="MV"
i=i+1;cd(i)="9999";id(i)=43599;lt(i)= 00.68;ln(i)=  73.15;el(i)=   2;lnm(i)="GAN                      ";st(i)="  ";ct(i)="MV"
i=i+1;cd(i)="9999";id(i)=44212;lt(i)= 49.97;ln(i)=  92.08;el(i)= 936;lnm(i)="ULAN-GOM                 ";st(i)="  ";ct(i)="MN"
i=i+1;cd(i)="9999";id(i)=44231;lt(i)= 49.63;ln(i)= 100.17;el(i)=1288;lnm(i)="MUREN                    ";st(i)="  ";ct(i)="MN"
i=i+1;cd(i)="9999";id(i)=44259;lt(i)= 48.07;ln(i)= 114.50;el(i)= 756;lnm(i)="CHOIBALSAN               ";st(i)="  ";ct(i)="MN"
i=i+1;cd(i)="9999";id(i)=44277;lt(i)= 46.40;ln(i)=  96.25;el(i)=2147;lnm(i)="ALTAI                    ";st(i)="  ";ct(i)="MN"
i=i+1;cd(i)="9999";id(i)=44288;lt(i)= 46.27;ln(i)= 102.78;el(i)=1832;lnm(i)="ARBAIHER                 ";st(i)="  ";ct(i)="MN"
i=i+1;cd(i)="9999";id(i)=44292;lt(i)= 47.93;ln(i)= 106.98;el(i)=1313;lnm(i)="ULAN-BATOR               ";st(i)="  ";ct(i)="MN"
i=i+1;cd(i)="9999";id(i)=44373;lt(i)= 43.58;ln(i)= 104.42;el(i)=1470;lnm(i)="DALANZADGAD              ";st(i)="  ";ct(i)="MN"
i=i+1;cd(i)="9999";id(i)=45004;lt(i)= 22.32;ln(i)= 114.17;el(i)=  66;lnm(i)="KING'S PARK              ";st(i)="  ";ct(i)="HK"
i=i+1;cd(i)="9999";id(i)=46750;lt(i)= 22.67;ln(i)= 120.45;el(i)=  24;lnm(i)="PINGTUNG SOUTH           ";st(i)="  ";ct(i)="TW"
i=i+1;cd(i)="9999";id(i)=46780;lt(i)= 22.68;ln(i)= 121.50;el(i)= 280;lnm(i)="LU-TAU                   ";st(i)="  ";ct(i)="TW"
i=i+1;cd(i)="9999";id(i)=47058;lt(i)= 39.03;ln(i)= 125.78;el(i)=  38;lnm(i)="PYONGYANG                ";st(i)="  ";ct(i)="KO"
i=i+1;cd(i)="9999";id(i)=47090;lt(i)= 38.25;ln(i)= 128.57;el(i)=  18;lnm(i)="SOKCHO                   ";st(i)="  ";ct(i)="KO"
i=i+1;cd(i)="9999";id(i)=47102;lt(i)= 38.03;ln(i)= 127.15;el(i)=  70;lnm(i)="CHONGONG-NI              ";st(i)="  ";ct(i)="KO"
i=i+1;cd(i)="RKSO";id(i)=47122;lt(i)= 37.10;ln(i)= 127.03;el(i)=  12;lnm(i)="OSAN (US/KOR-AFB)        ";st(i)="  ";ct(i)="KO"
i=i+1;cd(i)="9999";id(i)=47138;lt(i)= 36.03;ln(i)= 129.38;el(i)=   6;lnm(i)="POHANG                   ";st(i)="  ";ct(i)="KO"
i=i+1;cd(i)="RKJJ";id(i)=47158;lt(i)= 35.12;ln(i)= 126.82;el(i)=  13;lnm(i)="KWANGJU (KOR-AFB)        ";st(i)="  ";ct(i)="KO"
i=i+1;cd(i)="9999";id(i)=47169;lt(i)= 34.68;ln(i)= 125.45;el(i)=  83;lnm(i)="HEUKSANDO                ";st(i)="  ";ct(i)="KO"
i=i+1;cd(i)="9999";id(i)=47185;lt(i)= 33.28;ln(i)= 126.17;el(i)=  73;lnm(i)="CHEJU                    ";st(i)="  ";ct(i)="KO"
i=i+1;cd(i)="9999";id(i)=47261;lt(i)= 34.55;ln(i)= 126.58;el(i)=  14;lnm(i)="HAEMAN                   ";st(i)="  ";ct(i)="KO"
i=i+1;cd(i)="9999";id(i)=47401;lt(i)= 45.42;ln(i)= 141.68;el(i)=  11;lnm(i)="WAKKANAI                 ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47412;lt(i)= 43.05;ln(i)= 141.33;el(i)=  19;lnm(i)="SAPPORO                  ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="RJCS";id(i)=47418;lt(i)= 42.97;ln(i)= 144.40;el(i)=  37;lnm(i)="KUSHIRO                  ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47420;lt(i)= 43.33;ln(i)= 145.58;el(i)=  26;lnm(i)="NEMURO                   ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="RJSM";id(i)=47580;lt(i)= 40.68;ln(i)= 141.38;el(i)=  39;lnm(i)="MISAWA                   ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47582;lt(i)= 39.72;ln(i)= 140.10;el(i)=  21;lnm(i)="AKITA                    ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47590;lt(i)= 38.27;ln(i)= 140.90;el(i)=  43;lnm(i)="SENDAI                   ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47600;lt(i)= 37.38;ln(i)= 136.90;el(i)=  14;lnm(i)="WAJIMA                   ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47646;lt(i)= 36.05;ln(i)= 140.13;el(i)=  31;lnm(i)="TATENO                   ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47678;lt(i)= 33.12;ln(i)= 139.78;el(i)=  80;lnm(i)="HACHIJOJIMA/OMURE        ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="RJNH";id(i)=47681;lt(i)= 34.73;ln(i)= 137.67;el(i)=  48;lnm(i)="HAMAMATSU AB             ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47741;lt(i)= 35.46;ln(i)= 133.06;el(i)=  22;lnm(i)="MATSUE                   ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47744;lt(i)= 35.43;ln(i)= 133.35;el(i)=   8;lnm(i)="YONAGO                   ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47778;lt(i)= 33.45;ln(i)= 135.77;el(i)=  75;lnm(i)="SHIONOMISAKI             ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47807;lt(i)= 33.58;ln(i)= 130.38;el(i)=  14;lnm(i)="FUKUOKA                  ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47827;lt(i)= 31.63;ln(i)= 130.58;el(i)=  31;lnm(i)="KAGOSHIMA/YOSHINO        ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="RJOS";id(i)=47881;lt(i)= 34.13;ln(i)= 134.62;el(i)=  11;lnm(i)="TOKUSHIMA(JMSDF/CV)      ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47909;lt(i)= 28.38;ln(i)= 129.55;el(i)= 295;lnm(i)="NAZE/FUNCHATOGE          ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="ROIG";id(i)=47918;lt(i)= 24.33;ln(i)= 124.17;el(i)=   7;lnm(i)="ISHIGAKIJIMA ISLAND      ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=47936;lt(i)= 26.20;ln(i)= 127.68;el(i)=  27;lnm(i)="NAHA AIRPORT             ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="ROMD";id(i)=47945;lt(i)= 25.83;ln(i)= 131.23;el(i)=  20;lnm(i)="MINAMIDAITOJIMA ISLAND   ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="RJAO";id(i)=47971;lt(i)= 27.08;ln(i)= 142.18;el(i)=   8;lnm(i)="CHICHIJIMA ISLAND        ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="RJAW";id(i)=47981;lt(i)= 24.78;ln(i)= 141.33;el(i)= 116;lnm(i)="IWOJIMA (JMSDF)          ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="RJAM";id(i)=47991;lt(i)= 24.30;ln(i)= 153.97;el(i)=   9;lnm(i)="MINAMITORISHIMA          ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="VBRM";id(i)=48042;lt(i)= 21.98;ln(i)=  96.10;el(i)=  76;lnm(i)="MANDALAY                 ";st(i)="  ";ct(i)="MM"
i=i+1;cd(i)="9999";id(i)=48097;lt(i)= 16.77;ln(i)=  96.17;el(i)=  15;lnm(i)="RANGOON                  ";st(i)="  ";ct(i)="BM"
i=i+1;cd(i)="VTCC";id(i)=48327;lt(i)= 18.78;ln(i)=  98.98;el(i)= 314;lnm(i)="CHIANG MAI (CIV/AFB)     ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="9999";id(i)=48354;lt(i)= 17.38;ln(i)= 102.80;el(i)= 178;lnm(i)="UDON THANI               ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="9999";id(i)=48378;lt(i)= 16.82;ln(i)= 100.27;el(i)=  44;lnm(i)="PHITSANULOK              ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="VTUU";id(i)=48407;lt(i)= 15.25;ln(i)= 104.87;el(i)= 127;lnm(i)="UBON/RATCHATHANI         ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="9999";id(i)=48431;lt(i)= 14.97;ln(i)= 102.08;el(i)= 188;lnm(i)="NAKHON RATCHASIMA        ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="9999";id(i)=48453;lt(i)= 13.67;ln(i)= 100.60;el(i)=   3;lnm(i)="BANGNA AGROMET           ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="9999";id(i)=48455;lt(i)= 13.73;ln(i)= 100.50;el(i)=  20;lnm(i)="BANGKOK                  ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="9999";id(i)=48477;lt(i)= 12.68;ln(i)= 100.98;el(i)=  16;lnm(i)="SATTAHIP                 ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="9999";id(i)=48480;lt(i)= 12.60;ln(i)= 102.12;el(i)=   5;lnm(i)="CHANTHABURI              ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="9999";id(i)=48500;lt(i)= 11.80;ln(i)=  99.80;el(i)=   5;lnm(i)="PRACHUAP KHIRIKHAN       ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="9999";id(i)=48551;lt(i)=  9.12;ln(i)=  99.35;el(i)=  10;lnm(i)="SURAT THANI              ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="VTSP";id(i)=48565;lt(i)=  8.12;ln(i)=  98.32;el(i)=  10;lnm(i)="PHUKET INTL AIRPORT      ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="VTSH";id(i)=48568;lt(i)=  7.20;ln(i)= 100.60;el(i)=   5;lnm(i)="SONGKHLA (THAI-NAVY)     ";st(i)="  ";ct(i)="TH"
i=i+1;cd(i)="WMKP";id(i)=48601;lt(i)=  5.30;ln(i)= 100.27;el(i)=   4;lnm(i)="PENANG/BAYAN LEPAS       ";st(i)="  ";ct(i)="MY"
i=i+1;cd(i)="WMKC";id(i)=48615;lt(i)=  6.17;ln(i)= 102.28;el(i)=   5;lnm(i)="KOTA BHARU/SULTAN P      ";st(i)="  ";ct(i)="MY"
i=i+1;cd(i)="9999";id(i)=48620;lt(i)=  4.22;ln(i)= 100.70;el(i)=   8;lnm(i)="SITIAWAN                 ";st(i)="  ";ct(i)="MS"
i=i+1;cd(i)="9999";id(i)=48650;lt(i)=  2.73;ln(i)= 101.70;el(i)=  17;lnm(i)="SEPANG                   ";st(i)="  ";ct(i)="MY"
i=i+1;cd(i)="WMKD";id(i)=48657;lt(i)=  3.78;ln(i)= 103.22;el(i)=  16;lnm(i)="KUANTAN (AFB)            ";st(i)="  ";ct(i)="MY"
i=i+1;cd(i)="WSSS";id(i)=48698;lt(i)=  1.37;ln(i)= 103.98;el(i)=  16;lnm(i)="SINGAPORE/CHANGI         ";st(i)="  ";ct(i)="SG"
i=i+1;cd(i)="VVNB";id(i)=48811;lt(i)= 21.40;ln(i)= 113.02;el(i)= 472;lnm(i)="DIEN BIEN PHU            ";st(i)="  ";ct(i)="VN"
i=i+1;cd(i)="VVNB";id(i)=48820;lt(i)= 21.02;ln(i)= 105.80;el(i)=   6;lnm(i)="HANOI/NOIBAI INTL        ";st(i)="  ";ct(i)="VN"
i=i+1;cd(i)="9999";id(i)=48839;lt(i)= 20.13;ln(i)= 107.72;el(i)=  56;lnm(i)="BACH LONGVI HP           ";st(i)="  ";ct(i)="VN"
i=i+1;cd(i)="VVVH";id(i)=48845;lt(i)= 18.70;ln(i)= 105.66;el(i)=   6;lnm(i)="VINH                     ";st(i)="  ";ct(i)="VS"
i=i+1;cd(i)="VVDN";id(i)=48855;lt(i)= 16.03;ln(i)= 108.18;el(i)=   7;lnm(i)="DA NANG INTL ARPT        ";st(i)="  ";ct(i)="VN"
i=i+1;cd(i)="9999";id(i)=48870;lt(i)= 13.77;ln(i)= 109.22;el(i)=   5;lnm(i)="QUI NHON                 ";st(i)="  ";ct(i)="VS"
i=i+1;cd(i)="9999";id(i)=48887;lt(i)= 10.93;ln(i)= 108.10;el(i)=   5;lnm(i)="PHAN THIET               ";st(i)="  ";ct(i)="VS"
i=i+1;cd(i)="VVTS";id(i)=48900;lt(i)= 10.82;ln(i)= 106.67;el(i)=  19;lnm(i)="HO CHI MINH/TANSONN      ";st(i)="  ";ct(i)="VN"
i=i+1;cd(i)="9999";id(i)=48914;lt(i)=  9.18;ln(i)= 105.15;el(i)=   2;lnm(i)="CA MAU                   ";st(i)="  ";ct(i)="VS"
i=i+1;cd(i)="9999";id(i)=50527;lt(i)= 49.22;ln(i)= 119.75;el(i)= 611;lnm(i)="HAILAR                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=50557;lt(i)= 49.17;ln(i)= 125.23;el(i)= 243;lnm(i)="NENJIANG                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=50774;lt(i)= 47.72;ln(i)= 128.90;el(i)= 232;lnm(i)="YICHUN                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=50953;lt(i)= 45.75;ln(i)= 126.77;el(i)= 143;lnm(i)="HARBIN                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=51076;lt(i)= 47.73;ln(i)=  88.08;el(i)= 737;lnm(i)="ALTAY                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=51156;lt(i)= 46.78;ln(i)=  85.72;el(i)=1294;lnm(i)="HOBOG SAIR               ";st(i)="  ";ct(i)="CI"
i=i+1;cd(i)="ZWYN";id(i)=51431;lt(i)= 43.95;ln(i)=  81.33;el(i)= 663;lnm(i)="YINING                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=51463;lt(i)= 43.78;ln(i)=  87.62;el(i)= 919;lnm(i)="URUMQI                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=51644;lt(i)= 41.72;ln(i)=  82.95;el(i)=1100;lnm(i)="KUQA                     ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZWSH";id(i)=51709;lt(i)= 39.47;ln(i)=  75.98;el(i)=1291;lnm(i)="KASHI                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=51777;lt(i)= 39.03;ln(i)=  88.17;el(i)= 889;lnm(i)="RUOQIANG                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZWTN";id(i)=51828;lt(i)= 37.13;ln(i)=  79.93;el(i)=1375;lnm(i)="HOTAN                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=51839;lt(i)= 37.07;ln(i)=  82.77;el(i)=1410;lnm(i)="MINFENG/NIYA             ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZWHM";id(i)=52203;lt(i)= 42.82;ln(i)=  93.52;el(i)= 739;lnm(i)="HAMI                     ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=52267;lt(i)= 41.95;ln(i)= 101.07;el(i)= 941;lnm(i)="EJIN QI                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=52323;lt(i)= 41.80;ln(i)=  97.03;el(i)=1770;lnm(i)="MAZONG SHAN (MOUNT)      ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=52418;lt(i)= 40.15;ln(i)=  94.68;el(i)=1140;lnm(i)="DUNHUANG                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZLJQ";id(i)=52533;lt(i)= 39.77;ln(i)=  98.48;el(i)=1478;lnm(i)="JIUQUAN/SUZHOU           ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=52681;lt(i)= 38.63;ln(i)= 103.08;el(i)=1367;lnm(i)="MINQIN                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=52818;lt(i)= 36.42;ln(i)=  94.90;el(i)=2809;lnm(i)="GOLMUD                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=52836;lt(i)= 36.30;ln(i)=  98.10;el(i)=3192;lnm(i)="DULAN/QAGAN US           ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZLXN";id(i)=52866;lt(i)= 36.62;ln(i)= 101.77;el(i)=2262;lnm(i)="XINING                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=52983;lt(i)= 35.87;ln(i)= 104.15;el(i)=1875;lnm(i)="YU ZHONG                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=53068;lt(i)= 43.65;ln(i)= 112.00;el(i)= 966;lnm(i)="ERENHOT                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZBHH";id(i)=53463;lt(i)= 40.82;ln(i)= 111.68;el(i)=1065;lnm(i)="HOHHOT                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=53513;lt(i)= 40.77;ln(i)= 107.40;el(i)=1041;lnm(i)="LINHE                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZLIC";id(i)=53614;lt(i)= 38.48;ln(i)= 106.22;el(i)=1112;lnm(i)="YINCHUAN                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZBYN";id(i)=53772;lt(i)= 37.78;ln(i)= 112.55;el(i)= 779;lnm(i)="TAIYUAN/WUSU             ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZLYA";id(i)=53845;lt(i)= 36.60;ln(i)= 109.50;el(i)= 959;lnm(i)="YAN AN                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=53915;lt(i)= 35.55;ln(i)= 106.67;el(i)=1348;lnm(i)="PINGLIANG                ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=54102;lt(i)= 43.95;ln(i)= 116.07;el(i)= 991;lnm(i)="XILIN HOT/ABAGNAR        ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=54135;lt(i)= 43.60;ln(i)= 122.27;el(i)= 180;lnm(i)="TONGLIAO                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZYCC";id(i)=54161;lt(i)= 43.90;ln(i)= 125.22;el(i)= 238;lnm(i)="CHANGCHUN                ";st(i)="  ";ct(i)="CN" 
i=i+1;cd(i)="9999";id(i)=54218;lt(i)= 42.27;ln(i)= 118.97;el(i)= 572;lnm(i)="CHIFENG/ULANHAD          ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=54292;lt(i)= 42.88;ln(i)= 129.47;el(i)= 178;lnm(i)="YANJI                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=54337;lt(i)= 41.12;ln(i)= 121.07;el(i)=  70;lnm(i)="JINZHOU                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZYYY";id(i)=54342;lt(i)= 41.77;ln(i)= 123.43;el(i)=  43;lnm(i)="SHENYANG/DONGTA          ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=54374;lt(i)= 41.72;ln(i)= 126.92;el(i)= 333;lnm(i)="LINJIANG                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=54497;lt(i)= 40.08;ln(i)= 124.33;el(i)=  14;lnm(i)="DANDONG                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZBAA";id(i)=54511;lt(i)= 39.93;ln(i)= 116.28;el(i)=  55;lnm(i)="BEIJING/PEKING           ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZYTL";id(i)=54662;lt(i)= 38.90;ln(i)= 121.63;el(i)=  97;lnm(i)="DALIAN/DAIREN/LUDA       ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=54727;lt(i)= 36.70;ln(i)= 117.55;el(i)= 123;lnm(i)="ZHANGQUI                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZSTN";id(i)=54823;lt(i)= 36.68;ln(i)= 116.98;el(i)=  58;lnm(i)="JINAN/TSINAN             ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZSQD";id(i)=54857;lt(i)= 36.07;ln(i)= 120.33;el(i)=  77;lnm(i)="QINGDAO/TSINGTAO         ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=55299;lt(i)= 31.48;ln(i)=  92.07;el(i)=4508;lnm(i)="NAGQU                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZULS";id(i)=55591;lt(i)= 29.67;ln(i)=  91.13;el(i)=3650;lnm(i)="LHASA                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=56029;lt(i)= 33.02;ln(i)=  97.02;el(i)=3682;lnm(i)="YUSHU                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=56080;lt(i)= 34.97;ln(i)= 102.90;el(i)=2910;lnm(i)="HEZUO                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=56137;lt(i)= 31.15;ln(i)=  97.17;el(i)=3307;lnm(i)="QAMDO                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=56146;lt(i)= 31.62;ln(i)= 100.00;el(i)=3394;lnm(i)="GARZE                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=56187;lt(i)= 30.70;ln(i)= 103.83;el(i)= 541;lnm(i)="WENJIANG                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZUUU";id(i)=56294;lt(i)= 30.67;ln(i)= 104.02;el(i)= 508;lnm(i)="CHENGDU                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=56571;lt(i)= 27.90;ln(i)= 102.27;el(i)=1599;lnm(i)="XICHANG                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=56691;lt(i)= 26.87;ln(i)= 104.28;el(i)=2236;lnm(i)="WEINING                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=56739;lt(i)= 25.03;ln(i)=  98.48;el(i)=1649;lnm(i)="TENGCHONG                ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZPPP";id(i)=56778;lt(i)= 25.02;ln(i)= 102.68;el(i)=1892;lnm(i)="KUNMING/WUJIABA          ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=56964;lt(i)= 22.77;ln(i)= 100.98;el(i)=1303;lnm(i)="SIMAO                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=56985;lt(i)= 23.38;ln(i)= 103.38;el(i)=1302;lnm(i)="MENGZI                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZLSN";id(i)=57036;lt(i)= 34.30;ln(i)= 108.93;el(i)= 398;lnm(i)="XI'AN                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=57067;lt(i)= 34.00;ln(i)= 111.02;el(i)= 569;lnm(i)="LU SHIH                  ";st(i)="  ";ct(i)="CI"
i=i+1;cd(i)="ZHCC";id(i)=57083;lt(i)= 34.72;ln(i)= 113.65;el(i)= 111;lnm(i)="ZHENGZHOU                ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=57127;lt(i)= 33.07;ln(i)= 107.03;el(i)= 509;lnm(i)="HANZHONG                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=57131;lt(i)= 34.43;ln(i)= 108.97;el(i)= 387;lnm(i)="JINJHE                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=57178;lt(i)= 33.02;ln(i)= 112.53;el(i)= 131;lnm(i)="NANYANG                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=57447;lt(i)= 30.28;ln(i)= 109.47;el(i)= 458;lnm(i)="ENSHI                    ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=57461;lt(i)= 30.70;ln(i)= 111.28;el(i)= 134;lnm(i)="YICHANG                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZHHH";id(i)=57494;lt(i)= 30.62;ln(i)= 114.13;el(i)=  23;lnm(i)="WUHAN/NANHU              ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZUCK";id(i)=57516;lt(i)= 29.52;ln(i)= 106.48;el(i)= 351;lnm(i)="CHONGQING/CHUNGKING      ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZGCS";id(i)=57679;lt(i)= 28.20;ln(i)= 113.08;el(i)=  46;lnm(i)="CHANGSHA/DATUOPU         ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=57749;lt(i)= 27.57;ln(i)= 110.00;el(i)= 261;lnm(i)="HUAIHUA                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZUGY";id(i)=57816;lt(i)= 26.58;ln(i)= 106.72;el(i)=1074;lnm(i)="GUIYANG                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZGKL";id(i)=57957;lt(i)= 25.33;ln(i)= 110.30;el(i)= 166;lnm(i)="GUILIN                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=57972;lt(i)= 25.82;ln(i)= 113.02;el(i)= 185;lnm(i)="CHENZHOU                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZSGZ";id(i)=57993;lt(i)= 25.85;ln(i)= 114.95;el(i)= 125;lnm(i)="GANZHOU                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=58027;lt(i)= 34.28;ln(i)= 117.15;el(i)=  42;lnm(i)="XUZHOU                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=58150;lt(i)= 33.77;ln(i)= 120.25;el(i)=   7;lnm(i)="SHEYANG/HEDE             ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=58203;lt(i)= 32.93;ln(i)= 115.83;el(i)=  39;lnm(i)="FUYANG                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZSNJ";id(i)=58238;lt(i)= 32.00;ln(i)= 118.80;el(i)=  12;lnm(i)="NANJING/NANKING          ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=58362;lt(i)= 31.40;ln(i)= 121.47;el(i)=   4;lnm(i)="SHANGHAI                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=58424;lt(i)= 30.52;ln(i)= 117.03;el(i)=  20;lnm(i)="ANQING                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZSHC";id(i)=58457;lt(i)= 30.23;ln(i)= 120.17;el(i)=  43;lnm(i)="HANGZHOU/JIANQIAO        ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZSCN";id(i)=58606;lt(i)= 28.60;ln(i)= 115.92;el(i)=  50;lnm(i)="NANCHANG                 ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=58633;lt(i)= 28.97;ln(i)= 118.87;el(i)=  71;lnm(i)="QU XIAN                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=58665;lt(i)= 28.65;ln(i)= 120.08;el(i)=   9;lnm(i)="LUQIAO                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=58725;lt(i)= 27.33;ln(i)= 117.43;el(i)= 192;lnm(i)="SHAOWU                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZSFZ";id(i)=58847;lt(i)= 26.08;ln(i)= 119.28;el(i)=  85;lnm(i)="FUZHOU                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=58968;lt(i)= 25.03;ln(i)= 121.53;el(i)=   9;lnm(i)="TAIPEI                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=59134;lt(i)= 24.45;ln(i)= 118.07;el(i)=  63;lnm(i)="XIA-MEN                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=59211;lt(i)= 23.92;ln(i)= 106.62;el(i)= 242;lnm(i)="BOSE                     ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=59265;lt(i)= 23.48;ln(i)= 111.30;el(i)= 120;lnm(i)="WUZHOU                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=59280;lt(i)= 23.67;ln(i)= 113.05;el(i)=  19;lnm(i)="PING YUAN                ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZGOW";id(i)=59316;lt(i)= 23.40;ln(i)= 116.68;el(i)=   3;lnm(i)="SHANTOU                  ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZGNN";id(i)=59431;lt(i)= 22.82;ln(i)= 108.35;el(i)=  73;lnm(i)="NANNING/WUXU             ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZGBH";id(i)=59644;lt(i)= 21.48;ln(i)= 109.10;el(i)=  16;lnm(i)="FUCHENG/BEIHAI AIRPORT   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=59663;lt(i)= 21.83;ln(i)= 111.97;el(i)=  91;lnm(i)="YANGJIANG                ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="ZGHK";id(i)=59758;lt(i)= 20.03;ln(i)= 110.35;el(i)=  15;lnm(i)="HAIKOU                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=59948;lt(i)= 18.22;ln(i)= 109.58;el(i)= 420;lnm(i)="YULIN/YAXIAN SANYA       ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=59981;lt(i)= 16.83;ln(i)= 112.33;el(i)=   5;lnm(i)="XISHA ISLAND             ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=60010;lt(i)= 28.30;ln(i)= -16.50;el(i)=2368;lnm(i)="IZANA MOUNTAIN TOP       ";st(i)="  ";ct(i)="ES"
i=i+1;cd(i)="9999";id(i)=60018;lt(i)= 28.31;ln(i)= -16.25;el(i)= 120;lnm(i)="VALLE DE GUIMAR TENERIFE ";st(i)="  ";ct(i)="ES"
i=i+1;cd(i)="9999";id(i)=60020;lt(i)= 28.47;ln(i)= -16.25;el(i)=  36;lnm(i)="SANTA CRUZ TENERIFE      ";st(i)="  ";ct(i)="ES"
i=i+1;cd(i)="GMMC";id(i)=60155;lt(i)= 33.57;ln(i)=  -7.67;el(i)=  62;lnm(i)="CASABLANCA/ANFA          ";st(i)="  ";ct(i)="MA"
i=i+1;cd(i)="9999";id(i)=60191;lt(i)= 32.37;ln(i)=  -6.40;el(i)= 468;lnm(i)="BENI-MELLAL              ";st(i)="  ";ct(i)="MA"
i=i+1;cd(i)="GMAD";id(i)=60252;lt(i)= 30.33;ln(i)=  -9.42;el(i)=  74;lnm(i)="AL MASSIRA MC            ";st(i)="  ";ct(i)="MA"
i=i+1;cd(i)="DAAG";id(i)=60390;lt(i)= 36.72;ln(i)=   3.25;el(i)=  25;lnm(i)="DAR-EL-BEIDA/HOUARI      ";st(i)="  ";ct(i)="DZ"
i=i+1;cd(i)="DAOO";id(i)=60490;lt(i)= 35.63;ln(i)=  -0.60;el(i)=  90;lnm(i)="ORAN/ES SENIA            ";st(i)="  ";ct(i)="DZ"
i=i+1;cd(i)="DAUB";id(i)=60525;lt(i)= 34.80;ln(i)=   5.73;el(i)=  87;lnm(i)="BISKRA                   ";st(i)="  ";ct(i)="DZ"
i=i+1;cd(i)="DAAY";id(i)=60549;lt(i)= 34.93;ln(i)=  -0.43;el(i)=1149;lnm(i)="MECHERIA                 ";st(i)="  ";ct(i)="DZ"
i=i+1;cd(i)="DAOR";id(i)=60571;lt(i)= 31.62;ln(i)=  -2.23;el(i)= 773;lnm(i)="BECHAR/OUAKDA            ";st(i)="  ";ct(i)="DZ"
i=i+1;cd(i)="DAUU";id(i)=60580;lt(i)= 31.92;ln(i)=   5.40;el(i)= 141;lnm(i)="OUARGLA                  ";st(i)="  ";ct(i)="DZ"
i=i+1;cd(i)="9999";id(i)=60630;lt(i)= 27.20;ln(i)=   2.47;el(i)= 293;lnm(i)="IN SALAH                 ";st(i)="  ";ct(i)="DZ"
i=i+1;cd(i)="DAOF";id(i)=60656;lt(i)= 27.67;ln(i)=  -8.13;el(i)= 431;lnm(i)="TINDOUF                  ";st(i)="  ";ct(i)="DZ"
i=i+1;cd(i)="9999";id(i)=60686;lt(i)= 21.33;ln(i)=   0.95;el(i)= 399;lnm(i)="BORDJ-BADJ-MOKHTAR       ";st(i)="  ";ct(i)="DZ"
i=i+1;cd(i)="9999";id(i)=60680;lt(i)= 22.78;ln(i)=   5.52;el(i)=1378;lnm(i)="TAMANRASSET              ";st(i)="  ";ct(i)="DZ"
i=i+1;cd(i)="DTTA";id(i)=60715;lt(i)= 36.83;ln(i)=  10.23;el(i)=   4;lnm(i)="TUNIS/CARTHAGE           ";st(i)="  ";ct(i)="TN"
i=i+1;cd(i)="DTTX";id(i)=60750;lt(i)= 34.72;ln(i)=  10.68;el(i)=  23;lnm(i)="SFAX/EL-MAOU             ";st(i)="  ";ct(i)="TN"
i=i+1;cd(i)="DTTZ";id(i)=60760;lt(i)= 33.92;ln(i)=   8.17;el(i)=  93;lnm(i)="TOZEUR/NEFTA             ";st(i)="  ";ct(i)="TN"
i=i+1;cd(i)="9999";id(i)=61017;lt(i)= 18.68;ln(i)=  12.92;el(i)= 357;lnm(i)="BILMA                    ";st(i)="  ";ct(i)="NR"
i=i+1;cd(i)="DRZA";id(i)=61024;lt(i)= 16.97;ln(i)=   7.98;el(i)= 502;lnm(i)="AGADEZ SOUTH (MIL)       ";st(i)="  ";ct(i)="NE"
i=i+1;cd(i)="DRRR";id(i)=61052;lt(i)= 13.48;ln(i)=   2.17;el(i)= 227;lnm(i)="NIAMEY-AERO              ";st(i)="  ";ct(i)="NE"
i=i+1;cd(i)="DRRB";id(i)=61075;lt(i)= 13.80;ln(i)=   5.25;el(i)= 273;lnm(i)="BIRNI-N'KONNI            ";st(i)="  ";ct(i)="NE"
i=i+1;cd(i)="9999";id(i)=61090;lt(i)= 13.78;ln(i)=   8.98;el(i)= 460;lnm(i)="ZINDER                   ";st(i)="  ";ct(i)="NR"
i=i+1;cd(i)="GATS";id(i)=61202;lt(i)= 20.20;ln(i)=   0.98;el(i)= 491;lnm(i)="TESSALIT                 ";st(i)="  ";ct(i)="ML"
i=i+1;cd(i)="9999";id(i)=61214;lt(i)= 18.43;ln(i)=   1.35;el(i)= 458;lnm(i)="KIDAL                    ";st(i)="  ";ct(i)="MI"
i=i+1;cd(i)="GATB";id(i)=61223;lt(i)= 16.72;ln(i)=  -3.00;el(i)= 264;lnm(i)="TOMBOUCTOU/TIMBUKTU      ";st(i)="  ";ct(i)="ML"
i=i+1;cd(i)="9999";id(i)=61226;lt(i)= 16.27;ln(i)=  -0.05;el(i)= 265;lnm(i)="GAO                      ";st(i)="  ";ct(i)="MI"
i=i+1;cd(i)="GANK";id(i)=61233;lt(i)= 15.17;ln(i)=  -7.28;el(i)= 265;lnm(i)="NARA                     ";st(i)="  ";ct(i)="ML"
i=i+1;cd(i)="9999";id(i)=61265;lt(i)= 14.52;ln(i)=  -4.10;el(i)= 276;lnm(i)="MOPTI                    ";st(i)="  ";ct(i)="MI"
i=i+1;cd(i)="GABS";id(i)=61291;lt(i)= 12.53;ln(i)=  -7.95;el(i)= 381;lnm(i)="BAMAKO/SENOU (MIL)       ";st(i)="  ";ct(i)="ML"
i=i+1;cd(i)="9999";id(i)=61404;lt(i)= 22.75;ln(i)= -12.48;el(i)= 343;lnm(i)="ZOUERATE                 ";st(i)="  ";ct(i)="MT"
i=i+1;cd(i)="GQPP";id(i)=61415;lt(i)= 20.93;ln(i)= -17.03;el(i)=   3;lnm(i)="NOUADHIBOU               ";st(i)="  ";ct(i)="MR"
i=i+1;cd(i)="GQNN";id(i)=61442;lt(i)= 18.10;ln(i)= -15.95;el(i)=   3;lnm(i)="NOUAKCHOTT               ";st(i)="  ";ct(i)="MR"
i=i+1;cd(i)="9999";id(i)=61600;lt(i)= 16.05;ln(i)= -16.45;el(i)=   4;lnm(i)="SAINT-LOUIS              ";st(i)="  ";ct(i)="SG"
i=i+1;cd(i)="GOOY";id(i)=61641;lt(i)= 14.73;ln(i)= -17.50;el(i)=  24;lnm(i)="DAKAR/YOFF               ";st(i)="  ";ct(i)="SN"
i=i+1;cd(i)="GOTT";id(i)=61687;lt(i)= 13.76;ln(i)= -13.68;el(i)=  50;lnm(i)="TAMBACOUNDA              ";st(i)="  ";ct(i)="SN"
i=i+1;cd(i)="GOGG";id(i)=61695;lt(i)= 12.55;ln(i)= -16.27;el(i)=  23;lnm(i)="ZIGUINCHOR               ";st(i)="  ";ct(i)="SN"
i=i+1;cd(i)="9999";id(i)=61831;lt(i)=  9.57;ln(i)= -13.62;el(i)=  49;lnm(i)="CONAKRY                  ";st(i)="  ";ct(i)="GN"
i=i+1;cd(i)="9999";id(i)=61901;lt(i)= -5.97;ln(i)=  -5.70;el(i)= 436;lnm(i)="ST. HELENA ISLAND        ";st(i)="  ";ct(i)="HE"
i=i+1;cd(i)="FHAW";id(i)=61902;lt(i)= -7.97;ln(i)= -14.40;el(i)=  79;lnm(i)="WIDE AWAKE FIELD         ";st(i)="  ";ct(i)="HE"
i=i+1;cd(i)="FJDG";id(i)=61967;lt(i)= -7.30;ln(i)=  72.40;el(i)=   3;lnm(i)="DIEGO GARCIA             ";st(i)="  ";ct(i)="IO"
i=i+1;cd(i)="9999";id(i)=61976;lt(i)=-15.88;ln(i)=  54.52;el(i)=  13;lnm(i)="SERGE-FROLOW/TROMELIN    ";st(i)="  ";ct(i)="RE"
i=i+1;cd(i)="FMEE";id(i)=61980;lt(i)=-20.88;ln(i)=  55.52;el(i)=  20;lnm(i)="SAINT DENIS              ";st(i)="  ";ct(i)="RE"
i=i+1;cd(i)="9999";id(i)=61995;lt(i)=-20.30;ln(i)=  57.50;el(i)= 425;lnm(i)="VACOAS                   ";st(i)="  ";ct(i)="MU"
i=i+1;cd(i)="9999";id(i)=61996;lt(i)=-37.80;ln(i)=  77.53;el(i)=  29;lnm(i)="MARTIN DE VIVIES         ";st(i)="  ";ct(i)="FR"
i=i+1;cd(i)="9999";id(i)=61998;lt(i)=-49.35;ln(i)=  70.25;el(i)=  30;lnm(i)="PORT-AUX-FRANCAIS        ";st(i)="  ";ct(i)="FR"
i=i+1;cd(i)="HLLT";id(i)=62010;lt(i)= 32.67;ln(i)=  13.15;el(i)=  81;lnm(i)="TRIPOLI INTL ARPT        ";st(i)="  ";ct(i)="LY"
i=i+1;cd(i)="HLLB";id(i)=62053;lt(i)= 32.08;ln(i)=  20.27;el(i)= 132;lnm(i)="BENINA/BENGHAZI          ";st(i)="  ";ct(i)="LY"
i=i+1;cd(i)="HLLS";id(i)=62124;lt(i)= 27.02;ln(i)=  14.43;el(i)= 432;lnm(i)="SEBHA (AUT)              ";st(i)="  ";ct(i)="LY"
i=i+1;cd(i)="HEMM";id(i)=62306;lt(i)= 31.33;ln(i)=  27.22;el(i)=  30;lnm(i)="MERSA MATRUH (MIL)       ";st(i)="  ";ct(i)="EG"
i=i+1;cd(i)="HEAR";id(i)=62337;lt(i)= 31.08;ln(i)=  33.75;el(i)=  32;lnm(i)="EL ARISH                 ";st(i)="  ";ct(i)="EG"
i=i+1;cd(i)="9999";id(i)=62378;lt(i)= 29.87;ln(i)=  31.33;el(i)= 141;lnm(i)="HELWAN                   ";st(i)="  ";ct(i)="EG"
i=i+1;cd(i)="HESN";id(i)=62414;lt(i)= 23.97;ln(i)=  32.78;el(i)= 194;lnm(i)="ASWAN (CIV/MIL)          ";st(i)="  ";ct(i)="EG"
i=i+1;cd(i)="9999";id(i)=62423;lt(i)= 27.05;ln(i)=  27.97;el(i)=  90;lnm(i)="FARAFRA                  ";st(i)="  ";ct(i)="EG"
i=i+1;cd(i)="HSSS";id(i)=62721;lt(i)= 15.60;ln(i)=  32.55;el(i)= 380;lnm(i)="KHARTOUM (CIV/MIL)       ";st(i)="  ";ct(i)="SD"
i=i+1;cd(i)="9999";id(i)=62403;lt(i)= 26.20;ln(i)=  32.75;el(i)=  96;lnm(i)="SOUTH OF VALLEY UNIV     ";st(i)="  ";ct(i)="EG"
i=i+1;cd(i)="9999";id(i)=62423;lt(i)= 27.05;ln(i)=  27.97;el(i)=  92;lnm(i)="FARAFRA (OASIS)          ";st(i)="  ";ct(i)="EG"
i=i+1;cd(i)="HAAB";id(i)=63450;lt(i)=  8.98;ln(i)=  38.80;el(i)=2355;lnm(i)="ADDIS ABABA/BOLE         ";st(i)="  ";ct(i)="ET"
i=i+1;cd(i)="HUEN";id(i)=63705;lt(i)=  0.05;ln(i)=  32.45;el(i)=1155;lnm(i)="ENTEBBE INTL ARPT        ";st(i)="  ";ct(i)="UG"
i=i+1;cd(i)="HKNC";id(i)=63741;lt(i)= -1.30;ln(i)=  36.75;el(i)=1798;lnm(i)="NAIROBI/DAGORETTI        ";st(i)="  ";ct(i)="KE"
i=i+1;cd(i)="9999";id(i)=63894;lt(i)= -6.88;ln(i)=  39.20;el(i)=  55;lnm(i)="DAR ES SALAAM            ";st(i)="  ";ct(i)="TN"
i=i+1;cd(i)="FSSS";id(i)=63985;lt(i)= -4.68;ln(i)=  55.53;el(i)=   4;lnm(i)="SEYCHELLES INTL          ";st(i)="  ";ct(i)="SC"
i=i+1;cd(i)="FCPP";id(i)=64400;lt(i)= -4.82;ln(i)=  11.90;el(i)=  17;lnm(i)="POINTE-NOIRE             ";st(i)="  ";ct(i)="CG"
i=i+1;cd(i)="FCBB";id(i)=64450;lt(i)= -4.25;ln(i)=  15.25;el(i)= 316;lnm(i)="BRAZZAVILLE/MAYA-MAYA    ";st(i)="  ";ct(i)="CG"
i=i+1;cd(i)="FCOU";id(i)=64458;lt(i)=  1.62;ln(i)=  16.05;el(i)= 352;lnm(i)="OUESSO                   ";st(i)="  ";ct(i)="CG"
i=i+1;cd(i)="FOOL";id(i)=64500;lt(i)=  0.45;ln(i)=   9.42;el(i)=  15;lnm(i)="LIBREVILLE/LEON MBA      ";st(i)="  ";ct(i)="GA"
i=i+1;cd(i)="FEFF";id(i)=64650;lt(i)=  4.40;ln(i)=  18.52;el(i)= 366;lnm(i)="BANGUI/M'POKO (MIL)      ";st(i)="  ";ct(i)="CF"
i=i+1;cd(i)="9999";id(i)=64665;lt(i)=  4.32;ln(i)=  21.18;el(i)= 406;lnm(i)="MOBAYE                   ";st(i)="  ";ct(i)="CF"
i=i+1;cd(i)="FTTJ";id(i)=64700;lt(i)= 12.13;ln(i)=  15.03;el(i)= 295;lnm(i)="NDJAMENA (CIV/MIL)       ";st(i)="  ";ct(i)="TD"
i=i+1;cd(i)="9999";id(i)=64750;lt(i)=  9.15;ln(i)=  18.38;el(i)= 365;lnm(i)="SARH                     ";st(i)="  ";ct(i)="CD"
i=i+1;cd(i)="9999";id(i)=64870;lt(i)=  7.35;ln(i)=  13.57;el(i)=1114;lnm(i)="NGAOUNDERE               ";st(i)="  ";ct(i)="CM"
i=i+1;cd(i)="FKKD";id(i)=64910;lt(i)=  4.00;ln(i)=   9.70;el(i)=   9;lnm(i)="DOUALA (CIV/MIL)         ";st(i)="  ";ct(i)="CM"
i=i+1;cd(i)="9999";id(i)=65046;lt(i)= 12.05;ln(i)=   8.53;el(i)= 476;lnm(i)="KANO                     ";st(i)="  ";ct(i)="NI"
i=i+1;cd(i)="DNAA";id(i)=65125;lt(i)=  9.25;ln(i)=   7.00;el(i)= 344;lnm(i)="ABUJA                    ";st(i)="  ";ct(i)="NI"
i=i+1;cd(i)="9999";id(i)=65257;lt(i)=  6.48;ln(i)=   7.55;el(i)= 142;lnm(i)="ENUGU                    ";st(i)="  ";ct(i)="NI"
i=i+1;cd(i)="DBBP";id(i)=65330;lt(i)=  9.35;ln(i)=   2.62;el(i)= 393;lnm(i)="PARAKOU                  ";st(i)="  ";ct(i)="BJ"
i=i+1;cd(i)="9999";id(i)=65306;lt(i)= 11.13;ln(i)=   2.93;el(i)= 292;lnm(i)="KANDI                    ";st(i)="  ";ct(i)="BJ"
i=i+1;cd(i)="9999";id(i)=65344;lt(i)=  6.35;ln(i)=   2.38;el(i)=   6;lnm(i)="COTONOU                  ";st(i)="  ";ct(i)="BJ"
i=i+1;cd(i)="9999";id(i)=65387;lt(i)=  6.17;ln(i)=   1.25;el(i)=  21;lnm(i)="LOME                     ";st(i)="  ";ct(i)="TG"
i=i+1;cd(i)="DGLE";id(i)=65418;lt(i)=  9.50;ln(i)=  -0.85;el(i)= 173;lnm(i)="TAMALE                   ";st(i)="  ";ct(i)="GH"
i=i+1;cd(i)="DHHH";id(i)=65503;lt(i)= 12.35;ln(i)=  -1.52;el(i)= 306;lnm(i)="OUAGADOUGOU              ";st(i)="  ";ct(i)="BF"
i=i+1;cd(i)="9999";id(i)=65510;lt(i)= 11.17;ln(i)=  -4.30;el(i)= 460;lnm(i)="BOBO-DIOULASSO           ";st(i)="  ";ct(i)="HV"
i=i+1;cd(i)="DIMN";id(i)=65548;lt(i)=  7.38;ln(i)=  -7.52;el(i)= 340;lnm(i)="MAN                      ";st(i)="  ";ct(i)="CI"
i=i+1;cd(i)="DIAP";id(i)=65578;lt(i)=  5.25;ln(i)=  -3.93;el(i)=   8;lnm(i)="ABIDJAN/PORT BOUET       ";st(i)="  ";ct(i)="CI"
i=i+1;cd(i)="9999";id(i)=67002;lt(i)=-11.53;ln(i)=  43.27;el(i)=  29;lnm(i)="MORONI/HAHAYA            ";st(i)="  ";ct(i)="IC"
i=i+1;cd(i)="9999";id(i)=67027;lt(i)=-15.67;ln(i)=  46.35;el(i)=  26;lnm(i)="MAJUNGA                  ";st(i)="  ";ct(i)="MG"
i=i+1;cd(i)="FMMI";id(i)=67083;lt(i)=-18.80;ln(i)=  47.48;el(i)=1276;lnm(i)="ANTANANARIVO/IVATO       ";st(i)="  ";ct(i)="MG"
i=i+1;cd(i)="9999";id(i)=67095;lt(i)=-18.12;ln(i)=  49.40;el(i)=   5;lnm(i)="TAMATAVE                 ";st(i)="  ";ct(i)="MG"
i=i+1;cd(i)="FMSD";id(i)=67197;lt(i)=-25.03;ln(i)=  46.95;el(i)=   9;lnm(i)="FT. DAUPHIN/TOLAGNA      ";st(i)="  ";ct(i)="MG"
i=i+1;cd(i)="FQNP";id(i)=67237;lt(i)=-15.09;ln(i)=  39.28;el(i)= 441;lnm(i)="NAMPULA                  ";st(i)="  ";ct(i)="MZ"
i=i+1;cd(i)="FLLC";id(i)=67666;lt(i)=-15.42;ln(i)=  28.32;el(i)=1280;lnm(i)="LUSAKA CITY AIRPORT      ";st(i)="  ";ct(i)="ZB"
i=i+1;cd(i)="9999";id(i)=67774;lt(i)=-17.83;ln(i)=  31.02;el(i)=1472;lnm(i)="HARARE/BELVEDERE         ";st(i)="  ";ct(i)="ZW"
i=i+1;cd(i)="9999";id(i)=67964;lt(i)=-20.15;ln(i)=  28.62;el(i)=1344;lnm(i)="BULAWAYO/GOETZ           ";st(i)="  ";ct(i)="ZW"
i=i+1;cd(i)="FBMN";id(i)=68032;lt(i)=-19.98;ln(i)=  23.42;el(i)= 900;lnm(i)="MAUN                     ";st(i)="  ";ct(i)="BW"
i=i+1;cd(i)="FBLT";id(i)=68040;lt(i)=-21.42;ln(i)=  25.60;el(i)= 985;lnm(i)="LETLHAKANE               ";st(i)="  ";ct(i)="BW"
i=i+1;cd(i)="9999";id(i)=68098;lt(i)=-22.96;ln(i)=  14.66;el(i)= 150;lnm(i)="WALVIS BAY AIRPORT       ";st(i)="  ";ct(i)="NM"
i=i+1;cd(i)="FAWW";id(i)=68110;lt(i)=-22.57;ln(i)=  17.10;el(i)=1725;lnm(i)="WINDHOEK/EROS(SAAF)      ";st(i)="  ";ct(i)="NA"
i=i+1;cd(i)="FAPB";id(i)=68174;lt(i)=-23.87;ln(i)=  29.45;el(i)=1222;lnm(i)="PIETERSBURG (SAAF)       ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="FBSK";id(i)=68240;lt(i)=-24.22;ln(i)=  25.92;el(i)=1005;lnm(i)="SIR SERETSE KHAMA        ";st(i)="  ";ct(i)="BW"
i=i+1;cd(i)="FAIR";id(i)=68263;lt(i)=-25.92;ln(i)=  28.22;el(i)=1500;lnm(i)="PRETORIA/IRENE           ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="9999";id(i)=68312;lt(i)=-26.53;ln(i)=  18.12;el(i)=1073;lnm(i)="KEETMENSHOOP             ";st(i)="  ";ct(i)="NM"
i=i+1;cd(i)="FBTS";id(i)=68328;lt(i)=-26.05;ln(i)=  22.45;el(i)=1000;lnm(i)="TSABONG                  ";st(i)="  ";ct(i)="BW"
i=i+1;cd(i)="FAUP";id(i)=68424;lt(i)=-28.40;ln(i)=  21.27;el(i)= 836;lnm(i)="UPINGTON/PIERRE VAN      ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="FABL";id(i)=68442;lt(i)=-29.10;ln(i)=  26.30;el(i)=1348;lnm(i)="BLOEMFONTEIN/HERTZOG     ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="FABM";id(i)=68461;lt(i)=-28.25;ln(i)=  28.33;el(i)=1682;lnm(i)="BETHLEHEM AIRPORT        ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="FASB";id(i)=68512;lt(i)=-29.67;ln(i)=  17.87;el(i)=1006;lnm(i)="SPRINGBOK                ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="9999";id(i)=68538;lt(i)=-30.67;ln(i)=  24.02;el(i)=1287;lnm(i)="DEAAR(UA)                ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="FADN";id(i)=68588;lt(i)=-29.97;ln(i)=  30.95;el(i)=   8;lnm(i)="DURBAN/LOUIS BOTHA       ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="FALE";id(i)=68592;lt(i)=-29.61;ln(i)=  31.12;el(i)= 109;lnm(i)="KING SHAKA INTL ARPT     ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="FACT";id(i)=68816;lt(i)=-33.98;ln(i)=  18.60;el(i)=  42;lnm(i)="CAPETOWN/DF MALAN        ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="FAPE";id(i)=68842;lt(i)=-33.98;ln(i)=  25.60;el(i)=  60;lnm(i)="PORT ELIZABETH           ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="FAGE";id(i)=68906;lt(i)=-40.35;ln(i)=  -9.88;el(i)=  54;lnm(i)="GOUGH ISLAND             ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="FAME";id(i)=68994;lt(i)=-46.88;ln(i)=  37.87;el(i)=  22;lnm(i)="MARION ISLAND            ";st(i)="  ";ct(i)="ZA"
i=i+1;cd(i)="KSTA";id(i)=69990;lt(i)= 38.99;ln(i)= -77.50;el(i)=  88;lnm(i)="STERLING                 ";st(i)=" V";ct(i)="US"
i=i+1;cd(i)=" SYA";id(i)=70414;lt(i)= 52.72;ln(i)= 174.10;el(i)= 30 ;lnm(i)="SHEMYA                   ";st(i)=" A";ct(i)="US"
i=i+1;cd(i)="9999";id(i)=71561;lt(i)= 44.23;ln(i)= -79.78;el(i)= 251;lnm(i)="EGBERT                   ";st(i)="  ";ct(i)="CN"
i=i+1;cd(i)="9999";id(i)=74004;lt(i)= 32.50;ln(i)=-114.00;el(i)= 231;lnm(i)="YUMA PROVING GROUND      ";st(i)=" A";ct(i)="US"
i=i+1;cd(i)="KNID";id(i)=74612;lt(i)= 35.68;ln(i)=-117.68;el(i)= 665;lnm(i)="CHINA LAKE NAF           ";st(i)="  ";ct(i)="US"
i=i+1;cd(i)="9999";id(i)=78807;lt(i)=  8.97;ln(i)= -79.66;el(i)=   9;lnm(i)="ALBROOK AFB/BALBOA       ";st(i)="  ";ct(i)="PA"
i=i+1;cd(i)="9999";id(i)=78955;lt(i)= 13.15;ln(i)= -59.62;el(i)= 113;lnm(i)="CARRIBBEAN MET INSTITUTE ";st(i)="  ";ct(i)="BA"
i=i+1;cd(i)="9999";id(i)=80001;lt(i)= 12.58;ln(i)= -81.70;el(i)=   4;lnm(i)="SESQUICENTENARIO         ";st(i)="  ";ct(i)="CO"
i=i+1;cd(i)="SKRH";id(i)=80035;lt(i)= 11.53;ln(i)= -72.93;el(i)=   4;lnm(i)="RIOHACHA/ALMIRANTE       ";st(i)="  ";ct(i)="CO"
i=i+1;cd(i)="SKBO";id(i)=80222;lt(i)=  4.70;ln(i)= -74.13;el(i)=2548;lnm(i)="BOGOTA/ELDORADO          ";st(i)="  ";ct(i)="CO"
i=i+1;cd(i)="9999";id(i)=80241;lt(i)=  4.55;ln(i)= -70.92;el(i)= 167;lnm(i)="GAVIOTAS                 ";st(i)="  ";ct(i)="CO"
i=i+1;cd(i)="SKLT";id(i)=80398;lt(i)= -4.17;ln(i)= -69.95;el(i)=  84;lnm(i)="LETICIA/VASQUEZ COBO     ";st(i)="  ";ct(i)="CO"
i=i+1;cd(i)="SVBS";id(i)=80413;lt(i)= 10.25;ln(i)= -67.65;el(i)= 437;lnm(i)="MARACAY                  ";st(i)="  ";ct(i)="VN"
i=i+1;cd(i)="9999";id(i)=80422;lt(i)= 10.67;ln(i)= -63.25;el(i)=  10;lnm(i)="CARUPANO                 ";st(i)="  ";ct(i)="VN"
i=i+1;cd(i)="SVCB";id(i)=80444;lt(i)=  8.15;ln(i)= -63.55;el(i)=   8;lnm(i)="CIUDAD BOLIVAR           ";st(i)="  ";ct(i)="VN"
i=i+1;cd(i)="SVSA";id(i)=80447;lt(i)=  7.85;ln(i)= -72.45;el(i)= 378;lnm(i)="SAN ANTONIO TACHIRA      ";st(i)="  ";ct(i)="VN" 
i=i+1;cd(i)="SVSR";id(i)=80450;lt(i)=  7.90;ln(i)= -67.42;el(i)=  48;lnm(i)="SAN FERNANDO DE APURE    ";st(i)="  ";ct(i)="VN"
i=i+1;cd(i)="SVSE";id(i)=80462;lt(i)=  4.60;ln(i)= -61.12;el(i)= 907;lnm(i)="SANTA ELENA UAIREN       ";st(i)="  ";ct(i)="VN"
i=i+1;cd(i)="SOCA";id(i)=81405;lt(i)=  4.83;ln(i)= -52.37;el(i)=   9;lnm(i)="CAYENNE/ROCHAMBEAU       ";st(i)="  ";ct(i)="GF"
i=i+1;cd(i)="9999";id(i)=81729;lt(i)= -3.29;ln(i)= -60.63;el(i)=  21;lnm(i)="MANACAPURU               ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82107;lt(i)= -0.12;ln(i)= -67.07;el(i)=  79;lnm(i)="SAO GABRIEL CACH         ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBBV";id(i)=82022;lt(i)=  2.83;ln(i)= -60.70;el(i)= 140;lnm(i)="BOA VISTA (CIV/MIL)      ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBTS";id(i)=82025;lt(i)= -2.22;ln(i)= -55.93;el(i)= 325;lnm(i)="TIRIOS                   ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBTS";id(i)=82026;lt(i)=  2.22;ln(i)= -55.95;el(i)= 326;lnm(i)="TIRIOS                   ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82099;lt(i)=  0.05;ln(i)= -51.07;el(i)=  17;lnm(i)="MACAPA AIRPORT           ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBBE";id(i)=82193;lt(i)= -1.38;ln(i)= -48.48;el(i)=  16;lnm(i)="BELEM/VAL DE CAES        ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82244;lt(i)= -2.43;ln(i)= -54.72;el(i)=  72;lnm(i)="SANTAREM                 ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82276;lt(i)= -2.32;ln(i)= -44.42;el(i)=  50;lnm(i)="ALCANTARA                ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBSL";id(i)=82281;lt(i)= -2.60;ln(i)= -44.23;el(i)=  53;lnm(i)="SAO LUIS/MARECHAL        ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBMN";id(i)=82332;lt(i)= -3.15;ln(i)= -59.98;el(i)=  84;lnm(i)="MANAUS/PONTA PELADA      ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82397;lt(i)= -3.73;ln(i)= -38.55;el(i)=  19;lnm(i)="FORTALEZA                ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBFN";id(i)=82400;lt(i)= -3.85;ln(i)= -32.42;el(i)=  56;lnm(i)="FERNANDO DE NORONHA      ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82411;lt(i)= -4.23;ln(i)= -69.92;el(i)= 120;lnm(i)="TABATINGA                ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82532;lt(i)= -5.83;ln(i)= -61.28;el(i)=  48;lnm(i)="MANICORE                 ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBNT";id(i)=82599;lt(i)= -5.92;ln(i)= -35.25;el(i)=  52;lnm(i)="NATAL/AUGUSTO SEVER      ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82678;lt(i)= -6.77;ln(i)= -43.02;el(i)= 138;lnm(i)="FLORIANO                 ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82705;lt(i)= -7.62;ln(i)= -72.67;el(i)= 170;lnm(i)="CRUZERIO DO SUL          ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBPV";id(i)=82824;lt(i)= -8.77;ln(i)= -63.92;el(i)=  88;lnm(i)="PORTO VELHO(CV/MIL)      ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82900;lt(i)= -8.05;ln(i)= -34.92;el(i)=   7;lnm(i)="RECIFE/CURADO            ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBAT";id(i)=82965;lt(i)= -9.87;ln(i)= -56.10;el(i)= 288;lnm(i)="ALTA FLORESTA            ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82983;lt(i)= -9.38;ln(i)= -40.48;el(i)= 370;lnm(i)="PETROLINA                ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBVH";id(i)=83208;lt(i)=-12.73;ln(i)= -60.13;el(i)= 652;lnm(i)="VILHENA                  ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=83229;lt(i)=-13.02;ln(i)= -38.52;el(i)=  51;lnm(i)="SALVADOR                 ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBLP";id(i)=83288;lt(i)=-13.27;ln(i)= -43.42;el(i)= 458;lnm(i)="BOM JESUS DA LAPA        ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBCY";id(i)=83362;lt(i)=-15.65;ln(i)= -56.10;el(i)= 182;lnm(i)="CUIABA/MARECHAL          ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBBR";id(i)=83378;lt(i)=-15.87;ln(i)= -47.93;el(i)=1061;lnm(i)="BRASILIA (CIV/MIL)       ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=83498;lt(i)=-17.73;ln(i)= -39.25;el(i)=   3;lnm(i)="CARAVELAS                ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBUL";id(i)=83525;lt(i)=-18.87;ln(i)= -48.22;el(i)= 922;lnm(i)="UBERLANDIA               ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=83554;lt(i)=-19.00;ln(i)= -57.67;el(i)= 142;lnm(i)="CORUMBA (AEROPORTO)      ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=83566;lt(i)=-19.62;ln(i)= -43.57;el(i)= 827;lnm(i)="CONFIS INTNL ARPT        ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBCG";id(i)=83612;lt(i)=-20.47;ln(i)= -54.67;el(i)= 556;lnm(i)="CAMPO GRANDE INTL        ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBVT";id(i)=83649;lt(i)=-20.26;ln(i)= -40.28;el(i)=   4;lnm(i)="VITORIA AEROPORTO        ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=83650;lt(i)=-20.50;ln(i)= -29.32;el(i)=   5;lnm(i)="TRINDADE ISLAND          ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=83708;lt(i)=-22.78;ln(i)= -45.20;el(i)= 537;lnm(i)="GUARATINGUETA            ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBGL";id(i)=83746;lt(i)=-22.82;ln(i)= -43.25;el(i)=   6;lnm(i)="GALEAO/RIO(CIV/MIL)      ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82765;lt(i)= -7.33;ln(i)= -47.47;el(i)= 212;lnm(i)="CAROLINA                 ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=82917;lt(i)=-10.00;ln(i)= -67.80;el(i)= 143;lnm(i)="RIO BRANCO               ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBLO";id(i)=83768;lt(i)=-23.33;ln(i)= -51.12;el(i)= 570;lnm(i)="LONDRINA AIRPORT         ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBMT";id(i)=83779;lt(i)=-23.52;ln(i)= -46.63;el(i)= 722;lnm(i)="MARTE (CIV/MIL)          ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBFI";id(i)=83827;lt(i)=-25.52;ln(i)= -54.58;el(i)= 243;lnm(i)="FOZ DO IGUACU ARPT       ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBCT";id(i)=83840;lt(i)=-25.52;ln(i)= -49.17;el(i)= 908;lnm(i)="CURITIBA/AFONSO PEN      ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="9999";id(i)=83899;lt(i)=-27.67;ln(i)= -48.55;el(i)=   5;lnm(i)="FLORIANOPOLIS (AEROPORTO)";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBUG";id(i)=83928;lt(i)=-29.78;ln(i)= -57.03;el(i)=  74;lnm(i)="URUGUAIANA/RUBEM         ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBSM";id(i)=83937;lt(i)=-29.72;ln(i)= -53.70;el(i)=  85;lnm(i)="SANTA MARIA AEROPORT     ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SBPA";id(i)=83971;lt(i)=-30.00;ln(i)= -51.18;el(i)=   3;lnm(i)="PORTO ALEGRE/SALGAD      ";st(i)="  ";ct(i)="BR"
i=i+1;cd(i)="SEST";id(i)=84008;lt(i)= -0.90;ln(i)= -89.60;el(i)=   6;lnm(i)="SAN CRISTOBAL ISL        ";st(i)="  ";ct(i)="EC"
i=i+1;cd(i)="9999";id(i)=84203;lt(i)= -2.15;ln(i)= -79.88;el(i)=   6;lnm(i)="GUAYAQUIL AEROPUERTO     ";st(i)="  ";ct(i)="EC"
i=i+1;cd(i)="9999";id(i)=84378;lt(i)= -3.73;ln(i)= -73.25;el(i)= 117;lnm(i)="MORONA                   ";st(i)="  ";ct(i)="PE"
i=i+1;cd(i)="9999";id(i)=84416;lt(i)= -5.18;ln(i)= -80.60;el(i)=  52;lnm(i)="PIURA GRUP7              ";st(i)="  ";ct(i)="PE"
i=i+1;cd(i)="SPIM";id(i)=84628;lt(i)=-12.00;ln(i)= -77.12;el(i)=  13;lnm(i)="LIMA/JORGE CHAVEZ        ";st(i)="  ";ct(i)="PE"
i=i+1;cd(i)="9999";id(i)=84629;lt(i)=-12.15;ln(i)= -77.00;el(i)=  80;lnm(i)="LAS PALMAS               ";st(i)="  ";ct(i)="PE"
i=i+1;cd(i)="SPTU";id(i)=84658;lt(i)=-12.62;ln(i)= -69.20;el(i)= 266;lnm(i)="PUERTO MALDONADO         ";st(i)="  ";ct(i)="PE"
i=i+1;cd(i)="9999";id(i)=84659;lt(i)=-12.63;ln(i)= -69.23;el(i)= 200;lnm(i)="PUERTO MALDONANDO BAMAL  ";st(i)="  ";ct(i)="PE"
i=i+1;cd(i)="SCFA";id(i)=85442;lt(i)=-23.43;ln(i)= -70.43;el(i)= 120;lnm(i)="ANTOFAGASTA/CERRO        ";st(i)="  ";ct(i)="CL"
i=i+1;cd(i)="SCIP";id(i)=85469;lt(i)=-27.15;ln(i)= -09.42;el(i)=  47;lnm(i)="EASTER ISLAND            ";st(i)="  ";ct(i)="CL"
i=i+1;cd(i)="SCEL";id(i)=85574;lt(i)=-33.38;ln(i)= -70.78;el(i)= 476;lnm(i)="PUDAHUEL/ARTURO          ";st(i)="  ";ct(i)="CL"
i=i+1;cd(i)="SCSN";id(i)=85586;lt(i)=-33.65;ln(i)= -71.62;el(i)=  75;lnm(i)="SANTO DOMINGO            ";st(i)="  ";ct(i)="CL"
i=i+1;cd(i)="SCTE";id(i)=85799;lt(i)=-41.42;ln(i)= -73.08;el(i)=  86;lnm(i)="PUERTO MONTT/TEPUAL      ";st(i)="  ";ct(i)="CL"
i=i+1;cd(i)="SCCI";id(i)=85934;lt(i)=-53.00;ln(i)= -70.85;el(i)=  37;lnm(i)="PUNTA ARENAS             ";st(i)="  ";ct(i)="CL"
i=i+1;cd(i)="9999";id(i)=86218;lt(i)=-25.28;ln(i)= -57.63;el(i)=  83;lnm(i)="SILVIO PETTIROSSI AIRPORT";st(i)="  ";ct(i)="PY"
i=i+1;cd(i)="SASA";id(i)=87047;lt(i)=-24.85;ln(i)= -65.48;el(i)=1216;lnm(i)="SALTA AIRPORT            ";st(i)="  ";ct(i)="AR"
i=i+1;cd(i)="SARE";id(i)=87155;lt(i)=-27.45;ln(i)= -59.05;el(i)=  52;lnm(i)="RESISTENCIA AIRPORT      ";st(i)="  ";ct(i)="AR"
i=i+1;cd(i)="SACO";id(i)=87344;lt(i)=-31.32;ln(i)= -64.22;el(i)= 474;lnm(i)="CORDOBA AIRPORT          ";st(i)="  ";ct(i)="AR"
i=i+1;cd(i)="SAME";id(i)=87418;lt(i)=-32.83;ln(i)= -68.78;el(i)= 704;lnm(i)="MENDOZA/EL PLUMERIL      ";st(i)="  ";ct(i)="AR"
i=i+1;cd(i)="SAEZ";id(i)=87576;lt(i)=-34.82;ln(i)= -58.53;el(i)=  20;lnm(i)="BUENOS AIRES/EZEIZA      ";st(i)="  ";ct(i)="AR"
i=i+1;cd(i)="SAZR";id(i)=87623;lt(i)=-36.57;ln(i)= -64.27;el(i)= 191;lnm(i)="SANTA ROSA AIRPORT       ";st(i)="  ";ct(i)="AR"
i=i+1;cd(i)="SAZN";id(i)=87715;lt(i)=-38.95;ln(i)= -68.13;el(i)= 271;lnm(i)="NEUQUEN AIRPORT          ";st(i)="  ";ct(i)="AR"
i=i+1;cd(i)="SAVC";id(i)=87860;lt(i)=-45.78;ln(i)= -67.50;el(i)=  46;lnm(i)="COMODORO RIVADAVIA       ";st(i)="  ";ct(i)="AR"
i=i+1;cd(i)="EGYP";id(i)=88889;lt(i)=-51.82;ln(i)= -58.45;el(i)=  73;lnm(i)="MOUNT PLEASANT ARPT      ";st(i)="  ";ct(i)="GS"
i=i+1;cd(i)="9999";id(i)=89002;lt(i)=-70.67;ln(i)=  -8.25;el(i)=  40;lnm(i)="VON-NEUMAYER G-BASE      ";st(i)="  ";ct(i)="DE"
i=i+1;cd(i)="9999";id(i)=89009;lt(i)=-90.00;ln(i)=  -0.00;el(i)=2830;lnm(i)="AMUNDSEN-SCOTT           ";st(i)="  ";ct(i)="US"
i=i+1;cd(i)="9999";id(i)=89022;lt(i)=-75.50;ln(i)= -26.65;el(i)=  30;lnm(i)="HALLEY BRI-BASE          ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=89055;lt(i)=-64.23;ln(i)= -56.72;el(i)= 198;lnm(i)="VICECOMODORO MARAM       ";st(i)="  ";ct(i)="AR"
i=i+1;cd(i)="9999";id(i)=89062;lt(i)=-67.57;ln(i)= -68.13;el(i)=  16;lnm(i)="ROTHERA PT BRI-BASE      ";st(i)="  ";ct(i)="GB"
i=i+1;cd(i)="9999";id(i)=89512;lt(i)=-70.77;ln(i)=  11.83;el(i)= 102;lnm(i)="NOVOLAZAREVSKAJA         ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=89532;lt(i)=-69.00;ln(i)=  39.58;el(i)=  21;lnm(i)="SYOWA JAPAN-BASE         ";st(i)="  ";ct(i)="JP"
i=i+1;cd(i)="9999";id(i)=89564;lt(i)=-67.60;ln(i)=  62.88;el(i)=  16;lnm(i)="MAWSON AUS-BASE          ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=89571;lt(i)=-68.57;ln(i)=  77.95;el(i)=  13;lnm(i)="DAVIS AUS-BASE           ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=89592;lt(i)=-66.55;ln(i)=  93.02;el(i)=  30;lnm(i)="MIRNYJ SOVIET-BASE       ";st(i)="  ";ct(i)="RU"
i=i+1;cd(i)="9999";id(i)=89611;lt(i)=-66.28;ln(i)= 110.52;el(i)=  41;lnm(i)="CASEY AUS-BASE           ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=89642;lt(i)=-66.67;ln(i)= 140.02;el(i)=  43;lnm(i)="DUMONT D'URVILLE         ";st(i)="  ";ct(i)="FR"
i=i+1;cd(i)="9999";id(i)=89662;lt(i)=-74.70;ln(i)= 164.10;el(i)=  80;lnm(i)="BAIA TERRA NOVA          ";st(i)="  ";ct(i)="IT"
i=i+1;cd(i)="9999";id(i)=89664;lt(i)=-77.85;ln(i)= 166.67;el(i)=  34;lnm(i)="MCMURDO USA-BASE         ";st(i)="  ";ct(i)="US"
i=i+1;cd(i)="PGAC";id(i)=91212;lt(i)= 13.55;ln(i)= 144.83;el(i)= 111;lnm(i)="GUAM,MARIANA IS          ";st(i)="  ";ct(i)="UM"
i=i+1;cd(i)="PTKK";id(i)=91334;lt(i)=  7.47;ln(i)= 151.85;el(i)=   2;lnm(i)="TRUK INTL/MOEN ISL       ";st(i)="  ";ct(i)="FM"
i=i+1;cd(i)="PTPN";id(i)=91348;lt(i)=  6.97;ln(i)= 158.22;el(i)=  46;lnm(i)="PONAPE ISLAND            ";st(i)="  ";ct(i)="FM"
i=i+1;cd(i)="9999";id(i)=91364;lt(i)=  9.40;ln(i)= 167.47;el(i)=   4;lnm(i)="ROI-NAMUR                ";st(i)="  ";ct(i)="MH"
i=i+1;cd(i)="PKWA";id(i)=91366;lt(i)=  8.73;ln(i)= 167.73;el(i)=   8;lnm(i)="KWAJALEIN/BUCHOLZ        ";st(i)="  ";ct(i)="MH"
i=i+1;cd(i)="PKMJ";id(i)=91376;lt(i)=  7.08;ln(i)= 171.38;el(i)=   3;lnm(i)="MAJURO/MARSHALL ISL      ";st(i)="  ";ct(i)="MH"
i=i+1;cd(i)="PTRO";id(i)=91408;lt(i)=  7.33;ln(i)= 134.48;el(i)=  33;lnm(i)="KOROR/PALAU ISLAND       ";st(i)="  ";ct(i)="FM"
i=i+1;cd(i)="PTYA";id(i)=91413;lt(i)=  9.48;ln(i)= 138.08;el(i)=  17;lnm(i)="YAP ISLAND               ";st(i)="  ";ct(i)="FM"
i=i+1;cd(i)="9999";id(i)=91517;lt(i)= -9.42;ln(i)= 159.97;el(i)=  56;lnm(i)="HONIARA                  ";st(i)="  ";ct(i)="SB"
i=i+1;cd(i)="9999";id(i)=91532;lt(i)= -0.53;ln(i)= 166.92;el(i)=   0;lnm(i)="NAURU ARC-2 (AUS)        ";st(i)="  ";ct(i)="SB"
i=i+1;cd(i)="NWWN";id(i)=91592;lt(i)=-22.27;ln(i)= 166.45;el(i)=  72;lnm(i)="NOUMEA                   ";st(i)="  ";ct(i)="NC"
i=i+1;cd(i)="NGTA";id(i)=91610;lt(i)=  1.35;ln(i)= 172.92;el(i)=   4;lnm(i)="TARAWA/BONRIKI INTL      ";st(i)="  ";ct(i)="KI"
i=i+1;cd(i)="NGFU";id(i)=91643;lt(i)= -8.52;ln(i)= 179.22;el(i)=   2;lnm(i)="FUNAFUTI INTL ARPT       ";st(i)="  ";ct(i)="TV"
i=i+1;cd(i)="NFFN";id(i)=91680;lt(i)=-17.75;ln(i)= 177.45;el(i)=  18;lnm(i)="NANDI/NADI INTL          ";st(i)="  ";ct(i)="FJ"
i=i+1;cd(i)="NTSU";id(i)=91765;lt(i)=-14.33;ln(i)=-170.72;el(i)=   3;lnm(i)="PAGO PAGO INTL ARPT      ";st(i)="  ";ct(i)="UM"
i=i+1;cd(i)="9999";id(i)=91801;lt(i)= -9.00;ln(i)=-158.05;el(i)=   1;lnm(i)="PENRHYN ISLAND           ";st(i)="  ";ct(i)="CK"
i=i+1;cd(i)="NCRG";id(i)=91843;lt(i)=-21.20;ln(i)=-159.82;el(i)=   7;lnm(i)="AVARUA/RAROTONGA IL      ";st(i)="  ";ct(i)="CK"
i=i+1;cd(i)="9999";id(i)=91925;lt(i)= -9.80;ln(i)=-139.03;el(i)=  52;lnm(i)="ATUONA                   ";st(i)="  ";ct(i)="PF"
i=i+1;cd(i)="NTAA";id(i)=91938;lt(i)=-17.55;ln(i)=-149.62;el(i)=   2;lnm(i)="TAHITI ISLAND/FAAA       ";st(i)="  ";ct(i)="PF"
i=i+1;cd(i)="9999";id(i)=91943;lt(i)=-14.48;ln(i)=-145.03;el(i)=   3;lnm(i)="TAKAROA ATOLL            ";st(i)="  ";ct(i)="PF"
i=i+1;cd(i)="9999";id(i)=91948;lt(i)=-23.13;ln(i)=-134.97;el(i)=  89;lnm(i)="RIKITEA                  ";st(i)="  ";ct(i)="PF"
i=i+1;cd(i)="NTAT";id(i)=91954;lt(i)=-23.35;ln(i)=-149.48;el(i)=   3;lnm(i)="TUBUAI ISLAND            ";st(i)="  ";ct(i)="PF"
i=i+1;cd(i)="9999";id(i)=91958;lt(i)=-27.62;ln(i)=-144.33;el(i)=   2;lnm(i)="RAPA ISLAND              ";st(i)="  ";ct(i)="PF"
i=i+1;cd(i)="AYPY";id(i)=92035;lt(i)= -9.43;ln(i)= 147.20;el(i)=  49;lnm(i)="PORT MORESBY W.O.        ";st(i)="  ";ct(i)="PG"
i=i+1;cd(i)="9999";id(i)=92044;lt(i)= -2.07;ln(i)= 147.43;el(i)=   5;lnm(i)="MOMOTE W.O.              ";st(i)="  ";ct(i)="PG"
i=i+1;cd(i)="NZWP";id(i)=93112;lt(i)=-36.78;ln(i)= 174.63;el(i)=  27;lnm(i)="WHENUAPAI (NZ-AFB)       ";st(i)="  ";ct(i)="NZ"
i=i+1;cd(i)="9999";id(i)=93291;lt(i)=-38.67;ln(i)= 177.98;el(i)=   8;lnm(i)="GISBORNE                 ";st(i)="  ";ct(i)="NZ"
i=i+1;cd(i)="9999";id(i)=93308;lt(i)=-39.02;ln(i)= 174.18;el(i)=  36;lnm(i)="NEW PLYMOUTH             ";st(i)="  ";ct(i)="NZ"
i=i+1;cd(i)="NZPP";id(i)=93417;lt(i)=-40.90;ln(i)= 174.98;el(i)=  12;lnm(i)="PARAPARAUMU AERO         ";st(i)="  ";ct(i)="NZ"
i=i+1;cd(i)="9999";id(i)=93614;lt(i)=-42.72;ln(i)= 170.98;el(i)=  40;lnm(i)="HOKITIKA                 ";st(i)="  ";ct(i)="NZ"
i=i+1;cd(i)="NZNV";id(i)=93844;lt(i)=-46.42;ln(i)= 168.33;el(i)=   1;lnm(i)="INVERCARGILL AERO        ";st(i)="  ";ct(i)="NZ"
i=i+1;cd(i)="NZCI";id(i)=93986;lt(i)=-43.95;ln(i)=-176.57;el(i)=  48;lnm(i)="CHATHAM ISL/TUUTA        ";st(i)="  ";ct(i)="NZ"
i=i+1;cd(i)="NZRN";id(i)=93997;lt(i)=-29.25;ln(i)=-177.92;el(i)=  49;lnm(i)="RAOUL ISL/KERMADEC       ";st(i)="  ";ct(i)="NZ"
i=i+1;cd(i)="YPDN";id(i)=94120;lt(i)=-12.40;ln(i)= 130.87;el(i)=  30;lnm(i)="DARWIN (CIV/MIL)         ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YDGV";id(i)=94150;lt(i)=-12.27;ln(i)= 136.82;el(i)=  54;lnm(i)="GOVE AIRPORT             ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94170;lt(i)=-12.63;ln(i)= 141.90;el(i)=   0;lnm(i)="WEIPA AERO               ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YBRM";id(i)=94203;lt(i)=-17.95;ln(i)= 122.22;el(i)=   9;lnm(i)="BROOME AIRPORT           ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94212;lt(i)=-18.23;ln(i)= 127.67;el(i)= 424;lnm(i)="HALLS CREEK              ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94238;lt(i)=-19.63;ln(i)= 134.18;el(i)= 377;lnm(i)="TENNANT CREEK AIRPORT    ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YBTL";id(i)=94294;lt(i)=-19.25;ln(i)= 146.75;el(i)=   6;lnm(i)="TOWNSVILLE(CIV/MIL)      ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94299;lt(i)=-16.30;ln(i)= 149.98;el(i)=   9;lnm(i)="WILLIS ISLAND            ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94300;lt(i)=-24.88;ln(i)= 113.67;el(i)=   8;lnm(i)="CARNARVON AIRPORT        ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YPLM";id(i)=94302;lt(i)=-22.23;ln(i)= 114.08;el(i)=   6;lnm(i)="LEARMOUTH                ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YPPD";id(i)=94312;lt(i)=-20.37;ln(i)= 118.62;el(i)=   6;lnm(i)="PORT HEDLAND ARRT        ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YBAS";id(i)=94326;lt(i)=-23.80;ln(i)= 133.90;el(i)= 541;lnm(i)="ALICE SPRINGS ARPT       ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YBMA";id(i)=94332;lt(i)=-20.67;ln(i)= 139.48;el(i)= 344;lnm(i)="MOUNT ISA AIRPORT        ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94346;lt(i)=-23.47;ln(i)= 144.25;el(i)= 191;lnm(i)="LONGREACH                ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94367;lt(i)=-21.12;ln(i)= 149.22;el(i)=  33;lnm(i)="MACKAY MO                ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YBRK";id(i)=94374;lt(i)=-23.38;ln(i)= 150.47;el(i)=  14;lnm(i)="ROCKHAMPTON AIRPORT      ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YPGN";id(i)=94403;lt(i)=-28.78;ln(i)= 114.70;el(i)=  34;lnm(i)="GERALDTON AIRPORT        ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YPMR";id(i)=94430;lt(i)=-26.60;ln(i)= 118.53;el(i)= 518;lnm(i)="MEEKATHARRA AIRPORT      ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94461;lt(i)=-25.03;ln(i)= 128.28;el(i)= 599;lnm(i)="GILES MET STATION        ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YBCV";id(i)=94510;lt(i)=-26.40;ln(i)= 146.27;el(i)= 304;lnm(i)="CHARLEVILLE ARPT         ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YBBN";id(i)=94578;lt(i)=-27.38;ln(i)= 153.10;el(i)=   5;lnm(i)="BRISBANE INTL ARPT       ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YPPH";id(i)=94610;lt(i)=-31.93;ln(i)= 115.95;el(i)=  29;lnm(i)="PERTH INTL/BELMONT       ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YPKG";id(i)=94637;lt(i)=-30.77;ln(i)= 121.45;el(i)= 360;lnm(i)="KALGOORLIE/BOULDER       ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94638;lt(i)=-33.82;ln(i)= 121.88;el(i)=  26;lnm(i)="ESPERANCE                ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94647;lt(i)=-31.67;ln(i)= 128.88;el(i)=  99;lnm(i)="EUCLA                    ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94653;lt(i)=-32.13;ln(i)= 133.70;el(i)=  22;lnm(i)="CENUDA                   ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YPWR";id(i)=94659;lt(i)=-31.13;ln(i)= 136.82;el(i)= 167;lnm(i)="WOOMERA (AUS-AFB)        ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YPAD";id(i)=94672;lt(i)=-34.93;ln(i)= 138.52;el(i)=   4;lnm(i)="ADELAIDE INTL ARPT       ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94693;lt(i)=-34.23;ln(i)= 142.08;el(i)=  52;lnm(i)="MILDURA AIRPORT          ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94711;lt(i)=-31.48;ln(i)= 145.82;el(i)= 265;lnm(i)="COBAR                    ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94750;lt(i)=-34.95;ln(i)= 150.53;el(i)= 110;lnm(i)="NOWRA RAN AIR STATION    ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YSSY";id(i)=94767;lt(i)=-33.95;ln(i)= 151.18;el(i)=   3;lnm(i)="SYDNEY INTL AIRPORT      ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YSWM";id(i)=94776;lt(i)=-32.78;ln(i)= 151.82;el(i)=   8;lnm(i)="WILLIAMTOWN(AUS-AB)      ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=94791;lt(i)=-30.32;ln(i)= 153.12;el(i)=   6;lnm(i)="COFFS HARBOUR MO         ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YPAL";id(i)=94802;lt(i)=-34.93;ln(i)= 117.80;el(i)=  69;lnm(i)="ALBANY AIRPORT           ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YMMG";id(i)=94821;lt(i)=-37.73;ln(i)= 140.78;el(i)=  69;lnm(i)="MOUNT GAMBIER ARPT       ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YMML";id(i)=94866;lt(i)=-37.67;ln(i)= 144.83;el(i)= 141;lnm(i)="MELBOURNE INTL ARPT      ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YSWG";id(i)=94910;lt(i)=-35.15;ln(i)= 147.45;el(i)= 213;lnm(i)="WAGGA WAGGA(CV/MIL)      ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YMHB";id(i)=94975;lt(i)=-42.83;ln(i)= 147.48;el(i)=  27;lnm(i)="HOBART AIRPORT           ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="ASLH";id(i)=94995;lt(i)=-31.53;ln(i)= 159.07;el(i)=   6;lnm(i)="LORD HOWE ISLAND         ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YSNF";id(i)=94996;lt(i)=-29.03;ln(i)= 167.93;el(i)= 109;lnm(i)="NORFOLK ISLAND ARPT      ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="YMMQ";id(i)=94998;lt(i)=-54.48;ln(i)= 158.93;el(i)=   6;lnm(i)="MACQUARIE ISLAND         ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=95527;lt(i)=-29.50;ln(i)= 149.83;el(i)= 214;lnm(i)="MOREE AMO                ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=96009;lt(i)=  5.13;ln(i)=  97.17;el(i)=  87;lnm(i)="LHOSEUMAWE               ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96011;lt(i)=  5.52;ln(i)=  95.42;el(i)=  21;lnm(i)="BANDA ATJEH/BLANGBINTANG ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96015;lt(i)=  4.25;ln(i)=  96.12;el(i)=  90;lnm(i)="MEULABOH/CUT NYAK DHIEN  ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WIMM";id(i)=96035;lt(i)=  3.57;ln(i)=  98.68;el(i)=  25;lnm(i)="MEDAN/POLONIA (MIL)      ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96073;lt(i)=  1.55;ln(i)=  98.88;el(i)=   3;lnm(i)="SIBOLGA/PINANGSORE       ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96075;lt(i)=  1.50;ln(i)=  97.63;el(i)=   6;lnm(i)="GUNUNG SITOLI            ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96091;lt(i)=  0.92;ln(i)= 104.53;el(i)=  18;lnm(i)="TANJUNGPINANG            ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96109;lt(i)=  0.47;ln(i)= 101.45;el(i)=  31;lnm(i)="PAKANBARU                ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96145;lt(i)=  3.20;ln(i)= 106.25;el(i)=   3;lnm(i)="TAREMPA                  ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WION";id(i)=96147;lt(i)=  3.95;ln(i)= 108.38;el(i)=   2;lnm(i)="RANAI (mil/civ)          ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WIMG";id(i)=96163;lt(i)= -0.88;ln(i)= 100.35;el(i)=   3;lnm(i)="PADANG/TABING            ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96171;lt(i)= -0.43;ln(i)= 102.45;el(i)=  19;lnm(i)="RENGAT DJAPURA           ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96179;lt(i)= -0.48;ln(i)= 104.58;el(i)=  31;lnm(i)="SINGKEP/DABO             ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96195;lt(i)= -1.63;ln(i)= 103.65;el(i)=  25;lnm(i)="DIAMBI PAALMERAH         ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96221;lt(i)= -2.90;ln(i)= 104.70;el(i)=  10;lnm(i)="PALEMBANG TALANGBETUTU   ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WIKK";id(i)=96237;lt(i)= -2.17;ln(i)= 106.13;el(i)=  33;lnm(i)="PANGKAL-PINANG           ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96249;lt(i)= -2.75;ln(i)= 107.75;el(i)=  44;lnm(i)="TANDJUNGPANDAN BULUHTUMBA";st(i)="NG";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96253;lt(i)= -3.87;ln(i)= 102.33;el(i)=  16;lnm(i)="BENGKULU PADANGKEMILING  ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96295;lt(i)= -5.27;ln(i)= 105.18;el(i)=  61;lnm(i)="TELUKBETUNG BRANTI       ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WBGI";id(i)=96315;lt(i)=  4.93;ln(i)= 114.93;el(i)=  15;lnm(i)="BRUNEI AIRPORT           ";st(i)="  ";ct(i)="BN"
i=i+1;cd(i)="WBGG";id(i)=96413;lt(i)=  1.48;ln(i)= 110.33;el(i)=  27;lnm(i)="KUCHING (CIV/MIL)        ";st(i)="  ";ct(i)="MY"
i=i+1;cd(i)="WBGB";id(i)=96441;lt(i)=  3.20;ln(i)= 113.03;el(i)=   5;lnm(i)="BINTULU/KALIMANTAN       ";st(i)="  ";ct(i)="MY"
i=i+1;cd(i)="WBKK";id(i)=96471;lt(i)=  5.93;ln(i)= 116.05;el(i)=   3;lnm(i)="KOTA KINABALU INTL       ";st(i)="  ";ct(i)="MY"
i=i+1;cd(i)="WBKW";id(i)=96481;lt(i)=  4.27;ln(i)= 117.88;el(i)=  20;lnm(i)="TAWAU/KALIMANTAN IL      ";st(i)="  ";ct(i)="MY"
i=i+1;cd(i)="WBKS";id(i)=96491;lt(i)=  5.90;ln(i)= 118.07;el(i)=  13;lnm(i)="SANDAKAN/KALIMANTAN      ";st(i)="  ";ct(i)="MY"
i=i+1;cd(i)="9999";id(i)=96509;lt(i)=  3.33;ln(i)= 117.57;el(i)=   6;lnm(i)="TANTKAN/DJUWATA          ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96535;lt(i)= -1.70;ln(i)= 109.30;el(i)=  15;lnm(i)="PALOH                    ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96581;lt(i)= -0.15;ln(i)= 109.40;el(i)=   3;lnm(i)="PONTIANAK                ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96595;lt(i)= -0.95;ln(i)= 114.90;el(i)=  60;lnm(i)="MUARATEWE                ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96607;lt(i)= -0.62;ln(i)= 117.15;el(i)= 230;lnm(i)="SAMARINDA                ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96633;lt(i)= -1.27;ln(i)= 116.90;el(i)=   3;lnm(i)="BALIKPAPAN SEPINGGAN     ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96645;lt(i)= -2.70;ln(i)= 110.70;el(i)=  25;lnm(i)="PANGKALAN                ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96655;lt(i)= -1.00;ln(i)= 114.00;el(i)=  27;lnm(i)="PALANGKA RAYA            ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96685;lt(i)= -3.45;ln(i)= 114.75;el(i)=  20;lnm(i)="BANDJARMASIN ULIN        ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96737;lt(i)= -6.12;ln(i)= 106.13;el(i)=  40;lnm(i)="SERANG                   ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96739;lt(i)= -6.23;ln(i)= 106.65;el(i)=  46;lnm(i)="CURUG                    ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WIII";id(i)=96749;lt(i)= -6.12;ln(i)= 106.65;el(i)=   8;lnm(i)="SOEKARNO-HATTA INTL      ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96791;lt(i)= -6.75;ln(i)= 108.27;el(i)=  50;lnm(i)="JATIWANGI                ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96797;lt(i)= -6.85;ln(i)= 109.15;el(i)=  10;lnm(i)="TEGAL                    ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96801;lt(i)= -7.37;ln(i)= 108.25;el(i)= 366;lnm(i)="TASIKMALAYA              ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96805;lt(i)= -7.73;ln(i)= 109.02;el(i)=   6;lnm(i)="TJILATJAP                ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96839;lt(i)= -6.98;ln(i)= 110.38;el(i)=   3;lnm(i)="SEMARANG KALIBANTENG     ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96925;lt(i)= -5.85;ln(i)= 112.63;el(i)=   3;lnm(i)="SANGKAPURA               ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WRSJ";id(i)=96935;lt(i)= -7.37;ln(i)= 112.77;el(i)=   3;lnm(i)="SURABAYA/JUANDA MIL      ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96973;lt(i)= -7.05;ln(i)= 113.97;el(i)=   3;lnm(i)="KALIANGET                ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=96987;lt(i)= -8.22;ln(i)= 114.38;el(i)=   5;lnm(i)="BANYUWANGI               ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="YPCC";id(i)=96996;lt(i)=-12.18;ln(i)=  96.82;el(i)=   3;lnm(i)="COCOS ISL INTL ARPT      ";st(i)="  ";ct(i)="AU"
i=i+1;cd(i)="9999";id(i)=97008;lt(i)=  3.58;ln(i)= 125.47;el(i)=  38;lnm(i)="TAHUNA/NAHA              ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WAMM";id(i)=97014;lt(i)=  1.53;ln(i)= 124.92;el(i)=  80;lnm(i)="MENADO/SAM RATULANG      ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97028;lt(i)=  1.02;ln(i)= 120.80;el(i)=   2;lnm(i)="TOLI-TOLI/LALOS          ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97048;lt(i)=  0.52;ln(i)= 123.07;el(i)=   2;lnm(i)="GORONTALO                ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WAML";id(i)=97072;lt(i)= -0.68;ln(i)= 119.73;el(i)=   6;lnm(i)="PALU/MUTIARA             ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97086;lt(i)= -0.90;ln(i)= 122.78;el(i)=   2;lnm(i)="LUWUK                    ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97096;lt(i)= -1.38;ln(i)= 120.73;el(i)=   2;lnm(i)="POSO/KASIGUNCU           ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97120;lt(i)= -2.50;ln(i)= 119.00;el(i)=   8;lnm(i)="MAJENE                   ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WAAA";id(i)=97180;lt(i)= -5.07;ln(i)= 119.55;el(i)=  14;lnm(i)="HASANUDDIN/UJUNG AB      ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97192;lt(i)= -5.47;ln(i)= 122.62;el(i)=   2;lnm(i)="BAU BAU                  ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97230;lt(i)= -8.75;ln(i)= 115.17;el(i)=   1;lnm(i)="DENPASAR                 ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97240;lt(i)= -8.57;ln(i)= 116.07;el(i)=   3;lnm(i)="AMPENAN                  ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97260;lt(i)= -8.52;ln(i)= 117.42;el(i)=   3;lnm(i)="SUMBAWA BESAR            ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97270;lt(i)= -8.55;ln(i)= 118.70;el(i)=   2;lnm(i)="BIMA/MOHAMMED SALAH      ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97300;lt(i)= -8.50;ln(i)= 122.23;el(i)=  32;lnm(i)="MAUMERE                  ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97340;lt(i)= -9.67;ln(i)= 120.33;el(i)=  12;lnm(i)="MAINGAPU MAUHAU          ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WRKK";id(i)=97372;lt(i)=-10.17;ln(i)= 123.67;el(i)= 108;lnm(i)="KUPANG/EL TARI           ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97430;lt(i)=  0.78;ln(i)= 127.38;el(i)=  23;lnm(i)="TERNATE                  ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97460;lt(i)= -1.62;ln(i)= 124.55;el(i)=   3;lnm(i)="LABUHA/TALIABU           ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97502;lt(i)= -0.93;ln(i)= 131.12;el(i)=   2;lnm(i)="SORONG                   ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97530;lt(i)= -0.88;ln(i)= 134.05;el(i)=   3;lnm(i)="MANOKWARI                ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WABB";id(i)=97560;lt(i)= -1.18;ln(i)= 136.12;el(i)=  11;lnm(i)="BIAK/FRANS KAISIEPO      ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97600;lt(i)= -2.08;ln(i)= 126.00;el(i)=   2;lnm(i)="SANANA                   ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97686;lt(i)= -4.07;ln(i)= 138.95;el(i)=1660;lnm(i)="WAMENA                   ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97690;lt(i)= -2.57;ln(i)= 140.48;el(i)=  99;lnm(i)="SENTANI                  ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WAPP";id(i)=97724;lt(i)= -3.70;ln(i)= 128.08;el(i)=  12;lnm(i)="AMBON/PATTIMURA          ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97748;lt(i)= -3.88;ln(i)= 130.90;el(i)=   3;lnm(i)="GESER                    ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97760;lt(i)= -3.67;ln(i)= 133.75;el(i)=   3;lnm(i)="KAIMANA, IRIAN BARAT     ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97810;lt(i)= -5.68;ln(i)= 132.75;el(i)=  12;lnm(i)="TUAL/DUMATUBUN           ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="9999";id(i)=97900;lt(i)= -7.98;ln(i)= 131.30;el(i)=  24;lnm(i)="SAUMLAKI                 ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="WAKK";id(i)=97980;lt(i)= -8.47;ln(i)= 140.38;el(i)=   3;lnm(i)="MERAUKE/MOPAH            ";st(i)="  ";ct(i)="ID"
i=i+1;cd(i)="RPLI";id(i)=98223;lt(i)= 18.18;ln(i)= 120.53;el(i)=   5;lnm(i)="LAOAG INTL(PH-ARMY)      ";st(i)="  ";ct(i)="PH"
i=i+1;cd(i)="RPUB";id(i)=98328;lt(i)= 16.37;ln(i)= 120.62;el(i)=1501;lnm(i)="BAGUIO                   ";st(i)="  ";ct(i)="PH"
i=i+1;cd(i)="9999";id(i)=98433;lt(i)= 14.50;ln(i)= 121.35;el(i)= 614;lnm(i)="TANAY                    ";st(i)="  ";ct(i)="PH"
i=i+1;cd(i)="RPMP";id(i)=98444;lt(i)= 13.13;ln(i)= 123.73;el(i)=  17;lnm(i)="LEGAZPI/LUZON ISL        ";st(i)="  ";ct(i)="PH"
i=i+1;cd(i)="RPVP";id(i)=98618;lt(i)=  9.74;ln(i)= 118.76;el(i)=  15;lnm(i)="PUERTO PRINCESA          ";st(i)="  ";ct(i)="PH"
i=i+1;cd(i)="RPMT";id(i)=98646;lt(i)= 10.30;ln(i)= 123.97;el(i)=  24;lnm(i)="MACTAN INTL(CIV/AF)      ";st(i)="  ";ct(i)="PH"
i=i+1;cd(i)="9999";id(i)=98747;lt(i)=  8.41;ln(i)= 124.61;el(i)= 188;lnm(i)="LUMBIA AIRPORT           ";st(i)="  ";ct(i)="PH"
i=i+1;cd(i)="RPMD";id(i)=98753;lt(i)=  7.12;ln(i)= 125.65;el(i)=  18;lnm(i)="DAVAO/FRANCISCO BAN      ";st(i)="  ";ct(i)="PH"

      num_RadSnd_Stat = i
      allocate(MR_Snd_cd(num_RadSnd_Stat))   ! Station code
      allocate(MR_Snd_id(num_RadSnd_Stat))   ! WMO station ID
      allocate(MR_Snd_lt(num_RadSnd_Stat))   ! Station latitude
      allocate(MR_Snd_ln(num_RadSnd_Stat))   ! Station longitude
      allocate(MR_Snd_el(num_RadSnd_Stat))   ! Station elevation
      allocate(MR_Snd_lnm(num_RadSnd_Stat))  ! Station long name
      allocate(MR_Snd_st(num_RadSnd_Stat))   ! Station state
      allocate(MR_Snd_ct(num_RadSnd_Stat))   ! Station country

      ! Now copy to MR module arrays
      MR_Snd_cd(1:num_RadSnd_Stat) = cd(1:num_RadSnd_Stat)
      MR_Snd_id(1:num_RadSnd_Stat) = id(1:num_RadSnd_Stat)
      MR_Snd_lt(1:num_RadSnd_Stat) = lt(1:num_RadSnd_Stat)
      MR_Snd_ln(1:num_RadSnd_Stat) = ln(1:num_RadSnd_Stat)
      MR_Snd_el(1:num_RadSnd_Stat) = el(1:num_RadSnd_Stat)
      MR_Snd_lnm(1:num_RadSnd_Stat)= lnm(1:num_RadSnd_Stat)
      MR_Snd_st(1:num_RadSnd_Stat) = st(1:num_RadSnd_Stat)
      MR_Snd_ct(1:num_RadSnd_Stat) = ct(1:num_RadSnd_Stat)

      end subroutine MR_Load_Radiosonde_Station_Data

!##############################################################################
!
!     MR_Get_Radiosonde_Station_Coord
!
!     This subroutine take a station code and returns the coordinates (lon,lat,elv)
!
!##############################################################################

      subroutine MR_Get_Radiosonde_Station_Coord(StatID,StatIdx,Stat_lon,Stat_lat,Stat_elv)

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_error,verbosity_info,stdin,&
         MR_Snd_cd,MR_Snd_lnm,MR_Snd_el,MR_Snd_lt,MR_Snd_ln,MR_Snd_id,&
         num_RadSnd_Stat

      ! This module requires Fortran 2003 or later
      use iso_fortran_env, only : &
         input_unit

      implicit none

      integer, parameter :: sp             = 4

      integer,      intent(in)  :: StatID
      integer,      intent(out) :: StatIdx
      real(kind=sp),intent(out) :: Stat_lon,Stat_lat,Stat_elv

      integer :: i
      logical :: FoundStation
      integer :: iostatus
      character(len=120) :: iomessage
      character(len=80) :: linebuffer080
      real(kind=sp)     :: tmp_sp

      integer :: io                           ! Index for output streams

      StatIdx = -1
      FoundStation = .false.
      do i=1,num_RadSnd_Stat
        if(StatID.eq.MR_Snd_id(i))then
          Stat_lon = real(MR_Snd_ln(i),kind=sp)
          Stat_lat = real(MR_Snd_lt(i),kind=sp)
          Stat_elv = real(MR_Snd_el(i),kind=sp)
          FoundStation = .true.
          StatIdx = i
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)" Radiosonde station: ",MR_Snd_cd(i),StatID,MR_Snd_lnm(i)
            write(outlog(io),'(a14,f10.3)')"  longitude = ",Stat_lon
            write(outlog(io),'(a14,f10.3)')"  latitude  = ",Stat_lat
            write(outlog(io),'(a14,f10.3)')"  elevation = ",Stat_elv
          endif;enddo
          exit
        endif
      enddo

      if(.not.FoundStation)then
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Radiosonde stations not found for ID ",StatID
          write(outlog(io),*)"  Please enter the station longitude (degrees) : "
        endif;enddo
        read(input_unit,*,iostat=iostatus,iomsg=iomessage) linebuffer080
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) tmp_sp
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: Expecting a real value for station longitude"
            write(errlog(io),*)"          You entered :",linebuffer080
            write(errlog(io),*)'MR System Message: '
            write(errlog(io),*)iomessage
          endif;enddo
          stop 1
        endif
        if(tmp_sp.lt.180.0_sp.or.tmp_sp.gt.360.0_sp)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: longitude must be between -180 and 360."
            write(errlog(io),*)"          You entered :",tmp_sp
          endif;enddo
          stop 1
        endif
        Stat_lon = tmp_sp

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)&
             "  Please enter the station latitude  (degrees) : "
        endif;enddo
        read(stdin,*,iostat=iostatus,iomsg=iomessage) linebuffer080
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) tmp_sp
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then          
            write(errlog(io),*)"MR ERROR: Expecting a real value for station latitude"
            write(errlog(io),*)"          You entered :",linebuffer080
            write(errlog(io),*)'MR System Message: '
            write(errlog(io),*)iomessage
          endif;enddo
          stop 1
        endif
        if(tmp_sp.lt.-90.0_sp.or.tmp_sp.gt.90.0_sp)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then                    
            write(errlog(io),*)"MR ERROR: latitude must be between -90 and 90."
            write(errlog(io),*)"          You entered :",tmp_sp
          endif;enddo
          stop 1
        endif
        Stat_lat = tmp_sp

        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)&
            "  Please enter the station elevation (m a.s.l.): "
        endif;enddo
        read(stdin,*,iostat=iostatus,iomsg=iomessage) linebuffer080
        read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) tmp_sp
        if(iostatus.ne.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then                    
            write(errlog(io),*)"MR ERROR: Expecting a real value for station elevation"
            write(errlog(io),*)"          You entered :",linebuffer080
            write(errlog(io),*)'MR System Message: '
            write(errlog(io),*)iomessage
          endif;enddo
          stop 1
        endif
        if(tmp_sp.lt.0.0_sp.or.tmp_sp.gt.8000.0_sp)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then                    
            write(errlog(io),*)"MR ERROR: elevation must be between 0 and 8000."
            write(errlog(io),*)"          You entered :",tmp_sp
          endif;enddo
          stop 1
        endif
        Stat_elv = tmp_sp

      endif

      end subroutine MR_Get_Radiosonde_Station_Coord

!##############################################################################
!
!     MR_Get_Radiosonde_Stations_InDomain
!
!     This subroutine compares the radiosonde station list with the computational
!     domain and lists the potentially available stations.
!
!##############################################################################

      subroutine MR_Get_Radiosonde_Stations_InDomain

      use MetReader,       only : &
         MR_nio,VB,outlog,verbosity_info,&
         MR_Snd_ct,MR_Snd_st,MR_Snd_lnm,MR_Snd_id,MR_Snd_ln,MR_Snd_lt,Snd_idx,&
         num_RadSnd_Stat,MR_nSnd_Locs,x_comp_sp,y_comp_sp,nx_comp,ny_comp

      implicit none

      integer, parameter :: sp             = 4

      integer :: nSnd_InDomain
      integer :: i,ii
      real(kind=sp)     :: de,dn
      real(kind=sp)     :: latLL,lonLL
      real(kind=sp)     :: latUR,lonUR

      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Radiosonde stations currently loaded:"
        do i=1,MR_nSnd_Locs
          ii = Snd_idx(i)
          write(outlog(io),*)MR_Snd_id(ii),trim(adjustl(MR_Snd_lnm(ii))),',',MR_Snd_st(ii),',',MR_Snd_ct(ii)
        enddo
      endif;enddo

      de = x_comp_sp(2)-x_comp_sp(1)
      dn = y_comp_sp(2)-y_comp_sp(1)
      lonLL = x_comp_sp(1)      -de
      lonUR = x_comp_sp(nx_comp)+de
      latLL = y_comp_sp(1)      -dn
      latUR = y_comp_sp(ny_comp)+dn
      if(lonLL.gt.180.0_sp)then
        lonLL = lonLL - 360.0_sp
        lonUR = lonUR - 360.0_sp
      endif

      nSnd_InDomain = 0
      do i=1,num_RadSnd_Stat
        if(latLL.lt.MR_Snd_lt(i).and.&
           latUR.gt.MR_Snd_lt(i).and.&
           lonLL.lt.MR_Snd_ln(i).and.&
           lonUR.gt.MR_Snd_ln(i))then
          if(nSnd_InDomain.eq.0)then
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)" Stations found in domain:"
            endif;enddo
            nSnd_InDomain = nSnd_InDomain + 1
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)&
                 MR_Snd_id(i),trim(adjustl(MR_Snd_lnm(i))),',',MR_Snd_st(i),',',MR_Snd_ct(i)
            endif;enddo
          endif
        endif
      enddo

      end subroutine MR_Get_Radiosonde_Stations_InDomain

!##############################################################################
!
!     MR_Calc_Snd_Weights_NearNeigh
!
!     This subroutine calculates the interpolation weights on the computational
!     grid from a sparse set of radiosonde staions using nearest neighbor.
!     First, the distance from each comp point to each station must be
!     calculated.
!     This subroutine is obsolete since the same effect can be achieved using
!     MR_Calc_Snd_Weights_InvDistWeight with a max station of 1 or appoximated
!     by using a high exponent.
!
!##############################################################################

!      subroutine MR_Calc_Snd_Weights_NearNeigh
!
!      use MetReader,       only : &
!         nx_comp,ny_comp,MR_nSnd_Locs,x_comp_sp,y_comp_sp,MR_Snd_ln,MR_Snd_lt,&
!         MR_Snd2Comp_map_wgt,MR_Snd2Comp_map_idx,Snd_idx
!
!      implicit none
!
!      integer, parameter :: sp             = 4
!      integer, parameter :: dp             = 8
!
!      integer :: i,j,k,idx
!      real(kind=dp)     :: lat1,lon1
!      real(kind=dp)     :: lat2,lon2
!
!      real(kind=dp),dimension(:,:,:),allocatable :: MR_Snd2Comp_dist
!
!      INTERFACE
!        real(kind=4) function MR_Haversine(lon1,lat1,lon2,lat2)
!          real(kind=4)  :: lon1
!          real(kind=4)  :: lat1
!          real(kind=4)  :: lon2
!          real(kind=4)  :: lat2
!        end function MR_Haversine
!      END INTERFACE
!
!      allocate(MR_Snd2Comp_dist(nx_comp,ny_comp,MR_nSnd_Locs))
!
!      MR_Snd2Comp_map_wgt = 0.0_sp
!
!      do i=1,nx_comp
!        lon1 = real(x_comp_sp(i),kind=dp)
!        do j=1,ny_comp
!          lat1 = real(y_comp_sp(j),kind=dp)
!          do k=1,MR_nSnd_Locs
!            lon2 = real(MR_Snd_ln(Snd_idx(k)),kind=8)+360.0_dp
!            lat2 = real(MR_Snd_lt(Snd_idx(k)),kind=8)
!            MR_Snd2Comp_dist(i,j,k) = MR_Haversine(lon1,lat1,lon2,lat2)
!            MR_Snd2Comp_map_idx(i,j,k) = k
!          enddo
!          idx = minloc(MR_Snd2Comp_dist(i,j,:),1)
!          MR_Snd2Comp_map_wgt(i,j,idx) = 1.0_sp
!        enddo
!      enddo
!
!      end subroutine MR_Calc_Snd_Weights_NearNeigh

!##############################################################################
!
!     MR_Calc_Snd_Weights_InvDistWeight
!
!     pext  = positive expo of inv dist func (default=4.0)
!     nstat = num of stations to consider (default=4; All=0; NearNeigh=1)
!
!     Sets MR_Snd2Comp_map_wgt
!
!     This subroutine calculates the interpolation weights on the computational
!     grid from a sparse set of radiosonde staions using an inverse distance.
!     First, the distance from each comp point to each station must be
!     calculated.
!
!##############################################################################

      subroutine MR_Calc_Snd_Weights_InvDistWeight(pexp,nstat)

      use MetReader,       only : &
         MR_Snd_lt,y_comp_sp,x_comp_sp,nx_comp,ny_comp,MR_nSnd_Locs,MR_Snd2Comp_map_idx,&
         Snd_idx,MR_Snd_ln,MR_Snd_lt,MR_Snd2Comp_map_wgt

      implicit none

      integer, parameter :: sp               = 4
      !integer, parameter :: dp               = 8

      real(kind=sp), intent(in)           :: pexp  ! positive expo of inv dist func
      integer,       intent(in), optional :: nstat ! num of stations to consider

      integer       :: i,j,k,l,idx
      real(kind=sp) :: lat1,lon1
      real(kind=sp) :: lat2,lon2
      real(kind=sp) :: norm
      real(kind=sp) :: dist
      real(kind=sp) :: lenscale
      integer :: maxstations
      !logical :: inStencil
      !real(kind=sp) :: domwgt
      logical,dimension(:),allocatable :: StatInStencil

      real(kind=sp),dimension(:,:,:),allocatable :: MR_Snd2Comp_dist

      INTERFACE
        real(kind=4) function MR_Haversine(lon1,lat1,lon2,lat2)
          real(kind=4)  :: lon1
          real(kind=4)  :: lat1
          real(kind=4)  :: lon2
          real(kind=4)  :: lat2
        end function MR_Haversine
      END INTERFACE

      ! First, get a domain length scale so that our weighting terms do not approach
      ! machine precision.  Using the diaganol distance
      lenscale = MR_Haversine(x_comp_sp(1),y_comp_sp(1),x_comp_sp(nx_comp),y_comp_sp(ny_comp))

      ! Figure out how many stations to use in the inv.dist. weighting
      if(present(nstat))then
        if(nstat.eq.0)then
          maxstations = MR_nSnd_Locs
        else
          maxstations = nstat
        endif
      else
        maxstations = MR_nSnd_Locs
      endif
      allocate(StatInStencil(MR_nSnd_Locs))

      allocate(MR_Snd2Comp_dist(nx_comp,ny_comp,MR_nSnd_Locs))
      MR_Snd2Comp_map_wgt = 0.0_sp
      MR_Snd2Comp_map_idx = 0

      ! Loop over all nodes on the computational grid
      do i=1,nx_comp
        lon1 = real(x_comp_sp(i),kind=sp)
        do j=1,ny_comp
          lat1 = real(y_comp_sp(j),kind=sp)
          ! Loop over all the sonde locations and get distances to the comp node
          do k=1,MR_nSnd_Locs
            lon2 = MR_Snd_ln(Snd_idx(k))+360.0_sp
            lat2 = MR_Snd_lt(Snd_idx(k))
            MR_Snd2Comp_dist(i,j,k) = MR_Haversine(lon1,lat1,lon2,lat2)/lenscale
            MR_Snd2Comp_map_wgt(i,j,k) = 1.0_sp/(MR_Snd2Comp_dist(i,j,k)**pexp)
            !norm = norm + MR_Snd2Comp_map_wgt(i,j,k)
            !MR_Snd2Comp_map_idx(i,j,k) = k
          enddo
          ! For this comp node, only use the maxstations closest sonde sites
          StatInStencil(1:MR_nSnd_Locs) = .false.
          idx = 0
          do l=1,maxstations
            ! Find the distance to the furthest station, then look through
            ! all the stations to find the closest not yet accounted for
            dist = maxval(MR_Snd2Comp_dist(i,j,1:MR_nSnd_Locs))
            ! Looking for the next closest station
            do k=1,MR_nSnd_Locs
              if(StatInStencil(k)) cycle ! skip to next k if this is already in stencil
              if(MR_Snd2Comp_dist(i,j,k).le.dist)then
                dist = MR_Snd2Comp_dist(i,j,k)
                idx = k
              endif
            enddo
            ! Log this station and mark it as part of the stencil
            MR_Snd2Comp_map_idx(i,j,l) = idx
            StatInStencil(idx) = .true.
          enddo

          ! Zero-out any sonde points that are not part of the stencil
!          do k=1,MR_nSnd_Locs
!            inStencil = .false.
!            do l=1,maxstations
!              if(MR_Snd2Comp_map_idx(i,j,l).eq.k)inStencil=.true.
!            enddo
!            if(.not.inStencil)MR_Snd2Comp_map_wgt(i,j,k) = 0.0_sp
!          enddo

          ! Zero-out any sonde points that do not contribute more that 0.1%
!          domwgt = maxval(MR_Snd2Comp_map_wgt(i,j,:))
!          do k=1,MR_nSnd_Locs
!            if(MR_Snd2Comp_map_wgt(i,j,k)/domwgt.lt.0.001_sp)MR_Snd2Comp_map_wgt(i,j,k)=0.0_sp
!          enddo

          ! Now find the normalization factor
          norm = 0.0_sp
          do l=1,maxstations
            idx = MR_Snd2Comp_map_idx(i,j,l)
            norm = norm + MR_Snd2Comp_map_wgt(i,j,idx)
          enddo
          ! Finally, normalize the requested weights and zero the others
          do k=1,MR_nSnd_Locs
            if(StatInStencil(k))then
              MR_Snd2Comp_map_wgt(i,j,k) = MR_Snd2Comp_map_wgt(i,j,k)/norm
            else
              MR_Snd2Comp_map_wgt(i,j,k) = 0.0_sp
            endif
          enddo
        enddo
      enddo

      end subroutine MR_Calc_Snd_Weights_InvDistWeight

!##############################################################################
!
!    MR_Haversine
!
!    Calculates the distance in km given to lat/lon points.  Assumes a spherical
!    Earth and uses the Haversine function.
!
!##############################################################################

      function MR_Haversine(lon1,lat1,lon2,lat2)

      use MetReader,       only : &
         MR_RAD_EARTH,MR_DEG2RAD

      implicit none

      real(kind=4)  :: MR_Haversine
      real(kind=4)  :: lon1,lat1,lon2,lat2
      real(kind=4)  :: lam1,psi1,lam2,psi2
      real(kind=4)  :: dpsi,dlam
      real(kind=4)  :: a,c
      real(kind=4)  :: fac1,fac2

      ! From https://www.movable-type.co.uk/scripts/latlong.html
      lam1 = lon1*MR_DEG2RAD
      psi1 = lat1*MR_DEG2RAD
      lam2 = lon2*MR_DEG2RAD
      psi2 = lat2*MR_DEG2RAD
      dpsi = psi2-psi1
      dlam = lam2-lam1

      a = sin(0.5_4*dpsi)*sin(0.5_4*dpsi) + &
          sin(0.5_4*dlam)*sin(0.5_4*dlam) * cos(psi1)*cos(psi2)
      fac1 = sqrt(      a)
      fac2 = sqrt(1.0_4-a)
      c = 2.0_4*atan2( fac1, fac2 )

      MR_Haversine = MR_RAD_EARTH * c

      end function MR_Haversine

!##############################################################################
