      module MetReader

      ! This file is the output of the simple script
!get_version.sh 
!  echo -n "      character(len=40),parameter,public :: MR_GitComID ='" > version.h
!  git log -n 1 | grep commit | cut -f 2 -d' ' | tr -d $'\n' >> version.h
!  echo -n "'" >> version.h
!    which set the variable containing the git commit ID : MR_GitComID
#include "version.h"

      integer, parameter,private :: sp        = 4 ! single precision
      integer, parameter,private :: dp        = 8 ! double precision

      integer, parameter         :: MR_VERB = 1

      integer, parameter :: MR_MAXVARS        = 50 ! Maximum number of variables in fixed arrays

      real(kind=dp), parameter :: MR_EPS_SMALL  = 1.0e-7_dp  ! Small number
      real(kind=dp), parameter :: MR_EPS_TINY   = 1.0e-12_dp ! Very small number

      real(kind=sp), parameter :: RAD_EARTH_MET = 6371.229_sp ! Radius of Earth in km
      real(kind=sp), parameter :: DEG2RAD_MET   = 1.7453292519943295e-2_sp
      real(kind=sp), parameter :: MR_MIN_DZ     = 1.0e-4_sp   ! used in reassigning z for low-level negative GPH

      integer,public :: MR_global_essential    = 6
      integer,public :: MR_global_production   = 6
      integer,public :: MR_global_debug        = 6
      integer,public :: MR_global_info         = 6
      integer,public :: MR_global_log          = 9
      integer,public :: MR_global_error        = 0

      logical        :: MR_useCompGrid         = .true. ! Reset this to .false. if you only need the Met grid
      logical        :: MR_useCompTime         = .true. ! Reset this to .false. if you only need the time of the file
      logical        :: MR_useCompH            = .true.
      logical        :: MR_useCompP            = .true.

      integer,public :: MR_iwind       !     MR_IWIND specifies the type of wind input to the model:
                             !   MR_IWIND=1 if a 1-D wind sounding is use, 
                             !           =2 if a 3-D grid is read from a ASCII file.
                             !           =3 if a single, multistep 3-D file is used
                             !           =4 if multiple 3-D NetCDF files are used
                             !           =5 if multiple file with multiple steps are used
      integer,public :: MR_iwindformat !      MR_iwindformat specifies the NWP product
                                       !  0 Custom format based on template
                                       !  1 ASCII profile
                                       !  2 Radiosonde data
                                       !  3 NARR3D 32 km North America files (32 km) :: ds608.0
                                       !  4 NAM Regional North America 221 (32 km)
                                       !  5 NAM AK 216  (45 km)
                                       !  6 NAM Regional 104 (90 km)
                                       !  7 NAM CONUS 212 (40 km)
                                       !  8 NAM CONUS 218 (12 km)
                                       !  9 NAM CONUS 227 (5.08 km)
                                       ! 10 NAM AK 242 (11.25 km)
                                       ! 11 NAM HI 196 (2.5 km)
                                       ! 12 NAM AK 198 (5.953 km)
                                       ! 13 NAM AK 91 (2.976 km)
                                       ! 14 NAM CONUS 1227 (3.0 km)
                                       ! 20 GFS 0.5
                                       ! 21 GFS 1.0
                                       ! 22 GFS 0.25
                                       ! 23 NCEP / DOE reanalysis 2.5 degree files  :: ds091.0
                                       ! 24 NASA-MERRA-2 reanalysis 0.625/0.5 degree files
                                       ! 25 NCEP/NCAR reanalysis 2.5 degree files   :: ds090.0  iwind=4 or 5
                                       ! 26 JRA-55                                  :: ds628.0  iwind=5
                                       ! 27 NOAA-CIRES reanalysis 2.0 degree files  :: ds131.2,3  iwind=5
                                       ! 28 ECMWF Interim Reanalysis (ERA-Interim)  :: ds627.0  requires catted grib files
                                       ! 29 ECMWF ERA5                              :: ds633.0  iwind=5
                                       ! 30 ECMWF 20-Century (ERA-20C)              :: ds626.0  iwind=5
                                       ! 32 Air Force Weather Agency subcenter = 0
                                       ! 33 CCSM3.0 Community Atmosphere Model (CAM)
                                       ! 40 NASA-GEOS Cp
                                       ! 41 NASA-GEOS Np
                                       ! 50 WRF - output
      logical,public :: MR_Use_RDA     = .true.
      integer,public :: MR_iversion    = -1 ! version of the product; sometimes the grid used changes with newer
                                            ! versions. -1 indicates latest version, but can be overridden here
      integer,public :: MR_iGridCode   !   MR_iGridCode specifies the NCEP grid described in:
                                    !   http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html
      integer,public :: MR_idataFormat !   Specifies the data model used
                                    !    =1 ASCII
                                    !    =2 netcdf
                                    !    =3 grib1 or grib2

        ! These variables describe the full list of windfiles read
      integer                                       ,public :: MR_iwindfiles           ! number of files provided
#ifdef USEPOINTERS
!      character(len=130), allocatable,dimension(:)  ,public :: fc_windfilename
      character(len=130), pointer    ,dimension(:)  ,public :: MR_windfiles            ! name of file
#else
      character(len=130), allocatable,dimension(:)  ,public :: MR_windfiles            ! name of file
#endif
      character(len=130)                            ,public :: MR_iw5_root
      character(len=42)                             ,public :: MR_iw5_prefix
      character(len=24)                             ,public :: MR_iw5_suffix1  ! e.g. 1912060100_1912063021.nc
      character(len=24)                             ,public :: MR_iw5_suffix2
      real(kind=dp)                                 ,public :: MR_iw5_hours_per_file

#if USEPOINTERS      
      real(kind=dp)     , pointer,dimension(:)  ,    public :: MR_windfile_starthour
      real(kind=dp)     , pointer,dimension(:,:),    public :: MR_windfile_stephour
#else
      real(kind=dp)     , allocatable,dimension(:)  ,public :: MR_windfile_starthour   ! start hour of the file
      real(kind=dp)     , allocatable,dimension(:,:),public :: MR_windfile_stephour    ! offset hours of step
#endif

      integer           , allocatable,dimension(:)  ,public :: MR_windfiles_nt_fullmet ! number of steps in files
      character(len=80)                             ,public :: MR_iwf_template         ! name of the template file
      logical           , allocatable,dimension(:)  ,public :: MR_windfiles_Have_GRIB_index
      character(len=130), allocatable,dimension(:)  ,public :: MR_windfiles_GRIB_index ! name of grib index file

        ! These variables are a list of the same data above, but specific to the simulation
        ! duration
      integer                                     ,public :: MR_MetSteps_Total
      integer                                     ,public :: MR_iMetStep_Now
      character(len=130), allocatable,dimension(:),public :: MR_MetStep_File
      integer           , allocatable,dimension(:),public :: MR_MetStep_findex
      integer           , allocatable,dimension(:),public :: MR_MetStep_tindex
      real(kind=dp)     , allocatable,dimension(:),public :: MR_MetStep_Hour_since_baseyear
      real(kind=dp)     , allocatable,dimension(:),public :: MR_MetStep_Interval
      integer           , allocatable,dimension(:),public :: MR_MetStep_year
      integer           , allocatable,dimension(:),public :: MR_MetStep_month
      integer           , allocatable,dimension(:),public :: MR_MetStep_day
      integer           , allocatable,dimension(:),public :: MR_MetStep_DOY
      real(kind=dp)     , allocatable,dimension(:),public :: MR_MetStep_Hour_Of_Day

      logical                                             :: MR_Reannalysis   = .false.
      integer           ,dimension(:), allocatable,public :: MR_iwind5_year

      !    Native grid of Met file using Pressure as vertical coordinate
#ifdef USEPOINTERS
      integer      ,dimension(:,:)  ,pointer, public :: MR_dum2d_met_int   => null() ! Used for categorical variables
      real(kind=sp),dimension(:,:)  ,pointer, public :: MR_dum2d_met       => null() ! Used for surface variables
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d_metP      => null()
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d2_metP     => null()
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_geoH_metP_last  => null() ! These are needed for compH interpolation
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_geoH_metP_next  => null() 

      real(kind=sp),dimension(:,:,:),pointer, public :: MR_vx_metP_last  => null() !These might need to be stored to avoid a 
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_vx_metP_next  => null() !second reading
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_vy_metP_last  => null() 
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_vy_metP_next  => null() 
#else
      integer      ,dimension(:,:)  ,allocatable, public :: MR_dum2d_met_int   ! Used for categorical variables
      real(kind=sp),dimension(:,:)  ,allocatable, public :: MR_dum2d_met       ! Used for surface variables
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_metP
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d2_metP
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_geoH_metP_last  ! These are needed for compH interpolation
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_geoH_metP_next

      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_vx_metP_last  !These might need to be stored to avoid a 
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_vx_metP_next  !second reading
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_vy_metP_last
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_vy_metP_next
#endif
      logical, public :: MR_Save_Velocities = .false.

      real(kind=sp), public :: MR_Max_geoH_metP_predicted   ! Highest expected height in Met file, based on pressure
      real(kind=sp) :: Max_geoH_metP_last
      real(kind=sp) :: Max_geoH_metP_next
      real(kind=sp) :: Max_geoH_metP
      real(kind=sp) :: Suppl_H
      !Parameter iHeightHandler specifies what to do if the maximum height
      !of the simulation region exceeds the maximum height in the wind files.
      !If iHeightHandler = 1, stop the program if the plume height exceeds mesoscale height
      !                    2, wind direction at levels above the highest node 
      !                       equal that of the highest node.  Temperatures in the
      !                       upper nodes don't change between 11 and 20 km; above
      !                       20 km they increase by 2 C/km, as in the Standard
      !                       atmosphere.  A warning is written to the log file.
      integer       :: MR_iHeightHandler = 2

        ! The following are used for sonde data
      integer :: MR_nSnd_Locs      = 1  ! Number of Sonde locations
      integer :: MR_Snd_nt_fullmet = 1  ! Number of times at the Sonde locations
      integer :: MR_Snd_nvars      = 5  ! Number of Sonde variables (P,H,U,V,T,+user-specified)
      logical :: Snd_Have_PT    = .false.
      logical :: Snd_Have_Coord = .false.  ! If the 1-d data have the optional projection params
                                           ! then it will be used, otherwise vel will be relative
                                           ! to comp grid.
      real(kind=sp),dimension(:,:,:,:),allocatable, public :: MR_SndVars_metP   ! (MR_nSnd_Locs,MR_Snd_nt_fullmet,MR_Snd_nvars,300)
      integer      ,dimension(:,:)    ,allocatable, public :: MR_Snd_np_fullmet ! Number of pressure values for each location/time
      integer      ,dimension(:)      ,allocatable, public :: MR_SndVarsID      ! Lists which vars are in which columns of  MR_SndVars_metP
      real(kind=sp),dimension(:,:,:)  ,allocatable, public :: MR_Snd2Comp_tri_map_wgt ! weights of nearby sondes for every comp point
      integer      ,dimension(:,:,:)  ,allocatable, public :: MR_Snd2Comp_tri_map_idx ! sonde index of weights

      !    Native grid of Met file using Height as vertical coordinate
      !    (resampled onto z-gridpoints of computational grid)
#ifdef USEPOINTERS
      real(kind=sp),dimension(:,:,:),    pointer, public :: MR_dum3d_metH => null()
#else
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_metH
#endif
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_u_ER_metP ! For the cases where Met is proj and comp
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_v_ER_metP !  different we need to rotate so these
                                                                         !  store Earth-Relative velocities on MetP

      !    Computations grid
#ifdef USEPOINTERS
      integer      ,dimension(:,:)  ,pointer, public :: MR_dum2d_comp_int => null() ! Used for categorical variables
      real(kind=sp),dimension(:,:)  ,pointer, public :: MR_dum2d_comp     => null() ! Used for surface variables
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d_compP    => null() ! 
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d_compP_2  => null() ! 
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d_compH    => null() ! 
      real(kind=sp),dimension(:,:,:),pointer, public :: MR_dum3d_compH_2  => null() ! Used only when a vector field
                                                                               ! rotation is needed
#else
      integer      ,dimension(:,:)  ,allocatable, public :: MR_dum2d_comp_int  ! Used for categorical variables
      real(kind=sp),dimension(:,:)  ,allocatable, public :: MR_dum2d_comp      ! Used for surface variables
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_compP     !
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_compP_2   !
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_compH     !
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_compH_2   ! Used only when a vector field
                                                                               ! rotation is needed
#endif

      integer       :: nt_fullmet    ! length of t of file
      integer       :: nx_fullmet    ! length of x or lon of full met grid
      integer       :: ny_fullmet    ! length of y or lat of full met grid
      integer       :: np_fullmet    ! length of pressure of full met grid for H,U,V
      integer       :: np_fullmet_pad = 1 ! We might need to pad the top of the pressure grid.
      integer       :: neta_fullmet  ! Only used by WRF
      real(kind=sp) :: Pressure_Conv_Fac = 100.0_sp  ! factor for converting from Met file units to Pa

#ifdef USEPOINTERS
      real(kind=sp),dimension(:),   pointer :: x_fullmet_sp    => null() ! x-coordinates of full met grid
      real(kind=sp),dimension(:),   pointer :: y_fullmet_sp    => null() ! y-coordinates of full met grid
      real(kind=sp),dimension(:),   pointer :: p_fullmet_sp    => null() ! z-coordinates of full met grid for H
      real(kind=sp),dimension(:,:), pointer :: levs_fullmet_sp => null() ! This hold each of the numbered level coordinates:
                                                                         !    i.e. isobaric, isobaric1, isobaric2, but also
                                                                         !    height_above_ground, depth_below_surface_layer, etc.
                                                                         !    p_fullmet_sp,p_fullmet_[Vz,RH]sp are copies of one
                                                                         !    of the slices
      integer,dimension(:), pointer :: nlevs_fullmet                     ! length of z coordinate
      integer,dimension(:), pointer :: levs_code                         ! code indicating how to map to the GPH grid
                                                                         !   0 = no mapping (not a pressure coordinate)
                                                                         !   1 = one-to-one mapping (U,V)
                                                                         !   2 = upper truncation (missing upper levels)
                                                                         !   3 = interpolation (missing mid-levels)
                                                                         !   4 = more levels than GPH grid
#else
      real(kind=sp),dimension(:),  allocatable :: x_fullmet_sp    ! x-coordinates of full met grid
      real(kind=sp),dimension(:),  allocatable :: y_fullmet_sp    ! y-coordinates of full met grid
      real(kind=sp),dimension(:),  allocatable :: p_fullmet_sp    ! z-coordinates of full met grid for H
      real(kind=sp),dimension(:,:),allocatable :: levs_fullmet_sp ! This hold each of the numbered level coordinates:
                                                                  !    i.e. isobaric, isobaric1, isobaric2, but also
                                                                  !    height_above_ground, depth_below_surface_layer, etc.
                                                                  !    p_fullmet_sp,p_fullmet_[Vz,RH]sp are copies of one
                                                                  !    of the slices
      integer,dimension(:), allocatable :: nlevs_fullmet
      integer,dimension(:), allocatable :: levs_code                     ! code indicating how to map to the GPH grid
                                                                         !   0 = no mapping (not a pressure coordinate)
                                                                         !   1 = one-to-one mapping (U,V)
                                                                         !   2 = upper truncation (missing upper levels)
                                                                         !   3 = interpolation (missing mid-levels)
                                                                         !   4 = more levels than GPH grid
#endif

      logical       :: IsLatLon_MetGrid
      logical       :: IsGlobal_MetGrid     = .false.  ! Not all Lon/Lat grids are periodic
      logical       :: IsLatLon_CompGrid
      logical       :: IsPeriodic_CompGrid  = .false.
      logical       :: UseFullMetGrid       = .false.  ! This is the special case where the comp grid
                                                       ! equals the Met grid

      integer       :: nx_submet ! length of x or lon of sub-grid
      integer       :: ny_submet ! length of y or lat of sub-grid
#ifdef USEPOINTERS
      real(kind=sp),dimension(:), pointer:: x_submet_sp => null() ! x-coordinates of met sub-grid
      real(kind=sp),dimension(:), pointer:: y_submet_sp => null() ! y-coordinates of met sub-grid
      real(kind=sp),dimension(:), pointer:: z_approx    => null() ! approximate altidue from STD Atmos and press (in km)
#else
      real(kind=sp),dimension(:), allocatable :: x_submet_sp ! x-coordinates of met sub-grid
      real(kind=sp),dimension(:), allocatable :: y_submet_sp ! y-coordinates of met sub-grid
      real(kind=sp),dimension(:), allocatable :: z_approx ! approximate altitude from STD Atmos and press (in km)
#endif

      logical :: x_inverted     = .false.
      logical :: y_inverted     = .false. ! Some LatLon grids start at the North Pole and increment down
      logical :: z_inverted     = .false. ! Some grids give top pressure first
      logical :: y_pad_North    = .false. ! Some computational grids will require values above the top met point
      logical :: y_pad_South    = .false. !   
      logical :: isGridRelative = .true.  ! Most windfiles, whether Lat/Lon or projected, give
                                          ! velocities relative to the grid of the file.  Some (NARR)
                                          ! give velocities relative to earth coordinates and need to
                                          ! be rotated
      real(kind=dp),dimension(:,:) ,allocatable :: theta_Met  ! Earth to grid
      real(kind=dp),dimension(:,:) ,allocatable :: theta_Comp ! Earth to grid

      ! Met copies of projection variables, used for proj call on Met Grid
      integer      :: Met_iprojflag
      real(kind=8) :: Met_Re
      real(kind=8) :: Met_k0
      real(kind=8) :: Met_phi0    != 90.0_ip        ! latitude of projection point
      real(kind=8) :: Met_lam0 != -135.0_ip   ! longitude of projection point
      real(kind=8) :: Met_lam1,Met_phi1
      real(kind=8) :: Met_lam2,Met_phi2

      integer      :: Comp_iprojflag
      real(kind=8) :: Comp_Re
      real(kind=8) :: Comp_k0
      real(kind=8) :: Comp_phi0    != 90.0_ip        ! latitude of projection point
      real(kind=8) :: Comp_lam0 != -135.0_ip   ! longitude of projection point
      real(kind=8) :: Comp_lam1,Comp_phi1
      real(kind=8) :: Comp_lam2,Comp_phi2

      integer      :: Map_Case

      ! Some geometry terms
#ifdef USEPOINTERS
      real(kind=sp),dimension(:,:)     ,pointer :: rdphi_MetP_sp           => null()
      real(kind=sp),dimension(:,:,:)   ,pointer :: rdlambda_MetP_sp        => null()
      real(kind=sp),dimension(:)       ,pointer :: MR_dx_met, MR_dx_submet => null()
      real(kind=sp),dimension(:)       ,pointer :: MR_dy_met, MR_dy_submet => null()
#else
      real(kind=sp),dimension(:,:)     ,allocatable,private :: rdphi_MetP_sp
      real(kind=sp),dimension(:,:,:)   ,allocatable,private :: rdlambda_MetP_sp
      real(kind=sp),dimension(:)       ,allocatable         :: MR_dx_met, MR_dx_submet
      real(kind=sp),dimension(:)       ,allocatable         :: MR_dy_met, MR_dy_submet
#endif

      logical       :: IsRegular_MetGrid                ! True if the grid-spacing is uniform
      real(kind=sp) :: dx_met_const, dy_met_const

        ! There are some computational grid variables we might need, so make local copies
      integer       :: MR_BaseYear            = 1900    ! This should be reset in calling program
      logical       :: MR_useLeap             = .true.  ! This too
      integer       :: MR_Comp_StartYear
      integer       :: MR_Comp_StartMonth
      integer       :: MR_Comp_StartDay
      real(kind=dp) :: MR_Comp_StartHour      = 0.0_dp  ! Note that these must be double-precision to
      real(kind=dp) :: MR_Comp_Time_in_hours  = 0.0_dp  ! be passed correctly to HoursSince
      integer       :: nx_comp
      integer       :: ny_comp
      integer       :: nz_comp
      real(kind=sp) :: dx_comp,dy_comp
      real(kind=sp) :: MaxZ_comp_sp
#ifdef USEPOINTERS
      real(kind=sp),dimension(:),    pointer :: x_comp_sp => null() ! x-coordinates of computational grid
      real(kind=sp),dimension(:),    pointer :: y_comp_sp => null() ! y-coordinates of computational grid
      real(kind=sp),dimension(:),    pointer :: z_comp_sp => null() ! z-coordinates of computational grid
      real(kind=sp),dimension(:,:),  pointer :: CompPoint_X_on_Met_sp   => null() ! x-coord (on Met grid) of comp point
      real(kind=sp),dimension(:,:),  pointer :: CompPoint_Y_on_Met_sp   => null() ! y-coord (on Met grid) of comp point
      integer      ,dimension(:,:,:),pointer :: CompPoint_on_subMet_idx => null() ! index on met sub-grid of comp point
      real(kind=sp),dimension(:,:,:),pointer :: bilin_map_wgt => null()
#else
      real(kind=sp),dimension(:),    allocatable :: x_comp_sp ! x-coordinates of computational grid
      real(kind=sp),dimension(:),    allocatable :: y_comp_sp ! y-coordinates of computational grid
      real(kind=sp),dimension(:),    allocatable :: z_comp_sp ! z-coordinates of computational grid
      real(kind=sp),dimension(:,:),  allocatable :: CompPoint_X_on_Met_sp   ! x-coord (on Met grid) of comp point
      real(kind=sp),dimension(:,:),  allocatable :: CompPoint_Y_on_Met_sp   ! y-coord (on Met grid) of comp point
      integer      ,dimension(:,:,:),allocatable :: CompPoint_on_subMet_idx ! index on met sub-grid of comp point
      real(kind=sp),dimension(:,:,:),allocatable :: bilin_map_wgt
#endif
      ! Here are a few variables needed for sigma-altitude coordinates
      logical        :: MR_use_SigmaAlt   = .false.
      integer        :: MR_ZScaling_ID    = 0  ! = 0 for no scaling (i.e. s = z)
                                               ! = 1 for altitude shifting (s=z=zsurf)
                                               ! = 2 for sigma-altitude (s=(z-surf)/(top-surf))
      real(kind=sp)  :: MR_ztop
#ifdef USEPOINTERS
      real(kind=sp),dimension(:),    pointer :: s_comp_sp => null() ! s-coordinates (scaled z) of computational grid
      real(kind=sp),dimension(:,:),  pointer :: MR_zsurf  => null() ! surface elevation in km
      real(kind=sp),dimension(:,:),  pointer :: MR_jacob  => null() ! Jacobian of trans. = MR_ztop-MR_zsurf
#else
      real(kind=sp),dimension(:),    allocatable :: s_comp_sp ! s-coordinates (scaled z) of computational grid
      real(kind=sp),dimension(:,:),  allocatable :: MR_zsurf  ! surface elevation in km
      real(kind=sp),dimension(:,:),  allocatable :: MR_jacob  ! Jacobian of trans. = MR_ztop-MR_zsurf
#endif

      logical                        :: FoundFillVAttr = .false.
      real(kind=sp)                  :: fill_value_sp  = -9999.0_sp  ! Initialize to the commonly used fill value
      character(len=30),dimension(9) :: Met_dim_names      ! name of dimension
      logical          ,dimension(9) :: Met_dim_IsAvailable
      real(kind=sp)    ,dimension(9) :: Met_dim_fac  = 1.0_sp
      ! Here is the list of variables that can be read.  Each iwindformat will
      ! have just a sub-set available with specific names.  For now, allocate
      ! space for 50 variable names
        ! Mechanical / State variables
        !   1 = Geopotential Height
        !   2 = Vx                               : m/s
        !   3 = Vy                               : m/s
        !   4 = Vz                               : m/s
        !   5 = Temperature                      : K
        !   6 = pressure (when 3d for eta grids) : Pa
        !   7 = PVV (Pressure Vertical Vel.)     : Pa/s
        ! Surface
        !  10 = Planetary Boundary Layer Height
        !  11 = U @ 10m
        !  12 = V @ 10m
        !  13 = Friction velocity
        !  15 = Snow cover
        !  16 = Soil moisture
        !  17 = Surface Roughness
        !  18 = Wind gust speed
        !  19 = surface temperature
        ! Atmospheric Structure
        !  20 = pressure at lower cloud base
        !  21 = pressure at lower cloud top
        !  22 = temperature at lower cloud top
        !  23 = Total Cloud cover
        !  24 = Cloud cover (low)
        !  25 = Cloud cover (convective)
        ! Moisture
        !  30 = Rel. Hum
        !  31 = QV (specific humidity)
        !  32 = QL (liquid)
        !  33 = QI (ice)
        ! Precipitation
        !  40 = Categorical rain
        !  41 = Categorical snow
        !  42 = Categorical frozen rain
        !  43 = Categorical ice
        !  44 = Precipitation rate large-scale (liquid)
        !  45 = Precipitation rate convective  (liquid)
        !  46 = Precipitation rate large-scale (ice)
        !  47 = Precipitation rate convective  (ice)
      logical          ,dimension(MR_MAXVARS)   :: Met_var_IsAvailable      ! true if iwf contains the var
      character(len=71),dimension(MR_MAXVARS)   :: Met_var_NC_names          ! name in the file
      character(len=71),dimension(MR_MAXVARS)   :: Met_var_GRIB_names        ! name in the file
      character(len=5) ,dimension(MR_MAXVARS)   :: Met_var_WMO_names        ! WMO version of the name
      integer          ,dimension(MR_MAXVARS)   :: Met_var_ndim             ! 
      integer          ,dimension(MR_MAXVARS)   :: Met_var_zdim_idx         ! The index of this coordinate (used in Met_var_nlevs) 
      integer          ,dimension(MR_MAXVARS)   :: Met_var_zdim_ncid        ! The dimID of the dimension in the nc file
      integer                                   :: nlev_coords_detected = 0
      integer          ,dimension(MR_MAXVARS,4) :: Met_var_GRIB2_DPcPnSt    ! Grib2 files have variables identified by
                                                                            ! discpln,param_cat,param_num,surf_class
      integer          ,dimension(MR_MAXVARS)   :: Met_var_GRIB1_Param      !   indicatorOfParameter
      integer          ,dimension(MR_MAXVARS)   :: Met_var_GRIB1_Table      !   table2Version
      character(len=3) ,dimension(MR_MAXVARS)   :: Met_var_GRIB1_St         ! level type (pl, src, 116 etc)
      integer                                   :: MR_GRIB_Version  = 0

      !logical          ,dimension(MR_MAXVARS)   :: Met_var_IsFloat  ! true if kind=4 otherwise false
      real(kind=sp)    ,dimension(MR_MAXVARS)   :: Met_var_conversion_factor

      integer          ,dimension(MR_MAXVARS)   :: Met_var_nlevs

        ! Variables needed by netcdf reader
      real(kind=sp) :: iwf25_scale_facs(MR_MAXVARS)
      real(kind=sp) :: iwf25_offsets(MR_MAXVARS)
      real(kind=sp) :: x_in_iwf25_sp(192)
      real(kind=sp) :: y_in_iwf25_sp(94)
        ! Here is the mapping for bilinear weighting coefficients (amap) and
        ! indices (imap) from the 1.875-deg 2d grid to the 2.5-deg 
#ifdef USEPOINTERS
      real(kind=sp),dimension(:,:,:)   ,pointer :: amap_iwf25 => null()
      integer      ,dimension(:,:,:)   ,pointer :: imap_iwf25 => null()
#else
      real(kind=sp),dimension(:,:,:)   ,allocatable :: amap_iwf25
      integer      ,dimension(:,:,:)   ,allocatable :: imap_iwf25
#endif
      integer(kind=sp),dimension(:,:,:),allocatable :: tmpsurf2d_short

      integer :: istart, iend
      integer :: jstart, jend
      integer :: ilhalf_fm_l, ilhalf_fm_r
      integer :: irhalf_fm_l, irhalf_fm_r
      integer :: ilhalf_nx, irhalf_nx
      logical :: wrapgrid

      real(kind=sp)   ,dimension(:)       ,allocatable :: temp1d_sp
      real(kind=sp)   ,dimension(:,:,:)   ,allocatable :: temp2d_sp
      real(kind=sp)   ,dimension(:,:,:,:) ,allocatable :: temp3d_sp
      integer(kind=sp),dimension(:,:,:)   ,allocatable :: temp2d_int
      integer(kind=sp),dimension(:,:,:)   ,allocatable :: temp2d_short
      integer(kind=sp),dimension(:,:,:,:) ,allocatable :: temp3d_short

      real(kind=4),dimension(:,:),allocatable :: Met_Proj_lat
      real(kind=4),dimension(:,:),allocatable :: Met_Proj_lon

      ! Status variables for error-checking
      logical :: Check_prereq_conditions            = .true.
      logical :: CALLED_MR_Allocate_FullMetFileList = .false.
      logical :: CALLED_MR_Read_Met_DimVars         = .false.
      logical :: CALLED_MR_Set_CompProjection       = .false.
      logical :: CALLED_MR_Initialize_Met_Grids     = .false.
      logical :: CALLED_MR_Set_Met_Times            = .false.

      contains

!##############################################################################
!
!     MR_Reset_Memory
!
!     This subroutine reinitializes MetReader by deallocating all MR variables.
!     This is useful if a program needs to use multiple types of wind files.
!
!##############################################################################

      subroutine MR_Reset_Memory

      implicit none

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"-------------------------------------------------------"
        write(MR_global_production,*)"-------- Resetting all MetReader Memory ---------------"
        write(MR_global_production,*)"-------------------------------------------------------"
      endif

#ifndef USEPOINTERS
       if(allocated(MR_windfile_starthour         ))deallocate(MR_windfile_starthour)
       if(allocated(MR_windfile_stephour          ))deallocate(MR_windfile_stephour)
#endif       
       if(allocated(MR_windfiles_nt_fullmet       ))deallocate(MR_windfiles_nt_fullmet)
       if(allocated(MR_windfiles_Have_GRIB_index  ))deallocate(MR_windfiles_Have_GRIB_index)
       if(allocated(MR_windfiles_GRIB_index       ))deallocate(MR_windfiles_GRIB_index)
       if(allocated(MR_MetStep_File               ))deallocate(MR_MetStep_File)
       if(allocated(MR_MetStep_findex             ))deallocate(MR_MetStep_findex)
       if(allocated(MR_MetStep_tindex             ))deallocate(MR_MetStep_tindex)
       if(allocated(MR_MetStep_Hour_since_baseyear))deallocate(MR_MetStep_Hour_since_baseyear)
       if(allocated(MR_MetStep_Interval           ))deallocate(MR_MetStep_Interval)
       if(allocated(MR_MetStep_year               ))deallocate(MR_MetStep_year)
       if(allocated(MR_MetStep_month              ))deallocate(MR_MetStep_month)
       if(allocated(MR_MetStep_day                ))deallocate(MR_MetStep_day)
       if(allocated(MR_MetStep_DOY                ))deallocate(MR_MetStep_DOY)
       if(allocated(MR_MetStep_Hour_Of_Day        ))deallocate(MR_MetStep_Hour_Of_Day)
       if(allocated(MR_iwind5_year                ))deallocate(MR_iwind5_year)

#ifdef USEPOINTERS
       if(associated(MR_windfiles                  ))deallocate(MR_windfiles)
       if(associated(MR_dum2d_met_int              ))deallocate(MR_dum2d_met_int)
       if(associated(MR_dum2d_met                  ))deallocate(MR_dum2d_met)
       if(associated(MR_dum3d_metP                 ))deallocate(MR_dum3d_metP)
       if(associated(MR_dum3d2_metP                ))deallocate(MR_dum3d2_metP)
       if(associated(MR_geoH_metP_last             ))deallocate(MR_geoH_metP_last)
       if(associated(MR_geoH_metP_next             ))deallocate(MR_geoH_metP_next)
       if(associated(MR_vx_metP_last               ))deallocate(MR_vx_metP_last)
       if(associated(MR_vx_metP_next               ))deallocate(MR_vx_metP_next)
       if(associated(MR_vy_metP_last               ))deallocate(MR_vy_metP_last)
       if(associated(MR_vy_metP_next               ))deallocate(MR_vy_metP_next)
       if(associated(MR_dum3d_metH                 ))deallocate(MR_dum3d_metH)
       if(associated(MR_dum2d_comp_int             ))deallocate(MR_dum2d_comp_int)
       if(associated(MR_dum2d_comp                 ))deallocate(MR_dum2d_comp)
       if(associated(MR_dum3d_compP                ))deallocate(MR_dum3d_compP)
       if(associated(MR_dum3d_compP_2              ))deallocate(MR_dum3d_compP_2)
       if(associated(MR_dum3d_compH                ))deallocate(MR_dum3d_compH)
       if(associated(MR_dum3d_compH_2              ))deallocate(MR_dum3d_compH_2)
       if(associated(x_fullmet_sp                  ))deallocate(x_fullmet_sp)
       if(associated(y_fullmet_sp                  ))deallocate(y_fullmet_sp)
       if(associated(p_fullmet_sp                  ))deallocate(p_fullmet_sp)
       if(associated(x_submet_sp                   ))deallocate(x_submet_sp)
       if(associated(y_submet_sp                   ))deallocate(y_submet_sp)
       if(associated(levs_fullmet_sp               ))deallocate(levs_fullmet_sp)
       if(associated(nlevs_fullmet                 ))deallocate(nlevs_fullmet)
       if(associated(levs_code                     ))deallocate(levs_code)
       if(associated(z_approx                      ))deallocate(z_approx)
       if(associated(rdphi_MetP_sp                 ))deallocate(rdphi_MetP_sp)
       if(associated(rdlambda_MetP_sp              ))deallocate(rdlambda_MetP_sp)
       if(associated(MR_dx_met                     ))deallocate(MR_dx_met)
       if(associated(MR_dx_submet                  ))deallocate(MR_dx_submet)
       if(associated(MR_dy_met                     ))deallocate(MR_dy_met)
       if(associated(MR_dy_submet                  ))deallocate(MR_dy_submet)
       if(associated(x_comp_sp                     ))deallocate(x_comp_sp)
       if(associated(y_comp_sp                     ))deallocate(y_comp_sp)
       if(associated(z_comp_sp                     ))deallocate(z_comp_sp)
       if(associated(CompPoint_X_on_Met_sp         ))deallocate(CompPoint_X_on_Met_sp)
       if(associated(CompPoint_Y_on_Met_sp         ))deallocate(CompPoint_Y_on_Met_sp)
       if(associated(CompPoint_on_subMet_idx       ))deallocate(CompPoint_on_subMet_idx)
       if(associated(bilin_map_wgt                 ))deallocate(bilin_map_wgt)
       if(associated(amap_iwf25                    ))deallocate(amap_iwf25)
       if(associated(imap_iwf25                    ))deallocate(imap_iwf25)
       if(associated(s_comp_sp                     ))deallocate(s_comp_sp)
       if(associated(MR_zsurf                      ))deallocate(MR_zsurf)
       if(associated(MR_jacob                      ))deallocate(MR_jacob)
#else
       if(allocated(MR_windfiles                  ))deallocate(MR_windfiles)
       if(allocated(MR_dum2d_met_int              ))deallocate(MR_dum2d_met_int)
       if(allocated(MR_dum2d_met                  ))deallocate(MR_dum2d_met)
       if(allocated(MR_dum3d_metP                 ))deallocate(MR_dum3d_metP)
       if(allocated(MR_dum3d2_metP                ))deallocate(MR_dum3d2_metP)
       if(allocated(MR_geoH_metP_last             ))deallocate(MR_geoH_metP_last)
       if(allocated(MR_geoH_metP_next             ))deallocate(MR_geoH_metP_next)
       if(allocated(MR_vx_metP_last               ))deallocate(MR_vx_metP_last)
       if(allocated(MR_vx_metP_next               ))deallocate(MR_vx_metP_next)
       if(allocated(MR_vy_metP_last               ))deallocate(MR_vy_metP_last)
       if(allocated(MR_vy_metP_next               ))deallocate(MR_vy_metP_next)
       if(allocated(MR_dum3d_metH                 ))deallocate(MR_dum3d_metH)
       if(allocated(MR_dum2d_comp_int             ))deallocate(MR_dum2d_comp_int)
       if(allocated(MR_dum2d_comp                 ))deallocate(MR_dum2d_comp)
       if(allocated(MR_dum3d_compP                ))deallocate(MR_dum3d_compP)
       if(allocated(MR_dum3d_compP_2              ))deallocate(MR_dum3d_compP_2)
       if(allocated(MR_dum3d_compH                ))deallocate(MR_dum3d_compH)
       if(allocated(MR_dum3d_compH_2              ))deallocate(MR_dum3d_compH_2)
       if(allocated(x_fullmet_sp                  ))deallocate(x_fullmet_sp)
       if(allocated(y_fullmet_sp                  ))deallocate(y_fullmet_sp)
       if(allocated(p_fullmet_sp                  ))deallocate(p_fullmet_sp)
       if(allocated(x_submet_sp                   ))deallocate(x_submet_sp)
       if(allocated(y_submet_sp                   ))deallocate(y_submet_sp)
       if(allocated(levs_fullmet_sp               ))deallocate(levs_fullmet_sp)
       if(allocated(nlevs_fullmet                 ))deallocate(nlevs_fullmet)
       if(allocated(levs_code                     ))deallocate(levs_code)
       if(allocated(z_approx                      ))deallocate(z_approx)
       if(allocated(rdphi_MetP_sp                 ))deallocate(rdphi_MetP_sp)
       if(allocated(rdlambda_MetP_sp              ))deallocate(rdlambda_MetP_sp)
       if(allocated(MR_dx_met                     ))deallocate(MR_dx_met)
       if(allocated(MR_dx_submet                  ))deallocate(MR_dx_submet)
       if(allocated(MR_dy_met                     ))deallocate(MR_dy_met)
       if(allocated(MR_dy_submet                  ))deallocate(MR_dy_submet)
       if(allocated(x_comp_sp                     ))deallocate(x_comp_sp)
       if(allocated(y_comp_sp                     ))deallocate(y_comp_sp)
       if(allocated(z_comp_sp                     ))deallocate(z_comp_sp)
       if(allocated(CompPoint_X_on_Met_sp         ))deallocate(CompPoint_X_on_Met_sp)
       if(allocated(CompPoint_Y_on_Met_sp         ))deallocate(CompPoint_Y_on_Met_sp)
       if(allocated(CompPoint_on_subMet_idx       ))deallocate(CompPoint_on_subMet_idx)
       if(allocated(bilin_map_wgt                 ))deallocate(bilin_map_wgt)
       if(allocated(amap_iwf25                    ))deallocate(amap_iwf25)
       if(allocated(imap_iwf25                    ))deallocate(imap_iwf25)
       if(allocated(s_comp_sp                     ))deallocate(s_comp_sp)
       if(allocated(MR_zsurf                      ))deallocate(MR_zsurf)
       if(allocated(MR_jacob                      ))deallocate(MR_jacob)
#endif
       if(allocated(MR_u_ER_metP                  ))deallocate(MR_u_ER_metP)
       if(allocated(MR_v_ER_metP                  ))deallocate(MR_v_ER_metP)
       if(allocated(theta_Met                     ))deallocate(theta_Met)
       if(allocated(theta_Comp                    ))deallocate(theta_Comp)
       if(allocated(tmpsurf2d_short               ))deallocate(tmpsurf2d_short)
       if(allocated(temp1d_sp                     ))deallocate(temp1d_sp)
       if(allocated(temp2d_sp                     ))deallocate(temp2d_sp)
       if(allocated(temp3d_sp                     ))deallocate(temp3d_sp)
       if(allocated(temp2d_int                    ))deallocate(temp2d_int)
       if(allocated(temp2d_short                  ))deallocate(temp2d_short)
       if(allocated(temp3d_short                  ))deallocate(temp3d_short)
       if(allocated(Met_Proj_lat                  ))deallocate(Met_Proj_lat)
       if(allocated(Met_Proj_lon                  ))deallocate(Met_Proj_lon)

       if(allocated(MR_SndVars_metP               ))deallocate(MR_SndVars_metP)
       if(allocated(MR_SndVarsID                  ))deallocate(MR_SndVarsID)
       if(allocated(MR_Snd_np_fullmet             ))deallocate(MR_Snd_np_fullmet)
       if(allocated(MR_Snd2Comp_tri_map_wgt       ))deallocate(MR_Snd2Comp_tri_map_wgt)
       if(allocated(MR_Snd2Comp_tri_map_idx       ))deallocate(MR_Snd2Comp_tri_map_idx)

       nlev_coords_detected = 0

      end subroutine MR_Reset_Memory

!##############################################################################
!
!     MR_Allocate_FullMetFileList
!
!     This subroutine allocates the list of windfiles and does some
!     error-checking based on iwind and iwindfiles.
!
!     From the calling program, this is called once the information about
!     the NWP files is available (e.g. the iwind, iwindformat, grid ID, data
!     format, and number of windfiles.
!
!     Variables allocated:
!        MR_windfiles(MR_iwindfiles)
!        MR_windfiles_nt_fullmet(MR_iwindfiles)
!       and possibly:
!        MR_windfiles_Have_GRIB_index(MR_iwindfiles)
!        MR_windfiles_GRIB_index(MR_iwindfiles)
!
!     The next step is for the calling program to populate MR_windfiles
!
!##############################################################################

      subroutine MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)

      implicit none

      integer,intent(in)      :: iw
      integer,intent(in)      :: iwf
      integer,intent(in)      :: igrid
      integer,intent(in)      :: idf
      integer,intent(in)      :: iwfiles

      integer :: i

      MR_iwind              = iw
      MR_iwindformat        = iwf
      MR_iGridCode          = igrid
      MR_idataFormat        = idf
      MR_iwindfiles         = iwfiles

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Allocate_FullMetFileList                           ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif
      if ((MR_iwind.ne.1).and.(MR_iwind.ne.2).and. &
          (MR_iwind.ne.3).and.(MR_iwind.ne.4).and. &
          (MR_iwind.ne.5)) then
        write(MR_global_error,*)'MR_iwind must be between 1 and 5. Program stopped'
        write(MR_global_error,*)' MR_iwind = ',MR_iwind
        write(MR_global_error,*)' MR_IWIND OPTIONS:'
        write(MR_global_error,*)' MR_iwind = 1 read from a 1-D wind sounding'
        write(MR_global_error,*)'            2 read from 3D gridded ASCII files'
        write(MR_global_error,*)'            3 read from a single, multistep file'
        write(MR_global_error,*)'            4 read from multiple files'
        write(MR_global_error,*)'            5 read variables from separate files'
        stop 1
      endif

      ! Initialize the dimension and variable arrays.  Select slots in these arrays will be
      ! overwritten from the calls in the case block below
      Met_dim_IsAvailable=.false.
      Met_var_IsAvailable=.false.
      Met_var_IsAvailable(1:MR_MAXVARS)       = .false.
      Met_var_zdim_idx(1:MR_MAXVARS)          = 0
      Met_var_zdim_ncid(1:MR_MAXVARS)         = 0
      Met_var_GRIB2_DPcPnSt(1:MR_MAXVARS,1:4) = 0
      Met_var_GRIB1_Param(1:MR_MAXVARS)       = 0
      Met_var_GRIB1_St(1:MR_MAXVARS)          = ""
      Met_var_conversion_factor(1:MR_MAXVARS) = 1.0_sp
      Met_var_nlevs(1:MR_MAXVARS)             = 0
      isGridRelative = .true.

      !--------------------------------
      ! Mechanical / State variables
      !--------------------------------
      !   Note: All default names are that assigned by netcdf-java in the grib-to-nc conversion
      !         Command used :: java -Xmx2048m -classpath ~/ncj/netcdfAll-4.5.jar ucar.nc2.dataset.NetcdfDataset \
      !                          -in ${GribFile} -out ${NCFile} -IsLargeFile
      ! 
      !         Alternatively, one could use ncl_convert2nc
        !  Geopotential Height  (m^2/s^2)
      Met_var_NC_names(1)          = "Geopotential_height_isobaric"
      Met_var_GRIB_names(1)        = "gh"
      Met_var_WMO_names(1)         = "HGT"
      Met_var_GRIB2_DPcPnSt(1,1:4) = (/0, 3, 5, 100/)
      Met_var_GRIB1_Param(1)       = 7
      Met_var_GRIB1_St(1)          = "pl"
      Met_var_ndim(1)              = 4
        ! Velocity component in x (or E) direction  (m/s)
      Met_var_NC_names(2)          = "u-component_of_wind_isobaric"
      Met_var_GRIB_names(2)        = "u"
      Met_var_WMO_names(2)         = "UGRD"
      Met_var_GRIB2_DPcPnSt(2,1:4) = (/0, 2, 2, 100/)
      Met_var_GRIB1_Param(2)       = 33
      Met_var_GRIB1_St(2)          = "pl"
      Met_var_ndim(2)              = 4
        ! Velocity component in y (or N) direction (m/s)
      Met_var_NC_names(3)          = "v-component_of_wind_isobaric"
      Met_var_GRIB_names(3)        = "v"
      Met_var_WMO_names(3)         = "VGRD"
      Met_var_GRIB2_DPcPnSt(3,1:4) = (/0, 2, 3, 100/)
      Met_var_GRIB1_Param(3)       = 34
      Met_var_GRIB1_St(3)          = "pl"
      Met_var_ndim(3)              = 4
        ! Velocity component in z direction  (Pa/s)
        !  This actually returns Vz with VPP*dp/dz with a FD approx of dp/dz
      Met_var_NC_names(4)          = "Vertical_velocity_pressure_isobaric"
      Met_var_GRIB_names(4)        = "w"
      Met_var_WMO_names(4)         = "VVEL"
      Met_var_GRIB2_DPcPnSt(4,1:4) = (/0, 2, 8, 100/)
      Met_var_GRIB1_Param(4)       = 39
      Met_var_GRIB1_St(4)          = "pl"
      Met_var_ndim(4)              = 4
        ! Temperature  (K)
      Met_var_NC_names(5)          = "Temperature_isobaric"
      Met_var_GRIB_names(5)        = "t"
      Met_var_WMO_names(5)         = "TMP"
      Met_var_GRIB2_DPcPnSt(5,1:4) = (/0, 0, 0, 100/)
      Met_var_GRIB1_Param(5)       = 11
      Met_var_GRIB1_St(5)          = "pl"
      Met_var_ndim(5)              = 4
        ! Velocity component in z direction  (Pa/s)
        !  This returns the true vert. vel. pres.
      Met_var_NC_names(7)          = "Vertical_velocity_pressure_isobaric"
      Met_var_GRIB_names(7)        = "w"
      Met_var_WMO_names(7)         = "VVEL"
      Met_var_GRIB2_DPcPnSt(7,1:4) = (/0, 2, 8, 100/)
      Met_var_GRIB1_Param(7)       = 39
      Met_var_GRIB1_St(7)          = "pl"
      Met_var_ndim(7)              = 4

      !--------------------------------
      ! Surface
      !--------------------------------
        ! Height of planetary boundary layer  (m)
      Met_var_NC_names(10)          = "Planetary_Boundary_Layer_Height_surface"
      Met_var_GRIB_names(10)        = "hpbl"
      Met_var_WMO_names(10)         = "HPBL"
      Met_var_GRIB2_DPcPnSt(10,1:4) = (/0, 3, 196, 1/)
      Met_var_GRIB1_Param(10)       = 221
      Met_var_GRIB1_St(10)          = "1"
      Met_var_ndim(10)              = 3
        ! Velocity component in x (or E) direction at 10 m above ground surface  (m/s)
      Met_var_NC_names(11)          = "u-component_of_wind_height_above_ground"
      Met_var_GRIB_names(11)        = "u"
      Met_var_WMO_names(11)         = "UGRD"
      Met_var_GRIB2_DPcPnSt(11,1:4) = (/0, 2, 2, 103/)
      Met_var_GRIB1_Param(11)       = 33
      Met_var_GRIB1_St(11)          = "105"
      Met_var_ndim(11)              = 4
        ! Velocity component in y (or N) direction at 10 m above ground surface  (m/s)
      Met_var_NC_names(12)          = "v-component_of_wind_height_above_ground"
      Met_var_GRIB_names(12)        = "v"
      Met_var_WMO_names(12)         = "VGRD"
      Met_var_GRIB2_DPcPnSt(12,1:4) = (/0, 2, 3, 103/)
      Met_var_GRIB1_Param(12)       = 34
      Met_var_GRIB1_St(12)          = "105"
      Met_var_ndim(12)              = 4
        ! Friction velocity  (m/s)
      Met_var_NC_names(13)          = "Frictional_Velocity_surface"
      Met_var_GRIB_names(13)        = "fricv"
      Met_var_WMO_names(13)         = "FRICV"
      Met_var_GRIB2_DPcPnSt(13,1:4) = (/0, 2, 197, 1/)
      Met_var_GRIB1_Param(13)       = 253
      Met_var_GRIB1_St(13)          = "1"
      Met_var_ndim(13)              = 3
        ! Snow depth  (m)
      Met_var_NC_names(15)          = "Snow_depth_surface"
      Met_var_GRIB_names(15)        = "sd"
      Met_var_WMO_names(15)         = "SNOD"
      Met_var_GRIB2_DPcPnSt(15,1:4) = (/0, 1, 11, 1/)
      Met_var_GRIB1_Param(15)       = 66
      Met_var_GRIB1_St(15)          = "1"
      Met_var_ndim(15)              = 3
        ! Soil Moisture  (fraction)
      Met_var_NC_names(16)          = "Volumetric_Soil_Moisture_Content_depth_below_surface_layer"
      Met_var_GRIB_names(16)        = "soilw"
      Met_var_WMO_names(16)         = "SOILW"
      Met_var_GRIB2_DPcPnSt(16,1:4) = (/2, 0, 192, 106/)
      Met_var_GRIB1_Param(16)       = 0
      Met_var_GRIB1_St(16)          = ""
      Met_var_ndim(16)              = 4
        ! Surface roughness  (m)
      Met_var_NC_names(17)          = "Surface_roughness_surface"
      Met_var_GRIB_names(17)        = "sr"
      Met_var_WMO_names(17)         = "SFCR"
      Met_var_GRIB2_DPcPnSt(17,1:4) = (/2, 0, 1, 1/)
      Met_var_GRIB1_Param(17)       = 0
      Met_var_GRIB1_St(17)          = ""
      Met_var_ndim(17)              = 3
        ! Wind gust speed  (m/s)
      Met_var_NC_names(18)          = "Wind_speed_gust_surface"
      Met_var_GRIB_names(18)        = "gust"
      Met_var_WMO_names(18)         = "GUST"
      Met_var_GRIB2_DPcPnSt(18,1:4) = (/0, 2, 22, 1/)
      Met_var_GRIB1_Param(18)       = 0
      Met_var_GRIB1_St(18)          = ""
      Met_var_ndim(18)              = 3
        ! Surface temperature  (K)
      Met_var_NC_names(19)          = "Temperature_surface"
      Met_var_WMO_names(19)         = ""
      Met_var_GRIB2_DPcPnSt(19,1:4) = (/0, 0, 0, 0/)
      Met_var_GRIB1_Param(19)       = 11
      Met_var_GRIB1_St(19)          = "1"
      Met_var_ndim(19)              = 3
      !--------------------------------
      ! Atmospheric Structure
      !--------------------------------
        ! Pressure at base of lower cloud level  (Pa)
      Met_var_NC_names(20)          = "Pressure_cloud_base"
      Met_var_GRIB_names(20)        = "pres"
      Met_var_WMO_names(20)         = "PRES"
      Met_var_GRIB2_DPcPnSt(20,1:4) = (/0, 3, 0, 2/)
      Met_var_GRIB1_Param(20)       = 1
      Met_var_GRIB1_St(20)          = "2"
      Met_var_ndim(20)              = 3
        ! Pressure at top of lower cloud level  (Pa)
      Met_var_NC_names(21)          = "Pressure_cloud_tops"
      Met_var_GRIB_names(21)        = "pres"
      Met_var_WMO_names(21)         = "PRES"
      Met_var_GRIB2_DPcPnSt(21,1:4) = (/0, 3, 0, 3/)
      Met_var_GRIB1_Param(21)       = 1
      Met_var_GRIB1_St(21)          = "3"
      Met_var_ndim(21)              = 3
        ! Temperature at the top of lower cloud level  (K)
      Met_var_NC_names(22)          = "Temperature_cloud_tops"
      Met_var_GRIB_names(22)        = "t"
      Met_var_WMO_names(22)         = "TMP"
      Met_var_GRIB2_DPcPnSt(22,1:4) = (/0, 0, 0, 3/)
      Met_var_GRIB1_Param(22)       = 11
      Met_var_GRIB1_St(22)          = "3"
      Met_var_ndim(22)              = 3
        ! Total cloud cover  (%)
      Met_var_NC_names(23)          = "Total_cloud_cover_entire_atmosphere"
      Met_var_GRIB_names(23)        = "tcc"
      Met_var_WMO_names(23)         = "TCDC"
      Met_var_GRIB2_DPcPnSt(23,1:4) = (/0, 6, 1, 200/)
      Met_var_GRIB1_Param(23)       = 71
      Met_var_GRIB1_St(23)          = "200"
      Met_var_ndim(23)              = 3
        ! Cloud cover of lower cloud level  (%)
      Met_var_NC_names(24)          = "Low_cloud_cover_low_cloud"
      Met_var_GRIB_names(24)        = "lcc"
      Met_var_WMO_names(24)         = "LCDC"
      Met_var_GRIB2_DPcPnSt(24,1:4) = (/0, 6, 3, 214/)
      Met_var_GRIB1_Param(24)       = 0
      Met_var_GRIB1_St(24)          = ""
      Met_var_ndim(24)              = 3
      !--------------------------------
      ! Moisture
      !--------------------------------
        ! Relative humidity  (%)
      Met_var_NC_names(30)          = "Relative_humidity_isobaric"
      Met_var_GRIB_names(30)        = "r"
      Met_var_WMO_names(30)         = "RH"
      Met_var_GRIB2_DPcPnSt(30,1:4) = (/0, 1, 1, 100/)
      Met_var_GRIB1_Param(30)       = 52
      Met_var_GRIB1_St(30)          = "pl" 
      Met_var_ndim(30)              = 4
        ! Specific humidity  (kg/kg)
      Met_var_NC_names(31)          = "Specific_humidity_isobaric"
      Met_var_GRIB_names(31)        = "q"
      Met_var_WMO_names(31)         = "SPFH"
      Met_var_GRIB2_DPcPnSt(31,1:4) = (/0, 1, 0, 100/)
      Met_var_GRIB1_Param(31)       = 51
      Met_var_GRIB1_St(31)          = "pl"
      Met_var_ndim(31)              = 4
        ! Cloud water mixing ratio  (kg/kg)
      Met_var_NC_names(32)          = "Cloud_mixing_ratio_isobaric"
      Met_var_GRIB_names(32)        = "clwmr"
      Met_var_WMO_names(32)         = "CLWMR"
      Met_var_GRIB2_DPcPnSt(32,1:4) = (/0, 1, 22, 100/)
      Met_var_ndim(32)              = 4
        ! Snow mixing ratio  (kg/kg)
      Met_var_NC_names(33)          = "Snow_mixing_ratio_isobaric"
      Met_var_GRIB_names(33)        = "snmr"
      Met_var_WMO_names(33)         = "SNMR"
      Met_var_GRIB2_DPcPnSt(33,1:4) = (/0, 1, 25, 100/)
      Met_var_GRIB1_Param(33)       = 0
      Met_var_GRIB1_St(33)          = ""
      Met_var_ndim(33)              = 4
      !--------------------------------
      ! Precipitation
      !--------------------------------
        ! Categorical rain at ground surface  (0/1 no/yes)
      Met_var_NC_names(40)          = "Categorical_Rain_surface"
      Met_var_GRIB_names(40)        = "crain"
      Met_var_WMO_names(40)         = "CRAIN"
      Met_var_GRIB2_DPcPnSt(40,1:4) = (/0, 1, 192, 1/)
      Met_var_GRIB1_Param(40)       = 140
      Met_var_GRIB1_St(40)          = "1"
      Met_var_ndim(40)              = 3
        ! Categorical snow at ground surface  (0/1 no/yes)
      Met_var_NC_names(41)          = "Categorical_Snow_surface"
      Met_var_GRIB_names(41)        = "csnow"
      Met_var_WMO_names(41)         = "CSNOW"
      Met_var_GRIB2_DPcPnSt(41,1:4) = (/0, 1, 195, 1/)
      Met_var_GRIB1_Param(41)       = 143
      Met_var_GRIB1_St(41)          = "1"
      Met_var_ndim(41)              = 3
        ! Categorical freezing rain at ground surface  (0/1 no/yes)
      Met_var_NC_names(42)          = "Categorical_Freezing_Rain_surface"
      Met_var_GRIB_names(42)        = "cfrzr"
      Met_var_WMO_names(42)         = "CFRZR"
      Met_var_GRIB2_DPcPnSt(42,1:4) = (/0, 1, 193, 1/)
      Met_var_GRIB1_Param(42)       = 141
      Met_var_GRIB1_St(42)          = "1"
      Met_var_ndim(42)              = 3
        ! Categorical ice pellets at ground surface  (0/1 no/yes)
      Met_var_NC_names(43)          = "Categorical_Ice_Pellets_surface"
      Met_var_GRIB_names(43)        = "cicep"
      Met_var_WMO_names(43)         = "CICEP"
      Met_var_GRIB2_DPcPnSt(43,1:4) = (/0, 1, 194, 1/)
      Met_var_GRIB1_Param(43)       = 142
      Met_var_GRIB1_St(43)          = "1"
      Met_var_ndim(43)              = 3
        ! Precipitation rate at surface  (kg/m2/s)
      Met_var_NC_names(44)          = "Precipitation_rate_surface"
      Met_var_GRIB_names(44)        = "prate"
      Met_var_WMO_names(44)         = "PRATE"
      Met_var_GRIB2_DPcPnSt(44,1:4) = (/0, 1, 7, 1/)
      Met_var_GRIB1_Param(44)       = 59
      Met_var_GRIB1_St(44)          = "pl"
      Met_var_ndim(44)              = 3
        ! Convective liquid precipitation rate at surface  (kg/m2/s)
      Met_var_NC_names(45)          = "Precip.rate convective  (liquid)"
      Met_var_WMO_names(45)         = "CPRAT"
      Met_var_GRIB1_Param(45)       = 0
      Met_var_GRIB1_St(45)          = "cprat"
      Met_var_ndim(45)              = 3
        ! Large-scale precipitation rate at surface  (kg/m2/s)
      Met_var_NC_names(46)          = "Precip.rate large-scale (ice)"
      Met_var_WMO_names(46)         = ""
      Met_var_GRIB1_Param(46)       = 0
      Met_var_GRIB1_St(46)          = ""
      Met_var_ndim(46)              = 3
        ! Convective frozen precipitation rate at surface for (kg/m2/s)
      Met_var_NC_names(47)          = "Precip.rate convective  (ice)"
      Met_var_WMO_names(47)         = ""
      Met_var_GRIB1_Param(47)       = 0
      Met_var_GRIB1_St(47)          = ""
      Met_var_ndim(47)              = 3

      if (MR_iwindformat.eq.0) then
          ! Custom format based on template
        MR_Reannalysis = .false.
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "Custom format based on template"

          ! This expects that MR_iwf_template has been filled by the calling program
        call MR_Read_Met_Template

      elseif (MR_iwindformat.eq.1) then
          ! ASCII profile
        MR_Reannalysis = .false.
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "ASCII profile"
      elseif (MR_iwindformat.eq.2) then
          ! Radiosonde data
        MR_Reannalysis = .false.
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "Radiosonde data"
      elseif (MR_iwindformat.eq.3) then 
        ! NARR3D NAM221 32 km North America files
          ! https://rda.ucar.edu/datasets/ds608.0/
          !   merged_AWIP32.2018062000(.nc)
          ! Note that winds are "earth-relative" and must be rotated!
          !   See  http://www.emc.ncep.noaa.gov/mmb/rreanl/faq.html#eta-winds

        if(MR_iversion.eq.-1)MR_iversion = 0 ! the lastest version of NARR data is v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NARR3D NAM221 32 km North America files with Re=6367.47"

        MR_iGridCode = 1221  ! This is almost NAM221, but uses a diff Re
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode) 
        MR_Reannalysis = .true.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "reftime"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric2"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "y"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "x"

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.  ! Geopotential Height
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="u_wind_isobaric"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="v_wind_isobaric"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="Pressure_vertical_velocity_isobaric"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="Temp_isobaric"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="Pressure_vertical_velocity_isobaric"
        ! Surface
        Met_var_IsAvailable(10)=.true.; Met_var_NC_names(10)="Planetary_boundary_layer_height_surface"
        Met_var_IsAvailable(11)=.true.; Met_var_NC_names(11)="u_wind_height_above_ground"
        Met_var_IsAvailable(12)=.true.; Met_var_NC_names(12)="v_wind_height_above_ground"
        Met_var_IsAvailable(13)=.true.; Met_var_NC_names(13)="Surface_friction_velocity_surface"
        Met_var_IsAvailable(15)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        ! Moisture
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.; Met_var_NC_names(32)="Cloud_water_isobaric"
        Met_var_IsAvailable(33)=.true.; Met_var_NC_names(33)="Ice_mixing_ratio_isobaric"
        ! Precipitation
        Met_var_IsAvailable(40)=.true.; Met_var_NC_names(40)="Categorical_rain_yes1_no0_surface"
        Met_var_IsAvailable(41)=.true.; Met_var_NC_names(41)="Categorical_snow_yes1_no0_surface"
        Met_var_IsAvailable(42)=.true.; Met_var_NC_names(42)="Categorical_freezing_rain_yes1_no0_surface"
        Met_var_IsAvailable(43)=.true.; Met_var_NC_names(43)="Categorical_ice_pellets_yes1_no0_surface"
        Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp
      elseif (MR_iwindformat.eq.4) then
        ! NAM Regional North America 221 32 km North America files
          ! http://motherlode.ucar.edu/native/conduit/data/nccf/com/nam/prod/
          !   nam.t00z.awip3200.tm00.grib2

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NAM221 32 km North America files with Re=6371.229"

        MR_iGridCode = 221
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric3"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "y"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "x"

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(16)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(45)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.5) then
        ! NAM216 AK 45km
          ! http://motherlode.ucar.edu/native/conduit/data/nccf/com/nam/prod/
          !  nam.t00z.awipak00.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NAM216 AK 45km"

        MR_iGridCode = 216
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time1"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric1"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "y"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "x"

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.; Met_var_NC_names(10)="Planetary_Boundary_Layer_Height"
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.; Met_var_NC_names(13)="Frictional_Velocity"
        Met_var_IsAvailable(14)=.true.; Met_var_NC_names(16)="Volumetric_Soil_Moisture_Content"
        ! Atmospheric Structure
        Met_var_IsAvailable(23)=.true.; Met_var_NC_names(23)="Total_cloud_cover"
        ! Moisture
        Met_var_IsAvailable(30)=.true.;
        Met_var_IsAvailable(31)=.true.;
        Met_var_IsAvailable(32)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.; Met_var_NC_names(40)="Categorical_Rain"
        Met_var_IsAvailable(41)=.true.; Met_var_NC_names(41)="Categorical_Snow"
        Met_var_IsAvailable(42)=.true.; Met_var_NC_names(42)="Categorical_Freezing_Rain"
        Met_var_IsAvailable(43)=.true.; Met_var_NC_names(43)="Categorical_Ice_Pellets"
        Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="Large_scale_precipitation_non-convective"
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_precipitation"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.6) then
        ! NAM Regional 90 km grid 104
        !  ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod
        !    nam.t00z.grbgrd00.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NAM Regional 90 km grid 104"

        MR_iGridCode = 104
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "y"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "x"

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(16)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(23)=.true.; Met_var_NC_names(23)="Total_cloud_cover_entire_atmosphere_0_Hour_Average"
        Met_var_IsAvailable(24)=.true.; Met_var_NC_names(24)="Convective_cloud_cover_entire_atmosphere_0_Hour_Average"
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="Large_scale_precipitation_non-convective_surface_0_Hour_Accumulation"
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_precipitation_surface_0_Hour_Accumulation"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.7) then
        ! CONUS 212 40km
          !  ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/
          !    nam.t00z.awp21100.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "CONUS 212 40km"

        MR_iGridCode = 212
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time1"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric3"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "y"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "x"

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(23)=.true.; Met_var_NC_names(23)="Total_cloud_cover"
        Met_var_IsAvailable(24)=.true.; Met_var_NC_names(24)="Convective_cloud_cover"
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_Precipitation_Rate_surface"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.8) then
        ! CONUS 218 (12km)
          ! ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/
          !   nam.t00z.awphys00.grb2.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "CONUS 218 (12km)"

        MR_iGridCode = 218
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric1"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "y"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "x"

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.9) then
        ! CONUS 227 (5.08 km)
          !  http://motherlode.ucar.edu/native/conduit/data/nccf/com/nam/prod/
          !    nam.t00z.conusnest.hiresf00.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "CONUS 227 (5.1 km)"

        MR_iGridCode = 227
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time1"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric2"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "y"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "x"

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.
        Met_var_IsAvailable(17)=.true.
        Met_var_IsAvailable(18)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.10)then
        ! NAM 242 11.25 km AK
          !  http://motherlode.ucar.edu/native/conduit/data/nccf/com/nam/prod/
          !    nam.t00z.awipak00.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NAM 242 11.25 km AK"

        MR_iGridCode = 242
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric2"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "y"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "x"

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="Total_precipitation_surface_0_Hour_Accumulation"
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_precipitation_surface_0_Hour_Accumulation"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.11)then
        ! NAM 196 2.5 km HI
          !  ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/
          !    nam.t00z.hawaiinest.hiresf00.tm0

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  " NAM 196 2.5 km HI"

        MR_iGridCode = 196
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric2"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "y"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "x"

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.; Met_var_NC_names(15)="Snow_Cover_surface"
        Met_var_IsAvailable(16)=.true.; Met_var_NC_names(16)="Soil_moisture_content_depth_below_surface_layer"
        Met_var_IsAvailable(17)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.12)then
        ! NAM 198 5.953 km AK
          !  ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/ (this is actually now 91)
          !    nam.t00z.alaskanest.hiresf00.tm00

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NAM 198 5.953 km AK"

        MR_iGridCode = 198
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric2"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "y"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "x"

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.; Met_var_NC_names(16)="Soil_moisture_content_depth_below_surface_layer"
        Met_var_IsAvailable(16)=.true.
        Met_var_IsAvailable(17)=.true.
        Met_var_IsAvailable(18)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_precipitation_surface_0_Hour_Accumulation"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.13)then  ! NAM 91 2.976 km AK
        ! NAM 91 2.976 km AK
          !  ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/ (this is actually now 91)
          !    nam.t00z.alaskanest.hiresf00.tm00
          !
          ! Note: the dimension names given below are those generated by netcdf-java 4.5
          !       acting on the truncated grib files generated by get_nam91.sh which 
          !       uses get_inv.pl to get just the grib layers needed.
          !       This is relavent because the numbering of the isobaric dimensions
          !       appears to be in the order they are processed in the grib file

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NAM 91 2.976 km AK"

        MR_iGridCode = 91
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "y"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "x"

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.
        Met_var_IsAvailable(17)=.true.
        Met_var_IsAvailable(18)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.14) then
        ! CONUS 1227 (3.0 km)
          !  
          !  nam.t00z.conusnest.hiresf00.tm00.grib2.nc

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "CONUS no grid ID (3.0 km)"

        MR_iGridCode = 1227
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time1"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric2"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "y"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "x"

        ! Mechanical / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        Met_var_IsAvailable(13)=.true.
        Met_var_IsAvailable(15)=.true.
        Met_var_IsAvailable(16)=.true.
        Met_var_IsAvailable(17)=.true.
        Met_var_IsAvailable(18)=.true.
        ! Atmospheric Structure
        Met_var_IsAvailable(20)=.true.
        Met_var_IsAvailable(21)=.true.
        Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.
        Met_var_IsAvailable(41)=.true.
        Met_var_IsAvailable(42)=.true.
        Met_var_IsAvailable(43)=.true.
        Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.20)then
        ! GFS 0.5 deg
          !  http://www.nco.ncep.noaa.gov/pmb/products/gfs/
          !  http://motherlode.ucar.edu/native/conduit/data/nccf/com/gfs/prod/
          !    

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "GFS 0.5"

        MR_iGridCode = 4
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric3"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        !Met_var_IsAvailable(10)=.true.; Met_var_GRIB2_DPcPnSt(10,1:4)=(/0, 3, 18, 1/)
        !Met_var_IsAvailable(11)=.true.
        !Met_var_IsAvailable(12)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(31)=.true.
        Met_var_IsAvailable(32)=.true.
        Met_var_IsAvailable(33)=.true.
        ! Precipitation
        !Met_var_IsAvailable(40)=.true.; Met_var_NC_names(40)="Categorical_Rain"
        !Met_var_IsAvailable(41)=.true.; Met_var_NC_names(41)="Categorical_Snow"
        !Met_var_IsAvailable(42)=.true.; Met_var_NC_names(42)="Categorical_Freezing_Rain"
        !Met_var_IsAvailable(43)=.true.; Met_var_NC_names(43)="Categorical_Ice_Pellets"
        !Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="Precipitation_rate"
        !Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_Precipitation_Rate"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.21)then
        ! GFS 1.0 deg
          ! http://www.nco.ncep.noaa.gov/pmb/products/gfs/

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "GFS 0.5-degree"

        MR_iGridCode = 3
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric3"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.; Met_var_GRIB2_DPcPnSt(10,1:4)=(/0, 3, 18, 1/)
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(32)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.22)then
        ! GFS 0.25 deg
          ! http://www.nco.ncep.noaa.gov/pmb/products/gfs/
          ! http://motherlode.ucar.edu/native/conduit/data/nccf/com/gfs/prod/
          !  

        if(MR_iversion.eq.-1)MR_iversion = 0 ! forecasts are all v.0
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "GFS 0.25"

        MR_iGridCode = 193
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric3"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Surface
        Met_var_IsAvailable(10)=.true.; Met_var_GRIB2_DPcPnSt(10,1:4)=(/0, 3, 18, 1/)
        Met_var_IsAvailable(11)=.true.
        Met_var_IsAvailable(12)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        Met_var_IsAvailable(32)=.true.
        ! Precipitation
        Met_var_IsAvailable(40)=.true.; Met_var_NC_names(40)="Categorical_Rain"
        Met_var_IsAvailable(41)=.true.; Met_var_NC_names(41)="Categorical_Snow"
        Met_var_IsAvailable(42)=.true.; Met_var_NC_names(42)="Categorical_Freezing_Rain"
        Met_var_IsAvailable(43)=.true.; Met_var_NC_names(43)="Categorical_Ice_Pellets"
        Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="Precipitation_rate"
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="Convective_Precipitation_Rate"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.23)then
         ! NCEP / DOE reanalysis 2.5 degree files 
         ! https://rda.ucar.edu/datasets/ds091.0

        if(MR_iversion.eq.-1)MR_iversion = 0 ! latest is 0, but deprecated
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NCEP / DOE reanalysis 2.5 degree files"

        MR_iGridCode = 2
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

        Met_var_GRIB1_Table(1:MR_MAXVARS) = 132

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="Vertical_velocity_isobaric"
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="Vertical_velocity_isobaric"
        ! Moisture
        Met_var_IsAvailable(30)=.true.
        ! Precipitation
        Met_var_IsAvailable(44)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.24)then
         ! NASA-MERRA reanalysis 0.625 x 0.5 degree files 

        if(MR_iversion.eq.-1)MR_iversion = 2 ! v1 finished in 2016
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NASA-MERRA-2 reanalysis 0.625/0.5 degree files"

        MR_iGridCode = 1024
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "lev"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"

        Met_dim_fac(1) = 1.0_sp/60.0_sp

        ! Momentum / State variables
        !   Available in MERRA2_400.inst3_3d_asm_Np.YYYYMMDD.nc4
        !    from https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I3NPASM.5.12.4
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="H"
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="OMEGA"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="OMEGA"
        ! Moisture
        !   Available in MERRA2_400.inst3_3d_asm_Np.YYYYMMDD.nc4
        !    from https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I3NPASM.5.12.4
        Met_var_IsAvailable(30)=.true.; Met_var_NC_names(30)="RH"          ! float percent
        Met_var_IsAvailable(31)=.true.; Met_var_NC_names(31)="QV"          ! float cloud liquid water mixing ratio kg/kg
        Met_var_IsAvailable(32)=.true.; Met_var_NC_names(32)="QL"          ! float cloud liquid water mixing ratio kg/kg
        Met_var_IsAvailable(33)=.true.; Met_var_NC_names(33)="QI"          ! float cloud ice mixing ratio kg/kg

        fill_value_sp = 1.0e15_sp

      elseif (MR_iwindformat.eq.25)then
         ! NCEP/NCAR reanalysis 2.5 degree files 

        if(MR_iversion.eq.-1)MR_iversion = 0 !
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NCEP/NCAR reanalysis 2.5 degree files"

        MR_iGridCode = 2
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

        if(MR_iwind.eq.4)then
           ! https://rda.ucar.edu/datasets/ds090.0

          Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
          Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric"
          Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
          Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"
  
          ! Momentum / State variables
          Met_var_IsAvailable(1)=.true.
          Met_var_IsAvailable(2)=.true.
          Met_var_IsAvailable(3)=.true.
          Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="Pressure_vertical_velocity_isobaric"
          Met_var_IsAvailable(5)=.true.
          Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="Pressure_vertical_velocity_isobaric"
          ! Moisture
          !Met_var_IsAvailable(30)=.true.
  
          fill_value_sp = -9999.0_sp

        elseif(MR_iwind.eq.5)then
          MR_iw5_hours_per_file = 8760.0_dp

          Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
          Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "level"
          Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
          Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"
  
          ! Momentum / State variables
          Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="hgt"        ! short m^2/s^2 (32066.f,1.f)
          Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="uwnd"       ! short m/s (202.66f,0.01f)
          Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="vwnd"       ! short m/s (202.66f,0.01f)
          Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="omega"      ! short Pa/s (29.765f,0.001f)
          Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="air"        ! short K (477.66f,0.01f)
          Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="omega"      ! short Pa/s (29.765f,0.001f)
          ! Atmospheric Structure
          !Met_var_IsAvailable(20)=.true.
          !Met_var_IsAvailable(21)=.true.
          ! Moisture
          !Met_var_IsAvailable(30)=.true.; Met_var_NC_names(30)="rhum"      ! short  (302.66f,0.01f)
          !Met_var_IsAvailable(31)=.true.; Met_var_NC_names(31)="shum"      ! short SpecHum ~ mixing ratio kg/kg(0.032666f,1.e-06f)
          !Met_var_IsAvailable(32)=.true.; Met_var_NC_names(32)="shum"      ! short should really be QL (liquid)
          ! Precipitation
          !Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="prate"     ! short surface precipitation rate (kg/m2/s) (0.0032765f,1.e-07f)
          !Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="cprat"     ! short surface convective precip kg/m2/s (0.0031765f,1.e-07f)
  
          fill_value_sp = -9999.0_sp

        else
          write(MR_global_error,*)"MR ERROR : NCEP Reannalysis provided only as iwind 4 or iwind 5"
          stop 1
        endif

      elseif (MR_iwindformat.eq.26)then
         ! JRA-55 reanalysis 1.25 degree files 
         ! https://rda.ucar.edu/datasets/ds628.0/

        if(MR_iversion.eq.-1)MR_iversion = 0 ! 
        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "JRA-55 reanalysis 1.25 degree files"

        MR_iGridCode = 45
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.
        MR_iw5_hours_per_file = 672.0_dp

        Met_var_GRIB1_Table(1:MR_MAXVARS) = 200

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "initial_time0_hours"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "lv_ISBL1"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "g0_lat_2"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "g0_lon_3"

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="HGT_GDS0_ISBL"
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="UGRD_GDS0_ISBL"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="VGRD_GDS0_ISBL"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="VVEL_GDS0_ISBL"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="TMP_GDS0_ISBL"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="VVEL_GDS0_ISBL"
        ! Moisture
        Met_var_IsAvailable(31)=.true.; Met_var_GRIB_names(31)= "SPFH_GDS0_ISBL"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.27)then
         ! NOAA-CIRES reanalysis 2.0 degree files 
         ! https://rda.ucar.edu/datasets/ds131.2/
         ! https://www.esrl.noaa.gov/psd/data/gridded/data.20thC_ReanV2c.pressure.html

        if(MR_iversion.eq.-1)then
          MR_iversion = 3 ! This is the latest
        endif
        if(MR_iversion.eq.2)then

          write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                    "NOAA-CIRES reanalysis 2.0 degree files"
  
          MR_iGridCode = 1027
          call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
          MR_Reannalysis = .true.
          MR_iw5_hours_per_file = 8760.0_dp
  
          Met_var_GRIB1_Table(1:MR_MAXVARS) = 2
  
          if(MR_Use_RDA)then
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Version 2c (https://rda.ucar.edu/datasets/ds131.2/)
            ! Note: these must be converted from grib using
            !   ncl_convert2nc pgrbanl_mean_1912_VVEL_pres.grib -L
            ! netcdf-java does not seem to work on these
            Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "initial_time0_hours"
            Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "lv_ISBL1"
            Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "g0_lat_2"
            Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "g0_lon_3"
    
            ! Momentum / State variables
            Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="HGT_GDS0_ISBL"
            Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U_GRD_GDS0_ISBL"
            Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V_GRD_GDS0_ISBL"
            Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="V_VEL_GDS0_ISBL"
            Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="TMP_GDS0_ISBL"
            Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="V_VEL_GDS0_ISBL"
          else 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Version 2c (https://www.esrl.noaa.gov/psd/data/gridded/data.20thC_ReanV2c.pressure.html)
            Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
            Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "level"
            Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
            Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"
    
            ! Momentum / State variables
            Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="hgt"
            Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="uwnd"
            Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="vwnd"
            Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="omega"
            Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="air"
            Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="omega"
          endif
          fill_value_sp = 1.0e+20_sp
        elseif(MR_iversion.eq.3)then
          ! NOAA-CIRES-DOE Twentieth Century Reanalysis Version 3 
          ! https://rda.ucar.edu/datasets/ds131.3/
          ! Note: this is also avalable from
          ! https://psl.noaa.gov/data/gridded/data.20thC_ReanV3.html, but might
          ! need format verification
          write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                    "NOAA-CIRES-DOE reanalysis v3"

          MR_iGridCode = 1027
          call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
          MR_Reannalysis = .true.
          MR_iw5_hours_per_file = 8760.0_dp
  
          Met_var_GRIB1_Table(1:MR_MAXVARS) = 2
  
          Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
          Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaricInhPa"
          Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "latitude"
          Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "longitude"
  
          ! Momentum / State variables
          Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="gh"
          Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="u"
          Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="v"
          Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="w"
          Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="t"
          Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="w"
          fill_value_sp = 9999.0_dp
        endif

      elseif (MR_iwindformat.eq.28)then
         ! ECMWF Interim Reanalysis (ERA-Interim)
         ! https://rda.ucar.edu/datasets/ds627.0

        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "ECMWF Interim Reanalysis (ERA-Interim)"

        MR_iGridCode = 170
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

        Met_var_GRIB1_Table(1:MR_MAXVARS) = 128

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="Geopotential_isobaric"
                                       Met_var_GRIB1_Param(1)=129
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U_component_of_wind_isobaric"
                                       Met_var_GRIB1_Param(2)=131
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V_component_of_wind_isobaric"
                                       Met_var_GRIB1_Param(3)=132
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="Vertical_velocity_isobaric"
                                       Met_var_GRIB1_Param(4)=135
        Met_var_IsAvailable(5)=.true.
                                       Met_var_GRIB1_Param(5)=130
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="Vertical_velocity_isobaric"
                                       Met_var_GRIB1_Param(7)=135
        Met_var_IsAvailable(30)=.true.
                                       Met_var_GRIB1_Param(30)=157
        Met_var_IsAvailable(31)=.true.
                                       Met_var_GRIB1_Param(31)=133

        fill_value_sp = -9999.0_sp
        Met_var_conversion_factor(1) = 1.0_sp/9.81_sp

      elseif (MR_iwindformat.eq.29)then
         ! ECMWF ERA5
         ! https://rda.ucar.edu/datasets/ds630.0
         ! Note: files are provided as one variable per file

        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "ECMWF ERA5 reanalysis"

        MR_iGridCode = 1029
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.
        MR_iw5_hours_per_file = 24.0_dp

        Met_var_GRIB1_Table(1:MR_MAXVARS) = 128

        !Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        !Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "level"
        !Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "latitude"
        !Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "longitude"

        !Modified for ERA5 windfiles downloaded from ds633.0 ERA5 files
        !downloaded from https://rda.ucar.edu on May 13, 2020
        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "initial_time0_hours"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "lv_ISBL1"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "g0_lat_2"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "g0_lon_3"
        Met_dim_IsAvailable(5)=.true.; Met_dim_names(5) = "ncl_strlen_0"

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="Z_GDS0_ISBL" ! e5.oper.an.pl.128_129_z.ll025sc.1991061500_1991061523.nc
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U_GDS0_ISBL" ! e5.oper.an.pl.128_129_u.ll025uv.1991061500_1991061523.nc
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V_GDS0_ISBL" ! e5.oper.an.pl.128_129_v.ll025uv.1991061500_1991061523.nc
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="W_GDS0_ISBL" ! e5.oper.an.pl.128_129_w.ll025sc.1991061500_1991061523.nc
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T_GDS0_ISBL" ! e5.oper.an.pl.128_129_t.ll025sc.1991061500_1991061523.nc
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="W_GDS0_ISBL" ! e5.oper.an.pl.128_129_w.ll025sc.19910
        !Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="Z" ! e5.oper.an.pl.128_129_z.regn320sc.2018062000_2018062023.nc
        !Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U" ! e5.oper.an.pl.128_131_u.regn320uv.2018062000_2018062023.nc
        !Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V" ! e5.oper.an.pl.128_132_v.regn320uv.2018062000_2018062023.nc
        !Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="W" ! e5.oper.an.pl.128_135_w.regn320sc.2018062000_2018062023.nc
        !Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T" ! e5.oper.an.pl.128_130_t.regn320sc.2018062000_2018062023.nc
        ! Atmospheric Structure
        !Met_var_IsAvailable(23)=.true.; Met_var_NC_names(23)="CC" ! e5.oper.an.pl.128_248_cc.regn320sc.2018062000_2018062023.nc
        ! Moisture
        !Met_var_IsAvailable(30)=.true.; Met_var_NC_names(30)="R"   ! e5.oper.an.pl.128_157_r.regn320sc.2018062000_2018062023.nc
        !Met_var_IsAvailable(31)=.true.; Met_var_NC_names(31)="Q"   ! e5.oper.an.pl.128_133_q.regn320sc.2018062000_2018062023.nc
        !Met_var_IsAvailable(32)=.true.; Met_var_NC_names(32)="CLWC"! e5.oper.an.pl.128_246_clwc.regn320sc.2018062000_2018062023.nc
        !Met_var_IsAvailable(32)=.true.; Met_var_NC_names(32)="CIWC"! e5.oper.an.pl.128_247_ciwc.regn320sc.2018062000_2018062023.nc
        ! Precipitation
        !Met_var_IsAvailable(32)=.true.; Met_var_NC_names(32)="CRWC"! e5.oper.an.pl.128_075_crwc.regn320sc.2018062000_2018062023.nc
        !Met_var_IsAvailable(32)=.true.; Met_var_NC_names(32)="CSWC"! e5.oper.an.pl.128_076_cswc.regn320sc.2018062000_2018062023.nc

        fill_value_sp = -9999.0_sp
        Met_var_conversion_factor(1) = 1.0_sp/9.81_sp

      elseif (MR_iwindformat.eq.30)then
         ! ECMWF ERA 20C
         ! https://rda.ucar.edu/datasets/ds626.0
         ! Note: files are provided as one variable per file

        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "ECMWF ERA20C reanalysis"

        MR_iGridCode = 1030
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.
        MR_iw5_hours_per_file = 672.0_dp

        Met_var_GRIB1_Table(1:MR_MAXVARS) = 128


        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "initial_time0_hours"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "lv_ISBL1"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "g4_lat_2"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "g4_lon_3"

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="Z_GDS4_ISBL" ! e5.oper.an.pl.128_129_z.regn320sc.2018062000_2018062023.nc
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U_GDS4_ISBL" ! e5.oper.an.pl.128_131_u.regn320uv.2018062000_2018062023.nc
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V_GDS4_ISBL" ! e5.oper.an.pl.128_132_v.regn320uv.2018062000_2018062023.nc
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="W_GDS4_ISBL" ! e5.oper.an.pl.128_135_w.regn320sc.2018062000_2018062023.nc
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T_GDS4_ISBL" ! e5.oper.an.pl.128_130_t.regn320sc.2018062000_2018062023.nc
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="W_GDS4_ISBL" ! e5.oper.an.pl.128_135_w.regn320sc.201

        fill_value_sp = -9999.0_sp

        Met_var_conversion_factor(1) = 1.0_sp/9.81_sp

      elseif (MR_iwindformat.eq.32)then
         ! Air Force Weather Agency subcenter = 0
         ! GALWEM

        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "Air Force Weather Agency subcenter = 0"

        MR_iGridCode = 1032
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "isobaric"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.
        Met_var_IsAvailable(2)=.true.
        Met_var_IsAvailable(3)=.true.
        Met_var_IsAvailable(4)=.true.
        Met_var_IsAvailable(5)=.true.
        Met_var_IsAvailable(7)=.true.
        ! Moisture
        Met_var_IsAvailable(30)=.true.

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.33)then
         ! CCSM3.0 Community Atmosphere Model (CAM)
         ! http://www.cesm.ucar.edu/models/atm-cam/
         ! peleoclimate monthly averages

        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "CCSM3.0 Community Atmosphere Model (CAM)"

        MR_iGridCode = 1033
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "lev"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"

        Met_dim_fac(1) = 24.0

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="Z3"      ! float m
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="OMEGA"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="OMEGA"
        ! Moisture
        Met_var_IsAvailable(30)=.true.; Met_var_NC_names(30)="RELHUM"
        Met_var_IsAvailable(31)=.true.; Met_var_NC_names(31)="Q"

        fill_value_sp = -9999.0_sp

      elseif (MR_iwindformat.eq.40)then
         ! NASA-GEOS Cp

        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NASA-GEOS Cp"

        MR_iGridCode = 1040
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "lev"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"

        Met_dim_fac(1) = 1.0_sp/60.0_sp

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="H"
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="OMEGA"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="OMEGA"

        fill_value_sp = 1.0e15_sp

      elseif (MR_iwindformat.eq.41)then
         ! NASA-GEOS Np

        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NASA-GEOS Np"

        MR_iGridCode = 1041
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .false.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "lev"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "lat"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "lon"

        Met_dim_fac(1) = 1.0_sp/60.0_sp

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="H"
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="OMEGA"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="OMEGA"

        fill_value_sp = 1.0e15_sp

      elseif (MR_iwindformat.eq.50)then
         ! WRF - output

        write(MR_global_info,*)"  NWP format to be used = ",MR_iwindformat,&
                  "WRF"

        MR_iGridCode = 1050
        call MR_Set_Met_NCEPGeoGrid(MR_iGridCode)
        MR_Reannalysis = .true.

        Met_dim_IsAvailable(1)=.true.; Met_dim_names(1) = "Time"
        Met_dim_IsAvailable(2)=.true.; Met_dim_names(2) = "bottom_top"
        Met_dim_IsAvailable(3)=.true.; Met_dim_names(3) = "south_north"
        Met_dim_IsAvailable(4)=.true.; Met_dim_names(4) = "west_east"
        ! for pressure, read "P"  :: perturbation pressure
        !               and  "PB" :: base pressure

        ! for geopotential, read "PH"  :: perturbation geopotential
        !                   and  "PHB" :: base-state geopotential

        ! Momentum / State variables
        Met_var_IsAvailable(1)=.true.; Met_var_NC_names(1)="PHB"
        Met_var_IsAvailable(2)=.true.; Met_var_NC_names(2)="U"
        Met_var_IsAvailable(3)=.true.; Met_var_NC_names(3)="V"
        Met_var_IsAvailable(4)=.true.; Met_var_NC_names(4)="W"
        Met_var_IsAvailable(5)=.true.; Met_var_NC_names(5)="T"      ! float K perturbation potential temperature (theta-t0)
        Met_var_IsAvailable(6)=.true.; Met_var_NC_names(6)="PB"
        Met_var_IsAvailable(7)=.true.; Met_var_NC_names(7)="W"

        ! Surface
        Met_var_IsAvailable(10)=.true.; Met_var_NC_names(10)="PBLH"
        Met_var_IsAvailable(11)=.true.; Met_var_NC_names(11)="U10"
                                        Met_var_ndim(11)=3
        Met_var_IsAvailable(12)=.true.; Met_var_NC_names(12)="V10"
                                        Met_var_ndim(12)=3
        Met_var_IsAvailable(13)=.true.; Met_var_NC_names(13)="UST"
        Met_var_IsAvailable(15)=.true.; Met_var_NC_names(15)="SNOWH"
        Met_var_IsAvailable(16)=.true.; Met_var_NC_names(16)="SMOIS" !Soil moisture m3 m-3
        ! Moisture
        Met_var_IsAvailable(31)=.true.; Met_var_NC_names(31)="QVAPOR" !QV (specific humidity)
        ! Precipitation
        Met_var_IsAvailable(44)=.true.; Met_var_NC_names(44)="RAINC"   ! ACCUMULATED TOTAL CUMULUS PRECIPITATION in mm
        Met_var_IsAvailable(45)=.true.; Met_var_NC_names(45)="RAINNC"  ! ACCUMULATED TOTAL GRID SCALE PRECIPITATION in mm

        fill_value_sp = 1.0e20_sp

        Met_var_conversion_factor(1) = 1.0_sp/9.81_sp

      else
        write(MR_global_error,*)'MR ERROR : MR_iwindformat not supported'
        write(MR_global_error,*)'           MR_iwindformat=',MR_iwindformat,&
                                     '. Program stopped in MetReader.f90'
        stop 1
      endif

      write(MR_global_info,*)"                grid ID = ",MR_iGridCode
      select case (MR_idataFormat)
      case(1)
        write(MR_global_info,*)"            data format = ",MR_idataFormat,"ASCII"
      case(2)
        write(MR_global_info,*)"            data format = ",MR_idataFormat,"NETCDF"
#ifndef USENETCDF
        write(MR_global_error,*)"MR ERROR: Met files are netcdf, but MetReader was not"
        write(MR_global_error,*)"          compiled with netcdf support."
        write(MR_global_error,*)"          Please recompile MetReader with netcdf support"
        stop 1
#endif
      case(3)
        write(MR_global_info,*)"            data format = ",MR_idataFormat,"GRIB"
#ifndef USEGRIB
        write(MR_global_error,*)"MR ERROR: Met files are grib, but MetReader was not"
        write(MR_global_error,*)"          compiled with grib support."
        write(MR_global_error,*)"          Please recompile MetReader with grib support"
        stop 1
#endif
      case default
        write(MR_global_error,*)"MR ERROR:  MR_idataFormat not 1:4."
        write(MR_global_error,*)"           Exiting in MR_Allocate_FullMetFileList"
        stop 1
      end select

      write(MR_global_info,*)"     Allocating space for ",MR_iwindfiles,"file(s)."

      if (MR_iwind.eq.5)then
        write(MR_global_info,*)"For iwf=5:: one variable per file,"
        write(MR_global_info,*)"Only the directory should be listed. The remaining"
        write(MR_global_info,*)"path is hard-coded."
        ! Reset MR_iwindfiles to 2: only one "file" will be read,
        ! but we need this to be 2 to accommodate runs that might span two years/months
        ! For ERA5 files, data is provided in daily files so MR_iwindfiles might be > 2
        if(MR_Comp_Time_in_hours.gt.0.0)then
          if(MR_iwindformat.eq.25)then  ! NCEP 2.5 degree
            MR_iwindfiles = ceiling(MR_Comp_Time_in_hours/MR_iw5_hours_per_file) +1  ! These are yearly files
          elseif(MR_iwindformat.eq.26)then ! 
            MR_iwindfiles = ceiling(MR_Comp_Time_in_hours/MR_iw5_hours_per_file) +1 ! Monthly files
          elseif(MR_iwindformat.eq.27)then
            MR_iwindfiles = ceiling(MR_Comp_Time_in_hours/MR_iw5_hours_per_file) +1 ! yearly files
          elseif(MR_iwindformat.eq.29)then
            MR_iwindfiles = ceiling(MR_Comp_Time_in_hours/MR_iw5_hours_per_file) +1 ! daily files
          elseif(MR_iwindformat.eq.30)then
            MR_iwindfiles = ceiling(MR_Comp_Time_in_hours/MR_iw5_hours_per_file) +1 ! monthly
          endif
        else
          MR_iwindfiles = 2
        endif
      endif

      allocate (MR_windfiles(MR_iwindfiles))
      do i=1,MR_iwindfiles
        write(MR_windfiles(i),'(130x)')
      enddo
      allocate (MR_windfiles_nt_fullmet(MR_iwindfiles))
      MR_windfiles_nt_fullmet(:)=0
      if(MR_idataFormat.eq.3)then
        allocate (MR_windfiles_Have_GRIB_index(MR_iwindfiles))
          ! This will be reset to true if the index files are found
        MR_windfiles_Have_GRIB_index = .false.
        allocate (MR_windfiles_GRIB_index(MR_iwindfiles))
        do i=1,MR_iwindfiles
          write(MR_windfiles_GRIB_index(i),'(130x)')
        enddo
      endif

      if(MR_iwind.eq.1)then
        ! For the 1d profile or radiosonde case, igrid is used for the number of sonde
        ! locations.  If it is not provided, set it to one
        if(igrid.eq.0)then
          MR_iGridCode = 1
          MR_nSnd_Locs = MR_iGridCode
        else
          MR_iGridCode = igrid
          MR_nSnd_Locs = MR_iGridCode
        endif
          ! Now make sure that the number of windfiles is a multiple of the number of locations
          ! since this will be the number of timesteps
        if(mod(MR_iwindfiles,MR_nSnd_Locs).eq.0)then
          MR_Snd_nt_fullmet = MR_iwindfiles / MR_nSnd_Locs
        else
          write(MR_global_error,*)"MR ERROR:  The grid code for 1d ascii input is interpreted to be"
          write(MR_global_error,*)"           the number of sonde locations.  Each group of sondes can"
          write(MR_global_error,*)"           be repeated, correspoding to multiple timesteps.  The"
          write(MR_global_error,*)"           number of windfiles must be a multiple of the number of"
          write(MR_global_error,*)"           locations"
          write(MR_global_error,*)"                   MR_iwind = ",MR_iwind
          write(MR_global_error,*)"             MR_iwindformat = ",MR_iwindformat
          write(MR_global_error,*)"               MR_iGridCode = ",MR_iGridCode
          write(MR_global_error,*)"             MR_idataFormat = ",MR_idataFormat
          write(MR_global_error,*)"               MR_nSnd_Locs = ",MR_iGridCode
          write(MR_global_error,*)"          MR_Snd_nt_fullmet = ",real(MR_iwindfiles)/real(MR_nSnd_Locs)
          stop 1
        endif
      endif

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      CALLED_MR_Allocate_FullMetFileList = .true.

      return

      end subroutine MR_Allocate_FullMetFileList


!##############################################################################
!
!     MR_Read_Met_DimVars
!
!     This subroutine expects that MR_windfiles has been filled by the calling
!     program, and checks for their existence.
!
!     Depending on the NWP data format, a subroutine is called to read the spatial
!     structure of the NWP files, followed by one for the temporal span
!      e.g. MR_Read_Met_DimVars_netcdf
!           MR_Read_Met_Times_netcdf
!
!     From the calling program, this is called once the names of the NWP files
!     is specified.  If a custom netcdf template file is to be used (iwf=0), then
!     the variable MR_iwf_template must also be filled so that it can be read
!     from MR_Read_Met_DimVars_[].  
!
!     After this subroutine completes, the following variables will be set:
!       All the projection parameters of NWP grid
!       The lengths of all the dimensions of the file
!       p_fullmet_sp (converted to Pa)
!       x_fullmet_sp, y_fullmet_sp
!       IsLatLon_MetGrid, IsGlobal_MetGrid, IsRegular_MetGrid 
!
!     The next step is for the calling program to specify the projection parameters
!     of the computational grid (i.e. the grid MetReader should be returing values to)
!
!##############################################################################

      subroutine MR_Read_Met_DimVars(iy)

      implicit none

      integer, optional,intent(in) :: iy  ! Note: this is only needed for iw=5
                                          !       since we need to know how many
                                          !       metsteps to allocate

      integer      :: i
      logical      :: IsThere
      character(len=130) :: tmp_str
      character(len=130) :: infile
      integer            :: ivar
      real(kind=8)       :: inhour

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Read_Met_DimVars                                   ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .false., &  ! CALLED_MR_Read_Met_DimVars
                                 .false., &  ! CALLED_MR_Set_CompProjection
                                 .false., &  ! CALLED_MR_Initialize_Met_Grids
                                 .false.)    ! CALLED_MR_Set_Met_Times

#ifdef USEPOINTERS
      ! Need to add a check here for the case where MR_windfiles is a pointer, but has not been filled
#else
      if(.not.allocated(MR_windfiles))then
        write(MR_global_error,*)"MR ERROR:  The list of windfile names, MR_windfiles, has not been"
        write(MR_global_error,*)"           allocated.  The calling program must allocate this"
        write(MR_global_error,*)"           array via the subroutine MR_Allocate_FullMetFileList."
        stop 1
      endif
#endif

      if(MR_iwind.eq.5)then
        ! For iwind=5 files (NCEP 2.5 degree reanalysis, NOAA, etc. ), only the directory
        ! was read into slot 1 of MR_windfiles(:).  We need to copy to all other slots
        ! to make sure we don't throw an error
        do i=1,MR_iwindfiles
          MR_windfiles(i)   = MR_windfiles(1)
        enddo
        if(present(iy)) then
          ! This is needed at this point for allocating the number
          ! of steps per file (this depends on the year), but this
          ! is reset in MR_Set_Met_Times when the actual start time
          ! is given
          MR_Comp_StartYear = iy
        else
          write(MR_global_info,*)"MR WARNING: If the iw=5 (NCEP Reannalysis, NOAA Reannalysis, "
          write(MR_global_info,*)"            etc.) are used, then MR_Read_Met_DimVars"
          write(MR_global_info,*)"            should be called with a start year.  This is needed"
          write(MR_global_info,*)"            to allocate the correct number of time steps per file."
          write(MR_global_info,*)"            Setting MR_Comp_StartYear to 2018 for a non-leap year."
          write(MR_global_info,*)"            However, MR_Comp_StartYear will be checked later to"
          write(MR_global_info,*)"            verify that the start year is not a leap year."
          write(MR_global_info,*)"            If there is an inconsistancy, the program will stop."
          write(MR_global_info,*)"            If MR_Comp_StartYear is changed to a leap year outside"
          write(MR_global_info,*)"            of MetReader, then the results will be incorrect."
          MR_Comp_StartYear = 2018
        endif

      endif

      ! Verify that the first windfile has been changed from its initialized value
      tmp_str = MR_windfiles(1)
      if(tmp_str(1:15).eq.'              ')then
        write(MR_global_error,*)"MR ERROR: The array MR_windfiles appears to not be set."
        write(MR_global_error,*)"          Please have the calling program write the names of"
        write(MR_global_error,*)"          the windfiles to MR_windfiles(1:MR_iwindfiles)"
        write(MR_global_error,*)" "
        write(MR_global_error,*)"          Contents of the first slot is -",MR_windfiles(1),"-"
        stop 1
      endif

      ! Check the existence of the wind files
      write(MR_global_info,*)"  Verifying existance of windfiles:"
      ! Note iwf=5 cases have the number of windfiles (MR_iwindfiles)
      ! modified in MR_Allocate_FullMetFileList to be the number of
      ! anticipated files based on the length of the simulation and the number
      ! of steps per file.
      if(MR_iwind.ne.5)then
        do i=1,MR_iwindfiles
          inquire( file=adjustl(trim(MR_windfiles(i))), exist=IsThere )
          write(MR_global_info,*)"     ",i,adjustl(trim(MR_windfiles(i))),IsThere
          if(.not.IsThere)then
            write(MR_global_error,*)"MR ERROR: Could not find windfile ",i
            write(MR_global_error,*)"          ",adjustl(trim(MR_windfiles(i)))
            stop 1
          endif
        enddo
      else
        MR_iw5_root = MR_windfiles(1)
        write(*,*)MR_iw5_root,MR_iwindfiles
        do i = 1,MR_iwindfiles-1
          inhour = MR_Comp_StartHour + (i-1)*MR_iw5_hours_per_file
          do ivar = 1,5
            call MR_Set_iwind5_filenames(inhour,ivar,infile)
            inquire( file=adjustl(trim(infile)), exist=IsThere )
            write(MR_global_info,*)" ",i,trim(adjustl(trim(infile))),IsThere
            if(.not.IsThere)then
              write(MR_global_error,*)"MR ERROR: Could not find windfile ",i
              write(MR_global_error,*)"          ",adjustl(trim(infile))
              stop 1
            endif
          enddo
        enddo
      endif

      ! Now set up the full spatial and temporal grids
      select case (MR_iwind)
      case(1)   ! if we're using a 1-D wind sounding
        call MR_Read_Met_DimVars_ASCII_1d
      case(2)
        !call Set_Read_Met_DimVars_ASCII_3d
      case (3:5)
         ! call routine to populate variable lists and dimension axies based on
         ! MR_iwindformat; also do a bit of error checking by reading netcdf files
        if(MR_idataFormat.eq.2)then
#ifdef USENETCDF
          call MR_Read_Met_Times_netcdf
          call MR_Read_Met_DimVars_netcdf
#endif
        elseif(MR_idataFormat.eq.3)then
#ifdef USEGRIB
          call MR_Read_Met_Times_GRIB
          call MR_Read_Met_DimVars_GRIB
#endif
        else
          write(MR_global_error,*)"MR ERROR: Unknown MR_idataFormat:",MR_idataFormat
          stop 1
        endif
      case default
        write(MR_global_error,*)"MR ERROR:  MR_iwind not 1-4."
        write(MR_global_error,*)"           Exiting in MR_Read_Met_DimVars"
        stop 1
      end select

      CALLED_MR_Read_Met_DimVars = .true.

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      end subroutine MR_Read_Met_DimVars

!##############################################################################


!##############################################################################
!
!     MR_Set_CompProjection
!
!     This subroutine only sets the projection parameters for the computational
!     grid for MetReader.  The library needs to know these parameters in order
!     to know how to map the met grid onto the comp grid and whether or not to
!     rotate wind vectors.
!
!     This subroutine can be called from the calling program at any time, but
!     must be called before MR_Initialize_Met_Grids
!
!##############################################################################

      subroutine MR_Set_CompProjection(LL_flag,ipf,lam0,phi0,phi1,phi2,ko,Re)

      implicit none

      logical     ,intent(in) :: LL_flag
      integer     ,intent(in) :: ipf
      real(kind=8),intent(in) :: lam0,phi0
      real(kind=8),intent(in) :: phi1
      real(kind=8),intent(in) :: phi2
      real(kind=8),intent(in) :: ko
      real(kind=8),intent(in) :: Re

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Set_CompProjection                                 ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      isLatLon_CompGrid = LL_flag

      Comp_lam0 = 0.0_8
      Comp_phi0 = 0.0_8
      Comp_phi1 = 0.0_8
      Comp_phi2 = 0.0_8
      Comp_k0   = 0.0_8
      Comp_Re   = 0.0_8

      if(.not.isLatLon_CompGrid)then
        Comp_iprojflag = ipf 
        if(Comp_iprojflag.eq.1)then
          ! Polar stereographic
          Comp_lam0    = lam0
          Comp_phi0    = phi0
          Comp_phi1    = phi1
          Comp_k0      = ko
          Comp_Re      = Re
        elseif(Comp_iprojflag.eq.2)then
          ! Albers Equal Area
          Comp_lam0    = lam0
          Comp_phi0    = phi0
          Comp_phi1    = phi1
          Comp_phi2    = phi2
        elseif(Comp_iprojflag.eq.3)then
          ! UTM
          write(MR_global_error,*)"MR ERROR : WARNING: UTM not yet verified"
          stop 1
        elseif(Comp_iprojflag.eq.4)then
          ! Lambert conformal conic (NARR, NAM218, NAM221)
          Comp_lam0    = lam0
          Comp_phi0    = phi0
          Comp_phi1    = phi1
          Comp_phi2    = phi2
          Comp_Re      = Re
        elseif(Comp_iprojflag.eq.5)then
          ! Mercator (NAM196)
          Comp_lam0    = lam0
          Comp_phi0    = phi0
          Comp_Re      = Re
        endif
      endif

      !to do: perform a sanity check on these projection parameters
      ! take a lon/lat, project, then inverse project and check

      CALLED_MR_Set_CompProjection = .true.

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      return

      end subroutine MR_Set_CompProjection


!##############################################################################
!
!     MR_Initialize_Met_Grids
!
!     This subroutine sets up local copies of the computational grid and
!     evaluates the full NWP grid for the subgrid needed, through MR_Set_MetComp_Grids
!     for NWP cases and MR_Set_MetComp_Grids_1dascii for 1-d or sonde cases.
!
!     From the calling program, this must be called after all of:
!       MR_Allocate_FullMetFileList
!       MR_Read_Met_DimVars
!       MR_Set_CompProjection
!
!     In some cases, it is useful to have velocity values saved within MetReader.
!     For example, velocities might be read to calculate diffusivities, stored
!     locally, then reused later for advection.  If the calling program sets the
!     variable MR_Save_Velocities=.true. , then space is allocated for the local
!     copies of velocity.
!
!     Takes as input :: specs of computational grid defined in the calling program
!     Sets: n[x,y,z]_comp and [x,y,z]_comp_sp :: local variables holding computational grid info
!           n[t,x,y,p]_met     :: sets the size of the dimensions of the sub-met grid
!           [x,y,p]_met_sp     :: arrays holding dimension values of the sub-met grid
!           rdphi_MetP_sp      :: length scale along y (in meters)
!           rdlambda_MetP_sp   :: length scale along x (in meters)
!
!##############################################################################

      subroutine MR_Initialize_Met_Grids(nx,ny,nz, &
                                      dumx_sp,dumy_sp,dumz_sp,periodic)

      implicit none

      integer      ,intent(in) :: nx,ny,nz
      real(kind=sp),intent(in) :: dumx_sp(nx)
      real(kind=sp),intent(in) :: dumy_sp(ny)
      real(kind=sp),intent(in) :: dumz_sp(nz)
      logical      ,intent(in) :: periodic

      integer       :: i,j,k

      !                allocates *_Met_P for subset of Met grid on pressure levels
      !                          *_Met_H for subset of Met grid on height levels
      !                          *_comp_H for data regridded to computational grid

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Initialize_Met_Grids                               ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .false., &  ! CALLED_MR_Initialize_Met_Grids
                                 .false.)    ! CALLED_MR_Set_Met_Times
      nx_comp = nx
      ny_comp = ny
      nz_comp = nz
      allocate(x_comp_sp(nx_comp))
      allocate(y_comp_sp(ny_comp))
      allocate(z_comp_sp(nz_comp))
      if(MR_useCompGrid.eqv..false.)then
        ! This is the case where we will not be interpolating to a computational grid, but want access
        ! to the full met grid.  All parameters to this subroutine should have been dummy values
        x_comp_sp(1)  = x_fullmet_sp(1)
        x_comp_sp(nx) = x_fullmet_sp(nx_fullmet)
        y_comp_sp(1)  = min(y_fullmet_sp(1),y_fullmet_sp(ny_fullmet))
        y_comp_sp(ny) = max(y_fullmet_sp(1),y_fullmet_sp(ny_fullmet))
        z_comp_sp(1)  = 0.0_sp
        z_comp_sp(nz) = 1.1_sp*MR_Max_geoH_metP_predicted
      else
        x_comp_sp = dumx_sp
        y_comp_sp = dumy_sp
        z_comp_sp = dumz_sp
      endif
      ! Do some error-checking on these grids to make sure they are
      ! strictly increasing
      do i=1,nx-1
        if(x_comp_sp(i).gt.x_comp_sp(i+1))then
          write(MR_global_error,*)"MR ERROR:  x_comp not strictly increasing"
          write(MR_global_error,*)s_comp_sp
          stop 1
        endif
      enddo

      do j=1,ny-1
        if(y_comp_sp(j).gt.y_comp_sp(j+1))then
          write(MR_global_error,*)"MR ERROR:  y_comp not strictly increasing"
          write(MR_global_error,*)y_comp_sp
          stop 1
        endif
      enddo

      dx_comp = x_comp_sp(2) - x_comp_sp(1)
      dy_comp = abs(y_comp_sp(2) - y_comp_sp(1))
      MaxZ_comp_sp = maxval(z_comp_sp)

      IsPeriodic_CompGrid = periodic

      select case (MR_iwind)
      case(1)   ! if we're using a 1-D wind sounding
        call MR_Set_MetComp_Grids_ASCII_1d
        call MR_Set_Comp2Met_Map
      case(2)
        !call MR_Set_MetComp_Grids_ASCII_3d
      case (3:5)
        ! Now that we have the full grids defined in MR_Read_Met_DimVars_netcdf,
        ! calculate the subgrid needed for the simulation
        call MR_Set_MetComp_Grids

        if(MR_Save_Velocities)then
          write(MR_global_info,*)"Velocities will be saved on the metP grid"
          allocate(MR_vx_metP_last(nx_submet,ny_submet,np_fullmet));MR_vx_metP_last(:,:,:)=0.0_sp
          allocate(MR_vx_metP_next(nx_submet,ny_submet,np_fullmet));MR_vx_metP_next(:,:,:)=0.0_sp
          allocate(MR_vy_metP_last(nx_submet,ny_submet,np_fullmet));MR_vy_metP_last(:,:,:)=0.0_sp
          allocate(MR_vy_metP_next(nx_submet,ny_submet,np_fullmet));MR_vy_metP_next(:,:,:)=0.0_sp
        else
          write(MR_global_info,*)"Velocities not saved"
        endif

      case default
        write(MR_global_error,*)"MR ERROR:  MR_iwind not 1:5."
        write(MR_global_error,*)"           Exiting in MR_Initialize_Met_Grids"
        stop 1
      end select

      ! For LatLon Met grids, calculate some additional geometry terms
      !  These are currently only needed for calculating DelMetP_Dx and
      !  DelMetP_Dy
      if(IsLatLon_MetGrid)then
        allocate(rdphi_MetP_sp(ny_submet,np_fullmet))             ;   rdphi_MetP_sp(:,:)=0.0_sp
        allocate(rdlambda_MetP_sp(nx_submet,ny_submet,np_fullmet));rdlambda_MetP_sp(:,:,:)=0.0_sp
        do k=1,np_fullmet
          if(IsRegular_MetGrid)then
            ! length scale along y (in meters)
            rdphi_MetP_sp(:,k) = dy_met_const*DEG2RAD_MET * (RAD_EARTH_MET+z_approx(k))*1000.0_sp
            do j=1,ny_submet
              ! length scale along x (in meters)
              rdlambda_MetP_sp(:,j,k) =(RAD_EARTH_MET+z_approx(k))*1000.0_sp * &
                                      cos(DEG2RAD_MET*(y_submet_sp(j)-0.5_sp*dy_met_const)) * &
                                      dx_met_const*DEG2RAD_MET
            enddo
          else
            do i=1,nx_submet
              do j=1,ny_submet
                ! length scale along y (in meters)
                rdphi_MetP_sp(:,k) = MR_dy_submet(j)*DEG2RAD_MET * (RAD_EARTH_MET+z_approx(k))*1000.0_sp
                ! length scale along x (in meters)
                rdlambda_MetP_sp(i,j,k) =(RAD_EARTH_MET+z_approx(k))*1000.0_sp * &
                                        cos(DEG2RAD_MET*(y_submet_sp(j)-0.5_sp*MR_dy_submet(j))) * &
                                        MR_dx_submet(i)*DEG2RAD_MET
              enddo
            enddo
          endif
        enddo
      endif

      CALLED_MR_Initialize_Met_Grids = .true.

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      end subroutine MR_Initialize_Met_Grids

!##############################################################################
!
!     MR_Set_SigmaAlt_Scaling

!##############################################################################

      subroutine MR_Set_SigmaAlt_Scaling(nx,ny,nz, &
                                         dum_sp, dumz_sp, dumxy1_sp, &
                                         dum_int)

      implicit none

      integer      ,intent(in) :: nx,ny,nz
      real(kind=sp),intent(in) :: dum_sp
      real(kind=sp),intent(in) :: dumz_sp(nz)
      real(kind=sp),intent(in) :: dumxy1_sp(nx,ny)
      integer      ,intent(in) :: dum_int

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Set_SigmaAlt_Scaling                               ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      MR_use_SigmaAlt   = .true.
      MR_ztop           = dum_sp
      MR_ZScaling_ID    = dum_int
      allocate(s_comp_sp(nz));     s_comp_sp(1:nz) = dumz_sp(1:nz)
      allocate(MR_zsurf(nx,ny));   MR_zsurf(1:nx,1:ny) = dumxy1_sp(1:nx,1:ny)
      allocate(MR_jacob(nx,ny))
      if (MR_ZScaling_ID.eq.2)then
        MR_jacob(1:nx,1:ny) = MR_ztop - MR_zsurf(1:nx,1:ny)
      else
        MR_jacob(1:nx,1:ny) = 1.0_sp
      endif

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      end subroutine MR_Set_SigmaAlt_Scaling

!##############################################################################
!
!     MR_Set_Met_Times
!
!     This subroutine opens each file of the MR_windfile list, determines the hour
!     of each time step and logs the files/time-steps needed to bracket the
!     simulation duration.
!
!     Sets: MetStep_File(i=1:MR_MetSteps_Total)   :: contains file name for met step i
!           MetStep_findex(i=1:MR_MetSteps_Total) :: contains the index of the file in
!                                            the MR_windfiles(:) list
!           MetStep_tindex(i=1:MR_MetSteps_Total) :: contains the index within the file
!                                            for met step i
!
!##############################################################################

      subroutine MR_Set_Met_Times(eStartHour,Duration)

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer, parameter :: NT_MAXOUT = 40

      real(kind=8),intent(in) :: eStartHour
      real(kind=8),intent(in) :: Duration

      integer :: i
      integer :: iw
      !real(kind=8)       :: HS_HourOfDay
      !integer            :: HS_YearOfEvent
      !integer            :: HS_MonthOfEvent
      !integer            :: HS_DayOfEvent
      !integer            :: HS_DayOfYear
      integer            :: iwstep
      integer            :: istep
      real(kind=8)       :: stephour
      logical :: Found_First_Step = .false.
      logical :: Found_Last_Step  = .false.
      integer :: nMetSteps_Comp   = 0
      real(kind=8) :: StepInterval
      real(kind=8) :: met_t1,met_t2,met_dt1
      logical      :: prestep, poststep

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
        character (len=13) function HS_yyyymmddhhmm_since(HoursSince,byear,useLeaps)
          real(kind=8)          :: HoursSince
          integer               :: byear
          logical               :: useLeaps
        end function HS_yyyymmddhhmm_since

      END INTERFACE

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Set_Met_Times                                      ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .false., &  ! CALLED_MR_Set_CompProjection
                                 .false., &  ! CALLED_MR_Initialize_Met_Grids
                                 .false.)    ! CALLED_MR_Set_Met_Times

      if(MR_useCompTime.eqv..false.)then
        ! This is for the case where only one file is provided and steps will be accessed
        ! through the file step index.  e.g. used for evaluating windfiles
        ! Here, we assign the start hour to the first step and the duration with the last step
        MR_Comp_StartHour     = MR_windfile_starthour(1)+MR_windfile_stephour(1,1)
        MR_Comp_Time_in_hours = MR_windfile_stephour(1,nt_fullmet)
      else
        MR_Comp_StartHour        = eStartHour
        MR_Comp_Time_in_hours    = Duration
      endif
      MR_Comp_StartYear        = HS_YearOfEvent(MR_Comp_StartHour,MR_BaseYear,MR_useLeap)
      MR_Comp_StartMonth       = HS_MonthOfEvent(MR_Comp_StartHour,MR_BaseYear,MR_useLeap)
      MR_Comp_StartDay         = HS_DayOfEvent(MR_Comp_StartHour,MR_BaseYear,MR_useLeap)

      ! Now we want to loop through all the steps of each windfile and find the
      ! file/step that immediately preceeds the time needed (MR_Comp_StartHour)
      !  First, define MR_Comp_StartHour if we are running a forecast run

      ! Now is a good time to check to make sure the MR_Comp_StartHour is within
      !  the range of data
      ! First we need the step interval for the first step
      ! This is for checking pre- and post-steps
      if(MR_iwind.eq.1)then
        if(MR_iwindformat.eq.1)then
          ! 1-D ASCII
          if(MR_Snd_nt_fullmet.eq.1)then
            StepInterval = 1000.0_8 ! some large number
          else
            met_t1 = MR_windfile_starthour(1)
            met_t2 = MR_windfile_starthour(MR_nSnd_Locs+1)
            StepInterval = met_t2 - met_t1
          endif
        elseif(MR_iwindformat.eq.2)then
          ! 1-D Radiosonde
          StepInterval = 12.0_8
        endif
      else
        ! For all other cases, just check the arrays read in MR_Read_Met_DimVars
        met_t1  = MR_windfile_starthour(1)+MR_windfile_stephour(1,1)
        if(nt_fullmet.gt.1)then
            ! If there are multiple steps per file, use the second step
          met_t2 = MR_windfile_starthour(1)+MR_windfile_stephour(1,2)
        elseif(MR_iwindfiles.gt.1)then
            ! If only 1 step/file, use the next file (if present)
          met_t2 = MR_windfile_starthour(1)+MR_windfile_stephour(2,1)
        else
          !
          write(MR_global_info,*)"WARNING: Only one time step available."
          write(MR_global_info,*)"         Setting interval step to 1000.0"
          met_t2 = met_t1 + 1000.0
        endif
        StepInterval = met_t2 - met_t1
      endif

      !   Checking if a prestep is needed
      if(MR_iwind.eq.1.and.MR_iwindformat.eq.1)then
        ! For the ASCII profile cases (not the radiosonde), hours are given as offset
        ! from the start time.  So reset all hours to relative to eStartHour
        write(MR_global_info,*)"Note:  The hours value in 1d ascii profiles are interpreted as offset"
        write(MR_global_info,*)"       hours.  Shifting the file time to reflect the requested start"
        write(MR_global_info,*)"       time plus offset."
        do i=1,MR_Snd_nt_fullmet
          MR_windfile_starthour(i) = MR_windfile_starthour(i) + MR_Comp_StartHour
        enddo
      endif
      met_t1  = MR_windfile_starthour(1)+MR_windfile_stephour(1,1)
      met_dt1 = StepInterval
      if(MR_Comp_StartHour.lt.met_t1)then
        ! Start time requested is before that available, check if we can extrapolate
        if(MR_Comp_StartHour.ge.met_t1-met_dt1)then
          write(MR_global_info,*)"WARNING: Start time is before the first time step"
          write(MR_global_info,*)"         However, it is within one interval so we will"
          write(MR_global_info,*)"         apply the value at MR_windfile_starthour(1)."
          write(MR_global_info,*)"         Using a time step interval of ",StepInterval
          if(StepInterval.gt.720.0_dp)then
            if(MR_iwind.eq.1.and.MR_iwindformat.eq.1.and.nt_fullmet.eq.1)then
              write(MR_global_info,*)"         Note: This interval is set high since only one"
              write(MR_global_info,*)"               sonde file was provided."
            endif
          endif
          prestep = .true.
        else
          write(MR_global_error,*)"MR ERROR: Start time is before the first available data and"
          write(MR_global_error,*)"       cannot be extrapolated."
          stop 1
        endif
      else
        write(MR_global_info,*)"Start hour of simulation is in range of NWP data"
        prestep  = .false.
      endif
      !   Checking if a poststep is needed
      met_t1  = MR_windfile_starthour(MR_iwindfiles)+MR_windfile_stephour(MR_iwindfiles,nt_fullmet)
      met_dt1 = StepInterval
      if(MR_Comp_StartHour+MR_Comp_Time_in_hours.ge.met_t1)then
        ! Start time requested is after that available, check if we can extrapolate
        if(MR_Comp_StartHour+MR_Comp_Time_in_hours.le.met_t1+met_dt1)then
          write(MR_global_info,*)"WARNING: End time is at or after the last time step."
          write(MR_global_info,*)"         However, it is within one interval so we will"
          write(MR_global_info,*)"         apply the value at MR_windfile_starthour(MR_iwindfiles)"
          if(StepInterval.gt.720.0_dp)then
            if(MR_iwind.eq.1.and.MR_iwindformat.eq.1.and.nt_fullmet.eq.1)then
              write(MR_global_info,*)"         Note: This interval is set high since only one"
              write(MR_global_info,*)"               sonde file was provided."
            endif
          endif
          poststep = .true.
        else
          write(MR_global_error,*)"MR ERROR: End time is after the last available data and"
          write(MR_global_error,*)"       cannot be extrapolated."
          write(MR_global_error,*)"  MR_Comp_StartHour    = ",real(MR_Comp_StartHour,kind=4),&
                                  HS_yyyymmddhhmm_since(MR_Comp_StartHour,MR_BaseYear,MR_useLeap)
          write(MR_global_error,*)"  MR_Comp_Time_in_hours= ",real(MR_Comp_Time_in_hours,kind=4),&
                                  HS_yyyymmddhhmm_since(MR_Comp_Time_in_hours,MR_BaseYear,MR_useLeap)
          write(MR_global_error,*)"  met_t1               = ",real(met_t1,kind=4),&
                                  HS_yyyymmddhhmm_since(met_t1,MR_BaseYear,MR_useLeap)
          write(MR_global_error,*)"  met_dt1              = ",real(met_dt1,kind=4),&
                                  HS_yyyymmddhhmm_since(met_dt1,MR_BaseYear,MR_useLeap)
          stop 1
        endif
      else
        write(MR_global_info,*)"End hour of simulation is in range of NWP data"
        poststep  = .false.
      endif

      ! Loop through all the files and steps and count how many are needed to bracket the time
      ! needed (MR_Comp_StartHour -> MR_Comp_StartHour+MR_Comp_Time_in_hours)
      ! Note: If prestep or poststep is needed, istep will be incremented accordingly
      ! Once we know the number of steps needed, we will allocate space, then fill the variables
      ! with just the step info needed
      write(MR_global_info,*)'    File num  | Step in file  |  stephour   |   SimStartHour   | nMetStep | Note'
      if(prestep)then
        Found_First_Step = .true.
        istep = 1
        stephour = MR_windfile_starthour(1) + MR_windfile_stephour(1,1) - StepInterval
        write(MR_global_info,150)1,1,stephour,MR_Comp_StartHour,istep,"Prestep before MR_Comp_StartHour "
      else
        Found_First_Step = .false.
        istep = 0
      endif
      Found_Last_Step  = .false.
      do iw = 1,MR_iwindfiles
        do iwstep = 1,MR_windfiles_nt_fullmet(iw)
          stephour = MR_windfile_starthour(iw) + MR_windfile_stephour(iw,iwstep)
          ! Unless the start time is before this step hour, reset the index istep to 1
          !  Otherwise, increment index
          if(stephour.lt.MR_Comp_StartHour)then
            Found_First_Step = .false.
            istep = 1
            if(MR_windfiles_nt_fullmet(iw).lt.NT_MAXOUT)then
              ! This suppresses outputting the windfile step info for the files that contain a year's
              ! worth of data (NCEP 2.5, etc)
              write(MR_global_info,150)                    &
                iw,iwstep,stephour,MR_Comp_StartHour,istep,&
                "Before or at MR_Comp_StartHour."
            endif
          else
            !  Otherwise, increment index if we are still in the needed time bracket
            if(.not.Found_Last_Step)then
              Found_First_Step = .true.
              istep = istep + 1
              if(MR_windfiles_nt_fullmet(iw).lt.NT_MAXOUT)then
                write(MR_global_info,150)                    &
                  iw,iwstep,stephour,MR_Comp_StartHour,istep,&
                  "After MR_Comp_StartHour "
              endif
            endif
          endif

          ! Check if the current step hour is exclusively after the end of the simulation
          if(stephour.gt.MR_Comp_StartHour+MR_Comp_Time_in_hours)then
            if(.not.Found_Last_Step)nMetSteps_Comp=istep
            Found_Last_Step = .true.
          endif
          if(Found_Last_Step &
             .and.MR_windfiles_nt_fullmet(iw).lt.NT_MAXOUT)then ! This condition suppresses output for NCEP 2.5
              write(MR_global_info,150)                    &
                iw,iwstep,stephour,MR_Comp_StartHour,istep,&
                "At or after END OF SIM  "
          endif
        enddo
      enddo

      ! If we went through all the steps and didn't find the last step, and if poststep=T, then
      ! increment istep and set nMetSteps_Comp
      if(.not.Found_Last_Step)then
        if(poststep)then
          nMetSteps_Comp = istep+1
          stephour = MR_windfile_starthour(MR_iwindfiles) + &
                       MR_windfile_stephour(MR_iwindfiles,MR_windfiles_nt_fullmet(MR_iwindfiles)) + &
                       StepInterval
          write(MR_global_info,150)MR_iwindfiles,1,stephour,MR_Comp_StartHour,nMetSteps_Comp,"Poststep after END OF SIM  "

        else
          write(MR_global_error,*)"MR ERROR:  Something is wrong.  Could not find the last MetStep needed for"
          write(MR_global_error,*)"           the simulation."
          stop 1
        endif
      endif
 150  format(8x,i3,9x,i4,4x,f15.2,2x,f13.2,10x,i4,8x,a25)

      ! We now have the number of steps needed for the computation
      ! Allocate the lists
      MR_MetSteps_Total = nMetSteps_Comp
      write(MR_global_info,*)"MR: Allocating space for ",MR_MetSteps_Total,"step(s)"
      if(prestep.or.poststep)then
        write(MR_global_info,*)"         Including:"
        if(prestep) write(MR_global_info,*)"             1 prestep"
        if(poststep)write(MR_global_info,*)"             1 poststep"
      endif
      if (.not. allocated(MR_MetStep_File)) then
        allocate(MR_MetStep_File(MR_MetSteps_Total))               ;MR_MetStep_File(:)=''
        allocate(MR_MetStep_findex(MR_MetSteps_Total))             ;MR_MetStep_findex(:)=0
        allocate(MR_MetStep_tindex(MR_MetSteps_Total))             ;MR_MetStep_tindex(:)=0
        allocate(MR_MetStep_Hour_since_baseyear(MR_MetSteps_Total));MR_MetStep_Hour_since_baseyear(:)  =0.0_sp  
        allocate(MR_MetStep_Interval(MR_MetSteps_Total))           ;MR_MetStep_Interval(:)=0  
        allocate(MR_MetStep_year(MR_MetSteps_Total))               ;MR_MetStep_year(:)=0  
        allocate(MR_MetStep_month(MR_MetSteps_Total))              ;MR_MetStep_month(:)=  0  
        allocate(MR_MetStep_day(MR_MetSteps_Total))                ;MR_MetStep_day(:)=0  
        allocate(MR_MetStep_DOY(MR_MetSteps_Total))                ;MR_MetStep_DOY(:)=0  
        allocate(MR_MetStep_Hour_Of_Day(MR_MetSteps_Total))        ;MR_MetStep_Hour_Of_Day(:)=0.0_sp
        allocate(MR_iwind5_year(MR_MetSteps_Total))                ;MR_iwind5_year(:)=0
      endif

      ! Finally, we need to loop through the steps exactly as above, but this time populate 
      ! the lists just allocated
      if(prestep)then
        Found_First_Step = .true.
        istep = 1                                                         ! pre-step is one interval back
        stephour = MR_windfile_starthour(1) + MR_windfile_stephour(1,1) - StepInterval
        MR_MetStep_File(istep)            = trim(adjustl(MR_windfiles(1)))
        MR_MetStep_findex(istep)          = 1
        MR_MetStep_tindex(istep)          = 1
        MR_MetStep_Hour_since_baseyear(istep) = stephour
        MR_MetStep_year(istep)            = HS_YearOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
        MR_MetStep_month(istep)           = HS_MonthOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
        MR_MetStep_day(istep)             = HS_DayOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
        MR_MetStep_DOY(istep)             = HS_DayOfYear(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
        MR_MetStep_Hour_Of_Day(istep)     = real(HS_HourOfDay(real(stephour,kind=dp),MR_BaseYear,MR_useLeap),kind=sp)
        MR_iwind5_year(istep)             = MR_Comp_StartYear
      else
        Found_First_Step = .false.
        istep = 0
      endif
      Found_Last_Step  = .false.
      do iw = 1,MR_iwindfiles
        do iwstep = 1,MR_windfiles_nt_fullmet(iw)
          stephour = MR_windfile_starthour(iw) + MR_windfile_stephour(iw,iwstep)
          ! Unless the start time is before this step hour, reset the index istep to 1
          !  Otherwise, increment index
          if(stephour.lt.MR_Comp_StartHour)then
            Found_First_Step = .false.
            istep = 1
          else
            !  Otherwise, increment index if we are still in the needed time bracket
            if(.not.Found_Last_Step)then
              Found_First_Step = .true.
              istep = istep + 1
            endif
          endif

          if(.not.Found_Last_Step)then
            MR_MetStep_File(istep)            = trim(adjustl(MR_windfiles(iw)))
            MR_MetStep_findex(istep)          = iw
            MR_MetStep_tindex(istep)          = iwstep
            MR_MetStep_Hour_since_baseyear(istep) = stephour
            MR_MetStep_year(istep)            = HS_YearOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
            MR_MetStep_month(istep)           = HS_MonthOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
            MR_MetStep_day(istep)             = HS_DayOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
            MR_MetStep_DOY(istep)             = HS_DayOfYear(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
            MR_MetStep_Hour_Of_Day(istep)     = real(HS_HourOfDay(real(stephour,kind=dp),MR_BaseYear,MR_useLeap),kind=sp)
            MR_iwind5_year(istep)             = MR_MetStep_year(istep)
          endif

          ! Check if the current step hour is exclusively after the end of the simulation
          if(stephour.gt.MR_Comp_StartHour+MR_Comp_Time_in_hours)then
            if(.not.Found_Last_Step)nMetSteps_Comp=istep
            Found_Last_Step = .true.
          endif
        enddo
      enddo
      if(.not.Found_Last_Step)then
        if(poststep)then
          istep = istep+1                                      ! post-step is one interval later
          stephour = MR_MetStep_Hour_since_baseyear(istep-1) + StepInterval
          MR_MetStep_File(istep)            = MR_MetStep_File(istep-1)
          MR_MetStep_findex(istep)          = MR_MetStep_findex(istep-1)
          MR_MetStep_tindex(istep)          = MR_MetStep_tindex(istep-1)
          MR_MetStep_Hour_since_baseyear(istep) = stephour
          MR_MetStep_year(istep)            = HS_YearOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
          MR_MetStep_month(istep)           = HS_MonthOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
          MR_MetStep_day(istep)             = HS_DayOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
          MR_MetStep_DOY(istep)             = HS_DayOfYear(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
          MR_MetStep_Hour_Of_Day(istep)     = real(HS_HourOfDay(real(stephour,kind=dp),MR_BaseYear,MR_useLeap),kind=sp)
          MR_iwind5_year(istep)             = MR_MetStep_year(istep)
        else
          write(MR_global_error,*)"MR ERROR:  Something is wrong.  Could not find the last MetStep needed for"
          write(MR_global_error,*)"           the simulation."
          stop 1
        endif
      endif

      do istep = 1,MR_MetSteps_Total-1
        MR_MetStep_Interval(istep) = MR_MetStep_Hour_since_baseyear(istep+1) - &
                                     MR_MetStep_Hour_since_baseyear(istep)
      enddo
      MR_MetStep_Interval(MR_MetSteps_Total)=MR_MetStep_Interval(MR_MetSteps_Total-1)

      write(MR_global_info,*)"Comp start = ",real(MR_Comp_StartHour,kind=4)
      write(MR_global_info,*)"Comp end   = ",real(MR_Comp_StartHour+MR_Comp_Time_in_hours,kind=4)

      do istep = 1,MR_MetSteps_Total
        iw       = MR_MetStep_findex(istep)
        iwstep   = MR_MetStep_tindex(istep)
        stephour = MR_MetStep_Hour_since_baseyear(istep)
      enddo

!     MAKE SURE THE WIND MODEL TIME WINDOW COVERS THE ENTIRE SIMULATION TIME
      write(MR_global_info,99)
99    format(/,4x,'Making sure the mesoscale model time covers the simulation time . . . ')
      if (MR_MetStep_Hour_since_baseyear(1).gt.MR_Comp_StartHour) then
        write(MR_global_error,*)"MR ERROR:  Wind data starts after SimStartHour"
        write(MR_global_error,*)"WindHour(MR_iwindfiles) = ",MR_MetStep_Hour_since_baseyear(MR_MetSteps_Total)
        write(MR_global_error,*)"SimStartHour            = ",MR_Comp_StartHour
        write(MR_global_error,*)"Simtime_in_hours        = ",MR_Comp_Time_in_hours
        write(MR_global_error,*)"  All steps:"
        do i = 1,MR_MetSteps_Total
          write(MR_global_error,*)"    ",i,MR_MetStep_Hour_since_baseyear(i)
        enddo
        !write(MR_global_log ,101)
!101     format(4x,'Error.  The model starts after the first eruption. Program stopped.')
        stop 1
      elseif (MR_MetStep_Hour_since_baseyear(MR_MetSteps_Total).lt.(MR_Comp_StartHour+MR_Comp_Time_in_hours)) then
        write(MR_global_error,*)"MR ERROR:  Last met time is earlier than simulation time"
        write(MR_global_error,*)"MR_MetSteps_Total     = ",MR_MetSteps_Total
        write(MR_global_error,*)"Last time step       = ",MR_MetStep_Hour_since_baseyear(MR_MetSteps_Total)
        write(MR_global_error,*)"SimStartHour         = ",MR_Comp_StartHour
        write(MR_global_error,*)"Simtime_in_hours     = ",MR_Comp_Time_in_hours
        do i = 1,MR_MetSteps_Total
          write(MR_global_error,*)"    ",i,MR_MetStep_Hour_since_baseyear(i)
        enddo
        write(MR_global_error ,102)
102     format(4x,'Error.  The model ends before the end of the last eruption.  Program stopped.')
        stop 1
      endif
      write(MR_global_info,103)
103   format(4x,'Good.  It does.',/)

      CALLED_MR_Set_Met_Times = .true.

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      end subroutine MR_Set_Met_Times


!##############################################################################
!
!     MR_Read_HGT_arrays
!
!     This subroutine basically does the same thing as Read_3d_MetP_Variable
!     but specifically for the variable HGT (ivar=1).  This is needed because the
!     more general Read_3d_MetP_Variable just returns data in dum3d_metP and
!     expects the calling program to save the results.  The HGT variables,
!     geoH_metP_last and geoH_metP_next need to persist locally, however,  since they are
!     used in QC calculations and are needed in converting PresVertVel (from Pa
!     s to m/s).
!     The values extracted are just on the needed subgrid of the full met
!     grid on pressure coordinates.
!
!     Takes as input :: istep :: specified the met step
!
!     Sets: geoH_metP_last, geoH_metP_next
!
!##############################################################################

      subroutine MR_Read_HGT_arrays(istep,reset_first_time)

      implicit none

      integer,intent(in)           :: istep
      logical, optional,intent(in) :: reset_first_time

      integer :: ivar
      logical,save :: first_time = .true.
      integer :: k

      if(present(reset_first_time)) then
        first_time = .true.
      endif

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Read_HGT_arrays                                    ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      write(MR_global_production,*)"     Reading HGT array for istep = ",istep

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .false., &  ! CALLED_MR_Set_CompProjection
                                 .false., &  ! CALLED_MR_Initialize_Met_Grids
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)
      ivar = 1 ! HGT
      if(first_time)then
        if(MR_idataFormat.eq.1)then
          call MR_Read_MetP_Variable_ASCII_1d(ivar,istep)
        elseif(MR_idataFormat.eq.2)then
#ifdef USENETCDF
          call MR_Read_MetP_Variable_netcdf(ivar,istep)
#endif
        elseif(MR_idataFormat.eq.3)then
#ifdef USEGRIB
          call MR_Read_MetP_Variable_GRIB(ivar,istep)
#endif
        endif
        MR_geoH_metP_last(:,:,:) = MR_dum3d_metP(:,:,:)
        first_time = .false.
        Max_geoH_metP_last = maxval(MR_geoH_metP_last(:,:,np_fullmet))
      else
        MR_geoH_metP_last(:,:,:) = MR_geoH_metP_next(:,:,:)
        Max_geoH_metP_last = Max_geoH_metP_next
        Max_geoH_metP_next = maxval(MR_geoH_metP_next(:,:,np_fullmet))
      endif
      if(MR_idataFormat.eq.1)then
          call MR_Read_MetP_Variable_ASCII_1d(ivar,istep+1)
      elseif(MR_idataFormat.eq.2)then
#ifdef USENETCDF
        call MR_Read_MetP_Variable_netcdf(ivar,istep+1)
#endif
      elseif(MR_idataFormat.eq.3)then
#ifdef USEGRIB
        call MR_Read_MetP_Variable_GRIB(ivar,istep+1)
#endif
      endif

      MR_geoH_metP_next(:,:,:) = MR_dum3d_metP(:,:,:)

      Max_geoH_metP_next = maxval(MR_geoH_metP_next(:,:,np_fullmet))
        ! Now determine the maximum value provided between last and next steps
      Max_geoH_metP = max(Max_geoH_metP_last,Max_geoH_metP_next)
        ! and compare with the maximum value needed by the computational grid
      if(Max_geoH_metP.lt.z_comp_sp(nz_comp))then
        ! use highest needed point (rounded up to 5-km increment) for the
        ! padding height
        Suppl_H = real(ceiling(z_comp_sp(nz_comp)/5.0_sp),kind=sp) * 5.0_sp
      else
        ! otherwise use the height from the met files to determine padding
        Suppl_H = real(ceiling(Max_geoH_metP/5.0_sp),kind=sp) * 5.0_sp
      endif
      if(MR_iHeightHandler.eq.1)then
        if(z_comp_sp(nz_comp).gt.Max_geoH_metP)then
          write(MR_global_error,*)"MR ERROR : Computational grid extends higher than met grid"
          write(MR_global_error,*)"           iHeightHandler = 1:  exiting"
          stop 1
        endif
      endif

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      end subroutine MR_Read_HGT_arrays


!##############################################################################
!
!     MR_Read_3d_MetP_Variable
!
!     This subroutine extracts the variable Met_var_NC_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid on pressure coordinates.
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Sets: dum3d_metP
!##############################################################################

      subroutine MR_Read_3d_MetP_Variable(ivar,istep)

      implicit none

      integer,intent(in)        :: ivar
      integer,intent(in)        :: istep

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Read_3d_MetP_Variable                              ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .false., &  ! CALLED_MR_Set_CompProjection
                                 .false., &  ! CALLED_MR_Initialize_Met_Grids
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      select case (MR_iwind)
      case(1)   ! if we're using a 1-D wind sounding
        call MR_Read_MetP_Variable_ASCII_1d(ivar,istep)
      case(2)
        !call Read_3d_ASCII
      case (3:5)
        if(MR_idataFormat.eq.2)then
#ifdef USENETCDF
          call MR_Read_MetP_Variable_netcdf(ivar,istep)
#endif
        elseif(MR_idataFormat.eq.3)then
#ifdef USEGRIB
          call MR_Read_MetP_Variable_GRIB(ivar,istep)
#endif
        endif

      case default

      end select

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      end subroutine MR_Read_3d_MetP_Variable

!##############################################################################
!
!     MR_Read_3d_MetH_Variable
!
!     This subroutine extracts the variable Met_var_NC_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on height coordinates.
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Reads : MR_dum3d_metP
!     Sets  : MR_dum3d_metH
!##############################################################################

      subroutine MR_Read_3d_MetH_Variable(ivar,istep)

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer,intent(in)        :: ivar
      integer,intent(in)        :: istep

      real(kind=sp),dimension(:),allocatable :: z_col_metP
      real(kind=sp),dimension(:),allocatable :: var_col_metP
      real(kind=sp),dimension(:),allocatable :: var_col_metH
      integer :: i,j,k
      integer :: kc,knext
      integer :: np_fully_padded
      real(kind=sp),dimension(:),allocatable :: dumVertCoord_sp

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Read_3d_MetH_Variable                              ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      np_fully_padded = np_fullmet+1+np_fullmet_pad

      ! First get the variable on the native pressure coordinate
      call MR_Read_3d_MetP_Variable(ivar,istep)

      ! Now remap from
      !    nx_submet,ny_submet,np_fullmet to
      !    nx_submet,ny_submet,nz_comp
      ! Allocate 1-d arrays, padding one layer to the bottom and np_fullmet_pad
      ! layers on the top
      allocate(  z_col_metP(np_fully_padded));      z_col_metP(:) = 0.0_sp
      allocate(var_col_metP(np_fully_padded));    var_col_metP(:) = 0.0_sp
      allocate(var_col_metH(nz_comp))        ;    var_col_metH(:) = 0.0_sp
      allocate(dumVertCoord_sp(nz_comp))     ; dumVertCoord_sp(:) = 0.0_sp

!     CREATE 1-D ARRAYS IN P, AND REGRID THEM INTO 1-D ARRAYS IN z
      do i=1,nx_submet
        do j=1,ny_submet
            ! copy the column of z values for this i,j
            ! Note: These are the height values from the windfile plus an extra
            ! point at the surface and an extra point above the wind grid.
          if(istep.eq.MR_iMetStep_Now)then
            z_col_metP(2:np_fullmet+1)  = MR_geoH_metP_last(i,j,1:np_fullmet)
          else
            z_col_metP(2:np_fullmet+1)  = MR_geoH_metP_next(i,j,1:np_fullmet)
          endif
          ! Loop through the nodes and reset any that are negative or non-ascending in p
          ! This is usually never invoked, but would only affect potentially the bottom few nodes
          ! The intent is to maintain a monotonic GPH, while getting irrelevant nodes out of the
          ! way.
          do k=1,np_fullmet
            if(z_col_metP(k).le.real(k-1,kind=sp)*MR_MIN_DZ) z_col_metP(k)=real(k-1,kind=sp)*MR_MIN_DZ
          enddo
          var_col_metP(1)              = MR_dum3d_metP(i,j,1)
          var_col_metP(2:np_fullmet+1) = MR_dum3d_metP(i,j,1:np_fullmet)
          ! Fix occasional erroneous pressure levels that are non-increasing by
          ! resetting height to average of the neighbors in z, or, if more than
          ! one z value is lower than a previous one, reset all such values to
          ! increase linearly with k
          do k=2,np_fullmet
            if (z_col_metP(k)<z_col_metP(k-1))then
              knext=k       !find the first value of z_col_metP that's above z_col_metP(k-1)
              do while ((z_col_metP(knext)<z_col_metP(k-1)).and. &
                        (knext.lt.np_fullmet+1))
                 knext=knext+1
              enddo
              do kc=k,knext-1    !correct all intervening values of z_col_metP
                 z_col_metP(kc)=z_col_metP(k-1)+(real(kc-(k-1),kind=sp)/real(knext-(k-1),kind=sp))* &
                                              (z_col_metP(knext)-z_col_metP(k-1))
              enddo
            endif
          enddo
          if(np_fullmet_pad.gt.0)then
            ! Assign top node
            if(MR_iHeightHandler.eq.1)then
              z_col_metP(np_fullmet+2) = z_col_metP(np_fullmet+1) + 10.0_sp
              write(MR_global_info,*)z_col_metP(np_fullmet+1:np_fullmet+2)
            elseif(MR_iHeightHandler.eq.2)then
              z_col_metP(np_fullmet+2) = Suppl_H
              if(ivar.ne.5)then
                ! For all variables except temperature (ivar=5), just copy highest
                ! met node
                var_col_metP(np_fullmet+2) = var_col_metP(np_fullmet+1)
              else
                ! For temperature, use the 1976 US Standard Atmosphere
                var_col_metP(np_fullmet+2) = MR_Temp_US_StdAtm(Suppl_H)
              endif
            endif
          endif

          if (MR_use_SigmaAlt)then
              ! this recovers the real-world z coordinate from the sigma level and topography
            dumVertCoord_sp(1:nz_comp) = MR_zsurf(i,j) + &
                                           s_comp_sp(1:nz_comp) * MR_jacob(i,j)
          else
            dumVertCoord_sp(1:nz_comp) = z_comp_sp(1:nz_comp)
          endif

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !   Interpolate these values to a regular grid with
          !   spacing equal to the simulation grid
          call MR_Regrid_P2H_linear(np_fullmet+2, z_col_metP,       var_col_metP, & 
                                    nz_comp,      dumVertCoord_sp,  var_col_metH)

          MR_dum3d_metH(i,j,:) = var_col_metH

        enddo ! j
      enddo  ! i

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      return

      end subroutine MR_Read_3d_MetH_Variable

!##############################################################################
!
!     MR_Read_3d_Met_Variable_to_CompGrid
!      Note: this should be named MR_Read_3d_Met_Variable_to_CompH since
!            these is now an analogous subroutine mapping MetP to CompP
!
!     This subroutine extracts the variable Met_var_NC_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on to the computational grid.  This is done by calling:
!       Read_3d_MetH_Variable                   (gets variable on Met_x, Met_y, Comp_z)
!        -> Read_3d_MetP_Variable               (gets variable on native Met subgrid)
!            -> Read_3d_MetP_Variable_[format]  (direct read of variable in
!                                                whatever format : nc,grib1/2,ascii)
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Sets  : MR_dum3d_compH
!               MR_dum3d_metH and MR_dum3d_metP are also filled in the course of
!               generating MR_dum3d_compH
!
!##############################################################################

      subroutine MR_Read_3d_Met_Variable_to_CompGrid(ivar,istep,IsNext)

      implicit none

      integer,intent(in)           :: ivar
      integer,intent(in)           :: istep
      logical, optional,intent(in) :: IsNext

      integer             :: i,j,k
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Read_3d_Met_Variable_to_CompGrid                   ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

        ! First get the variable on the height coordinate
      call MR_Read_3d_MetH_Variable(ivar,istep)

      if(MR_Save_Velocities)then
        ! Note: this flag for saving the velocity values is useful in special
        ! cased such where the velocities might be read and used for a local
        ! calculation, but then can be used later.  Variable diffusivity uses
        ! this.
        if(present(IsNext)) then
          ! MR_dum3d_metP still contains the variable just read
          if(ivar.eq.2)then
            if(IsNext)then
              MR_vx_metP_last = MR_vx_metP_next
              MR_vx_metP_next = MR_dum3d_metP
            else
              MR_vx_metP_last = MR_dum3d_metP
            endif
          elseif(ivar.eq.3)then
            if(IsNext)then
              MR_vy_metP_last = MR_vy_metP_next
              MR_vy_metP_next = MR_dum3d_metP
            else
              MR_vy_metP_last = MR_dum3d_metP
            endif
          endif
        endif
      endif
        ! Now we have MR_dum3d_metH; interpolate onto computational grid
        !  Since MR_dum3d_metH and MR_dum3d_compH have the same z-coordinate, we only
        !  need to do 2d regridding on each k-slice
      allocate(tmp_regrid2d_sp(nx_comp,ny_comp)); tmp_regrid2d_sp(:,:)=0.0_sp

      do k=1,nz_comp
        if(MR_iwindformat.eq.1.or.MR_iwindformat.eq.2)then
          !NOTE: This will not work for multi-site sonde data
          MR_dum3d_compH(:,:,k) = MR_dum3d_metH(1,1,k)
          cycle
        endif
  
        call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum3d_metH(1:nx_submet,1:ny_submet,k),       &
                                nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
  
        do i = 1,nx_comp
          do j = 1,ny_comp
            if(isnan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
          enddo
        enddo
        MR_dum3d_compH(:,:,k) = tmp_regrid2d_sp(:,:)
      enddo

      if (MR_use_SigmaAlt) then
        ! If we are using sigma-altitude coordinates, then MR_dum3d_compH is returned on the
        ! s-grid, but we might need to scale the variable
        if(ivar.eq.4)then
          ! vertical velocity is scaled by jacobian
          do i=1,nx_comp
            do j=1,ny_comp
              MR_dum3d_compH(i,j,:)=MR_dum3d_compH(i,j,:)/MR_jacob(i,j)
            enddo
          enddo
        endif
      endif

      deallocate(tmp_regrid2d_sp)

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      return

      end subroutine MR_Read_3d_Met_Variable_to_CompGrid

!##############################################################################
!
!     MR_Read_3d_Met_Variable_to_CompP
!
!     This subroutine extracts the variable Met_var_NC_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on to the computational grid along pressure levels.
!     This is done by calling:
!       Read_3d_MetP_Variable                   (gets variable on Met_x,Met_y, Met_p)
!        -> Read_3d_MetP_Variable               (gets variable on native Met subgrid)
!            -> Read_3d_MetP_Variable_[format]  (direct read of variable in
!                                                whatever format :
!                                                nc,grib1/2,ascii)
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Sets  : MR_dum3d_compP
!               MR_dum3d_metP is also filled in the course of
!               generating MR_dum3d_compP
!
!##############################################################################

      subroutine MR_Read_3d_Met_Variable_to_CompP(ivar,istep,IsNext)

      implicit none

      integer,intent(in)           :: ivar
      integer,intent(in)           :: istep
      logical, optional,intent(in) :: IsNext

      integer             :: i,j,k
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Read_3d_Met_Variable_to_CompP                      ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times(this check is needed)

      if(.not.MR_useCompP)then
        write(MR_global_error,*)"MR ERROR: Trying to run a CompP subroutine"
        write(MR_global_error,*)"          with MR_useCompP = .false."
        stop 0
      endif

        ! First get the variable on the pressure levels
      call MR_Read_3d_MetP_Variable(ivar,istep)

      if(MR_Save_Velocities)then
        ! Note: this flag for saving the velocity values is useful in special
        ! cased such where the velocities might be read and used for a local
        ! calculation, but then can be used later.  Variable diffusivity uses
        ! this.
        if(present(IsNext)) then
          ! MR_dum3d_metP still contains the variable just read
          if(ivar.eq.2)then
            if(IsNext)then
              MR_vx_metP_last = MR_vx_metP_next
              MR_vx_metP_next = MR_dum3d_metP
            else
              MR_vx_metP_last = MR_dum3d_metP
            endif
          elseif(ivar.eq.3)then
            if(IsNext)then
              MR_vy_metP_last = MR_vy_metP_next
              MR_vy_metP_next = MR_dum3d_metP
            else
              MR_vy_metP_last = MR_dum3d_metP
            endif
          endif
        endif
      endif
        ! Now we have MR_dum3d_metP; interpolate onto computational grid
        !  Since MR_dum3d_metP and MR_dum3d_compP have the same
        !  pressure level, we only
        !  need to do 2d regridding on each p-slice
      allocate(tmp_regrid2d_sp(nx_comp,ny_comp));tmp_regrid2d_sp(:,:)=0.0_sp

      do k=1,np_fullmet
        if(MR_iwindformat.eq.1.or.MR_iwindformat.eq.2)then
          !NOTE: This will not work for multi-site sonde data
          MR_dum3d_compP(:,:,k) = MR_dum3d_metP(1,1,k)
          cycle
        endif

        call MR_Regrid_Met2Comp(nx_submet,ny_submet,MR_dum3d_metP(1:nx_submet,1:ny_submet,k),       &
                                nx_comp,  ny_comp,tmp_regrid2d_sp(1:nx_comp,1:ny_comp))

        do i = 1,nx_comp
          do j = 1,ny_comp
            if(isnan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
          enddo
        enddo
        MR_dum3d_compP(:,:,k) = tmp_regrid2d_sp(:,:)
      enddo

      deallocate(tmp_regrid2d_sp)

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      return

      end subroutine MR_Read_3d_Met_Variable_to_CompP

!##############################################################################
!
!     MR_Read_2d_Met_Variable
!
!     This subroutine extracts the variable Met_var_NC_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on height coordinates.
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Sets  : MR_dum2d_met or MR_dum2d_met_int
!##############################################################################

      subroutine MR_Read_2d_Met_Variable(ivar,istep)

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer,intent(in)        :: ivar
      integer,intent(in)        :: istep   

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Read_2d_Met_Variable                               ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .false., &  ! CALLED_MR_Set_CompProjection
                                 .false., &  ! CALLED_MR_Initialize_Met_Grids
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      select case (MR_iwind)
      case(1)   ! if we're using a 1-D wind sounding
        !call Read_1d_windfile
      case(2)
        !call Read_2d_ASCII
      case (3:5)
        if(MR_idataFormat.eq.2)then
#ifdef USENETCDF
          call MR_Read_MetP_Variable_netcdf(ivar,istep)
#endif
        elseif(MR_idataFormat.eq.3)then
#ifdef USEGRIB
          call MR_Read_MetP_Variable_GRIB(ivar,istep)
#endif
        endif
      case default

      end select

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      end subroutine MR_Read_2d_Met_Variable

!##############################################################################
!
!     MR_Read_2d_Met_Variable_to_CompGrid
!
!     This subroutine extracts the variable Met_var_NC_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on to the computational grid.  This is done by calling:
!       Read_2d_Met_Variable                   (gets variable on Met_x, Met_y)
!          -> Read_Met_Variable_[format]       (direct read of variable in
!                                                whatever format :
!                                                nc,grib1/2,ascii)
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Sets  : MR_dum2d_comp
!               MR_dum2d_met and MR_dum2d_met are also filled in the course of
!               generating MR_dum2d_comp
!
!##############################################################################

      subroutine MR_Read_2d_Met_Variable_to_CompGrid(ivar,istep)

      implicit none

      integer,intent(in)        :: ivar
      integer,intent(in)        :: istep

      integer             :: i,j
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Read_2d_Met_Variable_to_CompGrid                   ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      call MR_Read_2d_Met_Variable(ivar,istep)

      allocate(tmp_regrid2d_sp(nx_comp,ny_comp)); tmp_regrid2d_sp(:,:)=0.0_sp

        call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum2d_met(1:nx_submet,1:ny_submet),       &
                                nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))

      do i = 1,nx_comp
        do j = 1,ny_comp
          if(isnan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
        enddo
      enddo
      MR_dum2d_comp(:,:) = tmp_regrid2d_sp(:,:)

      deallocate(tmp_regrid2d_sp)

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      return

      end subroutine MR_Read_2d_Met_Variable_to_CompGrid

!##############################################################################
!
!     MR_Rotate_UV_GR2ER_Met
!
!     This subroutine reads U and V on the metP grid, rotates to EarthRelative
!
!     Takes as input :: istep   :: specified the met step
!                       SetComp :: (optional) logical flag to regrid to comp grid
!
!     Sets  : MR_dum3d_compH    holds U
!             MR_dum3d_compH_2  holds V
!
!##############################################################################

      subroutine MR_Rotate_UV_GR2ER_Met(istep,SetComp)

      implicit none

      integer,intent(in)  :: istep
      logical,optional,intent(in) :: SetComp

      integer             :: i,j,k

      real(kind=sp) :: vx_old,vy_old
      real(kind=sp) :: vx_new,vy_new
      real(kind=sp) :: rotang

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Rotate_UV_GR2ER_Met                                ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      ! Rotate wind vectors on the MetP grid for NARR cases with Map_Case = 2 (both projected)
      ! or rotate to LL on MetP grid for Map_Case = 4 (Met=proj, Comp=LL)
      !write(MR_global_info,*)" Rotating Earth-relative projected Met winds to grid-relative"
      !write(MR_global_info,*)"  or rotating projected grid-relative Met winds to Earth-relative"

      !if(.not.MR_Save_Velocities)then
      !  write(MR_global_info,*)"MR WARNING: Velocities not saved"
      !endif

      call MR_Read_3d_MetP_Variable(2,istep)
        MR_u_ER_metP = MR_dum3d_metP
      call MR_Read_3d_MetP_Variable(3,istep)
        MR_v_ER_metP = MR_dum3d_metP
      ! Using these earth relative velocities, rotate to grid relative and
      ! copy to MR_dum3d_metP
      do i=1,nx_submet
        do j=1,ny_submet
          ! The angle theta for the Earth to Grid conversion was
          ! precalculated in Set_MetComp_Grids_netcdf
          rotang = real(theta_Met(i,j),kind=sp)
          do k=1,np_fullmet
            vx_old = MR_u_ER_metP(i,j,k)
            vy_old = MR_v_ER_metP(i,j,k)
          ! Project vx_old onto grid
            vx_new = vx_old * cos(rotang)
            vy_new = vx_old * sin(rotang)
          ! Add projection of vel_old
            vx_new = vx_new - vy_old * sin(rotang)
            vy_new = vy_new + vy_old * cos(rotang)
            MR_u_ER_metP(i,j,k) = vx_new
            MR_v_ER_metP(i,j,k) = vy_new
          enddo
        enddo
      enddo

      if(present(SetComp)) then
        if(SetComp)then
          if(MR_useCompH)then
            ! copy v to main workspace and regrid
            MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_fullmet) = &
              MR_v_ER_metP(1:nx_submet,1:ny_submet,1:np_fullmet)
            call MR_Regrid_MetP_to_CompGrid(istep)
            MR_dum3d_compH_2(1:nx_comp,1:ny_comp,1:nz_comp) = &
              MR_dum3d_compH(1:nx_comp,1:ny_comp,1:nz_comp)
  
            ! copy u to main workspace and regrid
            MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_fullmet) = &
              MR_u_ER_metP(1:nx_submet,1:ny_submet,1:np_fullmet)
            call MR_Regrid_MetP_to_CompGrid(istep)
            !MR_dum3d_compH is set by this regrid call
          endif

          if(MR_useCompP)then
            ! copy v to main workspace and regrid
            MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_fullmet) = &
              MR_v_ER_metP(1:nx_submet,1:ny_submet,1:np_fullmet)
            call MR_Regrid_MetP_to_CompGrid(istep)
            MR_dum3d_compP_2(1:nx_comp,1:ny_comp,1:np_fullmet) = &
              MR_dum3d_compP(1:nx_comp,1:ny_comp,1:np_fullmet)
            ! copy u to main workspace and regrid
            MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_fullmet) = &
              MR_u_ER_metP(1:nx_submet,1:ny_submet,1:np_fullmet)
            call MR_Regrid_MetP_to_CompGrid(istep)
            !MR_dum3d_compP is set by this regrid call
          endif

        endif
      endif

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      end subroutine MR_Rotate_UV_GR2ER_Met


!##############################################################################
!
!     MR_Rotate_UV_ER2GR_Comp
!
!     This subroutine is needed when we have a Earth-relative Met grid (either
!     natively or through MR_Rotate_UV_GR2ER_Met) and need to rotate that into a
!     projected Comp grid.  This corresponds to:
!       Map_Case = 3 : Met=LL,   Comp=proj
!       Map_Case = 5 : Met=proj, Comp=proj (different)
!     This subroutine reads U and V on the metP projected grid,
!     then resamples onto the projected compH grid, rotates to grid relative and
!     returns to calling program
!
!     Takes as input :: istep :: specified the met step
!
!     Sets  : MR_dum3d_compH    holds U
!             MR_dum3d_compH_2  holds V
!
!##############################################################################

      subroutine MR_Rotate_UV_ER2GR_Comp(istep)

      implicit none

      integer,intent(in)  :: istep

      integer             :: i,j,k

      real(kind=sp) :: vx_old,vy_old
      real(kind=sp) :: vx_new,vy_new
      real(kind=sp) :: rotang1

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Rotate_UV_ER2GR_Met                                ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      if(Map_Case.eq.3.or.Map_Case.eq.4)then
        ! Met grid is natively LL and Comp grid is projected
        !  - or -
        ! Met grid is projected and comp grid is LL
        ! Fill the y velocities to spare array (comp_2)
        call MR_Read_3d_MetP_Variable(3,istep) ! This fills MR_dum3d_metP
        call MR_Regrid_MetP_to_CompGrid(istep) ! Takes MR_dum3d_metP and fills MR_dum3d_compH
          MR_dum3d_compH_2 = MR_dum3d_compH

        call MR_Read_3d_MetP_Variable(2,istep) 
        call MR_Regrid_MetP_to_CompGrid(istep)
        ! Now compH and compH_2 have vx and vy
      elseif(Map_Case.eq.5)then
        ! If the velocity components are grid relative, then
        ! assume MR_u_ER_metP and MR_v_ER_metP have already been set by
        ! MR_Rotate_UV_GR2ER_Met
        MR_dum3d_metP = MR_v_ER_metP
        call MR_Regrid_MetP_to_CompGrid(istep)
          MR_dum3d_compH_2 = MR_dum3d_compH
        MR_dum3d_metP = MR_u_ER_metP
        call MR_Regrid_MetP_to_CompGrid(istep)
        ! Now compH and compH_2 have vx and vy
      else
        write(MR_global_error,*)"Should not be calling this subroutine unless Map_Case=3 or 5"
        stop 1
      endif

      ! Now loop through the comp points and rotate the vectors in place
      do i=1,nx_comp
        do j=1,ny_comp
          ! The angle theta for the Earth to Grid conversion was
          ! precalculated in Set_MetComp_Grids_netcdf
          rotang1 = real(theta_Comp(i,j),kind=sp)

          do k=1,nz_comp
            vx_old = MR_dum3d_compH(i,j,k)
            vy_old = MR_dum3d_compH_2(i,j,k)

            vx_new = vx_old * cos(rotang1) - vy_old * sin(rotang1)
            vy_new = vx_old * sin(rotang1) + vy_old * cos(rotang1)

            MR_dum3d_compH(i,j,k) = vx_new
            MR_dum3d_compH_2(i,j,k) = vy_new
          enddo
        enddo
      enddo

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      end subroutine MR_Rotate_UV_ER2GR_Comp

!##############################################################################
!
!     MR_Regrid_MetP_to_CompGrid
!
!     This subroutine expects the calling program to populate MR_dum3d_metP.
!     This is regridded onto MR_dum3d_metH by the subroutine Regrid_MetP_to_MetH.
!     MR_dum3d_metH is then regridded onto MR_dum3d_compH
!
!     Takes as input :: istep :: specified the met step
!     Sets  : MR_dum3d_compH and MR_dum3d_compP
!
!##############################################################################

      subroutine MR_Regrid_MetP_to_CompGrid(istep)

      implicit none

      integer,intent(in)  :: istep

      integer             :: i,j,k
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Regrid_MetP_to_CompGrid                            ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif
      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      if(MR_useCompH)then
        ! convert MR_dum3d_MetP to MR_dum3d_metH
        call MR_Regrid_MetP_to_MetH(istep)
      endif

        ! Now we have MR_dum3d_metH; interpolate onto computational grid
        !  Since MR_dum3d_metH and MR_dum3d_compH have the same z-coordinate, we only
        !  need to do 2d regridding on each k-slice or each p-slice
      allocate(tmp_regrid2d_sp(nx_comp,ny_comp));tmp_regrid2d_sp(:,:)=0.0_sp
      if(MR_useCompH)then
        ! Regridding to CompH
        do k=1,nz_comp
          if(MR_iwindformat.eq.1.or.MR_iwindformat.eq.2)then
            !NOTE: This will not work for multi-site sonde data
            MR_dum3d_compH(1:nx_comp,1:ny_comp,k) = MR_dum3d_metH(1,1,k)
            cycle
          endif
          call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum3d_metH(1:nx_submet,1:ny_submet,k),       &
                                  nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
          do i = 1,nx_comp
            do j = 1,ny_comp
              if(isnan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
            enddo
          enddo
          MR_dum3d_compH(:,:,k) = tmp_regrid2d_sp(:,:)
        enddo
      endif

      if(MR_useCompP)then
        ! Regridding to CompP
        do k=1,np_fullmet
          if(MR_iwindformat.eq.1.or.MR_iwindformat.eq.2)then
            !NOTE: This will not work for multi-site sonde data
            MR_dum3d_compP(:,:,k) = MR_dum3d_metP(1,1,k)
            cycle
          endif
          call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum3d_metP(1:nx_submet,1:ny_submet,k),       &
                                  nx_comp,  ny_comp, tmp_regrid2d_sp(1:nx_comp,1:ny_comp))
          do i = 1,nx_comp
            do j = 1,ny_comp
              if(isnan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
            enddo
          enddo
          MR_dum3d_compP(:,:,k) = tmp_regrid2d_sp(:,:)
        enddo
      endif

      deallocate(tmp_regrid2d_sp)

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      return

      end subroutine MR_Regrid_MetP_to_CompGrid

!##############################################################################
!
!     MR_Regrid_MetP_to_MetH
!
!     This subroutine expects the calling program to populate MR_dum3d_metP.
!     This is regridded onto MR_dum3d_metH .
!
!     Takes as input :: istep :: specified the met step
!     Sets  : MR_dum3d_metH
!
!##############################################################################

      subroutine MR_Regrid_MetP_to_MetH(istep)

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer,intent(in)  :: istep

      real(kind=sp),dimension(:),allocatable :: z_col_metP
      real(kind=sp),dimension(:),allocatable :: var_col_metP
      real(kind=sp),dimension(:),allocatable :: var_col_metH
      integer :: i,j,k
      integer :: kc,knext

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Regrid_MetP_to_MetH                                ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

      ! Now remap from
      !    nx_submet,ny_submet,np_fullmet to
      !    nx_submet,ny_submet,nz_comp

      allocate(  z_col_metP(np_fullmet+2));   z_col_metP(:)=0.0_sp
      allocate(var_col_metP(np_fullmet+2)); var_col_metP(:)=0.0_sp
      allocate(var_col_metH(nz_comp))     ; var_col_metH(:)=0.0_sp

!     CREATE 1-D ARRAYS IN P, AND REGRID THEM INTO 1-D ARRAYS IN z
      do i=1,nx_submet
        do j=1,ny_submet
            ! copy the column of z values for this i,j
            ! Note: These are the height values from the windfile plus an extra
            ! point at the surface and an extra point above the wind grid.
          if(istep.eq.MR_iMetStep_Now)then
            z_col_metP(2:np_fullmet+1)  = MR_geoH_metP_last(i,j,1:np_fullmet)
          else
            z_col_metP(2:np_fullmet+1)  = MR_geoH_metP_next(i,j,1:np_fullmet)
          endif
          ! Loop through the nodes and reset any that are negative or non-acsending in p
          ! This is usually never invoked, but would only affect potentially the bottom few nodes
          ! The intent is to maintain a monotonic GPH, while getting irrelevant nodes out of the
          ! way. 
          do k=1,np_fullmet
            if(z_col_metP(k).le.real(k,kind=sp)*MR_MIN_DZ) z_col_metP(k)=real(k,kind=sp)*MR_MIN_DZ
          enddo
          var_col_metP(1)              = 0.0_sp
          var_col_metP(2:np_fullmet+1) = MR_dum3d_metP(i,j,1:np_fullmet)
          ! Fix occasional erroneous pressure levels that are non-increasing by
          ! resetting height to average of the neighbors in z, or, if more than
          ! one z value is lower than a previous one, reset all such values to
          ! increase linearly with k
          do k=2,np_fullmet
            if (z_col_metP(k)<z_col_metP(k-1))then
              write(MR_global_info,*)"MR_WARNING: Need to adjust value in MR_Regrid_MetP_to_MetH"
              write(MR_global_info,*)"            that are not increasing in height with decreasing p"
              write(MR_global_info,*)"i = ",i
              write(MR_global_info,*)"j = ",j
              write(MR_global_info,*)"k = ",k
              write(MR_global_info,*)z_col_metP
              knext=k       !find the first value of z_col_metP that's above z_col_metP(k-1)
              do while (z_col_metP(knext)<z_col_metP(k-1))
                 knext=knext+1
              enddo
              do kc=k,knext-1    !correct all intervening values of z_col_metP
                 z_col_metP(kc)=z_col_metP(k-1)+(real(kc-(k-1),kind=sp)/real(knext-(k-1),kind=sp))* &
                                              (z_col_metP(knext)-z_col_metP(k-1))
              enddo
            endif
          enddo
            ! Assign top node (Should these indicies be +1?)
          z_col_metP(np_fullmet+2) = max(z_comp_sp(nz_comp)+1.0_sp,z_col_metP(np_fullmet+1)+1.0_sp)
          var_col_metP(np_fullmet+2) = var_col_metP(np_fullmet+1)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !   Interpolate these values to a regular grid with
          !   spacing equal to the simulation grid
          call MR_Regrid_P2H_linear(np_fullmet+2, z_col_metP,  var_col_metP, & !vx
                                    nz_comp,       z_comp_sp,  var_col_metH)
          MR_dum3d_metH(i,j,:) = var_col_metH
        enddo ! j
      enddo  ! i

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      return

      end subroutine MR_Regrid_MetP_to_MetH


!##############################################################################
!
!     MR_Regrid_Met2d_to_Comp2d
!
!     This subroutine expects the calling program to populate MR_dum2d_met.
!     MR_dum2d_met is then regridded onto MR_dum2d_comp
!
!     Takes as input :: istep :: specified the met step
!     Sets  : MR_dum3d_compH
!
!##############################################################################

      subroutine MR_Regrid_Met2d_to_Comp2d

      implicit none

      integer             :: i,j
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Regrid_Met2d_to_Comp2d                             ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! Check prerequisites
      if(Check_prereq_conditions.eqv..true.) &
      call MR_Check_Prerequsites(.true.,  &  ! CALLED_MR_Allocate_FullMetFileList    (this check is needed)
                                 .true.,  &  ! CALLED_MR_Read_Met_DimVars            (this check is needed)
                                 .true.,  &  ! CALLED_MR_Set_CompProjection          (this check is needed)
                                 .true.,  &  ! CALLED_MR_Initialize_Met_Grids        (this check is needed)
                                 .true.)     ! CALLED_MR_Set_Met_Times               (this check is needed)

        ! Now we have MR_dum3d_metH; interpolate onto computational grid
        !  Since MR_dum3d_metH and MR_dum3d_compH have the same z-coordinate, we only
        !  need to do 2d regridding on each k-slice
      allocate(tmp_regrid2d_sp(nx_comp,ny_comp));tmp_regrid2d_sp(:,:)=0.0_sp

        call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum2d_met(1:nx_submet,1:ny_submet),       &
                                nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))

        do i = 1,nx_comp
          do j = 1,ny_comp
            if(isnan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
          enddo
        enddo
        MR_dum2d_comp(:,:) = tmp_regrid2d_sp(:,:)

      deallocate(tmp_regrid2d_sp)

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      return

      end subroutine MR_Regrid_Met2d_to_Comp2d


!##############################################################################
!
!     MR_DelMetP_Dx
!
!     Calculated the x derivative of the variable in MR_dum3d_MetP.
!     This subroutine expects the calling program to populate MR_dum3d2_metP.
!
!     Sets  : MR_dum3d_metH
!
!##############################################################################

      subroutine MR_DelMetP_Dx

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer :: i
      integer :: lside,rside
      real(kind=sp) :: dx_fac
      real(kind=sp) :: KM_2_M

      KM_2_M = 100.0_sp

      do i=1,nx_submet
        if(i.eq.1)then
          if(IsPeriodic_CompGrid)then
            ! Two-sided derivative wrapping around
            lside  = nx_submet
            rside  = i+1
            dx_fac = 2.0_sp
          else
            ! One-sided derivative on left
            lside  = i+0
            rside  = i+1
            dx_fac = 1.0_sp
          endif
        elseif(i.eq.nx_submet)then
          if(IsPeriodic_CompGrid)then
            ! Two-sided derivative wrapping around
            lside  = i-1
            rside  = 1
            dx_fac = 2.0_sp
          else
            ! One-sided derivative on right
            lside  = i-1
            rside  = i+0
            dx_fac = 1.0_sp
          endif
        else
          ! Two-sided derivative
          lside  = i-1
          rside  = i+1
          dx_fac = 2.0_sp
        endif
        if(IsLatLon_MetGrid)then
          if(IsRegular_MetGrid)then
            MR_dum3d2_metP(i,:,:) = (MR_dum3d_metP(lside,:,:)  + &
                                     MR_dum3d_metP(rside,:,:)) / &
                                    (rdlambda_MetP_sp(i,:,:)   * &
                                    dx_fac*KM_2_M)
          else
            MR_dum3d2_metP(i,:,:) = (MR_dum3d_metP(lside,:,:)  + &
                                     MR_dum3d_metP(rside,:,:)) / &
                                    (0.5_sp*(rdlambda_MetP_sp(i,:,:)+rdlambda_MetP_sp(rside,:,:)) * &
                                    dx_fac*KM_2_M)
          endif
        else
          if(IsRegular_MetGrid)then
            MR_dum3d2_metP(i,:,:) = (MR_dum3d_metP(lside,:,:)  + &
                                     MR_dum3d_metP(rside,:,:)) / &
                                    (dx_fac*MR_dx_submet(i)*KM_2_M)
          else
            write(MR_global_error,*)"Need to fix DelMetP_Dx for non-regular grids."
            stop 1
          endif
        endif
      enddo


      end subroutine MR_DelMetP_Dx


!##############################################################################
!
!     MR_DelMetP_Dy
!
!     Calculated the y derivative of the variable in MR_dum3d_MetP.
!     This subroutine expects the calling program to populate MR_dum3d2_metP.
!
!     Sets  : MR_dum3d_metH
!
!##############################################################################

      subroutine MR_DelMetP_Dy

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer :: i,j
      integer :: lside,rside
      real(kind=sp) :: dy_fac
      real(kind=sp) :: KM_2_M

      KM_2_M = 100.0_sp

      do j=1,ny_submet
        if(j.eq.1)then
          ! One-sided derivative on left
          lside=j+0
          rside=j+1
          dy_fac = 1.0_sp
        elseif(j.eq.ny_submet)then
          ! One-sided derivative on right
          lside=j-1
          rside=j+0
          dy_fac = 1.0_sp
        else
          ! Two-sided derivative
          lside=j-1
          rside=j+1
          dy_fac = 2.0_sp
        endif
        if(IsLatLon_MetGrid)then
          if(IsRegular_MetGrid)then
            do i=1,nx_submet
              MR_dum3d2_metP(i,j,:) = (MR_dum3d_metP(i,lside,:)  + &
                                       MR_dum3d_metP(i,rside,:)) / &
                                      (rdphi_MetP_sp(j,:)   * &
                                       dy_fac*KM_2_M)
            enddo
          else
            do i=1,nx_submet
              MR_dum3d2_metP(i,j,:) = (MR_dum3d_metP(i,lside,:)  + &
                                       MR_dum3d_metP(i,rside,:)) / &
                                      (0.5_sp*(rdphi_MetP_sp(j,:)+rdphi_MetP_sp(rside,:))   * &
                                       dy_fac*KM_2_M)
            enddo
          endif
        else
          if(IsRegular_MetGrid)then
            MR_dum3d2_metP(:,j,:) = (MR_dum3d_metP(:,lside,:)  + &
                                     MR_dum3d_metP(:,rside,:)) / &
                                    (dy_fac*MR_dy_submet(j)*KM_2_M)
          else
            write(MR_global_error,*)"Need to fix DelMetP_Dy for non-regular grids."
            stop 1
          endif
        endif
      enddo

      end subroutine MR_DelMetP_Dy

!##############################################################################
!
!  Returns the temperature (in K) given the height in km
!
!  From http://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
!
!##############################################################################

      function MR_Temp_US_StdAtm(zin)

      implicit none

      real(kind=sp) :: MR_Temp_US_StdAtm
      real(kind=sp) :: zin

      integer,parameter :: MAXATMOSNODES = 16
      real(kind=sp),dimension(MAXATMOSNODES) :: US_StdAtm_znodes
      real(kind=sp),dimension(MAXATMOSNODES) :: US_StdAtm_Tnodes

      real(kind=sp) :: frac
      real(kind=sp) :: Delta_temp
      integer :: k,kk

      US_StdAtm_znodes = (/ 0.0_sp, 11.0_sp, 20.0_sp, 32.0_sp, &
                           47.0_sp, 51.0_sp, 71.0_sp, 84.852_sp, &
                           90.0_sp,  95.0_sp, 100.0_sp, 105.0_sp,   &   ! These are extensions
                          110.0_sp, 115.0_sp, 120.0_sp, 125.0_sp/)      ! using mean CIRA ref atmos.
      US_StdAtm_Tnodes = (/288.15_sp, 216.65_sp, 216.65_sp, 228.65_sp, &
                           270.65_sp, 270.65_sp, 214.65_sp, 186.946_sp, &
                           183.80_sp, 190.30_sp, 203.50_sp, 228.00_sp,  &
                           265.50_sp, 317.10_sp, 380.60_sp, 452.30_sp/)

      if(zin.le.US_StdAtm_znodes(1))then
        MR_Temp_US_StdAtm = US_StdAtm_Tnodes(1)
      elseif(zin.ge.US_StdAtm_znodes(MAXATMOSNODES))then
        MR_Temp_US_StdAtm = US_StdAtm_Tnodes(MAXATMOSNODES)
        write(MR_global_error,*)"Stopping in Temp_US_StdAtm"
        stop 1
      else
        ! interpolate
        ! Start from the top since we assume the requested point is near the top
        kk = 1
        do k = MAXATMOSNODES-1,1,-1
          if(zin.ge.US_StdAtm_znodes(k).and.zin.lt.US_StdAtm_znodes(k+1)) kk = k
        enddo
        frac = (zin-US_StdAtm_znodes(kk)) / &
               (US_StdAtm_znodes(kk+1)-US_StdAtm_znodes(kk))
        Delta_temp = US_StdAtm_Tnodes(kk+1)-US_StdAtm_Tnodes(kk)
        MR_Temp_US_StdAtm = US_StdAtm_Tnodes(kk) + Delta_temp * frac
      endif

      return

      end function MR_Temp_US_StdAtm

!##############################################################################
!
!  Returns the height in km given the pressure in hPa
!
!  From http://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
!       supplemented by doi=10.1.1.455.9475
!
!##############################################################################

      function MR_Z_US_StdAtm(pin)

      implicit none

      real(kind=sp) :: MR_Z_US_StdAtm  ! in km
      real(kind=sp) :: pin          ! in mb

      integer,parameter :: MAXATMOSNODES = 16
      real(kind=sp),dimension(MAXATMOSNODES) :: US_StdAtm_znodes
      real(kind=sp),dimension(MAXATMOSNODES) :: US_StdAtm_pnodes

      real(kind=sp) :: frac
      real(kind=sp) :: Delta_z
      integer :: k,kk

      US_StdAtm_znodes = (/ 0.0_sp,  11.0_sp,  20.0_sp,  32.0_sp,   &
                           47.0_sp,  51.0_sp,  71.0_sp,  84.852_sp, &
                           90.0_sp,  95.0_sp, 100.0_sp, 105.0_sp,   &   ! These are extensions
                          110.0_sp, 115.0_sp, 120.0_sp, 125.0_sp/)      ! using mean CIRA ref atmos.
      US_StdAtm_pnodes = (/1013.25_sp, 226.321_sp, 54.7489_sp, 8.68019_sp, &
                           1.10906_sp, 6.69389e-1_sp, 3.95642e-2_sp, 3.68501e-3_sp, &
                           1.795e-3_sp,7.345e-4_sp,3.090e-4_sp,1.422e-4_sp, &
                           7.362e-5_sp,4.236e-5_sp,2.667e-5_sp,2.666e-5_sp/)

      if(pin.ge.US_StdAtm_pnodes(1))then
        MR_Z_US_StdAtm = US_StdAtm_pnodes(1)
      elseif(pin.le.US_StdAtm_pnodes(MAXATMOSNODES))then
        MR_Z_US_StdAtm = US_StdAtm_pnodes(MAXATMOSNODES)
        write(MR_global_error,*)"MR ERROR: Pressure provided is below lowest value of US StdAtm."
        write(MR_global_error,*)"Stopping in MR_Z_US_StdAtm"
        stop 1
      else
        ! interpolate
        ! Start from the top since we assume the requested point is near the top
        kk = 1
        do k = MAXATMOSNODES-1,1,-1
          if(pin.le.US_StdAtm_pnodes(k).and.pin.gt.US_StdAtm_pnodes(k+1)) kk = k
        enddo
        frac = (pin-US_StdAtm_pnodes(kk)) / &
               (US_StdAtm_pnodes(kk+1)-US_StdAtm_pnodes(kk))
        Delta_z = US_StdAtm_znodes(kk+1)-US_StdAtm_znodes(kk)
        ! Note: this should really be an exponential fit, not linear
        MR_Z_US_StdAtm = US_StdAtm_znodes(kk) + Delta_z * frac
      endif

      return

      end function MR_Z_US_StdAtm


!##############################################################################
!
!  Returns the pressure in hPa given the height in km
!
!  From http://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
!       supplemented by doi=10.1.1.455.9475
!
!##############################################################################

      function MR_Pres_US_StdAtm(zin)

      implicit none

      real(kind=sp) :: MR_Pres_US_StdAtm
      real(kind=sp) :: zin

      integer,parameter :: MAXATMOSNODES = 16
      real(kind=sp),dimension(MAXATMOSNODES) :: US_StdAtm_znodes
      real(kind=sp),dimension(MAXATMOSNODES) :: US_StdAtm_pnodes

      real(kind=sp) :: pres0 = 1013.25_sp
      real(kind=sp) :: skinz   = 7.0_sp

      US_StdAtm_znodes = (/ 0.0_sp,  11.0_sp,  20.0_sp,  32.0_sp,   &
                           47.0_sp,  51.0_sp,  71.0_sp,  84.852_sp, &
                           90.0_sp,  95.0_sp, 100.0_sp, 105.0_sp,   &   ! These are extensions
                          110.0_sp, 115.0_sp, 120.0_sp, 125.0_sp/)      ! using mean CIRA ref atmos.

      US_StdAtm_pnodes = (/1013.25_sp, 226.321_sp, 54.7489_sp, 8.68019_sp, &
                           1.10906_sp, 6.69389e-1_sp, 3.95642e-2_sp, 3.68501e-3_sp, &
                           1.795e-3_sp,7.345e-4_sp,3.090e-4_sp,1.422e-4_sp, &
                           7.362e-5_sp,4.236e-5_sp,2.667e-5_sp,2.666e-5_sp/)

        ! Eq 1.8 of Wallace and Hobbs
        !  This starts to fall apart around 90 km or so
      MR_Pres_US_StdAtm = pres0 * exp(-zin/skinz)

      return

      end function MR_Pres_US_StdAtm

!##############################################################################

!##############################################################################
!
!    MR_QC_3dvar
!
!##############################################################################

      subroutine MR_QC_3dvar(              &
                 ivar,                     &
                 nx_max,ny_max,nz1_max,    &
                 z_array_sp,               &
                 nz2_max,                  &
                 dum_array_sp,             &
                 fill_val_sp,              &
                 bc_low_sp, bc_high_sp)


      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer               ,intent(in)    :: ivar
      integer               ,intent(in)    :: nx_max,ny_max,nz1_max
      real(kind=sp)         ,intent(in)    :: z_array_sp(nx_max,ny_max,nz1_max)
      integer               ,intent(in)    :: nz2_max
      real(kind=sp)         ,intent(inout) :: dum_array_sp(nx_max,ny_max,nz1_max)
      real(kind=sp)         ,intent(in)    :: fill_val_sp
      real(kind=sp),optional,intent(in)    :: bc_low_sp
      real(kind=sp),optional,intent(in)    :: bc_high_sp

      logical,dimension(nz1_max) :: IsFillValue
      logical,dimension(nz1_max) :: InterpolateLev
      logical,dimension(nz1_max) :: ExtrapolateLev

      integer ::  i,j,k,kk,kkk,klow,khigh
      integer :: idx
      real(kind=sp) :: pp,pv,p1,p2,fac


      InterpolateLev(:) = .false.
      ExtrapolateLev(:) = .false.

      ! First check if this variable is on the same pressure grid as GPH
      idx = Met_var_zdim_idx(ivar)  ! index of pressure level for this variable
      if(levs_code(idx).eq.1)then
        ! one-to-one mapping: no special action needed
      else if (levs_code(idx).eq.2)then
        ! truncated grid: this will be addressed at the end of this subroutine
      else if (levs_code(idx).eq.3)then
        ! interpolated:
        ! We need the result to live on the same pressure levels as the GPH array, so
        ! for each level of the primary pressure coordinate, we need to find if a level
        ! of the variable can be copied or must be interpolated
        !
        ! Loop from the top of the variable array (because nz2_max<nz1_max)
        InterpolateLev(:) = .true.
        ExtrapolateLev(:) = .false.
        do k=nz2_max,1,-1  ! these are the p-levels the variable lives on (nz2_max<nz1_max)
          pv = levs_fullmet_sp(idx,k)  ! Pressure level in question

          do kk=nz1_max,1,-1           ! Looping over destination pressure levels
            pp = levs_fullmet_sp(1,kk)
            if(abs(pv-pp).lt.MR_EPS_SMALL)then
              ! pv matched with pp
              InterpolateLev(kk) = .false.
              ExtrapolateLev(kk) = .false.
              ! if the slices are at the same index, do nothing
              if (k.ne.kk)then
                ! Move the slice to the needed position and flag original slice with NaN/Fill_Value
                dum_array_sp(1:nx_max,1:ny_max,kk) = dum_array_sp(1:nx_max,1:ny_max,k)
                dum_array_sp(1:nx_max,1:ny_max,k)  = fill_val_sp
              endif
            elseif(levs_fullmet_sp(idx,1)-levs_fullmet_sp(1,kk).lt.MR_EPS_SMALL)then
                ! Check if this kk level of the master grid is less than the first level of the var grid
                ! pressure level is smaller (higher altitude) than top of GPH
              ExtrapolateLev(kk) = .true.
              InterpolateLev(kk) = .false.
            elseif(levs_fullmet_sp(idx,nz2_max)-levs_fullmet_sp(1,kk).gt.MR_EPS_SMALL)then
                ! This pressure level is greater (lower altitude) than bottom of GPH
              ExtrapolateLev(kk) = .true.
              InterpolateLev(kk) = .false.
            endif
          enddo

        enddo

        ! Now go to each of the layers we need to calculate, assuming we have the end layers
        do k=2,nz1_max-1
          if (InterpolateLev(k))then
            klow = 1
            do kk = k,1,-1
              if (.not.InterpolateLev(kk))then
                klow = kk
                exit
              endif
            enddo
            do kk = k,nz1_max
              if (.not.InterpolateLev(kk))then
                khigh = kk
                exit
              endif
            enddo
            pp = levs_fullmet_sp(1,k)     ! This is the pressure level we need to populate
            p1 = levs_fullmet_sp(1,klow)  ! pressure at low-altitude end of bracket
            p2 = levs_fullmet_sp(1,khigh) ! pressure at high-altitude end of bracket

            fac = (log(pp/p1))/(log(p2/p1))
            dum_array_sp(1:nx_max,1:ny_max,k) = dum_array_sp(1:nx_max,1:ny_max,klow) + &
                                                fac * (dum_array_sp(1:nx_max,1:ny_max,khigh) - &
                                                       dum_array_sp(1:nx_max,1:ny_max,klow))
          endif
        enddo
      else
        ! Should not be here!
        ! We assume that the variable lives on the GPH grid or a subset
        write(MR_global_error,*)'MR ERROR: Variable has more levels than GPH'
        write(MR_global_log  ,*)'MR ERROR: Variable has more levels than GPH'
        stop 1
      endif

      do i=1,nx_max
        do j=1,ny_max
          ! find all fill values in the column
          IsFillValue = .false.
          do k=1,nz2_max
            if(isnan(dum_array_sp(i,j,k)).or.      &  ! Some windfiles use NaN's for fill values
               dum_array_sp(i,j,k).eq.fill_val_sp) &  ! Others have a specific number for fill
                   IsFillValue(k) = .true.
          enddo

          ! Set lower BC if requested and needed
          if(IsFillValue(1).and.present(bc_low_sp)) &
              dum_array_sp(i,j,1) = bc_low_sp

          if(nz2_max.gt.1)then
            ! Set upper BC if requested and needed
            if(IsFillValue(nz2_max).and.present(bc_high_sp)) &
                dum_array_sp(i,j,nz2_max) = bc_high_sp

            ! Now find the lowest non-Fill value
            do klow=1,nz2_max
              if(.not.IsFillValue(klow))exit
            enddo
            ! And the highest non-Fill value
            do khigh=nz2_max,1,-1
              if(.not.IsFillValue(khigh))exit
            enddo

            ! Set bottom values to the lowest real number
            !if(klow.gt.1.and..not.present(bc_low_sp))then
            if(klow.gt.1)then
              do k = 1,klow-1
                dum_array_sp(i,j,k) = dum_array_sp(i,j,klow)
              enddo
            endif

            ! Set top values to the highest real number
            !if(khigh.lt.nz2_max.and..not.present(bc_high_sp))then
            if(khigh.lt.nz2_max)then
              do k = nz2_max,khigh+1,-1
                dum_array_sp(i,j,k) = dum_array_sp(i,j,khigh)
              enddo
            endif

            ! Now check if there are any intermediate FillValues and
            ! interpolate
            do k=klow+1,khigh-1
              if(IsFillValue(k))then
                 ! linearly interpolate in z
                   ! Find first number bove
                 do kk = k+1,khigh,1
                   if(.not.IsFillValue(kk))exit
                 enddo
                   ! Find first number beneath
                 do kkk = max(k-1,1),klow,-1
                   if(.not.IsFillValue(kkk))exit
                 enddo
                 dum_array_sp(i,j,k) = dum_array_sp(i,j,kk) + &
                       (dum_array_sp(i,j,kk)-dum_array_sp(i,j,kkk)) * &
                       (z_array_sp(i,j,kk)-z_array_sp(i,j,kkk))/ &
                       (z_array_sp(i,j,kk)-z_array_sp(i,j,kkk))
              endif

            enddo
          endif
        enddo
      enddo

      if(levs_code(idx).eq.2)then
        ! For truncated grids, assign truncated (high altitude) values to the top bc, if
        ! present, or copy the upper-most valid value if not
        if(present(bc_high_sp))then
          do k=nz2_max+1,nz1_max
            dum_array_sp(1:nx_max,1:ny_max,k) = bc_high_sp
          enddo
        else
          do k=nz2_max+1,nz1_max
            dum_array_sp(1:nx_max,1:ny_max,k) = dum_array_sp(1:nx_max,1:ny_max,nz2_max)
          enddo
        endif
      elseif(levs_code(idx).eq.3)then
        ! Some interpolated grids are also truncated
        do k=1,nz1_max
          if(k.gt.nz1_max/2)then
            if(ExtrapolateLev(k))then
              if(present(bc_high_sp))then
                dum_array_sp(1:nx_max,1:ny_max,k) = bc_high_sp
              else
                dum_array_sp(1:nx_max,1:ny_max,k) = dum_array_sp(1:nx_max,1:ny_max,nz2_max)
              endif
            endif
          else
            if(ExtrapolateLev(k))then
              if(present(bc_low_sp))then
                dum_array_sp(1:nx_max,1:ny_max,k) = bc_low_sp
              else
                dum_array_sp(1:nx_max,1:ny_max,k) = dum_array_sp(1:nx_max,1:ny_max,1)
              endif
            endif
          endif
        enddo
      endif

      end subroutine MR_QC_3dvar

!##############################################################################
!
!    MR_Check_Prerequsites
!
!    Subroutine called at various points in MetReader to verify that the
!    prerequisite steps have been completed.  The five parameters are logicals
!    specifying whether or not a particular test is required.
!
!##############################################################################

      subroutine MR_Check_Prerequsites(test_allocate,    &
                                       test_dimvars,     &
                                       test_compproj,    &
                                       test_initmetgrid, &
                                       test_setmettimes)

      implicit none

      logical ,intent(in) :: test_allocate    ! Check if MR_Allocate_FullMetFileList was called
      logical ,intent(in) :: test_dimvars     ! Check if MR_Read_Met_DimVars
      logical ,intent(in) :: test_compproj    ! Check if MR_Set_CompProjection
      logical ,intent(in) :: test_initmetgrid ! Check if MR_Initialize_Met_Grids
      logical ,intent(in) :: test_setmettimes ! Check if MR_Set_Met_Times

      logical :: failtest_allocate    = .false.
      logical :: failtest_dimvars     = .false.
      logical :: failtest_compproj    = .false.
      logical :: failtest_initmetgrid = .false.
      logical :: failtest_setmettimes = .false.

      ! Check prerequisites
      if(test_allocate.eqv..true.)then
        if(CALLED_MR_Allocate_FullMetFileList.eqv..false.)then
          write(MR_global_error,*)"MR ERROR:  The subroutine MR_Allocate_FullMetFileList has not be called."
          write(MR_global_error,*)"           This must first be called to allocate the space for the"
          write(MR_global_error,*)"           windfiles.  The calling format is:"
          write(MR_global_error,*)"             MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)"
          write(MR_global_error,*)"               where iw      = windfile class (1-5)"
          write(MR_global_error,*)"                     iwf     = windfile format number (1-51)"
          write(MR_global_error,*)"                     igrid   = 2d grid code used (if known)"
          write(MR_global_error,*)"                     idf     = data format ID (ascii, netcdf, grib, etc)"
          write(MR_global_error,*)"                     iwfiles = number of windfiles to be read"
          write(MR_global_error,*)" "
          write(MR_global_error,*)"           The array MR_windfiles(1:iwfiles) must then be filled within"
          write(MR_global_error,*)"           the calling program."
          failtest_allocate = .true.
        endif
      endif

      if(test_dimvars.eqv..true.)then
        if(CALLED_MR_Read_Met_DimVars.eqv..false.)then
          write(MR_global_error,*)"MR ERROR:  The subroutine MR_Read_Met_DimVars has not be called."
          write(MR_global_error,*)"           This must first be called to determine the size (x,y,z,t) of"
          write(MR_global_error,*)"           the windfiles.  The calling format is:"
          write(MR_global_error,*)"             MR_Read_Met_DimVars(optional_argument = iyear)"
          write(MR_global_error,*)"               where the argument iyear determines nt for leap/non-leap years"
          failtest_dimvars = .true.
        endif
      endif

      if(test_compproj.eqv..true.)then
        if(CALLED_MR_Set_CompProjection.eqv..false.)then
          write(MR_global_error,*)"MR ERROR:  The subroutine MR_Set_CompProjection has not be called."
          write(MR_global_error,*)"           This must first be called to determine the subgrid needed"
          write(MR_global_error,*)"           from the windfiles.  The calling format is:"
          write(MR_global_error,*)"             MR_Set_CompProjection(LL_flag,ipf,lam0,phi0,phi1,phi2,ko,Re)"
          write(MR_global_error,*)"               where LL_flag = 0,1 (for projected or lon/lat, respectively)"
          write(MR_global_error,*)"                     ipf     = 1-5 (projection code)"
          write(MR_global_error,*)"                     lam0    = central longitude"
          write(MR_global_error,*)"                     phi0    = first latitude for projection"
          write(MR_global_error,*)"                     phi1    = second latitude for projection"
          write(MR_global_error,*)"                     phi2    = third latitude for projection"
          write(MR_global_error,*)"                     ko      = scaling term (<=1.0)"
          write(MR_global_error,*)"                     Re      = radius of earth in km"
          failtest_compproj = .true.
        endif
      endif

      if(test_initmetgrid.eqv..true.)then
        if(CALLED_MR_Initialize_Met_Grids.eqv..false.)then
          write(MR_global_error,*)"MR ERROR:  The subroutine MR_Initialize_Met_Grids has not be called."
          write(MR_global_error,*)"           This must first be called to allocate the computational grid needed"
          write(MR_global_error,*)"           and to generate the mapping from the windfiles.  The calling format is:"
          write(MR_global_error,*)"             MR_Initialize_Met_Grids(nx,ny,nz,dumx_sp,dumy_sp,dumz_sp,periodic)"
          write(MR_global_error,*)"               where nx,ny,nz = size of comp. grid in x,y,z"
          write(MR_global_error,*)"                     dumx_sp  = cell-centered grid points in x"
          write(MR_global_error,*)"                     dumy_sp  = cell-centered grid points in y"
          write(MR_global_error,*)"                     dumz_sp  = cell-centered grid points in z"
          write(MR_global_error,*)"                     periodic = true or false for periodic grids"
          failtest_initmetgrid = .true.
        endif
      endif

      if(test_setmettimes.eqv..true.)then
        if(CALLED_MR_Set_Met_Times.eqv..false.)then
          write(MR_global_error,*)"MR ERROR:  The subroutine MR_Set_Met_Times has not be called."
          write(MR_global_error,*)"           This must first be called to trim the list of windfiles to just"
          write(MR_global_error,*)"           that needed by the simulation.  The calling format is:"
          write(MR_global_error,*)"             MR_Set_Met_Times(eStartHour,Duration)"
          write(MR_global_error,*)"               where eStartHour = start time needed in hours since basetime"
          write(MR_global_error,*)"                     Duration   = length of time needed"
          failtest_setmettimes = .true.
        endif
      endif

      if(failtest_allocate    .or.&
         failtest_dimvars     .or.&
         failtest_compproj    .or.&
         failtest_initmetgrid .or.&
         failtest_setmettimes)then
        stop 1
      endif

      end subroutine MR_Check_Prerequsites

!##############################################################################

      end module MetReader

