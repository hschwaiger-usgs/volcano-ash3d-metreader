!##############################################################################
!
!     MR_Set_Met_NCEPGeoGrid
!
!     This subroutine reads the variable and dimension IDs, fills the
!     coordinate dimension variables, and calculates xLL_meso, xUR_meso,
!     yLL_meso, and yUR_meso.
!
!     Allocated the dummy arrays for storing met data on met and computational
!     grids:  dum2d_met(nx_submet,ny_submet)
!             dum3d_metP(nx_submet,ny_submet,np_fullmet)
!             dum3d_metH(nx_submet,ny_submet,nz_comp)
!             dum2d_comp(nx_comp,ny_comp)
!             dum3d_compH(nx_comp,ny_comp,nz_comp)
!
!     Note: NCEP grids described here:
!              http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html
!
!##############################################################################

      subroutine MR_Set_Met_NCEPGeoGrid(igrid)

      use MetReader

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer,intent(in) :: igrid

      write(MR_global_production,*)"--------------------------------------------------------------------------------"
      write(MR_global_production,*)"----------                          MR_Set_Met_NCEPGeoGrid            ----------"
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      if(igrid.eq.1227)then
        ! CONUS 3.0-km Lambert Conformal
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID218
        !        LambertConformal_Projection:grid_mapping_name = "lambert_conformal_conic" ;
        !        LambertConformal_Projection:latitude_of_projection_origin = 38.5.;
        !        LambertConformal_Projection:longitude_of_central_meridian = 265.;
        !        LambertConformal_Projection:standard_parallel = 38.5 ;
        !        LambertConformal_Projection:earth_radius = 6371229. ;
        ! proj +proj=lcc +lon_0=262.5 +lat_0=38.5 +lat_1=38.5 +lat_2=38.5 +R=6371.229
        !   226.541 12.190
        !     -4226.108 -832.6978
        !   310.615 57.290
        !      3246.974 4372.859
        !  latlonflag  = 0         : projected grid
        !  projflag    = 4         : Lambert conformal conic
        !  lam0     = 265.0     : longitude of projection point
        !  phi0        =  25.0     : latitude of projection point
        !  phi1        =  25.0     : latitude of cone intersection 1
        !  phi2        =  25.0     : latitude of cone intersection 2
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 4 262.5 38.5 38.5 38.5 6371.229    #Proj flags and params  

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

        Met_iprojflag     = 4
        Met_lam0          =  262.5_8
        Met_phi0          =  38.5_8
        Met_phi1          =  38.5_8
        Met_phi2          =  38.5_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

      elseif(igrid.eq.1221)then
        ! NAM 32-km Lambert Conformal used by NARR (used Met_Re=6367.470, not 6371.229)
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID221
        !        LambertConformal_Projection:grid_mapping_name = "lambert_conformal_conic" ;
        !        LambertConformal_Projection:latitude_of_projection_origin = 50. ;
        !        LambertConformal_Projection:longitude_of_central_meridian = 253. ;
        !        LambertConformal_Projection:standard_parallel = 50. ;
        !    Reported NARR Reanal    Lambert_Conformal:GRIB_param_grid_radius_spherical_earth = 6367.47 ;
        ! proj +proj=lcc +lon_0=-107.0 +lat_0=50.0 +lat_1=50.0 +lat_2=50.0 +R=6367.47
        !   214.50  1.0
        !     -5629.34  -4609.85
        !   357.43 46.352
        !      5661.26  4344.51
        !  latlonflag  = 0         : projected grid
        !  projflag    = 4         : Lambert conformal conic
        !  lam0     = 265.0     : longitude of projection point
        !  phi0        =  25.0     : latitude of projection point
        !  phi1        =  25.0     : latitude of cone intersection 1
        !  phi2        =  25.0     : latitude of cone intersection 2
        !  radius      =  6367.47 : earth radius for spherical earth
        ! 0 4 -107.0 50.0 50.0 50.0 6367.47    #Proj flags and params  

        IsLatLon_MetGrid  = .false.   
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .false.
        MR_Reannalysis    = .true.

        Met_iprojflag     = 4
        Met_lam0          =  -107.0_8
        Met_phi0          =  50.0_8
        Met_phi1          =  50.0_8
        Met_phi2          =  50.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6367.470_8

      elseif(igrid.eq.1050)then
         ! Not an NCEP grid
         !  This grid is for the WRF runs (must be read from file)

        IsLatLon_MetGrid  = .true.   ! this might be reset
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

      elseif(igrid.eq.1041)then
         ! Not an NCEP grid
         !  This grid is for the NASA Np

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

      elseif(igrid.eq.1040.or.igrid.eq.1024)then
         ! Not an NCEP grid
         !  This grid is for the NASA GEOS-5 Cp or MERRA-2

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

      elseif(igrid.eq.1033)then
         ! Not an NCEP grid
         !  This grid is for the CAM files

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

      elseif(igrid.eq.1032)then
         ! Not an NCEP grid
         !  This grid is for the AFWA files

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

      elseif(igrid.eq.1030)then
         ! Not an NCEP grid
         !  This grid is for the ECMWF ERA-20c

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .false.
        isGridRelative    = .true.

      elseif(igrid.eq.1029)then
         ! Not an NCEP grid
         !  This grid is for the ECMWF ERA5

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .false.
        isGridRelative    = .true.

      elseif(igrid.eq.1027)then
         ! Not an NCEP grid
         !  This grid is for the NOAA Reanalysis
        if(MR_iversion.eq.2)then
          !v2
          IsLatLon_MetGrid  = .true.
          IsGlobal_MetGrid  = .true.
          IsRegular_MetGrid = .true.
          isGridRelative    = .true.
        elseif(MR_iversion.eq.3)then
          !v3
          IsLatLon_MetGrid  = .true.
          IsGlobal_MetGrid  = .true.
          IsRegular_MetGrid = .false.
          isGridRelative    = .true.
        endif
      elseif(igrid.eq.2)then
       ! Used by NCEP DOE reanalysis, NCEP-1
       !  http://www.nco.ncep.noaa.gov/pmb/docs/on388/grids/grid002.gif

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

      elseif(igrid.eq.3)then
        ! Used by GFS forecast
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID3
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/grids/grid003.gif

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

      elseif(igrid.eq.4)then
        ! Used by GFS forecast
         !  http://www.nco.ncep.noaa.gov/pmb/docs/on388/grids/grid003.gif

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

      elseif(igrid.eq.45)then
        ! Used by JMA 55
          !  http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID45
          !  http://www.nco.ncep.noaa.gov/pmb/docs/on388/grids/grid045.gif

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

      elseif(igrid.eq.91)then
        ! NAM 3-km Polar Stereographic
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID91
        !   1822145-point (1649x1105) N. Hemisphere Polar Stereographic grid
        !   oriented 150W; 
        !        Polar_Stereographic:grid_mapping_name = "polar_stereographic" ;
        !        Polar_Stereographic:longitude_of_projection_origin = 210. ;
        !        Polar_Stereographic:straight_vertical_longitude_from_pole = 210.;
        !        Polar_Stereographic:scale_factor_at_projection_origin = 0.933 ;
        !        Polar_Stereographic:latitude_of_projection_origin = 90. ;
        !        Polar_Stereographic:earth_shape = "Earth spherical with radius of 6371229.0 m" ;
        !        Polar_Stereographic:GRIB_param_Dx = 2976.0  ;
        !        Polar_Stereographic:GRIB_param_Dy = 2976.0 ;
        !        Polar_Stereographic:GRIB_param_GDSkey = -1633822368 ;
        !        Polar_Stereographic:GRIB_param_La1 = 40.53 ;
        !        Polar_Stereographic:GRIB_param_LaD = 60. ;
        !        Polar_Stereographic:GRIB_param_Lo1 = 181.42899 ;
        !        Polar_Stereographic:GRIB_param_LoV = 210. ;
        !        Polar_Stereographic:GRIB_param_NpProj = "true" ;
        !        Polar_Stereographic:GRIB_param_Nx = 1649 ;
        !        Polar_Stereographic:GRIB_param_Ny = 1105 ;
        !        Polar_Stereographic:GRIB_param_ProjFlag = 0 ;
        !        Polar_Stereographic:GRIB_param_Quasi = "false" ;
        !        Polar_Stereographic:GRIB_param_ResCompFlag = 8 ;
        !        Polar_Stereographic:GRIB_param_VectorComponentFlag = "gridRelative" ;
        !        Polar_Stereographic:GRIB_param_Winds = "Relative" ;
        !        Polar_Stereographic:GRIB_param_grid_name = "Polar_Stereographic" ;
        !        Polar_Stereographic:GRIB_param_grid_radius_spherical_earth = 6371229. ;
        !        Polar_Stereographic:GRIB_param_grid_shape = "Earth spherical with radius of 6371229.0 m" ;
        !        Polar_Stereographic:GRIB_param_grid_shape_code = 6 ;
        !        Polar_Stereographic:GRIB_param_grid_type = 20 ;
        !        Polar_Stereographic:GRIB_param_grid_units = "m" ;
        ! proj +proj=stere  +lon_0=210  +lat_0=90 +k_0=0.933 +R=6371.229
        !   181.42899 40.5301
        !      -2619.36159134661 -4810.03724324973
        !   266.3082 63.9757
        !       2285.91081099714 -1523.98097371848
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 1         : polar stereographic projection
        !  lam0     = -150.0    : longitude of projection point
        !  phi0        =  90.0     : latitude of projection point
        !  k0          =  0.933    : scale factor at projection point
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 1 -150.0 90.0 0.933 6371.229    #Proj flags and params

        IsLatLon_MetGrid  = .false.   
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

        Met_iprojflag     = 1
        Met_lam0          = -150.0_8
        Met_phi0          =  90.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

      elseif(igrid.eq.104)then
        ! NAM 90-km Polar Stereographic
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID104
        !   16170-point (147x110) N. Hemisphere Polar Stereographic grid oriented
        !   105W; pole at (75.5,109.5). (NGM Super C grid)
        !   90.75464 km at 60N
        !        PolarStereographic_Projection:grid_mapping_name = "stereographic" ;
        !        PolarStereographic_Projection:longitude_of_projection_origin = 255. ;
        !        PolarStereographic_Projection:latitude_of_projection_origin = 90. ;
        !        PolarStereographic_Projection:scale_factor_at_projection_origin = 0.933012701892219 ;
        !        PolarStereographic_Projection:earth_radius = 6371229. ;
        ! proj +proj=stere  +lon_0=255  +lat_0=90 +k_0=0.933 +R=6371.229
        !
        ! -6761.21 -9846.821
        ! 
        !  6489.02 45.47379
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 1         : polar stereographic projection
        !  lam0     = -105.0    : longitude of projection point
        !  phi0        =  90.0     : latitude of projection point
        !  k0          =  0.933    : scale factor at projection point
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 1 -105.0 90.0 0.933 6371.229    #Proj flags and params

        IsLatLon_MetGrid  = .false.   
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

        Met_iprojflag     = 1
        Met_lam0          = -105.0_8
        Met_phi0          =  90.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

      elseif(igrid.eq.170)then
        ! Global Gaussian Lat/Lon T170
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID170
        ! This is used by the ERA-Itrm data

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .false.
        isGridRelative    = .true.

      elseif(igrid.eq.182)then
        ! HI N.Pacific 
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID182

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

      elseif(igrid.eq.193)then
       ! Used by GFS forecast (0.25)
       !  http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID193

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

      elseif(igrid.eq.196)then
        ! HI 2.5-km Mercator
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID196
        !   72225-point (321x225) Mercator
        !        Mercator:grid_mapping_name = "mercator" ;
        !        Mercator:standard_parallel = 20. ;
        !        Mercator:longitude_of_projection_origin = 198.475006103516 ;
        !        Mercator:earth_shape = "Earth spherical with radius of 6371229.0 m" ;
        !        Mercator:GRIB_param_BasicAngle = 0 ;
        !        Mercator:GRIB_param_Dx = 2500. ;
        !        Mercator:GRIB_param_Dy = 2500. ;
        !        Mercator:GRIB_param_GDSkey = -1286480248 ;
        !        Mercator:GRIB_param_La1 = 18.073 ;
        !        Mercator:GRIB_param_La2 = 23.088 ;
        !        Mercator:GRIB_param_LaD = 20. ;
        !        Mercator:GRIB_param_Lo1 = 198.475 ;
        !        Mercator:GRIB_param_Lo2 = 206.13101 ;
        !        Mercator:GRIB_param_Nx = 321 ;
        !        Mercator:GRIB_param_Ny = 225 ;
        !        Mercator:GRIB_param_Quasi = "false" ;
        !        Mercator:GRIB_param_ResCompFlag = 56 ;
        !        Mercator:GRIB_param_VectorComponentFlag = "gridRelative" ;
        !        Mercator:GRIB_param_Winds = "Relative" ;
        !        Mercator:GRIB_param_grid_name = "Mercator" ;
        !        Mercator:GRIB_param_grid_radius_spherical_earth = 6371229. ;
        !        Mercator:GRIB_param_grid_shape = "Earth spherical with radius of 6371229.0 m" ;
        !        Mercator:GRIB_param_grid_shape_code = 6 ;
        !        Mercator:GRIB_param_grid_type = 10 ;
        !        Mercator:GRIB_param_grid_units = "m" ;
        ! proj +proj=merc  +lat_ts=20.0 +lon_0=198.475 +R=6371.229
        ! 198.475 18.073
        !   0.00    1920.62
        ! 206.131 23.088
        !   800.00  2480.60
        ! 0 5 198.475 20.0 0.933 6371.229    #Proj flags and params

        IsLatLon_MetGrid  = .false.   
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

        Met_iprojflag     = 5
        Met_lam0          = 198.475_8
        Met_phi0          =  20.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

      elseif(igrid.eq.198)then
        ! NAM 6-km Polar Stereographic
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID198
        !   456225-point (825x553) N. Hemisphere Polar Stereographic grid
        !   oriented 150W; 
        !        Polar_Stereographic:grid_mapping_name = "polar_stereographic" ;
        !        Polar_Stereographic:longitude_of_projection_origin = 210. ;
        !        Polar_Stereographic:straight_vertical_longitude_from_pole = 210.;
        !        Polar_Stereographic:scale_factor_at_projection_origin = 0.933 ;
        !        Polar_Stereographic:latitude_of_projection_origin = 90. ;
        !        Polar_Stereographic:earth_shape = "Earth spherical with radius of 6371229.0 m" ;
        !        Polar_Stereographic:GRIB_param_Dx = 5953.0005 ;
        !        Polar_Stereographic:GRIB_param_Dy = 5953.0005 ;
        !        Polar_Stereographic:GRIB_param_GDSkey = -1633822368 ;
        !        Polar_Stereographic:GRIB_param_La1 = 40.53 ;
        !        Polar_Stereographic:GRIB_param_LaD = 60. ;
        !        Polar_Stereographic:GRIB_param_Lo1 = 181.42899 ;
        !        Polar_Stereographic:GRIB_param_LoV = 210. ;
        !        Polar_Stereographic:GRIB_param_NpProj = "true" ;
        !        Polar_Stereographic:GRIB_param_Nx = 825 ;
        !        Polar_Stereographic:GRIB_param_Ny = 553 ;
        !        Polar_Stereographic:GRIB_param_ProjFlag = 0 ;
        !        Polar_Stereographic:GRIB_param_Quasi = "false" ;
        !        Polar_Stereographic:GRIB_param_ResCompFlag = 8 ;
        !        Polar_Stereographic:GRIB_param_VectorComponentFlag = "gridRelative" ;
        !        Polar_Stereographic:GRIB_param_Winds = "Relative" ;
        !        Polar_Stereographic:GRIB_param_grid_name = "Polar_Stereographic" ;
        !        Polar_Stereographic:GRIB_param_grid_radius_spherical_earth = 6371229. ;
        !        Polar_Stereographic:GRIB_param_grid_shape = "Earth spherical with radius of 6371229.0 m" ;
        !        Polar_Stereographic:GRIB_param_grid_shape_code = 6 ;
        !        Polar_Stereographic:GRIB_param_grid_type = 20 ;
        !        Polar_Stereographic:GRIB_param_grid_units = "m" ;
        ! proj +proj=stere  +lon_0=210  +lat_0=90 +k_0=0.933 +R=6371.229
        !   181.42899 40.5301
        !      -2619.36159134661 -4810.03724324973
        !   266.3082 63.9757
        !       2285.91081099714 -1523.98097371848
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 1         : polar stereographic projection
        !  lam0     = -150.0    : longitude of projection point
        !  phi0        =  90.0     : latitude of projection point
        !  k0          =  0.933    : scale factor at projection point
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 1 -150.0 90.0 0.933 6371.229    #Proj flags and params

        IsLatLon_MetGrid  = .false.   
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

        Met_iprojflag     = 1
        Met_lam0          = -150.0_8
        Met_phi0          =  90.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

      elseif(igrid.eq.212)then
        ! CONUS 40-km Lambert Conformal
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID212
        !        LambertConformal_Projection:grid_mapping_name = "lambert_conformal_conic" ;
        !        LambertConformal_Projection:latitude_of_projection_origin = 25.;
        !        LambertConformal_Projection:longitude_of_central_meridian = 265.;
        !        LambertConformal_Projection:standard_parallel = 25. ;
        !        LambertConformal_Projection:earth_radius = 6371229. ;
        ! proj +proj=lcc +lon_0=265.0 +lat_0=25.0 +lat_1=25.0 +lat_2=25.0 +R=6371.229
        !   226.541 12.190
        !     -4226.108 -832.6978
        !   310.615 57.290
        !      3250.731 4368.582
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 4         : Lambert conformal conic
        !  lam0     = 265.0     : longitude of projection point
        !  phi0        =  25.0     : latitude of projection point
        !  phi1        =  25.0     : latitude of cone intersection 1
        !  phi2        =  25.0     : latitude of cone intersection 2
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 4 265.0 25.0 25.0 25.0 6371.229    #Proj flags and params                               

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

        Met_iprojflag     = 4
        Met_lam0          =  265.0_8
        Met_phi0          =  25.0_8
        Met_phi1          =  25.0_8
        Met_phi2          =  25.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

      elseif(igrid.eq.215)then
        ! CONUS 20-km Lambert Conformal
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID215
        !        LambertConformal_Projection:grid_mapping_name = "lambert_conformal_conic" ;
        !        LambertConformal_Projection:latitude_of_projection_origin = 25.;
        !        LambertConformal_Projection:longitude_of_central_meridian = 265.;
        !        LambertConformal_Projection:standard_parallel = 25. ;
        !        LambertConformal_Projection:earth_radius = 6371229. ;
        ! proj +proj=lcc +lon_0=265.0 +lat_0=25.0 +lat_1=25.0 +lat_2=25.0 +R=6371.229
        !   226.541 12.190
        !     -4226.108 -832.6978
        !   310.615 57.290
        !      3250.916 4368.71
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 4         : Lambert conformal conic
        !  lam0     = 265.0     : longitude of projection point
        !  phi0        =  25.0     : latitude of projection point
        !  phi1        =  25.0     : latitude of cone intersection 1
        !  phi2        =  25.0     : latitude of cone intersection 2
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 4 265.0 25.0 25.0 25.0 6371.229    #Proj flags and params  

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

        Met_iprojflag     = 4
        Met_lam0          =  265.0_8
        Met_phi0          =  25.0_8
        Met_phi1          =  25.0_8
        Met_phi2          =  25.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

      elseif(igrid.eq.216)then
        ! NAM 45-km Polar Stereographic
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID216
        !        PolarStereographic_Projection:grid_mapping_name = "stereographic" ;
        !        PolarStereographic_Projection:longitude_of_projection_origin = 225. ;
        !        PolarStereographic_Projection:latitude_of_projection_origin = 90. ;
        !        PolarStereographic_Projection:scale_factor_at_projection_origin = 0.933012701892219 ;
        !        PolarStereographic_Projection:earth_radius = 6371229. ;
        ! proj +proj=stere  +lon_0=225  +lat_0=90 +k_0=0.933 +R=6371.229
        ! 187.0 30.0
        !   -4225.928 -5408.941
        ! 297.15 70.111
        !    1984.072 -638.9415
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 1         : polar stereographic projection
        !  lam0     = -135.0    : longitude of projection point
        !  phi0        =  90.0     : latitude of projection point
        !  k0          =  0.933    : scale factor at projection point
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 1 -135.0 90.0 0.933 6371.229    #Proj flags and params

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

        Met_iprojflag     = 1
        Met_lam0          = -135.0_8
        Met_phi0          =  90.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

      elseif(igrid.eq.218)then
        ! CONUS 12-km Lambert Conformal
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID218
        !        LambertConformal_Projection:grid_mapping_name = "lambert_conformal_conic" ;
        !        LambertConformal_Projection:latitude_of_projection_origin = 25.;
        !        LambertConformal_Projection:longitude_of_central_meridian = 265.;
        !        LambertConformal_Projection:standard_parallel = 25. ;
        !        LambertConformal_Projection:earth_radius = 6371229. ;
        ! proj +proj=lcc +lon_0=265.0 +lat_0=25.0 +lat_1=25.0 +lat_2=25.0 +R=6371.229
        !   226.541 12.190
        !     -4226.108 -832.6978
        !   310.615 57.290
        !      3246.974 4372.859
        !  latlonflag  = 0         : projected grid
        !  projflag    = 4         : Lambert conformal conic
        !  lam0     = 265.0     : longitude of projection point
        !  phi0        =  25.0     : latitude of projection point
        !  phi1        =  25.0     : latitude of cone intersection 1
        !  phi2        =  25.0     : latitude of cone intersection 2
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 4 265.0 25.0 25.0 25.0 6371.229    #Proj flags and params  

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

        Met_iprojflag     = 4
        Met_lam0          =  265.0_8
        Met_phi0          =  25.0_8
        Met_phi1          =  25.0_8
        Met_phi2          =  25.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

      elseif(igrid.eq.221)then
        ! NAM 32-km Lambert Conformal
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID221
        !        LambertConformal_Projection:grid_mapping_name = "lambert_conformal_conic" ;
        !        LambertConformal_Projection:latitude_of_projection_origin = 50. ;
        !        LambertConformal_Projection:longitude_of_central_meridian = 253. ;
        !        LambertConformal_Projection:standard_parallel = 50. ;
        !    NCEP FC        LambertConformal_Projection:earth_radius = 6371229. ;
        ! Note: the NARR grid should use 1221
        ! proj +proj=lcc +lon_0=-107.0 +lat_0=50.0 +lat_1=50.0 +lat_2=50.0 +R=6371.229
        !   214.50  1.0
        !     -5632.668 -4612.566
        !   357.43 46.352
        !      5664.457 4347.222
        !  latlonflag  = 0         : projected grid
        !  projflag    = 4         : Lambert conformal conic
        !  lam0     = 265.0     : longitude of projection point
        !  phi0        =  25.0     : latitude of projection point
        !  phi1        =  25.0     : latitude of cone intersection 1
        !  phi2        =  25.0     : latitude of cone intersection 2
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 4 -107.0 50.0 50.0 50.0 6371.229    #Proj flags and params  

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

        Met_iprojflag     = 4
        Met_lam0          =  -107.0_8
        Met_phi0          =  50.0_8
        Met_phi1          =  50.0_8
        Met_phi2          =  50.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8
!
      elseif(igrid.eq.227)then
        ! CONUS 5.079-km Lambert Conformal
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID218
        !        LambertConformal_Projection:grid_mapping_name = "lambert_conformal_conic" ;
        !        LambertConformal_Projection:latitude_of_projection_origin = 25.0.;
        !        LambertConformal_Projection:longitude_of_central_meridian = 265.;
        !        LambertConformal_Projection:standard_parallel = 25.0 ;
        !        LambertConformal_Projection:earth_radius = 6371229. ;
        ! proj +proj=lcc +lon_0=265.0 +lat_0=25.0 +lat_1=25.0 +lat_2=25.0 +R=6371.229
        !   226.541 12.190
        !     -4226.11  -832.70
        !   310.615 57.290
        !      3250.81  4368.72
        !  latlonflag  = 0         : projected grid
        !  projflag    = 4         : Lambert conformal conic
        !  lam0     = 265.0     : longitude of projection point
        !  phi0        =  25.0     : latitude of projection point
        !  phi1        =  25.0     : latitude of cone intersection 1
        !  phi2        =  25.0     : latitude of cone intersection 2
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 4 265.0 25.0 25.0 25.0 6371.229    #Proj flags and params  

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

        Met_iprojflag     = 4
        Met_lam0          =  265.0_8
        Met_phi0          =  25.0_8
        Met_phi1          =  25.0_8
        Met_phi2          =  25.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

      elseif(igrid.eq.242)then
        ! NAM 11.25-km Polar Stereographic
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID242
        !        PolarStereographic_Projection:grid_mapping_name = "stereographic" ;
        !        PolarStereographic_Projection:longitude_of_projection_origin = 225. ;
        !        PolarStereographic_Projection:latitude_of_projection_origin = 90. ;
        !        PolarStereographic_Projection:scale_factor_at_projection_origin = 0.933012701892219 ;
        !        PolarStereographic_Projection:earth_radius = 6371229. ;
        ! proj +proj=stere  +lon_0=225  +lat_0=90 +k_0=0.933 +R=6371.229
        ! 187.0 30.0
        !   -4225.928 -5408.941
        ! 297.15 70.111
        !    1984.072 -638.9415
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 1         : polar stereographic projection
        !  lam0     = -135.0    : longitude of projection point
        !  phi0        =  90.0     : latitude of projection point
        !  k0          =  0.933    : scale factor at projection point
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 1 -135.0 90.0 0.933 6371.229    #Proj flags and params

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        isGridRelative    = .true.

        Met_iprojflag     = 1
        Met_lam0          = -135.0_8
        Met_phi0          =  90.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

      else
        write(MR_global_info,*)"MR ERROR: MR_Set_Met_NCEPGeoGrid called with invalid code."
        stop 1
      endif
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Set_Met_NCEPGeoGrid

!##############################################################################
!
!     MR_Set_MetComp_Grids
!
!     This subroutine evaluates the full NWP grid for the subgrid needed.
!     If the NWP projection and the computational grid projection differ,
!     then the mapping from each computational grid point onto the NWP grid
!     is calculated through MR_Set_Comp2Met_Map.
!
!     Sets: 
!           n[t,x,y,p]_met     :: sets the size of the dimensions of the sub-met grid
!           [x,y,p]_met_sp     :: arrays holding dimension values of the sub-met grid
!           MR_dum2d_met(nx_submet,ny_submet)
!           MR_dum3d_metP(nx_submet,ny_submet,np_fullmet)
!           MR_dum3d_metH(nx_submet,ny_submet,nz_comp)
!           MR_dum2d_comp(nx_comp,ny_comp)
!           MR_dum3d_compP(nx_comp,ny_comp,np_fullmet)
!           MR_dum3d_compH(nx_comp,ny_comp,nz_comp)
!           CompPoint_on_subMet_idx
!           bilin_map_wgt
!           CompPoint_X_on_Met_sp
!           CompPoint_Y_on_Met_sp
!
!##############################################################################

      subroutine MR_Set_MetComp_Grids

      use MetReader
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer :: i, j
      integer :: ii,jj
      integer :: isubmet, jsubmet

      real(kind=sp) :: xLL,yLL
      real(kind=sp) :: xUR,yUR

      real(kind=sp) :: dely_sp,x_loc,y_loc,xc_sp,yc_sp,xfrac_sp,yfrac_sp
      integer :: ilat,ilon,ix1,ix2,iy1,iy2

      real(kind=dp) :: ptlon,ptlat,xin,yin,de_x,de_y

      real(kind=sp) :: xc,xfrac,yc,yfrac,px,py
      real(kind=sp) :: x_start_sub,y_start_sub

      real(kind=dp) :: xout,yout
      logical       :: cond1, cond2, cond3
      integer       :: nx_tmp

      INTERFACE
        subroutine MR_Set_Comp2Met_Map
        end subroutine
      END INTERFACE

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------                 MR_Set_MetComp_Grids                       ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      call MR_Set_Comp2Met_Map

      ! Now calculate the indices of the subgrid containing the
      ! computational grid.  Note, the subgrid of the wind file
      ! will be much coarser than the computational grid
      if(Map_Case.eq.1.or. & !  Both Comp Grid and Met grids are Lat/Lon
         Map_Case.eq.2)then  !  Both Comp Grid and Met grids are the same projection
        xLL = x_comp_sp(1)
        yLL = y_comp_sp(1)
        xUR = x_comp_sp(nx_comp)
        yUR = y_comp_sp(ny_comp)
      else
        write(MR_global_info,*)"Met and comp grids differ:"
        write(MR_global_info,2504)
        write(MR_global_info,2505)x_comp_sp(1),&
                     y_comp_sp(1),&
                     CompPoint_X_on_Met_sp(1,1),&
                     CompPoint_Y_on_Met_sp(1,1)
        write(MR_global_info,2505)x_comp_sp(nx_comp),&
                     y_comp_sp(1),&
                     CompPoint_X_on_Met_sp(nx_comp,1),&
                     CompPoint_Y_on_Met_sp(nx_comp,1)
        write(MR_global_info,2505)x_comp_sp(nx_comp),&
                     y_comp_sp(ny_comp),&
                     CompPoint_X_on_Met_sp(nx_comp,ny_comp),&
                     CompPoint_Y_on_Met_sp(nx_comp,ny_comp)
        write(MR_global_info,2505)x_comp_sp(1),&
                     y_comp_sp(ny_comp),&
                     CompPoint_X_on_Met_sp(1,ny_comp),&
                     CompPoint_Y_on_Met_sp(1,ny_comp)
        write(MR_global_info,*)" "

          ! This the branch for when Met and Comp grids differ
        xLL = minval(CompPoint_X_on_Met_sp(:,:))
        yLL = minval(CompPoint_Y_on_Met_sp(:,:))
        xUR = maxval(CompPoint_X_on_Met_sp(:,:))
        yUR = maxval(CompPoint_Y_on_Met_sp(:,:))
 2504   format(8x,'Comp grid corner',28x,'Met grid corner')
 2505   format(4x,'(',f10.4,',',f10.4,')',8x,'--->',8x,'(',f10.4,',',f10.4,')')
      endif

      if(IsLatLon_MetGrid)then
        if(xLL.gt.x_fullmet_sp(nx_fullmet).and.x_fullmet_sp(nx_fullmet).le.180.0_sp)then
          ! If the comp grid starts in the western hemisphere (xLL>180) and if
          ! the global Met grid only extends to 180, then shift the comp grid
          ! into the domain of the met grid
          xLL=xLL-360.0_sp  ! This should only be true western hemisphere (xLL>180)
          xUR=xUR-360.0_sp  ! cases using MERRA
          x_comp_sp = x_comp_sp-360.0_sp
        endif
      endif

      write(MR_global_info,*)"Region of Met grid required by comp grid (in Met coordinates):"

      write(MR_global_info,2501)
      write(MR_global_info,2502)xLL,yUR,xUR,yUR
      write(MR_global_info,2503)
      write(MR_global_info,2503)
      write(MR_global_info,2502)xLL,yLL,xUR,yLL
      write(MR_global_info,2501)

 2501 format(4x,'----------------------------------------------------------------------')
 2502 format(4x,'| (',f10.4,',',f10.4,')',20x,'(',f10.4,',',f10.4,') |')
 2503 format(4x,'|                                                                    |')

      if(IsPeriodic_CompGrid)then
          ! If the domain is periodic, use the whole x-range of the wind file
          ! including a periodic mapping at either end
        nx_submet = nx_fullmet
        istart = 1
        iend = nx_fullmet
        write(MR_global_info,*) "Computational domain is periodic"
      else
        if(x_fullmet_sp(1).le.xLL)then
          ! Make sure the start of the comp grid is not below the domain of the
          ! met files
          istart = 1
          do i = 1,nx_fullmet
            ! For the start index, we assign the lower node of the interval
            ! Note: cond1 is not satisfied when xLL.eq.x_fullmet_sp(1) so we
            !       must initialize istart to 1
            cond1 = x_fullmet_sp(i  ).lt.xLL
            cond2 = x_fullmet_sp(i+1).ge.xLL
            if(cond1.and.cond2) istart = i
          enddo
        else
          write(MR_global_info,*)"MR ERROR: xLL < x_fullmet_sp(1)"
          write(MR_global_info,*)"     x_fullmet_sp(1) = ",x_fullmet_sp(1)
          write(MR_global_info,*)"     xLL             = ",xLL
          stop 1
        endif
        iend = 1
        do i = 1,nx_fullmet
          ! For the end index, we assign the upper node of the interval
          cond1 = x_fullmet_sp(i  ).lt.xUR
          cond2 = x_fullmet_sp(i+1).ge.xUR
          if(cond1.and.cond2) iend = i+1
        enddo
        if(iend.eq.1)then
          if(IsGlobal_MetGrid)then
          ! If iend was not assigned, then the wrap back to the beginning
            iend = nx_fullmet
            do i = 1,nx_fullmet
              ! For the end index, we assign the upper node of the interval
              cond1 = x_fullmet_sp(i  ).lt.xUR-360.0_sp
              cond2 = x_fullmet_sp(i+1).ge.xUR-360.0_sp
              if(cond1.and.cond2) iend = nx_fullmet+i+1
            enddo
          else
            write(MR_global_info,*)"MR ERROR: could not find iend"
            stop 1
          endif
        endif
        nx_submet = iend-istart+1
        write(MR_global_info,*) "Domain is NOT periodic"
      endif

      !SEE if COMPUTATIONAL REGION STRADDLES THE BREAK IN THE WIND FILE
      !  (EITHER THE PRIME OR ANTI-MERIDIAN)
      if(iend.le.nx_fullmet)then        !yes
        wrapgrid = .false.
        write(MR_global_info,*)"Comp grid maps within a contiguous region of the Met grid"
        write(MR_global_info,*)"           wrapgrid = ",wrapgrid
        write(MR_global_info,*)"Met Sub grid specifications:"
        write(MR_global_info,*)"             istart = ",istart
        write(MR_global_info,*)"               iend = ",iend
        write(MR_global_info,*)"          nx_submet = ",nx_submet
        write(MR_global_info,*)"         xsubMetMin = ",x_fullmet_sp(istart)
        write(MR_global_info,*)"                xLL = ",xLL
        write(MR_global_info,*)"                xUR = ",xUR
        write(MR_global_info,*)"         xsubMetMax = ",x_fullmet_sp(iend)
      else                            !no
        if(IsGlobal_MetGrid)then
          wrapgrid = .true.

          ilhalf_fm_l = istart                        ! start index of left half on full met grid
          ilhalf_fm_r = nx_fullmet                    ! end index of left half on full met grid
          ilhalf_nx   = ilhalf_fm_r - ilhalf_fm_l +1  ! width of left half
          irhalf_fm_l = 1                             ! start index of right half on full met grid
          irhalf_fm_r = nx_submet - ilhalf_nx         ! end index of right half on full met grid
          irhalf_nx   = irhalf_fm_r - irhalf_fm_l +1  ! width of right half

          write(MR_global_info,*)"Comp grid span beyond the upper end of the Met grid"
          write(MR_global_info,*)"           wrapgrid = ",wrapgrid
          write(MR_global_info,*)"Met Sub grid specifications:"
          write(MR_global_info,*)"        ilhalf_fm_l = ",ilhalf_fm_l  ! start index of left half on full met grid
          write(MR_global_info,*)"        ilhalf_fm_r = ",ilhalf_fm_r  ! end index of left half on full met grid
          write(MR_global_info,*)"          ilhalf_nx = ",ilhalf_nx    ! width of left half
          write(MR_global_info,*)"        irhalf_fm_l = ",irhalf_fm_l  ! start index of right half on full met grid
          write(MR_global_info,*)"        irhalf_fm_r = ",irhalf_fm_r  ! end index of right half on full met grid
          write(MR_global_info,*)"          irhalf_nx = ",irhalf_nx    ! width of right half

          write(MR_global_info,*)"          nx_submet = ",nx_submet
          write(MR_global_info,*)"ilhalf_nx+irhalf_nx = ",ilhalf_nx+irhalf_nx
          write(MR_global_info,*)"         xsubMetMin = ",x_fullmet_sp(ilhalf_fm_l)
          write(MR_global_info,*)"                xLL = ",xLL
          write(MR_global_info,*)"                xUR = ",xUR
          write(MR_global_info,*)"         xsubMetMax = ",x_fullmet_sp(irhalf_fm_r)
        else
          write(MR_global_info,*)"MR ERROR: Comp grid requirements extend beyond Met grid"
          write(MR_global_info,*)"                xLL = ",xLL
          write(MR_global_info,*)"                xUR = ",xUR
          write(MR_global_info,*)"             istart = ",istart
          write(MR_global_info,*)"               iend = ",iend
          write(MR_global_info,*)"         xsubMetMin = ",x_fullmet_sp(1),x_fullmet_sp(istart)
          write(MR_global_info,*)"         xsubMetMax = ",x_fullmet_sp(iend),x_fullmet_sp(nx_fullmet)
          stop 1
        endif
      endif
      write(MR_global_info,*)"-------------"

      !SEE IF THE MODEL DOMAIN EXTENDS NORTH OR SOUTH OF THE MESOSCALE DOMAIN
      If(UseFullMetGrid)then
          ! This is the special case where the comp grid equals the Met grid
        jstart = 1
        jend = ny_fullmet
      else
        if(y_inverted)then
          ! Find start index
          if(y_fullmet_sp(1).ge.yUR)then
              ! This is the normal case where the UR of the comp grid is within
              ! the lat values of the wind file
            jstart = 1
            do j = 1,ny_fullmet-1
              ! For the start index, we assign the lower node of the interval
              ! Note: cond1 is not satisfied when yUR.eq.y_fullmet_sp(1) so we
              !       must initialize jstart to 1
              cond1 = y_fullmet_sp(j  ).gt.yUR
              cond2 = y_fullmet_sp(j+1).le.yUR
              if(cond1.and.cond2) jstart = j
            enddo
          elseif(IsGlobal_MetGrid)then
              ! There are some special cases where the met grid is global, but do not
              ! have values at the poles (e.g. ERA and NAVGEMHA).  There are occasional
              ! instances where we need values between the extreme lat value and the pole
            jstart = 1
            y_pad_North = .true.
          else
             write(MR_global_info,*)"MR ERROR: yUR > y_fullmet_sp(1)"
             write(MR_global_info,*)"     y_fullmet_sp(1).gt.yUR", &
                        y_fullmet_sp(1),yUR
             stop 1
          endif

          ! Find end index
          if(y_fullmet_sp(ny_fullmet).le.yLL)then
              ! Again, this is the normal case where the LL of the comp grid is within
              ! the lat values of the wind file
            jend = 1
            do j = 1,ny_fullmet-1
              ! For the end index, we assign the lower node of the interval
              cond1 = y_fullmet_sp(j  ).gt.yLL
              cond2 = y_fullmet_sp(j+1).le.yLL
              if(cond1.and.cond2) jend = j + 1
            enddo
          elseif(IsGlobal_MetGrid)then
              ! Here is the same special case as above, but for the southern boundary
            jend = ny_fullmet
            y_pad_South = .true.
          else
             write(MR_global_info,*)"MR ERROR: y_fullmet_sp(ny_fullmet).lt.yLL",&
                               y_fullmet_sp(ny_fullmet),yLL
             stop 1
          endif

        else ! .not.y_inverted
          ! y values go from - to +
          if(y_fullmet_sp(1).le.yLL)then
            jstart = 1
            do j = 1,ny_fullmet-1
              ! For the start index, we assign the lower node of the interval
              ! Note: cond1 is not satisfied when yLL.eq.y_fullmet_sp(1) so we
              !       must initialize jstart to 1
              cond1 = y_fullmet_sp(j  ).lt.yLL
              cond2 = y_fullmet_sp(j+1).ge.yLL
              if(cond1.and.cond2) jstart = j
            enddo
          elseif(IsGlobal_MetGrid)then
              ! There are some special cases where the met grid is global, but does not
              ! have values at the poles (e.g. ERA and NAVGEMHA).  There are occasional
              ! instances where we need values between the extreme lat value and the pole
            jstart = 1
            y_pad_North = .true.
          else
            write(MR_global_info,*)"MR ERROR: yLL < y_fullmet_sp(1)"
            write(MR_global_info,*)"y_fullmet_sp(1) ",y_fullmet_sp(1)
            write(MR_global_info,*)"yLL",yLL
            write(MR_global_info,*)"y_fullmet_sp(1).lt.yLL",y_fullmet_sp(1).lt.yLL
            write(MR_global_info,*)"y_fullmet_sp(1)-yLL",y_fullmet_sp(1)-yLL
            stop 1
          endif
          if(y_fullmet_sp(ny_fullmet).ge.yUR)then
            jend = 1
            do j = 1,ny_fullmet-1
              ! For the end index, we assign the upper node of the interval
              cond1 = y_fullmet_sp(j  ).lt.yUR
              cond2 = y_fullmet_sp(j+1).ge.yUR
              if(cond1.and.cond2) jend = j + 1
            enddo
          elseif(IsGlobal_MetGrid)then
              ! Here is the same special case as above, but for the southern boundary
            jend = ny_fullmet
            y_pad_South = .true.
          else
            write(MR_global_info,*)"MR ERROR: yUR > y_fullmet_sp(ny_fullmet)"
            write(MR_global_info,*)"y_fullmet_sp(my_fullmet)",y_fullmet_sp(ny_fullmet)
            write(MR_global_info,*)"yUr",yUr
            stop 1
          endif
        endif
      endif

      ! Calculate size of arrays that will hold the relevant section of
      ! the mesoscale model
      ny_submet = jend-jstart+1
      write(MR_global_info,*)"-------------"
      write(MR_global_info,*)"             jstart =" ,jstart
      write(MR_global_info,*)"               jend =" ,jend
      write(MR_global_info,*)"          ny_submet =" ,ny_submet
      write(MR_global_info,*)"         ysubMetMin =" ,y_fullmet_sp(jstart)
      write(MR_global_info,*)"                yLL =" ,yLL
      write(MR_global_info,*)"                yUR =" ,yUR
      write(MR_global_info,*)"         ysubMetMax =",y_fullmet_sp(jend)
      write(MR_global_info,*)"-------------"

      if(IsPeriodic_CompGrid)then
        allocate( x_submet_sp(0:nx_submet+1))
        allocate( MR_dx_submet(0:nx_submet+1))
      else
        allocate( x_submet_sp(1:nx_submet))
        allocate( MR_dx_submet(1:nx_submet))
      endif
      allocate( y_submet_sp(1:ny_submet))
      allocate( MR_dy_submet(1:ny_submet))

      ! POPULATE THE X AND Y ARRAYS (Z ARRAY WILL BE POPULATED DIFFERENTLY
      ! AT EACH POINT, IN READ4DWINDARRAY.F90)
      if (wrapgrid) then
        x_submet_sp(          1:ilhalf_nx) = x_fullmet_sp(ilhalf_fm_l:ilhalf_fm_l+ilhalf_nx-1)
        x_submet_sp(ilhalf_nx+1:nx_submet) = x_fullmet_sp(irhalf_fm_l:irhalf_fm_l+irhalf_nx-1) + 360.0_sp
        MR_dx_submet(          1:ilhalf_nx) = MR_dx_met(ilhalf_fm_l:ilhalf_fm_l+ilhalf_nx-1)
        MR_dx_submet(ilhalf_nx+1:nx_submet) = MR_dx_met(irhalf_fm_l:irhalf_fm_l+irhalf_nx-1) + 360.0_sp
        if(IsPeriodic_CompGrid)then
          x_submet_sp(0)           = x_fullmet_sp(nx_submet  ) - 360.0_sp
          x_submet_sp(nx_submet+1) = x_fullmet_sp(nx_submet+1) + 360.0_sp
          MR_dx_submet(0)           = MR_dx_met(nx_submet  )
          MR_dx_submet(nx_submet+1) = MR_dx_met(nx_submet+1)
        endif
      else
        x_submet_sp(1:nx_submet) = x_fullmet_sp(istart:iend)
        MR_dx_submet(1:nx_submet) = MR_dx_met(istart:iend)
        if(IsPeriodic_CompGrid)then
          x_submet_sp(0)           = x_fullmet_sp(nx_submet) - 360.0_sp
          x_submet_sp(nx_submet+1) = x_fullmet_sp(1     ) + 360.0_sp
          MR_dx_submet(0)           = MR_dx_met(nx_submet)
          MR_dx_submet(nx_submet+1) = MR_dx_met(1     )
        endif
      endif
      do j=jstart,jend
        if(y_inverted)then
          y_submet_sp(jend-j+1) = y_fullmet_sp(j)
          MR_dy_submet(jend-j+1) = -MR_dy_met(j) ! Note that we need to negated dy so that
                                                 ! MR_dy_submet is always +
        else
          y_submet_sp(j-jstart+1) = y_fullmet_sp(j)
          MR_dy_submet(j-jstart+1) = MR_dy_met(j)
        endif
      enddo

      ! Set up for interpolation if needed
      write(MR_global_info,*)" Calculating mapping of comp "
      x_start_sub = x_submet_sp(1)
      y_start_sub = y_submet_sp(1)
      do i=1,nx_comp
        do j=1,ny_comp
          px = CompPoint_X_on_Met_sp(i,j)
          py = CompPoint_Y_on_Met_sp(i,j)
          if(IsLatLon_MetGrid)then
            ! When the met grid is lat/lon, we need to make sure that px
            ! map to the correct domain of the met file (-180->180 v.s. 0->360)
            if(IsGlobal_MetGrid)then
              if(IsPeriodic_CompGrid)then
                ! For global Lon/Lat Met data, allow points up to
                ! x_fullmet_sp(nx_fullmet)+dx
                if(px.ge.x_fullmet_sp(1)+360.0_sp)then
                  px=px-360.0_sp
                endif
              else
                ! For global Lon/Lat, but not a periodic comp grid
                if(px.gt.x_submet_sp(nx_submet))then
                  px=px-360.0_sp
                endif
                if(px.lt.x_submet_sp(1))then
                  px=px+360.0_sp
                endif
              endif
            else
              ! For non-global Met data, require values to be strictly within
              ! limits of the SUB-grid (i.e. might be >360)
              if(px.gt.x_submet_sp(nx_submet))then
                px=px-360.0_sp
              endif
              if(px.lt.x_submet_sp(1))then
                px=px+360.0_sp
              endif
            endif
          endif
          if(.not.IsPeriodic_CompGrid)then
            if(px.lt.x_start_sub.or.px.gt.x_submet_sp(nx_submet))then
              write(MR_global_info,*)"MR ERROR: Comp point maps out of sub_Met in x."
              write(MR_global_info,*)"Comp i,j, x      :",i,j,px
              write(MR_global_info,*)"sub_Met xmin,xmax:",x_start_sub,x_submet_sp(nx_submet)
              stop 1
            endif
            if((py.lt.y_start_sub           .and..not.y_pad_South).or.&
               (py.gt.y_submet_sp(ny_submet).and..not.y_pad_North))then
              write(MR_global_info,*)"MR ERROR: Comp point maps out of sub_Met in y."
              write(MR_global_info,*)"Comp i,j, y      :",i,j,px,py
              write(MR_global_info,*)"sub_Met ymin,ymax:",y_start_sub,y_submet_sp(ny_submet)
              stop 1
            endif
          endif

          ! Get the sub_Met index of LL corner of cell containing px,py
          isubmet = 1
          if(IsPeriodic_CompGrid)then
            nx_tmp = nx_submet
          else
            nx_tmp = nx_submet-1
          endif
          do ii = 1,nx_tmp
              ! Set interval inclusive of lower node
            cond1 = px.ge.x_submet_sp(ii  )
            cond2 = px.lt.x_submet_sp(ii+1)
            cond3 = px.le.x_submet_sp(ii+1)
            if(ii.lt.nx_tmp)then
              if(cond1.and.cond2)then
                isubmet = ii
                exit
              endif
            else ! This is when ii = nx_tmp
              if(cond1.and.cond3)then
                isubmet = ii
                exit
              endif
            endif
          enddo
          CompPoint_on_subMet_idx(i,j,1) = isubmet

            ! Check if the point is within the upper and lower bounds
          if(py.lt.y_submet_sp(ny_submet).and.py.ge.y_submet_sp(1))then
            jsubmet = 1
            do jj = 1,ny_submet-1
              ! Set interval inclusive of lower node
              cond1 = py.ge.y_submet_sp(jj  )
              cond2 = py.lt.y_submet_sp(jj+1)
              if(cond1.and.cond2)then
                jsubmet = jj
                exit
              endif
            enddo
          elseif(abs(py-y_submet_sp(ny_submet)).lt.1.0e-7_sp)then
              ! This is to fix the occasional instances where the top comp point
              ! maps almost right on the top submet point
            jsubmet = 1
            do jj = 1,ny_submet-1
              ! Set interval inclusive of lower node
              cond1 = py.ge.y_submet_sp(jj  )
              cond2 = py.le.y_submet_sp(jj+1)
              if(cond1.and.cond2)then
                jsubmet = jj
                exit
              endif
            enddo
          elseif(py.gt.y_submet_sp(ny_submet).and.IsGlobal_MetGrid.and.y_pad_North)then
            jsubmet = ny_submet
          elseif(py.lt.y_submet_sp(1)        .and.IsGlobal_MetGrid.and.y_pad_South)then
            jsubmet = 1
          endif
          CompPoint_on_subMet_idx(i,j,2) = jsubmet

          ! Get fractional position of comp point in met cell
          xfrac=(px-x_submet_sp(isubmet))/MR_dx_submet(isubmet)
          if(py.gt.y_submet_sp(ny_submet).and.IsGlobal_MetGrid.and.y_pad_North)then
              ! If comp point is above all met points
            yfrac=(py- y_submet_sp(ny_submet) ) /  MR_dy_submet(ny_submet)
          elseif(py.lt.y_submet_sp(1).and.IsGlobal_MetGrid.and.y_pad_South)then
              ! If comp point is below all met points
            yfrac=(py- (y_submet_sp(1)-abs(MR_dy_submet(1))) ) / abs(MR_dy_submet(1))
          else
              ! Normal case where comp point is strictly within the met grid
            yfrac=(py-y_submet_sp(jsubmet))/MR_dy_submet(jsubmet)
          endif
          xc = 1.0_sp-xfrac
          yc = 1.0_sp-yfrac

          if(xfrac.gt.1.0_sp.or.xfrac.lt.0.0_sp.or.&
             yfrac.gt.1.0_sp.or.yfrac.lt.0.0_sp)then
            ! The point is mapping outside the expected cell
            write(MR_global_error,*)"MR ERROR : Error calculating Met to Comp mapping."
            write(MR_global_error,*)"Comp point : ",i,j,x_comp_sp(i),y_comp_sp(j)
            write(MR_global_error,*)"Coord on Met: ",CompPoint_X_on_Met_sp(i,j),CompPoint_Y_on_Met_sp(i,j)
            write(MR_global_error,*)"Index on subMet: ",isubmet,jsubmet
            !write(MR_global_error,*)MR_dx_submet(:)
            write(MR_global_error,*)"fractional pos.: ",xfrac,yfrac
            stop 1
          endif

          bilin_map_wgt(i,j,1)=xc*yc
          bilin_map_wgt(i,j,2)=xfrac*yc
          bilin_map_wgt(i,j,3)=xfrac*yfrac
          bilin_map_wgt(i,j,4)=yfrac*xc

        enddo
      enddo
      ! Now we have all the dimension sizes needed for allocating all grids

      ! We might need to rotate the wind vectors on the met grid in place if
      !   we need to convert ER to GR for the same grid or
      !   we need ER vectors from a projected Met grid
      if(.not.isGridRelative.or. &  ! We are dealing with NARR data
               Map_Case.eq.4.or. &  ! Met Grid is projected and Comp grid is Lat/Lon
               Map_Case.eq.5)then   ! Met Grid and Comp grids have different projections
        write(MR_global_info,*)"  Setting up arrays for rotating vectors on Met grid."
        allocate(MR_u_ER_metP(nx_submet,ny_submet,np_fullmet))
        allocate(MR_v_ER_metP(nx_submet,ny_submet,np_fullmet))
        allocate(theta_Met(nx_submet,ny_submet))  ! This holds the angle between the projected
                                                  ! Met grid and the earth grid; used for rotating
                                                  ! grid velocities to Earth-Relative, or in the
                                                  ! special NARR case, rotating ER to GR

          !Tot_Bytes_on_Heap = Tot_Bytes_on_Heap + ip*(nc_nxmax*nc_nymax)
        do i=1,nx_submet
          do j=1,ny_submet
              ! Get lon/lat of point in question
            xin = real(x_submet_sp(i),kind=dp)
            yin = real(y_submet_sp(j),kind=dp)
            call PJ_proj_inv(xin,yin, Met_iprojflag, &
                          Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                           ptlon,ptlat)
              ! Get projected coordinate of de at the current point
            call PJ_proj_for(ptlon+1.0_dp/60.0_dp,ptlat, Met_iprojflag, &
                       Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                       xout,yout)
            de_x = xout-xin
            de_y = yout-yin
              ! Now recover the angle between de and x (of map grid)
            theta_Met(i,j) = atan(de_y/de_x)
          enddo
        enddo
      endif

      ! theta_Met ensures that we have Met data that is Earth-Relative, even if the
      ! underlying wind data are projected.  We might, however, need to rotate these
      ! ER values to a projected computational grid.  So we set up another rotation
      ! to map ER values that were interpolated onto a computational grid to GR
      if(Map_Case.eq.3.or. & ! Met is Lat/Lon, but Comp is projected
         Map_Case.eq.4.or. & ! Met is projected, but Comp is Lat/Lon
         Map_Case.eq.5)then  ! Met Grid and Comp grids have different projections
        write(MR_global_info,*)"  Setting up arrays for rotating vectors on comp grid."
        if(MR_useCompH)allocate(MR_dum3d_compH_2(nx_comp,ny_comp,nz_comp))
        if(MR_useCompP)allocate(MR_dum3d_compP_2(nx_comp,ny_comp,np_fullmet))
        allocate(theta_Comp(nx_comp,ny_comp))
        if(Map_Case.eq.3.or.Map_Case.eq.5)then
          ! We only need to calculate inverse projections if the comp grid is projected.
          ! If met and comp grids differ, first get Met grid winds as Earth-relative
          ! Note: This is only needed if Met grid is projected, ie for Map_Case = 4 or 5
          do i=1,nx_comp
            do j=1,ny_comp
                ! Get lon/lat of point in question
              xin = real(x_comp_sp(i),kind=dp)
              yin = real(y_comp_sp(j),kind=dp)
              call PJ_proj_inv(xin,yin, Comp_iprojflag, &
                            Comp_lam0,Comp_phi0,Comp_phi1,Comp_phi2,Comp_k0,Comp_Re, &
                             ptlon,ptlat)
                ! Get projected coordinate of de at the current point
              call PJ_proj_for(ptlon+1.0_dp/60.0_dp,ptlat, Comp_iprojflag, &
                         Comp_lam0,Comp_phi0,Comp_phi1,Comp_phi2,Comp_k0,Comp_Re, &
                         xout,yout)
              de_x = xout-xin
              de_y = yout-yin
                ! Again, we want the angle between de and x (but now, of comp grid)
              theta_Comp(i,j) = atan(de_y/de_x)
            enddo
          enddo
        endif
      endif

      allocate(MR_dum2d_met_int(nx_submet,ny_submet))
      allocate(MR_dum2d_met(nx_submet,ny_submet))
      allocate(MR_dum3d_metP(nx_submet,ny_submet,np_fullmet))
      allocate(MR_dum3d2_metP(nx_submet,ny_submet,np_fullmet))
      allocate(MR_dum3d_metH(nx_submet,ny_submet,nz_comp))
      allocate(MR_dum2d_comp_int(nx_comp,ny_comp))
      allocate(MR_dum2d_comp(nx_comp,ny_comp))
      allocate(MR_dum3d_compP(nx_comp,ny_comp,np_fullmet))
      allocate(MR_dum3d_compH(nx_comp,ny_comp,nz_comp))
      !  The only 3d met data that persists locally is Geopotential Height.
      !  This is needed for interpolating the other variables onto a Cartesian
      !  grid.
      allocate(MR_geoH_metP_last(nx_submet,ny_submet,np_fullmet))
      allocate(MR_geoH_metP_next(nx_submet,ny_submet,np_fullmet))

      if(MR_iwindformat.eq.25)then
        ! The following is only needed if we are reading 2d variables from NCEP
        ! wind files, but we will set up the grids regardless as a precaution
          ! Here are the weights starting with LL then counter-clockwise
        allocate(amap_iwf25(nx_submet,ny_submet,4))
        amap_iwf25 = 0.0_sp
        !Tot_Bytes_on_Heap = Tot_Bytes_on_Heap + sp*(as_nxmax*as_nymax*4)
          ! Here is the x1,x2,y1,y2 indices
        allocate(imap_iwf25(nx_submet,ny_submet,4))
        imap_iwf25 = 0
        !Tot_Bytes_on_Heap = Tot_Bytes_on_Heap + sp*(as_nxmax*as_nymax*4)

        ! These should be read directly, but for now, just hardcode it.        
        do i = 1,192
          x_in_iwf25_sp(i)=(i-1)*1.875_sp
        enddo
        y_in_iwf25_sp(1:94) = (/ &
         88.542_sp,  86.6531_sp,  84.7532_sp,  82.8508_sp,  80.9473_sp,   79.0435_sp,  77.1394_sp, 75.2351_sp, &
        73.3307_sp,  71.4262_sp,  69.5217_sp,  67.6171_sp,  65.7125_sp,   63.8079_sp,  61.9033_sp, 59.9986_sp, &
        58.0939_sp,  56.1893_sp,  54.2846_sp,  52.3799_sp,  50.4752_sp,   48.5705_sp,  46.6658_sp, 44.7611_sp, &
        42.8564_sp,  40.9517_sp,   39.047_sp,  37.1422_sp,  35.2375_sp,   33.3328_sp,  31.4281_sp, 29.5234_sp, &
        27.6186_sp,  25.7139_sp,  23.8092_sp,  21.9044_sp,  19.9997_sp,    18.095_sp,  16.1902_sp, 14.2855_sp, &
        12.3808_sp, 10.47604_sp,  8.57131_sp,  6.66657_sp,  4.76184_sp,    2.8571_sp, 0.952368_sp, &
      -0.952368_sp,  -2.8571_sp, -4.76184_sp, -6.66657_sp, -8.57131_sp, -10.47604_sp, -12.3808_sp, &
       -14.2855_sp, -16.1902_sp,  -18.095_sp, -19.9997_sp, -21.9044_sp,  -23.8092_sp, -25.7139_sp, &
       -27.6186_sp, -29.5234_sp, -31.4281_sp, -33.3328_sp, -35.2375_sp,  -37.1422_sp,  -39.047_sp, &
       -40.9517_sp, -42.8564_sp, -44.7611_sp, -46.6658_sp, -48.5705_sp,  -50.4752_sp, -52.3799_sp, &
       -54.2846_sp, -56.1893_sp, -58.0939_sp, -59.9986_sp, -61.9033_sp,  -63.8079_sp, -65.7125_sp, &
       -67.6171_sp, -69.5217_sp, -71.4262_sp, -73.3307_sp, -75.2351_sp,  -77.1394_sp, -79.0435_sp, &
       -80.9473_sp, -82.8508_sp, -84.7532_sp, -86.6531_sp,  -88.542_sp /)

        do ilon = 1,nx_submet
          x_loc = max(0.0_sp,x_submet_sp(ilon))
          do i = 1,191
            if(max(0.0_sp,x_in_iwf25_sp(i)).le.x_loc.and.x_in_iwf25_sp(i+1).gt.x_loc)then
              ix1 = i
              ix2 = i+1
              xfrac_sp = (x_loc - x_in_iwf25_sp(i))/1.875_sp
              xc_sp = 1.0_sp - xfrac_sp
              exit ! leave do loop
            endif
          enddo

          do ilat = 1,ny_submet
            y_loc = y_submet_sp(ilat)
            do j = 94,2,-1
              if(y_in_iwf25_sp(j).le.y_loc.and.y_in_iwf25_sp(j-1).gt.y_loc)then
                iy1 = j
                iy2 = j-1
                dely_sp = y_in_iwf25_sp(j-1)-y_in_iwf25_sp(j)
                yfrac_sp = (y_loc - y_in_iwf25_sp(j))/dely_sp
                yc_sp = 1.0_sp - yfrac_sp
                exit ! leave do loop
              endif
            enddo

            imap_iwf25(ilon,ilat,1)=ix1
            imap_iwf25(ilon,ilat,2)=ix2
            imap_iwf25(ilon,ilat,3)=iy1
            imap_iwf25(ilon,ilat,4)=iy2
            amap_iwf25(ilon,ilat,1)=xc_sp*yc_sp
            amap_iwf25(ilon,ilat,2)=xfrac_sp*yc_sp
            amap_iwf25(ilon,ilat,3)=xfrac_sp*yfrac_sp
            amap_iwf25(ilon,ilat,4)=yfrac_sp*xc_sp

          enddo
        enddo

      endif
      !write(MR_global_info,*)"Finished MR_Set_MetComp_Grids"
      !write(MR_global_info,*)" "
      !write(MR_global_production,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Set_MetComp_Grids

!##############################################################################


!##############################################################################
!
!     MR_Set_Comp2Met_Map
!
!     This subroutine generates the mapping between the computational grid and
!     the Met grid.  
!
!     Allocates the mapping arrays:
!
!        CompPoint_on_subMet_idx
!        bilin_map_wgt
!        CompPoint_X_on_Met_sp
!        CompPoint_Y_on_Met_sp
!
!##############################################################################


      subroutine MR_Set_Comp2Met_Map

      use MetReader
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision
      real(kind=dp), parameter :: tol = 1.0e-3_dp

      real(kind=dp) :: dum1,dum2,dum3,dum4,dum5
      integer :: i,j
      real(kind=dp) :: x_in ,y_in
      real(kind=dp) :: x_out,y_out

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------                          MR_Set_Comp2Met_Map               ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! We now have the full definition of the Met grid and the Comp grid
      ! Figure out if we need to do any remapping
      !Possibilities are:
      !(1) Both Comp Grid and Met grids are Lat/Lon
      !(2) Both Comp Grid and Met grids are the same projection
      !(3) Met Grid is Lat/Lon and Comp grid is projected
      !(4) Met Grid is projected and Comp grid is Lat/Lon
      !(5) Met Grid and Comp grids have different projections
      !      For 3,4,5, we need to map the comp points to the Met grid
      !      Set up CompPoint_Metx, CompPoint_Mety (these are the xy (lon/lat) of each comp point)
      !             CompPoint_Meti, CompPoint_Metj (these are the i,j indices of the sub-Met grid)
      !             CompPoint_Met_Wgt (weights given to the four surrounding points
      if(IsLatLon_MetGrid.and.IsLatLon_CompGrid)then
        ! Both Comp and Met are in Lat/Lon 
        Map_Case = 1
      elseif(IsLatLon_MetGrid.and..not.IsLatLon_CompGrid)then
        ! Met is Lat/Lon, but Comp is projected
        Map_Case = 3
      elseif(.not.IsLatLon_MetGrid.and.IsLatLon_CompGrid)then
        ! Met is projected, but Comp is Lat/Lon
        Map_Case = 4
      else
        ! Both Met and Comp are projected.
        ! Test if the projections are the same.
        if(Met_iprojflag.ne.Comp_iprojflag)then
          ! Met and Comp are completely different projection types
          Map_Case = 5
        else
          ! Projections are the same type, test individual parameters
          if(Comp_iprojflag.eq.0)then
            ! Both Comp and Met are non-geographic, Cartesian grids
            Map_Case = 2
          elseif(Comp_iprojflag.eq.1)then
            ! Polar stereographic
            dum1 = abs(Comp_lam0 - Met_lam0)
            dum2 = abs(Comp_phi0 - Met_phi0)
            dum3 = abs(Comp_k0 - Met_k0)
            dum4 = abs(Comp_Re - Met_Re)
            if(dum1.gt.tol.and.&
               dum2.gt.tol.and.&
               dum3.gt.tol.and.&
               dum4.gt.tol)then
              Map_Case = 5
            else
              Map_Case = 2
            endif
          elseif(Comp_iprojflag.eq.2)then
            ! Albers Equal Area
            dum1 = abs(Comp_lam0 - Met_lam0)
            dum2 = abs(Comp_phi0 - Met_phi0)
            dum3 = abs(Comp_phi1 - Met_phi1)
            dum4 = abs(Comp_phi2 - Met_phi2)
            if(dum1.gt.tol.and.&
               dum2.gt.tol.and.&
               dum3.gt.tol.and.&
               dum4.gt.tol)then
              Map_Case = 5
            else
              Map_Case = 2
            endif
          elseif(Comp_iprojflag.eq.3)then
            ! UTM
            stop 1
          elseif(Comp_iprojflag.eq.4)then
            ! Lambert conformal conic (NARR, NAM218, NAM221)
            dum1 = abs(Comp_lam0 - Met_lam0)
            dum2 = abs(Comp_phi0 - Met_phi0)
            dum3 = abs(Comp_phi1 - Met_phi1)
            dum4 = abs(Comp_phi2 - Met_phi2)
            dum5 = abs(Comp_Re - Met_Re)
            if(dum1.gt.tol.and.&
               dum2.gt.tol.and.&
               dum3.gt.tol.and.&
               dum4.gt.tol.and.&
               dum5.gt.tol)then
              Map_Case = 5
            else
              Map_Case = 2
            endif
          elseif(Comp_iprojflag.eq.5)then
            ! Mercator (NAM196)
            dum1 = abs(Comp_lam0 - Met_lam0)
            dum2 = abs(Comp_phi0 - Met_phi0)
            dum3 = abs(Comp_Re - Met_Re)
            if(dum1.gt.tol.and.&
               dum2.gt.tol.and.&
               dum3.gt.tol)then
              Map_Case = 5
            else
              Map_Case = 2
            endif
          endif !Comp_iprojflag
        endif !Met_iprojflag.ne.Comp_iprojflag
      endif ! Met and Comp projected

      if(Map_Case.eq.1)then
        write(MR_global_info,*)"Map_Case = ",Map_Case
        write(MR_global_info,*)"  Both Comp Grid and Met grids are Lat/Lon"
      elseif(Map_Case.eq.2)then
        write(MR_global_info,*)"Map_Case = ",Map_Case
        write(MR_global_info,*)"  Both Comp Grid and Met grids are the same projection"
        write(MR_global_info,*)"   Met_iprojflag",Met_iprojflag
        write(MR_global_info,*)"   Met_lam0",real(Met_lam0,kind=sp)
        write(MR_global_info,*)"   Met_phi0",real(Met_phi0,kind=sp)
        write(MR_global_info,*)"   Met_phi1",real(Met_phi1,kind=sp)
        write(MR_global_info,*)"   Met_phi2",real(Met_phi2,kind=sp)
        write(MR_global_info,*)"   Met_Re  ",real(Met_Re,kind=sp)
        write(MR_global_info,*)"   Met_k0  ",real(Met_k0,kind=sp)
        write(MR_global_info,*)"   Comp_iprojflag",Comp_iprojflag
        write(MR_global_info,*)"   Comp_lam0",real(Comp_lam0,kind=sp)
        write(MR_global_info,*)"   Comp_phi0",real(Comp_phi0,kind=sp)
        write(MR_global_info,*)"   Comp_phi1",real(Comp_phi1,kind=sp)
        write(MR_global_info,*)"   Comp_phi2",real(Comp_phi2,kind=sp)
        write(MR_global_info,*)"   Comp_Re  ",real(Comp_Re,kind=sp)
        write(MR_global_info,*)"   Comp_k0  ",real(Comp_k0,kind=sp)
      elseif(Map_Case.eq.3)then
        write(MR_global_info,*)"Map_Case = ",Map_Case
        write(MR_global_info,*)"  Met Grid is Lat/Lon and Comp grid is projected"
        write(MR_global_info,*)"   Comp_iprojflag",Comp_iprojflag
        write(MR_global_info,*)"   Comp_lam0",real(Comp_lam0,kind=sp)
        write(MR_global_info,*)"   Comp_phi0",real(Comp_phi0,kind=sp)
        write(MR_global_info,*)"   Comp_phi1",real(Comp_phi1,kind=sp)
        write(MR_global_info,*)"   Comp_phi2",real(Comp_phi2,kind=sp)
        write(MR_global_info,*)"   Comp_Re  ",real(Comp_Re,kind=sp)
        write(MR_global_info,*)"   Comp_k0  ",real(Comp_k0,kind=sp)
      elseif(Map_Case.eq.4)then
        write(MR_global_info,*)"Map_Case = ",Map_Case
        write(MR_global_info,*)"  Met Grid is projected and Comp grid is Lat/Lon"
        write(MR_global_info,*)"   Met_iprojflag",Met_iprojflag
        write(MR_global_info,*)"   Met_lam0",real(Met_lam0,kind=sp)
        write(MR_global_info,*)"   Met_phi0",real(Met_phi0,kind=sp)
        write(MR_global_info,*)"   Met_phi1",real(Met_phi1,kind=sp)
        write(MR_global_info,*)"   Met_phi2",real(Met_phi2,kind=sp)
        write(MR_global_info,*)"   Met_Re  ",real(Met_Re,kind=sp)
        write(MR_global_info,*)"   Met_k0  ",real(Met_k0,kind=sp)
      elseif(Map_Case.eq.5)then
        write(MR_global_info,*)"Map_Case = ",Map_Case
        write(MR_global_info,*)"  Met Grid and Comp grids have different projections"
        write(MR_global_info,*)"   Met_iprojflag",Met_iprojflag
        write(MR_global_info,*)"   Met_lam0",real(Met_lam0,kind=sp)
        write(MR_global_info,*)"   Met_phi0",real(Met_phi0,kind=sp)
        write(MR_global_info,*)"   Met_phi1",real(Met_phi1,kind=sp)
        write(MR_global_info,*)"   Met_phi2",real(Met_phi2,kind=sp)
        write(MR_global_info,*)"   Met_Re  ",real(Met_Re,kind=sp)
        write(MR_global_info,*)"   Met_k0  ",real(Met_k0,kind=sp)
        write(MR_global_info,*)"   Comp_iprojflag",Comp_iprojflag
        write(MR_global_info,*)"   Comp_lam0",real(Comp_lam0,kind=sp)
        write(MR_global_info,*)"   Comp_phi0",real(Comp_phi0,kind=sp)
        write(MR_global_info,*)"   Comp_phi1",real(Comp_phi1,kind=sp)
        write(MR_global_info,*)"   Comp_phi2",real(Comp_phi2,kind=sp)
        write(MR_global_info,*)"   Comp_Re  ",real(Comp_Re,kind=sp)
        write(MR_global_info,*)"   Comp_k0  ",real(Comp_k0,kind=sp)
      endif

      allocate(CompPoint_on_subMet_idx(nx_comp,ny_comp,2))
      allocate(bilin_map_wgt(nx_comp,ny_comp,4))
      allocate(CompPoint_X_on_Met_sp(nx_comp,ny_comp))
      allocate(CompPoint_Y_on_Met_sp(nx_comp,ny_comp))

      if(Map_Case.eq.1.or.Map_Case.eq.2)then
        ! Map and Comp are on same grid
        do i=1,nx_comp
          x_in = x_comp_sp(i)
          do j=1,ny_comp
            y_in = y_comp_sp(j)
            CompPoint_X_on_Met_sp(i,j) = real(x_in,kind=sp)
            CompPoint_Y_on_Met_sp(i,j) = real(y_in,kind=sp)
          enddo
        enddo
      elseif(Map_Case.eq.3)then
          ! We just need to map the projected comp grid to the Lon/Lat Met grid
        do i=1,nx_comp
          x_in = x_comp_sp(i)
          do j=1,ny_comp
            y_in = y_comp_sp(j)
            call PJ_proj_inv(x_in, y_in, Comp_iprojflag, &
                           Comp_lam0,Comp_phi0,Comp_phi1,Comp_phi2,Comp_k0,Comp_Re, &
                           x_out,y_out)
            CompPoint_X_on_Met_sp(i,j) = real(x_out,kind=sp)
            CompPoint_Y_on_Met_sp(i,j) = real(y_out,kind=sp)
          enddo
        enddo
      elseif(Map_Case.eq.4)then
          ! We just need to map the Lon/Lat comp grid to the projected Met grid
        do i=1,nx_comp
          x_in = x_comp_sp(i)
          do j=1,ny_comp
            y_in = y_comp_sp(j)
            call PJ_proj_for(x_in, y_in, Met_iprojflag, &
                           Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                           x_out,y_out)
            CompPoint_X_on_Met_sp(i,j) = real(x_out,kind=sp)
            CompPoint_Y_on_Met_sp(i,j) = real(y_out,kind=sp)
          enddo
        enddo
      elseif(Map_Case.eq.5)then
          ! Here, we need to map the projected comp grid to a Lon/Lat grid, then
          ! map to the projected Met grid
        do i=1,nx_comp
          do j=1,ny_comp
            x_in = x_comp_sp(i)
            y_in = y_comp_sp(j)
            call PJ_proj_inv(x_in, y_in, Comp_iprojflag, &
                           Comp_lam0,Comp_phi0,Comp_phi1,Comp_phi2,Comp_k0,Comp_Re, &
                           x_out,y_out)
            x_in = x_out
            y_in = y_out
            call PJ_proj_for(x_in, y_in, Met_iprojflag, &
                           Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                           x_out,y_out)
            CompPoint_X_on_Met_sp(i,j) = real(x_out,kind=sp)
            CompPoint_Y_on_Met_sp(i,j) = real(y_out,kind=sp)
          enddo
        enddo
      endif ! Map_Case

      !write(MR_global_production,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Set_Comp2Met_Map

!##############################################################################
!
!     MR_Regrid_Met2Comp
!
!     This subroutine does the 2-D regridding using the mapping
!     arrays filled in MR_Set_Comp2Met_Map
!
!##############################################################################


      subroutine MR_Regrid_Met2Comp(nx1,ny1,wrk_met,nx2,ny2,wrk_comp)

      use MetReader

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer                         ,intent(in)  :: nx1,ny1
      real(kind=sp),dimension(nx1,ny1),intent(in)  :: wrk_met
      integer                         ,intent(in)  :: nx2,ny2
      real(kind=sp),dimension(nx2,ny2),intent(out) :: wrk_comp

      integer :: i,j,ii,jj
      integer :: nx_max
      real(kind=sp) :: a1,a2,a3,a4
      real(kind=sp) :: tmp

      real(kind=sp),dimension(:,:),allocatable :: wrk_loc

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Regrid_Met2Comp                                    ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      if(IsPeriodic_CompGrid)then
        nx_max = nx1+1
        allocate(wrk_loc(0:nx1+1,0:ny1+1))
        if(y_pad_North)then
          tmp = sum(wrk_met(1:nx1,ny1))/real(nx1,kind=sp)
          wrk_loc(0:nx1+1,ny1+1) = tmp
        else
          wrk_loc(0:nx1+1,ny1+1) = 0.0_sp
        endif
        if(y_pad_South)then
          tmp = sum(wrk_met(1:nx1,1))/real(nx1,kind=sp)
          wrk_loc(0:nx1+1,0) = tmp
        else
          wrk_loc(0:nx1+1,0) = 0.0_sp
        endif
        wrk_loc(1:nx1,1:ny1) = wrk_met(1:nx1,1:ny1) 
        wrk_loc(0    ,1:ny1) = wrk_met(  nx1,1:ny1)
        wrk_loc(nx1+1,1:ny1) = wrk_met(1    ,1:ny1)
      else
        nx_max = nx1
        allocate(wrk_loc(nx1,0:ny1+1))
        if(y_pad_North)then
          tmp = sum(wrk_met(1:nx1,ny1))/real(nx1,kind=sp)
          wrk_loc(1:nx1,ny1+1) = tmp
        else
          wrk_loc(1:nx1,ny1+1) = 0.0_sp
        endif
        if(y_pad_South)then
          tmp = sum(wrk_met(1:nx1,1))/real(nx1,kind=sp)
          wrk_loc(1:nx1,0) = tmp
        else
          wrk_loc(1:nx1,0) = 0.0_sp
        endif
        wrk_loc(1:nx1,1:ny1) = wrk_met(1:nx1,1:ny1)
      endif

      ! Loop over all comp points
      do i = 1,nx2
        do j = 1,ny2
          ! Get the Met cell id this comp point maps to
          ii = CompPoint_on_subMet_idx(i,j,1)
          jj = CompPoint_on_subMet_idx(i,j,2)
          if(ii.lt.1.or.ii.gt.nx_max-1)then
            write(MR_global_info,*)"MR ERROR: ii maps out of grid: ",ii
            stop 1
          endif
          if(jj.lt.0.or.jj.gt.ny1)then
            write(MR_global_info,*)"MR ERROR: jj maps out of grid: ",jj,ny1
            stop 1
          endif
          ! Look up this comp points weights
          a1 = bilin_map_wgt(i,j,1)
          a2 = bilin_map_wgt(i,j,2)
          a3 = bilin_map_wgt(i,j,3)
          a4 = bilin_map_wgt(i,j,4)
          ! Now interpolate from Met to Comp
          wrk_comp(i,j) = a1*wrk_loc(ii  ,jj  ) + &
                          a2*wrk_loc(ii+1,jj  ) + &
                          a3*wrk_loc(ii+1,jj+1) + &
                          a4*wrk_loc(ii  ,jj+1)
        enddo
      enddo

      end subroutine MR_Regrid_Met2Comp

!##############################################################################
!
!     MR_Regrid_P2H_linear
!
!     This subroutine does the 1-D regridding (replaces rgrd1d) linearly
!     interpolating in z.
!
!##############################################################################


      subroutine MR_Regrid_P2H_linear(nzm,z_met ,var_met, &
                                      nzc,z_comp,var_comp)

      use MetReader

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer                     ,intent(in)  :: nzm
      real(kind=sp),dimension(nzm),intent(in)  :: z_met
      real(kind=sp),dimension(nzm),intent(in)  :: var_met
      integer                     ,intent(in)  :: nzc
      real(kind=sp),dimension(nzc),intent(in)  :: z_comp
      real(kind=sp),dimension(nzc),intent(out) :: var_comp

      integer :: km,kc,km_interv
      real(kind=sp) :: a1,a2,dz,z1

      logical :: found_interv

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------      MR_Regrid_P2H_linear                                  ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      var_comp = -9999.0_sp
      ! Loop over all comp points
      km_interv = 1
      do kc = 1,nzc
        found_interv = .false.
        z1 = z_comp(kc)
        ! For each comp point, check which met interval it is in, starting from
        ! the last interval found
        do km = km_interv,nzm-1
        !do km = 1,nzm-1
          if(z1.ge.z_met(km).and.z1.le.z_met(km+1))then
            found_interv = .true.
            km_interv = km
            dz = z_met(km_interv+1)-z_met(km_interv)
            a1 = (z1-z_met(km_interv))/dz
            a2 = 1.0_sp-a1
            var_comp(kc) = var_met(km_interv  ) * a2 + &
                           var_met(km_interv+1) * a1
            exit
          else
            cycle
          endif
        enddo
        ! Check that interval was found
        if(.not.found_interv)then
          write(MR_global_info,*)"MR ERROR:  Did not find interval in vertical 1-D interpolation."
          write(MR_global_info,*)"z_met = "
          write(MR_global_info,*)z_met
          write(MR_global_info,*)" "
          write(MR_global_info,*)"z_comp = "
          write(MR_global_info,*)z_comp
          stop 1
        endif
      enddo

      return

      end subroutine MR_Regrid_P2H_linear

!##############################################################################
!
!     MR_Read_Met_Template
!
!     This subroutine reads an auxiliary file that specifies the custom windfile
!     structure
!
!##############################################################################


      subroutine MR_Read_Met_Template

      use MetReader
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      logical      :: IsThere
      character(len=130) :: lllinebuffer
      character(len=80)  :: Met_projection_line

      integer :: ioerr
      character(len=3)   :: useLeap_str
      integer :: ndims_custom, nvars_custom
      integer :: i
      character     :: dv_char
      integer       :: dimID,varID,zindx,vndim
      real(kind=sp) :: fac
      character(len=30) :: dname
      character(len=71) :: vname
      character(len=5)  :: vname_WMO
      real(kind=8)      :: StepInterval

      write(MR_global_info,*)"Inside MR_Read_Met_Template"
      inquire( file=adjustl(trim(MR_iwf_template)), exist=IsThere )
      if(.not.IsThere)then
        write(MR_global_info,*)"MR ERROR: Could not find NWP template file ",&
                   adjustl(trim(MR_iwf_template))
        write(MR_global_info,*)"          Make sure the calling program sets MR_iwf_template"
        write(MR_global_info,*)"          and that it is linked to the cwd."
        stop 1
      endif

      open(unit=27,file=adjustl(trim(MR_iwf_template)),status='unknown')
      read(27,'(a130)')lllinebuffer
      Met_projection_line = lllinebuffer(1:80)
      call PJ_Set_Proj_Params(Met_projection_line)

      if(PJ_ilatlonflag.eq.0)then
        IsLatLon_MetGrid = .false.
        IsGlobal_MetGrid = .false.
      else
        IsLatLon_MetGrid = .true.
        if(PJ_iprojflag.eq.0)then
          IsGlobal_MetGrid = .false.
        else
          IsGlobal_MetGrid = .true.
        endif
      endif
      Met_iprojflag = PJ_iprojflag
      Met_lam0      = PJ_lam0
      Met_lam1      = PJ_lam1
      Met_lam2      = PJ_lam2
      Met_phi0      = PJ_phi0
      Met_phi1      = PJ_phi1
      Met_phi2      = PJ_phi2
      Met_k0        = PJ_k0
      Met_Re        = PJ_radius_earth

      read(27,'(a130)')lllinebuffer
      read(lllinebuffer,*,iostat=ioerr)StepInterval,useLeap_str

      if (ioerr.eq.0)then
        ! Two values read, process useLeap_str to determine T or F
        if(useLeap_str(1:1).eq.'F'.or.useLeap_str(1:1).eq.'f')then
          MR_useLeap = .false.
          write(MR_global_info,*)"This windfile template specifies that leap years are NOT to"
          write(MR_global_info,*)"be used.  Resetting MR_useLeap = .false."
        endif
      else
        ! default will be whatever is set in the host program or in
        ! MetReader.f90 if the calling program doesn't specify
        read(lllinebuffer,*,iostat=ioerr)StepInterval
      endif
      read(27,'(a130)')lllinebuffer
      read(lllinebuffer,*,err=2002)ndims_custom,nvars_custom

      write(MR_global_info,*)"  Reading dimensions: ",ndims_custom
      do i = 1,ndims_custom
        read(27,'(a130)')lllinebuffer
        read(lllinebuffer,1501)dv_char,dimID,fac,dname
        if(dv_char.ne.'d')then
          write(MR_global_info,*)"MR ERROR : Trying to read variable into dimension"
          write(MR_global_info,*)"dv_char = ",dv_char
          write(MR_global_info,*)"dimID   = ",dimID
          write(MR_global_info,*)"fac     = ",fac
          write(MR_global_info,*)"dname   = ",dname
          stop 1
        endif
        if(dimID.le.9)then
          Met_dim_IsAvailable(dimID) = .true.
          Met_dim_names(dimID)       = adjustl(trim(dname))
          Met_dim_fac(i)             = fac
          write(MR_global_info,*)dimID,' ',Met_dim_names(dimID)
        else
          write(MR_global_info,*)"MR ERROR: dimID too large",dimID
          stop 1
        endif
      enddo
      write(MR_global_info,*)"  Reading variables: ",nvars_custom
      do i = 1,nvars_custom
        read(27,'(a130)')lllinebuffer
        read(lllinebuffer,1511)dv_char,vndim,zindx,varID, &
                                fac,vname_WMO,vname
        if(dv_char.ne.'v')then
          write(MR_global_info,*)"MR ERROR : Trying to read dimension into variable"
          write(MR_global_info,*)"dv_char   = ",dv_char
          write(MR_global_info,*)"vndim     = ",vndim
          write(MR_global_info,*)"zindx     = ",zindx
          write(MR_global_info,*)"varID     = ",varID
          write(MR_global_info,*)"fac       = ",fac
          write(MR_global_info,*)"vname_WMO = ",vname_WMO
          write(MR_global_info,*)"vname     = ",vname
          stop 1
        endif

        if(varID.le.MR_MAXVARS)then
          Met_var_IsAvailable(varID)       = .true.
          Met_var_NC_names(varID)          = adjustl(trim(vname))
          Met_var_WMO_names(varID)         = adjustl(trim(vname_WMO))
          Met_var_ndim(varID)              = vndim
          Met_var_zdim_idx(varID)          = zindx
          Met_var_conversion_factor(varID) = fac
          write(MR_global_info,*)varID,Met_var_WMO_names(varID),' ',Met_var_NC_names(varID)
        else
          write(MR_global_info,*)"MR ERROR: varID too large",varID
          stop 1
        endif
      enddo

      close(27)

      return

1501  format(a1,i9     ,f9.2,a30)
1511  format(a1,i3,i3,i3,f9.2,a7,a71)

!2001  write(MR_global_info,*)  'error reading ForecastInterval'
!      write(MR_global_info,*)lllinebuffer
!      stop 1
2002  write(MR_global_info,*)  'error reading number of custom dims and vars.'
      write(MR_global_info,*)lllinebuffer
      stop 1

      end subroutine MR_Read_Met_Template

