      subroutine MR_Set_Gen_Index_GRIB(grib_file)

      use eccodes

      use, intrinsic :: iso_fortran_env, only : &
         real32,real64,input_unit,output_unit,error_unit

      implicit none
      !implicit none (type, external)

        ! These single and double precision parameters should be 4 and 8
      integer, parameter :: sp = real32  ! selected_real_kind( 6,   37) ! single precision
      integer, parameter :: dp = real64  ! selected_real_kind(15,  307) ! double precision

      character(len=130),intent(in)  :: grib_file

      integer            :: MR_global_info
      integer            :: MR_global_error

      character(len=130)  :: index_file

      integer            :: idx
      integer            :: iret
      integer            :: i,j,k,l,t
      integer            :: count1
        ! Used for keys
      integer(kind=sp) ,dimension(:),allocatable :: level_idx
      integer(kind=sp)  :: levelSize

      integer(kind=sp) ,dimension(:),allocatable :: discipline_idx
      integer(kind=sp) ,dimension(:),allocatable :: parameterCategory_idx
      integer(kind=sp) ,dimension(:),allocatable :: parameterNumber_idx
      integer(kind=sp) ,dimension(:),allocatable :: forecastTime_idx
      integer(kind=sp)  :: disciplineSize
      integer(kind=sp)  :: parameterCategorySize
      integer(kind=sp)  :: parameterNumberSize
      integer(kind=sp)  :: forecastTimeSize

      ! Message identifier.
      integer          :: igrib
      integer          :: gribver
      logical          :: Got_Version

      ! Initialization
      count1           = 0
      MR_global_info   = 6
      MR_global_error  = 0
      Got_Version      = .false.

      ! First, open grib file and determine if we have a grib1 or grib2 file
      call codes_open_file(idx,trim(adjustl(grib_file)),'r')
      call codes_new_from_file(idx,igrib,CODES_PRODUCT_GRIB,iret)
      call codes_get(igrib,'editionNumber',gribver)
      call codes_close_file(idx)

      if(gribver == 1.or.gribver == 2)then
        write(MR_global_info,*)"Grib version detected = ",gribver
        Got_Version = .true.
      else
        write(MR_global_error,*)&
          'MR ERROR : Could not detect grib version.'
        write(MR_global_error,*)'           version must be 1 or 2.'
        stop 0
      endif

      index_file = trim(adjustl(grib_file)) // ".index"
      write(MR_global_info,*)"Generating index file: ",index_file

      if(gribver == 1)then
        ! For grib1 files, we use the parameter numbers
        !!!  grib_ls -p indicatorOfParameter,date,time,level
          ! create an index from a grib file using some keys
        call codes_index_create(idx,trim(adjustl(grib_file)),&
              'indicatorOfParameter,level')
        call codes_grib_multi_support_on()

          ! get the number of distinct values of all the keys in the index
        call codes_index_get_size(idx,'indicatorOfParameter',parameterNumberSize)
        call codes_index_get_size(idx,'level',levelSize)

          ! allocate the arry to contain the list of distinct values
        allocate(parameterNumber_idx(parameterNumberSize))
        allocate(level_idx(levelSize))

          ! get the list of distinct key values from the index
        call codes_index_get(idx,'indicatorOfParameter',parameterNumber_idx)
        call codes_index_get(idx,'level',level_idx)

        count1=0
        do l=1,parameterNumberSize
          call codes_index_select(idx,'indicatorOfParameter',parameterNumber_idx(l))

          do i=1,levelSize
            call codes_index_select(idx,'level',level_idx(i))
            call codes_new_from_index(idx,igrib, iret)
            do while (iret /= GRIB_END_OF_INDEX)
              count1=count1+1
              call codes_release(igrib)
              call codes_new_from_index(idx,igrib, iret)
            enddo
            call codes_release(igrib)

          enddo ! loop on level
        enddo ! loop on ParamN

        call codes_index_write(idx,trim(adjustl(index_file)))
        call codes_index_release(idx)

      elseif(gribver == 2)then
        ! For grib2 files, we use the message discipline, parameter category, parameter number,
        ! and layer type (pressure level, surface, tropopause, etc.) as keys for builing the index

        !----------------------------------------------------------------------
        ! discipline
        ! 0 Meteorological products
        ! 1 Hydrological products
        ! 2 Land surface products
        !10 Oceanographic products

        ! parameterCategory
        ! 0 Temperature
        ! 1 Moisture
        ! 2 Momentum
        ! 3 Mass
        ! 4 Short-wave Radiation
        ! 5 Long-wave Radiation
        ! 6 Cloud
        ! 7 Thermodynamic Stability indices

        ! parameterNumber
        !  This is a sub-category of the discipline and parameterCategory

        ! typeOfFirstFixedSurface
        !  1 Ground or water surface
        !  2 Cloud base level
        !  3 Level of cloud tops
        !100 Isobaric surface  (Pa)
        !103 Specified height level above ground  (m)
        !106 Depth below land surface  (m)
        !200 Unknown code table entry
        ! http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-1.shtml
        !    shortName
        !            discipline
        !              parameterCategory
        !                  parameterNumber
        ! 1  gh      0 3   5 100   HGT   Geopotential_height_isobaric
        ! 2  u       0 2   2 100  UGRD   u-component_of_wind_isobaric
        ! 3  v       0 2   3 100  VGRD   v-component_of_wind_isobaric
        ! 4  w       0 2   8 100  VVEL   Vertical_velocity_pressure_isobaric
        ! 5  t       0 0   0 100   TMP   Temperature_isobaric
        ! 10 hpbl    0 3 196   1  HPBL   Planetary_Boundary_Layer_Height_surface
        ! 11 u       0 2   2 103  UGRD   u-component_of_wind_height_above_ground
        ! 12 v       0 2   3 103  VGRD   v-component_of_wind_height_above_ground
        ! 13 fricv   0 2 197   1 FRICV   Frictional_Velocity_surface
        ! 15 sd      0 1  11   1  SNOD   Snow_depth_surface
        ! 16 soilw   2 0 192 106 SOILW
        ! Volumetric_Soil_Moisture_Content_depth_below_surface_layer
        ! 17 sr      2 0   1   1  SFCR   Surface_roughness_surface
        ! 18 gust    0 2  22   1  GUST   Wind_speed_gust_surface
        ! 20 pres    0 3   0   2  PRES   Pressure_cloud_base
        ! 21 pres    0 3   0   3  PRES   Pressure_cloud_tops
        ! 23 tcc     0 6   1 200  TCDC   Total_cloud_cover_entire_atmosphere
        ! 30 r       0 1   1 100    RH   Relative_humidity_isobaric
        ! 31 q       0 1   0 100  SPFH   Specific_humidity_isobaric
        ! 32 clwmr   0 1  22 100 CLWMR   Cloud_mixing_ratio_isobaric
        ! 33 snmr    0 1  25 100  SNMR   Snow_mixing_ratio_isobaric
        ! 40 crain   0 1 192   1 CRAIN   Categorical_Rain_surface
        ! 41 csnow   0 1 195   1 CSNOW   Categorical_Snow_surface
        ! 42 cfrzr   0 1 193   1 CFRZR   Categorical_Freezing_Rain_surface
        ! 43 cicep   0 1 194   1 CICEP   Categorical_Ice_Pellets_surface
        ! 44 prate   0 1   7   1 PRATE   Precipitation_rate_surface
        ! 45                   1 CPRAT   Convective_Precipitation_Rate_surface
        ! 45 tp      0 1   8   APCP   Total_precipitation_surface_0_Hour_Accumulation
        ! 46                  ACPCP
        ! Convective_precipitation_surface_0_Hour_Accumulation
        ! 47 ncpcp   0 1   9  NCPCP
        ! Large-scale_precipitation_non-convective_surface_0_Hour_Accumulation
        !----------------------------------------------------------------------

          ! create an index from a grib file using some keys
        call codes_index_create(idx,trim(adjustl(grib_file)),&
              'discipline,parameterCategory,parameterNumber,scaledValueOfFirstFixedSurface,forecastTime')
        call codes_grib_multi_support_on()

          ! get the number of distinct values of all the keys in the index
        call codes_index_get_size(idx,'discipline',disciplineSize)
        call codes_index_get_size(idx,'parameterCategory',parameterCategorySize)
        call codes_index_get_size(idx,'parameterNumber',parameterNumberSize)
        call codes_index_get_size(idx,'scaledValueOfFirstFixedSurface',levelSize)
        call codes_index_get_size(idx,'forecastTime',forecastTimeSize)

          ! allocate the arry to contain the list of distinct values
        allocate(discipline_idx(disciplineSize))
        allocate(parameterCategory_idx(parameterCategorySize))
        allocate(parameterNumber_idx(parameterNumberSize))
        allocate(level_idx(levelSize))
        allocate(forecastTime_idx(forecastTimeSize))

          ! get the list of distinct key values from the index
        call codes_index_get(idx,'discipline',discipline_idx)
        call codes_index_get(idx,'parameterCategory',parameterCategory_idx)
        call codes_index_get(idx,'parameterNumber',parameterNumber_idx)
        call codes_index_get(idx,'scaledValueOfFirstFixedSurface',level_idx)
        call codes_index_get(idx,'forecastTime',forecastTime_idx)

        count1=0
        do l=1,disciplineSize
          call codes_index_select(idx,'discipline',discipline_idx(l))
          do j=1,parameterCategorySize
            call codes_index_select(idx,'parameterCategory',parameterCategory_idx(j))
            do k=1,parameterNumberSize
              call codes_index_select(idx,'parameterNumber',parameterNumber_idx(k))
              do i=1,levelSize
                call codes_index_select(idx,&
                  'scaledValueOfFirstFixedSurface',level_idx(i))
                do t=1,forecastTimeSize
                     !  forecastTime_idx(t))
                  call codes_index_select(idx,'forecastTime',&
                    forecastTime_idx(t))
                call codes_new_from_index(idx,igrib, iret)
                do while (iret /= GRIB_END_OF_INDEX)
                   count1=count1+1
                   call codes_release(igrib)
                   call codes_new_from_index(idx,igrib, iret)
                enddo
                call codes_release(igrib)
                enddo ! loop on forecastTime
              enddo ! loop on level
            enddo ! loop on parameterNumber
          enddo ! loop on parameterCategory
        enddo ! loop on discipline

        call codes_index_write(idx,trim(adjustl(index_file)))
        call codes_index_release(idx)

      endif  ! grib version number

      end subroutine MR_Set_Gen_Index_GRIB

