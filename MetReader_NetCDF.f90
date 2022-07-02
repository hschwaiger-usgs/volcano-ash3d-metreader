!##############################################################################
!##############################################################################
!##############################################################################


!##############################################################################
!
!     MR_Read_Met_DimVars_netcdf
!
!     Called once from MR_Read_Met_DimVars 
!
!     This subroutine reads the variable and dimension IDs, and fills the
!     coordinate dimension variables
!
!     After this subroutine completes, the following variables will be set:
!       All the projection parameters of NWP grid
!       Met_dim_names, Met_var_GRIB_names, Met_var_conversion_factor,
!       Met_var_IsAvailable
!       The lengths of all the dimensions of the file
!       p_fullmet_sp (converted to Pa)
!       x_fullmet_sp, y_fullmet_sp
!       IsLatLon_MetGrid, IsGlobal_MetGrid, IsRegular_MetGrid 
!
!##############################################################################

      subroutine MR_Read_Met_DimVars_netcdf

      use MetReader
      use netcdf
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer :: i, k
      real(kind=sp) :: xLL_fullmet
      real(kind=sp) :: yLL_fullmet
      real(kind=sp) :: xUR_fullmet
      real(kind=sp) :: yUR_fullmet

      integer :: ncid
      integer :: ivar,in_var_id,var_ndims
      integer :: i_dim,iivar
      integer :: xdim_id, ydim_id, tdim_id
      integer :: dimlen,maxdimlen
      integer :: nSTAT
      character(len=NF90_MAX_NAME)  :: invar,dimname
      integer :: var_xtype,var_id,idx
      integer :: xtype, length, attnum
      real(kind=dp), parameter :: tol = 1.0e-3_dp
      real(kind=dp) :: x_start,y_start

      integer,dimension(:),allocatable :: var_dimIDs
      logical :: FoundOldDim
      logical :: IsPressureDimension
      real(kind=dp),dimension(:), allocatable :: dum1d_dp
      real(kind=sp),dimension(:), allocatable :: dum1d_sp
      character(len=31)  :: ustring
      logical :: IsTruncatedDim
      character(len=130)   :: infile

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------                MR_Read_Met_DimVars_netcdf                  ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      if(MR_iwind.eq.5)then
        ! For the case where variables are in different files, we will just hard-code the
        ! grid geometry

        if (MR_iwindformat.eq.25) then
           ! NCEP/NCAR reanalysis 2.5 degree files 
           ! https://rda.ucar.edu/datasets/ds090.0
          maxdimlen            = 17
          nlev_coords_detected = 3
          allocate(nlevs_fullmet(nlev_coords_detected));nlevs_fullmet(:)=0
          nlevs_fullmet(1) = 17
          nlevs_fullmet(2) = 12
          nlevs_fullmet(3) = 8
          allocate(levs_code(nlev_coords_detected));levs_code(:)=0
          levs_code(1) = 1
          levs_code(2) = 2
          levs_code(3) = 2
          allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen));levs_fullmet_sp(:,:)=0
          np_fullmet = 17
          allocate(p_fullmet_sp(np_fullmet))
          levs_fullmet_sp(:,:) = 0.0_sp
          p_fullmet_sp(1:np_fullmet) = &
            (/1000.0_sp, 925.0_sp, 850.0_sp, 700.0_sp, 600.0_sp, &
               500.0_sp, 400.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, &
               150.0_sp, 100.0_sp,  70.0_sp,  50.0_sp, 30.0_sp, &
                20.0_sp,  10.0_sp /)
          z_inverted = .false.
          levs_fullmet_sp(1,1:17) = p_fullmet_sp(1:17)*Pressure_Conv_Fac
          levs_fullmet_sp(2,1:12) = p_fullmet_sp(1:12)*Pressure_Conv_Fac
          levs_fullmet_sp(3,1:8 ) = p_fullmet_sp(1:8 )*Pressure_Conv_Fac
          Met_var_zdim_idx( 1) = 1
          Met_var_zdim_idx( 2) = 1
          Met_var_zdim_idx( 3) = 1
          Met_var_zdim_idx( 4) = 2
          Met_var_zdim_idx( 5) = 1
          Met_var_zdim_idx( 7) = 2
          Met_var_zdim_idx(30) = 3

          IsLatLon_MetGrid  = .true.
          IsGlobal_MetGrid  = .true.
          IsRegular_MetGrid = .true.
          nx_fullmet = 144
          ny_fullmet = 73
          allocate(x_fullmet_sp(0:nx_fullmet+1))
          allocate(y_fullmet_sp(ny_fullmet))
          allocate(MR_dx_met(nx_fullmet))
          allocate(MR_dy_met(ny_fullmet))
          dx_met_const = 2.5_sp
          dy_met_const = 2.5_sp
          x_start =   0.0_dp
          y_start =  90.0_dp
          do i = 0,nx_fullmet+1
            x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
          enddo
          x_inverted = .false.
          do i = 1,ny_fullmet
            y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
          enddo
          y_inverted = .true.
          do i = 1,nx_fullmet
            MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
          enddo
          do i = 1,ny_fullmet-1
            MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
          enddo
          MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

          iwf25_scale_facs = 0.0_sp
          iwf25_offsets    = 0.0_sp
          iwf25_scale_facs(1)  = 1.0_sp    ; iwf25_offsets(1)  = 32066.0_sp   ! hgt
          iwf25_scale_facs(2)  = 0.01_sp   ; iwf25_offsets(2)  = 202.66_sp    ! uwnd
          iwf25_scale_facs(3)  = 0.01_sp   ; iwf25_offsets(3)  = 202.66_sp    ! vwnd
          iwf25_scale_facs(4)  = 0.001_sp  ; iwf25_offsets(4)  = 29.765_sp    ! omega
          iwf25_scale_facs(5)  = 0.01_sp   ; iwf25_offsets(5)  = 477.66_sp    ! air (temperature)
          iwf25_scale_facs(6)  = 1.0_sp    ; iwf25_offsets(6)  = 0.0_sp       ! level
          iwf25_scale_facs(20) = 10.0_sp   ; iwf25_offsets(20) = 327650.0_sp  ! pres (lcb)
          iwf25_scale_facs(21) = 10.0_sp   ; iwf25_offsets(21) = 327650.0_sp  ! pres (lct)
          iwf25_scale_facs(30) = 0.01_sp   ; iwf25_offsets(30) = 302.66_sp    ! rhum
          iwf25_scale_facs(31) = 1.0e-6_sp ; iwf25_offsets(31) = 0.032666_sp  ! shum
          iwf25_scale_facs(32) = 1.0e-6_sp ; iwf25_offsets(32) = 0.032666_sp  ! shum
          iwf25_scale_facs(44) = 1.0e-7_sp ; iwf25_offsets(44) = 0.0032765_sp ! prate
          iwf25_scale_facs(45) = 1.0e-7_sp ; iwf25_offsets(45) = 0.0031765_sp ! cprat

        elseif(MR_iwindformat.eq.26)then
         ! JRA-55 reanalysis 1.25 degree files 
         !  https://rda.ucar.edu/datasets/ds628.0/

          maxdimlen            = 37
          nlev_coords_detected = 2
          allocate(nlevs_fullmet(nlev_coords_detected))
          nlevs_fullmet(1) = 37
          nlevs_fullmet(2) = 27
          allocate(levs_code(nlev_coords_detected))
          levs_code(1) = 1
          levs_code(2) = 2
          allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen))
          np_fullmet = 37
          allocate(p_fullmet_sp(np_fullmet))
          levs_fullmet_sp(:,:) = 0.0_sp
          p_fullmet_sp(1:np_fullmet) = &
            (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
               875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
               750.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, 550.0_sp, &
               500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, 300.0_sp, &
               250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
               125.0_sp, 100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp, &
                20.0_sp,  10.0_sp,   7.0_sp,   5.0_sp,   3.0_sp, &
                 2.0_sp,   1.0_sp /)
          z_inverted = .true.

          levs_fullmet_sp(1,1:37) = p_fullmet_sp(1:37)*Pressure_Conv_Fac
          levs_fullmet_sp(2,1:27) = p_fullmet_sp(1:27)*Pressure_Conv_Fac
          Met_var_zdim_idx( 1) = 1
          Met_var_zdim_idx( 2) = 1
          Met_var_zdim_idx( 3) = 1
          Met_var_zdim_idx( 4) = 1
          Met_var_zdim_idx( 5) = 1
          Met_var_zdim_idx( 7) = 1
          Met_var_zdim_idx(30) = 2

          IsLatLon_MetGrid  = .true.
          IsGlobal_MetGrid  = .true.
          IsRegular_MetGrid = .true.
          nx_fullmet = 288
          ny_fullmet = 145
          allocate(x_fullmet_sp(0:nx_fullmet+1))
          allocate(y_fullmet_sp(ny_fullmet))
          allocate(MR_dx_met(nx_fullmet))
          allocate(MR_dy_met(ny_fullmet))
          dx_met_const = 1.25_sp
          dy_met_const = 1.25_sp
          x_start =   0.0_dp
          y_start =  90.0_dp
          do i = 0,nx_fullmet+1
            x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
          enddo
          x_inverted = .false.
          do i = 1,ny_fullmet
            y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
          enddo
          y_inverted = .true.
          do i = 1,nx_fullmet
            MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
          enddo
          do i = 1,ny_fullmet-1
            MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
          enddo
          MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

        elseif(MR_iwindformat.eq.27)then
          if(MR_iversion.eq.2)then
            !  NOAA-CIRES reanalysis 2.0 degree files  :: ds131.2
            !  Note that v.2 is an older version.
            !  Since v.2 has a regular grid, we just specify it here
            maxdimlen            = 24
            nlev_coords_detected = 2
            allocate(nlevs_fullmet(nlev_coords_detected))
            nlevs_fullmet(1) = 24
            nlevs_fullmet(2) = 19
            allocate(levs_code(nlev_coords_detected))
            levs_code(1) = 1
            levs_code(2) = 2
            allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen))
            np_fullmet = 24
            allocate(p_fullmet_sp(np_fullmet))
            levs_fullmet_sp(:,:) = 0.0_sp
            p_fullmet_sp(1:np_fullmet) = &
              (/1000.0_sp, 950.0_sp, 900.0_sp, 850.0_sp, 800.0_sp, &
                 750.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, 550.0_sp, &
                 500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, 300.0_sp, &
                 250.0_sp, 200.0_sp, 150.0_sp, 100.0_sp,  70.0_sp, &
                  50.0_sp,  30.0_sp,  20.0_sp,  10.0_sp /)
            if(MR_Use_RDA)then
              ! Files from rda.ucar.edu/datasets/ds131.2/ converted from grib with ncl_convert2nc
              ! are top down
              z_inverted = .true.
            else
              ! Files from www.esrl.noaa.gov are bottom up
              z_inverted = .false.
            endif
            levs_fullmet_sp(1,1:24) = p_fullmet_sp(1:24)*Pressure_Conv_Fac
            levs_fullmet_sp(2,1:19) = p_fullmet_sp(1:19)*Pressure_Conv_Fac
            Met_var_zdim_idx( 1) = 1
            Met_var_zdim_idx( 2) = 1
            Met_var_zdim_idx( 3) = 1
            Met_var_zdim_idx( 4) = 2
            Met_var_zdim_idx( 5) = 1
            Met_var_zdim_idx( 7) = 2
  
            IsLatLon_MetGrid  = .true.
            IsGlobal_MetGrid  = .true.
            IsRegular_MetGrid = .true.
            nx_fullmet = 180
            ny_fullmet = 91
            allocate(x_fullmet_sp(0:nx_fullmet+1))
            allocate(y_fullmet_sp(ny_fullmet))
            allocate(MR_dx_met(nx_fullmet))
            allocate(MR_dy_met(ny_fullmet))
            dx_met_const = 2.0_sp
            dy_met_const = 2.0_sp
            x_start =   0.0_dp
            y_start =  90.0_dp
            do i = 0,nx_fullmet+1
              x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
            enddo
            x_inverted = .false.
            do i = 1,ny_fullmet
              y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
            enddo
            y_inverted = .true.
            do i = 1,nx_fullmet
              MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
            enddo
            do i = 1,ny_fullmet-1
              MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
            enddo
            MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)
          elseif(MR_iversion.eq.3)then
            !  NOAA-CIRES reanalysis  :: ds131.3
            !  In v.3, we have a Gaussian grid so we read in the values in
            !  MR_Read_Met_Times_netcdf
            maxdimlen            = 28
            nlev_coords_detected = 2
            allocate(nlevs_fullmet(nlev_coords_detected))
            nlevs_fullmet(1) = 28
            nlevs_fullmet(2) = 21
            allocate(levs_code(nlev_coords_detected))
            levs_code(1) = 1
            levs_code(2) = 2
            allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen))
            np_fullmet = 28
            allocate(p_fullmet_sp(np_fullmet))
            levs_fullmet_sp(:,:) = 0.0_sp
            p_fullmet_sp(1:np_fullmet) = &
              (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
                 850.0_sp, 800.0_sp, 750.0_sp, 700.0_sp, 650.0_sp, &
                 600.0_sp, 550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, &
                 350.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, &
                 100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp,  20.0_sp, &
                  10.0_sp,   5.0_sp,   1.0_sp /)
              ! Files from www.esrl.noaa.gov are bottom up
              z_inverted = .false.
            levs_fullmet_sp(1,1:28) = p_fullmet_sp(1:28)*Pressure_Conv_Fac
            levs_fullmet_sp(2,1:21) = p_fullmet_sp(1:21)*Pressure_Conv_Fac
            Met_var_zdim_idx( 1) = 1
            Met_var_zdim_idx( 2) = 1
            Met_var_zdim_idx( 3) = 1
            Met_var_zdim_idx( 4) = 2
            Met_var_zdim_idx( 5) = 1
            Met_var_zdim_idx( 7) = 2
  
            IsLatLon_MetGrid  = .true.
            IsGlobal_MetGrid  = .true.
            IsRegular_MetGrid = .false.
            allocate(MR_dx_met(nx_fullmet))
            allocate(MR_dy_met(ny_fullmet))
            do i = 1,nx_fullmet
              MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
            enddo
            do i = 1,ny_fullmet-1
              MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
            enddo
            MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)
          endif

        elseif(MR_iwindformat.eq.29)then
          ! ECMWF ERA5
          maxdimlen            = 37
          nlev_coords_detected = 1
          allocate(nlevs_fullmet(nlev_coords_detected))
          nlevs_fullmet(1) = 37
          allocate(levs_code(nlev_coords_detected))
          levs_code(1) = 1
          allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen))
          np_fullmet = 37
          allocate(p_fullmet_sp(np_fullmet))
          levs_fullmet_sp(:,:) = 0.0_sp
          p_fullmet_sp(1:np_fullmet) = &
            (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
               875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
               750.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, 550.0_sp, &
               500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, 300.0_sp, &
               250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
               125.0_sp, 100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp, &
                20.0_sp,  10.0_sp,   7.0_sp,   5.0_sp,   3.0_sp, &
                 2.0_sp, 1.0_sp /)
          z_inverted = .true.
          levs_fullmet_sp(1,1:37) = p_fullmet_sp(1:37)*Pressure_Conv_Fac
          Met_var_zdim_idx( 1) = 1
          Met_var_zdim_idx( 2) = 1
          Met_var_zdim_idx( 3) = 1
          Met_var_zdim_idx( 4) = 1
          Met_var_zdim_idx( 5) = 1
          Met_var_zdim_idx( 7) = 1

          IsLatLon_MetGrid  = .true.
          IsGlobal_MetGrid  = .true.
          IsRegular_MetGrid = .false.
          allocate(MR_dx_met(nx_fullmet))
          allocate(MR_dy_met(ny_fullmet))
          do i = 1,nx_fullmet
            MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
          enddo
          do i = 1,ny_fullmet-1
            MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
          enddo
          MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

        elseif(MR_iwindformat.eq.30)then
          ! ECMWF ERA=20c
          maxdimlen            = 37
          nlev_coords_detected = 1
          allocate(nlevs_fullmet(nlev_coords_detected))
          nlevs_fullmet(1) = 37
          allocate(levs_code(nlev_coords_detected))
          levs_code(1) = 1
          allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen))
          np_fullmet = 37
          allocate(p_fullmet_sp(np_fullmet))
          levs_fullmet_sp(:,:) = 0.0_sp
          p_fullmet_sp(1:np_fullmet) = &
            (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
               875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
               750.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, 550.0_sp, &
               500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, 300.0_sp, &
               250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
               125.0_sp, 100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp, &
                20.0_sp,  10.0_sp,   7.0_sp,   5.0_sp,   3.0_sp, &
                 2.0_sp, 1.0_sp /)
          z_inverted = .true.
          levs_fullmet_sp(1,1:37) = p_fullmet_sp(1:37)*Pressure_Conv_Fac
          Met_var_zdim_idx( 1) = 1
          Met_var_zdim_idx( 2) = 1
          Met_var_zdim_idx( 3) = 1
          Met_var_zdim_idx( 4) = 1
          Met_var_zdim_idx( 5) = 1
          Met_var_zdim_idx( 7) = 1

          IsLatLon_MetGrid  = .true.
          IsGlobal_MetGrid  = .true.
          IsRegular_MetGrid = .false.

          allocate(MR_dx_met(nx_fullmet))
          allocate(MR_dy_met(ny_fullmet))
          do i = 1,nx_fullmet
            MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
          enddo
          do i = 1,ny_fullmet-1
            MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
          enddo
          MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)
        else
          stop 1
        endif

        idx = Met_var_zdim_idx(1)
        p_fullmet_sp(1:nlevs_fullmet(idx)) = levs_fullmet_sp(idx,1:nlevs_fullmet(idx))

      else  ! MR_iwind not equal to 5
        !-----------------------------------------------------------------------------------
        !  Start of all non-iwind=5 cases
        if(MR_iwindformat.eq.50)then
          ! WRF files have a special reader, but we still need to set up 

          call MR_Get_WRF_grid

        else  ! MR_iwindformat .ne. 50
          !---------------------------------------------------------------------------------
          ! Start of block for all non-iwind=5 and non-iwf=50
          ! This is where the Netcdf and Grib subroutines can be compared
          !
          ! Checking for dimension length and values for x,y,t,p
          !   Assume all files have the same format

          maxdimlen = 0
          infile = adjustl(trim(MR_windfiles(1)))
          nSTAT=nf90_open(adjustl(trim(infile)),NF90_NOWRITE, ncid)
          if(nSTAT.ne.NF90_NOERR)call MR_NC_check_status(nSTAT,1,"nf90_open")
          do ivar = 1,MR_MAXVARS
            if (.not.Met_var_IsAvailable(ivar)) cycle  ! Only look at variables that are available
            if (Met_var_ndim(ivar).ne.4) cycle         !  and only ones with a 'level' dimension
            invar = Met_var_NC_names(ivar)
            nSTAT = nf90_inq_varid(ncid,invar,in_var_id)  ! get the var_id for this named variable
            if(nSTAT.ne.NF90_NOERR)then
              call MR_NC_check_status(nSTAT,0,"inq_varid")
              write(MR_global_error,*)'  Cannot find variable ',invar
              write(MR_global_error,*)'  Setting Met_var_IsAvailable to .false.'
              Met_var_IsAvailable(ivar) = .false.
              cycle
            endif
            nSTAT = nf90_inquire_variable(ncid, in_var_id, invar, &
                      xtype = var_xtype, &
                      ndims = var_ndims)   ! get the number of dimensions
            if(nSTAT.ne.NF90_NOERR)call MR_NC_check_status(nSTAT,1,"inq_variable")
            if (var_ndims.ne.Met_var_ndim(ivar))then
              write(MR_global_error,*)'MR ERROR: The actual number of dimensions differs from'
              write(MR_global_error,*)'          what is expected'
              write(MR_global_error,*)'      Variable : ',ivar,Met_var_NC_names(ivar)
              write(MR_global_error,*)'      Expected : ',Met_var_ndim(ivar)
              write(MR_global_error,*)'      Found    : ',var_ndims
              stop 1
            endif
            allocate(var_dimIDs(var_ndims))
            nSTAT = nf90_inquire_variable(ncid, in_var_id, invar, &
                      dimids = var_dimIDs(:var_ndims))
            if(nSTAT.ne.NF90_NOERR)call MR_NC_check_status(nSTAT,1,"inq_variable")
    
            ! 3-d transient variables should be in the COORDS convention (time, level, y, x)
            !                                                                4      3  2  1
            ! if ivar = 1 (Geopotential Height), then take this opportunity to get all
            ! the info on the dimensions for x,y and time.  We will assume these
            ! are the same for all variables.  Pressure levels will be read
            ! separately since they do differ among variables.
            if(ivar.eq.1)then
              i_dim = 1  ! get x info
              nSTAT = nf90_inquire_dimension(ncid,var_dimIDs(i_dim), &
                           name =  dimname, &
                           len = dimlen)
              if(nSTAT.ne.NF90_NOERR)call MR_NC_check_status(nSTAT,1,"nf90_inquire_dimension X")
              if(index(dimname,Met_dim_names(4)).ne.0)then
                nx_fullmet = dimlen
                xdim_id    = var_dimIDs(i_dim)
              endif
              nSTAT = nf90_inq_varid(ncid,dimname,var_id) ! get the variable associated with this dim
              if(nSTAT.ne.NF90_NOERR)call MR_NC_check_status(nSTAT,1,"inq_variable X")
              ! Check if we need to read into a float or a double
              nSTAT = nf90_inquire_variable(ncid, var_id, dimname, xtype = var_xtype)
              if(nSTAT.ne.NF90_NOERR)call MR_NC_check_status(nSTAT,1,"nf90_inquire_variable X")
              allocate(x_fullmet_sp(0:nx_fullmet+1))
              if(var_xtype.eq.NF90_FLOAT)then
                allocate(dum1d_sp(dimlen))
                nSTAT = nf90_get_var(ncid,var_id,dum1d_sp, &
                       start = (/1/),count = (/dimlen/))
                if(nSTAT.ne.NF90_NOERR)call MR_NC_check_status(nSTAT,1,"get_var X flt")
                ! copy to local variable
                x_fullmet_sp(1:nx_fullmet) = dum1d_sp(1:nx_fullmet)
                deallocate(dum1d_sp)
              elseif(var_xtype.eq.NF90_DOUBLE)then
                allocate(dum1d_dp(dimlen))
                nSTAT = nf90_get_var(ncid,var_id,dum1d_dp, &
                       start = (/1/),count = (/dimlen/))
                if(nSTAT.ne.NF90_NOERR)call MR_NC_check_status(nSTAT,1,"get_var X dbl")
                ! copy to local variable
                x_fullmet_sp(1:nx_fullmet) = real(dum1d_dp(1:nx_fullmet),kind=sp)
                deallocate(dum1d_dp)
              else
                write(MR_global_error,*)'MR ERROR: Cannot recognize variable type for x'
                stop 1
              endif
              ! Check the units
              nSTAT = nf90_Inquire_Attribute(ncid, var_id,&
                                             "units",xtype, length, attnum)
              if(nSTAT.ne.NF90_NOERR)then
                call MR_NC_check_status(nSTAT,0,"nf90_Inquire_Attribute")
                write(MR_global_error,*)'MR WARNING: cannot file units ',dimname,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR WARNING: cannot file units ',dimname,nf90_strerror(nSTAT)
              else
                nSTAT = nf90_get_att(ncid, var_id,"units",ustring)
                call MR_NC_check_status(nSTAT,1,"nf90_get_att X units")
                if(index(ustring,'km').gt.0.or.&
                   index(ustring,'kilo').gt.0)then
                  ! This is a projected grid
                  IsLatLon_MetGrid  = .false.
                elseif(index(ustring,'deg').gt.0)then
                  ! This is a lon/lat grid
                  IsLatLon_MetGrid  = .true.
                else
                  write(MR_global_error,*)"MR ERROR: Cannot determine if the grid is lon/lat or projected"
                  stop 1
                endif
              endif
              ! Finally, check for orientation
              if(x_fullmet_sp(1).lt.x_fullmet_sp(2))then
                x_inverted = .false.
              else
                x_inverted = .true.
              endif
              allocate(MR_dx_met(nx_fullmet))
              x_start = x_fullmet_sp(1)
              dx_met_const = x_fullmet_sp(2)-x_fullmet_sp(1)
              x_fullmet_sp(0)            = x_fullmet_sp(1)          - dx_met_const
              x_fullmet_sp(nx_fullmet+1) = x_fullmet_sp(nx_fullmet) + dx_met_const
              if(abs(x_fullmet_sp(nx_fullmet+1)-360.0-x_fullmet_sp(1)).lt.0.1*dx_met_const)then
                IsGlobal_MetGrid = .true.
              else
                IsGlobal_MetGrid = .false.
              endif
              do i = 1,nx_fullmet
                MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
              enddo 
 
              i_dim = 2  ! get y info
              nSTAT = nf90_inquire_dimension(ncid,var_dimIDs(i_dim), &
                           name =  dimname, &
                           len = dimlen)
              call MR_NC_check_status(nSTAT,1,"nf90_inquire_dimension Y")
              if(index(dimname,Met_dim_names(3)).ne.0)then
                ny_fullmet = dimlen
                ydim_id    = var_dimIDs(i_dim)
              endif

              nSTAT = nf90_inq_varid(ncid,dimname,var_id) ! get the variable associated with this dim
              call MR_NC_check_status(nSTAT,1,"nf90_inq_varid Y")
              ! Check if we need to read into a float or a double
              nSTAT = nf90_inquire_variable(ncid, var_id, dimname, xtype = var_xtype)
              call MR_NC_check_status(nSTAT,1,"nf90_inquire_variable Y")
              allocate(y_fullmet_sp(ny_fullmet))
              if(var_xtype.eq.NF90_FLOAT)then
                allocate(dum1d_sp(dimlen))
                nSTAT = nf90_get_var(ncid,var_id,dum1d_sp, &
                       start = (/1/),count = (/dimlen/))
                call MR_NC_check_status(nSTAT,1,"nf90_get_var Y flt")
                ! copy to local variable
                y_fullmet_sp(1:ny_fullmet) = dum1d_sp(1:ny_fullmet)
                deallocate(dum1d_sp)
              elseif(var_xtype.eq.NF90_DOUBLE)then
                allocate(dum1d_dp(dimlen))
                nSTAT = nf90_get_var(ncid,var_id,dum1d_dp, &
                       start = (/1/),count = (/dimlen/))
                call MR_NC_check_status(nSTAT,1,"nf90_get_var Y dbl")
                ! copy to local variable
                y_fullmet_sp(1:ny_fullmet) = real(dum1d_dp(1:ny_fullmet),kind=sp)
                deallocate(dum1d_dp)
              else
                write(MR_global_error,*)'MR ERROR: Cannot recognize variable type for y'
                stop 1
              endif
              ! Check the units
              nSTAT = nf90_Inquire_Attribute(ncid, var_id,&
                                             "units",xtype, length, attnum)
              if(nSTAT.ne.NF90_NOERR)then
                call MR_NC_check_status(nSTAT,0,"nf90_Inquire_Attribute Y units")
                write(MR_global_error,*)'MR WARNING: cannot file units ',dimname,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR WARNING: cannot file units ',dimname,nf90_strerror(nSTAT)
              else
                nSTAT = nf90_get_att(ncid, var_id,"units",ustring)
                call MR_NC_check_status(nSTAT,1,"nf90_get_att Y units")
                if(index(ustring,'km').gt.0.or.&
                   index(ustring,'kilo').gt.0)then
                  ! This is a projected grid
                  IsLatLon_MetGrid  = .false.
                elseif(index(ustring,'deg').gt.0)then
                  ! This is a lon/lat grid
                  IsLatLon_MetGrid  = .true.
                else
                  write(MR_global_error,*)"MR ERROR: Cannot determine if the grid is lon/lat or projected"
                  stop 1
                endif
              endif
              ! check for orientation
              if(y_fullmet_sp(1).lt.y_fullmet_sp(2))then
                y_inverted = .false.
              else
                y_inverted = .true.
              endif
              allocate(MR_dy_met(ny_fullmet))
              y_start = y_fullmet_sp(1)
              dy_met_const = y_fullmet_sp(2)-y_fullmet_sp(1)
              do i = 1,ny_fullmet-1
                MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
              enddo
              MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)
              ! We need to check if this is a regular grid
              IsRegular_MetGrid = .true.
              do i = 1,nx_fullmet-1
                if(abs(MR_dx_met(i+1)-MR_dx_met(i)).gt.tol*MR_dx_met(i))then
                  IsRegular_MetGrid = .false.
                endif
              enddo
              do i = 1,ny_fullmet-1
                if(abs(MR_dy_met(i+1)-MR_dy_met(i)).gt.tol*MR_dy_met(i))then
                  IsRegular_MetGrid = .false.
                endif
              enddo

              i_dim = 4  ! get t info
              nSTAT = nf90_inquire_dimension(ncid,var_dimIDs(i_dim), &
                           name =  dimname, &
                           len = dimlen)
              call MR_NC_check_status(nSTAT,1,"nf90_inquire_dimension time")
              if(index(dimname,Met_dim_names(1)).ne.0)then
                nt_fullmet = dimlen
                tdim_id    = var_dimIDs(i_dim)
              endif
  
            endif ! ivar.eq.1

            ! Dimension variables for x,y,time were read while reading GPH, but
            ! we need the level information for all variables.
            ! Now checking level coordinates (pressure, height, depth);
            ! This is the third dimension for COORD compliant netcdf files.
            i_dim = 3
            nSTAT = nf90_inquire_dimension(ncid,var_dimIDs(i_dim), &
                         name =  dimname, & 
                         len = dimlen)
            call MR_NC_check_status(nSTAT,1,"nf90_inquire_dimension P")
            ! Here is the master list of the known level names. For each
            ! variable, we check if the 'level' coordinate is new or already
            ! noted in this file, building an indexed list
            nSTAT = nf90_inq_varid(ncid,dimname,var_id)
            if(nSTAT.eq.NF90_NOERR.and. &   ! This first condition excludes dims with no vars
               (index(dimname,'lev').ne.0.or.&
                index(dimname,'isobaric').ne.0.or.&
                index(dimname,'isobaricInhPa').ne.0.or.&
                index(dimname,'pressure').ne.0.or.&
                index(dimname,'height').ne.0.or.&
                index(dimname,'depth').ne.0.or.&
                index(dimname,'lv_ISBL1').ne.0.or.&
                index(dimname,'bottom_top').ne.0.or.&
                index(dimname,'bottom_top_stag').ne.0.or.&
                index(dimname,'soil_layers_stag').ne.0))then

              ! Log this level coordinate if it is the first
              if (nlev_coords_detected.eq.0)then
                nlev_coords_detected = nlev_coords_detected + 1
                Met_var_zdim_idx(ivar)  = nlev_coords_detected
                Met_var_zdim_ncid(ivar) = var_dimIDs(i_dim)
                maxdimlen = dimlen
              else
                ! Otherwise, check if this level coordinate has already been logged
                FoundOldDim = .false.
                do iivar = 1,ivar-1
                  if (Met_var_zdim_ncid(iivar).eq.var_dimIDs(i_dim))then
                    FoundOldDim = .true.
                    Met_var_zdim_idx(ivar)  = Met_var_zdim_idx(iivar)
                    Met_var_zdim_ncid(ivar) = var_dimIDs(i_dim)
                    exit
                  endif
                enddo
                if(.not.FoundOldDim)then
                  nlev_coords_detected = nlev_coords_detected + 1
                  Met_var_zdim_idx(ivar)  = nlev_coords_detected
                  Met_var_zdim_ncid(ivar) = var_dimIDs(i_dim)
                  if (maxdimlen.lt.dimlen) maxdimlen = dimlen
                endif
              endif
            else
              call MR_NC_check_status(nSTAT,0,"nf90_inq_varid")
              write(MR_global_error,*)'MR ERROR: level coordinate is not in pos. 3 for ',invar
              write(MR_global_error,*)'          Expected one of: lev, isobaric, pressure,'
              write(MR_global_error,*)'            height, depth, lv_ISBL1, bottom_top,'
              write(MR_global_error,*)'            bottom_top_stag, soil_layers_stag'
              write(MR_global_error,*)'          Instead, found: ',dimname
              write(MR_global_log  ,*)'MR ERROR: level coordinate is not in pos. 3 for ',invar
              stop 1
            endif
            ! tidy up
            deallocate(var_dimIDs)
          enddo ! ivar
  
          ! We have all the level dimension names and dim_ids; now we need to get the sizes
          allocate(nlevs_fullmet(nlev_coords_detected))
          allocate(levs_code(nlev_coords_detected))
          allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen))
          do ivar = 1,MR_MAXVARS
            ! Check if this variable has a z-dimension (pressure, height, depth, etc.)
            if(Met_var_zdim_ncid(ivar).gt.0)then
              ! log the length of the dimension for this level coordinate
              nSTAT = nf90_inquire_dimension(ncid,Met_var_zdim_ncid(ivar), &
                         name =  dimname, &
                         len = dimlen)
              call MR_NC_check_status(nSTAT,1,"nf90_inquire_dimension level")
              idx = Met_var_zdim_idx(ivar)
              nlevs_fullmet(idx) = dimlen
              ! Now inquire and populate the dimension variable info
              nSTAT = nf90_inq_varid(ncid,dimname,var_id)
              call MR_NC_check_status(nSTAT,1,"nf90_inq_varid level")
              ! Check if we need to read into a float or a double
              nSTAT = nf90_inquire_variable(ncid, var_id, dimname, xtype = var_xtype)
              call MR_NC_check_status(nSTAT,1,"nf90_inquire_variable level")
              if(var_xtype.eq.NF90_FLOAT)then
                allocate(dum1d_sp(dimlen))
                nSTAT = nf90_get_var(ncid,var_id,dum1d_sp, &
                       start = (/1/),count = (/dimlen/))
                call MR_NC_check_status(nSTAT,1,"nf90_get_var level flt")
                ! copy to local variable
                levs_fullmet_sp(idx,1:nlevs_fullmet(idx)) = dum1d_sp(1:nlevs_fullmet(idx))
                deallocate(dum1d_sp)
              elseif(var_xtype.eq.NF90_DOUBLE)then
                allocate(dum1d_dp(dimlen))
                nSTAT = nf90_get_var(ncid,var_id,dum1d_dp, &
                       start = (/1/),count = (/dimlen/))
                call MR_NC_check_status(nSTAT,1,"nf90_get_var level dbl")
                ! copy to local variable
                levs_fullmet_sp(idx,1:nlevs_fullmet(idx)) = real(dum1d_dp(1:nlevs_fullmet(idx)),kind=sp)
                deallocate(dum1d_dp)
              endif
              ! Check the units
              nSTAT = nf90_Inquire_Attribute(ncid, var_id,&
                                             "units",xtype, length, attnum)
              call MR_NC_check_status(nSTAT,0,"nf90_Inquire_Attribute level units")
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: cannot find dim units ',dimname,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: cannot find dim units ',dimname,nf90_strerror(nSTAT)
                IsPressureDimension = .false.
                stop 1
              else
                nSTAT = nf90_get_att(ncid, var_id,"units",ustring)
                call MR_NC_check_status(nSTAT,1,"nf90_get_att level units")

                ! Note: the variables below are single-valued, not arrays on ivar
                !       If a pressure is found, the assumption here is that all pressure coordinates will
                !       be given in the same units (hPa or Pa) and the same orientations (bot to top, or
                !       inverted.
                if(index(ustring,'Pa').gt.0.or.&
                   index(ustring,'millibar').gt.0)then
                  ! This is a pressure level
                  IsPressureDimension = .true.
                  if(index(ustring,'hPa').gt.0.or.&
                     index(ustring,'millibar').gt.0)then
                    Pressure_Conv_Fac = 100.0_sp
                  else
                    Pressure_Conv_Fac = 1.0_sp
                  endif
                elseif(index(ustring,'level').gt.0)then
                  ! this is a special case for CAM files which are on hybrid levels
                  Pressure_Conv_Fac = 100.0_sp
                  IsPressureDimension = .true.
                else
                  IsPressureDimension = .false.
                endif
              endif
              
              ! Finally, check for orientation
              if(IsPressureDimension.and. &     ! We are only concerned with orientation in pressure
                 nlevs_fullmet(idx).gt.1)then   ! Neglect single-valued pressure coordinates
                if(levs_fullmet_sp(idx,1).lt.levs_fullmet_sp(idx,2))then
                  z_inverted = .true.
                else
                  z_inverted = .false.
                endif
              endif
            endif
          enddo  ! ivar

          ! Close file
          nSTAT = nf90_close(ncid)
          call MR_NC_check_status(nSTAT,1,"nf90_close")
        endif ! else part of (MR_iwindformat.eq.50)
        !-----------------------------------------

        ! Now report on what we've found
        write(MR_global_production,*)" Found these levels"
        write(MR_global_production,*)"  VaribleID    LevelIdx       dimID      length"
        do ivar = 1,MR_MAXVARS
          if (Met_var_IsAvailable(ivar))then
            if(Met_var_zdim_idx(ivar).eq.0)then
              write(MR_global_production,*)ivar,Met_var_zdim_idx(ivar),&
                                           Met_var_zdim_ncid(ivar),0
            else
              write(MR_global_production,*)ivar,Met_var_zdim_idx(ivar),&
                                           Met_var_zdim_ncid(ivar),&
                                           nlevs_fullmet(Met_var_zdim_idx(ivar))
            endif
          endif
        enddo
#ifdef USEPOINTERS
        if(.not.associated(p_fullmet_sp))then
#else
        if(.not.allocated(p_fullmet_sp))then
#endif            
          ! Now invert if necessary and convert to Pa
          allocate(p_fullmet_sp(maxdimlen))
          do idx = 1,nlev_coords_detected
            if(z_inverted)then
              do i = 1,nlevs_fullmet(idx)
                p_fullmet_sp(nlevs_fullmet(idx)+1-i) = levs_fullmet_sp(idx,i)*Pressure_Conv_Fac
              enddo
            else
              p_fullmet_sp(1:nlevs_fullmet(idx)) = levs_fullmet_sp(idx,1:nlevs_fullmet(idx))*Pressure_Conv_Fac
            endif
            levs_fullmet_sp(idx,:) = 0.0_sp
            levs_fullmet_sp(idx,1:nlevs_fullmet(idx)) = p_fullmet_sp(1:nlevs_fullmet(idx))
          enddo
          deallocate(p_fullmet_sp)
        endif
 
        ! The native vertical coordinate used in MetReader is that used for GPH.
        ! Here we check each of the vertical coordinates and note how they
        ! compare with that of GPH, assigning the following codes:
        !  1 = same dim length: assume one-to-one mapping
        !  2 = fewer levels: truncated version of GPH (doesn't extend as high)
        !  3 = fewer levels: missing levels requiring interpolation
        !  4 = more values than that of GPH
        levs_code(1:nlev_coords_detected) = 0
        levs_code(1) = 1                       ! The first var checked (GPH) should have a one-to-one mapping
        ! Check how each of the pressure coordinates map onto the GPH grid
        if (nlev_coords_detected.gt.1)then
          ! Only bother if there are multiple pressure coordinates
          do idx = 2,nlev_coords_detected
            if (nlevs_fullmet(idx).gt.nlevs_fullmet(1))then
              ! This coordinate has more values than the GPH pressure coordinate
              levs_code(idx) = 4
            elseif (nlevs_fullmet(idx).lt.nlevs_fullmet(1))then
              ! It there are fewer levels, check if this is a truncated coordinate (code = 2)
              ! or one with missing levels that requires interpolation (code = 3)
              IsTruncatedDim = .true.
              do i=1,nlevs_fullmet(idx)
                if(abs(levs_fullmet_sp(idx,i)-levs_fullmet_sp(1,i)).gt.MR_EPS_SMALL)then
                  IsTruncatedDim = .false.
                  exit
                endif
              enddo
              if(IsTruncatedDim)then
                levs_code(idx) = 2
              else
                levs_code(idx) = 3
              endif
            else
              ! This coordinate has the same dimension as the GPH pressure coordinate.
              ! They are probably the same
              levs_code(idx) = 1
            endif
          enddo
        endif
  
        ! Now assign these levels to the working arrays
        np_fullmet    = nlevs_fullmet(Met_var_zdim_idx(1))  ! Assign fullmet the length of H,U,V
  
        ! Geopotential
#ifdef USEPOINTERS        
        if(.not.associated(p_fullmet_sp))  allocate(p_fullmet_sp(np_fullmet))
#else
        if(.not.allocated(p_fullmet_sp))   allocate(p_fullmet_sp(np_fullmet))
#endif
        idx = Met_var_zdim_idx(1)
        p_fullmet_sp(1:nlevs_fullmet(idx)) = levs_fullmet_sp(idx,1:nlevs_fullmet(idx))
  
      endif ! else part of (MR_iwind.eq.5)
      !---------------------------------------------------------------------------------

      if(MR_iwindformat.eq.0)then
        ! Template windfile (example for nam198)
        !  Need to populate
        call MR_Set_Met_Dims_Template_netcdf
      endif

      allocate(z_approx(np_fullmet))
      do k=1,np_fullmet
        ! Calculate heights for US Std Atmos while pressures are still in mbars
        ! or hPa
        z_approx(k) = MR_Z_US_StdAtm(p_fullmet_sp(k))
      enddo
      MR_Max_geoH_metP_predicted = z_approx(np_fullmet)

      write(MR_global_info,*)"Dimension info:"
      write(MR_global_info,*)"  record (time): ",nt_fullmet
      write(MR_global_info,*)"  level  (z)   : ",np_fullmet
      write(MR_global_info,*)"  y            : ",ny_fullmet
      write(MR_global_info,*)"  x            : ",nx_fullmet

      !************************************************************************
      ! assign boundaries of mesoscale model
      if(x_inverted)then
          ! I know of no windfiles with x-coordinate reversed
        xLL_fullmet = x_fullmet_sp(nx_fullmet)
        xUR_fullmet = x_fullmet_sp(1)
      else
        xLL_fullmet = x_fullmet_sp(1)
        xUR_fullmet = x_fullmet_sp(nx_fullmet)
      endif

      if(y_inverted)then
          ! Most lon/lat grids have y reversed
        yLL_fullmet = y_fullmet_sp(ny_fullmet)
        yUR_fullmet = y_fullmet_sp(1)
      else
          ! Projected grids have y not reversed
        yLL_fullmet = y_fullmet_sp(1)
        yUR_fullmet = y_fullmet_sp(ny_fullmet)
      endif

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      end subroutine MR_Read_Met_DimVars_netcdf

!##############################################################################


!##############################################################################
!
!     MR_Get_WRF_grid
!
!     Called once from MR_Read_Met_DimVars_netcdf
!
!     This subroutine reads the variable and dimension IDs, and fills the
!     coordinate dimension variables, just like MR_Read_Met_DimVars_netcdf,
!     but WRF files are a bit more complicated than most since the grid/projection can
!     vary and state variables (U,V,W,Gph) are each on grids staggered in the
!     relevant direction.
!
!     After this subroutine completes, the following variables will be set:
!       All the projection parameters of NWP grid
!       Met_dim_names, Met_var_NC_names, Met_var_conversion_factor, Met_var_IsAvailable
!       The lengths of all the dimensions of the file
!       p_fullmet_sp (converted to Pa)
!       x_fullmet_sp, y_fullmet_sp
!       IsLatLon_MetGrid, IsGlobal_MetGrid, IsRegular_MetGrid 
!
!##############################################################################

      subroutine MR_Get_WRF_grid

      use MetReader
      use netcdf
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer :: nSTAT
      integer :: ncid

      character(len = nf90_max_name) :: name_dum

      integer :: t_dim_id          = 0 ! x or lon
      integer :: x_dim_id          = 0 ! x or lon
      integer :: y_dim_id          = 0 ! y or lat
      integer :: z_dim_id          = 0 ! x or lon

      integer :: lon_var_id
      integer :: lat_var_id
      integer :: PB_var_id
      integer :: Ppert_var_id

      integer :: Map_Proj
      real(kind=sp) :: WRF_dx,WRF_dy
      real(kind=sp) :: Cen_Lat,Cen_Lon
      real(kind=sp) :: Moad_Cen_Lat
      real(kind=sp) :: Stand_Lon
      real(kind=sp) :: Truelat1,Truelat2
      real(kind=sp) :: Pole_Lat, Pole_Lon
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dum3d_sp
      real(kind=sp),dimension(:,:,:,:),allocatable :: dum4d_sp
      integer :: i,j

      real(kind=dp) :: x_start,y_start

      real(kind=dp) :: lat_in,lon_in

        ! MAP_PROJ - Model projection  1 = Lambert
        !                              2 = polar stereographic
        !                              3 = Mercator
        !                              6 = lat-lon

        ! First set spatial (x/y) grid
        ! Open first windfile and assume all grids are the same

      write(MR_global_info,*)"About to open first WRF file : ",MR_windfiles(1)
      nSTAT=nf90_open(adjustl(trim(MR_windfiles(1))),NF90_NOWRITE, ncid)
      call MR_NC_check_status(nSTAT,1,"nf90_open WRF file")
      
      Met_dim_names(1) = "Time"             ! time
      Met_dim_names(2) = "bottom_top"       ! pressure (24 levels 10 -> 1000)
      Met_dim_names(3) = "south_north"      ! y        (90.0 -> -90.0)
      Met_dim_names(4) = "west_east"        ! x        (0.0 -> 378.0)
      Met_dim_names(5) = "bottom_top_stag"  ! Stag in Z (pressure coordinate for Vz)
      Met_dim_names(6) = "west_east_stag"   ! Stag in x
      Met_dim_names(7) = "south_north_stag" ! Stag in y

      ! Get dim ids and sizes
      nSTAT = nf90_inq_dimid(ncid,Met_dim_names(1),t_dim_id)
      call MR_NC_check_status(nSTAT,1,"nf90_inq_dimid Time")
      nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,name=name_dum,len=nt_fullmet)
      call MR_NC_check_status(nSTAT,1,"nf90_Inquire_Dimension Time")
      nSTAT = nf90_inq_dimid(ncid,Met_dim_names(4),x_dim_id)
      call MR_NC_check_status(nSTAT,1,"nf90_inq_dimid X")
      nSTAT = nf90_Inquire_Dimension(ncid,x_dim_id,name=name_dum,len=nx_fullmet)
      call MR_NC_check_status(nSTAT,1,"nf90_Inquire_Dimension X")
      nSTAT = nf90_inq_dimid(ncid,Met_dim_names(3),y_dim_id)
      call MR_NC_check_status(nSTAT,1,"nf90_inq_dimid Y")
      nSTAT = nf90_Inquire_Dimension(ncid,y_dim_id,name=name_dum,len=ny_fullmet)
      call MR_NC_check_status(nSTAT,1,"nf90_Inquire_Dimension Y")
      nSTAT = nf90_inq_dimid(ncid,Met_dim_names(2),z_dim_id)
      call MR_NC_check_status(nSTAT,1,"nf90_inq_dimid level")
      nSTAT = nf90_Inquire_Dimension(ncid,z_dim_id,name=name_dum,len=neta_fullmet)
      call MR_NC_check_status(nSTAT,1,"nf90_Inquire_Dimension level")
      nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "MAP_PROJ", Map_Proj)
      call MR_NC_check_status(nSTAT,1,"nf90_get_att MAP_PROJ")
      if(Map_Proj.eq.1)then
         ! Lambert
         !  Mandatory parameters are given in module_map_utils.F line 297:
         !      ' truelat1, truelat2, lat1, lon1, knowni, knownj, stdlon, dx'
         !   TRUELAT1            = lat_1 # Cone intersects with the sphere
         !   TRUELAT2 (optional) = lat_2 # Cone intersects with the sphere
         !   STAND_LON           = lon_0 # Center Point
         !  Other global attributes:
         !   MOAD_CEN_LAT        = lat_0 # Center Point
         !   CEN_LAT
         !   CEN_LON
         !   POLE_LAT
         !   POLE_LON
         !proj +proj=lcc +lon_0=-175.0 +lat_0=55.0 +lat_1=50.0 +lat_2=60.0 +R=6371.229

        write(MR_global_info,*)"  WRF projection detected: Lambert Conformal"

        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "DX", WRF_dx)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att DX")
        ! Note: DY might not be present, in which case DY=DX
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "DY", WRF_dy)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att DY")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "MOAD_CEN_LAT", Moad_Cen_Lat)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att MOAD_CEN_LAT")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "STAND_LON", Stand_Lon)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att STAND_LON")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT1", Truelat1)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att TRUELAT1")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT2", Truelat2)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att TRUELAT2")

          ! convert dx, dy to km
        IsLatLon_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        dx_met_const = WRF_dx*1.0e-3_4
        dy_met_const = WRF_dy*1.0e-3_4

        ! Projected grids have Lon and Lat provided as 2d fields
        allocate(Met_Proj_lat(nx_fullmet,ny_fullmet))
        allocate(Met_Proj_lon(nx_fullmet,ny_fullmet))
        allocate(dum3d_sp(nx_fullmet,ny_fullmet,1))

        nSTAT = nf90_inq_varid(ncid,"XLONG",lon_var_id)
        call MR_NC_check_status(nSTAT,1,"nf90_inq_varid XLONG")
        nSTAT = nf90_inq_varid(ncid,"XLAT",lat_var_id)
        call MR_NC_check_status(nSTAT,1,"nf90_inq_varid XLAT")

        nSTAT = nf90_get_var(ncid,lon_var_id,dum3d_sp, &
               start = (/1,1,1/),count = (/nx_fullmet,ny_fullmet,1/))
        call MR_NC_check_status(nSTAT,1,"nf90_get_var XLONG")
        Met_Proj_lon(:,:) = dum3d_sp(:,:,1)

        nSTAT = nf90_get_var(ncid,lat_var_id,dum3d_sp, &
               start = (/1,1,1/),count = (/nx_fullmet,ny_fullmet,1/))
        call MR_NC_check_status(nSTAT,1,"nf90_get_var XLAT")
        Met_Proj_lat(:,:) = dum3d_sp(:,:,1)

        ! In the example WRF files, x and y projected values are not actually
        ! provided, so we recreate them here using the coordinates if the LL
        ! point of the Lon/Lat grid
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))

        lon_in = real(Met_Proj_lon(1,1),kind=8)
        lat_in = real(Met_Proj_lat(1,1),kind=8)

          ! Setting the projection parameters as libprojection.a expects
        Met_iprojflag = 4  
        Met_lam0   = real(Stand_Lon,kind=8)
        Met_phi0   = real(Moad_Cen_Lat,kind=8)
        Met_phi1   = real(Truelat1,kind=8)
        Met_phi2   = real(Truelat2,kind=8)
        Met_k0     = real(1.0,kind=8)
        Met_Re     = 6370.0_8
        call PJ_proj_for(lon_in,lat_in, &
                       Met_iprojflag,Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                       x_start,y_start)
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=4)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=4)
        enddo
        do i = 1,nx_fullmet-1
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        MR_dx_met(nx_fullmet)    = MR_dx_met(nx_fullmet-1)
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(Map_Proj.eq.2)then
        ! Polar Stereographic
         !  Mandatory parameters are given in module_map_utils.F line 308:
         !      ' truelat1, lat1, lon1, knowni, knownj, stdlon, dx'
         !      and in subroutine set_ps on line 661
         !   TRUELAT1            = lat_1 # Cone intersects with the sphere
         !   STAND_LON           = lon_0 # Center Point
         !  Other global attributes:
         !   MOAD_CEN_LAT        = lat_0 # Center Point
         !   CEN_LAT
         !   CEN_LON
         !   POLE_LAT
         !   POLE_LON
         !   TRUELAT2 (optional) = lat_2 # Cone intersects with the sphere
         !proj +proj=lcc +lon_0=-175.0 +lat_0=55.0 +lat_1=50.0 +lat_2=60.0 +R=6371.229

        write(MR_global_info,*)"  WRF projection detected: Polar Stereographic"

        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "DX", WRF_dx)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att DX")
        ! Note: DY might not be present, in which case DY=DX
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "DY", WRF_dy)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att DY")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "MOAD_CEN_LAT", Moad_Cen_Lat)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att MOAD_CEN_LAT")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "STAND_LON", Stand_Lon)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att STAND_LON")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT1", Truelat1)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att TRUELAT1")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "POLE_LAT", Pole_Lat)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att POLE_LAT")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "POLE_LON", Pole_Lon)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att POLE_LON")

          ! convert dx, dy to km
        IsLatLon_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        dx_met_const = WRF_dx*1.0e-3_4
        dy_met_const = WRF_dy*1.0e-3_4

        ! Projected grids have Lon and Lat provided as 2d fields
        allocate(Met_Proj_lat(nx_fullmet,ny_fullmet))
        allocate(Met_Proj_lon(nx_fullmet,ny_fullmet))
        allocate(dum3d_sp(nx_fullmet,ny_fullmet,1))

        nSTAT = nf90_inq_varid(ncid,"XLONG",lon_var_id)
        call MR_NC_check_status(nSTAT,1,"nf90_inq_varid XLONG")
        nSTAT = nf90_inq_varid(ncid,"XLAT",lat_var_id)
        call MR_NC_check_status(nSTAT,1,"nf90_inq_varid XLAT")

        nSTAT = nf90_get_var(ncid,lon_var_id,dum3d_sp, &
               start = (/1,1,1/),count = (/nx_fullmet,ny_fullmet,1/))
        call MR_NC_check_status(nSTAT,1,"nf90_get_var XLONG")
        Met_Proj_lon(:,:) = dum3d_sp(:,:,1)

        nSTAT = nf90_get_var(ncid,lat_var_id,dum3d_sp, &
               start = (/1,1,1/),count = (/nx_fullmet,ny_fullmet,1/))
        call MR_NC_check_status(nSTAT,1,"nf90_get_var XLAT")
        Met_Proj_lat(:,:) = dum3d_sp(:,:,1)

        ! In the example WRF files, x and y projected values are not actually
        ! provided, so we recreate them here using the coordinates if the LL
        ! point of the Lon/Lat grid
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))

        lon_in = real(Met_Proj_lon(1,1),kind=8)
        lat_in = real(Met_Proj_lat(1,1),kind=8)

          ! Setting the projection parameters as libprojection.a expects
        Met_iprojflag = 1
        Met_lam0   = real(Stand_Lon,kind=8)
        Met_phi0   = real(-Pole_Lat,kind=8) ! phi0 is the lat of projection origin
        Met_phi1   = real(Truelat1,kind=8)
        Met_phi2   = real(Truelat2,kind=8)
        !Met_k0     = real(0.972759,kind=8)
        Met_k0     = real(1.0,kind=8)
        Met_Re     = 6370.0_8
        call PJ_proj_for(lon_in,lat_in, &
                       Met_iprojflag,Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                       x_start,y_start)
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=4)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=4)
        enddo
        do i = 1,nx_fullmet-1
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        MR_dx_met(nx_fullmet)    = MR_dx_met(nx_fullmet-1)
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(Map_Proj.eq.3)then
        ! Mercator
         !   TRUELAT1            = lat_1 # Cone intersects with the sphere
         !  Other global attributes:
         !   MOAD_CEN_LAT        = lat_0 # Center Point
         !   STAND_LON           = lon_0 # Center Point
         !   CEN_LAT
         !   CEN_LON
         !   POLE_LAT
         !   POLE_LON
         !   TRUELAT2 (optional) = lat_2 # Cone intersects with the sphere


         !proj +proj=merc 

        write(MR_global_info,*)"  WRF projection detected: Mercator"

        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "DX", WRF_dx)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att DX")
        ! Note: DY might not be present, in which case DY=DX
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "DY", WRF_dy)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att DY")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "CEN_LAT", Cen_Lat)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att CEN_LAT")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "STAND_LON", Stand_Lon)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att STAND_LON")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT1", Truelat1)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att TRUELAT1")
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT2", Truelat2)
        call MR_NC_check_status(nSTAT,1,"nf90_get_att TRUELAT2")

          ! convert dx, dy to km
        IsLatLon_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        dx_met_const = WRF_dx*1.0e-3_4
        dy_met_const = WRF_dy*1.0e-3_4

        ! Projected grids have Lon and Lat provided as 2d fields
        allocate(Met_Proj_lat(nx_fullmet,ny_fullmet))
        allocate(Met_Proj_lon(nx_fullmet,ny_fullmet))
        allocate(dum3d_sp(nx_fullmet,ny_fullmet,1))

        nSTAT = nf90_inq_varid(ncid,"XLONG",lon_var_id)
        call MR_NC_check_status(nSTAT,1,"nf90_inq_varid XLONG")
        nSTAT = nf90_inq_varid(ncid,"XLAT",lat_var_id)
        call MR_NC_check_status(nSTAT,1,"nf90_inq_varid XLAT")

        nSTAT = nf90_get_var(ncid,lon_var_id,dum3d_sp, &
               start = (/1,1,1/),count = (/nx_fullmet,ny_fullmet,1/))
        call MR_NC_check_status(nSTAT,1,"nf90_get_var XLONG")
        Met_Proj_lon(:,:) = dum3d_sp(:,:,1)
        nSTAT = nf90_get_var(ncid,lat_var_id,dum3d_sp, &
               start = (/1,1,1/),count = (/nx_fullmet,ny_fullmet,1/))
        call MR_NC_check_status(nSTAT,1,"nf90_get_var XLAT")
        Met_Proj_lat(:,:) = dum3d_sp(:,:,1)

        ! In the example WRF files, x and y projected values are not actually
        ! provided, so we recreate them here using the coordinates if the LL
        ! point of the Lon/Lat grid
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))

        lon_in = real(Met_Proj_lon(1,1),kind=8)
        lat_in = real(Met_Proj_lat(1,1),kind=8)

          ! Setting the projection parameters as libprojection.a expects
        Met_iprojflag = 5
        Met_lam0   = real(Stand_Lon,kind=8)
        Met_phi0   = real(Cen_Lat,kind=8)
        Met_phi1   = real(Truelat1,kind=8)
        Met_phi2   = real(Truelat2,kind=8)
        Met_k0     = real(1.0,kind=8)
        Met_Re     = 6370.0_8
        call PJ_proj_for(lon_in,lat_in, &
                       Met_iprojflag,Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                       x_start,y_start)
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=4)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=4)
        enddo
        do i = 1,nx_fullmet-1
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        MR_dx_met(nx_fullmet)    = MR_dx_met(nx_fullmet-1)
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(Map_Proj.eq.6)then
        ! Lon-Lat or cylindrical equidistant
        !   pole_lat
        !   pole_lon
        !   stand_lon

        write(MR_global_info,*)"  WRF projection detected: Lon-Lat"

        write(MR_global_info,*)"WRF: MAP_PROJ=6 : Lon-Lat : Not implemented"
        stop 1
      else
        write(MR_global_info,*)&
         "The MAP_PROJ global attribute is either not present or is"
        write(MR_global_info,*)"not a recognized projection."
        stop 1
      endif

      ! Now setting up pressure coordinate
      ! WRF data are provided on eta levels instead of pressure level
      ! Fortunately, it will usually be adequate to pretend that they are on
      ! pressure levels since Ash3d interpolates onto a z-grid using GPH.
      ! The "atmosphere" module, will need access to physical pressure and will
      ! need a special case for WRF files
      np_fullmet    = neta_fullmet
      allocate(p_fullmet_sp(np_fullmet))
      allocate(dum4d_sp(nx_fullmet,ny_fullmet,np_fullmet,1))

      ! To populate a place-holder p_fullmet_sp, read the full pressure grid and
      ! copy a representative column
      nSTAT = nf90_inq_varid(ncid,"PB",PB_var_id)
      call MR_NC_check_status(nSTAT,1,"nf90_inq_varid PB")

      nSTAT = nf90_get_var(ncid,PB_var_id,dum4d_sp, &
             start = (/1,1,1,1/),count = (/nx_fullmet,ny_fullmet,neta_fullmet,1/))
      call MR_NC_check_status(nSTAT,1,"nf90_get_var PB")
      p_fullmet_sp(:) = dum4d_sp(1,1,:,1)
      nSTAT = nf90_inq_varid(ncid,"P",Ppert_var_id)
      call MR_NC_check_status(nSTAT,1,"nf90_inq_varid P")
      nSTAT = nf90_get_var(ncid,Ppert_var_id,dum4d_sp, &
             start = (/1,1,1,1/),count = (/nx_fullmet,ny_fullmet,neta_fullmet,1/))
      call MR_NC_check_status(nSTAT,1,"nf90_get_var P")
      p_fullmet_sp(:) = p_fullmet_sp(:) + dum4d_sp(1,1,:,1)
      MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)/100.0_sp) 
      !p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa

      x_inverted = .false.
      y_inverted = .false.
      z_inverted = .false.

      ! Close file
      nSTAT = nf90_close(ncid)
      call MR_NC_check_status(nSTAT,1,"nf90_close WRF file")

      nlev_coords_detected = 1
      allocate(nlevs_fullmet(nlev_coords_detected))
      nlevs_fullmet(1) = np_fullmet
      allocate(levs_code(nlev_coords_detected))
      levs_code(1) = 1
      allocate(levs_fullmet_sp(nlev_coords_detected,np_fullmet))
      levs_fullmet_sp(1,1:np_fullmet) = p_fullmet_sp(1:np_fullmet)

      Met_var_zdim_idx(:)  = 1
      !Met_var_zdim_ncid(ivar) = var_dimIDs(i_dim)

      end subroutine MR_Get_WRF_grid

!##############################################################################
!
!     MR_Read_Met_Times_netcdf
!
!     Called once from MR_Read_Met_DimVars 
!
!     This subroutine opens each netcdf file and determine the time of each
!     time step of each file in the number of hours since MR_BaseYear.
!     In most cases, the length of the time variable (nt_fullmet) will be 
!     read directly from the file and overwritten (is was set in MR_Read_Met_DimVars_netcdf
!     above).
!
!     After this subroutine completes, the following variables will be set:
!       MR_windfile_starthour(MR_iwindfiles)
!       MR_windfile_stephour(MR_iwindfiles,nt_fullmet)
!
!##############################################################################

      subroutine MR_Read_Met_Times_netcdf

      use MetReader
      use netcdf

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer :: iw,iws
      integer :: itstart_year,itstart_month
      integer :: itstart_day
      real(kind=sp) :: filestart_hour

      integer :: itstart_hour,itstart_min,itstart_sec

      integer :: nSTAT
      integer :: ncid
      integer :: time_var_id = 0
      integer :: gph_var_id  = 0
      integer :: reftime_var_id
      integer :: t_dim_id
      integer :: x_dim_id,y_dim_id,x_var_id,y_var_id
      integer :: reftimedimID
      integer :: var_ndims
      integer,dimension(:),allocatable :: var_dimIDs
      integer :: gph_ndims
      integer,dimension(NF90_MAX_VAR_DIMS) :: gph_DimIDs

      integer :: reftimedimlen
      real(kind=sp),dimension(:),allocatable :: filetime_in_sp
      character(len=19) :: Timestr_WRF

      integer            :: var_xtype
      character(len=NF90_MAX_NAME) :: invar
      character(len=NF90_MAX_NAME) :: indim
      integer            :: xtype, length, attnum
      character(len=31)  :: tstring2
      !real(kind=8)       :: HS_hours_since_baseyear !,HS_HourOfDay
      real(kind=8)       :: iwf_int,iwf_tot
      integer            :: iwstep
      logical            :: TimeHasUnitsAttr = .false.
      integer            :: i,ii
      real(kind=dp),dimension(:), allocatable :: dum1d_dp
      real(kind=sp),dimension(:), allocatable :: dum1d_sp
      integer(kind=4),dimension(:), allocatable :: dum1d_int4
      integer,dimension(8)  :: values
      integer               :: Current_Year,nt_tst
      character(len=130)    :: Z_infile
      !character(len=1)      :: answer          !for debugging
      !integer               :: HS_YearOfEvent
      !integer               :: HS_MonthOfEvent
      !integer               :: HS_DayOfEvent

      INTERFACE
        real(kind=8) function HS_hours_since_baseyear(iyear,imonth,iday,hours,byear,useLeaps)
          integer            :: iyear
          integer            :: imonth
          integer            :: iday
          real(kind=8)       :: hours
          integer            :: byear
          logical            :: useLeaps
        end function HS_hours_since_baseyear
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
      END INTERFACE

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------                MR_Read_Met_Times_netcdf                    ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      if(.not.Met_dim_IsAvailable(1))then
        write(MR_global_error,*)"MR ERROR: Time dimension is required and not listed"
        write(MR_global_error,*)"          in template windfile specification file."
        stop 1
      endif

      allocate(MR_windfile_starthour(MR_iwindfiles))
      if(MR_iwind.eq.5)then
        ! For iw=5, we need to know the start and end day,month,year
        MR_Comp_StartYear        = HS_YearOfEvent(MR_Comp_StartHour,MR_BaseYear,MR_useLeap)
        MR_Comp_StartMonth       = HS_MonthOfEvent(MR_Comp_StartHour,MR_BaseYear,MR_useLeap)
        MR_Comp_StartDay         = HS_DayOfEvent(MR_Comp_StartHour,MR_BaseYear,MR_useLeap)

        ! Here the branch for when MR_iwindformat = 25, 26, 27, 29, or 30
        ! First copy path read in to slot 2
        !MR_iw5_root = MR_windfiles(1)
 110    format(a50,a1,i4,a1)
        write(MR_windfiles(1),110)trim(ADJUSTL(MR_iw5_root)),'/', &
                                   MR_Comp_StartYear,'/'
        if(MR_iwindformat.eq.25)then
          iwf_int = 6.0_dp
          iwf_tot = MR_iw5_hours_per_file
        elseif(MR_iwindformat.eq.26)then
          iwf_int = 6.0_dp
          iwf_tot = MR_iw5_hours_per_file
        elseif(MR_iwindformat.eq.27)then
          if(MR_iversion.eq.2)then
            iwf_int = 6.0_dp
            iwf_tot = MR_iw5_hours_per_file
          elseif(MR_iversion.eq.3)then
            iwf_int = 3.0_dp
            iwf_tot = MR_iw5_hours_per_file
          endif
        elseif(MR_iwindformat.eq.29)then
          iwf_int = 1.0_dp
          iwf_tot = MR_iw5_hours_per_file
        elseif(MR_iwindformat.eq.30)then
          iwf_int = 3.0_dp
          iwf_tot = MR_iw5_hours_per_file
        endif
        nt_fullmet = 1
        do iw = 1,MR_iwindfiles
          if(MR_iwindformat.eq.25)then
            MR_windfile_starthour(iw) = &
              real(HS_hours_since_baseyear(MR_Comp_StartYear+(iw-1),1,&
                                           1,0.0_8,MR_BaseYear,MR_useLeap),kind=sp)
          elseif(MR_iwindformat.eq.26)then
            MR_windfile_starthour(iw) = &
              real(HS_hours_since_baseyear(MR_Comp_StartYear,MR_Comp_StartMonth+(iw-1),&
                                           1,0.0_8,MR_BaseYear,MR_useLeap),kind=sp)
          elseif(MR_iwindformat.eq.27)then
            MR_windfile_starthour(iw) = &
              real(HS_hours_since_baseyear(MR_Comp_StartYear+(iw-1),1,&
                                           1,0.0_8,MR_BaseYear,MR_useLeap),kind=sp)
          elseif(MR_iwindformat.eq.29)then
            MR_windfile_starthour(iw) = &
              real(HS_hours_since_baseyear(MR_Comp_StartYear,MR_Comp_StartMonth,&
                                           MR_Comp_StartDay+(iw-1),0.0_8,MR_BaseYear,MR_useLeap),kind=sp)
          elseif(MR_iwindformat.eq.30)then
            MR_windfile_starthour(iw) = &
              real(HS_hours_since_baseyear(MR_Comp_StartYear,MR_Comp_StartMonth+(iw-1),&
                                           1,0.0_8,MR_BaseYear,MR_useLeap),kind=sp)
          endif

          ! Building the name of the first windfile (for hgt) to inspect for nt
          call MR_Set_iwind5_filenames(MR_Comp_StartHour+(iw-1)*iwf_tot,1,Z_infile)
          nSTAT = nf90_open(trim(ADJUSTL(Z_infile)),NF90_NOWRITE,ncid)
          if(nSTAT.ne.NF90_NOERR)then
            call MR_NC_check_status(nSTAT,0,"nf90_open")
            if(iw.eq.1)then
              ! Do a hard stop if we can't even read the first file
              write(MR_global_error,*)'MR ERROR: nf90_open: ',nf90_strerror(nSTAT)
              write(MR_global_error,*)"    Could not open file: ",trim(ADJUSTL(Z_infile))
              write(MR_global_log  ,*)'MR ERROR: nf90_open: ',nf90_strerror(nSTAT)
              write(MR_global_error,*)'Exiting'
              stop 1
            else
              ! This is probably OK as long as the sim time is within the first file
              write(MR_global_error,*)'MR WARNING: nf90_open: ',nf90_strerror(nSTAT)
              write(MR_global_error,*)"    Could not open file: ",trim(ADJUSTL(Z_infile))
              write(MR_global_error,*)"    This should be OK if the previous file exits."
              exit
            endif
          endif
          nSTAT = nf90_inq_dimid(ncid,Met_dim_names(1),t_dim_id)
          call MR_NC_check_status(nSTAT,1,"nf90_inq_dimid time")
          nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,len=nt_tst)
          call MR_NC_check_status(nSTAT,1,"nf90_Inquire_Dimension time")
          if(iw.eq.1.and.(.not.IsRegular_MetGrid))then ! for iwf=5, this is 27,29, or 30
            ! Normally we would populate the x and y arrays in MR_Read_Met_DimVars_netcdf, but
            ! for Gaussian or otherwise irregular grids, it is easier to just read the grids
            ! directly.  We will do this now while we have the Geopotential Height file open.  
            nSTAT = nf90_inq_dimid(ncid,Met_dim_names(3),y_dim_id)
            call MR_NC_check_status(nSTAT,1,"nf90_inq_dimid Y (lat)")
            nSTAT = nf90_Inquire_Dimension(ncid,y_dim_id,len=ny_fullmet)
            call MR_NC_check_status(nSTAT,1,"nf90_Inquire_Dimension Y (lat)")
            allocate(y_fullmet_sp(ny_fullmet))
            allocate(dum1d_dp(ny_fullmet))
            nSTAT = nf90_inq_varid(ncid,Met_dim_names(3),y_var_id)
            call MR_NC_check_status(nSTAT,1,"nf90_inq_varid Y (lat)")
            nSTAT = nf90_get_var(ncid,y_var_id,dum1d_dp, &
                     start = (/1/),count = (/ny_fullmet/))
            call MR_NC_check_status(nSTAT,1,"nf90_get_var Y (lat)")
            y_fullmet_sp(1:ny_fullmet) = real(dum1d_dp(1:ny_fullmet),kind=sp)
            y_inverted = .true.
            deallocate(dum1d_dp)

            nSTAT = nf90_inq_dimid(ncid,Met_dim_names(4),x_dim_id)
            call MR_NC_check_status(nSTAT,1,"nf90_inq_dimid X (lon)")
            nSTAT = nf90_Inquire_Dimension(ncid,x_dim_id,len=nx_fullmet)
            call MR_NC_check_status(nSTAT,1,"nf90_Inquire_Dimension X (lon)")
            allocate(x_fullmet_sp(0:nx_fullmet+1))
            allocate(dum1d_dp(nx_fullmet))
            nSTAT = nf90_inq_varid(ncid,Met_dim_names(4),x_var_id)
            call MR_NC_check_status(nSTAT,1,"nf90_inq_varid X (lon)")
            nSTAT = nf90_get_var(ncid,x_var_id,dum1d_dp, &
                     start = (/1/),count = (/nx_fullmet/))
            call MR_NC_check_status(nSTAT,1,"nf90_get_var X (lon)")
            x_fullmet_sp(1:nx_fullmet) = real(dum1d_dp(1:nx_fullmet),kind=sp)
            x_inverted = .false.
            x_fullmet_sp(0)            = 2.0_sp*x_fullmet_sp(1)-x_fullmet_sp(2)
            x_fullmet_sp(nx_fullmet+1) = 2.0_sp*x_fullmet_sp(nx_fullmet) - &
                                                x_fullmet_sp(nx_fullmet-1)
            deallocate(dum1d_dp)

          endif ! .not.IsRegular_MetGrid
          nSTAT = nf90_close(ncid)
          call MR_NC_check_status(nSTAT,1,"nf90_close")

           ! Set this (and subsequent) iw to the current number of timesteps
          MR_windfiles_nt_fullmet(iw:) = nt_tst
          if(nt_tst.gt.nt_fullmet) nt_fullmet = nt_tst
        enddo
        
        if(MR_iwindformat.eq.25)then
          ! Getting current year
          call date_and_time(VALUES=values)
          Current_Year = values(1)
          if(MR_Comp_StartYear.lt.Current_Year.and.nt_tst.lt.nt_fullmet)then
            write(MR_global_info,*)"WARNING:  The NCEP files are for an archived year yet are incomplete."
            write(MR_global_info,*)"          To get the complete year, run the script "
            write(MR_global_info,*)"            autorun_scripts/get_NCEP_50YearReanalysis.sh",MR_Comp_StartYear
            write(MR_global_info,*)"          Steps available = ",nt_tst
            write(MR_global_info,*)"          Hours into year = ",(nt_tst-1)*6
          endif
        endif
        allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))

          ! the interval for both iwf25 and iwf27 is 6 hours
        do iwstep = 1,nt_fullmet
          MR_windfile_stephour(:,iwstep) = (iwstep-1)*iwf_int
        enddo
      else ! MR_iwind = 3 or 4
        if(MR_iwindformat.eq.50)then
          ! Branch for WRF files
          ! Loop through all the windfiles
          do iw = 1,MR_iwindfiles
            nSTAT = nf90_open(trim(ADJUSTL(MR_windfiles(iw))),NF90_NOWRITE,ncid)
            call MR_NC_check_status(nSTAT,0,"nf90_open WRF")
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: nf90_open to read header:', &
                             nf90_strerror(nSTAT)
              write(MR_global_error,*)'Could not open ',trim(ADJUSTL(MR_windfiles(iw)))
              write(MR_global_error,*)'Exiting'
              stop 1
            endif
            if(iw.eq.1)then
              ! Find the id of the time dimension
              nSTAT = nf90_inq_dimid(ncid,trim(ADJUSTL(Met_dim_names(1))),t_dim_id)
              call MR_NC_check_status(nSTAT,1,"nf90_inq_dimid time")
              ! Get length of time dimension and allocate MR_windfile_stephour
              nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,len=nt_fullmet)
              call MR_NC_check_status(nSTAT,1,"nf90_Inquire_Dimension time")
              write(MR_global_info,*)"  Assuming all NWP files have the same number of steps."
              write(MR_global_info,*)"   Allocating time arrays for ",MR_iwindfiles,"file(s)"
              write(MR_global_info,*)"                              ",nt_fullmet,"step(s) each"
              allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))
              MR_windfile_stephour(:,:) = 0.0_dp
              nSTAT = nf90_inq_varid(ncid,"Times",time_var_id)
              call MR_NC_check_status(nSTAT,1,"nf90_inq_varid Times")
              nSTAT = nf90_inquire_variable(ncid, time_var_id, invar, &
                  xtype = var_xtype)
              call MR_NC_check_status(nSTAT,1,"nf90_inquire_variable Times")
              allocate(filetime_in_sp(nt_fullmet))
              if(nt_fullmet.gt.1)then
                write(MR_global_error,*)"MR ERROR: Currently WRF files are expected to only have one"
                write(MR_global_error,*)"       timestep/file"
                stop 1
              endif
              filetime_in_sp = 0.0_sp
            endif
            nSTAT = nf90_get_var(ncid,time_var_id,Timestr_WRF,&
                           start = (/1,1/),       &
                           count = (/19,1/))
            call MR_NC_check_status(nSTAT,1,"nf90_get_var Times")

            nSTAT = nf90_close(ncid)
            call MR_NC_check_status(nSTAT,1,"nf90_close WRF")

            read(Timestr_WRF,121)itstart_year,itstart_month,itstart_day, &
                              itstart_hour,itstart_min,itstart_sec
            filestart_hour = real(itstart_hour,kind=sp) + &
                             real(itstart_min,kind=sp)/60.0_sp      + &
                             real(itstart_sec,kind=sp)/3600.0_sp
 121        format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2,1x)
            MR_windfiles_nt_fullmet(iw)=nt_fullmet
            MR_windfile_starthour(iw) = real(HS_hours_since_baseyear(itstart_year,itstart_month, &
                                         itstart_day,real(filestart_hour,kind=8),MR_BaseYear,MR_useLeap),kind=4)
            do iwstep = 1,nt_fullmet
              MR_windfile_stephour(iw,iwstep) = MR_windfile_stephour(iw,1) + filetime_in_sp(iwstep)
            enddo
          enddo
        else
          ! For all other formats, try to read the GRIB_orgReferenceTime string
          ! Loop through all the windfiles
          do iw = 1,MR_iwindfiles

            ! Each wind file needs a ref-time which in almost all cases is given
            ! in the 'units' attribute of the time variable
            write(MR_global_info,*)iw,trim(ADJUSTL(MR_windfiles(iw)))
            nSTAT = nf90_open(trim(ADJUSTL(MR_windfiles(iw))),NF90_NOWRITE,ncid)
            call MR_NC_check_status(nSTAT,0,"nf90_open")
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: nf90_open to read header:', &
                             nf90_strerror(nSTAT)
              write(MR_global_error,*)'Could not open ',trim(ADJUSTL(MR_windfiles(iw)))
              write(MR_global_error,*)'Exiting'
              stop 1
            endif

            ! Find the id of the time dimension
            !   First, search for the GPH variable and check its
            !   dimensions
            nSTAT = nf90_inq_varid(ncid,Met_var_NC_names(1),gph_var_id)
            call MR_NC_check_status(nSTAT,1,"nf90_inq_varid GPH")
            nSTAT = nf90_inquire_variable(ncid, gph_var_id, ndims = gph_ndims)
            call MR_NC_check_status(nSTAT,1,"nf90_inquire_variable GPH")
            nSTAT = nf90_inquire_variable(ncid, gph_var_id, dimids = gph_DimIDs(:gph_ndims))
            call MR_NC_check_status(nSTAT,1,"nf90_inquire_variable GPH")
            if(gph_ndims.lt.4)then
              write(MR_global_error,*)'MR ERROR:'
              write(MR_global_error,*)&
                'GPH variable does not have 4 dimensions.'
              write(MR_global_error,*)&
                'Expecting time, level, y, x'
              stop 1
            endif
            nSTAT = nf90_inquire_dimension(ncid,gph_DimIDs(1),indim)
            call MR_NC_check_status(nSTAT,1,"nf90_inquire_dimension X")
            Met_dim_names(4)= trim(adjustl(indim))
            nSTAT = nf90_inquire_dimension(ncid,gph_DimIDs(2),indim)
            call MR_NC_check_status(nSTAT,1,"nf90_inquire_dimension Y")
            Met_dim_names(3)= trim(adjustl(indim))
            nSTAT = nf90_inquire_dimension(ncid,gph_DimIDs(3),indim)
            call MR_NC_check_status(nSTAT,1,"nf90_inquire_dimension level")
            Met_dim_names(2)= trim(adjustl(indim))
            nSTAT = nf90_inquire_dimension(ncid,gph_DimIDs(4),indim)
            call MR_NC_check_status(nSTAT,1,"nf90_inquire_dimension Time")
            Met_dim_names(1)= trim(adjustl(indim))

            nSTAT = nf90_inq_dimid(ncid,trim(ADJUSTL(Met_dim_names(1))),t_dim_id)
            call MR_NC_check_status(nSTAT,1,"nf90_inq_dimid Time")
            if(iw.eq.1)then
              ! Get length of time dimension and allocate MR_windfile_stephour
              nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,len=nt_fullmet)
              call MR_NC_check_status(nSTAT,1,"nf90_inquire_dimension Time")
              write(MR_global_info,*)"  Assuming all NWP files have the same number of steps."
              write(MR_global_info,*)"   Allocating time arrays for ",MR_iwindfiles,"files"
              write(MR_global_info,*)"                              ",nt_fullmet,"step(s) each"
              allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))
            endif

            ! get variable id for time
            nSTAT = nf90_inq_varid(ncid,trim(ADJUSTL(Met_dim_names(1))),time_var_id)
            call MR_NC_check_status(nSTAT,1,"nf90_inq_varid Time")
            ! We need the reftime for this file, check time variable for 'units'
            nSTAT = nf90_Inquire_Attribute(ncid, time_var_id,&
                                           "units",xtype, length, attnum)
            if(nSTAT.eq.0)then
              TimeHasUnitsAttr = .true.
              nSTAT = nf90_get_att(ncid, time_var_id,"units",tstring2)
              call MR_NC_check_status(nSTAT,0,"nf90_get_att time units")
              reftimedimlen = 31  ! set the length to the full amount
            else
              ! if the call to check units above failed, issue a warning
              call MR_NC_check_status(nSTAT,0,"nf90_Inquire_Attribute")
              write(MR_global_log  ,*)'  Trying to read GRIB_orgReferenceTime'
              ! Try GRIB_orgReferenceTime
              nSTAT = nf90_Inquire_Attribute(ncid, time_var_id,&
                                             "GRIB_orgReferenceTime",xtype, length, attnum)
              call MR_NC_check_status(nSTAT,1,"nf90_Inquire_Attribute Time GRIB_orgReferenceTime")
              nSTAT = nf90_get_att(ncid, time_var_id,"GRIB_orgReferenceTime",tstring2)
              call MR_NC_check_status(nSTAT,1,"nf90_get_att Time GRIB_orgReferenceTime")

              TimeHasUnitsAttr = .true.
            endif

            if(TimeHasUnitsAttr)then
              if(index(tstring2,'since').ne.0)then
                ! loop over the first 26 characters looking for 'since'
                do i=1,26
                  ! try to parse
                  !  time:units = "Hour since 2016-01-11T00:00:00Z" ;
                  !  time:units = "days since 0001-01-01 00:00:00" ;
                  if(tstring2(i:i+5).eq.'since ')then
                    ii = i+6
                    read(tstring2(ii:31),103)itstart_year,itstart_month,itstart_day, &
                                      itstart_hour,itstart_min,itstart_sec
                    write(MR_global_info,2100)"Ref time = ",itstart_year,itstart_month,itstart_day, &
                                               itstart_hour,itstart_min,itstart_sec
                    filestart_hour = real(itstart_hour,kind=sp) + &
                                     real(itstart_min,kind=sp)/60.0_sp      + &
                                     real(itstart_sec,kind=sp)/3600.0_sp
                    exit
                  endif
                enddo
              endif
            endif

            if(index(tstring2,'since').eq.0)then
              ! Time variable does not have units attribute
              ! Try variable 'reftime'
              nSTAT = nf90_inq_varid(ncid,'reftime',reftime_var_id)
              call MR_NC_check_status(nSTAT,0,"nf90_inq_varid reftime")
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)"MR ERROR:  Could not read time:units or reftime"
                write(MR_global_error,*)"        Windfile start time is not defined."
                stop 1
              endif

              var_ndims = 1
              allocate(var_dimIDs(1))
              nSTAT = nf90_inquire_variable(ncid, reftime_var_id, invar, &
                        dimids = var_dimIDs(:var_ndims))
              call MR_NC_check_status(nSTAT,1,"nf90_inquire_variable reftime")
              reftimedimID = var_dimIDs(1)
              nSTAT = nf90_inquire_dimension(ncid,reftimedimID, &
                           len = reftimedimlen)
              call MR_NC_check_status(nSTAT,1,"nf90_inquire_dimension reftime")
              nSTAT = nf90_get_var(ncid,reftime_var_id,tstring2(:reftimedimlen))
              call MR_NC_check_status(nSTAT,0,"nf90_get_var reftime")
              if(nSTAT.ne.0)then
                write(MR_global_error,*)"MR ERROR:  Could not read reftime"
                write(MR_global_error,*)"        Windfile start time is not defined."
                stop 1
              endif
              if(index(tstring2,'20').ne.0.or.index(tstring2,'19').ne.0)then
                do i=1,reftimedimlen-1
                  if(tstring2(i:i+1).eq.'20'.or.tstring2(i:i+1).eq.'19')then
                    write(MR_global_info,*)"Found reference time: ",tstring2(i:reftimedimlen)
                    read(tstring2(i:reftimedimlen),103)itstart_year,itstart_month,itstart_day, &
                                      itstart_hour,itstart_min,itstart_sec
                    write(MR_global_info,2100)"Ref time = ",itstart_year,itstart_month,itstart_day, &
                                               itstart_hour,itstart_min,itstart_sec
                    filestart_hour = real(itstart_hour,kind=sp) + &
                                     real(itstart_min,kind=sp)/60.0_sp      + &
                                     real(itstart_sec,kind=sp)/3600.0_sp
                    exit
                  endif
                enddo
              endif
            endif
2100        format(20x,a11,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)
 103        format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)

            ! Assume we now have the parsed reftime

            ! Now get time data
            ! Check if we need to read into an int, float or a double
            nSTAT = nf90_inquire_variable(ncid, time_var_id, name = invar, xtype = var_xtype)
            call MR_NC_check_status(nSTAT,1,"nf90_inquire_variable reftime")
            if(var_xtype.eq.NF90_FLOAT)then
              allocate(dum1d_sp(nt_fullmet))
              nSTAT = nf90_get_var(ncid,time_var_id,dum1d_sp, &
                     start = (/1/),count = (/nt_fullmet/))
              call MR_NC_check_status(nSTAT,1,"nf90_get_var reftime flt")
              ! copy to local variable
              MR_windfile_stephour(iw,1:nt_fullmet) = dum1d_sp(1:nt_fullmet)* &
                                                           Met_dim_fac(1)
              deallocate(dum1d_sp)
            elseif(var_xtype.eq.NF90_DOUBLE)then
              allocate(dum1d_dp(nt_fullmet))
              nSTAT = nf90_get_var(ncid,time_var_id,dum1d_dp, &
                     start = (/1/),count = (/nt_fullmet/))
              call MR_NC_check_status(nSTAT,1,"nf90_get_var reftime dbl")
              ! copy to local variable
              MR_windfile_stephour(iw,1:nt_fullmet) = real(dum1d_dp(1:nt_fullmet),kind=4)* &
                                                           Met_dim_fac(1)
              deallocate(dum1d_dp)
            elseif(var_xtype.eq.NF90_INT)then
              allocate(dum1d_int4(nt_fullmet))
              nSTAT = nf90_get_var(ncid,time_var_id,dum1d_int4, &
                     start = (/1/),count = (/nt_fullmet/))
              call MR_NC_check_status(nSTAT,1,"nf90_get_var reftime int")
              ! copy to local variable
              MR_windfile_stephour(iw,1:nt_fullmet) = real(dum1d_int4(1:nt_fullmet),kind=4)* &
                                                           Met_dim_fac(1)
              deallocate(dum1d_int4)
            else
              write(MR_global_error,*)"MR ERROR: Unexpected time variable type ",Met_dim_names(i)
              stop 1
            endif

            nSTAT = nf90_close(ncid)
            call MR_NC_check_status(nSTAT,1,"nf90_close")

            MR_windfiles_nt_fullmet(iw) = nt_fullmet
            MR_windfile_starthour(iw) =  real(HS_hours_since_baseyear(itstart_year,itstart_month, &
                                           itstart_day,real(filestart_hour,kind=8),MR_BaseYear,MR_useLeap),kind=4)
          enddo
        endif  ! MR_iwindformat = 50 v.s. all others
      endif  ! MR_iwind = 5 v.s. 3/4
      ! Finished setting up the start time of each wind file in HoursSince : MR_windfile_starthour(iw)
      !  and the forecast (offset from start of file) for each step        : MR_windfile_stephour(iw,iwstep)

      if (MR_iwind.ne.5)then
        write(MR_global_info,*)"  File,  step,        Ref,     Offset,  HoursSince"
        do iw = 1,MR_iwindfiles
          do iws = 1,nt_fullmet
            write(MR_global_info,800)iw,iws,real(MR_windfile_starthour(iw),kind=4),&
                             real(MR_windfile_stephour(iw,iws),kind=4),&
                             real(MR_windfile_starthour(iw)+MR_windfile_stephour(iw,iws),kind=4)
          enddo
        enddo
      endif
 800  format(i7,i7,3f12.2)

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      end subroutine MR_Read_Met_Times_netcdf

!##############################################################################


!##############################################################################
!
!     MR_Set_iwind5_filenames
!
!     Called from MR_Read_Met_Times_netcdf
!
!     Sets the name of the iwind=5 file that contains the data for the requested
!     time as well as the next filename
!
!##############################################################################

      subroutine MR_Set_iwind5_filenames(inhour,ivar,infile)

      use MetReader

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      real(kind=8)      ,intent(in)  :: inhour
      integer           ,intent(in)  :: ivar
      character(len=130),intent(out) :: infile

      !integer               :: HS_YearOfEvent
      !integer               :: HS_MonthOfEvent
      !integer               :: HS_DayOfEvent
      !logical               :: HS_IsLeapYear
      integer,dimension(12) :: DaysInMonth
      integer               :: dum_i1,dum_i2,dum_i3

      integer :: thisYear,thisMonth,thisDay
      logical :: Use_ds630 = .false.

      INTERFACE
        logical function HS_IsLeapYear(iyear)
          integer            :: iyear
        end function HS_IsLeapYear
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
      END INTERFACE

      if(ivar.ne.1.and. &
         ivar.ne.2.and. &
         ivar.ne.3.and. &
         ivar.ne.4.and. &
         ivar.ne.5.and. &
         ivar.ne.7)then
        write(MR_global_error,*)"MR ERROR: iwind=5 only compatible with the following variables:"
        write(MR_global_error,*)"  ivar = 1 :: ",Met_var_NC_names(1)
        write(MR_global_error,*)"  ivar = 2 :: ",Met_var_NC_names(2)
        write(MR_global_error,*)"  ivar = 3 :: ",Met_var_NC_names(3)
        write(MR_global_error,*)"  ivar = 4 :: ",Met_var_NC_names(4)
        write(MR_global_error,*)"  ivar = 5 :: ",Met_var_NC_names(5)
        write(MR_global_error,*)"  ivar = 7 :: ",Met_var_NC_names(4)
      endif

      thisYear        = HS_YearOfEvent( inhour,MR_BaseYear,MR_useLeap)
      thisMonth       = HS_MonthOfEvent(inhour,MR_BaseYear,MR_useLeap)
      thisDay         = HS_DayOfEvent(  inhour,MR_BaseYear,MR_useLeap)

      if(HS_IsLeapYear(thisYear))then
        DaysInMonth=(/31,29,31,30,31,30,31,31,30,31,30,31 /)
      else
        DaysInMonth=(/31,28,31,30,31,30,31,31,30,31,30,31 /)
      endif
      if(MR_iwindformat.eq.25)then   ! Not needed for iwf=25
        ! YYYY/hgt.2018.nc
        dum_i1 = 1                                 ! Start day in file
        dum_i2 = DaysInMonth(thisMonth)   ! End day in file
        dum_i3 = 18                                ! End hour in file
        if(ivar.eq.1)then
          write(MR_iw5_prefix ,251)'hgt.'
        elseif(ivar.eq.2)then
          write(MR_iw5_prefix ,252)'uwnd.'
        elseif(ivar.eq.3)then
          write(MR_iw5_prefix ,253)'vwnd.'
        elseif(ivar.eq.4)then
          write(MR_iw5_prefix ,254)'omega.'
        elseif(ivar.eq.5)then
          write(MR_iw5_prefix ,255)'air.'
        elseif(ivar.eq.7)then
          write(MR_iw5_prefix ,254)'omega.'
        endif
        write(MR_iw5_suffix1,325)thisYear,'.nc'
        write(MR_iw5_suffix2,325)thisYear+1,'.nc'   ! Next file for iwf=25 is next year
        write(infile,425)trim(ADJUSTL(MR_iw5_root)),'/',thisYear,'/', &
                         trim(adjustl(MR_iw5_prefix)),   &
                         trim(adjustl(MR_iw5_suffix1))
 251    format(a4)
 252    format(a5)
 253    format(a5)
 254    format(a6)
 255    format(a4)
 325    format(i4,a3)
 425    format(a50,a1,i4,a1,a,a7)
      elseif(MR_iwindformat.eq.26)then
        ! YYYY/anl_p125.007_hgt.2018060100_2018063018.nc
        dum_i1 = 1                                 ! Start day in file
        dum_i2 = DaysInMonth(thisMonth)   ! End day in file
        dum_i3 = 18                                ! End hour in file
        if(ivar.eq.1)then
          write(MR_iw5_prefix ,261)'anl_p125.007_hgt.'
        elseif(ivar.eq.2)then
          write(MR_iw5_prefix ,262)'anl_p125.033_ugrd.'
        elseif(ivar.eq.3)then
          write(MR_iw5_prefix ,263)'anl_p125.034_vgrd.'
        elseif(ivar.eq.4)then
          write(MR_iw5_prefix ,264)'anl_p125.039_vvel.'
        elseif(ivar.eq.5)then
          write(MR_iw5_prefix ,265)'anl_p125.011_tmp.'
        elseif(ivar.eq.7)then
          write(MR_iw5_prefix ,264)'anl_p125.039_vvel.'
        endif
        write(MR_iw5_suffix1,326)thisYear,thisMonth,dum_i1,'00_',&
                                 thisYear,thisMonth,dum_i2,dum_i3,'.nc'
        write(MR_iw5_suffix2,326)thisYear,thisMonth,dum_i1,'00_',&
                                 thisYear,thisMonth,dum_i2,dum_i3,'.nc'

        write(infile,426)trim(ADJUSTL(MR_iw5_root)),'/',MR_Comp_StartYear,'/', &
                         trim(adjustl(MR_iw5_prefix)),   &
                         trim(adjustl(MR_iw5_suffix1))
 261    format(a17)
 262    format(a18)
 263    format(a18)
 264    format(a18)
 265    format(a17)
 326    format(i4,i0.2,i0.2,a3,i4,i0.2,i0.2,i0.2,a3)
 426    format(a50,a1,i4,a1,a,a24)
      elseif(MR_iwindformat.eq.27)then  ! Not needed for iwf=27
        if(MR_iversion.eq.2)then
          dum_i1 = 1                        ! Start day in file
          dum_i2 = DaysInMonth(thisMonth)   ! End day in file
          dum_i3 = 18                       ! End hour in file
          if(MR_Use_RDA)then
            if(ivar.eq.1)then
              write(MR_iw5_prefix ,271)'pgrbanl_mean_',thisYear,'_HGT'
            elseif(ivar.eq.2)then
              write(MR_iw5_prefix ,271)'pgrbanl_mean_',thisYear,'_UGRD'
            elseif(ivar.eq.3)then
              write(MR_iw5_prefix ,271)'pgrbanl_mean_',thisYear,'_VGRD'
            elseif(ivar.eq.4)then
              write(MR_iw5_prefix ,271)'pgrbanl_mean_',thisYear,'_VVEL'
            elseif(ivar.eq.5)then
              write(MR_iw5_prefix ,271)'pgrbanl_mean_',thisYear,'_TMP'
            elseif(ivar.eq.7)then
              write(MR_iw5_prefix ,271)'pgrbanl_mean_',thisYear,'_VVEL'
            endif
            write(MR_iw5_suffix1,272)'_pres.nc'
          else
            if(ivar.eq.1)then  ! These are the same as for iwf=25 so just use the format
                               ! statements above
              write(MR_iw5_prefix ,251)'hgt.'
            elseif(ivar.eq.2)then
              write(MR_iw5_prefix ,252)'uwnd.'
            elseif(ivar.eq.3)then
              write(MR_iw5_prefix ,253)'vwnd.'
            elseif(ivar.eq.4)then
              write(MR_iw5_prefix ,254)'omega.'
            elseif(ivar.eq.5)then
              write(MR_iw5_prefix ,255)'air.'
            elseif(ivar.eq.7)then
              write(MR_iw5_prefix ,254)'omega.'
            endif
            write(MR_iw5_suffix1,327)thisYear,'.nc'
          endif

          write(infile,427)trim(ADJUSTL(MR_iw5_root)),'/',thisYear,'/', &
                           trim(adjustl(MR_iw5_prefix)),   &
                           trim(adjustl(MR_iw5_suffix1))
        elseif(MR_iversion.eq.3)then
          dum_i1 = 1                        ! Start day in file
          dum_i2 = DaysInMonth(thisMonth)   ! End day in file
          dum_i3 = 18                       ! End hour in file
          if(ivar.eq.1)then
            write(MR_iw5_prefix ,273)'anl_mean_',thisYear,'_HGT'
          elseif(ivar.eq.2)then
            write(MR_iw5_prefix ,273)'anl_mean_',thisYear,'_UGRD'
          elseif(ivar.eq.3)then
            write(MR_iw5_prefix ,273)'anl_mean_',thisYear,'_VGRD'
          elseif(ivar.eq.4)then
            write(MR_iw5_prefix ,273)'anl_mean_',thisYear,'_VVEL'
          elseif(ivar.eq.5)then
            write(MR_iw5_prefix ,273)'anl_mean_',thisYear,'_TMP'
          elseif(ivar.eq.7)then
            write(MR_iw5_prefix ,273)'anl_mean_',thisYear,'_VVEL'
          endif
        write(MR_iw5_suffix1,272)'_pres.nc'

        write(infile,427)trim(ADJUSTL(MR_iw5_root)),'/',thisYear,'/', &
                         trim(adjustl(MR_iw5_prefix)),   &
                         trim(adjustl(MR_iw5_suffix1))
        endif

 271    format(a13,i4,a)
 272    format(a8)
 273    format(a9,i4,a)
 327    format(i4,a3)
 427    format(a50,a1,i4,a1,a,a)
      elseif(MR_iwindformat.eq.29)then
        dum_i1 = thisDay                  ! Start day in file
        dum_i2 = thisDay                  ! End day in file
        dum_i3 = 23                                ! End hour in file
        if (Use_ds630) then
          ! Dataset https://rda.ucar.edu/datasets/ds630.0/
          !  uses the following format:
          ! YYYY/e5.oper.an.pl.128_129_z.regn320sc.2018062000_2018062023.nc
          if(ivar.eq.1)then
            write(MR_iw5_prefix,291)'e5.oper.an.pl.128_129_z.regn320sc.'
          elseif(ivar.eq.2)then
            write(MR_iw5_prefix,291)'e5.oper.an.pl.128_131_u.regn320uv.'
          elseif(ivar.eq.3)then
            write(MR_iw5_prefix,291)'e5.oper.an.pl.128_132_v.regn320uv.'
          elseif(ivar.eq.4)then
            write(MR_iw5_prefix,291)'e5.oper.an.pl.128_135_w.regn320sc.'
          elseif(ivar.eq.5)then
            write(MR_iw5_prefix,291)'e5.oper.an.pl.128_130_t.regn320sc.'
          elseif(ivar.eq.7)then
            write(MR_iw5_prefix,291)'e5.oper.an.pl.128_135_w.regn320sc.'
          endif
        else
          !  whereas https://rda.ucar.edu/datasets/ds633.0
          !  uses the following format:
          ! YYYY/e5.oper.an.pl.128_129_z.ll025sc.2018062000_2018062023.nc
          if(ivar.eq.1)then
            write(MR_iw5_prefix ,292)'e5.oper.an.pl.128_129_z.ll025sc.'
          elseif(ivar.eq.2)then
            write(MR_iw5_prefix ,292)'e5.oper.an.pl.128_131_u.ll025uv.'
          elseif(ivar.eq.3)then
            write(MR_iw5_prefix ,292)'e5.oper.an.pl.128_132_v.ll025uv.'
          elseif(ivar.eq.4)then
            write(MR_iw5_prefix ,292)'e5.oper.an.pl.128_135_w.ll025sc.'
          elseif(ivar.eq.5)then
            write(MR_iw5_prefix ,292)'e5.oper.an.pl.128_130_t.ll025sc.'
          elseif(ivar.eq.7)then
            write(MR_iw5_prefix ,292)'e5.oper.an.pl.128_135_w.ll025sc.'
          endif
        endif
        write(MR_iw5_suffix1,329)thisYear,thisMonth,dum_i1,'00_',&
                                 thisYear,thisMonth,dum_i2,dum_i3,'.nc'
        write(infile,429)trim(ADJUSTL(MR_iw5_root)),'/',MR_Comp_StartYear,'/', &
                         trim(adjustl(MR_iw5_prefix)),   &
                         trim(adjustl(MR_iw5_suffix1))
 291    format(a34)
 292    format(a32)
 329    format(i4,i0.2,i0.2,a3,i4,i0.2,i0.2,i0.2,a3)
 429    format(a50,a1,i4,a1,a,a24)
      elseif(MR_iwindformat.eq.30)then
        ! YYYY/e20c.oper.an.pl.3hr.128_129_z.regn80sc.1912060100_1912063021.nc
        dum_i1 = 1                                 ! Start day in file
        dum_i2 = DaysInMonth(thisMonth)   ! End day in file
        dum_i3 = 21                                ! End hour in file
        if(ivar.eq.1)then
          write(MR_iw5_prefix ,230)'e20c.oper.an.pl.3hr.128_129_z.regn80sc.'
        elseif(ivar.eq.2)then
          write(MR_iw5_prefix ,230)'e20c.oper.an.pl.3hr.128_131_u.regn80uv.'
        elseif(ivar.eq.3)then
          write(MR_iw5_prefix ,230)'e20c.oper.an.pl.3hr.128_132_v.regn80uv.'
        elseif(ivar.eq.4)then
          write(MR_iw5_prefix ,230)'e20c.oper.an.pl.3hr.128_135_w.regn80sc.'
        elseif(ivar.eq.5)then
          write(MR_iw5_prefix ,230)'e20c.oper.an.pl.3hr.128_130_t.regn80sc.'
        elseif(ivar.eq.7)then
          write(MR_iw5_prefix ,230)'e20c.oper.an.pl.3hr.128_135_w.regn80sc.'
        endif
        write(MR_iw5_suffix1,330)thisYear,thisMonth,dum_i1,'00_',&
                                 thisYear,thisMonth,dum_i2,dum_i3,'.nc'
        write(infile,430)trim(ADJUSTL(MR_iw5_root)),'/',MR_Comp_StartYear,'/', &
                         trim(adjustl(MR_iw5_prefix)),   &
                         trim(adjustl(MR_iw5_suffix1))
 230    format(a39)
 330    format(i4,i0.2,i0.2,a3,i4,i0.2,i0.2,i0.2,a3)
 430    format(a50,a1,i4,a1,a,a24)
      endif

      end subroutine MR_Set_iwind5_filenames


!##############################################################################


!##############################################################################
!
!     MR_Set_Met_Dims_Template_netcdf
!
!     Called from MR_Set_MetComp_Grids_netcdf
!
!     Sets Met grid for Template windfiles
!
!##############################################################################


      subroutine MR_Set_Met_Dims_Template_netcdf

      use MetReader
      use netcdf
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision
      real(kind=sp), parameter :: tol = 1.0e-7_sp

      integer :: iw,i

      integer :: nSTAT
      integer :: ncid

      !integer            :: var_xtype
      integer            :: xtype
      integer            :: length
      integer            :: attnum
      character(len=40)  :: invar
      character(len=nf90_max_name) :: name_dum
      integer :: t_dim_id
      integer :: var_id
      real(kind=sp):: dum_sp
      integer :: ivar

      if(MR_VERB.ge.1)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------                MR_Set_Met_Dims_Template_netcdf             ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      ! To set up the grid, we assume that the grid is the same for all
      ! windfiles.  There is no checking if this is actually the case.
      ! Just read the first windfile.
      iw = 1
      write(MR_global_info,*)"Opening ",iw,trim(ADJUSTL(MR_windfiles(iw)))
      nSTAT = nf90_open(trim(ADJUSTL(MR_windfiles(iw))),NF90_NOWRITE,ncid)
      call MR_NC_check_status(nSTAT,0,"nf90_open")
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: nf90_open to read header:', nf90_strerror(nSTAT)
        write(MR_global_error,*)'Could not open ',trim(ADJUSTL(MR_windfiles(iw)))
        write(MR_global_error,*)'Exiting'
        stop 1
      else
        write(MR_global_info,*)"Opened ",trim(ADJUSTL(MR_windfiles(iw)))
      endif

      ! Get dim ids, sizes, and associated dimension variable for dims:
      !  1 = time
      !  2 = pressure used for state variables
      !  3 = y or lat
      !  4 = x or lon
      !  5 = pressure used for Vz
      !  6 = pressure used for RH or SH
      !  7 = height above ground
      !  8 = depth below surface
      !  9 = extra pressure dimension

      ! Time
      !  This will be repeated in MR_Set_Met_Times_netcdf where the time values
      !  are read, but for now, we just want the dimension size
      i = 1
      if(.not.Met_dim_IsAvailable(i))then
        write(MR_global_error,*)"MR ERROR: TIME dimension is required and not listed"
        write(MR_global_error,*)"          in template windfile specification file."
        stop 1
      endif
      nSTAT = nf90_inq_dimid(ncid,Met_dim_names(i),t_dim_id)
      call MR_NC_check_status(nSTAT,1,"nf90_inq_dimid Time")
      nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,name=name_dum,len=nt_fullmet)
      call MR_NC_check_status(nSTAT,1,"nf90_Inquire_Dimension Time")

      ! Need to get fill_value.  If we can't find it, fill_value will remain as initialized
      ! according to the MR_iwindformat (probably -9999.0_sp)
      ! We'll look in variables 1,2,3 (GPH, U, V) and assume if it is found, that the
      ! value is used throughout the file
      FoundFillVAttr = .false.
      do ivar=1,3
        nSTAT = nf90_inq_varid(ncid,Met_var_NC_names(ivar),var_id)
        call MR_NC_check_status(nSTAT,1,"nf90_inq_varid")
        nSTAT = nf90_Inquire_Attribute(ncid, var_id,&
                                       "_FillValue",xtype, length, attnum)
        call MR_NC_check_status(nSTAT,0,"nf90_Inquire_Attribute")
        if(nSTAT.eq.0)then
          FoundFillVAttr = .true.
          nSTAT = nf90_get_att(ncid, var_id,"_FillValue",dum_sp)
          call MR_NC_check_status(nSTAT,1,"nf90_get_att _FillValue")
          fill_value_sp = dum_sp
          write(MR_global_info,*)"    Found fill value",fill_value_sp
          exit
        endif
      enddo

      nSTAT = nf90_close(ncid)
      call MR_NC_check_status(nSTAT,0,"nf90_close")

      end subroutine MR_Set_Met_Dims_Template_netcdf

!##############################################################################
!
!     MR_Read_MetP_Variable_netcdf
!
!     Called from Read_HGT_arrays and once from Read_3d_MetP_Variable.
!
!     Sets MR_dum3d_metP, MR_dum2d_met, or MR_dum2d_met_int as appropriate
!
!##############################################################################

      subroutine MR_Read_MetP_Variable_netcdf(ivar,istep)

      use MetReader
      use netcdf

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer,intent(in) :: ivar
      integer,intent(in) :: istep

      integer :: iw,iwstep
      integer :: np_met_loc
      character(len=130) :: infile
      character(len=71)  :: invar
      character(len=40)  :: fileposstr

      integer :: ncid       = 0
      integer :: nSTAT      = 0
      integer :: in_var_id  = 0
      integer :: in_var_id1 = 0
      integer :: in_var_id2 = 0

      real(kind=sp) :: del_H,del_P,dpdz

      integer :: i,j,k,ii,jj,kk,kkk,itmp
      integer :: ict, ileft(2),iright(2)   !if wrapgrid=.true. ict=2 and left & iright have 2 values, otherwise 1
      integer :: iistart(2),iicount(2)     !if (wrapgrid), iistart(1)=istart, iistart(2)=1

      integer :: Dimension_of_Variable
      logical :: IsCategorical

      integer :: var_xtype
      integer :: NC_version

      real(kind=sp),dimension(:,:,:,:),allocatable :: dum3d_metP_aux
      real(kind=sp) :: theta,cofac

      real(kind=sp) :: Z_top, T_top
      real(kind=sp) :: pp
      integer       :: idx

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
        write(MR_global_production,*)"----------                MR_Read_MetP_Variable_netcdf                ----------"
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      if(.not.Met_var_IsAvailable(ivar))then
        write(MR_global_error,*)"MR ERROR:  Variable not available for this windfile"
        write(MR_global_error,*)"             ivar = ",ivar
        write(MR_global_error,*)"            vname = ",Met_var_NC_names(ivar)
        write(MR_global_error,*)"             iwf  = ",MR_iwindformat
        stop 1
      endif

      iw     = MR_MetStep_findex(istep)
      iwstep = MR_MetStep_tindex(istep)

      if(Met_var_NC_names(ivar).eq."")then
        write(MR_global_info,*)"Variable ",ivar," not available for MR_iwindformat = ",&
                  MR_iwindformat
        stop 1
      endif

      ! Get the dimension of the variable requested (either 2 or 3-D)
      if(ivar.eq.1 ) Dimension_of_Variable = 3 ! Geopotential Height
      if(ivar.eq.2 ) Dimension_of_Variable = 3 ! Vx
      if(ivar.eq.3 ) Dimension_of_Variable = 3 ! Vy
      if(ivar.eq.4 ) Dimension_of_Variable = 3 ! Vz
      if(ivar.eq.5 ) Dimension_of_Variable = 3 ! Temperature
      if(ivar.eq.6 ) Dimension_of_Variable = 3 ! Pressure (only for WRF or other eta-level files)
      if(ivar.eq.7 ) Dimension_of_Variable = 3 ! PVV

      if(ivar.eq.10) Dimension_of_Variable = 2 ! Planetary Boundary Layer Height
      if(ivar.eq.11) Dimension_of_Variable = 2 ! U @ 10m
      if(ivar.eq.12) Dimension_of_Variable = 2 ! V @ 10m
      if(ivar.eq.13) Dimension_of_Variable = 2 ! Friction velocity
      if(ivar.eq.14) Dimension_of_Variable = 2 ! Displacement Height
      if(ivar.eq.15) Dimension_of_Variable = 2 ! Snow cover
      if(ivar.eq.16) Dimension_of_Variable = 2 ! Soil moisture
      if(ivar.eq.17) Dimension_of_Variable = 2 ! Surface roughness
      if(ivar.eq.18) Dimension_of_Variable = 2 ! Wind_speed_gust_surface

      if(ivar.eq.20) Dimension_of_Variable = 2 ! pressure at lower cloud base
      if(ivar.eq.21) Dimension_of_Variable = 2 ! pressure at lower cloud top
      if(ivar.eq.22) Dimension_of_Variable = 2 ! temperature at lower cloud top
      if(ivar.eq.23) Dimension_of_Variable = 2 ! Total Cloud cover
      if(ivar.eq.24) Dimension_of_Variable = 2 ! Cloud cover (low)
      if(ivar.eq.25) Dimension_of_Variable = 2 ! Cloud cover (convective)

      if(ivar.eq.30) Dimension_of_Variable = 3 ! Rel. Hum
      if(ivar.eq.31) Dimension_of_Variable = 3 ! QV (specific humidity)
      if(ivar.eq.32) Dimension_of_Variable = 3 ! QL (liquid)
      if(ivar.eq.33) Dimension_of_Variable = 3 ! QI (ice)

      if(ivar.eq.40) Dimension_of_Variable = 2 ! Categorical rain
      if(ivar.eq.41) Dimension_of_Variable = 2 ! Categorical snow
      if(ivar.eq.42) Dimension_of_Variable = 2 ! Categorical frozen rain
      if(ivar.eq.43) Dimension_of_Variable = 2 ! Categorical ice
      if(MR_iwindformat.eq.24)then
          ! NASA MERRA has 3d precipitation
        if(ivar.eq.44) Dimension_of_Variable = 3 ! Precipitation rate large-scale (liquid)
        if(ivar.eq.45) Dimension_of_Variable = 3 ! Precipitation rate convective (liquid)
      else
          ! All other met files use surface precipitation
        if(ivar.eq.44) Dimension_of_Variable = 2 ! Precipitation rate large-scale (liquid)
        if(ivar.eq.45) Dimension_of_Variable = 2 ! Precipitation rate convective (liquid)
      endif
      if(ivar.eq.46) Dimension_of_Variable = 3 ! Precipitation rate large-scale (ice)
      if(ivar.eq.47) Dimension_of_Variable = 3 ! Precipitation rate convective (ice)

      if(ivar.eq.40.or.&
         ivar.eq.41.or.&
         ivar.eq.42.or.&
         ivar.eq.43)then
          ! Categorical variables are integers and need special interpolation
        IsCategorical = .true.
      else
          ! The default is to read floating point values
        IsCategorical = .false.
      endif

      if(MR_iwind.eq.5)then
          ! Files are hard-coded
        call MR_Set_iwind5_filenames(MR_MetStep_Hour_since_baseyear(istep),ivar,infile)
        infile = trim(adjustl(infile))
      else
          ! Files are provided directly by calling program, not hard-coded
        infile = trim(adjustl(MR_MetStep_File(istep)))
      endif
      np_met_loc = nlevs_fullmet(Met_var_zdim_idx(ivar))

      invar = Met_var_NC_names(ivar)

      write(fileposstr,'(a9,i4,a9,i4,a10,i4)')"  step = ",istep,&
                         ", file = ",iw,&
                         ", slice = ",iwstep
      write(MR_global_info,*)"Reading ",trim(adjustl(invar))," from file : ",&
                trim(adjustl(infile)),fileposstr
      nSTAT = nf90_open(trim(adjustl(infile)),NF90_NOWRITE,ncid)
      call MR_NC_check_status(nSTAT,0,"nf90_open")
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR open file:',infile,nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR open file:',infile,nf90_strerror(nSTAT)
        write(MR_global_error,*)trim(adjustl(infile)),NF90_NOWRITE,ncid,nSTAT
        write(MR_global_error,*)'Exiting'
        stop 1
      endif

      ! Load data variables for just the subgrid defined above
      if (wrapgrid) then
        ict        = 2
          ! index on the sub-met
        ileft(1)   = 1;         ileft(2)   = ilhalf_nx+1
        iright(1)  = ilhalf_nx; iright(2)  = nx_submet
          ! index on the full-met
        iistart(1) = ilhalf_fm_l; iistart(2) = irhalf_fm_l
        iicount(1) = ilhalf_nx  ; iicount(2) = irhalf_nx
      else
        ict        = 1
        ileft(1)   = 1
        iright(1)  = nx_submet
        iistart(1) = istart
        iicount(1) = nx_submet
      endif

      nSTAT = nf90_inq_varid(ncid,invar,in_var_id)
      call MR_NC_check_status(nSTAT,0,"nf90_inq_varid")
      nSTAT = nf90_inquire_variable(ncid, in_var_id, invar, &
                xtype = var_xtype)
      call MR_NC_check_status(nSTAT,0,"nf90_inquire_variable")

      ! Test for version 3 vs 4 NCEP files
      if(var_xtype.eq.NF90_FLOAT)then
        NC_version = 4
      elseif(var_xtype.eq.NF90_SHORT)then
        NC_version = 3
      endif

      if(Dimension_of_Variable.eq.3)then
        MR_dum3d_metP = 0.0_sp
        if(MR_iwindformat.ne.50)then
          allocate(temp3d_sp(nx_submet,ny_submet,np_met_loc,1))
        else
            ! For MR_iwindformat = 50 (WRF), we need an extra point in p
            ! Allocate auxiliary array
          if(ivar.eq.1)then
                ! Geopotential lives on z-staggered grid
            allocate(temp3d_sp(nx_submet,ny_submet,np_met_loc+1,1))
            allocate(dum3d_metP_aux(nx_submet,ny_submet,np_met_loc+1,1))
          elseif(ivar.eq.2)then
                ! U wind lives on x-staggered grid
            allocate(temp3d_sp(nx_submet+1,ny_submet,np_met_loc,1))
            allocate(dum3d_metP_aux(nx_submet+1,ny_submet,np_met_loc,1))
          elseif(ivar.eq.3)then
                ! V wind lives on y-staggered grid
            allocate(temp3d_sp(nx_submet,ny_submet+1,np_met_loc,1))
            allocate(dum3d_metP_aux(nx_submet,ny_submet+1,np_met_loc,1))
          elseif(ivar.eq.4)then
                ! W wind lives on z-staggered grid
            allocate(temp3d_sp(nx_submet,ny_submet,np_met_loc+1,1))
            allocate(dum3d_metP_aux(nx_submet,ny_submet,np_met_loc+1,1))
          elseif(ivar.eq.5)then
                ! Temp lives on non-staggered grid, but we need pres. base and pert.
            allocate(temp3d_sp(nx_submet,ny_submet,np_met_loc,1))
            allocate(dum3d_metP_aux(nx_submet,ny_submet,np_met_loc,1))
          elseif(ivar.eq.6)then
                ! pressure lives on non-staggered grid, but we need base and pert.
            allocate(temp3d_sp(nx_submet,ny_submet,np_met_loc,1))
            allocate(dum3d_metP_aux(nx_submet,ny_submet,np_met_loc,1))
          endif
          dum3d_metP_aux(:,:,:,:)=0.0_sp
        endif ! MR_iwindformat.ne.50
        temp3d_sp(:,:,:,:)=0.0_sp

        do i=1,ict        !read subgrid at current time step
          ! Branch on four cases: (1) iw=5, NCEP with NetCDFv3
          !                       (2) All other iw=5
          !                       (3) WRF iwf=50
          !                       (4) All other iw=3/4
          if(MR_iwind.eq.5.and.MR_iwindformat.eq.25.and.var_xtype.eq.NF90_SHORT)then
            ! NCEP reanalysis files are now NCv4 (stored as float), but the
            ! older version, NCv3 (stored as short) might still be around.
            if(i.eq.1)allocate(temp3d_short(nx_submet,ny_submet,np_met_loc,1))
            nSTAT = nf90_get_var(ncid,in_var_id,temp3d_short(ileft(i):iright(i),:,:,:), &
                     start = (/iistart(i),jstart,1,iwstep/),       &
                     count = (/iicount(i),ny_submet,np_met_loc,1/))
            call MR_NC_check_status(nSTAT,0,"nf90_get_var iw5 wf25 short")
          elseif(MR_iwind.eq.5.and.(MR_iwindformat.eq.25.or.&
                                    MR_iwindformat.eq.26.or.&
                                    MR_iwindformat.eq.27.or.&
                                    MR_iwindformat.eq.29.or.&
                                    MR_iwindformat.eq.30))then
            nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                     start = (/iistart(i),jstart,1,iwstep/),       &
                     count = (/iicount(i),ny_submet,np_met_loc,1/))
            call MR_NC_check_status(nSTAT,0,"nf90_get_var iw5")
          elseif(MR_iwindformat.eq.50)then
            ! Now read the data and convert if necessary
            if(ivar.eq.1)then
                ! Geopotential
                ! Get PHB
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc+1,1/))
              call MR_NC_check_status(nSTAT,0,"nf90_get_var PHB")
              write(MR_global_info,*)istep,"Reading ","PH"," from file : ",trim(adjustl(infile))
              nSTAT = nf90_inq_varid(ncid,"PH",in_var_id2)
              call MR_NC_check_status(nSTAT,0,"nf90_inq_varid PH")
              nSTAT = nf90_get_var(ncid,in_var_id,dum3d_metP_aux(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc+1,1/))
              call MR_NC_check_status(nSTAT,0,"nf90_get_var PH")
              temp3d_sp(:,:,:,:) = temp3d_sp(:,:,:,:) + dum3d_metP_aux(:,:,:,:)
              do kk=1,np_met_loc
                MR_dum3d_metP(:,:,kk) = 0.5_sp*(temp3d_sp(:,:,kk  ,1) + &
                                          temp3d_sp(:,:,kk+1,1))
              enddo

            elseif(ivar.eq.2)then
                ! U wind lives on x-staggered grid
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i)+1,:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i)+1,ny_submet,np_met_loc,1/))
              call MR_NC_check_status(nSTAT,0,"nf90_get_var U")
              do ii=1,iicount(i)
                MR_dum3d_metP(ii,:,:) = 0.5_sp*(temp3d_sp(ii  ,:,:,1) + &
                                          temp3d_sp(ii+1,:,:,1))
              enddo
            elseif(ivar.eq.3)then
                ! V wind lives on y-staggered grid
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet+1,np_met_loc,1/))
              call MR_NC_check_status(nSTAT,0,"nf90_get_var V")
              do jj=1,ny_submet
                MR_dum3d_metP(:,jj,:) = 0.5_sp*(temp3d_sp(:,jj,:,1) + &
                                          temp3d_sp(:,jj,:,1))
              enddo
            elseif(ivar.eq.4)then
                ! W wind lives on z-staggered grid
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc+1,1/))
              call MR_NC_check_status(nSTAT,0,"nf90_get_var W")
              do kk=1,np_met_loc
                MR_dum3d_metP(:,:,kk) = 0.5_sp*(temp3d_sp(:,:,kk  ,1) + &
                                          temp3d_sp(:,:,kk+1,1))
              enddo
            elseif(ivar.eq.5)then
                ! Temperature is actually stored as potential temperature
                ! perturbation: we convert to potential temperature via
                !  theta = (pot.temp.pert + 300.0)
                ! and convert to thermodynamic temperature via
                !  Temp = theta*(p/p_0)^kappa
                !    where p_0   = 1000.0mb (or 1.0e5 Pa)
                !    and   kappa = R/c_p = 0.2854 (for dry air)
                ! Temp lives on non-staggered grid, but we need pres. base and pert.
                ! First get PB (base pressure)
              write(MR_global_info,*)istep,"Reading ","PB"," from file : ",trim(adjustl(infile))
              nSTAT = nf90_inq_varid(ncid,"PB",in_var_id1)
              call MR_NC_check_status(nSTAT,0,"nf90_inq_varid PB")
              nSTAT = nf90_get_var(ncid,in_var_id1,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc,1/))
              call MR_NC_check_status(nSTAT,0,"nf90_get_var PB")
                ! Now get P (perturbation pressure)
              write(MR_global_info,*)istep,"Reading ","P"," from file : ",trim(adjustl(infile))
              nSTAT = nf90_inq_varid(ncid,"P",in_var_id2)
              call MR_NC_check_status(nSTAT,0,"nf90_inq_varid P")
              nSTAT = nf90_get_var(ncid,in_var_id,dum3d_metP_aux(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc,1/))
              call MR_NC_check_status(nSTAT,0,"nf90_get_var P")
                ! Now get total pressure in Pa
              MR_dum3d_metP(:,:,:) = temp3d_sp(:,:,:,1) + dum3d_metP_aux(:,:,:,1)

                ! Read potential temperature perturbation
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc,1/))
              call MR_NC_check_status(nSTAT,0,"nf90_get_var T")
              do ii=1,iicount(i)
                do jj=1,ny_submet
                  do kk=1,np_met_loc
                    theta = temp3d_sp(ii,jj,kk,1) + 300.0_sp
                    cofac = (MR_dum3d_metP(ii,jj,kk)*1.0e-5_sp)**(0.2854_sp)
                    MR_dum3d_metP(ii,jj,kk) = theta * cofac
                  enddo
                enddo
              enddo

            elseif(ivar.eq.6)then
                ! pressure
                ! Get PB (base pressure)
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc,1/))
              call MR_NC_check_status(nSTAT,0,"nf90_get_var PB")
                ! Get P (perturbation pressure)
              nSTAT = nf90_inq_varid(ncid,"P",in_var_id2)
              call MR_NC_check_status(nSTAT,0,"nf90_inq_varid P")
              nSTAT = nf90_get_var(ncid,in_var_id,dum3d_metP_aux(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc,1/))
              call MR_NC_check_status(nSTAT,0,"nf90_get_var P")
                ! Now get total pressure in Pa
              temp3d_sp(:,:,:,:) = temp3d_sp(:,:,:,:) + dum3d_metP_aux(:,:,:,:)
              MR_dum3d_metP(:,:,:) = temp3d_sp(:,:,:,1)
            else
            ! for any other 3d WRF variable, assume non-staggered grid

            endif

          else ! end of MR_iwind=5 and iwf=50 (WRF) sections

            ! for any other 3d variable (non-WRF, non-NCEP/2.5 reanalysis)
            nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                     start = (/iistart(i),jstart,1,iwstep/),       &
                     count = (/iicount(i),ny_submet,np_met_loc,1/))
            call MR_NC_check_status(nSTAT,0,"nf90_get_var")
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: get_var: ',nf90_strerror(nSTAT)
              write(MR_global_error,*)i
              write(MR_global_error,*)ileft(i),iright(i)
              write(MR_global_error,*)iistart(i),jstart,1,iwstep
              write(MR_global_error,*)iicount(i),ny_submet,np_met_loc,1
              write(MR_global_log  ,*)'MR ERROR: get_var: ',nf90_strerror(nSTAT)
              stop 1
            endif
          endif
        enddo
        !  At this point, we have temp3d_sp with the direct read from the file
        !  Next, we need to copy it to MR_dum3d_metP

        if(MR_iwindformat.ne.50)then
          do j=1,ny_submet
            itmp = ny_submet-j+1
            !reverse the j indices (since they increment from N to S)
            if(MR_iwindformat.eq.25)then
              if(var_xtype.eq.NF90_FLOAT)then
                  ! No scaling/offset needed
                MR_dum3d_metP(1:nx_submet,j,1:np_met_loc) = &
                                      real(temp3d_sp(:,itmp,:,1),kind=sp)
              elseif(var_xtype.eq.NF90_SHORT)then
                  ! NC_version = 3 :: need to scale and offset shorts to get float
                MR_dum3d_metP(1:nx_submet,j,1:np_met_loc) = &
                                      real(temp3d_short(1:nx_submet,itmp,1:np_met_loc,1),kind=sp) * &
                                    iwf25_scale_facs(ivar) + iwf25_offsets(ivar)
              endif
            else
              if(y_inverted)then
                MR_dum3d_metP(1:nx_submet,j,1:np_met_loc)  = temp3d_sp(1:nx_submet,itmp,1:np_met_loc,1)
              else
                MR_dum3d_metP(1:nx_submet,j,1:np_met_loc)  = temp3d_sp(1:nx_submet,j,1:np_met_loc,1)
              endif
            endif
          enddo
          if(z_inverted)then ! reverse the vertical coordinate
            temp3d_sp(1:nx_submet,1:ny_submet,1:np_met_loc,1) = MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_met_loc)
            do i=1,np_met_loc
              itmp = np_met_loc-i+1
              MR_dum3d_metP(1:nx_submet,1:ny_submet,i) = temp3d_sp(1:nx_submet,1:ny_submet,itmp,1)
            enddo
          endif
        endif !MR_iwindformat.eq.50, MR_iwindformat.eq.25, else

        if(MR_iwind.eq.5.and.MR_iwindformat.eq.25)then
          if(allocated(temp3d_short)) deallocate(temp3d_short)
        endif
        if(MR_iwindformat.eq.50)then
          deallocate(dum3d_metP_aux)
        endif
        deallocate(temp3d_sp)

      elseif(Dimension_of_Variable.eq.2)then
        if(IsCategorical)then
          allocate(temp2d_int(nx_submet,ny_submet,1));temp2d_int(:,:,:)=0
          do i=1,ict        !read subgrid at current time step
            if(MR_iwindformat.eq.25)then
              ! No categorical variables for MR_iwindformat = 25
            else
              nSTAT = nf90_get_var(ncid,in_var_id,temp2d_int(ileft(i):iright(i),:,:), &
                         start = (/iistart(i),jstart,iwstep/),       &
                         count = (/iicount(i),ny_submet,1/))
              call MR_NC_check_status(nSTAT,0,"nf90_get_var categorical")
              do j=1,ny_submet
                itmp = ny_submet-j+1
                if(y_inverted)then
                  MR_dum2d_met_int(1:nx_submet,j)  = temp2d_int(1:nx_submet,itmp,1)
                else
                  MR_dum2d_met_int(1:nx_submet,j)  = temp2d_int(1:nx_submet,j,1)
                endif
              enddo
            endif
          enddo
          deallocate(temp2d_int)
        else ! non-categorical 2-d variable
          allocate(temp2d_sp(nx_submet,ny_submet,1));temp2d_sp(:,:,:)=0.0_sp
          if(ivar.eq.11.or.ivar.eq.12)then
              ! Surface winds usually have a z coordinate as well
            allocate(temp3d_sp(nx_submet,ny_submet,1,1));temp3d_sp(:,:,:,:)=0.0_sp
          endif
  
          do i=1,ict        !read subgrid at current time step
            if(MR_iwindformat.eq.25)then
              allocate(tmpsurf2d_short(192,94,1))
              if(NC_version.eq.4)then
                nSTAT = nf90_get_var(ncid,in_var_id,temp2d_sp(:,:,1), &
                         start = (/1,1,iwstep/),       &
                         count = (/192,94,1/))
                call MR_NC_check_status(nSTAT,0,"nf90_get_var NCEP 2d non-cat ncv4")
                call MR_interp_iwf25_grid(nx_submet,ny_submet,tmpsurf2d_short,temp2d_sp,&
                                    iwf25_scale_facs(ivar),iwf25_offsets(ivar))
                MR_dum2d_met(1:nx_submet,:) = temp2d_sp(1:nx_submet,:,1)
              else
                nSTAT = nf90_get_var(ncid,in_var_id,tmpsurf2d_short(:,:,1), &
                         start = (/1,1,iwstep/),       &
                         count = (/192,94,1/))
                call MR_NC_check_status(nSTAT,0,"nf90_get_var NCEP 2d non-cat ncv3")
                call MR_interp_iwf25_grid(nx_submet,ny_submet,tmpsurf2d_short,temp2d_sp,&
                                    iwf25_scale_facs(ivar),iwf25_offsets(ivar))
                MR_dum2d_met(1:nx_submet,:) = temp2d_sp(1:nx_submet,:,1)
              endif
            else  ! Not MR_iwindformat.eq.25
              ! 2d variables for iwf .ne. 25
              if(ivar.eq.11.or.ivar.eq.12)then
                ! Surface velocities do have a z dimension
                nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                         start = (/iistart(i),jstart,1,iwstep/),       &
                         count = (/iicount(i),ny_submet,1,1/))
                call MR_NC_check_status(nSTAT,0,"nf90_get_var surface_vel")
                do j=1,ny_submet
                  itmp = ny_submet-j+1
                  if(y_inverted)then
                    MR_dum2d_met(1:nx_submet,j)  = temp3d_sp(1:nx_submet,itmp,1,1)
                  else
                    MR_dum2d_met(1:nx_submet,j)  = temp3d_sp(1:nx_submet,j,1,1)
                  endif
                enddo
              else  ! not MR_iwindformat.eq.25 and not ivar.eq.11.or.ivar.eq.12
                nSTAT = nf90_get_var(ncid,in_var_id,temp2d_sp(ileft(i):iright(i),:,:), &
                         start = (/iistart(i),jstart,iwstep/),       &
                         count = (/iicount(i),ny_submet,1/))
                call MR_NC_check_status(nSTAT,0,"nf90_get_var non-surface_vel")
                do j=1,ny_submet
                  itmp = ny_submet-j+1
                  if(y_inverted)then
                    MR_dum2d_met(1:nx_submet,j)  = temp2d_sp(1:nx_submet,itmp,1)
                  else
                    MR_dum2d_met(1:nx_submet,j)  = temp2d_sp(1:nx_submet,j,1)
                  endif
                enddo
              endif
            endif
          enddo
          deallocate(temp2d_sp)
          if(ivar.eq.11.or.ivar.eq.12) deallocate(temp3d_sp)
          if(MR_iwindformat.eq.25) deallocate(tmpsurf2d_short)
        endif ! IsCategorical
      endif ! Dimension_of_Variable.eq.2

      ! Quality control checks
      if(ivar.eq.1)then
        ! If this is filling HGT, then we need to do a special QC check
        do i=1,nx_submet
          do j=1,ny_submet
            do k=1,np_met_loc
              if(abs(MR_dum3d_metP(i,j,k)-fill_value_sp).lt.MR_EPS_SMALL.or.&
                     MR_dum3d_metP(i,j,k).lt.0.0_sp)then  ! also flag value as to be reset if it
                                                          ! maps below sea level
                 ! linearly interpolate in z
                 ! find the first non NaN above k
                 do kk = k+1,np_met_loc,1
                   if(abs(MR_dum3d_metP(i,j,kk)-fill_value_sp).gt.MR_EPS_SMALL.and.&
                          MR_dum3d_metP(i,j,kk).ge.0.0_sp)exit
                 enddo
                 if(kk.eq.np_met_loc+1)then
                   kk=np_met_loc
                   MR_dum3d_metP(i,j,kk) = 0.0_sp
                 endif
                 ! find the first non NaN below k if k!=1
                 do kkk = max(k-1,1),1,-1
                   if(abs(MR_dum3d_metP(i,j,kkk)-fill_value_sp).gt.MR_EPS_SMALL.and.&
                          MR_dum3d_metP(i,j,kkk).ge.0.0_sp)exit
                 enddo
                 if(kkk.eq.0)then
                   kkk=1
                   MR_dum3d_metP(i,j,kkk) = 0.0_sp
                 endif
                 MR_dum3d_metP(i,j,k) = MR_dum3d_metP(i,j,kkk) + &
                       (MR_dum3d_metP(i,j,kk)-MR_dum3d_metP(i,j,kkk)) * &
                       real(k-kkk,kind=sp)/real(kk-kkk,kind=sp)
              endif
            enddo
          enddo
        enddo
        ! convert m to km
        MR_dum3d_metP = MR_dum3d_metP / 1000.0_sp
      elseif(Dimension_of_Variable.eq.3)then
        ! Do QC checking of all other 3d variables
        if(ivar.eq.2.or.ivar.eq.3.or.ivar.eq.4.or.ivar.eq.7)then
          ! taper winds (vx,vy,vz) to zero at ground surface
          if(istep.eq.MR_iMetStep_Now)then
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp, &
                          bc_low_sp=0.0_sp)
          else
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_next,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp, &
                          bc_low_sp=0.0_sp)
          endif
        elseif(ivar.eq.5)then
          ! set ground and top-level conditions for temperature
          idx = Met_var_zdim_idx(ivar)
          pp = levs_fullmet_sp(idx,nlevs_fullmet(idx))/real(100.0,kind=sp)
          Z_top = MR_Z_US_StdAtm(pp)
          T_top = MR_Temp_US_StdAtm(Z_top)
          if(istep.eq.MR_iMetStep_Now)then
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp, &
                          bc_low_sp=293.0_sp, bc_high_sp=T_top)
          else
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_next,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp, &
                          bc_low_sp=293.0_sp, bc_high_sp=T_top)
          endif
        else
          ! For other variables, use the top and bottom non-fill values
          if(istep.eq.MR_iMetStep_Now)then
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp)
          else
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_next,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp)
          endif
        endif
      endif

      if(ivar.eq.4)then
        if(MR_iwindformat.ne.50)then
            ! For vertical velocity, we convert from pressure vertical velocity
            ! Pa s to m/s by dividing by pressure gradient.
            !  Note:  This seems problematic in the stratosphere where pressure
            !         levels are further apart.  It might be better to read PVV
            !         directly, and use dp/dz = -rho g where we get rho from the
            !         ideal gas law
          idx = Met_var_zdim_idx(ivar)
          do k=1,np_met_loc
            do i=1,nx_submet
              do j=1,ny_submet
                if(k.eq.1)then
                  ! Use one-sided gradients for bottom
                  del_p = levs_fullmet_sp(idx,2) - levs_fullmet_sp(idx,1)
                  !del_P = p_fullmet_Vz_sp(2)-p_fullmet_Vz_sp(1)
                  if(istep.eq.MR_iMetStep_Now)then
                    del_H = MR_geoH_metP_last(i,j,2) - MR_geoH_metP_last(i,j,1)
                  else
                    del_H = MR_geoH_metP_next(i,j,2) - MR_geoH_metP_next(i,j,1)
                  endif
                elseif(k.eq.np_met_loc)then
                  ! Use one-sided gradients for top
                  !del_P = p_fullmet_Vz_sp(np_met_loc) - &
                  !         p_fullmet_Vz_sp(np_met_loc-1)
                  del_p = levs_fullmet_sp(idx,nlevs_fullmet(idx)) - &
                          levs_fullmet_sp(idx,nlevs_fullmet(idx)-1)
                  if(istep.eq.MR_iMetStep_Now)then
                    del_H = MR_geoH_metP_last(i,j,np_met_loc) - &
                             MR_geoH_metP_last(i,j,np_met_loc-1)
                  else
                    del_H = MR_geoH_metP_next(i,j,np_met_loc) - &
                             MR_geoH_metP_next(i,j,np_met_loc-1)
                  endif
                else
                  ! otherwise, two-sided calculation
                  !del_P = p_fullmet_Vz_sp(k+1)-p_fullmet_Vz_sp(k-1)
                  del_p = levs_fullmet_sp(idx,k+1) - &
                          levs_fullmet_sp(idx,k-1)
                  if(istep.eq.MR_iMetStep_Now)then
                    del_H = MR_geoH_metP_last(i,j,k+1) - MR_geoH_metP_last(i,j,k-1)
                  else
                    del_H = MR_geoH_metP_next(i,j,k+1) - MR_geoH_metP_next(i,j,k-1)
                  endif
                endif
                del_H = del_H * 1000.0_sp ! convert to m
                if(abs(del_H).gt.MR_EPS_SMALL)then
                  dpdz  = del_P/del_H
                else
                  write(MR_global_error,*)'MR ERROR: failed to calculate dpdz'
                  write(MR_global_error,*)i,j,k,del_P,del_H
                  write(MR_global_error,*)MR_geoH_metP_last(i,j,:)
                  stop 1
                endif
                MR_dum3d_metP(i,j,k) = MR_dum3d_metP(i,j,k) / dpdz
              enddo ! j
            enddo ! i
          enddo ! k
        endif
      endif

      ! Close file
      nSTAT = nf90_close(ncid)
      call MR_NC_check_status(nSTAT,0,"nf90_close")

      MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_met_loc) =  &
      MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_met_loc) * Met_var_conversion_factor(ivar)

      if(MR_VERB.ge.2)then
        write(MR_global_production,*)"--------------------------------------------------------------------------------"
      endif

      end subroutine MR_Read_MetP_Variable_netcdf

!##############################################################################
!
!    MR_interp_iwf25_grid
!
!##############################################################################

      subroutine MR_interp_iwf25_grid(imax,jmax,invar,outvar,scale_fac,offset)

      use MetReader

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer         ,intent(in)  :: imax,jmax
      integer(kind=sp),intent(in)  :: invar(192,94,1)
      real(kind=sp)   ,intent(out) :: outvar(imax,jmax)
      real(kind=sp)   ,intent(in)  :: scale_fac
      real(kind=sp)   ,intent(in)  :: offset

      real(kind=sp)    :: a1,a2,a3,a4
      real(kind=sp)    :: v1,v2,v3,v4

      integer :: ilon,ilat

      do ilon = 1,imax
        do ilat = 1,jmax
          a1 = amap_iwf25(ilon,ilat,1)
          a2 = amap_iwf25(ilon,ilat,2)
          a3 = amap_iwf25(ilon,ilat,3)
          a4 = amap_iwf25(ilon,ilat,4)
          v1 = invar(imap_iwf25(ilon,ilat,1),imap_iwf25(ilon,ilat,3),1) &
                 * scale_fac + offset
          v2 = invar(imap_iwf25(ilon,ilat,2),imap_iwf25(ilon,ilat,3),1) &
                 * scale_fac + offset
          v3 = invar(imap_iwf25(ilon,ilat,2),imap_iwf25(ilon,ilat,4),1) &
                 * scale_fac + offset
          v4 = invar(imap_iwf25(ilon,ilat,1),imap_iwf25(ilon,ilat,4),1) &
                 * scale_fac + offset

          outvar(ilon,ilat) = a1*v1 + a2*v2 + a3*v3 + a4*v4
        enddo
      enddo

      end subroutine MR_interp_iwf25_grid

!##############################################################################
!
!    MR_NC_check_status
!
!    nSTAT   = error code returned from netcdf call
!    errcode = user-supplied return value on stopping of code
!    operation = string descriptor of function call causing error
!
!    Error-checking routine for NetCDF function calls.
!    Modeled after a subroutine posted at:
!    https://climate-cms.org/2018/10/12/create-netcdf.html
!
!##############################################################################

      subroutine MR_NC_check_status(nSTAT, errcode, operation)

      use MetReader
      use netcdf

      implicit none

      integer, intent(in) :: nSTAT
      integer, intent(in) :: errcode
      character(len=*), intent(in) :: operation

      character(len=9) :: severity

      if (errcode.eq.0)then
        severity = "WARNING: "
       else
        severity = "ERROR:   "
      endif

      if (nSTAT == nf90_noerr) return
      !write(MR_global_essential ,*)severity,errcode,operation,nf90_strerror(nSTAT)
      write(MR_global_log ,*)severity,errcode,operation,nf90_strerror(nSTAT)
      write(MR_global_error ,*)severity,errcode,operation,nf90_strerror(nSTAT)

      ! If user-supplied error code is 0, then consider this a warning,
      ! otherwise do a hard stop
      if (errcode.ne.0) stop

      end subroutine MR_NC_check_status


