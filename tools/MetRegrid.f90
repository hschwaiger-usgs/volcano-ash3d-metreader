!##############################################################################
!##############################################################################
!
! MetRegrid
!
! This is a stand-alone program that uses the MetReader interface to calculate
! regridded atmospheric data.  A control file is provided which defines the
! computational grid specifications, output time, output variable, output
! grid (could be CompH, CompP, MetH, or MetP), the output dimension (0,1, 2, or 3).
! If 0, then the x,y,z indexes is needed (point output)
! If 1, then 2 coordinate indexes are needed (line output)
!       e.g. a z profile at x=1,y=2 would be : 1 2 0
!            a y profile at x=5,z=1 would be : 5 0 1
! If 2, then 1 index along off-axis is needed
!       e.g. a xy grid at the surface would be : 0 0 1
!            a yz grid at x=5 would be         : 5 0 0
! If 3, then the full grid is written
!
!##############################################################################

      program MetRegrid

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_error,verbosity_info,&
         MR_Comp_StartHour,MR_Comp_Time_in_hours,MR_iMetStep_Now, &
         x_submet_sp,y_submet_sp,p_fullmet_sp,&
         MR_MetStep_Interval,Met_var_GRIB_names,Met_var_NC_names,&
         Met_var_IsAvailable,&
         MR_dum3d_compH,MR_dum3d_compP,MR_dum3d_metH,MR_dum3d_metP,&
         np_fullmet,nx_submet,ny_submet,&
           MR_Read_Met_DimVars,&
           MR_Set_Met_Times,&
           MR_Initialize_Met_Grids,&
           MR_Read_3d_MetP_Variable,&
           MR_Read_3d_Met_Variable_to_CompH,&
           MR_Read_3d_Met_Variable_to_CompP,&
           MR_Read_3d_MetH_Variable,&
           MR_Read_HGT_arrays,&
           MR_Allocate_FullMetFileList,&
           MR_MetStep_Hour_since_baseyear,&
           MR_Reset_Memory,&
           MR_FileIO_Error_Handler

      implicit none

      integer,parameter :: fid_outfile = 20

      ! These are the variables that must be set in the input file or command line
      logical             :: IsPeriodic
      logical             :: IsLatLon
      real(kind=4)        :: xLL,xUR
      real(kind=4)        :: yLL,yUR
      real(kind=4)        :: zbot,ztop
      real(kind=4)        :: gridwidth_x
      real(kind=4)        :: gridwidth_y
      real(kind=4)        :: gridwidth_z
      real(kind=4)        :: dx
      real(kind=4)        :: dy
      real(kind=4)        :: dz_const
      integer             :: inyear,inmonth,inday
      real(kind=8)        :: inhour
      real(kind=8)        :: tfrac
      integer             :: outgrid_ID
      integer             :: outvar_ID
      integer             :: outdim_ID
      integer             :: outformat       ! output format
      logical             :: IsH
      integer             :: iidx
      integer             :: jidx
      integer             :: kidx

      integer             :: nxmax,nymax,nzmax
      real(kind=4),dimension(:)    ,allocatable :: x_cc
      real(kind=4),dimension(:)    ,allocatable :: y_cc
      real(kind=4),dimension(:)    ,allocatable :: z_cc
      integer             :: i,j,k
      integer             :: nx,ny,nz
      character(len=13)   :: outfilename
      real(kind=4)        :: zout
      character(len=20)   :: filename_root

      real(kind=4),dimension(:,:,:),allocatable :: out3d_t1,out3d_t2,out3d
      
      integer :: io                           ! Index for output streams

      INTERFACE
      subroutine Read_ComdLine_InpFile(IsLatLon,xLL,yLL,zbot,&
                      gridwidth_x,gridwidth_y,gridwidth_z,&
                      dx,dy,dz_const,&
                      inyear,inmonth,inday,inhour,&
                      outgrid_ID,outvar_ID,outdim_ID,outformat,iidx,jidx,kidx)
           logical                  ,intent(out) :: IsLatLon
           real(kind=4)             ,intent(out) :: xLL
           real(kind=4)             ,intent(out) :: yLL
           real(kind=4)             ,intent(out) :: zbot
           real(kind=4)             ,intent(out) :: gridwidth_x
           real(kind=4)             ,intent(out) :: gridwidth_y
           real(kind=4)             ,intent(out) :: gridwidth_z
           real(kind=4)             ,intent(out) :: dx
           real(kind=4)             ,intent(out) :: dy
           real(kind=4)             ,intent(out) :: dz_const
           integer                  ,intent(out) :: inyear,inmonth,inday
           real(kind=8)             ,intent(out) :: inhour
           integer                  ,intent(out) :: outgrid_ID
           integer                  ,intent(out) :: outvar_ID
           integer                  ,intent(out) :: outdim_ID
           integer                  ,intent(out) :: outformat
           integer                  ,intent(out) :: iidx
           integer                  ,intent(out) :: jidx
           integer                  ,intent(out) :: kidx
        end subroutine Read_ComdLine_InpFile
        subroutine write_2D_ASCII(IsLatLon,nx,ny,xLL,yLL,dx,dy,OutVar,filename_root)
           logical                  ,intent(in) :: IsLatLon
           integer                  ,intent(in) :: nx,ny
           real(kind=4)             ,intent(in) :: xLL,yLL
           real(kind=4)             ,intent(in) :: dx,dy
           real(kind=4)             ,intent(in) :: OutVar(nx,ny)
           character(len=20)        ,intent(in) :: filename_root
        end subroutine
      END INTERFACE

      call Read_ComdLine_InpFile(IsLatLon,xLL,yLL,zbot,&
                      gridwidth_x,gridwidth_y,gridwidth_z,&
                      dx,dy,dz_const,&
                      inyear,inmonth,inday,inhour,&
                      outgrid_ID,outvar_ID,outdim_ID,outformat,iidx,jidx,kidx)

        ! Check for existance and compatibility with simulation time requirements
      call MR_Read_Met_DimVars(inyear)

      call MR_Set_Met_Times(MR_Comp_StartHour,MR_Comp_Time_in_hours)

      ! Now set up the computational grid
      IsPeriodic  = .false.
      if(IsLatLon)then
        if(xLL.lt.-360.0_4)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)"MR ERROR: xLL must be >-360 and <360"
            write(errlog(io),*)"          From control file, xLL = ",xLL
          endif;enddo
          stop 1
        elseif(xLL.lt.0.0_4)then
          ! shift negative lon to range 0->360
          xLL=xLL+360.0_4
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"Resetting to range 0->360"
          endif;enddo
        elseif(xLL.gt.360.0_4)then 
          do io=1,MR_nio;if(VB(io).le.verbosity_info)then
            write(outlog(io),*)"WARNING: xLL is greater than 360"
            write(outlog(io),*)"Resetting to range 0->360"
          endif;enddo
          xLL=xLL-360.0_4
        endif
        if(gridwidth_x.ge.360.0_4.or.&
           abs(gridwidth_x-360.0_4).lt.1.0e-5)then
          IsPeriodic  = .true.
          xLL         = 0.0_4
          gridwidth_x = 360.0_4
        endif
      endif

      xUR  = xLL  + gridwidth_x
      yUR  = yLL  + gridwidth_y
      ztop = zbot + gridwidth_z
      nxmax = ceiling((xUR-xLL)/dx)         !number of x nodes
      nymax = ceiling((yUR-yLL)/dy)         !number of y nodes
      nzmax = ceiling((ztop-zbot)/dz_const)
      allocate(x_cc(1:nxmax))
      allocate(y_cc(1:nymax))
      allocate(z_cc(1:nzmax))
      do i=1,nxmax
        x_cc(i) = xLL  + dx*real(i,kind=4) - dx*0.5_4
      enddo
      do j=1,nymax
        y_cc(j) = yLL  + dy*real(j,kind=4) - dy*0.5_4
      enddo
      do k=1,nzmax
        z_cc(k) = zbot + dz_const*real(k,kind=4) - dz_const*0.5_4
      enddo

      call MR_Initialize_Met_Grids(nxmax,nymax,nzmax,                    &
                              x_cc(1:nxmax),y_cc(1:nymax),z_cc(1:nzmax), &
                              IsPeriodic)

      ! Noting what we are planning to do to the stdout
      if(Met_var_IsAvailable(outvar_ID))then
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Planning to read variable:",&
                    trim(adjustl(Met_var_NC_names(outvar_ID)))," (",&
                    trim(adjustl(Met_var_GRIB_names(outvar_ID))),")"
        endif;enddo
      else
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"Variable is not available for this windfile"
          write(errlog(io),*)"  outvar_ID = ",outvar_ID
          write(errlog(io),*)"  var name  = ",trim(adjustl(Met_var_NC_names(outvar_ID)))," (",&
                                     trim(adjustl(Met_var_GRIB_names(outvar_ID))),")"
        endif;enddo
        stop 1
      endif

      ! Output grid (1,2,3,4 for CompH,CompP,MetH,MetP)
      if(outgrid_ID.eq.1)then
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Output will be interpolated onto the CompH grid"
        endif;enddo
        IsH = .true.
      elseif(outgrid_ID.eq.2)then
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Output will be interpolated onto the CompP grid"
        endif;enddo
        IsH = .false.
      elseif(outgrid_ID.eq.3)then
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Output will be interpolated onto the MetH grid"
        endif;enddo
        IsH = .true.
      elseif(outgrid_ID.eq.4)then
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Output will be interpolated onto the MetP grid"
        endif;enddo
        IsH = .false.
      else
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: outgrid_ID not recognized; expecting 1-4"
          write(errlog(io),*)"          outgrid_ID = ",outgrid_ID
        endif;enddo
        stop 1
      endif

      ! Now actually interpolate in time
      MR_iMetStep_Now = 1
      tfrac = (MR_Comp_StartHour-MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now))/ &
               MR_MetStep_Interval(MR_iMetStep_Now)
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Hour = ",inhour
        write(outlog(io),*)"Using a fractional step of ",tfrac
      endif;enddo
      if(outgrid_ID.eq.1)then
        ! We need CompH grids
        nx = nxmax
        ny = nymax
        nz = nzmax
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)" Output array is of size ",nx,ny,nz
        endif;enddo
        allocate(out3d_t1(nx,ny,nz))
        allocate(out3d_t2(nx,ny,nz))
        allocate(   out3d(nx,ny,nz))
        call MR_Read_HGT_arrays(MR_iMetStep_Now)
        call MR_Read_3d_Met_Variable_to_CompH(outvar_ID,MR_iMetStep_Now)
        out3d_t1(1:nx,1:ny,1:nz) = MR_dum3d_compH(1:nx,1:ny,1:nz)
        call MR_Read_3d_Met_Variable_to_CompH(outvar_ID,MR_iMetStep_Now+1)
        out3d_t2(1:nx,1:ny,1:nz) = MR_dum3d_compH(1:nx,1:ny,1:nz)
        out3d = out3d_t1 + (out3d_t2-out3d_t1)*real(tfrac,kind=4)
      elseif(outgrid_ID.eq.2)then
        ! We need CompP grids
        nx = nxmax
        ny = nymax
        nz = np_fullmet
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)" Output array is of size ",nx,ny,nz
        endif;enddo
        allocate(out3d_t1(nx,ny,nz))
        allocate(out3d_t2(nx,ny,nz))
        allocate(   out3d(nx,ny,nz))
        call MR_Read_3d_Met_Variable_to_CompP(outvar_ID,MR_iMetStep_Now)
        out3d_t1(1:nx,1:ny,1:nz) = MR_dum3d_compP(1:nx,1:ny,1:nz)
        call MR_Read_3d_Met_Variable_to_CompP(outvar_ID,MR_iMetStep_Now+1)
        out3d_t2(1:nx,1:ny,1:nz) = MR_dum3d_compP(1:nx,1:ny,1:nz)
        out3d = out3d_t1 + (out3d_t2-out3d_t1)*real(tfrac,kind=4)
      elseif(outgrid_ID.eq.3)then
        ! We need MetH grids
        nx = nx_submet
        ny = ny_submet
        nz = nzmax
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)" Output array is of size ",nx,ny,nz
        endif;enddo
        allocate(out3d_t1(nx,ny,nz))
        allocate(out3d_t2(nx,ny,nz))
        allocate(   out3d(nx,ny,nz))
        call MR_Read_HGT_arrays(MR_iMetStep_Now)
        call MR_Read_3d_MetH_Variable(outvar_ID,MR_iMetStep_Now)
        out3d_t1(1:nx,1:ny,1:nz) = MR_dum3d_metH(1:nx,1:ny,1:nz)
        call MR_Read_3d_MetH_Variable(outvar_ID,MR_iMetStep_Now+1)
        out3d_t2(1:nx,1:ny,1:nz) = MR_dum3d_metH(1:nx,1:ny,1:nz)
        out3d = out3d_t1 + (out3d_t2-out3d_t1)*real(tfrac,kind=4)
      elseif(outgrid_ID.eq.4)then
        ! We need MetP grids
        nx = nx_submet
        ny = ny_submet
        nz = np_fullmet
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)" Output array is of size ",nx,ny,nz
        endif;enddo
        allocate(out3d_t1(nx,ny,nz))
        allocate(out3d_t2(nx,ny,nz))
        allocate(   out3d(nx,ny,nz))
        call MR_Read_3d_MetP_Variable(outvar_ID,MR_iMetStep_Now)
        out3d_t1(1:nx,1:ny,1:nz) = MR_dum3d_metP(1:nx,1:ny,1:nz)
        call MR_Read_3d_MetP_Variable(outvar_ID,MR_iMetStep_Now+1)
        out3d_t2(1:nx,1:ny,1:nz) = MR_dum3d_metP(1:nx,1:ny,1:nz)
        out3d = out3d_t1 + (out3d_t2-out3d_t1)*real(tfrac,kind=4)
      else
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: outgrid_ID not recognized; expecting 1-4"
          write(errlog(io),*)"          outgrid_ID = ",outgrid_ID
        endif;enddo
        stop 1
      endif

      ! What are we writing out
200   format(a3,i0.2,a1,i1,a1,i1,a4)

      write(outfilename,200) 'out',outvar_ID,'_',outgrid_ID,'_',outdim_ID,'.dat'
      open(unit=fid_outfile,file=outfilename)
      ! First check the indexes for validity
      if(iidx.lt.0.or.jidx.lt.0.or.kidx.lt.0)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: Output index must be non-negative"
          write(errlog(io),*)"   iidx = ", iidx
          write(errlog(io),*)"   jidx = ", jidx 
          write(errlog(io),*)"   kidx = ", kidx
        endif;enddo
        stop 1
      elseif(iidx.gt.nx)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: Output index must not be greater than output array."
          write(errlog(io),*)"   nx   = ", nx
          write(errlog(io),*)"   iidx = ", iidx
        endif;enddo
        stop 1
      elseif(jidx.gt.ny)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: Output index must not be greater than output array."
          write(errlog(io),*)"   ny   = ", ny
          write(errlog(io),*)"   jidx = ", kidx
        endif;enddo
        stop 1
      elseif(kidx.gt.nz)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: Output index must not be greater than output array."
          write(errlog(io),*)"   nz   = ", nz
          write(errlog(io),*)"   kidx = ", kidx
        endif;enddo
        stop 1
      endif
210   format(3i8,4f12.3)
      if(outdim_ID.eq.0)then
        ! We need all three indexes to be non-zero for point output
        if(iidx.eq.0.or.jidx.eq.0.or.kidx.eq.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)" Unknown point",iidx,jidx,kidx
            write(errlog(io),*)" For 0-d output, each index must be specified"
          endif;enddo
          stop 1
        endif
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Writing point data at:"
          write(outlog(io),*)"    ", iidx
          write(outlog(io),*)"    ", jidx 
          write(outlog(io),*)"    ", kidx
        endif;enddo
        if(outgrid_ID.eq.1)then
          ! CompH
          write(fid_outfile,210)iidx,jidx,kidx,x_cc(iidx),y_cc(jidx),z_cc(kidx),out3d(iidx,jidx,kidx)
        elseif(outgrid_ID.eq.1)then
          ! CompP
          write(fid_outfile,210)iidx,jidx,kidx,x_cc(iidx),y_cc(jidx),p_fullmet_sp(kidx),out3d(iidx,jidx,kidx)
        elseif(outgrid_ID.eq.3)then
          ! MetH
          write(fid_outfile,210)iidx,jidx,kidx,x_submet_sp(iidx),y_submet_sp(jidx),z_cc(kidx),out3d(iidx,jidx,kidx)
        elseif(outgrid_ID.eq.4)then
          ! MetP
          write(fid_outfile,210)iidx,jidx,kidx,x_submet_sp(iidx),y_submet_sp(jidx),p_fullmet_sp(kidx),out3d(iidx,jidx,kidx)
        endif
      elseif(outdim_ID.eq.1)then
        ! We need two of three indexes to be non-zero for line output
        if((iidx.eq.0.and.(jidx.eq.0.or.kidx.eq.0)).or.&
           (jidx.eq.0.and.kidx.eq.0))then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)" Unknown line",iidx,jidx,kidx
            write(errlog(io),*)" For 1-d output, two of three indexes must be specified"
          endif;enddo
          stop 1
        endif
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Writing line data at:"
        endif;enddo
        if(iidx.eq.0)then
          ! line along x at jinx,kind
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.2)then
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    y = ",y_cc(jidx),jidx
              write(outlog(io),*)"    x range = ",x_cc(1),x_cc(nxmax)
            endif;enddo
          else
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    y = ",y_submet_sp(jidx),jidx
              write(outlog(io),*)"    x range = ",y_submet_sp(1),y_submet_sp(nymax)
            endif;enddo
          endif
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.3)then
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    z = ",z_cc(kidx),kidx
            endif;enddo
          else
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    p = ",p_fullmet_sp(kidx),kidx
            endif;enddo
          endif

          do i=1,nx
            if(outgrid_ID.eq.1)then
              ! CompH
              write(fid_outfile,210)i,jidx,kidx,x_cc(i),y_cc(jidx),z_cc(kidx),out3d(i,jidx,kidx)
            elseif(outgrid_ID.eq.2)then
              ! CompP
              write(fid_outfile,210)i,jidx,kidx,x_cc(i),y_cc(jidx),p_fullmet_sp(kidx),out3d(i,jidx,kidx)
            elseif(outgrid_ID.eq.3)then
              ! MetH
              write(fid_outfile,210)i,jidx,kidx,x_submet_sp(i),y_submet_sp(jidx),z_cc(kidx),out3d(i,jidx,kidx)
            elseif(outgrid_ID.eq.4)then
              ! MetP
              write(fid_outfile,210)i,jidx,kidx,x_submet_sp(i),y_submet_sp(jidx),p_fullmet_sp(kidx),out3d(i,jidx,kidx)
            endif
          enddo
        elseif(jidx.eq.0)then
          ! line along y at iinx,kidx
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.2)then
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    x = ",x_cc(iidx),iidx
              write(outlog(io),*)"    y range = ",y_cc(1),y_cc(nymax)
            endif;enddo
          else
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    x = ",x_submet_sp(iidx),iidx
              write(outlog(io),*)"    y range = ",y_submet_sp(1),y_submet_sp(nymax)
            endif;enddo
          endif
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.3)then
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    z = ",z_cc(kidx),kidx
            endif;enddo
          else
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    p = ",p_fullmet_sp(kidx),kidx
            endif;enddo
          endif
          do j=1,ny
            if(outgrid_ID.eq.1)then
              ! CompH
              write(fid_outfile,210)iidx,j,kidx,x_cc(iidx),y_cc(j),z_cc(kidx),out3d(iidx,j,kidx)
            elseif(outgrid_ID.eq.2)then
              ! CompP
              write(fid_outfile,210)iidx,j,kidx,x_cc(iidx),y_cc(j),p_fullmet_sp(kidx),out3d(iidx,j,kidx)
            elseif(outgrid_ID.eq.3)then
              ! MetH
              write(fid_outfile,210)iidx,j,kidx,x_submet_sp(iidx),y_submet_sp(j),z_cc(kidx),out3d(iidx,j,kidx)
            elseif(outgrid_ID.eq.4)then
              ! MetP
              write(fid_outfile,210)iidx,j,kidx,x_submet_sp(iidx),y_submet_sp(j),p_fullmet_sp(kidx),out3d(iidx,j,kidx)
            endif
          enddo
        elseif(kidx.eq.0)then
          ! line along z at iinx,jind
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.2)then
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    x = ",x_cc(iidx),iidx
              write(outlog(io),*)"    y = ",y_cc(jidx),jidx
            endif;enddo
          else
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    x = ",x_submet_sp(iidx),iidx
              write(outlog(io),*)"    y = ",y_submet_sp(jidx),jidx
            endif;enddo
          endif
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.3)then
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    z range = ",z_cc(1),z_cc(nzmax)
            endif;enddo
          else
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"    p range = ",p_fullmet_sp(1),p_fullmet_sp(nzmax)
            endif;enddo
          endif
          do k=1,nz
            if(outgrid_ID.eq.1)then
              ! CompH
              write(fid_outfile,210)iidx,jidx,k,x_cc(iidx),y_cc(jidx),z_cc(k),out3d(iidx,jidx,k)
            elseif(outgrid_ID.eq.2)then
              ! CompP
              write(fid_outfile,210)iidx,jidx,k,x_cc(iidx),y_cc(jidx),p_fullmet_sp(k),out3d(iidx,jidx,k)
            elseif(outgrid_ID.eq.3)then
              ! MetH
              write(fid_outfile,210)iidx,jidx,k,x_submet_sp(iidx),y_submet_sp(jidx),z_cc(k),out3d(iidx,jidx,k)
            elseif(outgrid_ID.eq.4)then
              ! MetP
              write(fid_outfile,210)iidx,jidx,k,x_submet_sp(iidx),y_submet_sp(jidx),p_fullmet_sp(k),out3d(iidx,jidx,k)
            endif
          enddo
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)" Unknown line",iidx,jidx,kidx
            write(errlog(io),*)" For 1-d output, the dimension along the line should have index=0"
          endif;enddo
          stop 1
        endif
      elseif(outdim_ID.eq.2)then
        ! We need one of three indexes to be non-zero for line output
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Writing slice data at:"
        endif;enddo
        if(iidx.eq.0.and.jidx.eq.0)then
          ! Map-view slices
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.3)then
            ! Vert Coord is Z
            zout = z_cc(kidx)
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"       xy slice data at z=",zout,kidx
          endif;enddo
          else
            ! Vert Coord is P
            zout = p_fullmet_sp(kidx)
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"       xp slice data at p=",zout,kidx
            endif;enddo
          endif
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.2)then
            ! Horz Coords are Comp
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"Writing data to out file:",nx,ny
            endif;enddo
            if(outformat.eq.0)then
              do i=1,nx
                do j=1,ny
                  write(fid_outfile,210)i,j,kidx,x_cc(i),y_cc(j),zout,out3d(i,j,kidx)
                enddo
              enddo
            elseif(outformat.eq.1)then
              filename_root = 'outvarComp'
              call write_2D_ASCII(.true.,nx,ny,x_cc(1),y_cc(1),dx,dy,out3d(1:nx,1:ny,kidx),filename_root)
            else
              ! Unknown output format
              stop 1
            endif
          else
            ! Horz Coords are Met
            if(outformat.eq.0)then
              do i=1,nx
                do j=1,ny
                  write(fid_outfile,210)i,j,kidx,x_submet_sp(i),y_submet_sp(j),zout,out3d(i,j,kidx)
                enddo
              enddo
            elseif(outformat.eq.1)then
              dx = x_submet_sp(2)-x_submet_sp(1)
              dy = y_submet_sp(2)-y_submet_sp(1)
              filename_root = 'outvarMet'
              call write_2D_ASCII(.true.,nx,ny,x_submet_sp(1),y_submet_sp(1),dx,dy,out3d(1:nx,1:ny,kidx),filename_root)
            else
              ! Unknown output format
              stop 1
            endif
          endif
        elseif(iidx.eq.0.and.kidx.eq.0)then
          ! east-west slices
          if(outgrid_ID.eq.1)then
            ! CompH
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"       xz slice data at y=",y_cc(jidx),jidx
            endif;enddo
            do i=1,nx
              do k=1,nz
                write(fid_outfile,210)i,jidx,k,x_cc(i),y_cc(jidx),z_cc(k),out3d(i,jidx,k)
              enddo
            enddo
          elseif(outgrid_ID.eq.2)then
            ! CompP
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"       xp slice data at y=",y_cc(jidx),jidx
            endif;enddo
            do i=1,nx
              do k=1,nz
                write(fid_outfile,210)i,jidx,k,x_cc(i),y_cc(jidx),p_fullmet_sp(k),out3d(i,jidx,k)
              enddo
            enddo
          elseif(outgrid_ID.eq.3)then
            ! MetH
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"       xz slice data at y=",y_submet_sp(jidx),jidx
            endif;enddo
            do i=1,nx
              do k=1,nz
                write(fid_outfile,210)i,jidx,k,x_submet_sp(i),y_submet_sp(jidx),z_cc(k),out3d(i,jidx,k)
              enddo
            enddo
          elseif(outgrid_ID.eq.4)then
            ! MetP
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"       xz slice data at y=",y_submet_sp(jidx),jidx
            endif;enddo
            do i=1,nx
              do k=1,nz
                write(fid_outfile,210)i,jidx,k,x_submet_sp(i),y_submet_sp(jidx),p_fullmet_sp(k),out3d(i,jidx,k)
              enddo
            enddo
          endif
        elseif(jidx.eq.0.and.kidx.eq.0)then
          ! north-south slices
          if(outgrid_ID.eq.1)then
            ! CompH
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"       yz slice data at x=",x_cc(iidx),iidx
            endif;enddo
            do j=1,ny
              do k=1,nz
                write(fid_outfile,210)iidx,j,k,x_cc(iidx),y_cc(j),z_cc(k),out3d(iidx,j,k)
              enddo
            enddo
          elseif(outgrid_ID.eq.2)then
            ! CompP
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"       yp slice data at x=",x_cc(iidx),iidx
            endif;enddo
            do j=1,ny
              do k=1,nz
                write(fid_outfile,210)iidx,j,k,x_cc(iidx),y_cc(j),p_fullmet_sp(k),out3d(iidx,j,k)
              enddo
            enddo
          elseif(outgrid_ID.eq.3)then
            ! MetH
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"       yz slice data at x=",x_submet_sp(iidx),iidx
            endif;enddo
            do j=1,ny
              do k=1,nz
                write(fid_outfile,210)iidx,j,k,x_submet_sp(iidx),y_submet_sp(j),z_cc(k),out3d(iidx,j,k)
              enddo
            enddo
          elseif(outgrid_ID.eq.4)then
            ! MetP
            do io=1,MR_nio;if(VB(io).le.verbosity_info)then
              write(outlog(io),*)"       yp slice data at x=",x_submet_sp(iidx),iidx
            endif;enddo
            do j=1,ny
              do k=1,nz
                write(fid_outfile,210)iidx,j,k,x_submet_sp(iidx),y_submet_sp(j),p_fullmet_sp(k),out3d(iidx,j,k)
              enddo
            enddo
          endif
        else
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)" Unknown slice",iidx,jidx,kidx
            write(errlog(io),*)" For 2-d output, the slicing dimension should have index=0"
          endif;enddo
          stop 1
        endif
      elseif(outdim_ID.eq.3)then
        ! We need none of three indexes to be non-zero for line output (all = 0)
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)"Writing volume data"
        endif;enddo
        do i=1,nx
          do j=1,ny
            do k=1,nz
              write(fid_outfile,210)i,j,k,x_submet_sp(i),y_submet_sp(j),p_fullmet_sp(k),out3d(i,j,k)
            enddo
          enddo
        enddo
      else
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(outlog(io),*)" Unknown data output dimension.",outdim_ID
        endif;enddo
        stop 1
      endif

      close(fid_outfile)

      call MR_Reset_Memory
      if(allocated(    x_cc)) deallocate(x_cc)
      if(allocated(    y_cc)) deallocate(y_cc)
      if(allocated(    z_cc)) deallocate(z_cc)
      if(allocated(out3d_t1)) deallocate(out3d_t1)
      if(allocated(out3d_t2)) deallocate(out3d_t2)
      if(allocated(out3d   )) deallocate(out3d)

      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Program ended normally."
      endif;enddo

      end program MetRegrid

!##############################################################################
!
!     Subroutines
!
!##############################################################################

!##############################################################################
!  Read_ComdLine_InpFile
!
!  This subroutine will parse the command-line options and set up the met and
!  computational grids.  This will be determined through a command file. All
!  aspects of the run will be set, except the wind-file names, which are set
!  in a separate subroutine (GetWindFile).
!
! Example control file:
! 1 1 -135.0 90.0 0.933 6371.229  ! Proj flags and params of computational grid
! -125.0  38.0   0.0              ! x,y,z of bottom LL corner of grid (km, or deg. if latlongflag=1)  
! 15.0    15.0  25.0              ! grid x-width, y-width, depth (km, or deg. if latlonflag=1)    
! 0.5     0.5    2.0              ! DX, DY, Dz of grid cells  (km, or deg.)
! 1                               ! iHeightHandler
! 1980 5 18 15.5                  ! YYYY MM DD HH.H
! 1                               ! Output grid (1,2,3,4 for CompH,CompP,MetH,MetP)
! 5                               ! Output Var (1=H,2=Vx,3=Vy,4=Vz,5=T)
! 2 [1]                           ! Output Dimension (1,2,or 3), format
! 0 0 1                           ! i,j,k of output grid with zeros for unused dims.
! 5 25 2 2                        ! iwind iwindformat igrid iformat
! 1                               ! number of windfiles
! NCEP
!##############################################################################

      subroutine Read_ComdLine_InpFile(IsLatLon,&
                      xLL,yLL,zbot,gridwidth_x,gridwidth_y,gridwidth_z,&
                      dx,dy,dz_const,&
                      inyear,inmonth,inday,inhour,&
                      outgrid_ID,outvar_ID,outdim_ID,outformat,iidx,jidx,kidx)

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_info,verbosity_error,&
         MR_windfiles,MR_BaseYear,MR_useLeap,MR_Comp_StartHour,MR_Comp_Time_in_hours,&
         MR_iHeightHandler,MR_MAXVARS,&
           MR_Set_CompProjection,&
           MR_Allocate_FullMetFileList,&
           MR_FileIO_Error_Handler

      use projection,      only : &
         PJ_iprojflag,PJ_k0,PJ_lam0,PJ_phi0,PJ_phi1,PJ_phi2,PJ_Re,&
           PJ_Set_Proj_Params

      implicit none

      integer,parameter :: fid_ctrlfile = 10

      ! These are the variables that must be set in the input file or command line

      logical                   ,intent(out) :: IsLatLon
      real(kind=4)              ,intent(out) :: xLL
      real(kind=4)              ,intent(out) :: yLL
      real(kind=4)              ,intent(out) :: zbot
      real(kind=4)              ,intent(out) :: gridwidth_x
      real(kind=4)              ,intent(out) :: gridwidth_y
      real(kind=4)              ,intent(out) :: gridwidth_z
      real(kind=4)              ,intent(out) :: dx
      real(kind=4)              ,intent(out) :: dy
      real(kind=4)              ,intent(out) :: dz_const
      integer                   ,intent(out) :: inyear,inmonth,inday
      real(kind=8)              ,intent(out) :: inhour
      integer                   ,intent(out) :: outgrid_ID
      integer                   ,intent(out) :: outvar_ID
      integer                   ,intent(out) :: outdim_ID
      integer                   ,intent(out) :: outformat
      integer                   ,intent(out) :: iidx
      integer                   ,intent(out) :: jidx
      integer                   ,intent(out) :: kidx

      integer             :: iw
      integer             :: iwf
      integer             :: igrid
      integer             :: idf
      integer             :: iwfiles

      integer             :: nargs
      integer             :: iostatus
      integer             :: inlen
      character(len=120)  :: iomessage
      character (len=100) :: arg

      integer :: i

      character (len=100):: infile
      logical            :: IsThere
      character(len=80)  :: linebuffer080
      character(len=130) :: linebuffer130
      character(len=80)  :: Comp_projection_line
      integer            :: ilatlonflag

      integer :: io                           ! Index for output streams

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
      END INTERFACE

      ! Test read command line arguments
      nargs = command_argument_count()
      if (nargs.ne.1) then
        ! We need one command-line argument (input file name)
        ! Write usage to stdout and exit
        call Print_Usage
      endif

      ! we're using a control file.
      do io=1,MR_nio;if(VB(io).le.verbosity_info)then
        write(outlog(io),*)"Reading control file"
      endif;enddo
      call get_command_argument(1, arg, length=inlen, status=iostatus)
      if(iostatus.ne.0)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR : Could not read first command-line argument"
          write(errlog(io),*)" arg = ",arg
        endif;enddo
        call Print_Usage
      endif

      read(arg,*,iostat=iostatus,iomsg=iomessage) infile
      linebuffer080 = arg(1:80)
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      ! Error-check infile
      inquire( file=infile, exist=IsThere )
      if(.not.IsThere)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"MR ERROR: Cannot find input file"
        endif;enddo
        stop 1
      endif
      open(unit=fid_ctrlfile,file=infile,status='old',err=1900)

      ! Line 1: projection parameters for computational grid
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      Comp_projection_line = linebuffer080
      read(Comp_projection_line,*,iostat=iostatus,iomsg=iomessage)ilatlonflag
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      if (ilatlonflag.eq.0) then
        ! expecting input variables to be in the same projection as
        ! specified by iprojflag and parameters
        IsLatLon          = .false.
      else
        ! expecting input variables to be in lat/lon
       IsLatLon          = .true.
      endif

      ! Line 2: Lower-left coordinate of computational grid (x,y,z)
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) xLL, yLL, zbot
      if(iostatus.ne.0)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)  'error reading xLL, yLL, zbot.'
          write(errlog(io),*)  'You entered: ',linebuffer080
          write(errlog(io),*)  'Program stopped.'
        endif;enddo
        call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        stop 1
      endif

      ! Line 3: Widths of computational grid
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) gridwidth_x, gridwidth_y, gridwidth_z
      if(iostatus.ne.0.or.&
         gridwidth_x.lt.0.0_8.or.gridwidth_y.lt.0.0_8.or.gridwidth_z.lt.0.0_8)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)  'error reading length, width and height of model domain.'
          write(errlog(io),*)  'You entered: ', linebuffer080
          write(errlog(io),*)  'Program stopped'
        endif;enddo
        call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        stop 1
      endif

      ! Line 4: dx, dy, dz of computational grid
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) dx, dy, dz_const
      if(iostatus.ne.0.or.&
         dx.lt.0.0_8.or.dy.lt.0.0_8.or.dz_const.lt.0.0_8)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)  'error reading dx, dy, or dz.'
          write(errlog(io),*)  'You entered: ',linebuffer080
          write(errlog(io),*)  'Program stopped.'
        endif;enddo
        call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        stop 1
      endif

      ! Line 5: iHeightHandler
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage) linebuffer080
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) MR_iHeightHandler
      if(iostatus.ne.0.or.&
         MR_iHeightHandler.lt.1.or.MR_iHeightHandler.gt.2)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)  'error reading iHeightHandler.'
          write(errlog(io),*)  'You gave: ', linebuffer080
          write(errlog(io),*)  'Program stopped.'
        endif;enddo
        call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        stop 1
      endif

      ! Line 6: YYYY MM DD HH.H
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) inyear,inmonth,inday,inhour
      if(iostatus.ne.0.or.&
         inmonth.lt.1.or.inmonth.gt.12.or.&
         inday.lt.1.or.inday.gt.31.or.&
         inhour.lt.0.0_8.or.inhour.gt.24.0_8)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)  'error reading inyear,inmonth,inday,inhour.'
          write(errlog(io),*)  'You gave: ', linebuffer080
          write(errlog(io),*)  'Program stopped.'
        endif;enddo
        call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        stop 1
      endif

      ! Line 7: Output grid code (1=CompH,2=CompP,3=MetH,4=MetP)
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) outgrid_ID
      if(iostatus.ne.0.or.&
         outgrid_ID.lt.1.or.outgrid_ID.gt.4)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)  'error reading outgrid_ID.'
          write(errlog(io),*)  'You gave: ', linebuffer080
          write(errlog(io),*)  'Program stopped.'
        endif;enddo
        call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        stop 1
      endif

      ! Line 8: Output var ID
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) outvar_ID
      if(iostatus.ne.0.or.&
         outvar_ID.lt.1.or.outvar_ID.gt.MR_MAXVARS)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)  'error reading outvar_ID.'
          write(errlog(io),*)  'You gave: ', linebuffer080
          write(errlog(io),*)  'Program stopped.'
        endif;enddo
        call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        stop 1
      endif

      ! Line 9: Output dim ID
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) outdim_ID
      if(iostatus.ne.0.or.&
         outdim_ID.lt.1.or.outdim_ID.gt.3)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)  'error reading outdim_ID.'
          write(errlog(io),*)  'You gave: ', linebuffer080
          write(errlog(io),*)  'Program stopped.'
        endif;enddo
        call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        stop 1
      endif
      ! Try to read the optional format code.  If not present then = 0
      !   0 = ASCII dump with i,j,k,x,y,z,var
      !   1 = ESRI ASCII (only works for 2d xy slices)
      ! Later, we might add other options such as kml, png, shapefile similar to Ash3d_PostProc
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) outdim_ID, outformat
      if(iostatus.ne.0.or.&
         outformat.lt.0.or.outformat.gt.1)then
        outformat = 0
      endif
      if(outdim_ID.eq.1.or.outdim_ID.eq.3)then
        ! For now; ESRI ASCII output requires 2-d output
        ! Reverting to full ASCII dump
        outformat = 0
      endif


      ! Line 10: i,j,k index with 0's for unused dimensions
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) iidx,jidx,kidx
      if(iostatus.ne.0.or.&
         iidx.lt.0.or.jidx.lt.0.or.kidx.lt.0)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)  'error reading iidx,jidx,kidx.'
          write(errlog(io),*)  'You gave: ', linebuffer080
          write(errlog(io),*)  'Program stopped.'
        endif;enddo
        call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        stop 1
      endif

      ! Line 11: Met specification=> iw, iwf, igrid, idf
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) iw,iwf,igrid,idf
      if(iostatus.ne.0.or.&
         iw.lt.1.or.iw.gt.5.or.&
         iwf.lt.0.or.iwf.gt.50.or.&
         idf.lt.1.or.idf.gt.5)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)  'error reading iw,iwf,igrid,idf.'
          write(errlog(io),*)  'You gave: ', linebuffer080
          write(errlog(io),*)  'wind format must be one of 1,2,3,4, or 5'
          write(errlog(io),*)  'wind product must be in range 1-50'
          write(errlog(io),*)  'wind data format must be one of 1,2, or 3'
          write(errlog(io),*)  'Program stopped.'
        endif;enddo
        call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        stop 1
      endif

      ! Line 12: number of met files
      read(fid_ctrlfile,'(a80)',iostat=iostatus,iomsg=iomessage)linebuffer080
      if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
      read(linebuffer080,*,iostat=iostatus,iomsg=iomessage) iwfiles
      if(iostatus.ne.0.or.&
         iwfiles.lt.1)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)  'error reading iwfiles.'
          write(errlog(io),*)  'You gave: ', linebuffer080
          write(errlog(io),*)  'Program stopped.'
        endif;enddo
        call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
        stop 1
      endif

      MR_Comp_StartHour = HS_hours_since_baseyear(inyear,inmonth,inday,inhour,&
                                                MR_BaseYear,MR_useLeap)
      MR_Comp_Time_in_hours = 0.0_8

      ! Set Projection Parameters
      call PJ_Set_Proj_Params(Comp_projection_line)
      call MR_Set_CompProjection(IsLatLon,PJ_iprojflag,PJ_lam0, &
                             PJ_phi0,PJ_phi1,PJ_phi2,       &
                             PJ_k0,PJ_Re)
      ! Set up memory for windfiles
      call MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)

      ! Line 13: starting to read the windfiles
      do i=1,iwfiles
        read(fid_ctrlfile,'(a130)',iostat=iostatus,iomsg=iomessage)linebuffer130
        if(iostatus.ne.0) call MR_FileIO_Error_Handler(iostatus,linebuffer130(1:80),iomessage)
        read(linebuffer130,'(a130)',iostat=iostatus,iomsg=iomessage)MR_windfiles(i)
        inquire( file=infile, exist=IsThere )
        if(iostatus.ne.0.or.&
           .not.IsThere)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)  'error reading MR_windfiles(i)'
            write(errlog(io),*)  'You gave: ', i, linebuffer080
            write(errlog(io),*)  'Inquire test= ',IsThere
            write(errlog(io),*)  'Program stopped.'
          endif;enddo
          call MR_FileIO_Error_Handler(iostatus,linebuffer080,iomessage)
          stop 1
        endif
        do io=1,MR_nio;if(VB(io).le.verbosity_info)then
          write(outlog(io),*)i,trim(adjustl(MR_windfiles(i)))
        endif;enddo
      enddo

      ! Successfully read control file; closing it
      close(fid_ctrlfile)

      return

!******************************************************************************
!     ERROR TRAPS

1900  do io=1,MR_nio;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)  'error: cannot find input file: ',infile
        write(errlog(io),*)  'Program stopped'
      endif;enddo
      stop 1

      end subroutine Read_ComdLine_InpFile

!##############################################################################
!
!  write_2D_ASCII
!
!  Arguments:
!    IsLatLon      = logical
!    nx            = x length of output array OutVar
!    ny            = y length of output array OutVar
!    xLL,yLL       = lower-left coordinates
!    dx,dy         = x and y spacing
!    OutVar        = 2-d array to be written to ASCII file
!    filename_root = root name of file (20 characters)
!
!  Subroutine that writes out 2-D arrays in ESRI ASCII raster format
!  This format can be post-processed with gmt converting to grid files with
!  gmt grdconvert out.dat=ef out.grd
!
!##############################################################################

      subroutine write_2D_ASCII(IsLatLon,nx,ny,xLL,yLL,dx,dy,OutVar,filename_root)

      use MetReader,       only : &
         MR_nio,VB,errlog,verbosity_error

      logical          ,intent(in) :: IsLatLon
      integer          ,intent(in) :: nx,ny
      real(kind=4)     ,intent(in) :: xLL,yLL
      real(kind=4)     ,intent(in) :: dx,dy
      real(kind=4)     ,intent(in) :: OutVar(nx,ny)
      character(len=20),intent(in) :: filename_root

      real(kind=4),parameter :: KM_2_M = 1000.0_4

      integer :: i,j
      character(len=50)  :: filename_out

      integer,parameter :: fid_outfile = 20

      integer :: io                           ! Index for output streams

      write(filename_out,*)trim(adjustl(filename_root)),'.dat'

      open(unit=fid_outfile,file=trim(adjustl(filename_out)), status='replace',action='write',err=2500)

      write(fid_outfile,3000) nx        ! write header values
      write(fid_outfile,3001) ny
      if (IsLatLon) then
        write(fid_outfile,3002) xLL
        write(fid_outfile,3003) yLL
        write(fid_outfile,3004) dx,dy
      else
        write(fid_outfile,3002) xLL*KM_2_M    ! convert xLL from km to meters so ArcMap can read it
        write(fid_outfile,3003) yLL*KM_2_M    ! same with yLL
        write(fid_outfile,3004) dx*KM_2_M,dy*KM_2_M    ! and with dx and dy
      endif
      write(fid_outfile,3005)'-9999.'

      ! Write out arrays
      do j=ny,1,-1
        write(fid_outfile,3006) (OutVar(i,j), i=1,nx)
        !write(fid_outfile,'(g0)') ''          ! make a blank line between rows
        write(fid_outfile,*)" "               ! make a blank line between rows
      enddo

      close(fid_outfile)

!     format statements
3000  format('NCOLS ',i5)
3001  format('NROWS ',i5)
3002  format('XLLCORNER ',f15.3)
3003  format('YLLCORNER ',f15.3)
3004  format('CELLSIZE ',2f15.3)
3005  format('NODATA_VALUE ',a6)
3006  format(10f18.6)
      return

!     Error traps
2500  do io=1,MR_nio;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error opening output file ASCII_output_file.txt.  Program stopped'
      endif;enddo
      stop 1

      end subroutine write_2D_ASCII

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
         MR_nio,VB,errlog,verbosity_error

      implicit none

      integer :: io                           ! Index for output streams

      do io=1,MR_nio;if(VB(io).le.verbosity_error)then
        write(errlog(io),*)"Too few command-line arguments:"
        write(errlog(io),*)"  Usage: MetRegrid command_file"
        write(errlog(io),*)"   Where the command file has the following format"
        write(errlog(io),*)"--------------------------------------------------"
        write(errlog(io),*)"1 1 -135.0 90.0 0.933 6371.229  ! Proj flags and params of computational grid"
        write(errlog(io),*)"-125.0  38.0   0.0              ! x,y,z of bottom LL corner of grid (km, or deg. if latlongflag=1)"
        write(errlog(io),*)"15.0    15.0  25.0              ! grid x-width, y-width, depth (km, or deg. if latlonflag=1)"
        write(errlog(io),*)"0.5     0.5    2.0              ! DX, DY, Dz of grid cells  (km, or deg.)"
        write(errlog(io),*)"1                               ! iHeightHandler"
        write(errlog(io),*)"1980 5 18 15.5                  ! YYYY MM DD HH.H"
        write(errlog(io),*)"1                               ! Output grid (1,2,3,4 for CompH,CompP,MetH,MetP)"
        write(errlog(io),*)"5                               ! Output Var (1=H,2=Vx,3=Vy,4=Vz,5=T)"
        write(errlog(io),*)"2                               ! Output Dimension (0,1,2,or 3)"
        write(errlog(io),*)"0 0 1                           ! i,j,k of output grid with zeros for unused dims."
        write(errlog(io),*)"5 25 2 2                        ! iwind iwindformat igrid iformat"
        write(errlog(io),*)"1                               ! number of windfiles"
        write(errlog(io),*)"NCEP"
        write(errlog(io),*)"--------------------------------------------------"
      endif;enddo
      stop 1

      end subroutine Print_Usage

!##############################################################################
!##############################################################################
