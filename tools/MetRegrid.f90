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

      use MetReader

      implicit none

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

      real(kind=4),dimension(:,:,:),allocatable :: out3d_t1,out3d_t2,out3d
      

      INTERFACE
      subroutine Read_ComdLine_InpFile(IsLatLon,xLL,yLL,zbot,&
                      gridwidth_x,gridwidth_y,gridwidth_z,&
                      dx,dy,dz_const,&
                      inyear,inmonth,inday,inhour,&
                      outgrid_ID,outvar_ID,outdim_ID,iidx,jidx,kidx)
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
           integer                  ,intent(out) :: iidx
           integer                  ,intent(out) :: jidx
           integer                  ,intent(out) :: kidx
        end subroutine Read_ComdLine_InpFile
      END INTERFACE

      call Read_ComdLine_InpFile(IsLatLon,xLL,yLL,zbot,&
                      gridwidth_x,gridwidth_y,gridwidth_z,&
                      dx,dy,dz_const,&
                      inyear,inmonth,inday,inhour,&
                      outgrid_ID,outvar_ID,outdim_ID,iidx,jidx,kidx)

        ! Check for existance and compatibility with simulation time requirements
      call MR_Read_Met_DimVars(inyear)

      call MR_Set_Met_Times(MR_Comp_StartHour,MR_Comp_Time_in_hours)

      ! Now set up the computational grid
      IsPeriodic  = .false.
      if(IsLatLon)then
        if(xLL.lt.-360.0_4)then
          write(MR_global_error,*)"ERROR: xLL must be >-360 and <360"
          write(MR_global_error,*)"       From control file, xLL = ",xLL
          stop 1
        elseif(xLL.lt.0.0_4)then
          ! shift negative lon to range 0->360
          xLL=xLL+360.0_4
          write(MR_global_info,*)"Resetting to range 0->360"
        elseif(xLL.gt.360.0_4)then 
          write(MR_global_info,*)"WARNING: xLL is greater than 360"
          write(MR_global_info,*)"Resetting to range 0->360"
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
        write(MR_global_info,*)"Planning to read variable:",&
                  trim(adjustl(Met_var_NC_names(outvar_ID)))," (",&
                  trim(adjustl(Met_var_GRIB_names(outvar_ID))),")"
      else
        write(MR_global_info,*)"Variable is not available for this windfile"
        write(MR_global_info,*)"  outvar_ID = ",outvar_ID
        write(MR_global_info,*)"  var name  = ",trim(adjustl(Met_var_NC_names(outvar_ID)))," (",&
                                   trim(adjustl(Met_var_GRIB_names(outvar_ID))),")"
        stop 1
      endif

      ! Output grid (1,2,3,4 for CompH,CompP,MetH,MetP)
      if(outgrid_ID.eq.1)then
        write(MR_global_info,*)"Output will be interpolated onto the CompH grid"
        IsH = .true.
      elseif(outgrid_ID.eq.2)then
        write(MR_global_info,*)"Output will be interpolated onto the CompP grid"
        IsH = .false.
      elseif(outgrid_ID.eq.3)then
        write(MR_global_info,*)"Output will be interpolated onto the MetH grid"
        IsH = .true.
      elseif(outgrid_ID.eq.4)then
        write(MR_global_info,*)"Output will be interpolated onto the MetP grid"
        IsH = .false.
      else
        write(MR_global_error,*)"ERROR: outgrid_ID not recognized; expecting 1-4"
        write(MR_global_error,*)"       outgrid_ID = ",outgrid_ID
        stop 1
      endif

      ! Now actually interpolate in time
      MR_iMetStep_Now = 1
      tfrac = (MR_Comp_StartHour-MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now))/ &
               MR_MetStep_Interval(MR_iMetStep_Now)
      write(MR_global_info,*)"Hour = ",inhour
      write(MR_global_info,*)"Using a fractional step of ",tfrac
      if(outgrid_ID.eq.1)then
        ! We need CompH grids
        nx = nxmax
        ny = nymax
        nz = nzmax
        write(MR_global_info,*)" Output array is of size ",nx,ny,nz
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
        write(MR_global_info,*)" Output array is of size ",nx,ny,nz
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
        write(MR_global_info,*)" Output array is of size ",nx,ny,nz
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
        write(MR_global_info,*)" Output array is of size ",nx,ny,nz
        allocate(out3d_t1(nx,ny,nz))
        allocate(out3d_t2(nx,ny,nz))
        allocate(   out3d(nx,ny,nz))
        call MR_Read_3d_MetP_Variable(outvar_ID,MR_iMetStep_Now)
        out3d_t1(1:nx,1:ny,1:nz) = MR_dum3d_metP(1:nx,1:ny,1:nz)
        call MR_Read_3d_MetP_Variable(outvar_ID,MR_iMetStep_Now+1)
        out3d_t2(1:nx,1:ny,1:nz) = MR_dum3d_metP(1:nx,1:ny,1:nz)
        out3d = out3d_t1 + (out3d_t2-out3d_t1)*real(tfrac,kind=4)
      else
        write(MR_global_error,*)"ERROR: outgrid_ID not recognized; expecting 1-4"
        write(MR_global_error,*)"       outgrid_ID = ",outgrid_ID
        stop 1
      endif

      ! What are we writing out
200   format(a3,i0.2,a1,i1,a1,i1,a4)

      write(outfilename,200) 'out',outvar_ID,'_',outgrid_ID,'_',outdim_ID,'.dat'
      open(unit=20,file=outfilename)
      ! First check the indexes for validity
      if(iidx.lt.0.or.jidx.lt.0.or.kidx.lt.0)then
        write(MR_global_error,*)"ERROR: Output index must be non-negative"
        write(MR_global_error,*)"   iidx = ", iidx
        write(MR_global_error,*)"   jidx = ", jidx 
        write(MR_global_error,*)"   kidx = ", kidx
        stop 1
      elseif(iidx.gt.nx)then
        write(MR_global_error,*)"ERROR: Output index must not be greater than output array."
        write(MR_global_error,*)"   nx   = ", nx
        write(MR_global_error,*)"   iidx = ", iidx
        stop 1
      elseif(jidx.gt.ny)then
        write(MR_global_error,*)"ERROR: Output index must not be greater than output array."
        write(MR_global_error,*)"   ny   = ", ny
        write(MR_global_error,*)"   jidx = ", kidx
        stop 1
      elseif(kidx.gt.nz)then
        write(MR_global_error,*)"ERROR: Output index must not be greater than output array."
        write(MR_global_error,*)"   nz   = ", nz
        write(MR_global_error,*)"   kidx = ", kidx
        stop 1
      endif
210   format(3i8,4f12.4)
      if(outdim_ID.eq.0)then
        ! We need all three indexes to be non-zero for point output
        if(iidx.eq.0.or.jidx.eq.0.or.kidx.eq.0)then
          write(MR_global_info,*)" Unknown point",iidx,jidx,kidx
          write(MR_global_info,*)" For 0-d output, each index must be specified"
          stop 1
        endif
        write(MR_global_info,*)"Writing point data at:"
        write(MR_global_info,*)"    ", iidx
        write(MR_global_info,*)"    ", jidx 
        write(MR_global_info,*)"    ", kidx
        if(outgrid_ID.eq.1)then
         ! CompH
         write(20,210)iidx,jidx,kidx,x_cc(iidx),y_cc(jidx),z_cc(kidx),out3d(iidx,jidx,kidx)
        elseif(outgrid_ID.eq.1)then
         ! CompP
         write(20,210)iidx,jidx,kidx,x_cc(iidx),y_cc(jidx),p_fullmet_sp(kidx),out3d(iidx,jidx,kidx)
        elseif(outgrid_ID.eq.3)then
         ! MetH
         write(20,210)iidx,jidx,kidx,x_submet_sp(iidx),y_submet_sp(jidx),z_cc(kidx),out3d(iidx,jidx,kidx)
        elseif(outgrid_ID.eq.4)then
         ! MetP
         write(20,210)iidx,jidx,kidx,x_submet_sp(iidx),y_submet_sp(jidx),p_fullmet_sp(kidx),out3d(iidx,jidx,kidx)
        endif
      elseif(outdim_ID.eq.1)then
        ! We need two of three indexes to be non-zero for line output
        if((iidx.eq.0.and.(jidx.eq.0.or.kidx.eq.0)).or.&
           (jidx.eq.0.and.kidx.eq.0))then
          write(MR_global_info,*)" Unknown line",iidx,jidx,kidx
          write(MR_global_info,*)" For 1-d output, two of three indexes must be specified"
          stop 1
        endif
        write(MR_global_info,*)"Writing line data at:"
        if(iidx.eq.0)then
          ! line along x at jinx,kind
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.2)then
            write(MR_global_info,*)"    y = ",y_cc(jidx),jidx
            write(MR_global_info,*)"    x range = ",x_cc(1),x_cc(nxmax)
          else
            write(MR_global_info,*)"    y = ",y_submet_sp(jidx),jidx
            write(MR_global_info,*)"    x range = ",y_submet_sp(1),y_submet_sp(nymax)
          endif
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.3)then
            write(MR_global_info,*)"    z = ",z_cc(kidx),kidx
          else
            write(MR_global_info,*)"    p = ",p_fullmet_sp(kidx),kidx
          endif

          do i=1,nx
            if(outgrid_ID.eq.1)then
              ! CompH
              write(20,210)i,jidx,kidx,x_cc(i),y_cc(jidx),z_cc(kidx),out3d(i,jidx,kidx)
            elseif(outgrid_ID.eq.2)then
              ! CompP
              write(20,210)i,jidx,kidx,x_cc(i),y_cc(jidx),p_fullmet_sp(kidx),out3d(i,jidx,kidx)
            elseif(outgrid_ID.eq.3)then
              ! MetH
              write(20,210)i,jidx,kidx,x_submet_sp(i),y_submet_sp(jidx),z_cc(kidx),out3d(i,jidx,kidx)
            elseif(outgrid_ID.eq.4)then
              ! MetP
              write(20,210)i,jidx,kidx,x_submet_sp(i),y_submet_sp(jidx),p_fullmet_sp(kidx),out3d(i,jidx,kidx)
            endif
          enddo
        elseif(jidx.eq.0)then
          ! line along y at iinx,kidx
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.2)then
            write(MR_global_info,*)"    x = ",x_cc(iidx),iidx
            write(MR_global_info,*)"    y range = ",y_cc(1),y_cc(nymax)
          else
            write(MR_global_info,*)"    x = ",x_submet_sp(iidx),iidx
            write(MR_global_info,*)"    y range = ",y_submet_sp(1),y_submet_sp(nymax)
          endif
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.3)then
            write(MR_global_info,*)"    z = ",z_cc(kidx),kidx
          else
            write(MR_global_info,*)"    p = ",p_fullmet_sp(kidx),kidx
          endif
          do j=1,ny
            if(outgrid_ID.eq.1)then
              ! CompH
              write(20,210)iidx,j,kidx,x_cc(iidx),y_cc(j),z_cc(kidx),out3d(iidx,j,kidx)
            elseif(outgrid_ID.eq.2)then
              ! CompP
              write(20,210)iidx,j,kidx,x_cc(iidx),y_cc(j),p_fullmet_sp(kidx),out3d(iidx,j,kidx)
            elseif(outgrid_ID.eq.3)then
              ! MetH
              write(20,210)iidx,j,kidx,x_submet_sp(iidx),y_submet_sp(j),z_cc(kidx),out3d(iidx,j,kidx)
            elseif(outgrid_ID.eq.4)then
              ! MetP
              write(20,210)iidx,j,kidx,x_submet_sp(iidx),y_submet_sp(j),p_fullmet_sp(kidx),out3d(iidx,j,kidx)
            endif
          enddo
        elseif(kidx.eq.0)then
          ! line along z at iinx,jind
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.2)then
            write(MR_global_info,*)"    x = ",x_cc(iidx),iidx
            write(MR_global_info,*)"    y = ",y_cc(jidx),jidx
          else
            write(MR_global_info,*)"    x = ",x_submet_sp(iidx),iidx
            write(MR_global_info,*)"    y = ",y_submet_sp(jidx),jidx
          endif
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.3)then
            write(MR_global_info,*)"    z range = ",z_cc(1),z_cc(nzmax)
          else
            write(MR_global_info,*)"    p range = ",p_fullmet_sp(1),p_fullmet_sp(nzmax)
          endif
          do k=1,nz
            if(outgrid_ID.eq.1)then
              ! CompH
              write(20,210)iidx,jidx,k,x_cc(iidx),y_cc(jidx),z_cc(k),out3d(iidx,jidx,k)
            elseif(outgrid_ID.eq.2)then
              ! CompP
              write(20,210)iidx,jidx,k,x_cc(iidx),y_cc(jidx),p_fullmet_sp(k),out3d(iidx,jidx,k)
            elseif(outgrid_ID.eq.3)then
              ! MetH
              write(20,210)iidx,jidx,k,x_submet_sp(iidx),y_submet_sp(jidx),z_cc(k),out3d(iidx,jidx,k)
            elseif(outgrid_ID.eq.4)then
              ! MetP
              write(20,210)iidx,jidx,k,x_submet_sp(iidx),y_submet_sp(jidx),p_fullmet_sp(k),out3d(iidx,jidx,k)
            endif
          enddo
        else
          write(MR_global_info,*)" Unknown line",iidx,jidx,kidx
          write(MR_global_info,*)" For 1-d output, the dimension along the line should have index=0"
          stop 1
        endif
      elseif(outdim_ID.eq.2)then
        write(MR_global_info,*)"Writing slice data at:"
        if(iidx.eq.0.and.jidx.eq.0)then
          ! Map-view slices
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.3)then
            ! Vert Coord is Z
            zout = z_cc(kidx)
            write(MR_global_info,*)"       xy slice data at z=",zout,kidx
          else
            ! Vert Coord is P
            zout = p_fullmet_sp(kidx)
            write(MR_global_info,*)"       xp slice data at p=",zout,kidx
          endif
          if(outgrid_ID.eq.1.or.outgrid_ID.eq.2)then
            ! Horz Coords are Comp
            write(MR_global_info,*)"Writing data to out file:",nx,ny
            do i=1,nx
              do j=1,ny
                write(20,210)i,j,kidx,x_cc(i),y_cc(j),zout,out3d(i,j,kidx)
              enddo
            enddo
          else
            ! Horz Coords are Met
            do i=1,nx
              do j=1,ny
                write(20,210)i,j,kidx,x_submet_sp(i),y_submet_sp(j),zout,out3d(i,j,kidx)
              enddo
            enddo
          endif
        elseif(iidx.eq.0.and.kidx.eq.0)then
          ! east-west slices
          if(outgrid_ID.eq.1)then
            ! CompH
            write(MR_global_info,*)"       xz slice data at y=",y_cc(jidx),jidx
            do i=1,nx
              do k=1,nz
                write(20,210)i,jidx,k,x_cc(i),y_cc(jidx),z_cc(k),out3d(i,jidx,k)
              enddo
            enddo
          elseif(outgrid_ID.eq.2)then
            ! CompP
            write(MR_global_info,*)"       xp slice data at y=",y_cc(jidx),jidx
            do i=1,nx
              do k=1,nz
                write(20,210)i,jidx,k,x_cc(i),y_cc(jidx),p_fullmet_sp(k),out3d(i,jidx,k)
              enddo
            enddo
          elseif(outgrid_ID.eq.3)then
            ! MetH
            write(MR_global_info,*)"       xz slice data at y=",y_submet_sp(jidx),jidx
            do i=1,nx
              do k=1,nz
                write(20,210)i,jidx,k,x_submet_sp(i),y_submet_sp(jidx),z_cc(k),out3d(i,jidx,k)
              enddo
            enddo
          elseif(outgrid_ID.eq.4)then
            ! MetP
            write(MR_global_info,*)"       xz slice data at y=",y_submet_sp(jidx),jidx
            do i=1,nx
              do k=1,nz
                write(20,210)i,jidx,k,x_submet_sp(i),y_submet_sp(jidx),p_fullmet_sp(k),out3d(i,jidx,k)
              enddo
            enddo
          endif
        elseif(jidx.eq.0.and.kidx.eq.0)then
          ! north-south slices
          if(outgrid_ID.eq.1)then
            ! CompH
            write(MR_global_info,*)"       yz slice data at x=",x_cc(iidx),iidx
            do j=1,ny
              do k=1,nz
                write(20,210)iidx,j,k,x_cc(iidx),y_cc(j),z_cc(k),out3d(iidx,j,k)
              enddo
            enddo
          elseif(outgrid_ID.eq.2)then
            ! CompP
            write(MR_global_info,*)"       yp slice data at x=",x_cc(iidx),iidx
            do j=1,ny
              do k=1,nz
                write(20,210)iidx,j,k,x_cc(iidx),y_cc(j),p_fullmet_sp(k),out3d(iidx,j,k)
              enddo
            enddo
          elseif(outgrid_ID.eq.3)then
            ! MetH
            write(MR_global_info,*)"       yz slice data at x=",x_submet_sp(iidx),iidx
            do j=1,ny
              do k=1,nz
                write(20,210)iidx,j,k,x_submet_sp(iidx),y_submet_sp(j),z_cc(k),out3d(iidx,j,k)
              enddo
            enddo
          elseif(outgrid_ID.eq.4)then
            ! MetP
            write(MR_global_info,*)"       yp slice data at x=",x_submet_sp(iidx),iidx
            do j=1,ny
              do k=1,nz
                write(20,210)iidx,j,k,x_submet_sp(iidx),y_submet_sp(j),p_fullmet_sp(k),out3d(iidx,j,k)
              enddo
            enddo
          endif
        else
          write(MR_global_info,*)" Unknown slice",iidx,jidx,kidx
          write(MR_global_info,*)" For 2-d output, the slicing dimension should have index=0"
          stop 1
        endif
      elseif(outdim_ID.eq.3)then
        write(MR_global_info,*)"Writing volume data"
        do i=1,nx
          do j=1,ny
            do k=1,nz
              write(20,210)i,j,k,x_submet_sp(i),y_submet_sp(j),p_fullmet_sp(k),out3d(i,j,k)
            enddo
          enddo
        enddo
      else
        write(MR_global_info,*)" Unknown data output dimension.",outdim_ID
        stop 1
      endif
      close(20)

      write(MR_global_info,*)"Program ended normally."

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
! 2                               ! Output Dimension (1,2,or 3)
! 0 0 1                           ! i,j,k of output grid with zeros for unused dims.
! 5 25 2 2                        ! iwind iwindformat igrid iformat
! 1                               ! number of windfiles
! NCEP
!##############################################################################

      subroutine Read_ComdLine_InpFile(IsLatLon,&
                      xLL,yLL,zbot,gridwidth_x,gridwidth_y,gridwidth_z,&
                      dx,dy,dz_const,&
                      inyear,inmonth,inday,inhour,&
                      outgrid_ID,outvar_ID,outdim_ID,iidx,jidx,kidx)

      use MetReader
      use projection

      implicit none

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
      integer                   ,intent(out) :: iidx
      integer                   ,intent(out) :: jidx
      integer                   ,intent(out) :: kidx

      integer             :: iw
      integer             :: iwf
      integer             :: igrid
      integer             :: idf
      integer             :: iwfiles

      integer             :: nargs
      integer             :: status
      character (len=100) :: arg

      integer :: i

      character (len=100):: infile
      logical            :: IsThere
      character(len=80)  :: linebuffer080
      character(len=130) :: linebuffer130
      character(len=80)  :: Comp_projection_line
      integer            :: ilatlonflag

      INTERFACE
        subroutine write_usage
        end subroutine write_usage
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
        call write_usage
      endif

      ! we're using a control file.
      write(MR_global_production,*)"Reading control file"
      call get_command_argument(1, arg, status)
      read(arg,*) infile
      inquire( file=infile, exist=IsThere )
      if(.not.IsThere)then
        write(MR_global_error,*)"ERROR: Cannot find input file"
        stop 1
      endif
      open(unit=10,file=infile,status='old',err=1900)

      ! Line 1: projection parameters for computational grid
      read(10,'(a80)') linebuffer080
      Comp_projection_line = linebuffer080
      read(Comp_projection_line,*)ilatlonflag
      if (ilatlonflag.eq.0) then
        ! expecting input variables to be in the same projection as
        ! specified by iprojflag and parameters
        IsLatLon          = .false.
      else
        ! expecting input variables to be in lat/lon
       IsLatLon          = .true.
      endif

      ! Line 2: Lower-left coordinate of computational grid (x,y,z)
      read(10,'(a80)') linebuffer080
      read(linebuffer080,*,err=1901) xLL, yLL, zbot

      ! Line 3: Widths of computational grid
      read(10,'(a80)') linebuffer080
      read(linebuffer080,*,err=1902) gridwidth_x, gridwidth_y, gridwidth_z

      ! Line 4: dx, dy, dz of computational grid
      read(10,'(a80)') linebuffer080
      read(linebuffer080,*,err=1904) dx, dy, dz_const

      ! Line 5: iHeightHandler
      read(10,'(a80)') linebuffer080
      read(linebuffer080,*,err=1905) MR_iHeightHandler

      ! Line 6: YYYY MM DD HH.H
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*) inyear,inmonth,inday,inhour

      ! Line 7: Output grid code (1=CompH,2=CompP,3=MetH,4=MetP)
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*) outgrid_ID

      ! Line 8: Output var ID
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*) outvar_ID

      ! Line 9: Output dim ID
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*) outdim_ID

      ! Line 10: i,j,k index with 0's for unused dimensions
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*) iidx,jidx,kidx

      ! Line 11: Met specification=> iw, iwf, igrid, idf
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*) iw,iwf,igrid,idf

      ! Line 12: number of met files
      read(10,'(a80)')linebuffer080
      read(linebuffer080,*) iwfiles

      MR_Comp_StartHour = HS_hours_since_baseyear(inyear,inmonth,inday,inhour,&
                                                MR_BaseYear,MR_useLeap)
      MR_Comp_Time_in_hours = 0.0_8

      ! Set Projection Parameters
      call PJ_Set_Proj_Params(Comp_projection_line)
      call MR_Set_CompProjection(IsLatLon,PJ_iprojflag,PJ_lam0, &
                             PJ_phi0,PJ_phi1,PJ_phi2,       &
                             PJ_k0,PJ_radius_earth)
      ! Set up memory for windfiles
      call MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)

      ! Line 13: starting to read the windfiles
      do i=1,iwfiles
        read(10,'(a130)')linebuffer130
        read(linebuffer130,'(a130)')MR_windfiles(i)
        write(MR_global_info,*)i,trim(adjustl(MR_windfiles(i)))
      enddo

      return

!******************************************************************************
!     ERROR TRAPS

      !ERROR TRAPS TO STDIN
1900  write(MR_global_error,*)  'error: cannot find input file: ',infile
      write(MR_global_error,*)  'Program stopped'
      stop 1

1901  write(MR_global_error,*)  'error reading xLL or yLL.'
      write(MR_global_error,*)  'You entered: ',linebuffer080
      write(MR_global_error,*)  'Program stopped.'
      stop 1

1902  write(MR_global_error,*)  'error reading width and height of model domain.'
      write(MR_global_error,*)  'You entered: ', linebuffer080
      write(MR_global_error,*)  'Program stopped'
      stop 1

1904  write(MR_global_error,*)  'error reading dx or dy.'
      write(MR_global_error,*)  'You entered: ',linebuffer080
      write(MR_global_error,*)  'Program stopped.'
      stop 1

1905  write(MR_global_error,*)  'error reading iHeightHandler.'
      write(MR_global_error,*)  'You gave: ', linebuffer080
      write(MR_global_error,*)  'Program stopped.'
      stop 1

      end subroutine Read_ComdLine_InpFile

!##############################################################################
!##############################################################################
!  write_usage
!
!  This subroutine is called if there is an error reading the command-line.
!  Expected usage is written to stdout and the program exits.
!
!##############################################################################

      subroutine write_usage

      use MetReader

      implicit none

        write(MR_global_info,*)"Too few command-line arguments:"
        write(MR_global_info,*)"  Usage: MetRegrid command_file"
        write(MR_global_info,*)"   Where the command file has the following format"
        write(MR_global_info,*)"--------------------------------------------------"
        write(MR_global_info,*)"1 1 -135.0 90.0 0.933 6371.229  ! Proj flags and params of computational grid"
        write(MR_global_info,*)"-125.0  38.0   0.0              ! x,y,z of bottom LL corner of grid (km, or deg. if latlongflag=1)"
        write(MR_global_info,*)"15.0    15.0  25.0              ! grid x-width, y-width, depth (km, or deg. if latlonflag=1)"
        write(MR_global_info,*)"0.5     0.5    2.0              ! DX, DY, Dz of grid cells  (km, or deg.)"
        write(MR_global_info,*)"1                               ! iHeightHandler"
        write(MR_global_info,*)"1980 5 18 15.5                  ! YYYY MM DD HH.H"
        write(MR_global_info,*)"1                               ! Output grid (1,2,3,4 for CompH,CompP,MetH,MetP)"
        write(MR_global_info,*)"5                               ! Output Var (1=H,2=Vx,3=Vy,4=Vz,5=T)"
        write(MR_global_info,*)"2                               ! Output Dimension (0,1,2,or 3)"
        write(MR_global_info,*)"0 0 1                           ! i,j,k of output grid with zeros for unused dims."
        write(MR_global_info,*)"5 25 2 2                        ! iwind iwindformat igrid iformat"
        write(MR_global_info,*)"1                               ! number of windfiles"
        write(MR_global_info,*)"NCEP"
        write(MR_global_info,*)"--------------------------------------------------"
        stop 1

      end subroutine write_usage

!##############################################################################
!##############################################################################
