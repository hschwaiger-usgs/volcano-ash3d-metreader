! MR_ASCII_check calculates the difference between two ascii grid files

      program MR_ASCII_check

      ! This module requires Fortran 2003 or later
      use iso_fortran_env, only : &
         input_unit,output_unit,error_unit

      use MetReader,       only : &
         MR_nio,VB,outlog,errlog,verbosity_error,verbosity_production

      implicit none

      integer, parameter :: ip             =  8 ! Internal precision
      integer, parameter :: fid_ascii2din  = 20 
      real(kind=ip), parameter :: EPS_TINY   = 1.0e-12_ip      ! Very small number

      integer           :: nargs
      integer           :: stat
      character(len=80) :: linebuffer080

      character(len=80) :: file_1,file_2
      logical :: IsThere1, IsThere2

      integer :: nx_1,nx_2
      integer :: ny_1,ny_2
      real(kind=ip) :: xll_1,xll_2
      real(kind=ip) :: yll_1,yll_2
      real(kind=ip) :: dx_1,dx_2
      real(kind=ip) :: dy_1,dy_2
      real(kind=ip) :: Fill_1,Fill_2
      real(kind=ip),dimension(:,:),allocatable :: XY_1, XY_2
      integer :: i,j, io
      real(kind=ip) :: err

      real(kind=ip) :: tmp_ip
      real(kind=ip) :: L2_toterror
      real(kind=ip) :: L2_tol = 1.0e-3
      integer           :: iostatus
      character(len=120):: iomessage

      MR_nio = 1  ! Turn off logging by setting output streams to stdout/stderr only

      nargs = command_argument_count()
      if (nargs.eq.0) then
        ! If no command-line arguments are given, then prompt user
        ! interactively for the two file names and the L2 tolerance
        do io=1,MR_nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter name of the first ESRI ASCII file:'
        endif;enddo
        read(input_unit,*) file_1
        do io=1,MR_nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter name of the second ESRI ASCII file:'
        endif;enddo
        read(input_unit,*) file_2
        do io=1,MR_nio;if(VB(io).le.verbosity_production)then
          write(outlog(io),*)'Enter the L2 tolerance (default is 1.0e-7):'
        endif;enddo
        read(input_unit,*) L2_tol
      elseif (nargs.eq.1.or.nargs.gt.3) then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)'ERROR: Too few command-line arguments.'
          write(errlog(io),*)'  Usage: Ash3d_ASCII_check file1 file2 (tol.)'
        endif;enddo
        stop 1
      else
        call get_command_argument(1, linebuffer080, status=stat)
        if(stat.gt.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 1'
          endif;enddo
          stop 1
        elseif (stat.lt.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Argument 1 has been truncated.'
            write(errlog(io),*)'       File name length is limited to 80 char.'
          endif;enddo
          stop 1
        endif
        file_1=trim(adjustl(linebuffer080))
        inquire( file=adjustl(trim(file_1)), exist=IsThere1 )

        call get_command_argument(2, linebuffer080, status=stat)
        if(stat.gt.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Could not parse argument 2'
          endif;enddo
          stop 1
        elseif (stat.lt.0)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Argument 2 has been truncated.'
            write(errlog(io),*)'       File name length is limited to 80 char.'
          endif;enddo
          stop 1
        endif
        file_2=trim(adjustl(linebuffer080))
        inquire( file=adjustl(trim(file_2)), exist=IsThere2 )

        if (.not.IsThere1.and..not.IsThere2)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Neither input files could be found'
          endif;enddo
          stop 1
        elseif (.not.IsThere1)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Input file 1 could not be found'
          endif;enddo
          stop 1
        elseif (.not.IsThere2)then
          do io=1,MR_nio;if(VB(io).le.verbosity_error)then
            write(errlog(io),*)'ERROR: Input file 2 could not be found'
          endif;enddo
          stop 1
        endif

        if (nargs.eq.3)then
          call get_command_argument(3, linebuffer080, status=stat)
          read(linebuffer080,*)tmp_ip
          if(stat.eq.0)then
            L2_tol = tmp_ip
          else
            do io=1,MR_nio;if(VB(io).le.verbosity_error)then
              write(errlog(io),*)'WARNING: argument 3 not read correctly.'
              write(errlog(io),*)'         Using default value.'
            endif;enddo
            stop 1
          endif
        endif
      endif

      ! Open and read ESRII ASCII file 1
      open(unit=fid_ascii2din,file=trim(adjustl(file_1)), status='old',action='read',err=2500)
      read(fid_ascii2din,3000,err=2600,iostat=iostatus,iomsg=iomessage) nx_1        ! read header values
      if(iostatus.ne.0)then
        ! We might have an empty file
        ! Issue warning and return
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) 'Error reading file ',trim(adjustl(file_1))
          write(errlog(io),*) 'Check for zero-length file.'
        endif;enddo
        stop 1
      endif
      read(fid_ascii2din,3001,err=2600,iostat=iostatus,iomsg=iomessage) ny_1
      allocate(XY_1(nx_1,ny_1))
      read(fid_ascii2din,3002,err=2600,iostat=iostatus,iomsg=iomessage) xll_1
      read(fid_ascii2din,3003,err=2600,iostat=iostatus,iomsg=iomessage) yll_1
      read(fid_ascii2din,3004,err=2600,iostat=iostatus,iomsg=iomessage) dx_1,dy_1
      read(fid_ascii2din,3005,err=2600,iostat=iostatus,iomsg=iomessage) Fill_1
      do j=ny_1,1,-1
        read(fid_ascii2din,3006,err=2600,iostat=iostatus,iomsg=iomessage) (XY_1(i,j), i=1,nx_1)
        read(fid_ascii2din,*,iostat=iostatus,iomsg=iomessage)
      enddo
      close(fid_ascii2din)

      ! Open and read ESRII ASCII file 2
      open(unit=fid_ascii2din,file=trim(adjustl(file_2)), status='old',action='read',err=2500)
      read(fid_ascii2din,3000,err=2600,iostat=iostatus,iomsg=iomessage) nx_2        ! read header values
      if(iostatus.ne.0)then
        ! We might have an empty file
        ! Issue warning and return
        do io=1,2;if(VB(io).le.verbosity_error)then
          write(errlog(io),*) 'Error reading file ',trim(adjustl(file_2))
          write(errlog(io),*) 'Check for zero-length file.'
        endif;enddo
        stop 1
      endif
      read(fid_ascii2din,3001,err=2600,iostat=iostatus,iomsg=iomessage) ny_2
      allocate(XY_2(nx_2,ny_2))
      read(fid_ascii2din,3002,err=2600,iostat=iostatus,iomsg=iomessage) xll_2
      read(fid_ascii2din,3003,err=2600,iostat=iostatus,iomsg=iomessage) yll_2
      read(fid_ascii2din,3004,err=2600,iostat=iostatus,iomsg=iomessage) dx_2,dy_2
      read(fid_ascii2din,3005,err=2600,iostat=iostatus,iomsg=iomessage) Fill_2
      do j=ny_2,1,-1
        read(fid_ascii2din,3006,err=2600,iostat=iostatus,iomsg=iomessage) (XY_2(i,j), i=1,nx_2)
        read(fid_ascii2din,*,iostat=iostatus,iomsg=iomessage)
      enddo
      close(fid_ascii2din)

!     format statements
3000  format(6x,i5)
3001  format(6x,i5)
3002  format(10x,f15.3)
3003  format(10x,f15.3)
3004  format(10x,2f15.3)
3005  format(13x,a6)
3006  format(10f18.6)

      if(nx_1.ne.nx_2)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"FAIL : nx differs"
        endif;enddo
        stop 1
      endif
      if(ny_1.ne.ny_2)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"FAIL : ny differs"
        endif;enddo
        stop 1
      endif
      if(abs(dx_1-dx_2).gt.EPS_TINY)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"FAIL : dx differs"
        endif;enddo
        stop 1
      endif
      if(abs(dy_1-dy_2).gt.EPS_TINY)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"FAIL : dy differs"
        endif;enddo
        stop 1
      endif
      if(abs(xll_1-xll_2).gt.EPS_TINY)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"FAIL : xll differs"
        endif;enddo
        stop 1
      endif
      if(abs(yll_1-yll_2).gt.EPS_TINY)then
        do io=1,MR_nio;if(VB(io).le.verbosity_error)then
          write(errlog(io),*)"FAIL : yll differs"
        endif;enddo
        stop 1
      endif

      L2_toterror = 0.0_ip
      do i=1,nx_1
        do j=1,ny_1
          err=abs(XY_1(i,j)-XY_2(i,j))
          L2_toterror = L2_toterror + err*err
        enddo
      enddo
      L2_toterror = sqrt(L2_toterror)
      L2_toterror = L2_toterror/nx_1/ny_1

      do io=1,MR_nio;if(VB(io).le.verbosity_error)then
        if(abs(L2_toterror).lt.L2_tol)then
          write(outlog(io),*)'PASS : ',L2_toterror
        else
          write(outlog(io),*)'FAIL : ',L2_toterror
        endif
      endif;enddo

      stop 0

2500  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error opening ASCII file. Program stopped'
      endif;enddo
      stop 1

2600  do io=1,2;if(VB(io).le.verbosity_error)then
        write(errlog(io),*) 'Error reading from ASCII file.'
      endif;enddo
      stop 1

      end program MR_ASCII_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
