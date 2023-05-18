      program gen_grib_index

      use MetReader
      use eccodes

      implicit none

      integer             :: nargs
      integer             :: status
      character (len=100) :: arg
      character(len=130)  :: grib_file

      INTERFACE
        subroutine MR_Set_Gen_Index_GRIB(grib_file)
          character(len=130),intent(in)  :: grib_file
        end subroutine MR_Set_Gen_Index_GRIB
      END INTERFACE

      nargs = command_argument_count()

      if (nargs.ne.1) then
        !do io=1,MR_nio;if(VB(io).le.verbosity_error)then
        !  write(errlog(io),*)"MR ERROR: no grib file given"
          write(6,*)"MR ERROR: no grib file given"
        !endif;enddo
        stop 1
      else
        call get_command_argument(1, arg, status)
        read(arg,*)grib_file
        call MR_Set_Gen_Index_GRIB(grib_file)
      endif

      end program
