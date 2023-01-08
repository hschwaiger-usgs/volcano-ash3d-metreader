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
        write(0,*)"MR ERROR: no grib file given"
        stop 1
      else
        call get_command_argument(1, arg, status)
        read(arg,*)grib_file
        call MR_Set_Gen_Index_GRIB(grib_file)
      endif

      end program
