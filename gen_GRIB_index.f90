!##############################################################################
!##############################################################################
!
! gen_grib_index
!
! This stand-alone program is a wrapper for the subroutine MR_Set_Gen_Index_GRIB.
! A grib2 filename as a single argument is expected. This program then generates
! a grib index file.
!
!##############################################################################

      program gen_grib_index

      use MetReader

      use eccodes

      use, intrinsic :: iso_fortran_env, only : &
         real32,real64,input_unit,output_unit,error_unit

      implicit none
      !implicit none (type, external)

      integer             :: nargs
      integer             :: iostatus
      character (len=100) :: arg
      character (len=130) :: grib_file

      INTERFACE
        subroutine MR_Set_Gen_Index_GRIB(grib_file)
          implicit none
          !implicit none (type, external)
          character(len=130),intent(in)  :: grib_file
        end subroutine MR_Set_Gen_Index_GRIB
      END INTERFACE

      nargs = command_argument_count()

      if (nargs /= 1) then
        !do io=1,MR_nio;if(MR_VB(io) <= verbosity_error)then
        !  write(errlog(io),*)"MR ERROR: no grib file given"
          write(input_unit,*)"MR ERROR: no grib file given"
        !endif;enddo
        stop 1
      else
        call get_command_argument(number=1, value=arg, status=iostatus)
        read(arg,*)grib_file
        call MR_Set_Gen_Index_GRIB(grib_file)
      endif

      end program gen_grib_index
