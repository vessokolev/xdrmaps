      program test

      use iso_c_binding
      use mod_cmd_line_d

      implicit none

      character(len=4096) :: grid_file
      integer(C_INT) :: dim_
      integer(C_INT) :: step_
      character(len=4096) :: path

      call check_input(grid_file,&
                       dim_,&
                       step_,&
                       path)


      end program test
