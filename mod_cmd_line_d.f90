      module mod_cmd_line_d
      !
      !########################################################
      !#                                                      #
      !#                        XDRMAPS                       #
      !#                                                      #
      !#    TOOL FOR PARSING TRR TRAJECTORIES AND CREATING    #
      !#     GRIDS CONTAINING THE SPATIAL DISTRIBUTION OF     #
      !#              SELECTED ATOMIC PROPERTIES              #
      !#                                                      #
      !########################################################
      !
      ! Version: 2017031400'
      ! Author: Veselin Kolev <vesso.kolev@gmail.com>'
      ! License: GPLv2'
      !
      ! This module parses the command line parameters before running
      ! dumping the maps.
      !
      use iso_c_binding,only:C_SHORT,C_INT
      !
      implicit none

      public :: check_input
      private :: valid_keys 
      private :: valid_num_args 
      private :: key_in_keys
      private :: print_header
      private :: print_help
      private :: check_if_string_is_int
      private :: print_error_message_01
      private :: print_error_message_02
      private :: print_error_message_03
      private :: print_error_message_04
      private :: print_error_message_05

      contains


      subroutine check_input(grid_file,&
                             dim_,&
                             step_,&
                             noise_red,&
                             path)
      !
      ! This subroutine checks the input parameters received when the
      ! program was invoked.
      !
      character(len=4096),intent(out) :: grid_file
      integer(C_INT),intent(out) :: dim_
      integer(C_INT),intent(out) :: step_
      logical,intent(out) :: noise_red
      character(len=4096),intent(out) :: path
      !
      ! Local variables:
      !
      character(len=4096) :: tmp
      character(len=4096) :: dim_string
      character(len=4096) :: step_str
      character(len=4096) :: noise_red_str
      character(len=2),dimension(5) :: keys
      integer(C_INT) :: i
      integer(C_INT) :: j
      integer(C_SHORT) :: num_args
      integer(C_INT) :: iostat_
      logical :: res
      !
      keys=(/'-g','-d','-s','-w','-o'/)
      !
      num_args = command_argument_count()
      !
      call print_header()
      !
      if (valid_keys(num_args,5,keys)) then

         do i=1,num_args
               if (.not.mod(i,2).eq.0) then
                  call get_command_argument(i,tmp)
                  if (trim(adjustl(tmp)).eq.'-g') then
                     call get_command_argument(i+1,grid_file)
                     grid_file=adjustl(grid_file)
                     inquire(file=trim(grid_file), exist=res)
                     if (.not.res) then
                        call print_help()
                        call print_error_message_01(grid_file)
                        call exit(1)
                     end if
                  end if
                  if (trim(adjustl(tmp)).eq.'-d') then
                     call get_command_argument(i+1,dim_string)
                     res=check_if_string_is_int(dim_string,dim_)
                     if (.not.res) then
                        call print_help()
                        call print_error_message_02(dim_string)
                        call exit(1)
                     end if
                  end if
                  if (trim(adjustl(tmp)).eq.'-s') then
                     call get_command_argument(i+1,step_str)
                     res=check_if_string_is_int(step_str,step_)
                     if (.not.res) then
                        call print_help()
                        call print_error_message_03(step_str)
                        call exit(1)
                     end if
                  end if
                  if (trim(adjustl(tmp)).eq.'-w') then
                     call get_command_argument(i+1,noise_red_str)
                     res=check_if_string_is_int(noise_red_str,j)
                     if (.not.res) then
                        call print_help()
                        call print_error_message_04(noise_red_str)
                        call exit(1)
                     else
                        if (j.eq.0) then
                           noise_red=.false.
                        else
                           if (j.eq.1) then
                              noise_red=.true.
                           else
                              call print_help()
                              call print_error_message_04(noise_red_str)
                              call exit(1)
                           end if
                        end if
                     end if
                  end if
                  if (trim(adjustl(tmp)).eq.'-o') then
                     call get_command_argument(i+1,path)
                     path=adjustl(path)
                     open(9283,file=trim(path)//&
                               '/check.3926865544354657',&
                               iostat=iostat_,status='unknown')
                     if (iostat_.gt.0) then
                        call print_help()
                        call print_error_message_05(path)
                        call exit(1)
                     else
                       close(9283,status='delete')
                     end if
                  end if
               end if
         end do
      else
         call print_help()
         print *,'FATAL ERROR: Wrong number of input parameters!'
         call exit(1)
      end if
      !
      end subroutine check_input


      function valid_keys(num_args,nkeys,keys) result(res)
      !
      ! This function validates each of the keys supplied through the
      ! command line.
      !
      ! Interface variables:
      !
      integer(C_SHORT) :: num_args
      integer(C_INT),intent(in) :: nkeys
      character(len=2),dimension(nkeys),intent(in) :: keys
      logical :: res
      !
      ! Local variables:
      !
      integer(C_INT) :: i
      integer(C_SHORT) :: counter
      character(len=4096) :: tmp
      !
      res=.false.
      counter=0
      !
      if (valid_num_args(num_args,nkeys)) then
         do i=1,num_args
            call get_command_argument(i,tmp)
            if (.not.mod(i,2).eq.0) then
               if (key_in_keys(tmp,nkeys,keys)) then
                  counter=counter+1
               end if
            end if
         end do
      end if
      !
      if ((counter.gt.0).and.(counter.eq.num_args/2)) then
         res=.true.
      end if
      !
      end function valid_keys


      function valid_num_args(num_args,nkeys) result(res)
      !
      ! Checks the number of pairs in the supplied through the command
      ! line parameters.
      !
      ! Interface variables:
      !
      integer(C_SHORT) :: num_args
      integer(C_INT),intent(in) :: nkeys
      !
      ! Local variables:
      !
      logical :: res
      !
      res=.false.
      !
      if (mod(num_args,2).eq.0) then
         res=.true.
         if ((num_args/2).ne.nkeys) then
            res=.false.
         end if
      end if
      !
      end function valid_num_args


      function key_in_keys(key,nkeys,keys) result(res)
      !
      ! This function checks if a given key string is among the key
      ! strings stored in the array keys.
      !
      ! Interface variables:
      !
      character(len=*),intent(in) :: key
      integer(C_INT),intent(in) :: nkeys
      character(len=2),dimension(nkeys),intent(in) :: keys
      logical :: res
      !
      ! Local variables:
      !
      integer(C_SHORT) :: i
      !
      res=.false.
      !
      if (len(trim(adjustl(key))).eq.2) then
         do i=1,nkeys
            if (trim(adjustl(key)).eq.keys(i)) then
               res=.true.
               exit
            end if
         end do
      end if
      !
      end function key_in_keys


      subroutine print_header()
      !
      ! This subroutine prints the program greeting message.
      !
      print *
      print *,'########################################################'
      print *,'#                                                      #'
      print *,'#                        XDRMAPS                       #'
      print *,'#                                                      #'
      print *,'#  TOOL FOR PARSING TRR/XTC TRAJECTORIES AND CREATING  #'
      print *,'#     GRIDS CONTAINING THE SPATIAL DISTRIBUTION OF     #'
      print *,'#              SELECTED ATOMIC PROPERTIES              #'
      print *,'#                                                      #'
      print *,'########################################################'
      print *
      print *,'VERSION: 2017031400'
      print *,'AUTHOR: Veselin Kolev <vesso.kolev@gmail.com>'
      print *,'LICENSE: GPLv2'
      print *
      print *,'The program implements the idea of the grid statistics'
      print *,'as described in details in the Supporting Information'
      print *,'material to the DOI 10.1021/acs.jpcb.5b06566. This tool'
      print *,'dumps the maps from the collected grid statistics into'
      print *,'a set of text files.'
      print *
      !
      end subroutine print_header


      subroutine print_help()
      !
      ! This subroutine only prints the help
      !
      ! Interface variables:
      !
      character(len=4096) :: tmp
      !
      call get_command_argument(0,tmp)
      !
      print *
      print *,'Invoke the program by supplying the parameters '//&
              'following this example:'
      print *
      print *,trim(adjustl(tmp))//' -g grid.h5 '//&
                                  ' -d 3 '//&
                                  ' -s 5 '//&
                                  ' -o /output/folder '
      print *
      print *,' where:'
      print *
      print *,'-g (grid.h5)'
      print *,'                    the HDF5 file containing the grid'
      print *,'-d (3)'
      print *,'                    the axis perpendicular to the maps'
      print *,'                    plane:'
      print *,'                    1 - x,'
      print *,'                    2 - y,'
      print *,'                    3 - z'
      print *,'                    in the grid statistics analysis'
      print *,'-s (5)'
      print *,'                    how many planes to unite in order to'
      print *,'                    construct each map - the number of'
      print *,'                    the planes is 2 times the number'
      print *,'                    given here + 1 (2*5+1) - that is the'
      print *,'                    width of the sliding window'
      print *,'-w (0)'
      print *,'                    if set to 1, perform a noise'
      print *,'                    reduction to make the maps easier to'
      print *,'                    read (set 0 here if the noise'
      print *,'                    reduction should be skip) - node'
      print *,'                    that the noise reduction changes'
      print *,'                    values of the investigated property'
      print *,'                    up so one must found a way how to'
      print *,'                    normalize the results after the'
      print *,'                    reduction'
      print *,'-o (/output/folder)'
      print *,'                    the folder where the maps will be'
      print *,'                    saved as text files'
      print *
      !
      end subroutine print_help


      subroutine print_error_message_01(grid_file)
      !
      ! This is the message that is going to be displayed in case the
      ! grid file does no exist or cannot be accessibe.
      !
      character(len=*),intent(in) :: grid_file
      !
      print *
      print *,'FATAL ERROR: The selected grid file:'
      print *
      print *,trim(adjustl(grid_file))
      print *
      print *,'does not exist or cannot be accessed!'
      print *
      !
      end subroutine print_error_message_01


      function check_if_string_is_int(param,int_) result(res)
      !
      ! Checks if the string "param" is an integer number.
      !
      ! Interface variables:
      !
      character(len=*),intent(in) :: param
      integer(C_INT),intent(out) :: int_
      logical :: res
      !
      ! Local variables:
      !
      integer(C_INT) :: i

      read(param,*,iostat=i) int_

      res=.true.

      if (i.ne.0) then
         res=.false.
      end if

      end function check_if_string_is_int


      subroutine print_error_message_02(dim_string)
      !
      ! This is the message that is going to be displayed in case the
      ! direction for getting the maps is not selected correctly.
      !
      character(len=*),intent(in) :: dim_string
      !
      print *
      print *,'FATAL ERROR: The requested direction for getting the '//&
              'maps:'
      print *
      print *,trim(adjustl(dim_string))
      print *
      print *,'is not one of 1, 2, 3!'
      print *
      !
      end subroutine print_error_message_02


      subroutine print_error_message_03(step)
      !
      ! This is the message that is going to be displayed in case the
      ! half-width of the sliding window is not selected propertly.
      !
      character(len=*),intent(in) :: step
      !
      print *
      print *,'FATAL ERROR: The selected number of planes '
      print *
      print *,trim(adjustl(step))
      print *
      print *,'required for getting the thickness of the sliding '//&
              'window is not the correct one. Note that this must '//&
              'an integer number and its value have to be lesser '//&
              'than the half of the size of the grid in the '//&
              'processed dimension.'
      print *
      !
      end subroutine print_error_message_03


      subroutine print_error_message_04(noise_red_str)
      !
      ! This is the message that is going to be displayed in case the
      ! half-width of the sliding window is not selected propertly.
      !
      character(len=*),intent(in) :: noise_red_str
      !
      print *
      print *,'FATAL ERROR: The selector supplied along with -w '
      print *
      print *,trim(adjustl(noise_red_str))
      print *
      print *,'does not specify whether to apply a noise reduction '//&
              'or not. Note that the possible values of that input '//&
              'parameter have to either 0 (skip the reduction) or '//&
              '1 (perform the reduction).'
      print *
      !
      end subroutine print_error_message_04


      subroutine print_error_message_05(path)
      !
      ! This is the message that is going to be displayed in case the
      ! selected folder for wtinting down the files with the maps does
      ! not exist or it is not accessible.
      !
      character(len=*),intent(in) :: path
      !
      print *
      print *,'FATAL ERROR: The selected output folder for storing '//&
              'the maps:'
      print *
      print *,trim(adjustl(path))
      print *
      print *,'does not exist or the permissions set upon it do not '//&
              'allow writing files there. If the file does not '//&
              'exist create it first and run this util again.'
      print *
      !
      end subroutine print_error_message_05


      end module mod_cmd_line_d
