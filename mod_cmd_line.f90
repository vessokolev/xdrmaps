      module mod_cmd_line
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
      ! the computations.
      !
      use iso_c_binding,only:C_SHORT,C_INT
      !
      implicit none

      public :: send_input_param_to_main
      private :: check_input
      private :: valid_keys 
      private :: valid_num_args 
      private :: key_in_keys
      private :: print_header
      private :: print_help
      private :: print_error_message_01
      private :: print_error_message_02
      private :: print_error_message_03
      private :: print_error_message_04
      private :: print_error_message_05
      private :: print_error_message_06

      contains

      subroutine send_input_param_to_main(intvl,&
                                          files,&
                                          factor,&
                                          extension)
      !
      ! The goal of this subroutine is to process the exit codes
      ! generated during the checking of the input parameters.
      !
      integer(C_INT),dimension(3),intent(out) :: intvl
      character(len=4096),dimension(4),intent(out) :: files
      character(len=1),intent(out) :: factor
      character(len=4),intent(out) :: extension
      !
      ! Local variables:
      !
      integer(C_SHORT) :: exit_code
      character(len=4096) :: int_str
      !
      ! Print the header. That shows the named of the program, the
      ! authors, the copyright, etc.
      !
      call print_header()
      !
      ! Check the command line input
      !
      call check_input(exit_code,&
                       intvl,&
                       int_str,&
                       files,&
                       factor,&
                       extension)
      !
      ! We presume all parameters are provided correctly by the user and
      ! we are going to do some additional processing only in case the
      ! exit_code value is not 0
      !
      if (exit_code.ne.0) then
         !
         ! Print the help every time an error is detected:
         !
         call print_help()
         !
         if (exit_code.eq.1) then
            !
            ! The trajectory file does not exist or it is not
            ! accessible.
            !
            call print_error_message_01(files(1))
         end if
         !
         if (exit_code.eq.2) then
            !
            ! The HDF5 with the atomic masses does not exist or it is
            ! not accessible.
            !
            call print_error_message_02(files(2))
         end if
         !
         if (exit_code.eq.3) then
            !
            ! The HDF5 with the presence of certain or all atoms does
            ! not exist or it is not accessible.
            !
            call print_error_message_03(files(3))
         end if
         !
         if (exit_code.eq.4) then
            !
            ! The HDF5 with the grid nodes cannot be created.
            !
            call print_error_message_04(files(4))
         end if
         !
         if (exit_code.eq.5) then
            !
            ! The number of nodes in each direction is not properly set.
            !
            call print_error_message_05(int_str)
         end if
         !
         if (exit_code.eq.6) then
            !
            ! The flag that defines which property of the atoms should
            ! be scanned is not properly set.
            !
            call print_error_message_06()
         end if
         !
         if (exit_code.eq.8) then
            !
            ! The flag indicates that either velocity of force
            ! distribution is requested on XTC trajectory instead based
            ! on TRR.
            !
            call print_error_message_07()
         end if
         !
         if (exit_code.eq.9) then
            !
            ! The flag indicates that either velocity of force
            ! distribution is requested on XTC trajectory instead based
            ! on TRR.
            !
            call print_error_message_08()
         end if
         !
         ! Terminate the program properly by supplying exit status to
         ! the invoking shell:
         !
         call exit(1)
      end if


      end subroutine send_input_param_to_main


      subroutine check_input(exit_code,&
                             intvl,&
                             int_str,&
                             files,&
                             factor,&
                             extension)
      !
      ! This subroutine checks the input parameters received when the
      ! program was invoked.
      !
      integer(C_SHORT),intent(out) :: exit_code
      integer(C_INT),dimension(3),intent(out) :: intvl
      character(len=4096),intent(out) :: int_str
      character(len=4096),dimension(4),intent(out) :: files
      character(len=1),intent(out) :: factor
      character(len=4),intent(out) :: extension
      !
      ! Local variables:
      !
      character(len=4096) :: tmp
      character(len=2),dimension(6) :: keys
      integer(C_INT) :: i
      integer(C_SHORT) :: num_args
      integer(C_INT) :: iostat_
      logical :: res
      !
      keys=(/'-f','-m','-p','-o','-i','-s'/)
      !
      num_args = command_argument_count()
      !
      if (valid_keys(num_args,6,keys)) then

         do i=1,num_args
               if (.not.mod(i,2).eq.0) then
                  call get_command_argument(i,tmp)
                  if (trim(adjustl(tmp)).eq.'-f') then
                     call get_command_argument(i+1,files(1))
                     files(1)=adjustl(files(1))
                     inquire(file=files(1), exist=res)
                     if (.not.res) then
                        exit_code=1
                        exit
                     end if
                  end if
                  if (trim(adjustl(tmp)).eq.'-m') then
                     call get_command_argument(i+1,files(2))
                     files(2)=adjustl(files(2))
                     inquire(file=files(2), exist=res)
                     if (.not.res) then
                        exit_code=2
                        exit
                     end if
                  end if
                  if (trim(adjustl(tmp)).eq.'-p') then
                     call get_command_argument(i+1,files(3))
                     files(3)=adjustl(files(3))
                     inquire(file=files(3), exist=res)
                     if (.not.res) then
                        exit_code=3
                        exit
                     end if
                  end if
                  if (trim(adjustl(tmp)).eq.'-o') then
                     call get_command_argument(i+1,files(4))
                     files(4)=adjustl(files(4))
                     open(9283,file=trim(adjustl(files(4))),&
                        iostat=iostat_,status='unknown')
                     if (iostat_.gt.0) then
                        exit_code=4
                        exit
                     else
                        close(9283,status='delete')
                     end if
                  end if
                  if (trim(adjustl(tmp)).eq.'-i') then
                     call get_command_argument(i+1,int_str)
                     read(int_str,*,iostat=iostat_) intvl(1),&
                                                    intvl(2),&
                                                    intvl(3)
                     if (iostat_.gt.0) then
                        exit_code=5
                        exit
                     end if
                  end if
                  if (trim(adjustl(tmp)).eq.'-s') then
                     call get_command_argument(i+1,tmp)
                     if (len(trim(adjustl(tmp))).eq.1) then
                        factor=trim(adjustl(tmp))
                        if (.not.((factor.eq.'M').or.&
                                  (factor.eq.'V').or.&
                                  (factor.eq.'F'))) then
                           exit_code=6
                           exit
                        end if
                     else
                        exit_code=6
                        exit
                     end if
                  end if
               end if
         end do
      else
         exit_code=7
      end if
      !
      if (exit_code.eq.0) then
         !
         ! Get the file extension of the trajectory file. Accept only
         ! one of ".xtc" or ".trr":
         !
         i=len(trim(files(1)))
         extension=files(1)(i-3:i)
         if ((extension.ne.'.xtc').and.(extension.ne.'.trr')) then
            exit_code=8
         end if
         !
         ! XTC files could be used only if the mass is a subject of
         ! investigation:
         !
         if (factor.eq.'V'.or.factor.eq.'F') then
            i=len(trim(files(1)))
            if (files(1)(i-3:i).eq.'.xtc') then
               exit_code=9
            end if
         end if
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
      print *,'material to the DOI 10.1021/acs.jpcb.5b06566'
      print *
      print *,'This program uses a source code from the following'
      print *,'external projects:'
      print *
      print *,'(1) xdrfile-1.1.4 : Copyright (c) 2009-2014, Erik'
      print *,'Lindahl & David van der Spoel All rights reserved'
      print *
      print *,'(2) xdrfort : XDR Fortran Interface with Wrappers'
      print *,'2014 (c) James W. Barnett <jbarnet4@tulane.edu>'
      print *,'https://github.com/wesbarnett/'
      print *,'Modified by Kai-Min Tu (2014) https://github.com/kmtu/'
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
      print *,trim(adjustl(tmp))//' -f trajectory.trr '//&
                                  ' -m mass.h5 '//&
                                  ' -p presence.h5 '//&
                                  ' -o output.h5 '//&
                                  ' -i "0 0 0" '//&
                                  ' -s M'
      print *
      print *,' where:'
      print *
      print *,'-f (trajectory.trr)'
      print *,'                    TRR trajectory file to analyze'
      print *,'-m (mass.h5)'
      print *,'                    HDF5 data set supplying the '
      print *,'                    atomic mass of the atoms involved'
      print *,'                    in the grid statistics analysis'
      print *,'-p (presence.h5)'
      print *,'                    HDF5 data set which specifies which'
      print *,'                    of the atoms in the topology should'
      print *,'                    participate in the grid statistics'
      print *,'-o (output.h5)'
      print *,'                    HDF5 data set containing the'
      print *,'                    collected grid statistics'
      print *,'-i ("0 0 0")'
      print *,'                    a vector of three integers which'
      print *,'                    specifies how many nodes should be'
      print *,'                    defined in each size of the box -'
      print *,'                    "N_x N_y N_z" (if any of the vector'
      print *,'                    "components is zero or negative, the'
      print *,'                    the program will optimize the number'
      print *,'                    of the nodes in each size of the'
      print *,'                    box)'
      print *,'-s (M)'
      print *,'                    a symbol which specifies which of'
      print *,'                    the atomic properties to investigate'
      print *,'                    (supply here M for the atomic mass,'
      print *,'                    V - for analyzing atomic velocity,'
      print *,'                    F - for analyzing the force upon the'
      print *,'                    atoms, note that to investigate the'
      print *,'                    velocity and the force they HAVE TO'
      print *,'                    be stored in the trajectory)'
      print *
      !
      end subroutine print_help


      subroutine print_error_message_01(traj_file)
      !
      ! This is the message that is going to be displayed in case the
      ! trajectory file does no exist or cannot be accessibe.
      !
      character(len=*),intent(in) :: traj_file
      !
      print *
      print *,'FATAL ERROR: The selected trajectory file:'
      print *
      print *,trim(adjustl(traj_file))
      print *
      print *,'does not exist or cannot be accessed!'
      print *
      !
      end subroutine print_error_message_01


      subroutine print_error_message_02(mass_h5_file)
      !
      ! This is the message that is going to be displayed in case the
      ! atomic mass HDF5 file does no exist or cannot be accessibe.
      !
      character(len=*),intent(in) :: mass_h5_file
      !
      print *
      print *,'FATAL ERROR: The selected HDF5 file that is '//&
              'supposed to supply the atomic massess:'
      print *
      print *,trim(adjustl(mass_h5_file))
      print *
      print *,'does not exist or cannot be accessed!'
      print *
      !
      end subroutine print_error_message_02


      subroutine print_error_message_03(presence_h5_file)
      !
      ! This is the message that is going to be displayed in case the
      ! presence HDF5 file does no exist or cannot be accessibe.
      !
      character(len=*),intent(in) :: presence_h5_file
      !
      print *
      print *,'FATAL ERROR: The selected HDF5 file that is '//&
              'supposed to supply the presence of the atoms in the '//&
              'grid statistics:'
      print *
      print *,trim(adjustl(presence_h5_file))
      print *
      print *,'does not exist or cannot be accessed!'
      print *
      !
      end subroutine print_error_message_03


      subroutine print_error_message_04(output_h5_file)
      !
      ! This is the message that is going to be displayed in case the
      ! presence HDF5 file does no exist or cannot be accessibe.
      !
      character(len=*),intent(in) :: output_h5_file
      !
      print *
      print *,'FATAL ERROR: The selected HDF5 output file which is '//&
              'supposed to contain the nodes of the grid:'
      print *
      print *,trim(adjustl(output_h5_file))
      print *
      print *,'cannot be created on the file system! Check if the '//&
              'target directory where the file should be created in'//&
              'to exists and the proper permissions are set on it. '//&
              'If the folder does not exist, create it first and '//&
              'then run the program again.'
      print *
      !
      end subroutine print_error_message_04


      subroutine print_error_message_05(int_str)
      !
      ! This is the message that is going to be displayed in case the
      ! number of nodes in each direction defined by the unit vectors in
      ! the box are not supplied.
      !
      character(len=4096),intent(in) :: int_str
      !
      print *
      print *,'FATAL ERROR: The amount of grid nodes in each '//&
              'direction of the box (where the directions are '//&
              'defined by the box unit vectors), supplied as:'
      print *
      print *,trim(adjustl(int_str))
      print *
      print *,'does not match the requirements. Some of the numbers '//&
              'specified are not recognized integers! Note that'//&
              'when specify the grid dimensions in each direction, '//&
              'that input must look like -i "ngrid_x ngrid_y '//&
              'ngrid_z" where ngrid_x, ngrid_y, and ngrid_z are '//&
              'integer numbers. If any of them is zero or negative '//&
              'the program will start an optimization to find these '//&
              'three numbers based on the best distance between the '//&
              'grid nodes, which later will result in map images '//&
              'with very high resolution.'
      print *
      !
      end subroutine print_error_message_05


      subroutine print_error_message_06()
      !
      ! This is the message that is going to be displayed in case the
      ! atomic property for investigation is not properly set by the
      ! user.
      !
      print *
      print *,'FATAL ERROR: The user requested an investigation of '//&
              'atomic property which cannot be processed by '//&
              'this program. It currently supports the '//&
              'investigation of the (i) the '//&
              'atomic mass, (ii) the velocity, (iii) the force, '//&
              'through the switches (specify only one of them):'
      print *
      print *,'-s M (this investigates the atomic mass)'
      print *,'-s V (this investigates the atomic velocity)'
      print *,'-s F (this investigates the force upon atoms)'
      print *
      !
      end subroutine print_error_message_06


      subroutine print_error_message_07()
      !
      ! This subroutine rises warning when the file extension of the
      ! trajectory file is not one of ".xtc" or ".trr".
      !
      print *
      print *,'FATAL ERROR: The file extension of the trajectory '//&
              'file is not one of ".xtc" and ".trr". Please, '//&
              'specify a correct GROMACS trajectory file after -f!'
      print *
      !
      end subroutine print_error_message_07


      subroutine print_error_message_08()
      !
      ! This subroutine explains why it is not possible to get the
      ! valocity or the forces from XTC trajectory file. That message
      ! must be displayed if the user is requested an investigation of
      ! the velocity and force but supplied XTC input trajectory file.
      !
      print *
      print *,'FATAL ERROR: The distribution of the atomic '//&
              'properties like velocity and force acting upon '//&
              'and atom, cannot be investigated based '//&
              'on XTC trajectories (files with extension *.xtc). '//&
              'That kind of investigation is possible only based '//&
              'on TRR trajectories (files with extension *.trr)! '//&
              'Run program again and provide TRR trajectory file '//&
              'after -f!'
      print *
      !
      end subroutine print_error_message_08




      end module mod_cmd_line
