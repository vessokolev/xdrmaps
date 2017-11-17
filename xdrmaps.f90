program xdrmaps
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
! Version: 2017111600'
! Author: Veselin Kolev <vesso.kolev@gmail.com>'
! License: GPLv2'
!
use iso_c_binding,only:C_INT,C_FLOAT
use hdf5
use mod_cmd_line
use mod_maps
!
implicit none
!
character(len=4096),dimension(4) :: files
character(len=4) :: extension
character(len=1) :: flag
integer(C_INT) :: i
integer(C_INT) :: code
integer(C_INT) :: num_frames
integer(C_INT) :: num_atoms
integer(C_INT),dimension(3) :: num_intervals
integer(C_INT),dimension(3) :: num_nodes
integer(C_INT),dimension(:,:,:),allocatable :: grid_f
real(C_FLOAT),dimension(3,3) :: box_average
real(C_FLOAT),dimension(3) :: step_average
real(C_FLOAT),dimension(3) :: corr_step
real(C_FLOAT),dimension(:),allocatable :: mass
real(C_FLOAT) :: step
real(C_FLOAT) :: time_s
real(C_FLOAT) :: time_f
real(C_FLOAT),dimension(:,:,:),allocatable :: grid
logical,dimension(:),allocatable :: presence
!
step=0.02_C_FLOAT
!
! Record the time of the beginning of the process of parsing the
! trajectory:
!
call cpu_time(time_s)
!
! Parse the command line parameters, check them and then supply
! them formatted to this scope of the program:
!
call send_input_param_to_main(num_nodes,files,flag,extension)
!
! Analyze the trajectory. Extract from there the average box size,
! the number of the atoms, and the number of the frames:
!
if (extension.eq.'.xtc') then
   !
   call analyze_ave_box_size_xtc(files(1),box_average,num_atoms,num_frames)
   !
else
   !
   call analyze_ave_box_size_trr(files(1),box_average,num_atoms,num_frames)
   !
end if
!
! Display the extracted information:
!
print *
print *,'Successfully finished parsing the trajectory. So far:'
print *,'Number of frames passed:',num_frames
print *,'Number of atoms detected:',num_atoms
print *,'Average box size:',(/(box_average(i,i),i=1,3)/)
print *
!
! Check if any of the requested number of nodes is less then or
! equal zero. If that is true start scanning the trajectory in
! order to get the optimum number of nodes, with respec to the
! paramater step as possible upper limit of the distance between
! two nodes in one direction.
!
if ((num_nodes(1).gt.0).and.(num_nodes(2).gt.0).and.(num_nodes(3).gt.0)) then
   !
   ! Get the number of interfvals each frame should be divided
   ! into:
   !
   num_intervals=num_nodes
   !
   ! 
   !
   if (any(num_intervals.lt.10)) then
      print *
      print *,'GENERAL ERROR: To few nodes. Rise up the number'
      print *,'of nodes.'
      print *
      call exit(1)
   end if
else
   !
   ! Get the number of interfvals each frame should be divided
   ! into:
   !
   call get_num_intervals_per_frame(box_average,step,num_intervals,corr_step)
   !
   num_nodes=num_intervals
   !
end if
!
! Compute the average step of division of the box volume as ratio
! between the average box size and the number of intervals:
!
step_average=(/(box_average(i,i),i=1,3)/)/num_intervals
!
! Read the presence file:
!
call read_atom_presence(files(3),'presence',presence)
!
! If the size of the presence array is not equal to the number of
! the atoms detected, that is a general error and the program
! cannot continue:
if (size(presence).ne.num_atoms) then
   print *
   print *,'GENERAL ERROR: The number of entries in the '//&
           'presence array is not equal to the number of '//&
           'the atoms in the frames of the trajectory'
   print *
   call exit(1)
end if
!
! Read the file with the atomic masses:
!
call read_atomic_mass(files(2),'mass',mass)
!
! If the size of the atomic mass array is not equal to the number
! of the atoms detected, that is a general error and the program
! cannot continue:
if (size(mass).ne.num_atoms) then
   print *
   print *,'GENERAL ERROR: The number of entries in the '//&
           'atomic mass array is not equal to the number of '//&
           'the atoms in the frames of the trajectory'
   print *
   call exit(1)
end if
!
! Create the grids:
!
print *,'Allocating the grid arrays ...'
allocate(grid(num_nodes(1),num_nodes(2),num_nodes(3)))
allocate(grid_f(num_nodes(1),num_nodes(2),num_nodes(3)))
print *,'[OK]'
!
! Starting the grid analysis:
!
if (flag.eq.'M') then
   !
   ! Processing the atomic mass:
   !
   if (extension.eq.'.xtc') then
      !
      call get_mass_grid_xtc(files(1),num_atoms,num_intervals,presence,mass,grid,grid_f)
      !
   else
      !
      call get_mass_grid_trr(files(1),num_atoms,num_intervals,presence,mass,grid,grid_f)
      !
   end if
   !
   ! Save the map:
   !
   code=1
   !
   call save_grid_as_hdf5((/(box_average(i,i),i=1,3)/),step_average,num_frames,num_nodes,&
                          num_atoms,files(4),code,grid,grid_f)
end if
!
if (flag.eq.'V') then
   !
   ! Processing the atomic velocities:
   !
   call get_vf_grid(files(1),num_atoms,num_intervals,presence,.true.,grid,grid_f)
   !
   ! Save the map:
   !
   code=2
   !
   call save_grid_as_hdf5((/(box_average(i,i),i=1,3)/),step_average,num_frames,num_nodes,&
                          num_atoms,files(4),code,grid,grid_f)
   !
end if
!
if (flag.eq.'F') then
   !
   ! Processing the force acting upon the atoms:
   !
   call get_vf_grid(files(1),num_atoms,num_intervals,presence,.false.,grid,grid_f)
   !
   ! Save the map:
   !
   code=3
   !
   call save_grid_as_hdf5((/(box_average(i,i),i=1,3)/),step_average,num_frames,&
                          num_nodes,num_atoms,files(4),code,grid,grid_f)
   !
end if   
!
! Get the current CPU time:
!
call cpu_time(time_f)
!
! Print out the total CPU time used for performing the
! derivation of the maps.
!
print *
print *,'Total CPU time for processing the trajectory:',&
        time_f-time_s,'seconds.'
print *
!
! Exit the program supplying a status code. Note that this is not
! a default behaviour of the FORTRAN programs.
!
call exit(0)
!
end program xdrmaps
