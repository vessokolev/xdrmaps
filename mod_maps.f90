module mod_maps
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
! Version: 2017111700'
! Author: Veselin Kolev <vesso.kolev@gmail.com>'
! License: GPLv2'
!
! This module processes the grid statistics and creates
! the maps.
!
use iso_c_binding,only:C_INT,C_FLOAT
use hdf5
use xdr, only: trrfile,xtcfile

implicit none

public :: analyze_ave_box_size_xtc
public :: analyze_ave_box_size_trr
public :: get_mass_grid_xtc
public :: get_mass_grid_trr
public :: get_num_intervals_per_frame
public :: get_vf_grid
public :: read_atom_presence
public :: read_atomic_mass
public :: save_grid_as_hdf5

private :: bubble_sort
private :: get_indx
private :: integer_to_string
private :: norm02

contains


subroutine bubble_sort(array_s,array_indx,array)
!
! This subroutine implements Bubble Sort method upon integer array
! in order to sort its elements in ascendant order. For more
! details regarding the algorithm visit:
! https://en.wikipedia.org/wiki/Bubble_sort
!
! In keeps the indexes of the original array to support checking
! the permutations occured during the sorting.
!
! Arguments:
!
! array_s - the size of the 1D array named "array".
!
!           NOTE: If outside the subroutine the array has size N
!                 but only first M elements are actual (the others
!                 might be zeroes), then set array_s=M.
!
! array - the 1D array of integers which elements are going to be
!         sorted.
!
! Interface variables:
!
integer(C_INT),intent(in) :: array_s
integer(C_INT),dimension(array_s),intent(inout) :: array_indx
real(C_FLOAT),dimension(array_s),intent(out) :: array
!
! Local variables:
!
integer(C_INT) :: i
real(C_FLOAT) :: dummy
integer(C_INT) :: dummy_indx
logical :: flag
!
! First, construct the index array (use the simpliest possible
! array constructor here):
!
array_indx=(/(i,i=1,array_s)/)

flag=.True.
do while(flag)
   flag=.False.
   do i=1,array_s-1
      if (array(i).gt.array(i+1)) then
         dummy=array(i+1)
         dummy_indx=array_indx(i+1)
         array(i+1)=array(i)
         array_indx(i+1)=array_indx(i)
         array(i)=dummy
         array_indx(i)=dummy_indx
         flag=.True.
      end if
   end do
end do

end subroutine bubble_sort


subroutine get_num_intervals_per_frame(box,min_step,num_intervals,corr_step)
!
! This subroutine compute the number of intervals in each
! direction computed according equation E2 in Supporting
! Information of DOI 10.1021/acs.jpcb.5b06566
!
! Interface variables:
!
real(C_FLOAT),dimension(3,3),intent(in) :: box
real(C_FLOAT),intent(in) :: min_step
integer(C_INT),dimension(3),intent(out) :: num_intervals
real(C_FLOAT),dimension(3),intent(out) :: corr_step
!
! Local variables:
!
integer(C_INT) :: i
integer(C_INT),dimension(3) :: box_diagonal_indx
real(C_FLOAT),dimension(3) :: box_diagonal

box_diagonal=(/(box(i,i),i=1,3)/)

call bubble_sort(3,box_diagonal_indx,box_diagonal)

num_intervals(box_diagonal_indx(1))=&
   int(floor(box_diagonal(1)/min_step))

num_intervals(box_diagonal_indx(2))=&
   nint(num_intervals(box_diagonal_indx(1))*&
   box_diagonal(2)/box_diagonal(1))

num_intervals(box_diagonal_indx(3))=&
   nint(num_intervals(box_diagonal_indx(1))*&
   box_diagonal(3)/box_diagonal(1))

do i=1,3
   corr_step(i)=box(i,i)/num_intervals(i)
end do

end subroutine get_num_intervals_per_frame


subroutine read_atom_presence(h5_file,dsetname,presence)
!
! This subroutine reads the array that is responsible to take in
! account which atom (by its serial number) should be include in
! the process of the grid statistics collection. The index number
! of the atom in the array+1 is equal to the actual serial number
! of the atom.
!
! Interface variables:
!
character(len=*),intent(in) :: h5_file
character(len=*),intent(in) :: dsetname
logical,allocatable :: presence(:)
!
! Local variables:
!
integer(HID_T) :: file_id
integer(HID_T) :: dset_id
integer(HID_T) :: dataspace
integer(HSIZE_T),allocatable :: presence_dims(:)
integer(HSIZE_T),allocatable :: max_presence_dims(:)
integer(C_INT) :: error
integer(C_INT) :: rankr
integer(C_INT) :: i
integer(C_INT),allocatable :: presence_int(:)
!
! Initialize the HDF5 interface:
!
call h5open_f(error)
!
if (error.ne.0) then
   print *
   print *,'FATAL ERROR: Cannot open the HDF5 interface!'
   print *,'Check if the application is compiled to a proper '//&
           'and locally existing HDF5 distribution!'
   print *
   call exit(1)
end if
!
! Open the HDF5 file:
!
call h5fopen_f(trim(adjustl(h5_file)),H5F_ACC_RDWR_F,file_id,error)
!
if (error.ne.0) then
   print *
   print *,'FATAL ERROR: Cannot open the HDF5 file:'
   print *
   print *,trim(adjustl(h5_file))
   print *
   print *,'Check if the file exists and the appropriate '//&
           'permissions are set on it.'
   print *
   call exit(1)
end if
!
! The name of the dataset is passed as an argument:
!
call h5dopen_f(file_id,trim(adjustl(dsetname)),dset_id,error)
!
if (error.ne.0) then
   print *
   print *,'FATAL ERROR: The requested HDF5 dataset:'
   print *
   print *,trim(adjustl(dsetname))
   print *
   print *,'is not presented in the file:'
   print *
   print *,trim(adjustl(h5_file))
   print *
   print *,'Provide the correct HDF5 file which contains this '//&
           'dataset!'
   print *
   call exit(1)
end if
!
! Get the dataspace id:
!
call h5dget_space_f(dset_id, dataspace, error)
!
if (error.ne.0) then
   print *
   print *,'FATAL ERROR: It is not possible to get the '//&
           'dataspace object of the data space whih name :'
   print *
   print *,trim(adjustl(dsetname))
   print *
   print *,'Check the consistence of the HDF5 file:'
   print *
   print *,trim(adjustl(h5_file))
   print *
   print *
   call exit(1)
end if
!
! Get the rank of the dataspace:
!
call h5sget_simple_extent_ndims_f(dataspace, rankr, error)
!
if (error.ne.0) then
   print *
   print *,'FATAL ERROR: It is not possible to get the '//&
           'rank of the data space whih name :'
   print *
   print *,trim(adjustl(dsetname))
   print *
   print *,'Check the consistence of the HDF5 file:'
   print *
   print *,trim(adjustl(h5_file))
   print *
   print *
   call exit(1)
end if
!
! Get the dimensions of the dataspace, but first allocate
! the corresponding dynamic arrays:
!
allocate(presence_dims(rankr))
allocate(max_presence_dims(rankr))
!
call h5sget_simple_extent_dims_f(dataspace,presence_dims,&
                                 max_presence_dims,error)
!
! No need to check the error here!!!
!
! Read the dataset. Before that allocate the dynamic array to
! transfer the dataset data into:
!
allocate(presence(presence_dims(1)))
allocate(presence_int(presence_dims(1)))
!
call h5dread_f(dset_id,H5T_NATIVE_INTEGER,presence_int,&
               presence_dims,error)
if (error.ne.0) then
   print *
   print *,'FATAL ERROR: Cannot read the content of the '//&
           'dataset whih name :'
   print *
   print *,trim(adjustl(dsetname))
   print *
   print *,'Check the consistence of the HDF5 file:'
   print *
   print *,trim(adjustl(h5_file))
   print *
   print *
   call exit(1)
end if
!
! To save memory assign the presence to a boolean array and
! dispose the integer one:
!
do i=1,size(presence_int)
   if (presence_int(i)==0) then
      presence(i)=.false.
   else
      presence(i)=.true.
   end if
end do
!
deallocate(presence_int)
!
! We do not need anymore the arrays to store the dimensions of the
! dataset space. Deallocate them:
!
deallocate(presence_dims)
deallocate(max_presence_dims)
!
! Close the HDF5 interface
!
call h5dclose_f(dset_id,error)
!
if (error.ne.0) then
   print *
   print *,'Cannot close the HDF5 interface after finishing '//&
           'with reading the data set with name:'
   print *
   print *,trim(adjustl(dsetname))
   print *
   print *,'Possible memory problem or problem with the '//&
           'the HDF5 implementation this application is '//&
           'linked to!'
   print *
   call exit(1)
end if
!
end subroutine read_atom_presence


subroutine read_atomic_mass(h5_file,dsetname, mass)
!
! This subroutine reads the array that contains the atomic mass of
! each atom presented in the topology. The index number of the
! element in the array + 1 is equal to the actual atom serial
! number of the atom in the topology.
!
! Interface variables:
!
character(len=*),intent(in) :: h5_file
character(len=*),intent(in) :: dsetname
real(C_FLOAT),dimension(:),allocatable :: mass
!
! Local variables:
!
integer(HID_T) :: file_id
integer(HID_T) :: dset_id
integer(HID_T) :: dataspace
integer(HSIZE_T),dimension(:),allocatable :: mass_dims
integer(HSIZE_T),dimension(:),allocatable :: max_mass_dims
integer(C_INT) :: error
integer(C_INT) :: rankr
!
! Initialize the HDF5 interface:
!
call h5open_f(error)
!
if (error.ne.0) then
   print *
   print *,'FATAL ERROR: Cannot open the HDF5 interface!'
   print *,'Check if the application is compiled to a proper '//&
           'and locally existing HDF5 distribution!'
   print *
   call exit(1)
end if
!
! Open the HDF5 file:
!
call h5fopen_f(trim(adjustl(h5_file)),H5F_ACC_RDWR_F,file_id,error)
!
if (error.ne.0) then
   print *
   print *,'FATAL ERROR: Cannot open the HDF5 file:'
   print *
   print *,trim(adjustl(h5_file))
   print *
   print *,'Check if the file exists and the appropriate '//&
           'permissions are set on it.'
   print *
   call exit(1)
end if
!
! The name of the dataset is passed as an argument:
!
call h5dopen_f(file_id,trim(adjustl(dsetname)),dset_id,error)
!
if (error.ne.0) then
   print *
   print *,'FATAL ERROR: The requested HDF5 dataset:'
   print *
   print *,trim(adjustl(dsetname))
   print *
   print *,'is not presented in the file:'
   print *
   print *,trim(adjustl(h5_file))
   print *
   print *,'Provide the correct HDF5 file which contains this '//&
           'dataset!'
   print *
   call exit(1)
end if
!
! Get the dataspace id:
!
call h5dget_space_f(dset_id, dataspace, error)
!
if (error.ne.0) then
   print *
   print *,'FATAL ERROR: It is not possible to get the '//&
           'dataspace object of the data space whih name :'
   print *
   print *,trim(adjustl(dsetname))
   print *
   print *,'Check the consistence of the HDF5 file:'
   print *
   print *,trim(adjustl(h5_file))
   print *
   print *
   call exit(1)
end if
!
! Get the rank of the dataspace:
!
call h5sget_simple_extent_ndims_f(dataspace, rankr, error)
!
if (error.ne.0) then
   print *
   print *,'FATAL ERROR: It is not possible to get the '//&
           'rank of the data space whih name :'
   print *
   print *,trim(adjustl(dsetname))
   print *
   print *,'Check the consistence of the HDF5 file:'
   print *
   print *,trim(adjustl(h5_file))
   print *
   print *
   call exit(1)
end if
!
! Get the dimensions of the dataspace, but first allocate
! the corresponding dynamic arrays:
!
allocate(mass_dims(rankr))
allocate(max_mass_dims(rankr))
!
call h5sget_simple_extent_dims_f(dataspace,mass_dims,max_mass_dims,error)
!
! No need to check the error here!!!
!
! Read the dataset. Before that allocate the dynamic array to
! transfer the dataset data into:
!
allocate(mass(mass_dims(1)))
!
call h5dread_f(dset_id,H5T_NATIVE_REAL,mass,mass_dims,error)
!
if (error.ne.0) then
   print *
   print *,'FATAL ERROR: Cannot read the content of the '//&
           'dataset whih name :'
   print *
   print *,trim(adjustl(dsetname))
   print *
   print *,'Check the consistence of the HDF5 file:'
   print *
   print *,trim(adjustl(h5_file))
   print *
   print *
   call exit(1)
end if
!
! We do not need anymore the arrays to store the dimensions of the
! dataset space. Deallocate them:
!
deallocate(mass_dims)
deallocate(max_mass_dims)
!
! Close the HDF5 interface
!
call h5dclose_f(dset_id,error)
!
if (error.ne.0) then
   print *
   print *,'Cannot close the HDF5 interface after finishing '//&
           'with reading the data set with name:'
   print *
   print *,trim(adjustl(dsetname))
   print *
   print *,'Possible memory problem or problem with the '//&
           'the HDF5 implementation this application is '//&
           'linked to!'
   print *
   call exit(1)
end if
!
end subroutine read_atomic_mass


subroutine save_grid_as_hdf5(box_average,step_average,num_frames,num_nodes,&
                             num_atoms,h5_file,code,grid,grid_f)
!
! Interface variables:
!
real(C_FLOAT),dimension(3),intent(in) :: box_average
real(C_FLOAT),dimension(3),intent(in) :: step_average
integer(C_INT),intent(in) :: num_frames
integer(C_INT),dimension(3),intent(in) :: num_nodes
integer(C_INT),intent(in) :: num_atoms
character(len=*),intent(in) :: h5_file
integer(C_INT),intent(in) :: code
real(C_FLOAT),dimension(:,:,:),intent(in) :: grid
integer(C_INT),dimension(:,:,:),intent(in) :: grid_f
!
! Local variables:
!
integer(HID_T) :: file_id
integer(HID_T) :: dset_id
integer(HID_T) :: dspace_id
integer(HSIZE_T),dimension(3) :: grid_s
integer(HSIZE_T),dimension(1) :: box_s
integer(HSIZE_T),dimension(1) :: num_frames_s
integer(C_INT) :: error
!
! First, get the dimensions of the arrays:
!
box_s=(/3/)
num_frames_s=(/1/)
grid_s=shape(grid)
!
! Open the HDF5 file containing the grid array.
!
call h5open_f(error)
!
! Create the HDF5 file.
!
call h5fcreate_f(trim(adjustl(h5_file)),H5F_ACC_TRUNC_F,file_id,error)
!
! Create a dataset for storing there the box dimensions
!
call h5screate_simple_f(1,box_s,dspace_id,error)
!
call h5dcreate_f(file_id,'box',H5T_NATIVE_REAL,dspace_id,dset_id,error)
!
call h5dwrite_f(dset_id,H5T_NATIVE_REAL,box_average,box_s,error)
!
call h5dclose_f(dset_id,error)
!
call h5sclose_f(dspace_id,error)
!
! Create a dataset for storing there the steps used to define the
! grid intervals
!
call h5screate_simple_f(1,box_s,dspace_id,error)
!
call h5dcreate_f(file_id,'step',H5T_NATIVE_REAL,dspace_id,dset_id,error)
!
call h5dwrite_f(dset_id,H5T_NATIVE_REAL,step_average,box_s,error)
!
call h5dclose_f(dset_id,error)
!
call h5sclose_f(dspace_id,error)
!
! Create a dataset for storing there the number of rames detected
! when parsing the trajectory:
!
call h5screate_simple_f(1,num_frames_s,dspace_id,error)
!
call h5dcreate_f(file_id,'num_frames',H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
!
call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,(/num_frames/),num_frames_s,error)
!
call h5dclose_f(dset_id,error)
!
call h5sclose_f(dspace_id,error)
!
! Create a dataset for storing there the code number of the grid
! content: 1 - mass, 2 - velocity, 3 - force.
!
call h5screate_simple_f(1,num_frames_s,dspace_id,error)
!
call h5dcreate_f(file_id,'code',H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
!
call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,(/code/),num_frames_s,error)
!
call h5dclose_f(dset_id,error)
!
call h5sclose_f(dspace_id,error)
!
! Create a dataset for storing there the number of the atoms
! in each and every frame of the trajectory:
!
call h5screate_simple_f(1,num_frames_s,dspace_id,error)
!
call h5dcreate_f(file_id,'num_atoms',H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
!
call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,(/num_atoms/),num_frames_s,error)
!
call h5dclose_f(dset_id,error)
!
call h5sclose_f(dspace_id,error)
!
! Create a dataset for storing there the number of nodes in each
! box size:
!
call h5screate_simple_f(1,box_s,dspace_id,error)
!
call h5dcreate_f(file_id,'num_nodes',H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
!
call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,num_nodes,box_s,error)
!
call h5dclose_f(dset_id,error)
!
call h5sclose_f(dspace_id,error)
!
! Create a dataset for storing the grid:
!
call h5screate_simple_f(3,grid_s,dspace_id,error)
!
call h5dcreate_f(file_id,'grid',H5T_NATIVE_REAL,dspace_id,dset_id,error)
!
call h5dwrite_f(dset_id,H5T_NATIVE_REAL,grid,grid_s,error)
!
call h5dclose_f(dset_id,error)
!
call h5sclose_f(dspace_id,error)
!
! Create a dataset for storing the grid frequency (that shows
! which nodes are more or less updated):
!
call h5screate_simple_f(3,grid_s,dspace_id,error)
!
call h5dcreate_f(file_id,'freq',H5T_NATIVE_INTEGER,dspace_id,dset_id,error)
!
call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,grid_f,grid_s,error)
!
call h5dclose_f(dset_id,error)
!
call h5sclose_f(dspace_id,error)
!
! Close the HDF5 file:
!
call h5fclose_f(file_id, error)
!
! Dispose the HDF5 interface.
!
call h5close_f(error)
!
end subroutine save_grid_as_hdf5


subroutine get_indx(box,corr_step,coord,indx)
!
! This subroutine get the index of the node in the grid with
! respect to the coordinates of the atom. It is an implementation
! of equation E6 in Supporting Information of
! DOI 10.1021/acs.jpcb.5b06566
!
! Note that the indxes given in eqiation E6 are are in C and
! Python style (that means the first index in the array in 0).
! Here, in the FORTRAN code the first index is 1! See how the
! variable indx is initializied bellow!
!
! Interface variables:
!
real(C_FLOAT),dimension(3,3),intent(in) :: box
real(C_FLOAT),dimension(3),intent(in) :: corr_step
real(C_FLOAT),dimension(3),intent(in) :: coord
integer(C_INT),dimension(3),intent(out) :: indx
!
! Local variables:
!
integer(C_INT) :: i
real(C_FLOAT),dimension(3) :: corr
!
indx=1
corr(:)=coord
!
! VERY IMPORTANT!!! Check if the coordinate of the atom
! are within the box. If not, correct them with respect to
! the PBC and box size!
!
do i=1,3
   !
   if (coord(i).lt.0.0_C_FLOAT) then
      corr(i)=coord(i)+box(i,i)
   end if
   !
   if (coord(i).ge.box(i,i)) then
      corr(i)=coord(i)-box(i,i)
   end if
   !
   indx(i)=int(corr(i)/corr_step(i))+1
   !
end do
!
end subroutine get_indx


subroutine analyze_ave_box_size_trr(trr_file,box_average,num_atoms,num_frames)
!
! This subroutine should be implemented only for trajectory
! derived from NPT simulations (or other simulations where the box
! volume oscilates). All the computations performed by this
! subroutine code are dedicated to solve the equation E1 in the
! Supporting Information of DOI: 10.1021/acs.jpcb.5b06566
!
! Interface variables:
!
character(len=*),intent(in) :: trr_file
real(C_FLOAT),dimension(3,3),intent(out) :: box_average
integer(C_INT),intent(out) :: num_atoms
integer(C_INT),intent(out) :: num_frames
!
! Local variables:
!
type(trrfile) :: trr
logical :: flag
character(len=50) :: tmp_1
integer(C_INT) :: tmp_1_len
!
! Initialize the output variables:
!
box_average=0.0_C_FLOAT
num_frames=0_C_INT
!
! Initialize the local logical variable flag:
!
flag=.true.
!
! Open the trr trajectory file:
!
call trr%init(trim(adjustl(trr_file)))
!
! Start reading the frames:
!
call trr%read
!
! And to that to the last frame of the trajectory (the variables
! trr%STAT gets <>0 when the end of the trajectory is reached):
!
do while (trr%STAT==0)
   !
   ! Use the high-rank operator for accumulating the sizes of the
   ! box for getting the average:
   !
   box_average=box_average+trr%box
   num_frames=num_frames+1
   !
   call integer_to_string(num_frames,tmp_1_len,tmp_1)
   !
   if (flag) then
      flag=.false.
      num_atoms=trr%NATOMS
   else
      if (trr%NATOMS.ne.num_atoms) then
         print *
         print *,'FATAL ERROR: The frame number '//&
                 tmp_1(1:tmp_1_len)//' at time',&
                 trr%time,'ps contains different '//&
                 'number of atoms then in the other frames'//&
                 '! The trajectory is either impropertly '//&
                 'manipulated or broken!'
         print *
         call exit(1)
      end if
   end if
   !
   print *,'Successfully read frame No. '//tmp_1(1:tmp_1_len)//&
           ' at t=',trr%time,'ps'
   print *
   !
   ! Read next frame:
   !
   call trr%read
   !
end do
!
! Cut the access to the trajectory file:
!
call trr%close
!
! Do average the box size:
!
box_average=box_average/num_frames
!
end subroutine analyze_ave_box_size_trr


subroutine analyze_ave_box_size_xtc(xtc_file,box_average,num_atoms,num_frames)
!
! This subroutine should be implemented only for trajectory
! derived from NPT simulations (or other simulations where the box
! volume oscilates). All the computations performed by this
! subroutine code are dedicated to solve the equation E1 in the
! Supporting Information of DOI: 10.1021/acs.jpcb.5b06566
!
! Interface variables:
!
character(len=*),intent(in) :: xtc_file
real(C_FLOAT),dimension(3,3),intent(out) :: box_average
integer(C_INT),intent(out) :: num_atoms
integer(C_INT),intent(out) :: num_frames
!
! Local variables:
!
type(xtcfile) :: xtc
logical :: flag
character(len=50) :: tmp_1
integer(C_INT) :: tmp_1_len
!
! Initialize the output variables:
!
box_average=0.0_C_FLOAT
num_frames=0_C_INT
!
! Initialize the local logical variable flag:
!
flag=.true.
!
! Open the trr trajectory file:
!
call xtc%init(trim(adjustl(xtc_file)))
!
! Start reading the frames:
!
call xtc%read
!
! And to that to the last frame of the trajectory (the variables
! trr%STAT gets <>0 when the end of the trajectory is reached):
!
do while (xtc%STAT==0)
   !
   ! Use the high-rank operator for accumulating the sizes of the
   ! box for getting the average:
   !
   box_average=box_average+xtc%box
   num_frames=num_frames+1
   !
   call integer_to_string(num_frames,tmp_1_len,tmp_1)
   !
   if (flag) then
      flag=.false.
      num_atoms=xtc%NATOMS
   else
      if (xtc%NATOMS.ne.num_atoms) then
         print *
         print *,'FATAL ERROR: The frame number '//&
                  tmp_1(1:tmp_1_len)//' at time',&
                  xtc%time,'ps contains different '//&
                 'number of atoms then in the other frames'//&
                 '! The trajectory is either impropertly '//&
                'manipulated or broken!'
         print *
         call exit(1)
      end if
   end if
   !
   print *,'Successfully read frame No. '//tmp_1(1:tmp_1_len)//&
           ' at t=',xtc%time,'ps'
   print *
   !
   ! Read next frame:
   !
   call xtc%read
   !
end do
!
! Cut the access to the trajectory file:
!
call xtc%close
!
! Do average the box size:
!
box_average=box_average/num_frames
!
end subroutine analyze_ave_box_size_xtc


subroutine get_mass_grid_trr(trr_file,num_atoms,num_intervals,presence,&
                             mass,mass_grid,mass_grid_f)
!
! Investigates the atomic mass distribution alongside the
! simulation box. This is the TRR version of the subroutine.
!
!
! Interface variables:
!
character(len=*),intent(in) :: trr_file
integer(C_INT),intent(in) :: num_atoms
integer(C_INT),dimension(3),intent(in) :: num_intervals
logical,dimension(:),intent(in) :: presence
real(C_FLOAT),dimension(:),intent(in) :: mass
real(C_FLOAT),dimension(:,:,:),intent(out) :: mass_grid
integer(C_INT),dimension(:,:,:),intent(out) :: mass_grid_f
!
! Local variables:
!
type(trrfile) :: trr
integer(C_INT) :: i
integer(C_INT) :: j
integer(C_INT) :: tmp_1_len
integer(C_INT) :: num_frames
integer(C_INT),dimension(3) :: indx
integer(C_INT),dimension(3) :: s
real(C_FLOAT),dimension(3) :: corr_step
character(len=50) :: tmp_1
!
! Initialize:
!
num_frames=0
mass_grid=0.0_C_FLOAT
mass_grid_f=0
s=shape(mass_grid)
!
! Open the trr trajectory file:
!
call trr%init(trim(adjustl(trr_file)))
!
! Start reading the frames:
!
call trr%read
!
! Do the reading until reach the end of the trajectory file:
!
do while (trr%STAT==0)
   !
   num_frames=num_frames+1
   !
   call integer_to_string(num_frames,tmp_1_len,tmp_1)
   !
   ! Now the new interval size should be computed. It is because
   ! in NPT simulations the box size oscillates. This is an
   ! implementation of equation E5 in Supplimentary Information of
   ! DOI 10.1021/acs.jpcb.5b06566
   !
   corr_step=(/(trr%box(i,i)/num_intervals(i),i=1,3)/)
   !
   do i=1,num_atoms
      !
      ! Process the i-th atom only it that is requested in the presence
      ! array:
      !
      if (presence(i)) then
         !
         ! Translate the atom coordinates into node indexes:
         !
         call get_indx(trr%box,corr_step,trr%pos(:,i),indx)
         !
         ! Sometimes (due to a bug of problematic compilation) some atoms
         ! might appear outside the box. Correct the indexes of the nodes
         ! with respect to that unusual atom coordinates:
         !
         do j=1,3
            if (indx(j)<1) then
               indx(j)=indx(j)+s(j)
            end if
            if (indx(j)>s(j)) then
               indx(j)=indx(j)-s(j)
            end if
         end do
         !
         ! Update the mass grid:
         !
         mass_grid(indx(1),indx(2),indx(3))=&
            mass_grid(indx(1),indx(2),indx(3))+mass(i)
         !
         ! Update the frequency grid:
         !
         mass_grid_f(indx(1),indx(2),indx(3))=&
            mass_grid_f(indx(1),indx(2),indx(3))+1
         !
         ! NOTE: The frequency grid is not always required when populating
         !       the nodes of the mass grid. It is given here just for the
         !       sake of completion.
         !
      end if
   end do
   !
   print *,'Successfully processed frame No. '//&
           tmp_1(1:tmp_1_len)//' at t=',trr%time,'ps'
   print *
   !
   ! Read next frame:
   !
   call trr%read
   !
end do
!
! Close the connection to the trajectory file:
!
call trr%close
!
end subroutine get_mass_grid_trr


subroutine get_mass_grid_xtc(xtc_file,num_atoms,num_intervals,presence,&
                             mass,mass_grid,mass_grid_f)
!
! Investigates the atomic mass distribution alongside the
! simulation box. This is the XTC version of the subroutine.
!
! Interface variables:
!
character(len=*),intent(in) :: xtc_file
integer(C_INT),intent(in) :: num_atoms
integer(C_INT),dimension(3),intent(in) :: num_intervals
logical,dimension(:),intent(in) :: presence
real(C_FLOAT),dimension(:),intent(in) :: mass
real(C_FLOAT),dimension(:,:,:),intent(out) :: mass_grid
integer(C_INT),dimension(:,:,:),intent(out) :: mass_grid_f
!
! Local variables:
!
type(xtcfile) :: xtc
integer(C_INT) :: i
integer(C_INT) :: j
integer(C_INT) :: tmp_1_len
integer(C_INT) :: num_frames
integer(C_INT),dimension(3) :: indx
integer(C_INT),dimension(3) :: s
real(C_FLOAT),dimension(3) :: corr_step
character(len=50) :: tmp_1
!
! Initialize:
!
num_frames=0
mass_grid=0.0_C_FLOAT
mass_grid_f=0
s=shape(mass_grid)
!
! Open the trr trajectory file:
!
call xtc%init(trim(adjustl(xtc_file)))
!
! Start reading the frames:
!
call xtc%read
!
! Do the reading until reach the end of the trajectory file:
!
do while (xtc%STAT==0)
   !
   num_frames=num_frames+1
   !
   call integer_to_string(num_frames,tmp_1_len,tmp_1)
   !
   ! Now the new interval size should be computed. It is because
   ! in NPT simulations the box size oscillates. This is an
   ! implementation of equation E5 in Supplimentary Information of
   ! DOI 10.1021/acs.jpcb.5b06566
   !
   corr_step=(/(xtc%box(i,i)/num_intervals(i),i=1,3)/)
   !
   do i=1,num_atoms
      !
      ! Process the i-th atom only it that is requested in the presence
      ! array:
      !
      if (presence(i)) then
         !
         ! Translate the atom coordinates into node indexes:
         !
         call get_indx(xtc%box,corr_step,xtc%pos(:,i),indx)
         !
         ! Sometimes (due to a bug of problematic compilation) some atoms
         ! might appear outside the box. Correct the indexes of the nodes
         ! with respect to that unusual atom coordinates:
         !
         do j=1,3
            if (indx(j)<1) then
               indx(j)=indx(j)+s(j)
            end if
            if (indx(j)>s(j)) then
               indx(j)=indx(j)-s(j)
            end if
         end do
         !
         ! Update the mass grid:
         !
         mass_grid(indx(1),indx(2),indx(3))=&
            mass_grid(indx(1),indx(2),indx(3))+mass(i)
         !
         ! Update the frequency grid:
         !
         mass_grid_f(indx(1),indx(2),indx(3))=&
            mass_grid_f(indx(1),indx(2),indx(3))+1
         !
         ! NOTE: The frequency grid is not always required when populating
         !       the nodes of the mass grid. It is given here just for the
         !       sake of completion.
         !
      end if
      !
   end do
   !
   print *,'Successfully processed frame No. '//&
           tmp_1(1:tmp_1_len)//' at t=',xtc%time,'ps'
   print *
   !
   ! Read next frame:
   !
   call xtc%read
   !
end do
!
! Close the connection to the trajectory file:
!
call xtc%close
!
end subroutine get_mass_grid_xtc


subroutine get_vf_grid(trr_file,num_atoms,num_intervals,presence,&
                       v_or_f,vf_grid,vf_grid_f)
!
! This subroutine scans the given trajectory and estimates the
! local distribution of the velocity of force.
!
! Interface variables:
!
character(len=*),intent(in) :: trr_file
integer(C_INT),intent(in) :: num_atoms
integer(C_INT),dimension(3),intent(in) :: num_intervals
logical,dimension(:),intent(in) :: presence
logical,intent(in) :: v_or_f
real(C_FLOAT),dimension(:,:,:),intent(out) :: vf_grid
integer(C_INT),dimension(:,:,:),intent(out) :: vf_grid_f
!
! Local variables:
!
type(trrfile) :: trr
integer(C_INT) :: i
integer(C_INT) :: j
integer(C_INT),dimension(3) :: indx
integer(C_INT) :: tmp_1_len
integer(C_INT) :: num_frames
integer(C_INT),dimension(3) :: s
real(C_FLOAT) :: tmp
real(C_FLOAT),dimension(3) :: corr_step
character(len=50) :: tmp_1
!
! Initialize:
!
num_frames=0
vf_grid=0.0_C_FLOAT
vf_grid_f=0
s=shape(vf_grid)
!
! Open the trr trajectory file:
!
call trr%init(trim(adjustl(trr_file)))
!
! Start reading the frames:
!
call trr%read
!
! Do the reading until reach the end of the trajectory file:
!
do while (trr%STAT==0)
   !
   num_frames=num_frames+1
   !
   call integer_to_string(num_frames,tmp_1_len,tmp_1)
   !
   ! Now the new interval size should be computed. It is because
   ! in NPT simulations the box size oscillates. This is an
   ! implementation of equation E5 in Supplimentary Information of
   ! DOI 10.1021/acs.jpcb.5b06566
   !
   corr_step=(/(trr%box(i,i)/num_intervals(i),i=1,3)/)
   !
   do i=1,num_atoms
      !
      ! Process the i-th atom only it that is requested in the presence
      ! array:
      !
      if (presence(i)) then
         !
         ! Translate the atom coordinates into node indexes:
         !
         call get_indx(trr%box,corr_step,trr%pos(:,i),indx)
         !
         ! Compute the amount to add to the grid nodes depending on which
         ! value is requested to be handled (force or velocity):
         !
         if (v_or_f) then
            tmp=norm02(trr%vel(:,i))
         else
            tmp=norm02(trr%force(:,i))
         end if
         !
         ! Sometimes (due to a bug of problematic compilation) some atoms
         ! might appear outside the box. Correct the indexes of the nodes
         ! with respect to that unusual atom coordinates:
         !
         do j=1,3
            if (indx(j)<1) then
               indx(j)=indx(j)+s(j)
            end if
            if (indx(j)>s(j)) then
               indx(j)=indx(j)-s(j)
            end if
         end do
         !
         ! Update the velocity/force grid:
         !
         vf_grid(indx(1),indx(2),indx(3))=&
            vf_grid(indx(1),indx(2),indx(3))+tmp
         !
         ! Update the frequency grid:
         !
         vf_grid_f(indx(1),indx(2),indx(3))=&
            vf_grid_f(indx(1),indx(2),indx(3))+1
      end if
   end do
   !
   print *,'Successfully processed frame No. '//&
           tmp_1(1:tmp_1_len)//' at t=',trr%time,'ps'
   print *
   !
   ! Read next frame:
   !
   call trr%read
   !
end do
!
! Close the connection to the trajectory file:
!
call trr%close
!
end subroutine get_vf_grid


subroutine integer_to_string(int_num,str_len,str)
!
! This subroutine prints integer number int_num into the string
! str and estimates the actual string length str_len.
!
! Interface variables:
!
integer(C_INT),intent(in) :: int_num
integer(C_INT),intent(out) :: str_len
character(len=50),intent(out) :: str

write(str,*) int_num

str=adjustl(str)

str_len=len(trim(str))

end subroutine integer_to_string


function norm02(vector) result(norm)
!
! This function computes the norm of a vector in 3D. Note that
! most of the compilers provide the function norm2 as an intrinsic
! one, so the goal of having norm02 implemented here is only to
! make it sure that the norm will be computed even if the compiler
! does not equipped with norm2.
!
! Interface variables:
!
real(C_FLOAT),dimension(3),intent(in) :: vector
real(C_FLOAT) :: norm
!
! Local variables:
!
integer(C_INT) :: i

norm=0.0_C_FLOAT

do i=1,3
   norm=norm+vector(i)**2
end do

norm=sqrt(norm)

end function norm02


end module mod_maps
