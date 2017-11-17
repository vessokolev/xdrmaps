module mod_maps_d
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
!
implicit none

public :: read_real_vector_set
public :: read_int_vector_set
public :: read_scalar_set
public :: read_grid
public :: read_freq_grid
public :: get_mass_slice
public :: get_vf_slice
public :: integer_to_string
public :: float_to_string
public :: write_map_file
public :: stencil
public :: get_slice_shape
public :: get_step_2D
public :: get_node_pbc


contains


subroutine read_real_vector_set(h5_file,dsetname,vect)
!
! This subroutine reads a real vector with 3 elements
! represented as HDF5 dataset.
!
! Interface variables:
!
character(len=*),intent(in) :: h5_file
character(len=*),intent(in) :: dsetname
real(C_FLOAT),dimension(3),intent(out) :: vect
!
! Local variables:
!
integer(HID_T) :: file_id
integer(HID_T) :: dset_id
integer(HID_T) :: dataspace
integer(HSIZE_T),dimension(3) :: dims
integer(HSIZE_T),dimension(3) :: max_dims
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
call h5dget_space_f(dset_id,dataspace,error)
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
call h5sget_simple_extent_ndims_f(dataspace,rankr,error)
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
! Get the dimensions of the dataspace:
!
call h5sget_simple_extent_dims_f(dataspace,dims,max_dims,error)
!
! No need to check the error here!!!
!
call h5dread_f(dset_id,H5T_NATIVE_REAL,vect,dims,error)
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
end subroutine read_real_vector_set


subroutine read_int_vector_set(h5_file,dsetname,vect)
!
! This subroutine reads an integer vector with 3 elements
! represented as HDF5 dataset.
!
! Interface variables:
!
character(len=*),intent(in) :: h5_file
character(len=*),intent(in) :: dsetname
integer(C_INT),dimension(3),intent(out) :: vect
!
! Local variables:
!
integer(HID_T) :: file_id
integer(HID_T) :: dset_id
integer(HID_T) :: dataspace
integer(HSIZE_T),dimension(3) :: dims
integer(HSIZE_T),dimension(3) :: max_dims
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
call h5dget_space_f(dset_id,dataspace,error)
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
call h5sget_simple_extent_ndims_f(dataspace,rankr,error)
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
! Get the dimensions of the dataspace:
!
call h5sget_simple_extent_dims_f(dataspace,dims,max_dims,error)
!
! No need to check the error here!!!
!
call h5dread_f(dset_id,H5T_NATIVE_INTEGER,vect,dims,error)
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
end subroutine read_int_vector_set


subroutine read_scalar_set(h5_file,dsetname,val)
!
! This subroutine reads a scalar represented as HDF5 dataset.
!
! Interface variables:
!
character(len=*),intent(in) :: h5_file
character(len=*),intent(in) :: dsetname
integer(C_INT),intent(out) :: val
!
! Local variables:
!
integer(HID_T) :: file_id
integer(HID_T) :: dset_id
integer(HID_T) :: dataspace
integer(HSIZE_T),dimension(1) :: dims
integer(HSIZE_T),dimension(1) :: max_dims
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
call h5dget_space_f(dset_id,dataspace,error)
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
call h5sget_simple_extent_ndims_f(dataspace,rankr,error)
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
! Get the dimensions of the dataspace:
!
call h5sget_simple_extent_dims_f(dataspace,dims,max_dims,error)
!
! No need to check the error here!!!
!
call h5dread_f(dset_id,H5T_NATIVE_INTEGER,val,dims,error)
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

end subroutine read_scalar_set


subroutine read_grid(h5_file,dsetname,grid)
!
! This subroutine reads the array that contains 3D array or real
! type.
!
! Interface variables:
!
character(len=*),intent(in) :: h5_file
character(len=*),intent(in) :: dsetname
real(C_FLOAT),dimension(:,:,:),allocatable,intent(out) :: grid
!
! Local variables:
!
integer(HID_T) :: file_id
integer(HID_T) :: dset_id
integer(HID_T) :: dataspace
integer(HSIZE_T),dimension(:),allocatable :: dims
integer(HSIZE_T),dimension(:),allocatable :: max_dims
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
call h5dget_space_f(dset_id,dataspace,error)
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
call h5sget_simple_extent_ndims_f(dataspace,rankr,error)
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
allocate(dims(rankr))
allocate(max_dims(rankr))
!
call h5sget_simple_extent_dims_f(dataspace,dims,max_dims,error)
!
! No need to check the error here!!!
!
! Read the dataset. Before that allocate the dynamic array to
! transfer the dataset data into:
!
allocate(grid(dims(1),dims(2),dims(3)))
!
call h5dread_f(dset_id,H5T_NATIVE_REAL,grid,dims,error)
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
deallocate(dims)
deallocate(max_dims)
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

end subroutine read_grid


subroutine read_freq_grid(h5_file,dsetname,grid)
!
! This subroutine reads the array that contains 3D array or
! integer type.
!
! Interface variables:
!
character(len=*),intent(in) :: h5_file
character(len=*),intent(in) :: dsetname
integer(C_INT),dimension(:,:,:),allocatable,intent(out) :: grid
!
! Local variables:
!
integer(HID_T) :: file_id
integer(HID_T) :: dset_id
integer(HID_T) :: dataspace
integer(HSIZE_T),dimension(:),allocatable :: dims
integer(HSIZE_T),dimension(:),allocatable :: max_dims
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
call h5dget_space_f(dset_id,dataspace,error)
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
call h5sget_simple_extent_ndims_f(dataspace,rankr,error)
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
allocate(dims(rankr))
allocate(max_dims(rankr))
!
call h5sget_simple_extent_dims_f(dataspace,dims,max_dims,error)
!
! No need to check the error here!!!
!
! Read the dataset. Before that allocate the dynamic array to
! transfer the dataset data into:
!
allocate(grid(dims(1),dims(2),dims(3)))
!
call h5dread_f(dset_id,H5T_NATIVE_INTEGER,grid,dims,error)
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
deallocate(dims)
deallocate(max_dims)
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
end subroutine read_freq_grid


subroutine get_mass_slice(grid_s,grid,sel_dim,current_node,slice_width,&
                          slice)
!
! Interface variables:
!
integer(C_INT),dimension(3),intent(in) :: grid_s
real(C_FLOAT),dimension(:,:,:),intent(in) :: grid
integer(C_INT),intent(in) :: sel_dim
integer(C_INT),intent(in) :: current_node
integer(C_INT),intent(in) :: slice_width
real(C_FLOAT),dimension(:,:),intent(out) :: slice
!
! Local variables:
!
integer(C_INT) :: i
integer(C_INT) :: tmp

slice=0.0_C_FLOAT

do i=current_node-slice_width,current_node+slice_width
   tmp=get_node_pbc(grid_s,sel_dim,i)
   if (sel_dim.eq.1) then
      slice=slice+grid(tmp,:,:)
   end if
   if (sel_dim.eq.2) then
      slice=slice+grid(:,tmp,:)
   end if
   if (sel_dim.eq.3) then
      slice=slice+grid(:,:,tmp)
   end if
end do

end subroutine get_mass_slice


subroutine get_vf_slice(grid_s,grid,freq,sel_dim,current_node,&
                        slice_width,slice_s,slice)
!
! Interface variables:
!
integer(C_INT),dimension(3),intent(in) :: grid_s
real(C_FLOAT),dimension(:,:,:),intent(in) :: grid
integer(C_INT),dimension(:,:,:),intent(in) :: freq
integer(C_INT),intent(in) :: sel_dim
integer(C_INT),intent(in) :: current_node
integer(C_INT),intent(in) :: slice_width
integer(C_INT),dimension(2),intent(in) :: slice_s
real(C_FLOAT),dimension(:,:),intent(out) :: slice
!
! Local variables:
!
integer(C_INT) :: i
integer(C_INT) :: j
integer(C_INT) :: tmp
integer(C_INT),dimension(slice_s(1),slice_s(2)) :: slice_f

slice=0.0_C_FLOAT
slice_f=0

do i=current_node-slice_width,current_node+slice_width
   tmp=get_node_pbc(grid_s,sel_dim,i)
   if (sel_dim.eq.1) then
      slice=slice+grid(tmp,:,:)
      slice_f=slice_f+freq(tmp,:,:)
   end if
   if (sel_dim.eq.2) then
      slice=slice+grid(:,tmp,:)
      slice_f=slice_f+freq(:,tmp,:)
   end if
   if (sel_dim.eq.3) then
      slice=slice+grid(:,:,tmp)
      slice_f=slice_f+freq(:,:,tmp)
   end if
end do

do i=1,slice_s(1)
   do j=1,slice_s(2)
      if (slice_f(i,j).ne.0) then
         slice(i,j)=slice(i,j)/slice_f(i,j)
      end if
   end do
end do

end subroutine get_vf_slice


function get_slice_shape(array_s,sel_dim) result(res)
!
! Computes the sizes of the square shape which collects a slice of
! the grid statistics to represent it as a map.
!
! Interface variables:
!
integer(C_INT),dimension(3),intent(in) :: array_s
integer(C_INT),intent(in) :: sel_dim
!
! This is the value returned by the function:
!
integer(C_INT),dimension(2) :: res
!
! Local variables:
!
integer(C_INT) :: i
integer(C_INT) :: j

res=(/0,0/)

j=1

do i=1,3
   if (i.ne.sel_dim) then
      res(j)=array_s(j)
      j=j+1
   end if
end do

end function get_slice_shape


function get_step_2D(step,sel_dim) result(res)
!
! Constructs an 2D array containing the step in each of the map
! dimensions.
!
! Interface variables:
!
real(C_FLOAT),dimension(3),intent(in) :: step
integer(C_INT),intent(in) :: sel_dim
!
! This is the value returned by the function:
!
real(C_FLOAT),dimension(2) :: res
!
! Local variables:
!
integer(C_INT) :: i
integer(C_INT) :: j

res=(/0,0/)

j=1

do i=1,3
   if (i.ne.sel_dim) then
      res(j)=step(j)
      j=j+1
   end if
end do

end function get_step_2D


function get_node_pbc(array_s,sel_dim,node) result(res)
!
! Computes the indexes of the nodes in the grid with respect to
! the Periodic Boundary Conditions (PBC) so the slice cannot go
! outside the box volume.
!
! Interface variables:
!
integer(C_INT),dimension(3),intent(in) :: array_s
integer(C_INT),intent(in) :: sel_dim
integer(C_INT),intent(in) :: node
!
! This is the value returned by the function:
!
integer(C_INT) :: res

res=node

if (node.lt.1) then
   res=node+array_s(sel_dim)-1
else
   if (node.ge.array_s(sel_dim)) then
      res=node-array_s(sel_dim)+1
   end if
end if

end function get_node_pbc


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


subroutine float_to_string(float_num,fmt_str,str_len,str)
!
! This subroutine prints the floating point number float_num into
! the string str and estimates the actual string length str_len.
! The format of the floating point number is supplied as fmt_str.
!
! Interface variables:
!
real(C_FLOAT),intent(in) :: float_num
character(len=*),intent(in) :: fmt_str
integer(C_INT),intent(out) :: str_len
character(len=50),intent(out) :: str

write(str,fmt=trim(adjustl(fmt_str))) float_num

str=adjustl(str)

str_len=len(trim(str))

end subroutine float_to_string


subroutine write_map_file(output_folder,map_file_prefix,map_number,&
                          slice_s,slice)
!
! Interface variables:
!
character(len=*),intent(in) :: output_folder
character(len=*),intent(in) :: map_file_prefix
integer(C_INT),intent(in) :: map_number
integer(C_INT),dimension(2),intent(in) :: slice_s
real(C_FLOAT),dimension(:,:),intent(in) :: slice
!
! Local variables:
!
character(len=4096) :: tmp
character(len=4096) :: fmt_
character(len=4096) :: file_name
integer(C_INT) :: i
integer(C_INT) :: j

write(tmp,fmt='(I10.10)') map_number

file_name=trim(adjustl(output_folder))//'/'//&
          trim(adjustl(map_file_prefix))//&
          trim(adjustl(tmp))//'.map'

write(tmp,fmt='(I10)') slice_s(2)

fmt_='('//trim(adjustl(tmp))//'E14.7)'

open(23432,file=trim(adjustl(file_name)),status='unknown')

do i=1,slice_s(1)
   write(23432,fmt=fmt_)(slice(i,j),j=1,slice_s(2))
end do

close(23432)

end subroutine write_map_file


subroutine stencil(slice_s,slice,stencil_nodes,w,slice_filtered)
!
! This subroutine removes the noise from the slice of distribution
! that represents 2D map, by implementing the stencil method
! defined with equation E17 in the Supporting Information of
! DOI: 10.1021/acs.jpcb.5b06566 
!
integer(C_INT),dimension(2),intent(in) :: slice_s
real(C_FLOAT),dimension(:,:),intent(in) :: slice
integer(C_INT),dimension(2,12),intent(in) :: stencil_nodes
real(C_FLOAT),dimension(4),intent(in) :: w
real(C_FLOAT),dimension(:,:),intent(out) :: slice_filtered
!
! Local variables:
!
integer(C_INT) :: i
integer(C_INT) :: j
integer(C_INT) :: p
integer(C_INT) :: q
integer(C_INT),dimension(2,12) :: node_c
!
! Initialize the slice_filtered (it is supposed to keep the
! noise filtered slice):
!
slice_filtered=0.0_C_FLOAT
!
! Do the noise reduction:
!
do i=1,slice_s(1)
   do j=1,slice_s(2)
      !
      ! Initialize the stencil nodes indexes (as if there is no
      ! PBC):
      !
      node_c(1,:)=stencil_nodes(1,:)+i
      node_c(2,:)=stencil_nodes(2,:)+j
      !
      ! Correct the node indexes of the stencil with respect to
      ! PBC:
      !
      do p=1,12
         do q=1,2
            if (node_c(q,p).gt.slice_s(q)) then
               node_c(q,p)=node_c(q,p)-slice_s(q)
            end if
            if (node_c(q,p).le.0) then
               node_c(q,p)=node_c(q,p)+slice_s(q)
            end if
         end do
      end do
      !
      ! Assign the nodes values to the filtered version of the
      ! slice:
      !
      slice_filtered(i,j)=w(1)*slice(i,j)+&
                          w(2)*(slice(node_c(1,1),node_c(2,1))+&
                                slice(node_c(1,2),node_c(2,2))+&
                                slice(node_c(1,3),node_c(2,3))+&
                                slice(node_c(1,4),node_c(2,4)))+&
                          w(3)*(slice(node_c(1,5),node_c(2,5))+&
                                slice(node_c(1,6),node_c(2,6))+&
                                slice(node_c(1,7),node_c(2,7))+&
                                slice(node_c(1,8),node_c(2,8)))+&
                          w(4)*(slice(node_c(1,9),node_c(2,9))+&
                                slice(node_c(1,10),node_c(2,10))+&
                                slice(node_c(1,11),node_c(2,11))+&
                                slice(node_c(1,12),node_c(2,12)))
   end do
end do
!
end subroutine stencil


end module mod_maps_d
