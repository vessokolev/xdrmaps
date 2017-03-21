      program xdrmaps_s

      use iso_c_binding,only:C_SHORT,C_INT,C_FLOAT
      use mod_maps_d
      use mod_cmd_line_d
      !
      implicit none
      !
      character(len=4096) :: grid_file
      character(len=4096) :: path
      character(len=4096) :: map_file_prefix
      integer(C_INT) :: i
      integer(C_INT) :: j
      integer(C_INT) :: k
      integer(C_INT) :: sel_dim
      integer(C_INT) :: code
      integer(C_INT) :: num_atoms
      integer(C_INT) :: num_frames
      integer(C_INT) :: slice_width
      integer(C_INT),dimension(3) :: num_nodes
      integer(C_INT),dimension(3) :: grid_s
      integer(C_INT),dimension(2) :: slice_s
      integer(C_INT),dimension(2,12) :: stencil_node
      integer(C_INT),dimension(:,:,:),allocatable :: freq
      real(C_FLOAT),dimension(3) :: box
      real(C_FLOAT),dimension(3) :: step
      real(C_FLOAT),dimension(2) :: step_2D
      real(C_FLOAT),dimension(4) :: w
      real(C_FLOAT),dimension(:,:,:),allocatable :: grid
      real(C_FLOAT),dimension(:,:),allocatable :: slice
      real(C_FLOAT),dimension(:,:),allocatable :: slice_filtered
      logical :: noise_red
      !
      call check_input(grid_file,&
                       sel_dim,&
                       slice_width,&
                       noise_red,&
                       path)
      !
      ! Define the weights used for the stencil. Here w are the weights
      ! defined in equation E17 in Supporting Information of
      ! DOI: 10.1021/acs.jpcb.5b06566 
      !
      w(1)=0.07692307692307693_C_FLOAT ! this is 1/13
      w(2)=0.30769230769230769_C_FLOAT ! this is 4/13
      w(3)=w(2)
      w(4)=w(2)
      !
      ! Define the stencil nodes. That nodes supplies equation E17 in
      ! Supporting Information of DOI: 10.1021/acs.jpcb.5b06566 with the
      ! correct index values of i and j with respect to PBC.
      !
      stencil_node(:,1)=(/-1,0/)
      stencil_node(:,2)=(/0,-1/)
      stencil_node(:,3)=(/1,0/)
      stencil_node(:,4)=(/0,1/)
      !
      stencil_node(:,5)=(/-1,-1/)
      stencil_node(:,6)=(/-1,1/)
      stencil_node(:,7)=(/1,-1/)
      stencil_node(:,8)=(/1,1/)
      !
      stencil_node(:,9)=(/-2,0/)
      stencil_node(:,10)=(/0,-2/)
      stencil_node(:,11)=(/2,0/)
      stencil_node(:,12)=(/0,2/)
      !
      ! Read the code number first:
      !
      call read_scalar_set(grid_file,&
                           'code',&
                           code)
      !
      ! Read first the dimension of the grid and freq data sets:
      !
      call read_int_vector_set(grid_file,&
                               'num_nodes',&
                               num_nodes)
      !
      ! Based on the information given by num_nodes array dynamically
      ! allocate the arrays for storing the grid array elements:
      !
      allocate(grid(num_nodes(1),num_nodes(2),num_nodes(3)))
      !
      ! Read the grid array:
      !
      call read_grid(grid_file,&
                     'grid',&
                     grid)
      if (code.gt.1) then
         !
         ! Based on the information given by num_nodes array dynamically
         ! allocate the arrays for storing the freq array elements:
         !
         allocate(freq(num_nodes(1),num_nodes(2),num_nodes(3)))
         !
         ! Read the frequency array:
         !
         call read_freq_grid(grid_file,&
                             'freq',&
                             freq)
         !
         ! Set the map file prefix accoring to the processed
         ! property:
         !
         if (code.eq.2) then
            map_file_prefix='velocity_'
         else
            map_file_prefix='force_'
         end if
      else
         !
         ! Read the number of the frames:
         !
         call read_scalar_set(grid_file,&
                              'num_frames',&
                              num_frames)
         !
         ! In this particular case the grid contains the distribution of
         ! the atomic mass. Divide the values stored in the dataset into
         ! the total number of the frames processed:
         !
         grid=grid/num_frames
         !
         ! Set the map file prefix accoring to the processed property:
         !
         map_file_prefix='mass_'
      end if
      !
      ! Read the step array:
      !
      call read_real_vector_set(grid_file,&
                                'step',&
                                step)

      !
      ! Define the slice array dimensions:
      !
      slice_s=get_slice_shape(num_nodes,&
                              sel_dim)
      
      allocate(slice(slice_s(1),slice_s(2)))
      allocate(slice_filtered(slice_s(1),slice_s(2)))

      step_2D=get_step_2D(step,&
                          sel_dim)

      if (code.eq.1) then

         do i=1,num_nodes(sel_dim)

            call get_mass_slice(num_nodes,&
                                grid,&
                                sel_dim,&
                                i,&
                                slice_width,&
                                slice_s,&
                                slice)

            if (noise_red) then
               call stencil(slice_s,&
                            slice,&
                            stencil_node,&
                            w,&
                            slice_filtered)
               call write_map_file(path,&
                                   map_file_prefix,&
                                   i,&
                                   slice_s,&
                                   slice_filtered,&
                                   step)
            else
               call write_map_file(path,&
                                   map_file_prefix,&
                                   i,&
                                   slice_s,&
                                   slice,&
                                   step)
            end if
            print *,'Successfully processed the map number',i

         end do

      end if

      if ((code.eq.2).or.(code.eq.3)) then

         do i=1,num_nodes(sel_dim)

            call get_vf_slice(num_nodes,&
                              grid,&
                              freq,&
                              sel_dim,&
                              i,&
                              slice_width,&
                              slice_s,&
                              slice)
            if (noise_red) then
               call stencil(slice_s,&
                            slice,&
                            stencil_node,&
                            w,&
                            slice_filtered)

               call write_map_file(path,&
                                   map_file_prefix,&
                                   i,&
                                   slice_s,&
                                   slice_filtered,&
                                   step)
            else
               call write_map_file(path,&
                                   map_file_prefix,&
                                   i,&
                                   slice_s,&
                                   slice,&
                                   step)
            end if
            print *,'Successfully processed the map number',i

         end do

      end if

      call exit(0)

      end program xdrmaps_s
