#!/usr/bin/env python3

########################################################
#                                                      #
#                        XDRMAPS                       #
#                                                      #
#    TOOL FOR PARSING TRR TRAJECTORIES AND CREATING    #
#     GRIDS CONTAINING THE SPATIAL DISTRIBUTION OF     #
#              SELECTED ATOMIC PROPERTIES              #
#                                                      #
########################################################

# Version: 2017111700'
# Author: Veselin Kolev <vesso.kolev@gmail.com>'
# License: GPLv2'
#
# Synopsis: This Python 3 code analyzes the content of the
#           grid maps stored in the output files produced
#           by running the tool xdrmaps. It displays the
#           information about the HDF5 data sets stored
#           there, so one can get information about: the
#           sizes of the simulation box, the number of
#           the nodes, the step corresponding to the
#           distance between each two neighboring nodes,
#           as well as the type of invesigated atomic
#           property.


import h5py
import sys

if len(sys.argv)<2:
   print('\nIncorrect number of input paramters!\n')
   print('Execute the scipt in the following manner:\n')
   print((sys.argv[0]+' maps.h5\n'))
else:

   import h5py  # This Python module can be installed by using
                # pip, as:
                # $ pip install h5py --user
   import numpy # This Python module can be installed by using
                # pip, as:
                # $ pip install numpy --user
   import os

   # Check if the input file exists:
   if not os.path.isfile(sys.argv[1]):
      print(('\nFATAL ERROR: The input file with maps '+sys.argv[1]+\
            ' does not exist!\n'))
      flag=False

   # Open the file containing the HDF5 datasets:
   try:

      h5_obj=h5py.File(sys.argv[1], "r")
   except:
      print('\nFATAL ERROR: The input file:\n')
      print((sys.argv[1]+'\n'))
      print('does not contain HDF5 structure!\n')
      sys.exit(1)

   # First, check if the data sets are presented by confirming that the datasets with the
   # following names are accessible in the HDF5 container:
   required_dset_names=['box','code','num_atoms','num_frames','num_nodes','step','grid','freq']

   read_dset_names=[i for i in list(h5_obj.keys())]

   # This kind of check is not the fastest one, but it would be improved in the future
   # versions of the program.
   counter=0
   for i in required_dset_names:
      if i in read_dset_names:
         counter+=1

   if counter!=len(required_dset_names):
     print('\nFATAL ERROR: The analyzed file:\n')
     print((sys.argv[1]+'\n'))
     print('does not contain the required HDF5 datasets!\n')
     print(('Required datasets: ',','.join(required_dset_names)))
     print(('Available datasets: ',','.join(read_dset_names)+'\n'))
     h5_obj.close()
     sys.exit(1)

   # Check the data type of each dataset:

   for i in required_dset_names:

      if i=='box':
         if type(h5_obj[i][:])!=numpy.ndarray:
            print(('\nFATAL ERROR: The dataset "'+i+'" is not an object of type '+\
                  'numpy.ndarray!'))
            sys.exit(1)
         else:
            if type(h5_obj[i][0])!=numpy.float32:
               print(('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                     'data type (numpy.float32).\n'))
               sys.exit(1)
            else:
               if h5_obj[i][:].shape!=(3,):
                  print(('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                        'shape (3).\n'))
                  sys.exit(1)

         box=h5_obj[i][:]


      if i=='code':
         if type(h5_obj[i][:])!=numpy.ndarray:
            print(('\nFATAL ERROR: The dataset "'+i+'" is not an object of type '+\
                  'numpy.ndarray!'))
            sys.exit(1)
         else:
            if type(h5_obj[i][0])!=numpy.int32:
               print(('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                     'data type (numpy.int32).\n'))
               sys.exit(1)
            else:
               if h5_obj[i][:].shape!=(1,):
                  print(('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                        'shape (1).\n'))
                  sys.exit(1)
         code=h5_obj[i][:]


      if i=='num_atoms':
         if type(h5_obj[i][:])!=numpy.ndarray:
            print(('\nFATAL ERROR: The dataset "'+i+'" is not an object of type '+\
                  'numpy.ndarray!'))
            sys.exit(1)
         else:
            if type(h5_obj[i][0])!=numpy.int32:
               print(('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                     'data type (numpy.int32).\n'))
               sys.exit(1)
            else:
               if h5_obj[i][:].shape!=(1,):
                  print(('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                        'shape (1).\n'))
                  sys.exit(1)

         num_atoms=h5_obj[i][:]


      if i=='num_frames':
         if type(h5_obj[i][:])!=numpy.ndarray:
            print(('\nFATAL ERROR: The dataset "'+i+'" is not an object of type '+\
                  'numpy.ndarray!'))
            sys.exit(1)
         else:
            if type(h5_obj[i][0])!=numpy.int32:
               print(('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                     'data type (numpy.int32).\n'))
               sys.exit(1)
            else:
               if h5_obj[i][:].shape!=(1,):
                  print(('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                        'shape (1).\n'))
                  sys.exit(1)

         num_frames=h5_obj[i][:]


      if i=='num_nodes':
         if type(h5_obj[i][:])!=numpy.ndarray:
            print(('\nFATAL ERROR: The dataset "'+i+'" is not an object of type '+\
                  'numpy.ndarray!'))
            sys.exit(1)
         else:
            if type(h5_obj[i][0])!=numpy.int32:
               print(('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                     'data type (numpy.int32).\n'))
               sys.exit(1)
            else:
               if h5_obj[i][:].shape!=(3,):
                  print(('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                        'shape (3).\n'))
                  sys.exit(1)

         num_nodes=h5_obj[i][:]


      if i=='step':
         if type(h5_obj[i][:])!=numpy.ndarray:
            print('\nFATAL ERROR: The dataset "'+i+'" is not an object of type '+\
                  'numpy.ndarray!')
            sys.exit(1)
         else:
            if type(h5_obj[i][0])!=numpy.float32:
               print('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                     'data type (numpy.float32).\n')
               sys.exit(1)
            else:
               if h5_obj[i][:].shape!=(3,):
                  print('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                        'shape (3).\n')
                  sys.exit(1)

         step=h5_obj[i][:]

      if i=='grid':
         if type(h5_obj[i][:])!=numpy.ndarray:
            print('\nFATAL ERROR: The dataset "'+i+'" is not an object of type '+\
                  'numpy.ndarray!')
            sys.exit(1)
         else:
            if type(h5_obj[i][0][0][0])!=numpy.float32:
               print('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                     'data type (numpy.float32).\n')
               sys.exit(1)
            else:
    
               if h5_obj[i][:].shape!=(num_nodes[0],num_nodes[1],num_nodes[2]):
                  print('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                        'shape ('+str(num_nodes[0])+','+str(num_nodes[1])+',',+str(num_nodes[2])+').\n')
                  sys.exit(1)


      if i=='freq':
         if type(h5_obj[i][:])!=numpy.ndarray:
            print('\nFATAL ERROR: The dataset "'+i+'" is not an object of type '+\
                  'numpy.ndarray!')
            sys.exit(1)
         else:
            if type(h5_obj[i][0][0][0])!=numpy.int32:
               print('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                     'data type (numpy.int32).\n')
               sys.exit(1)
            else:

               if h5_obj[i][:].shape!=(num_nodes[0],num_nodes[1],num_nodes[2]):
                  print('\nFATAL ERROR: The dataset "'+i+'" does not match the expected '+\
                        'shape ('+str(num_nodes[0])+','+str(num_nodes[1])+',',+str(num_nodes[2])+').\n')
                  sys.exit(1)

   # Print the summary:
   print('\n*********************************************')
   print('\nSummary of the content found in the map file:\n')
   print(sys.argv[1]+'\n')
   print('*********************************************\n')
   print('(*) Number of atoms in the processed topology:\n')
   print('    '+str(num_atoms[0])+'\n')
   print('(*) Number of frames processed:\n')
   print('    '+str(num_frames[0])+'\n')
   print('(*) Box size (averaged over the frames):\n')
   print('    in x: '+str(box[0])+' nm\n')
   print('    in y: '+str(box[1])+' nm\n')
   print('    in z: '+str(box[2])+' nm\n')
   print('(*) Number of nodes in each direction:\n')
   print('    in x: '+str(num_nodes[0])+'\n')
   print('    in y: '+str(num_nodes[1])+'\n')
   print('    in z: '+str(num_nodes[2])+'\n')
   print('(*) Average step size in each direction:\n')
   print('    in x: '+str(step[0])+' nm\n')
   print('    in y: '+str(step[1])+' nm\n')
   print('    in z: '+str(step[2])+' nm\n')
   print('(*) The grid contains the spatial distribution of:\n')
   if code[0]==1:
      print('    atomic mass\n')
      sys.exit(0)
   if code[0]==1:
      print('    atomic velocity\n')
      sys.exit(0)
   if code[0]==2:
      print('    force upon atoms\n')
      sys.exit(0)
   print('    unknown property code '+str(code[0]))
   
