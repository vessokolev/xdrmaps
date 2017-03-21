#!/usr/bin/env python2

########################################################
#                                                      #
#                        XDRMAPS                       #
#                                                      #
#    TOOL FOR PARSING TRR TRAJECTORIES AND CREATING    #
#     GRIDS CONTAINING THE SPATIAL DISTRIBUTION OF     #
#              SELECTED ATOMIC PROPERTIES              #
#                                                      #
########################################################

# Version: 2017031400'
# Author: Veselin Kolev <vesso.kolev@gmail.com>'
# License: GPLv2'
#
# Synopsis: This scripts reads a GRO file (in Gromos87)
#           format and selects the atoms to participate
#           in the grid statistics. It selects only the
#           atoms that have certain name and the same
#           time belonging to a certain residue name.


dataset_name='presence' # Do not change this name!

# ###################################################################
# #                                                                 #
# # The format of GRO (Gromos87) files is described in details at:  #
# #                                                                 #
# # http://manual.gromacs.org/current/online/gro.html               #
# #                                                                 #
# ###################################################################

# When reading the information about each of the atoms presented in a
# GRO (Gromos87) file in Python, one should keep into account the
# exact positions of the columns. It is because the selection process
# is entirely based on the information listed there. Note that the
# first two lines are not part of column structure. The first line is
# only a place where one can leave a comment. As for the second one
# it provides the number of the atoms listed in the investigated
# topology. All the lines bellow that line, except the last one (it
# contains the box unit vectors), and part of the column structure
# which presents the topology and they need to be read.
#
# The easiest way to read GRO file structure in Python is to assign
# each line content to a string variable and then check only those
# positions of the string which correspond to a certain column of
# interest. For instance, if the line content is assigned to the
# string variable "line", the residue id, residue name, and the atom
# name can be extracted as:
#
#  resid=line[0:5]
#  resname=line[5:10]
#  aname=line[10:15]


def select_certain_atoms_from_resdue(gro_file):
   # This is an example illustrating how to organize the selection
   # of atoms based on the residue name. Only the atoms belonging
   # to the residue with name SOL are selected bellow!

   atom_name_beginning_with='O2'
   residue_name='GMO'

   # Read the file efficiently whithout loading its content into the
   # memory.
   with open(gro_file,'r') as file_obj:

      # The counter takes into account the current line number in the
      # GRO file. Set it to 0 at the beginning.
      counter=0

      for line in file_obj:

         if counter==1:
            # This is the position of the second line in the GRO
            # file and it is supposed to contain the number of atoms
            # in the topology. Let's read that line and convert its
            # content to an inteleger number:
            numatoms=int(line)
            # Once the number of atoms is known, create and initialize
            # the array atom_presence and set all its elements 0, which
            # means that by default all atoms are excluded from the
            # selection. Later only the atoms which belong to the
            # residue with name SOL will be presented in atom_presence
            # array (the elements representing them in the presence
            # array will be set 1).
            atom_presence=[0 for i in xrange(numatoms)]
         else:
            if counter>1:
               # Parsing the lines describing the atoms in the topology.
               # Check only the name of the residue (the variable
               #"resname" bellow).
               resname=line[5:10]
               aname=line[10:15]
               if resname.split()[0]==residue_name and \
                  aname.split()[0][0:2]==atom_name_beginning_with:
                  # Pay attention how the counter variable value is
                  # shifted backward. That is deliberately done to
                  # compensate the accumulation of 2 when reading
                  # of the first two lines in GROM file, which contain
                  # only meta data but no atom description.
                  atom_presence[counter-2]=1
         counter+=1

   return atom_presence


# Parse the invoking command line parameters:

flag=False

import sys

if len(sys.argv)<3:

   print '\nIncorrect number of input paramters!\n'
   print 'Execute the scipt in the following manner:\n'
   print sys.argv[0]+' input.gro output.h5\n'

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
      flag=False
      print 'FATAL ERROR: The input file '+sys.argv[1]+\
            ' does not exist!'
      flag=False
   else:
      # Check if the output file could be created
      try:
         file_obj=open(sys.argv[2],'w')
         flag=True
      except:
         flag=False
      
      if flag:
         os.remove(sys.argv[2])
      else:
         print 'FATAL ERROR: The output file '+sys.argv[2]+\
               'cannot be created because either the target '+\
               'folder does not exists or the permissions '+\
               'set on that folder does not allow the '+\
               'requested file creation!'

if flag:

   # Parse the GRO file and get the atom presence array:
   atom_presence=select_certain_atoms_from_resdue(sys.argv[1])

   # Open the HDF5 file for writing
   h5py_obj=h5py.File(sys.argv[2],'w')

   # define the dataset size and data type
   dset=h5py_obj.create_dataset(dataset_name,(len(atom_presence),), \
                                dtype=numpy.int)

   # and transfer the content of the list atomic_mass into the
   # dataset:
   dset[:]=numpy.array(atom_presence)

   # Finally, close the HDF5 file and destroy the HDF5 object defined
   # in the memory:
   h5py_obj.close()

