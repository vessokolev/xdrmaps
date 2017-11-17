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
# Synopsis: This scripts reads a GRO file (in Gromos87
#           format) and assigns atomic mass to each of
#           the atom found there based on the first
#           symbol of the atom name. Note that the
#           atomic masses are supplied as Python
#           dictionary which have to be created prior to
#           the parsing of the GRO input file. The
#           result of the successful execution of the
#           script will be a HDF5 data set, stored in a
#           file. The name of the produced HDF5 datas
#           set is fixed to "mass" and MUST NOT be
#           changed.


# Python dictonary with the atomic masses (you migh need to modify that):
#
# List the atomic mass of each and every type of atom presented in the
# analyzed GRO file. Note the atom types are recognized by the first
# two symbols of the atom name listed in the atom name field of the
# GRO file. If you wanna change that default behaviour of the script,
# modify the code bellow accodringly.
atomic_mass_key={\
                  'H ':1.00794, \
                  'HA':1.00794, \
                  'HB':1.00794, \
                  'HD':1.00794, \
                  'HE':1.00794, \
                  'HG':1.00794, \
                  'HH':1.00794, \
                  'HN':1.00794, \
                  'HW':1.00794, \
                  'HZ':1.00794, \
                  'N ':14.00674, \
                  'ND':14.00674, \
                  'NE':14.00674, \
                  'NH':14.00674, \
                  'NZ':14.00674, \
                  'C ':12.0107, \
                  'CA':12.0107, \
                  'CB':12.0107, \
                  'CD':12.0107, \
                  'CE':12.0107, \
                  'CG':12.0107, \
                  'CH':12.0107, \
                  'CZ':12.0107, \
                  'O ':15.9994, \
                  'OD':15.9994, \
                  'OE':15.9994, \
                  'OG':15.9994, \
                  'OH':15.9994, \
                  'OT':15.9994, \
                  'OW':15.9994, \
                  'SD':32.065, \
                  'SG':32.065, \
                  'NA':22.9898,\
                  'CL':35.453
                }

dataset_name='mass' # DO NOT CHANGE THIS DATASET NAME!!!


###########################################
# DO NOT MODIFY THE CODE BELLOW THIS LINE #
#  UNLESS YOU KNOW EXACTLY WHAT YOU WANT  #
#           TO ACHIEVE BY THAT            #
###########################################


def get_masses_from_gro_file(gro_file,atomic_mass_key):
    # This function reads a GRO file in Gromos87 format, parses the
    # atom name column there, and reads the first symbol of each atom
    # name. Later that symbol is used as a key to receive the atomic
    # mass from the Python dictionary "atomic_mass_key" (see how it is
    # defined above). This function returns a list with the atomic
    # masses as 1D array. The position (index) of each atomic mass in
    # that lists matches the atom serial number of the corresponding
    # atom in GRO file + 1 (the atom serial numbers start at 1, while
    # the first index in the Python lists is 0). Note that the atom
    # serial number in the GRO file is not anymore the number listed
    # in the corresponding column. Nowadays the atom serial number is
    # taken as the consecutive number of the row the atom is described
    # into. The first two lines of the GRO files are excluded (the
    # first one contains comment, and the second one provides the
    # total number of the atoms).

    numbers='01234567890'

    gro_file_obj=open(gro_file,'r')

    comment=gro_file_obj.readline()

    numatoms=int(gro_file_obj.readline())

    atomic_mass=[0.0 for i in range(numatoms)]

    for i in range(numatoms):
       line=gro_file_obj.readline()
       aname=line[10:15]
       aname=aname.split()
       if aname[0][0] in numbers:
          print('FATAL ERROR: At least one atom name starts with number!')
          exit()
       else:
          if len(aname[0])==1:
             atomic_mass[i]=atomic_mass_key[aname[0][0]+' ']
          else:
             if not (aname[0][1] in numbers):
                atomic_mass[i]=atomic_mass_key[aname[0][0:2]]
             else:
                atomic_mass[i]=atomic_mass_key[aname[0][0]+' ']

    gro_file_obj.close()

    return atomic_mass

# Parse the invoking command line parameters:

flag=False

import sys

if len(sys.argv)<3:
   print('\nIncorrect number of input paramters!\n')
   print('Execute the scipt in the following manner:\n')
   print(sys.argv[0]+' input.gro output.h5\n')
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
      print('FATAL ERROR: The input file '+sys.argv[1]+\
            ' does not exist!')
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
         print('FATAL ERROR: The output file '+sys.argv[2]+\
               'cannot be created because either the target '+\
               'folder does not exists or the permissions '+\
               'set on that folder does not allow the '+\
               'requested file creation!')

if flag:

   # Parse the GRO file and get the atomic masses:
   atomic_mass=get_masses_from_gro_file(sys.argv[1],atomic_mass_key)

   # Open the HDF5 file for writing
   h5py_obj=h5py.File(sys.argv[2],'w')
   # define the dataset size and data type
   dset=h5py_obj.create_dataset(dataset_name,(len(atomic_mass),), \
                                dtype=numpy.float32)
   # and transfer the content of the list atomic_mass into the
   # dataset:
   dset[:]=numpy.array(atomic_mass)
   # Finally, close the HDF5 file and destroy the HDF5 object defined
   # in the memory:
   h5py_obj.close()

