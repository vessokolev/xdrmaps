## XDRMAPS

### Tools for analyzing GROMACS TRR or XTC trajectories, and obtaining the spatial distribution of selected atomic properties as a grid statistics.

#### Author: Veselin Kolev <vesso.kolev@gmail.com>
#### Released under GPLv2 licence.

#### Content:

#### 1. Introduction.
#### 2. How to download the source code of the project.
#### 3. Compilation of the source code.
#### 4. Preparation of the mass and presence arrays as HDF5 data sets.
#### 5. Collection of the grid statistics as 3D array  HDF5 data set.
#### 6. Storing the spatial distribution as a series of 2D maps.
#### 7. Displaying meta info about the content of the grid data set.
#### 8. Plotting the 2D maps with gnuplot.


_1. Introduction_

The xdrmaps package implements a protocol for collecting the spatial distribution of certain atomic properties like the atomic mass, the atomic velocity, or the aforce cting upon each atom, as an average by ensemble at a series of points across the volume box - nodes. By introducing nodes the sought spatial distribution is collected as a grid statistics - the simulation box is divided into sub-boxes by means of the positions of the nodesa and the investigated atom property in each of the boxes is assigned to the closest node. The xdrmaps tool reads the content of the simulation box by processing each frame stored in GROMACS TRR or XTC trajectories.

The mathematical formalism used to implement the grid statisitcs collection is described in details in section S6 of the "Supporting Information" material of DOI: 10.1021/acs.jpcb.5b06566:

Molecular Dynamics Investigation of Ion Sorption and Permeation in Desalination Membranes. V. Kolev, V. Freger, *J. Phys. Chem. B*, 2015, 119 (44), pp 14168–14179.

The "Supporting Information" material could be downloaded by following the link: http://pubs.acs.org/doi/suppl/10.1021/acs.jpcb.5b06566

If you cannot download that file send me an e-mail request.


_2. How to download the source code of the project._

The preferable method for obtaining the code is to use the tool ``git`` and clone the source tree locally to your file system:

```
git clone https://github.com/vessokolev/xdrmaps.git
```

You may download the source as a ZIP-archive by pressing the button "Clone or download" in the web-interface on GitHub, or by using wget:

```
wget https://github.com/vessokolev/xdrmaps/archive/master.zip
```


_3. Compilation of the source code._

You can compile the source code of xdrmaps by using any modern C and FORTRAN compilers under Linux. Refer to the files:


```
README.compile_gcc
README.compile_intel
README.compile_pgi

```

to obtain the minimum set of compiler options needed to compile fast executable code. Note that the compilation of the source code of xdrmaps project depends on the following external projects:

xdrfile-1.1.4: ftp://ftp.gromacs.org/contrib/xdrfile-1.1.4.tar.gz
xdrfort: https://github.com/kmtu/xdrfort.git

Please, strictly follow the instructions given in the corresponding README files listed above on how to combine the source code and proceed with the compilation. As a result two binary tools will be built - ``xdramaps`` and ``xdramps_d``.


_4. Preparation of the mass and presence arrays as HDF5 data sets._

The xdrmaps tool requires to be supplied with the atomic mass of each atom participating in the topology. That particular array have to be stored as HDF5 data set in a separate file. The document ``README.creating_atomic_mass_dset`` explains how to create the atomic mass array and how to save it as HDF5 file set on a file.

Having a presence array makes it possible to mark which of the atoms in the topology should be involved in the grid statisitcs collection process. It is analoque of ``*.ndx`` files implemented by GROMACS. The presence array have to be stored as HDF5 data set on a file. The procecedure for creating a presence array is described in the file ``README.creating_presence_dset``.

One need to supply both the mass array and at least the default presence array, to run the grid statistics.


_5. Collection of the grid statistics as 3D array  HDF5 data set._

The collection of the grid statistics is performed by running the tool ``xdrmaps``. Running it without input parameters shows the help on the possible command line options. Currently ``xdrmaps`` processes TRR and XTC trajectories. Note that the XTC trajectories can be involved only in the collection of the spatial distribution of the atomic mass. Also note that the results obtained by parsing TRR trajectories are more smooth and reliable then the ones collected by means of XTC trajectories.

Currently ``xdrmaps`` supports the collection of spatial distribution of the atomic mass, atomic velocity, and the force acting upon atoms. It is possible also to collect the spatial distribution of the atomic charge if the atomic mass saved as elements of the mass array is replaces by the partial atomic charges.

The grid statistics collected by ``xdrmaps`` is created as an array in the memory of the system and then stored as HDF5 data set.


_6. Storing the spatial distribution as a series of 2D maps._

The tool ``xdrmaps_d`` saves the collected by ``xdrmaps`` grid statistics as a series of two-dimensional maps (isolines). That method is based on slicing the grid array in arbitrary direction (x, y, or z). Once the maps are saved as files they can be plotted as 2D maps (see the explanation how to do that by employing gnuplot given next).

Run ``xdrmaps_d`` without input parameters to obtain a list of the possible starting options.

The tool ``xdramaps_d`` implements a stencil for reducing the noise level. Its shape and way of definint the nodes is described in Figure S4 in "Supporting Information" of DOI: 10.1021/acs.jpcb.5b06566.


_7. Displaying meta info about the content of the grid data set._

This package sipplies a Python 2 script which might help to print out an information about the saved HDF5 data set, obtained as a result of the grid statistics (that process is explained in 5). The script is named ``check_maps_dset.py`` is located inside the ``python/`` folder. It should be invoked as:

```
check_maps_dset.py maps.h5
```

The output may look like:

```
*********************************************

Summary of the content found in the map file:

maps.h5

*********************************************

(*) Number of atoms in the processed topology:

    102851

(*) Number of frames processed:

    10001

(*) Box size (averaged over the frames):

    in x: 10.1335 nm

    in y: 10.1335 nm

    in z: 10.1335 nm

(*) Number of nodes in each direction:

    in x: 506

    in y: 506

    in z: 506

(*) Average step size in each direction:

    in x: 0.0200267 nm

    in y: 0.0200267 nm

    in z: 0.0200267 nm

(*) The grid contains the spatial distribution of:

    force upon atoms
```


_8. Plotting the 2D maps with gnuplot._

8.1. Plotting the atomic mass.

To help explaining how the maps of the atomic mass should be visualized, the source tree of the project contains the following data sets:

```
mass_orig.dat
mass_nr.dat
```

and Gnuplot templates:

```
mass_orig.gnuplot
mass_nr.gnuplot
```

Here "orig" stands for "original" (i.e. not processed in order to remove the noise), and "nr" means "noise-reduction" - that is a map passed a noise reduction. To visualize the maps as EPS files the Gnuplot templates should be processed:

```
gnuplot mass_orig.gnuplot
gnuplot mass_nr.gnuplot
```

In case of successful processing the following EPS files will be created:

```
mass_orig.eps
mass_nr.eps
```

They might be replotted as PNG files with very high quality:

```
convert -flatten -density 300 mass_orig.eps mass_orig.png
convert -flatten -density 300 mass_nr.eps mass_nr.png
```

8.2. Plotting the atomic velocity.

To help explaining how the maps of the atomic velocity should be visualized, the source tree of the project contains the following data sets:

```
velocity_orig.dat
velocity_nr.dat
```

and Gnuplot templates:

```
velocity_orig.gnuplot
velocity_nr.gnuplot
```

Here "orig" stands for "original" (i.e. not processed in order to remove the noise), and "nr" means "noise-reduction" - that is a map passed a noise reduction. To visualize the maps as EPS files the Gnuplot templates should be processed:

```
gnuplot velocity_orig.gnuplot
gnuplot velocity_nr.gnuplot
```

In case of successful processing the following EPS files will be created:

```
velocity_orig.eps
velocity_nr.eps
```

They might be replotted as PNG files with very high quality:

```
convert -flatten -density 300 velocity_orig.eps velocity_orig.png
convert -flatten -density 300 velocity_nr.eps velocity_nr.png
```

8.3. Some notes on Gnuplot templates.

One must use in advance the script ```check_maps_dset.py``` as shown in 7. That shows the box sizes and the grid spacial parameter (that is the distance between two neighboring nodes, also step when divide the box into sub-boxes). The box size numbers are used to define the xrange and yrange definitions in the Gnuplot template:

```
set xrange [0:10.1335]
set yrange [0:10.1335]
```

while the grid spatial parameter becomes part of the plot definition:

```
plot 'mass_ng.dat' u ($1*0.0200267):($2*0.0200267):3 matrix with image
```

Note that the format of the files containing the maps does not include any of the abscissa or ordinate. Each line of the file contains only information about the values assigned to the "pixels" of each of the lines of the plotted image.
