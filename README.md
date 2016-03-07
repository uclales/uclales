name: Eckhard Dietze
date: Sun 13.09.2015

# README

This version of the UCLA-LES uses level-set-based front tracking to accurately represent sharp cloud boundaries instead of resolving them. The code is set up to simulate Bretherton et al.'s [1] smoke cloud and this case has been used to study the code's dependency on different choices in the numerics, with and without the Level Set Method. The results have been published in

> Dietze, E., Schmidt, H., Stevens, B., & Mellado, J. P. (2015). [Controlling entrainment in the smoke cloud using level set-based front tracking](https://www.schweizerbart.de/papers/metz/detail/23/84461/Controlling_entrainment_in_the_smoke_cloud_using_level_set_based_front_tracking). *Meteorologische Zeitschrift*, 23(6), 661–674. doi:10.1127/metz/2014/0595

This document describes how the code was set up and where the varied parameters can be set in order to reproduce the results.

An newer version of this code, configured to simulate the stratocumulus-topped atmospheric boundary layer, can be found in the [level-set-stbl](https://github.com/uclales/uclales/tree/level-set-stbl) branch.

For a more detailed description of how to use the UCLA-LES, please refer to the user documentation in ./doc/les_doc.pdf.

[1] Bretherton, C. S. et al. (1999). An intercomparison of radiatively-driven entrainment and turbulence in a smoke cloud, as simulated by different numerical models. *Quarterly Journal of the Royal Meteorological Society*, 125(554), 391–423.  doi:10.1002/qj.49712555402

## Set parameters

A set of 12 simulations was carried out with the following parameters varied:

1. Level set method (on, off)
2. Grid resolution (64 x 64 x 50, 128 x 128 x 100)
3. Scalar flux limiter (Minmod, Superbee, Monotonized Central)

The first two parameters are set in the main configuration file ./bin/NAMELIST. The last one is set  in the code directly in the scalar advection module ./src/advf.f90.

## NAMELIST file

### Level set method
The level set method is switched on or off by setting the `levset` parameter. The following settings are recognized:

	levset = 0	Level set is switched off.
	levset = 1	Level set is switched on.
	levset = 2	Like 1, but with additional output saved to the analysis file (<case>.nc)

### Grid resolution
The resolution is set by specifying values for `nxp`, `nyp`, `nzp`, and `deltax`, `deltay`, and `deltaz`. Note that `nxp` and `nyp` each include four ghost cells (two at each boundary) and that `nzp` includes one ghost cell below the ocean. Thus, for the computational domain of 3.2 km x 3.2 km x 1.25 km, we have

	nxp = 68
	nyp = 68
	nzp = 51
	deltax  = 50.0
	deltay  = 50.0
	deltaz  = 25.0
	
and 

	nxp = 132
	nyp = 132
	nzp = 101
	deltax  = 25.0
	deltay  = 25.0
	deltaz  = 12.5

In addition, by setting

	dzrat = 1.0
	
equidistant spacing is used in the vertical. Spacing in the horizontal directions is always equidistant.

## Scalar advection module advf.f90
### Scalar flux limiter
The flux limiter for scalar advection can be set at the beginning of ./src/advf.f90. via the `lmtr` variable. The following three settings were used:

	integer :: lmtr = 1 ! Minmod
	integer :: lmtr = 2 ! Superbee
	integer :: lmtr = 3 ! Monotonized Central
	
In addition, two simulations without the level set method and the limiter switched off (lmtr = 0)  were run.

## Running simulations

The code can be compiled on the terminal using make. For serial runs, type

	$ make seq

and for parallel runs, using the Message Passing Interface, type

	$ make mpi

After successful compilation, the code can be run serially

	$ nohup ./les.seq &

or in parallel

	$ nohup mpirun -n <number of MPI tasks> ./les.mpi &

In both cases, `nohup` is used to send the processes in the background and `&` is used to return to the terminal immediately after starting.

To specify the compiler you want to use (e.g gfortran or ifort), set it in the main makefile ./bin/Makefile by specifying

	COMPILER  := 'GFORTRAN'	

or

	COMPILER  := 'IFORT'	

Also, make sure that the `NCDF` variable in the selected compiler section of the makefile points to where NetCDF is installed on your system. E.g. on Ubuntu, using the system NetCDF package, specify

	NCDF = /usr	

NetCDF should be compiled with the same compiler that you are using to compile UCLA-LES.

