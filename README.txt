
            The AcCoRD Simulator
            (Actor-based Communication via Reaction-Diffusion)

This document is the README for AcCoRD v1.0 (2016-10-31)

# Introduction

Welcome to the AcCoRD Simulator (Actor-based Communication via Reaction-Diffusion). AcCoRD is a simulation tool for molecular communication. It is a hybrid microscopic/mesoscopic simulator that functions as a generic solver of stochastic reaction-diffusion but which has been developed from the perspective of communications analysis. The focus is to efficiently generate the statistics of molecule observations at some location (i.e., at a receiver).

Primary features include:
* As a hybrid solver, AcCoRD integrates two simulation models to define 3D environments with flexible local accuracy.  Microscopic regions define each molecule individually and evolve over discrete time steps. Mesoscopic regions count the number of molecules in disjoint virtual bins (called subvolumes) and evolve over time steps with continuous granularity.
* "Actors" are either active molecule sources (i.e., transmitters) or passive observers (i.e., receivers). Active actors can modulate binary sequence data and release molecules over specified intervals. Passive actors record the number of molecules observed and can also record molecule positions.
* Chemical reactions can be defined locally or globally. A general framework for defining reactions can accommodate reactions such as molecule degradation, enzyme kinetics, reversible or irreversible surface binding, ligand-receptor binding, transitions across boundary membranes, and simplified molecular crowding. Surface binding reactions include absorption (consumption), adsorption (sticking), and desorption (release).
* Independent realizations of a simulation can be repeated an arbitrary number of times (and on different computers) and then aggregated to determine the average behavior and channel statistics.
* AcCoRD's interface has been designed to be helpful to the novice user while providing meaningful output. AcCoRD also includes post-processing tools developed in MATLAB. These tools enable saving aggregated simulation output files, plotting receiver observations (either the time-varying behavior or empirical distributions at specified times), and visualizing the physical environment (either as still images or compiled into a video)

Quicklinks:
* Download latest stable version from http://adamnoel.ca/software/accord/download (user-friendly page) or https://github.com/adamjgnoel/AcCoRD/releases (developer page)
* Journal paper (draft): http://adamnoel.ca/software/accord/publications
* Usage documentation: http://adamnoel.ca/software/accord/user
* Sample configuration files: http://adamnoel.ca/software/accord/examples
* Old AcCoRD releases: https://github.com/adamjgnoel/AcCoRD/releases
* Main AcCoRD development page on Github: https://github.com/adamjgnoel/accord
* Maintenance of existing issues and planned features: https://github.com/adamjgnoel/AcCoRD/issues
* Github wiki (just a copy of this README): https://github.com/adamjgnoel/AcCoRD/wiki

## Audience

If you are a communications researcher who is interested in studying molecular communication, then you are the intended audience. However, anyone interested in the multi-scale simulation of reaction-diffusion systems might also find AcCoRD useful.

## Development Status

Development of AcCoRD is on-going and hosted on Github (https://github.com/adamjgnoel/AcCoRD).

Existing bugs and planned features are described on the Issues page: https://github.com/adamjgnoel/AcCoRD/issues

Most identified content on the Issues page are planned features.

# Background

AcCoRD development began in October 2014 with the motivation to develop a reaction-diffusion simulation tool that would be suitable for the communications research community. The goal was to integrate recent advances in stochastic reaction-diffusion models. The overall design of AcCoRD is very much inspired by Smoldyn (http://www.smoldyn.org/), which is also a microscopic reaction-diffusion solver (and which, coincidentally, also added mesoscopic regions in 2015). The key differences with Smoldyn relate to AcCoRD's communication focus. AcCoRD's intended use is to capture the dynamics of a communications system, i.e., with a transmitter and receiver. Transmitters are implemented as "active actors", which manually add molecules to the environment. Observations by receivers are recorded as "passive actors". The placement of actors is defined by the user, and the only information that AcCoRD saves is that related to actor behavior. Furthermore, simulations are designed to be run repeatedly (i.e., an arbitrary number of times, but could easily be tens to hundreds of thousands of times), so that accurate receiver statistics can be generated.

AcCoRD is developed in C (C99) for precise control over memory management and access to very fast algorithms for random number generation, linked list operations, and related computational tasks. Some utilities have been written (and will continue to be written) in MATLAB to facilitate post-processing.

A formal publication on the design of AcCoRD will be submitted in November 2016. For now, the general motivation and some of the implementation ideas were presented in the following conference papers:  
- A. Noel, K.C. Cheung, and R. Schober, On the Statistics of Reaction-Diffusion Simulations for Molecular Communication, in Proc. ACM NANOCOM 2015, Sep. 2015. DOI: http://dx.doi.org/10.1145/2800795.2800821 - the simulations in this paper were completed with AcCoRD v0.1, which was an early 2D build 
- A. Noel, K.C. Cheung, and R. Schober, Multi-Scale Stochastic Simulation for Diffusive Molecular Communication, in Proc. IEEE ICC 2015, pp. 2712--2718, Jun. 2015. DOI: http://dx.doi.org/10.1109/ICC.2015.7248471 - the simulations in this paper were completed with proof-of-concept MATLAB code 

Please contact the developer or consider following AcCoRD on Github (https://github.com/adamjgnoel/accord) for the latest updates.

# Feature Summary

AcCoRD currently has the following "mature" features:
* 3D multi-scale hybrid reaction-diffusion
* Environment described as union of cubes and spheres. Cubes can be microscopic (track individual molecules) or mesoscopic (assume uniform molecule density throughout the cube). Spheres must be microscopic but can be infinite in size (i.e., unbounded environment). Diffusion coefficients can be locally defined.
* Chemical reactions of 0th, 1st, and 2nd order (though 2nd order reactions in microscopic regime have a preliminary implementation that relies on specifying a binding radius)
* Surfaces can surround all or part of 3D regions and be sites for local chemical reactions, such as molecule generation or absorption. User can choose how reaction probabilities are determined for these reactions.
* Readable configuration and simulation output files. Configuration is described using the JSON interchange format.
* Utility functions in MATLAB for post-processing. Features include importing simulation output, generating figures and video of simulation progress, and plotting simulation data and distributions. Simulation environments can be plotted without running the simulation first.
* Flexible placement of molecule sources and observers, whose location can be defined as their own shapes or as the union of a set of regions. Molecule sources can also be defined as points.
* Molecule sources based on "Concentration Shift Keying" (CSK) modulation of independent symbols, with user-defined pulse widths and switches for stochastic or deterministic molecule release. Symbols can also be pre-defined. A second option, "Burst", does not modulate data and is analogous to CSK with 1 modulation bit that always has value 1.
* Molecule observers record number of specified types of molecules (and optionally their positions) over a subset of the simulation environment.
* Even in the microscopic regime, molecule creation and some chemical reactions can occur on a time scale smaller than the microscopic time step.

The following features have not been as extensively tested and are considered in "beta":

The following features have had very limited testing and are considered "alpha":
* Environments described entirely in 2D, with diffusion across a union of 2D rectangles

A complete version history can be found in the CHANGELOG.txt:
* latest version is in the Github code repository (https://github.com/adamjgnoel/AcCoRD/blob/master/CHANGELOG.txt)
* latest STABLE PUBLIC version is included with the corresponding download

# Known Issues and Limitations

This list is current as of v1.0. See https://github.com/adamjgnoel/AcCoRD/issues for the latest details.
* Full 2D simulations are mostly untested
* 2nd order chemical reactions in microscopic regions need a binding radius specified.
* There is an option to set whether actor is independent but the value is ignored since dependent actors have not yet been implemented
* 2 regions must share an outer face in order to be classified as neighbors, unless one region is the parent of the other. So, if a child region is up against a face of its parent, and that parent also has a parent but doesn't share the same face with its parent, then the parent's parent will not be recognized as a neighbor of the child. In such a case molecule transitions between the child and the "grandparent" are not possible.
* Surface regions must be microscopic and can only have microscopic regions as neighbors.
* A child region is always assumed to be a neighbor of its parent, which is not strictly true but should not create a problem during a simulation
* microscopic molecules in a spherical region may go a very small distance beyond the region boundary but are still believe to be inside the sphere. This is due to numerical underflow, and might mean that the potential chemical reactions for the molecule are inconsistent with its location.
* A mesoscopic subvolume adjacent to a microscopic region needs that region to cover the entirety of every face that is adjacent to that region (i.e., part of the adjacent face should not also be adjacent to some other region). Otherwise the hybrid interface will not have proper transitions.
* Actors that act at the same time will do so in a random order (this is by design). So, if an active actor is adding molecules at the exact time that a passive actor is supposed to observe them, it is unknown a priori which will act first. It will depend on whichever is currently ranked higher in the timer heap. This could be modified in the future to rank the actors in the order they are defined in the event that their timer values are the same.
* Molecules leaving a mesoscopic region and entering a microscopic region might be placed outside the valid environment boundary.
* Unbinding radii are only applied to second order microscopic reactions that have multiple products.


# Future Features

Development of AcCoRD is active and on-going. See https://github.com/adamjgnoel/AcCoRD/issues for the latest details. Here is a list of some planned end-user features (in no particular order):
* Add flow. This can first be added as steady and uniform, and then varying locally.
* Allow dependent actors, whose behavior depends on the observations of passive actors. This would allow implementation of more complex communication networks (e.g., relays)
* Add actors that track changes to molecule composition in their volumes
* Improve implementation of 2nd order chemical reactions in microscopic regions
* Enable surfaces to mesoscopic regions. Could be achieved via exclusion regions
* Write more MATLAB utility functions to perform common tasks with simulation output (e.g., bit error probabilities)
* Use different microscopic time steps in different microscopic regions
* Add checks on subvolume size relative to diffusion coefficients
* Add "tau leaping" in mesoscopic regime to improve scalability
* Add actor mobility, such that actor location can vary over time (either randomly or deterministically)
* Add more modulation schemes for active actors


# Installation

There are two installation options:

1) (BASIC) Use a pre-compiled binary from the releases page (https://github.com/adamjgnoel/AcCoRD/releases).
Debug and optimized versions exist for Windows, Debian/Ubuntu Linux, and RHEL/CentOS Linux. No further action is required before running the simulator, unless file permissions need to be modified (e.g., in Linux, type "chmod +x FILENAME" in a terminal to enable execution of FILENAME). It is recommended to use a "config" folder to store configuration files and a "results" folder to store simulation output. AcCoRD will detect these folders automatically if they are in the working directory or are subdirectories of the working directory's parent.

Windows binaries were compiled using GCC 4.9.3 on minGW and the C99 standard of C.
Debian/Ubuntu binaries were compiled using GCC 4.8.4 and the C99 standard of C.
RHEL/CentOS binaries were compiled using GCC 4.8.5 and the C99 standard of C.

It is possible that you may need C libraries on your system that correspond to the compiler used to generate the pre-compiled binaries. If you have trouble getting the pre-compiled binaries to work, then it might be easier to try the advanced installation.

2) (ADVANCED) Build from source code. This is needed if you want to use different compilation parameters. The src directory contains scripts for Windows, Debian/Ubuntu, and RHEL/CentOS builds (note that the linux build scripts are identical except for the binary filenames in order to distinguish between them). Run a script from the command line while in the "src" directory and the binary will be placed in the "bin" directory. GCC and standard C libraries are required.
* build_accord_opt_win.bat builds the optimized Windows version with executable accord_win.exe
* build_accord_debug_win.bat builds the debug Windows version with executable accord_win_debug.exe
* build_accord_opt_dub builds the optimized Debian/Ubuntu version with executable accord_dub.out
* build_accord_debug_dub builds the debug Debian/Ubuntu version with executable accord_dub_debug.out
* build_accord_opt_rc builds the optimized RHEL/CentOS version with executable accord_rc.out
* build_accord_debug_rc builds the debug RHEL/CentOS version with executable accord_rc_debug.out

If you are compiling on Windows, then minGW is recommended for the GCC compiler and related utilities http://www.mingw.org/

You may need to change the file permissions to execute the build script in Linux (e.g., chmod +x FILENAME).

Example calls for compiling the optimized version:
From Windows command line: build_accord_opt_win.bat
From shell in Debian/Ubuntu: ./build_accord_opt_dub


# Basic Usage

The code release page and root source directory both include "fake" config and output files (they include comments so they are not valid JSON files and cannot be used as-is). They are called "HOWTO_DEFINE_CONFIG.txt", "HOWTO_READ_SUMMARY_OUTPUT.txt", and "HOWTO_READ_OUTPUT.txt". The config subdirectory also includes a number of sample valid configuration files for you to start running and experimenting with.

AcCoRD is run from the command line. It takes 2 additional (optional) arguments when called:  
1) Configuration filename. Config file must be defined relative to one of 3 locations. AcCoRD will first search the current directory, then the "config" subdirectory, and finally the "../config/" directory. There are a number of sample configuration files provided to demonstrate AcCoRD functionality.  
2) Seed value for random number generator (will override value specified in config file). You can read more about the meaning of the seed value in HOWTO_DEFINE_CONFIG.txt

If there are no additional arguments, then a default config file is used ("accord_config_sample.txt").
If there is one additional argument, then it must be the configuration filename.

Sample call from Windows command prompt (where both the executable and the configuration file are in the current directory):
accord_win.exe myconfig.txt 2
Sample call from Ubuntu/Debian shell while in the AcCoRD "bin" directory while the configuration file is in "../config/":
./accord_dub.out myconfig.txt 2

AcCoRD will create 2 files. The main output file has the suffix '_SEEDX.txt', where X is the seed number used, and a summary file will have the suffix '_SEEDX_summary.txt'. These 2 files should be kept together, especially if you want to import them into MATLAB.

To import simulation output as a structure in MATLAB (version R2015a or newer recommended), use the accordImport function found in the matlab folder of the source code. Please note that the files in the "JSONLab" subdirectory are required for this function to work. The call is:
[data, config] = accordImport(FILENAME, SEEDRANGE, bWrite)

where
* FILENAME - the output filename to load (including the relative path, if applicable, but excluding the '_SEEDX.txt' suffix)
* SEEDRANGE - a vector specifying the seed values to import (each seed value corresponds to one pair of output files)
* bWrite - switch to write the "data" and "config" structures to a MATLAB mat-file named FILENAME_out.mat in the "matlab" directory.
* data - output structure of simulation data
* config - output structure with the parameters of the configuration file used to run the simulation

accordImport will use the summary file to learn the name of the original configuration file and will search for that file in the current directory, then in the "config", "../config/", "../", "../../", and "../../config" directories (and in that order). If a file with a matching configuration filename cannot be found, then "config" will be an empty array.

The "config" structure will contain all of the parameters specified in the original configuration file, using a format similar to that in the configuration file.
The "data" structure will contain the data from both the output and output summary files. The names of the structure members are similar to those used in the output files themselves, so reading the structure should be straightforward.
Examples:
* data.numRepeat -> total number of independent realizations (aggregated from all seed values).
* data.activeID(i) -> index (in the global actor list) of the ith actor in the active actor list. The global list includes both active and passive actors in no particular order.
* data.activeBits{i}(j,k) -> value of the kth bit sent by the ith active actor (indexed from the active actor list) in the jth simulation realization.
* data.passiveRecordMolID{i}(j) -> index (in the global molecule list) of the the jth molecule being observed by the ith recorded passive actor.
* data.passiveRecordCount{i}(j,k,l) -> number of molecules of the kth type observed by the ith recorded passive actor in the lth observation of the jth simulation realization.

If you have the original configuration file accessible and one of the passive actors recorded molecule positions, then you can watch the simulation progress by generating a video or sequence of figures. To get started, import the simulation output using accordImport, then copy and rename accordVideoMakerWrapper. This wrapper function creates all of the arguments needed for the accordVideoMaker function. You should modify the copy for your particular simulation. There are dozens of options that can be tweaked (practically any Matlab option for plotting figures, axes, patch objects, and markers, plus deciding what to display and at what times). When you run the copy of accordVideoMakerWrapper, accordVideoMaker will pause before making a video or sequence of figures so that you can manually tweak the camera and axes settings as desired.

To plot the simulation environment in MATLAB without needing to simulate (i.e., to plot just the regions or actors without showing molecules), then use a copy of accordEmptyEnvironmentWrapper. This function calls accordEmptyEnvironment, which loads the specified configuration file.

To plot passive actor observation data in MATLAB, then use a copy of accordPlotMakerWrapper, which calls accordPlotMaker. Complete details of the kinds of plots and how to define them can be found in the header of accordBuildObserverStruct. Options include time-varying signals, probability mass and cumulative distribution functions, and mutual information. Plots also have "Expected" versions that are based on Binominal, Poisson, or Gaussian distributions.


# Licensing

Main AcCoRD files are copyright 2016 by Adam Noel under the "New BSD" license. For full details of licensing, including use of third-party code, please see LICENSE.txt

# Credits

Developed and maintained by Adam Noel (http://www.adamnoel.ca)

Testing: M. Halimeh and T. Schwering

Supervision and Support: Profs. K. C. Cheung, A. Hafid, D. Makrakis, and R. Schober.