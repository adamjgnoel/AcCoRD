
            The AcCoRD Simulator
            (Actor-based Communication via Reaction-Diffusion)

This document is the README for AcCoRD v0.5 (public beta, 2016-04-15)

# Introduction

Welcome to the AcCoRD Simulator (Actor-based Communication via Reaction-Diffusion). AcCoRD is a simulation tool for molecular communication. It is a hybrid microscopic/mesoscopic simulator that functions as a generic solver of stochastic reaction-diffusion but which has been developed from the perspective of communications analysis. The focus is to efficiently generate the statistics of molecule observations at some location (i.e., at a receiver).

As of v0.5 (2016-04-15), this document serves as both the official README (https://github.com/adamjgnoel/AcCoRD/blob/master/README.txt) and the Github wiki page (https://github.com/adamjgnoel/AcCoRD/wiki)


## Audience

If you are a communications researcher who is interested in studying molecular communication, then you are the intended audience. However, anyone interested in the multi-scale simulation of reaction-diffusion systems might also find AcCoRD useful.

## Development Status

Development of AcCoRD is on-going and hosted on Github (https://github.com/adamjgnoel/AcCoRD).

The executables (for Windows, Debian/Ubuntu, and RHEL/CentOS) are available here: https://github.com/adamjgnoel/AcCoRD/releases

Existing bugs and planned features are described on the Issues page: https://github.com/adamjgnoel/AcCoRD/issues

Most identified "bugs" in the issues page are either:

1) Parameters ignored in configuration files because they haven't been implemented yet  
2) Cases where the simulator does not recognize errors in the configuration file. These cases have been greatly reduced as of v0.5.

# Background

AcCoRD development began in October 2014 with the motivation to develop a reaction-diffusion simulation tool that would be suitable for the communications research community. The goal was to integrate recent advances in stochastic reaction-diffusion models. The overall design of AcCoRD is very much inspired by Smoldyn (http://www.smoldyn.org/), which is also a microscopic reaction-diffusion solver (and which, coincidentally, has also recently added mesoscopic regions). The key differences with Smoldyn relate to AcCoRD's communication focus. AcCoRD's intended use is to capture the dynamics of a communications system, i.e., with a transmitter and receiver. Transmitters are implemented as "active actors", which manually add molecules to the environment. Observations by receivers are recorded as "passive actors". The placement of actors is defined by the user, and the only information that AcCoRD saves is that related to actor behavior. Furthermore, simulations are designed to be run repeatedly (i.e., an arbitrary number of times, but could easily be tens to hundreds of thousands of times), so that accurate receiver statistics can be generated.

AcCoRD is developed in C in order to have precise control over memory management and access to very fast algorithms for random number generation, linked list operations, and related computational tasks. Some utilities have been written (and will be written) in MATLAB to facilitate post-processing.

While a formal publication on the design of AcCoRD is in preparation, the general motivation and some of the implementation ideas were presented in the following conference papers:  
- A. Noel, K.C. Cheung, and R. Schober, On the Statistics of Reaction-Diffusion Simulations for Molecular Communication, in Proc. ACM NANOCOM 2015, Sep. 2015. DOI: http://dx.doi.org/10.1145/2800795.2800821 - the simulations in this paper were completed with proof-of-concept MATLAB code  
- A. Noel, K.C. Cheung, and R. Schober, Multi-Scale Stochastic Simulation for Diffusive Molecular Communication, in Proc. IEEE ICC 2015, pp. 2712--2718, Jun. 2015. DOI: http://dx.doi.org/10.1109/ICC.2015.7248471 - the simulations in this paper were completed with AcCoRD v0.1, which was an early 2D build

Additional documentation will be prepared as AcCoRD nears general public release. Until that time, please contact the developer or consider following AcCoRD on Github (https://github.com/adamjgnoel/accord) for the latest updates.

# Feature Summary

AcCoRD currently has the following "mature" features:
* 3D multi-scale hybrid reaction-diffusion
* Environment described as union of cubes and spheres. Cubes can be microscopic (track individual molecules) or mesoscopic (assume uniform molecule density throughout the cube). Spheres must be microscopic but can be infinite in size (i.e., unbounded environment).
* Chemical reactions of 0th, 1st, and 2nd order in mesoscopic regime. No 2nd order reactions permitted in microscopic regime (for now).
* Readable configuration and simulation output files. Configuration is described using the JSON interchange format. An import function can load output into MATLAB.
* Flexible placement of molecule sources and observers, whose location can be defined as their own shapes or as the union of a set of regions
* Molecule sources based on "Concentration Shift Keying" (CSK) modulation of independent symbols, with user-defined pulse widths and switches for stochastic or deterministic molecule release
* Molecule observers record number of specified types of molecules (and optionally their positions) over a subset of the simulation environment.
* Even in the microscopic regime, molecule creation and chemical reactions can occur on a time scale smaller than the microscopic time step.

The following features have not been extensively tested and are considered in "beta":
* Surfaces can surround all or part of 3D regions and be sites for local chemical reactions, such as molecule generation or absorption.
* Membrane surfaces control transitions between neighboring 3D regions

The following features have had very limited testing and are considered "alpha":
* Environments described entirely in 2D, with diffusion across a union of 2D rectangles

A complete version history can be found in CHANGELOG.txt

# Known Issues

This list is current as of v0.5. See https://github.com/adamjgnoel/AcCoRD/issues for the latest details.
* Simplified algorithms are used to handle molecule placement when transitioning between mesoscopic and microscopic regions
* Chemical reactions in microscopic regions must be 0th or 1st order.
* There is an option to set whether actor is independent but the value is ignored since dependent actors have not yet been implemented
* There is an option to set whether active actor has a message of random bits but the value is ignored so that the bits are always random (though they can be set to all 1s or all 0s by setting probability to 1 or 0, respectively)
* Active actor modulation scheme is presented as a configuration option but the only current valid choice is "CSK" (concentration shift keying)
* The bit stream of active actors will always be written to the output file, independent of the setting of "Is Actor Activity Recorded?". This is because active actor emission is currently always random.
* Point sources are not accommodated
* 2 regions must share an outer face in order to be classified as neighbors, unless one region is the parent of the other. So, if a child region is up against a face of its parent, and that parent also has a parent but doesn't share the same face with its parent, then the parent's parent will not be recognized as a neighbor of the child. In such a case molecule transitions between the child and the "grandparent" are not possible.
* A child region is always assumed to be a neighbor of its parent, which is not strictly true but should not create a problem during a simulation
* microscopic molecules in a spherical region may go a very small distance beyond the region boundary but are still believe to be inside the sphere. This is due to numerical underflow, and might mean that the potential chemical reactions for the molecule are inconsistent with its location.
* A mesoscopic subvolume adjacent to a microscopic region needs that region to cover the entirety of every face that is adjacent to that region (i.e., part of the adjacent face should not also be adjacent to some other region). Otherwise the hybrid interface will not have proper transitions.
* Actors that act at the same time will do so in a random order (this is by design). So, if an active actor is adding molecules at the exact time that a passive actor is supposed to observe them, it is unknown a priori which will act first. It will depend on whichever is currently ranked higher in the timer heap. This could be modified in the future to rank the actors in the order they are defined in the event that their timer values are the same.


# Future Features

Development of AcCoRD is active and on-going. See https://github.com/adamjgnoel/AcCoRD/issues for the latest details. Here is a list of planned end-user features in approximate order of priority:
* Display warnings if user defines active actor parameters for a passive actor (or vice versa)
* Allow dependent actors, whose behavior depends on the observations of passive actors. This would allow implementation of more complex communication networks (e.g., relays)
* Add actors that track changes to molecule composition in their volumes
* Add 2nd order chemical reactions in microscopic regions
* Improve accuracy of transitions between microscopic and mesoscopic regions
* Enable surfaces to mesoscopic regions
* Write more MATLAB utility functions to perform common tasks with simulation output (e.g., plot average behavior, simulation statistics, bit error probabilities, and video of simulation progression)
* Specify different diffusion coefficients for different regions
* Use different microscopic time steps in different microscopic regions
* Add checks on subvolume size relative to diffusion coefficients
* Add "tau leaping" in mesoscopic regime to improve scalability
* Allow user to specify the bits of an active actor transmission sequence
* Add steady uniform flow


# Installation

There are two installation options:

1) (BASIC) Use a pre-compiled binary from the releases page (https://github.com/adamjgnoel/AcCoRD/releases).
Debug and optimized versions exist for Windows, Debian/Ubuntu Linux, and RHEL/CentOS Linux. No further action is required before running the simulator, unless file permissions need to be modified (e.g., in Linux, type "chmod +x FILENAME" in a terminal to enable execution of FILENAME). It is recommended to use a "config" folder to store configuration files and a "results" folder to store simulation output. AcCoRD will detect these folders automatically if they are in the working directory or are subdirectories of the working directory's parent.

Windows binaries were compiled using GCC 4.9.3 on minGW and the C99 standard of C.
Debian/Ubuntu binaries were compiled using GCC 4.8.4 and the C99 standard of C.
RHEL/CentOS binaries were compiled using GCC 4.8.5 and the C99 standard of C.

It is possible that you may need C libraries on your system that correspond to the compiler used to generate the pre-compiled binaries. If you have trouble getting the pre-compiled binaries to work, then it might be easier to try the advanced installation.

2) (ADVANCED) Build from source code. This is needed if you want to use different compilation
parameters. The src directory contains scripts for Windows, Debian/Ubuntu, and RHEL/CentOS builds (note that the linux build scripts are identical except for the binary filenames in order to distinguish between them). Run a script from
the command line while in the "src" directory and the binary will be placed in the "bin" directory. GCC and standard C libraries are required.
* build_accord_opt_win.bat builds the optimized Windows version with executable accord_win.exe
* build_accord_debug_win.bat builds the debug Windows version with executable accord_win_debug.exe
* build_accord_opt_dub builds the optimized Debian/Ubuntu version with executable accord_dub.out
* build_accord_debug_dub builds the debug Debian/Ubuntu version with executable accord_dub_debug.out
* build_accord_opt_rc builds the optimized RHEL/CentOS version with executable accord_rc.out
* build_accord_debug_rc builds the debug RHEL/CentOS version with executable accord_rc_debug.out

If are compiling on Windows, then minGW is recommended for GCC http://www.mingw.org/

You may need to change the file permissions to execute the build script in Linux (e.g., chmod +x FILENAME).

Example calls for compiling the optimized version:
From Windows command line: build_accord_opt_win.bat
From shell in Debian/Ubuntu: ./build_accord_opt_dub


# Basic Usage

The root source directory includes "fake" config and output files (they include comments so they are not valid JSON files and cannot be used as-is). They are called "HOWTO_DEFINE_CONFIG.txt", "HOWTO_READ_SUMMARY_OUTPUT.txt", and "HOWTO_READ_OUTPUT.txt". The config subdirectory also includes a number of sample valid configuration files for you to start running and experimenting with.

AcCoRD is run from the command line. It takes 2 additional (optional) arguments when called:  
1) Configuration filename. Config file must be defined relative to one of 3 locations. AcCoRD will first search the current directory, then the "config" subdirectory, and finally the "../config/" directory. There are a number of sample configuration files provided to demonstrate AcCoRD functionality.  
2) Seed value for random number generator (will override value specified in config file). You can read more about the meaning of the seed value in HOWTO_DEFINE_CONFIG.txt

If there are no additional arguments, then a default config file is used ("accord_config_sample.txt").
If there is one additional argument, then it must be the configuration filename.

Sample call from Windows command prompt (where both the executable and the configuration file are in the current directory):
accord_win.exe myconfig.txt 2
Sample call from Linux shell while in the AcCoRD "bin" directory while the configuration file is in "../config/":
./accord myconfig.txt 2

AcCoRD will create 2 files. The main output file has the suffix '_SEEDX.txt', where X is the seed number used, and a summary file will have the suffix '_SEEDX_summary.txt'. These 2 files should be kept together, especially if you want to import them into MATLAB.

To import simulation output as a structure in MATLAB, use the accord_import function found in the matlab folder of the source code. Please note that the files in the "JSONLab" subdirectory are required for this function to work. The call is:
[data, config] = accord_import(FILENAME, SEEDRANGE, bWrite)

where
FILENAME - the output filename to load (including the relative path, if applicable, but excluding the '_SEEDX.txt' suffix)
SEEDRANGE - a vector specifying the seed values to import (each seed value corresponds to one pair of output files)
bWrite - switch to write the "data" and "config" structures to a MATLAB mat-file named CONFIG_NAME_out.mat in the "matlab" directory, where CONFIG_NAME is the name of the configuration file that was originally used to run the simulation.
data - output structure of simulation data
config - output structure with the parameters of the configuration file used to run the simulation

accord_import will use the summary file to learn the name of the original configuration file and will search for that file in the current directory, then in the "config", "../config/", "../", "../../", and "../../config" directories (and in that order). If a file with a matching configuration filename cannot be found, then "config" will be an empty array.

The "config" structure will contain all of the parameters specified in the original configuration file, using a format similar to that in the configuration file.
The "data" structure will contain the data from both the output and output summary files. The names of the structure members are similar to those used in the output files themselves, so reading the structure should be straightforward.
Examples:
* data.numRepeat -> total number of independent realizations (aggregated from all seed values).
* data.activeID(i) -> index (in the global actor list) of the ith actor in the active actor list. The global list includes both active and passive actors in no particular order.
* data.activeBits{i}(j,k) -> value of the kth bit sent by the ith active actor (indexed from the active actor list) in the jth simulation realization.
* data.passiveRecordMolID{i}(j) -> index (in the global molecule list) of the the jth molecule being observed by the ith recorded passive actor.
* data.passiveRecordCount{i}(j,k,l) -> number of molecules of the kth type observed by the ith recorded passive actor in the lth observation of the jth simulation realization.

A future release of AcCoRD will include more utilities for post-processing simulation results.


# Licensing

Main AcCoRD files are copyright 2016 by Adam Noel under the "New BSD" license. For full details of licensing, including use of third-party code, please see LICENSE.txt

# Credits

Developed and maintained by Adam Noel (http://www.adamnoel.ca)

Testing: M. Halimeh and T. Schwering

Supervision and Support: Profs. K. C. Cheung, A. Hafid, D. Makrakis, and R. Schober.