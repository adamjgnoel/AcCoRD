
            The AcCoRD Simulator
            (Actor-based Communication via Reaction-Diffusion)

This document is the README for AcCoRD v1.4.2 (2020-02-12)
Last updated 2024-11-26

# Introduction

Welcome to the AcCoRD Simulator (Actor-based Communication via Reaction-Diffusion). AcCoRD is a simulation tool for molecular communication. It is a hybrid microscopic/mesoscopic simulator that functions as a generic solver of stochastic advection-reaction-diffusion but which has been developed from the perspective of communications analysis. The focus is to efficiently generate the statistics of molecule observations at some location (i.e., at a receiver).

Quicklinks:
* User-friendly external webpage (includes links to feature list, download page, and instructions for use): https://www.engr.mun.ca/~adamnoel/accord.html
* Old AcCoRD releases: https://github.com/adamjgnoel/AcCoRD/releases
* Main AcCoRD development page on Github: https://github.com/adamjgnoel/accord
* Maintenance of existing issues and planned features: https://github.com/adamjgnoel/AcCoRD/issues
* Github wiki (just a copy of this README): https://github.com/adamjgnoel/AcCoRD/wiki

Primary features include:
* As a hybrid solver, AcCoRD integrates two simulation models to define 3D environments with flexible local accuracy.  Microscopic regions define each molecule individually and evolve over discrete time steps. Mesoscopic regions count the number of molecules in disjoint virtual bins (called subvolumes) and evolve over time steps with continuous granularity.
* "Actors" are either active molecule sources (i.e., transmitters) or passive observers (i.e., receivers). Active actors can modulate binary sequence data and release molecules over specified intervals. Passive actors record the number of molecules observed and can also record molecule positions.
* Chemical reactions can be defined locally or globally. A general framework for defining reactions can accommodate reactions such as molecule degradation, enzyme kinetics, reversible or irreversible surface binding, ligand-receptor binding, transitions across boundary membranes, and simplified molecular crowding. Surface binding reactions include absorption (consumption), adsorption (sticking), and desorption (release).
* Independent realizations of a simulation can be repeated an arbitrary number of times (and on different computers) and then aggregated to determine the average behavior and channel statistics.
* AcCoRD's interface has been designed to be helpful to the novice user while providing meaningful output. AcCoRD also includes post-processing tools developed in MATLAB. These tools enable saving aggregated simulation output files, plotting receiver observations (either the time-varying behavior or empirical distributions at specified times), and visualizing the physical environment (either as still images or compiled into a video)

## Audience

If you are a communications researcher who is interested in studying molecular communication, then you are the intended audience. However, anyone interested in the multi-scale simulation of reaction-diffusion systems might also find AcCoRD useful.

## Development Status

Development of AcCoRD is on-going and hosted on Github (https://github.com/adamjgnoel/AcCoRD).

Existing bugs and planned features are described on the Issues page: https://github.com/adamjgnoel/AcCoRD/issues

Most identified content on the Issues page are planned features.

# Publications

Publications that either describe or use AcCoRD can be found at this page: https://www.engr.mun.ca/~adamnoel/accord.html

# Known Issues and Limitations

This list is current as of v1.3. See https://github.com/adamjgnoel/AcCoRD/issues for the latest details.
* Full 2D simulations are mostly untested
* 2nd order chemical reactions in microscopic regions need a binding radius specified.
* There is an option to set whether actor is independent but the value is ignored since dependent actors have not yet been implemented
* 2 regions must share an outer face in order to be classified as neighbors, unless one region is the parent of the other. So, if a child region is up against a face of its parent, and that parent also has a parent but doesn't share the same face with its parent, then the parent's parent will not be recognized as a neighbor of the child. In such a case molecule transitions between the child and the "grandparent" are not possible.
* Surface regions must be microscopic and can only have microscopic regions as neighbors.
* A child region is always assumed to be a neighbor of its parent, which is not strictly true but should not create a problem during a simulation
* microscopic molecules in a spherical region may go a very small distance beyond the region boundary but are still believe to be inside the sphere. This is due to numerical underflow, and might mean that the potential chemical reactions for the molecule are inconsistent with its location.
* A mesoscopic subvolume adjacent to a microscopic region needs that region to cover the entirety of every face that is adjacent to that region (i.e., part of the adjacent face should not also be adjacent to some other region). Otherwise the hybrid interface will not have proper transitions.
* Actors that act at the same time will do so in a random order (this is by design). So, if an active actor is adding molecules at the exact time that a passive actor is supposed to observe them, it is unknown a priori which will act first. It will depend on whichever is currently ranked higher in the timer heap. This could be modified in the future to rank the actors in the order they are defined in the event that their timer values are the same.
* Unbinding radii are only applied to second order microscopic reactions that have multiple products.


# Future Features

Development of AcCoRD is on-going. See https://github.com/adamjgnoel/AcCoRD/issues for the latest details. Here is a list of some planned end-user features (in no particular order):
* Expand flow model. Currently, flow must be steady and uniform (within a given region).
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

Download and installation instructions can be found on this page: https://www.engr.mun.ca/~adamnoel/accord.html


# Basic Usage

Usage instructions can be found in the User's Manual that is included with the software download


# Licensing

Main AcCoRD files are copyright 2016-2020 by Adam Noel under the "New BSD" license. For full details of licensing, including use of third-party code, please see LICENSE.txt

# Credits

Developed and maintained by Adam Noel (https://www.engr.mun.ca/~adamnoel)

Testing: M. Halimeh (v0.2 to v0.4) and T. Schwering (v0.3 to v0.5)

Supervision and Support: K. C. Cheung (to v1.0), A. Hafid (v0.4 to v1.0), D. Makrakis (v0.4 to v1.1), and R. Schober  (to v1.0).
