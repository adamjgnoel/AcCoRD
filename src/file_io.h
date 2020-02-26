/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * file_io.h - interface with JSON configuration files
 *
 * Last revised for AcCoRD v1.4.2 (2020-02-12)
 *
 * Revision history:
 *
 * Revision v1.4.2 (2020-02-12)
 * - updated file output file directory creation to accommodate Mac OS
 *
 * Revision v1.4 (2018-08-06)
 * - added a priori monte carlo (APMC) absorption algorithm as a new surface
 * reaction type. Includes settings for how to define the a priori absorption
 * probability calculation and whether/how to apply a threshold to turn it off
 *
 * Revision v1.3 (2018-07-31)
 * - changed checks for the "Random Molecule Release Times?" "Slot Interval" parameters
 * for active actors, since whether they are needed depends on the values of other parameters.
 *
 * Revision v1.1 (2016-12-24)
 * - added flow parameters (based on either no flow or uniform flow).
 * Global settings apply to any subset of molecule types but can be modified for any
 * individual molecule type in any region.
 * - corrected default subvolume base size if is is not defined correctly
 * - corrected detection of missing surface reaction probability for surface reactions
 *
 * Revision v1.0 (2016-10-31)
 * - added warnings for defining passive actor parameters for active actors or
 * vice versa.
 * - added "BURST" modulation, which does not rely on bits and always releases molecules
 * in every action interval. Also is able to release multiple types of molecules.
 * - added simpler methods for defining region anchor coordinates and the number of
 * subvolumes along each dimension of rectangular regions.
 * - added local diffusion coefficients as region parameters. Any region can
 * over-ride the default diffusion coefficients defined for all molecules
 * - added specifying diffusion coefficient to use for surface reaction transition
 * probabilities. By default, the reactant's default diffusion coefficient is used
 * for absorption and membrane reactions, and the first product's default diffusion
 * coefficient is used for desorption reactions. Warning appears if coefficient is
 * defined for reaction types that cannot use it.
 *
 * Revision v0.7.0.1 (public beta, 2016-08-30)
 * - added measurement of simulation runtime to be written to simulation output
 * - fixed bug in writing index of active actor in simulation summary
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - added option for user to select small subvolume or big subvolume assumption
 * at hybrid interfaces. Only needed if there is at least 1 microscopic and 1
 * mesoscopic region defined.
 * - added option for user to choose maximum distance beyond which a micro to mesoscopic
 * substep transition will not be considered (either initial or final point should be beyond
 * this distance)
 * - added binding and unbinding radii as bimolecular reaction parameters (needed in
 * microscopic regime)
 * - made output of active actor data sequence a user option
 * - added option for user to define a constant active actor bit sequence
 * - added point active actors defined by their 3D coordinate.
 * - added check for positive (and not just non-negative) radii for spherical regions and
 * actors (i.e., 0 is not a valid radius)
 * - added warnings for unnecessary active actor parameters depending on values
 * of other active actor parameters
 *
 * Revision v0.5.1 (2016-05-06)
 * - added bReleaseProduct to chemical reaction. Applies to surface reactions
 * - added chemical reaction properties to define coupled reversible reactions
 * - added chemical reaction properties to configure absorbing, desorbing, and membrane
 * reactions, including how transition probabilities are calculated and how
 * desorbed molecules are placed. Membrane reactions must have no product molecules
 * specified (product is always the same as the reactant)
 * - shortened string used to indicate how actor location is defined. This avoids
 * a MATLAB warning when loaded simulation configuration files
 *
 * Revision v0.5 (2016-04-15)
 * - added ability to define location of actor by a list of regions
 * - modified check on number of subvolumes along each dimension of a rectangular region
 * - added type and surfaceType properties to region. Default values are REGION_NORMAL and
 * NO_SURFACE, respectively
 * - added bSurface and surfRxnType properties to chemical reaction. Default values are false and RXN_NORMAL, respectively.
 * - added stringWrite function to nest some calls to stringAllocate, strlen,
 * and strcpy
 * - removed NUM_DIM parameter from simulation spec
 * - removed upper limit on number of molecule types
 *
 * Revision v0.4.1
 * - added search for configuration file. First checks current directory, then
 * config subdirectory, and then config sibling directory
 * - added search for results folder in current directory and then as subdirectory
 * of parent. If it cannot be found, results folder is created in current directory
 * and output is placed there.
 * - improved use and format of error messages
 *
 * Revision v0.4
 * - modified use of unsigned long and unsigned long long to uint32_t and uint64_t
 * - added options to accommodate spherical regions and actors
 * - added restriction of chemical reactions to specific regions
 * - renamed region structure and array of structures to accommodate naming
 * of new boundary structure
 *
 * Revision v0.3.1.1
 * - corrected warning message when loading a configuration file to display the correct
 *   keys to quit or continue
 *
 * Revision v0.3.1
 * - check all JSON members for existence and correct format before reading them in
 * - header added
*/

#ifndef FILE_IO_H
#define FILE_IO_H

#include <stdio.h>
#include <sys/stat.h> // for stat(), S_ISDIR(), mkdir()
#include <stdlib.h> // for exit(), malloc
#include <inttypes.h> // for extended integer type macros
#include <stdbool.h> // for C++ bool naming, requires C99
#include <string.h> // for strcpy()
#include <time.h> // For time record keeping
#ifndef __linux__
	#include <direct.h> // for _mkdir() [Windows]
#endif // __linux__
#include "cJSON.h"
#include "region.h"
#include "actor.h"
#include "chem_rxn.h"
#include "micro_molecule.h" // for individual molecule definitions, operations
#include "actor_data.h" // for active actor binary data
#include "observations.h" // for observation structure (linked list)
#include "global_param.h" // for common global parameters

//
// Data type declarations
//

struct simSpec3D {
	char * OUTPUT_NAME;
	
	// Simulation Control
	unsigned int NUM_REPEAT;
	double TIME_FINAL;
	double DT_MICRO;
	uint32_t SEED;
	unsigned int MAX_UPDATES;
	bool B_HYBRID_SMALL_SUB;
	double MAX_HYBRID_DIST;
	
	// Environment
	double SUBVOL_BASE_SIZE;
	short NUM_REGIONS;
	struct spec_region3D * subvol_spec;
	short NUM_ACTORS;
	struct actorStructSpec3D * actorSpec;
	
	// Chemical Properties
	unsigned short NUM_MOL_TYPES;
	unsigned short MAX_RXNS;
	double * DIFF_COEF;
	unsigned short GLOBAL_FLOW_TYPE;
	double * GLOBAL_FLOW_VECTOR;
	bool * B_GLOBAL_MOL_FLOW;
	struct chem_rxn_struct * chem_rxn;
};

//
// Function Declarations
//

void loadConfig(const char * CONFIG_NAME,
	uint32_t customSEED,
	struct simSpec3D * curSpec);

void deleteConfig(struct simSpec3D curSpec);

void initializeOutput(FILE ** out,
	FILE ** outSummary,
	const char * CONFIG_NAME,
	const struct simSpec3D curSpec);

// Copy string (with memory allocation)
char * stringWrite(char * src);

// Allocate memory for a string
char * stringAllocate(long stringLength);

void printOneTextRealization(FILE * out,
	const struct simSpec3D curSpec,
	unsigned int curRepeat,
	ListObs3D observationArray[],
	short numPassiveRecord,
	short * passiveRecordID,
	short numActiveRecord,
	short * activeRecordID,
	const struct actorStruct3D actorCommonArray[],
	const struct actorActiveStruct3D actorActiveArray[],
	const struct actorPassiveStruct3D actorPassiveArray[],	
	uint32_t maxActiveBits[],
	uint32_t maxPassiveObs[]);
	
void printTextEnd(FILE * out,	
	short numActiveRecord,
	short numActorRecord,
	const struct actorStruct3D actorCommonArray[],
	const struct actorActiveStruct3D actorActiveArray[],
	const struct actorPassiveStruct3D actorPassiveArray[],
	short * passiveRecordID,
	short * activeRecordID,
	uint32_t maxActiveBits[],
	uint32_t maxPassiveObs[],
	double runTime);

#endif // FILE_IO_H