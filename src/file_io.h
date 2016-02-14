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
 * Last revised for AcCoRD v0.4
 *
 * Revision history:
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
#include <stdlib.h> // for exit(), malloc
#include <inttypes.h> // for extended integer type macros
#include <stdbool.h> // for C++ bool naming, requires C99
#include <string.h> // for strcpy()
#include <time.h> // For time record keeping
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
	
	// Environment
	unsigned short NUM_DIM;
	double SUBVOL_BASE_SIZE;
	short NUM_REGIONS;
	struct spec_region3D * subvol_spec;
	short NUM_ACTORS;
	struct actorStructSpec3D * actorSpec;
	
	// Chemical Properties
	unsigned short NUM_MOL_TYPES;
	unsigned short MAX_RXNS;
	double * DIFF_COEF;
	struct chem_rxn_struct * chem_rxn;
};

//
// Function Declarations
//

void loadConfig3D(const char * CONFIG_NAME,
	uint32_t customSEED,
	struct simSpec3D * curSpec);

void deleteConfig3D(struct simSpec3D curSpec);

void initializeOutput3D(FILE ** out,
	FILE ** outSummary,
	const char * CONFIG_NAME,
	const struct simSpec3D curSpec);

char * stringAllocate(long stringLength);

void printOneTextRealization3D(FILE * out,
	const struct simSpec3D curSpec,
	unsigned int curRepeat,
	ListObs3D observationArray[],
	short numActorRecord,
	short * actorRecordID,
	short NUM_ACTORS_ACTIVE,
	const struct actorStruct3D actorCommonArray[],
	const struct actorActiveStruct3D actorActiveArray[],
	const struct actorPassiveStruct3D actorPassiveArray[],	
	uint32_t maxActiveBits[],
	uint32_t maxPassiveObs[]);
	
void printTextEnd3D(FILE * out,	
	short NUM_ACTORS_ACTIVE,
	short numActorRecord,
	const struct actorStruct3D actorCommonArray[],
	const struct actorActiveStruct3D actorActiveArray[],
	const struct actorPassiveStruct3D actorPassiveArray[],
	short * actorRecordID,
	uint32_t maxActiveBits[],
	uint32_t maxPassiveObs[]);

#endif // FILE_IO_H