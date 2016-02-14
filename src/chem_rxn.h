/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * chem_rxn.h - structure for storing chemical reaction properties
 * Last revised for AcCoRD v0.4
 *
 * Revision history:
 *
 * Revision v0.4
 * - modified use of unsigned long and unsigned long long to uint32_t and uint64_t
 * - modified propensity updates to do a full re-calculation in order to avoid
 * numerical underflow problems
 * - removed deprecated debug function
 * - added restriction of chemical reactions to specific regions
 * - renamed region structure and array of structures to accommodate naming
 * of new boundary structure
 *
 * Revision v0.3.1
 * - header added
*/
#ifndef CHEM_RXN_H
#define CHEM_RXN_H

#include <math.h> // for exp()

#include "region.h"
#include "global_param.h"

/*
* Data Type Declarations
*/
//struct region; // Forward declaration (declared in region.h)

/* The chem_rxn structure defines the parameters relevant to the
 * firing of a single chemical reaction
*/
struct chem_rxn_struct { // Used to define a single chemical reaction
	// Indicate the indices of the reactants and the number of each
	// Only the first NUM_MOL_TYPES values will be defined
	uint32_t reactants[MAX_MOL_TYPES];
	
	// Indicate the indices of the products and the number of each
	// Only the first NUM_MOL_TYPES values will be defined
	uint32_t products[MAX_MOL_TYPES];
	
	// Base reaction rate k (units depends on order of reaction)
	double k;
	
	// Can the reaction take place anywhere by default?
	bool bEverywhere;
	
	// Number of regions that are exceptions to the location default
	short numRegionExceptions;
	
	// Names of regions that are exceptions to the location default
	// Length is numRegionExceptions
	char ** regionExceptionLabel;
	
	// TODO: Add parameters for reactions that take place across multiple
	// regions (i.e., surface interactions)
};

//
// Function Declarations
//

/* Allocate space for the regionArray members needed for chemical
 * reactions and initialize their values based on the chem_rxn
 * array of chem_rxn_struct
*/
void initialize_region_chem_rxn3D(const short NUM_REGIONS,
	struct region regionArray[],
	const unsigned short NUM_MOL_TYPES,
	const unsigned short MAX_RXNS,
	const struct chem_rxn_struct chem_rxn[]);

void delete_region_chem_rxn3D(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	struct region regionArray[]);

#endif // CHEM_RXN_H
