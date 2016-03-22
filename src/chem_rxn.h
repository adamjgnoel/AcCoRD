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
 * Last revised for AcCoRD v0.4.1
 *
 * Revision history:
 *
 * Revision v0.4.1
 * - improved use and format of error messages
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
	uint32_t * reactants;
	
	// Indicate the indices of the products and the number of each
	uint32_t * products;
	
	// Base reaction rate k (units depends on order of reaction)
	double k;
	
	// Is the reaction a surface reaction?
	// This will affect where a reaction will take place by default,
	// as indicated by bEverywhere
	bool bSurface;
	
	// Can the reaction take place anywhere by default?
	// Actual regions will depend on value of bSurface and whether a given
	// region is a normal region or a surface
	bool bEverywhere;
	
	// Number of regions that are exceptions to the location default
	short numRegionExceptions;
	
	// Names of regions that are exceptions to the location default
	// Length is numRegionExceptions
	char ** regionExceptionLabel;
	
	// TODO: Add parameters for reactions that take place across multiple
	// regions (i.e., surface interactions)
	
	// Type of surface reaction.
	// Affects how the reaction probability is calculated
	// Default is RXN_NORMAL, which determines reaction probability from
	// reaction rate as if it were a solution reaction
	short surfRxnType;
};

//
// Function Declarations
//

/* Allocate space for the regionArray members needed for chemical
 * reactions and initialize their values based on the chem_rxn
 * array of chem_rxn_struct
*/
void initializeRegionChemRxn(const short NUM_REGIONS,
	struct region regionArray[],
	const unsigned short NUM_MOL_TYPES,
	const unsigned short MAX_RXNS,
	const struct chem_rxn_struct chem_rxn[],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES]);

void deleteRegionChemRxn(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	struct region regionArray[]);

#endif // CHEM_RXN_H
