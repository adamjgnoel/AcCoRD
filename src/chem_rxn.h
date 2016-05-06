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
 *
 * Last revised for AcCoRD LATEST_VERSION
 *
 * Revision history:
 *
 * Revision LATEST_VERSION
 * - added label, bReversible, and labelCoupled so that reactions can be named
 * and coupled together
 * - added bReleaseProduct for surface reactions to indicate which products are
 * released from the surface
 * - updated reaction probabilities for surface reactions so that user has
 * choices for what calculation to use. Added types to store user choices
 * - adsorption, desorption, and membrane probability calculations are mostly based on
 * S.S. Andrews, "Accurate particle-based simulation of adsorption, desorption
 * and partial transmission" Physical Biology, vol. 6, p.046015, 2009
 * - constrained absorbing and desorbing reactions to one per type of molecule
 * at a given region. In many cases these reactions are now treated separately
 * from other types of 1st order reactions
 * - constrained membrane reactions to one inner and one outer reaction per type of
 * molecule at a given membrane region.
 *
 * Revision v0.5 (2016-04-15)
 * - removed limit on number of molecule types
 * - removed limit on number of products in a reaction
 * - added region label when giving errors about region initialization
 * - added bSurface, surfRxnType members to reaction structure to add specialized reactions
 * - implementations of some surface reactions
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
#include <complex.h> // for complex error function
#include <limits.h> // For USHRT_MAX

#include "region.h"
#include "global_param.h"
#include "cerf.h" // for complex error function
#include "defs.h" // definitions for cerf.h

/*
* Data Type Declarations
*/
//struct region; // Forward declaration (declared in region.h)

/* The chem_rxn structure defines the parameters relevant to the
 * firing of a single chemical reaction
*/
struct chem_rxn_struct { // Used to define a single chemical reaction
	// label is an optional string to name the reaction.
	// Needed if reaction is reversible so that we can identify the coupled reaction
	char * label;
	
	// Is reaction reversible?
	// Reversibility impacts the calculation of the reaction rates, especially
	// for surface transition reactions
	bool bReversible;
	
	// If reversible, what is the name of the coupled reaction?
	char * labelCoupled;
	
	// Indicate the indices of the reactants and the number of each
	// Length is the number of molecule types
	uint32_t * reactants;
	
	// Indicate the indices of the products and the number of each
	// Length is the number of molecule types
	uint32_t * products;
	
	// Base reaction rate k (units depends on order of reaction)
	double k;
	
	// Is the reaction a surface reaction?
	// This will affect where a reaction will take place by default,
	// as indicated by bEverywhere
	bool bSurface;
	
	// Are products of a surface reaction released from surface?
	// If true for given product, then the product molecule is placed in closest neighboring
	// region when it is created.
	// Length is the number of molecule types
	bool * bReleaseProduct;
	
	// Type of product release from surface
	// Default is PROD_PLACEMENT_LEAVE, which will leave the molecule next to
	// the surface.
	// Defined only for desorbing reactions
	short releaseType;
	
	// Can the reaction take place anywhere by default?
	// Actual regions will depend on value of bSurface and whether a given
	// region is a normal region or a surface
	bool bEverywhere;
	
	// Number of regions that are exceptions to the location default
	short numRegionExceptions;
	
	// Names of regions that are exceptions to the location default
	// Length is numRegionExceptions
	char ** regionExceptionLabel;
	
	// Type of surface reaction.
	// Affects how the reaction probability is calculated
	// Constraints are imposed on the number of surface transition reactions
	// Default is RXN_NORMAL, which determines reaction probability from
	// reaction rate as if it were a solution reaction
	short surfRxnType;
	
	// Type of surface reaction probability calculation
	// Default is RXN_PROB_NORMAL, which is inaccurate for surface transition reactions
	short rxnProbType;
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

// Calculate probability of desorption reaction for specified time step
bool calculateDesorptionProb(double * rxnProb,
	const short curRegion,
	const unsigned short curMolType,
	const unsigned short curRegionRxn,
	const double dt,
	const short NUM_REGIONS,
	const struct region regionArray[],
	const unsigned short NUM_MOL_TYPES,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES]);

// Calculate probability of absorption reaction for specified time step
double calculateAbsorptionProb(const short curRegion,
	const unsigned short curMolType,
	const unsigned short curRegionRxn,
	const double dt,
	const short NUM_REGIONS,
	const struct region regionArray[],
	const unsigned short NUM_MOL_TYPES,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES]);

// Calculate probability of membrane transition for specified time step
double calculateMembraneProb(const short curRegion,
	const unsigned short curMolType,
	const unsigned short curRegionRxn,
	const double dt,
	const short NUM_REGIONS,
	const struct region regionArray[],
	const unsigned short NUM_MOL_TYPES,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES]);

#endif // CHEM_RXN_H
