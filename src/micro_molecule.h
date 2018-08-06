/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * micro_molecule.h - 	linked list of individual molecules in same
 * 						microscopic region
 *
 * Last revised for AcCoRD v1.4 (2018-08-06)
 *
 * Revision history:
 *
 * Revision v1.4 (2018-08-06)
 * - added a priori monte carlo (APMC) absorption algorithm as a new surface
 * reaction type. Includes settings for how to define the a priori absorption
 * probability calculation and whether/how to apply a threshold to turn it off
 * - corrected missing assignment for first order reactions that don't have any
 * product molecules to place. Led to memory error when first molecule in a list
 * has to be removed
 *
 * Revision v1.2 (2018-05-30)
 * - fixed implementation of replication reactions, where a first order reactant produces
 * at least one copy of itself. If such a reactant reacted again within the same
 * microscopic time step, then the new molecule(s) previously went missing.
 *
 * Revision v1.1 (2016-12-24)
 * - simplified detection of whether molecules flow or diffuse in each region
 * - added uniform flow to the diffusion algorithm
 * - modified meso-to-micro hybrid transition algorithm when a molecule is placed
 * in the microscopic regime. Now, the trajectory of the molecule will be tracked
 * to make sure that it can reach its intended destination. Reflections are added
 * as necessary. Molecule is assumed to start from the middle of the mesoscopic
 * subvolume.
 *
 * Revision v1.0 (2016-10-31)
 * - added specifying diffusion coefficient that applies to specific surface
 * interaction reactions.
 * - moved mesoscopic structure fields from subvolume struct to meso subvolume struct
 *
 * Revision v0.7 (public beta, 2016-07-09)
 * - set microscopic partial time step to 0 when creating new molecule from meso
 * diffusion.
 * - corrected what region is indexed when new microscopic molecule diffuses away
 * from hybrid interface
 * - use molecule's partial time step when determining whether it entered and exited
 * mesoscopic regime within a time step
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - added check for molecules entering mesoscopic regime "during" a microscopic time step,
 * i.e., when molecule is in micro at start and end of diffusion step.
 * - added function for placing molecules in microscopic regime when they come from the
 * mesoscopic region
 * - modified random number generation. Now use PCG via a separate interface file.
 * - added bimolecular chemical reactions in microscopic regime. Transition of reactants
 * to reaction site is tested if reaction site is not in same region. Transition of
 * product molecules from reaction site is also tested if unbinding radius is used.
 * Molecules are compared with all potential reactants in current region and neighboring
 * regions where the reaction is also valid. A molecule cannot participate in more than
 * one bimolecular reaction is a single time step (as a reactant or product). Reaction
 * site and displacement of products depends on relative diffusion coefficients
 * - added check on total number of chemical reactions in a region before checking for
 * surface reactions during diffusion validation, since many surface reaction structure
 * members are not initialized unless there is at least one reaction at the surface
 * - corrected bug where a molecule is not correctly reflected off of a child region
 * surface
 *
 * Revision v0.5.1 (2016-05-06)
 * - updated first order reaction functions to account for surface reactions that
 * release products from the surface. This is done in a common function (for both old
 * and recent molecules). Placement of products depends on user configuration
 * - fixed membrane transition reactions. They can be reversible and should not have
 * product molecules explicitly defined.
 * - added calls to new functions to determine adsorption/desorption probabilities
 * for recent molecules
 * - corrected how molecules are locked to region boundary when they cross regions
 * - updating reaction probabilities for surface reactions so that user has
 * choices for what calculation to use.
 *
 * Revision v0.5 (2016-04-15)
 * - added surface reactions, including membrane transitions
 * - added switch to record all molecules in a region instead of just those
 * within some specified boundary
 * - corrected distance to end point when a molecule is "pushed" into a neighboring
 * region
 * - added fail check to while loop when a molecule is "pushed" into a 
 * neighboring region. Error will display if we did not end up in specified
 * region or one of its children.
 * - corrected molecule diffusion validation algorithm to reflect off of the correct
 * surface boundary when a reflection is needed
 *
 * Revision v0.4.1
 * - improved use and format of error messages
 *
 * Revision v0.4
 * - modified use of unsigned long and unsigned long long to uint32_t and uint64_t
 * - re-wrote diffusion validation so that molecule path is followed from its initial
 * position until it reaches the final position, hits a reflective boundary, or
 * is absorbed into a mesoscopic subvolume (whichever occurs first).
 * - renamed region structure and array of structures to accommodate naming
 * of new boundary structure
 *
 * Revision v0.3.1
 * - header added
 *
 * Created 2014-11-23
*/
#ifndef MICRO_MOLECULE_H
#define MICRO_MOLECULE_H

#include <stdio.h> // DEBUG: to create and edit files
#include <stdlib.h> // for exit(), malloc, free, NULL
#include <stdbool.h> // for C++ bool naming, requires C99
#include <limits.h> // For SHRT_MAX
#include <math.h> // For sqrt()
#include <complex.h> // for complex error function

#include "rand_accord.h" // For PRNGs
#include "region.h"
#include "meso.h"
#include "subvolume.h"
#include "chem_rxn.h"
#include "global_param.h" // for common global parameters
#include "cerf.h" // for complex error function
#include "defs.h" // definitions for cerf.h

// micro_molecule specific declarations

struct molecule_list3D {
	double x, y, z; // Coordinates of centre of molecule
	bool bNeedUpdate; 	// Indicate whether molecule can still move or react in current
						// time step
};

struct molecule_recent_list3D {
	double x, y, z; // Coordinates of centre of molecule
	double dt_partial; // Time between molecule creation and next micro time step
};

// General type declarations

typedef struct molecule_list3D ItemMol3D;
typedef struct molecule_recent_list3D ItemMolRecent3D;

typedef struct node3D{
	ItemMol3D item;
	struct node3D * next;
} NodeMol3D;

typedef struct nodeRecent3D{
	ItemMolRecent3D item;
	struct nodeRecent3D * next;
} NodeMolRecent3D;

typedef NodeMol3D * ListMol3D;
typedef NodeMolRecent3D * ListMolRecent3D;

// micro_molecule specific Prototypes

bool addMolecule(ListMol3D * p_list, double x, double y, double z);

bool addMoleculeRecent(ListMolRecent3D * p_list, double x, double y, double z, double dt_partial);

void moveMolecule(ItemMol3D * molecule, double x, double y, double z);

void moveMoleculeRecent(ItemMolRecent3D * molecule, double x, double y, double z);

void diffuseMolecules(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	ListMol3D p_list[NUM_REGIONS][NUM_MOL_TYPES],
	ListMolRecent3D p_listRecent[NUM_REGIONS][NUM_MOL_TYPES],
	const struct region regionArray[],
	struct mesoSubvolume3D mesoSubArray[],
	double sigma[NUM_REGIONS][NUM_MOL_TYPES],
	const struct chem_rxn_struct chem_rxn[],
	const double HYBRID_DIST_MAX,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES]);

void diffuseOneMolecule(ItemMol3D * molecule, double sigma);

// Move one molecule according to a flow vector
void flowTransportOneMolecule(ItemMol3D * molecule,
	const unsigned short flowType,
	double *flowConstant);

void diffuseOneMoleculeRecent(ItemMolRecent3D * molecule, double DIFF_COEF);

// Move one molecule according to a flow vector
void flowTransportOneMoleculeRecent(ItemMolRecent3D * molecule,
	const unsigned short flowType,
	double *flowVector);

// Did molecule enter mesoscopic region while diffusing?
bool bEnterMesoIndirect(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const struct region regionArray[],
	const short curType,
	const short curRegion,
	short * mesoRegion,
	const double oldPoint[3],
	const double newPoint[3],
	uint32_t * newSub,
	const double HYBRID_DIST_MAX,
	const double tLeft,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES]);

// Place a molecule entering microscopic region from a mesoscopic subvolume
bool placeInMicroFromMeso(const unsigned short curRegion,
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const unsigned short destRegion,
	uint32_t * newSub,
	const struct region regionArray[],
	const uint32_t curBoundSub,
	const bool bSmallSub,
	const unsigned short curMolType,
	ListMolRecent3D pRecentList[NUM_REGIONS][NUM_MOL_TYPES],
	const struct chem_rxn_struct chem_rxn[],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES]);

void rxnFirstOrder(const unsigned short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	unsigned short curRegion,
	ListMol3D p_list[NUM_REGIONS][NUM_MOL_TYPES],
	const struct region regionArray[],
	unsigned short curMolType,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	ListMolRecent3D pRecentList[NUM_REGIONS][NUM_MOL_TYPES]);

void rxnFirstOrderRecent(const unsigned short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	unsigned short curRegion,
	ListMolRecent3D pRecentList[NUM_REGIONS][NUM_MOL_TYPES],
	ListMol3D p_list[NUM_REGIONS][NUM_MOL_TYPES],
	const struct region regionArray[],
	unsigned short curMolType,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	bool bCheckCount,
	uint32_t numMolCheck[NUM_REGIONS][NUM_MOL_TYPES]);

// Place products of 1st order reaction
void rxnFirstOrderProductPlacement(const NodeMol3D * curMol,
	const NodeMolRecent3D * curMolRecent,
	const unsigned short curRxn,
	const unsigned short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	unsigned short curRegion,
	ListMol3D p_list[NUM_REGIONS][NUM_MOL_TYPES],
	ListMolRecent3D pRecentList[NUM_REGIONS][NUM_MOL_TYPES],
	const struct region regionArray[],
	unsigned short curMolType,	
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	const bool bRecent,
	bool * bProductIsReactant);

// If a molecule is the product of a surface reaction and it is supposed
// to be released from the surface, find the destination region
unsigned short findDestRegion(const double point[3],
	const unsigned short curRegion,
	const struct region regionArray[]);

// Check all second order reactions for microscopic regions
void rxnSecondOrder(const unsigned short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	ListMol3D p_list[NUM_REGIONS][NUM_MOL_TYPES],
	const struct region regionArray[],
	struct mesoSubvolume3D mesoSubArray[],
	const struct chem_rxn_struct chem_rxn[],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES]);

// Compare distance between 2 molecules with threshold
// Return true if distance is greater than threshold
bool moleculeSeparation(ItemMol3D * molecule1, ItemMol3D * molecule2, double threshSq);

void transferMolecules(ListMolRecent3D * molListRecent, ListMol3D * molList);

bool validateMolecule(double newPoint[3],
	double oldPoint[3],
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const short curRegion,
	short * newRegion,
	short * transRegion,
	bool * bPointChange,
	const struct region regionArray[],
	unsigned short molType,
	bool * bReaction,
	bool * bApmcRevert,
	bool bRecent,
	double dt,
	const struct chem_rxn_struct chem_rxn[],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	unsigned short * curRxn);

// Recursively follow a molecule's path through region boundaries from its diffusion
// start and end points
// Return whether molecule path had to be changed
bool followMolecule(const double startPoint[3],
	double endPoint[3],
	double lineVector[3],
	double lineLength,
	const short startRegion,
	short * endRegion,
	short * transRegion,
	bool * bPointChange,
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const struct region regionArray[],
	unsigned short molType,
	bool * bReaction,
	unsigned short * curRxn,
	bool * bApmcRevert,
	bool bRecent,
	double dt,
	const struct chem_rxn_struct chem_rxn[],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	unsigned int depth);

uint64_t countMolecules(ListMol3D * p_list,
	int obsType,
	double boundary[]);

uint64_t countMoleculesRecent(ListMolRecent3D * p_list,
	int obsType,
	double boundary[]);

uint64_t recordMolecules(ListMol3D * p_list,
	ListMol3D * recordList,
	int obsType,
	double boundary[],
	bool bRecordPos,
	bool bRecordAll);

uint64_t recordMoleculesRecent(ListMolRecent3D * p_list,
	ListMol3D * recordList,
	int obsType,
	double boundary[],
	bool bRecordPos,
	bool bRecordAll);

bool isMoleculeObserved(ItemMol3D * molecule,
	int obsType,
	double boundary[]);

bool isMoleculeObservedRecent(ItemMolRecent3D * molecule,
	int obsType,
	double boundary[]);

// General Prototypes

void initializeListMol(ListMol3D * p_list);

void initializeListMolRecent(ListMolRecent3D * p_list);

bool isListMol3DEmpty(const ListMol3D * p_list);

bool isListMol3DRecentEmpty(const ListMolRecent3D * p_list);

void emptyListMol(ListMol3D * p_list);

void emptyListMol3DRecent(ListMolRecent3D * p_list);


#endif // MICRO_MOLECULE_H
