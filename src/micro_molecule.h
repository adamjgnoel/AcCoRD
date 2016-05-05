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
 * Last revised for AcCoRD LATEST_VERSION
 *
 * Revision history:
 *
 * Revision LATEST_VERSION
 * - updated first order reaction functions to account for surface reactions that
 * release products from the surface. This is done in a common function (for both old
 * and recent molecules). Placement of products depends on user configuration
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

#include "randistrs.h" // For PRNGs
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
	bool bNeedUpdate; // Indicate whether molecule needs to be moved in current
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
	struct subvolume3D subvolArray[],
	double sigma[NUM_REGIONS][NUM_MOL_TYPES],
	const struct chem_rxn_struct chem_rxn[],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES]);

void diffuseOneMolecule(ItemMol3D * molecule, double sigma);

void diffuseOneMoleculeRecent(ItemMolRecent3D * molecule, double DIFF_COEF);

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
	const bool bRecent);

// If a molecule is the product of a surface reaction and it is supposed
// to be released from the surface, find the destination region
unsigned short findDestRegion(const double point[3],
	const unsigned short curRegion,
	const struct region regionArray[]);

void transferMolecules(ListMolRecent3D * molListRecent, ListMol3D * molList);

bool validateMolecule(double newPoint[3],
	double oldPoint[3],
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const short curRegion,
	short * newRegion,
	short * transRegion,
	const struct region regionArray[],
	short molType,
	bool * bReaction,
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
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const struct region regionArray[],
	short molType,
	bool * bReaction,
	unsigned short * curRxn,
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
