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
 * Last revised for AcCoRD v0.5
 *
 * Revision history:
 *
 * Revision v0.5
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
#include "randistrs.h" // For PRNGs
#include "region.h"
#include "meso.h"
#include "subvolume.h"
#include "global_param.h" // for common global parameters

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

bool addMolecule3D(ListMol3D * p_list, double x, double y, double z);

bool addMoleculeRecent3D(ListMolRecent3D * p_list, double x, double y, double z, double dt_partial);

void moveMolecule3D(ItemMol3D * molecule, double x, double y, double z);

void moveMoleculeRecent3D(ItemMolRecent3D * molecule, double x, double y, double z);

void diffuseMolecules3D(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	ListMol3D p_list[NUM_REGIONS][NUM_MOL_TYPES],
	ListMolRecent3D p_listRecent[NUM_REGIONS][NUM_MOL_TYPES],
	const struct region regionArray[],
	struct mesoSubvolume3D mesoSubArray[],
	struct subvolume3D subvolArray[],
	double sigma[NUM_REGIONS][NUM_MOL_TYPES],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES]);

void diffuseOneMolecule3D(ItemMol3D * molecule, double sigma);

void diffuseOneMoleculeRecent3D(ItemMolRecent3D * molecule, double DIFF_COEF);

void rxnFirstOrder3D(ListMol3D * p_list,
	const struct region regionArray,
	unsigned short curMolType,
	const unsigned short NUM_MOL_TYPES,
	ListMolRecent3D pRecentList[NUM_MOL_TYPES]);

void rxnFirstOrderRecent3D(const unsigned short NUM_MOL_TYPES,
	ListMolRecent3D pRecentList[NUM_MOL_TYPES],
	const struct region regionArray,
	unsigned short curMolType,
	bool bCheckCount,
	uint32_t numMolCheck[NUM_MOL_TYPES]);

void transferMolecules3D(ListMolRecent3D * molListRecent, ListMol3D * molList);

bool validateMolecule(double newPoint[3],
	double oldPoint[3],
	const short NUM_REGIONS,
	const short curRegion,
	short * newRegion,
	short * transRegion,
	const struct region regionArray[]);

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
	const struct region regionArray[],
	unsigned int depth);

uint64_t countMolecules3D(ListMol3D * p_list,
	int obsType,
	double boundary[]);

uint64_t countMoleculesRecent3D(ListMolRecent3D * p_list,
	int obsType,
	double boundary[]);

uint64_t recordMolecules3D(ListMol3D * p_list,
	ListMol3D * recordList,
	int obsType,
	double boundary[],
	bool bRecordPos);

uint64_t recordMoleculesRecent3D(ListMolRecent3D * p_list,
	ListMol3D * recordList,
	int obsType,
	double boundary[],
	bool bRecordPos);

bool isMoleculeObserved3D(ItemMol3D * molecule,
	int obsType,
	double boundary[]);

bool isMoleculeObservedRecent3D(ItemMolRecent3D * molecule,
	int obsType,
	double boundary[]);

// General Prototypes

void initializeListMol3D(ListMol3D * p_list);

void initializeListMolRecent3D(ListMolRecent3D * p_list);

bool isListMol3DEmpty(const ListMol3D * p_list);

bool isListMol3DRecentEmpty(const ListMolRecent3D * p_list);

void emptyListMol3D(ListMol3D * p_list);

void emptyListMol3DRecent(ListMolRecent3D * p_list);


#endif // MICRO_MOLECULE_H
