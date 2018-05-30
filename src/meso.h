/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * meso.h - heap of all mesoscopic subvolumes in simulation environment
 *
 * Last revised for AcCoRD v1.2 (2018-05-30)
 *
 * Revision history:
 *
 * Revision v1.2 (2018-05-30)
 * - corrected calculation of diffusion propensity when more than 1 molecule is added to a subvolume
 *
 * Revision v1.1 (2016-12-24)
 * - added uniform flow to the diffusion algorithm
 *
 * Revision v1.0 (2016-10-31)
 * - moved mesoscopic structure fields from subvolume struct to meso subvolume struct
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - modified random number generation. Now use PCG via a separate interface file.
 *
 * Revision v0.4.1
 * - improved use and format of error messages
 *
 * Revision v0.4
 * - modified use of unsigned long and unsigned long long to uint32_t and uint64_t
 * - removed deprecated debug functions
 * - renamed region structure and array of structures to accommodate naming
 * of new boundary structure
 *
 * Revision v0.3.1
 * - header added
*/
#ifndef MESO_H
#define MESO_H

#include "region.h"
#include "subvolume.h"
#include "rand_accord.h" // For PRNGs

//
// Constant definitions
//


//
// Subvolume data type declarations
//

/* The mesoSubvolume3D structure lists properties for a single mesoscopic subvolume,
* The entire simulation environment has only one array of these structures.
*/
struct mesoSubvolume3D {
	uint32_t subID; // ID of mesoscopic subvolume in list of all subvolumes
	
	// Total number of reactions (including diffusion) that molecules in this
	// subvolume can participate in
	unsigned short num_rxn;
	
	// Mesoscopic properties
	double totalProp; // Total reaction propensity
	double t_rxn; // Time of next reaction in this subvolume
	uint32_t heapID; // Index of subvolume in next reaction heap
	
	double * rxnProp;
	
	// Index of the first chemical reaction in rxnProp array. Used
	// to identify whether the current reaction is diffusion or chemical
	// and to act as a reference for finding elements corresponding
	// to chemical reactions
	unsigned short firstChemRxn;
	
	// Diffusion transition rate from a boundary subvolume to all
	// of its neighbours
	// (whether or not each neighbour is in a different region)
	// Only allocated if subvolume bBoundary == true or there is a flow
	// Size is NUM_MOL_TYPES x num_neigh.
	// Each element gives the diffusion rate to the corresponding
	// neighbour in the neighID array (in subvolume structure).
	double ** diffRateNeigh;
		
	// Number of each type of molecule in subvolume
	// Length NUM_MOL_TYPES
	uint64_t * num_mol;
	
	// FUTURE MEMBERS
};

//
// Function Declarations
//

// Allocate Array of pointers to Mesoscopic subvolume structures
void allocateMesoSubArray(const uint32_t numMesoSub,
	struct mesoSubvolume3D ** mesoSubArray);
	
// Allocate Heap (sorted mesoscopic subvolume IDs) and associated arrays
void allocateMesoHeapArray(const uint32_t numMesoSub,
	uint32_t ** heap_subvolID,
	uint32_t (**heap_childID)[2],
	bool (**b_heap_childValid)[2]);

// Free memory allocated to array of Mesoscopic subvolume structures
void deleteMesoSubArray(const uint32_t numMesoSub,
	struct mesoSubvolume3D mesoSubArray[],
	const struct subvolume3D subvolArray[],
	const unsigned short NUM_MOL_TYPES,
	const short NUM_REGIONS);

// Initialize Array of pointers to Mesoscopic subvolume structures
void initializeMesoSubArray(const uint32_t numMesoSub,
	const uint32_t numSub,
	struct mesoSubvolume3D mesoSubArray[],
	const struct subvolume3D subvolArray[],
	const double SUBVOL_BASE_SIZE,
	const unsigned short NUM_MOL_TYPES,
	const unsigned short MAX_RXNS,
	struct region regionArray[],
	const short NUM_REGIONS,
	uint32_t subCoorInd[numSub][3],
	unsigned short ** subNeighDir,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES]);
	
// Reset propensities and reaction times for all subvolumes
void resetMesoSubArray(const uint32_t numMesoSub,
	struct mesoSubvolume3D mesoSubArray[],
	const struct subvolume3D subvolArray[],
	const unsigned short NUM_MOL_TYPES,
	const unsigned short MAX_RXNS,
	const short NUM_REGIONS,
	struct region regionArray[]);

// Update propensities and next reaction time of subvolume
void updateMesoSub(const uint32_t curSub,
	bool bChemRxn,
	uint64_t numMolChange[],
	bool bMolAdd[],
	const unsigned short singleMolID,
	bool bAllowCorrectedTime,
	const uint32_t numMesoSub,
	struct mesoSubvolume3D mesoSubArray[],
	struct subvolume3D subvolArray[],
	double tCur,
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	struct region regionArray[]);

// Update propensities, next reaction times, and heap location of subvolumes
// along micro/meso interface.
// This function is meant to be called AFTER validateMolecules() in
// micro_molecule.c, which populates the bNeedUpdate and numMolFromMicro
// members of the array of region2D structures
void updateMesoSubBoundary(const uint32_t numSub,
	const uint32_t numMesoSub,
	struct mesoSubvolume3D mesoSubArray[],
	struct subvolume3D subvolArray[],
	struct region regionArray[],
	bool bTrue[],
	double tCur,
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	uint32_t heap_subvolID[],
	uint32_t heap_childID[][2],
	bool heap_childValid[][2]);

// Sum terms in reaction propensity vector
double updateTotalProp(const double rxnProp[],
	const unsigned short numChemRxn);

// Determine next relative subvolume reaction time
double mesoSubCalcTime(struct mesoSubvolume3D mesoSubArray[],
	const uint32_t ID);

// Heap Operations
void heapMesoFindChildren(const uint32_t numSub,
	uint32_t heap_childID[][2],
	bool heap_childValid[][2]);
	
void heapMesoDelete(const uint32_t numElements,
	uint32_t heap_subvolID[],
	uint32_t heap_childID[][2],
	bool heap_childValid[][2]);
	
void heapMesoBuild(const uint32_t numSub,
	struct mesoSubvolume3D mesoSubArray[],
	uint32_t heap_subvolID[],
	const unsigned int num_heap_levels,
	uint32_t heap_childID[][2],
	bool heap_childValid[][2]);
	
uint32_t heapMesoUpdate(const uint32_t numSub,
	struct mesoSubvolume3D mesoSubArray[],
	uint32_t heap_subvolID[],
	const uint32_t heapID,
	uint32_t heap_childID[][2],
	bool heap_childValid[][2]);
	
uint32_t heapMesoCompareDown(const uint32_t numSub,
	struct mesoSubvolume3D mesoSubArray[],
	uint32_t heap_subvolID[],
	const uint32_t parent,
	uint32_t heap_childID[][2],
	bool heap_childValid[][2]);
	
uint32_t heapMesoCompareUp(const uint32_t numSub,
	struct mesoSubvolume3D mesoSubArray[],
	uint32_t heap_subvolID[],
	const uint32_t child);
	
void heapMesoSwap(struct mesoSubvolume3D mesoSubArray[],
	uint32_t heap_subvolID[],
	uint32_t index1, uint32_t index2);
	
#endif // MESO_H
