/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * meso.c - heap of all mesoscopic subvolumes in simulation environment
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
#include <stdio.h>
#include <math.h> // for fabs(), INFINITY
#include <stdlib.h> // for exit(), malloc
#include <stdbool.h> // for C++ bool conventions
#include <inttypes.h> // for extended integer type macros
#include "meso.h"
#include "subvolume.h"
#include "region.h"

//
// "Private" Declarations
//

//
// Definitions
//

/* Allocate space for an array of mesoscopic subvolume
* structures.
*/
void allocateMesoSubArray(const uint32_t numMesoSub,
	struct mesoSubvolume3D ** mesoSubArray)
{
	
	*mesoSubArray = malloc(numMesoSub*sizeof(struct mesoSubvolume3D));
		
	if(*mesoSubArray == NULL){
		fprintf(stderr, "ERROR: Memory allocation for array of mesoscopic subvolume structures.\n");
		exit(EXIT_FAILURE);
	}
}

/* Allocate space for the heap of IDs pointing to the array of mesoscopic
 * subvolume structures, a 2D array listing the child IDs of each heap element,
 * and a 2D array of bools indicating whether each child is valid
*/
void allocateMesoHeapArray(const uint32_t numMesoSub,
	uint32_t ** heap_subvolID,
	uint32_t (**heap_childID)[2],
	bool (**b_heap_childValid)[2])
{
	
	*heap_subvolID = malloc(numMesoSub*sizeof(uint32_t));
	*heap_childID = malloc(numMesoSub*sizeof(uint32_t [2]));
	*b_heap_childValid = malloc(numMesoSub*sizeof(bool [2]));
		
	if(*heap_subvolID == NULL || *heap_childID == NULL || *b_heap_childValid == NULL)
	{
		fprintf(stderr, "ERROR: Memory allocation for heap of mesoscopic subvolume structures.\n");
		exit(EXIT_FAILURE);
	}
}

void deleteMesoSubArray(const uint32_t numMesoSub,
	struct mesoSubvolume3D mesoSubArray[],
	const struct subvolume3D subvolArray[],
	const unsigned short NUM_MOL_TYPES,
	const short NUM_REGIONS)
{
	uint32_t curMesoSub;
	uint32_t curSub = 0;
	uint32_t curMolType = 0;

	if(mesoSubArray == NULL)
		return;
	
	for(curMesoSub = 0; curMesoSub < numMesoSub; curMesoSub++)
	{
		if(mesoSubArray[curMesoSub].rxnProp != NULL)
			free(mesoSubArray[curMesoSub].rxnProp);
		if(mesoSubArray[curMesoSub].num_mol != NULL)
			free(mesoSubArray[curMesoSub].num_mol);
		curSub = mesoSubArray[curMesoSub].subID;
		if(NUM_REGIONS > 1 && mesoSubArray[curMesoSub].diffRateNeigh != NULL)
		{
			for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
			{
				if(mesoSubArray[curMesoSub].diffRateNeigh[curMolType] != NULL)
					free(mesoSubArray[curMesoSub].diffRateNeigh[curMolType]);
			}
			free(mesoSubArray[curMesoSub].diffRateNeigh);
		}
			
	}
		
	free(mesoSubArray);
}

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
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES])
{
	uint32_t curSub;
	uint32_t curMesoSub = 0;
	
	// Parameters for diffusion propensity calculation and storage
	unsigned short curMolType;
	short curRegion;
	double h_i, h_j; // Subvolume sizes (used for finding transition rates)
	double curSubBound[6]; // Boundary of current subvolume
	double curNeighBound[6]; // Boundary of prospective neighbor subvolume
	double boundOverlap[6]; // Overlap area of adjacent subvolumes
	bool bFlow; // Can molecules flow to neighboring subvolume?
	double flowVal; // Value of flow vector in direction of neighbor subvolume
	double flowProp; // Flow propensity
	double relOverlap; // Relative overlap of neighboring subvolume
	uint32_t curNeighID = 0; // Current Subvolume neighbour ID
	uint32_t neighID = 0; // Subvolume neighbour ID in master subvolume list
	short int neighRegion;
	double boundAdjError = SUBVOL_BASE_SIZE * SUB_ADJ_RESOLUTION;
	
	for(curSub=0; curSub < numSub; curSub++)
	{
		if(regionArray[subvolArray[curSub].regionID].spec.bMicro)
			continue; // Microscopic subvolumes are ignored
		
		if(curMesoSub < numMesoSub)
		{ // Failsafe so that we do not write to non-allocated memory
			mesoSubArray[curMesoSub].subID = curSub;
			mesoSubArray[curMesoSub].totalProp = 0.;
			mesoSubArray[curMesoSub].t_rxn = INFINITY;
			mesoSubArray[curMesoSub].heapID = 0;
			mesoSubArray[curMesoSub].firstChemRxn =
				subvolArray[curSub].num_neigh * NUM_MOL_TYPES;
			mesoSubArray[curMesoSub].rxnProp = 
				malloc((mesoSubArray[curMesoSub].firstChemRxn + MAX_RXNS)
				*sizeof(double));
			mesoSubArray[curMesoSub].num_mol = 
				malloc(NUM_MOL_TYPES*sizeof(uint64_t));
			if(mesoSubArray[curMesoSub].rxnProp == NULL ||
				mesoSubArray[curMesoSub].num_mol == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation to mesoscopic structure parameters for mesoscopic subvolume %" PRIu32 " (subvolume ID %" PRIu32 ").\n", curMesoSub, curSub);
				exit(EXIT_FAILURE);
			}
			curMesoSub++;
		} else
		{
			fprintf(stderr, "ERROR: Error: Tried to write to mesoscopic array beyond allocated memory. Intended mesoscopic subvolume %" PRIu32 " (subvolume ID %" PRIu32 ").\n", curMesoSub, curSub);
			exit(EXIT_FAILURE);
		}
	}
	
	if(NUM_REGIONS > 1 ||
		(numMesoSub > 0 && regionArray[subvolArray[0].regionID].bFlow))
	{
		// Determine transition rates out of mesoscopic subvolumes that are along
		// the boundary of their respective region, or that can have flow
		for(curMesoSub=0; curMesoSub < numMesoSub; curMesoSub++)
		{
			curSub = mesoSubArray[curMesoSub].subID;
			curRegion = subvolArray[curSub].regionID;
			if(!subvolArray[curSub].bBoundary &&
				!regionArray[subvolArray[curRegion].regionID].bFlow)
			{
				mesoSubArray[curMesoSub].diffRateNeigh = NULL;
				continue;  // This subvolume is not meso, not along region boundary,
						   // and has no flow
			}
			
			// Allocate memory to store transition rate for this subvolume
			mesoSubArray[curMesoSub].diffRateNeigh = 
				malloc(NUM_MOL_TYPES*sizeof(mesoSubArray[curMesoSub].diffRateNeigh));
			
			if(mesoSubArray[curMesoSub].diffRateNeigh == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation for mesoscopic diffusion rates for meso subvolume %" PRIu32 ", which is along the boundary or region %u.\n", curMesoSub, subvolArray[curSub].regionID);
				exit(EXIT_FAILURE);
			}
			
			for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
			{
				mesoSubArray[curMesoSub].diffRateNeigh[curMolType] =
					malloc(subvolArray[curSub].num_neigh*sizeof(double));
				if(mesoSubArray[curMesoSub].diffRateNeigh[curMolType] == NULL){
					fprintf(stderr, "ERROR: Memory allocation for mesoscopic diffusion rates for meso subvolume %" PRIu32 ", which is along the boundary or region %u.\n", curMesoSub, subvolArray[curSub].regionID);
					exit(EXIT_FAILURE);
				}
			}
			
			h_i = regionArray[curRegion].actualSubSize;
			
			findSubvolCoor(curSubBound, regionArray[curRegion], subCoorInd[curSub]);
			
			for(curNeighID = 0; curNeighID < subvolArray[curSub].num_neigh;
				curNeighID++)
			{
				// Find actual neighbor index and region
				neighID = subvolArray[curSub].neighID[curNeighID];
				neighRegion = subvolArray[neighID].regionID;
				if (curRegion == neighRegion)
				{ // Transition rate only depends on source volume and flow rate
					for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
					{
						mesoSubArray[curMesoSub].diffRateNeigh[curMolType][curNeighID] =
							DIFF_COEF[curRegion][curMolType]/h_i/h_i;
					}						
				} else
				{
					if(regionArray[curRegion].spec.type !=
						regionArray[neighRegion].spec.type)
					{
						// Regions are not of the same type (at least one is a surface)
						// Normal diffusion between these regions is not possible
						for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
						{
							// TODO: Correct this preliminary solution which prevents
							// any transition out of a meso subvolume to a surface.
							// Correct implementation would base transition rate
							// on corresponding chemical reaction probabilities
							mesoSubArray[curMesoSub].diffRateNeigh[curMolType][curNeighID] = 0.;
						}
						continue;
					}
					
					// TODO: Need to catch cases where a membrane lies in between two
					// subvolumes so that the diffusion rate can be adjusted properly
					
					if(regionArray[neighRegion].spec.bMicro)
						h_j = h_i;
					else
						h_j = regionArray[neighRegion].actualSubSize;
					
					// Determine overlap area
					if(regionArray[neighRegion].spec.shape == RECTANGULAR_BOX)
					{
						findSubvolCoor(curNeighBound, regionArray[neighRegion],
							subCoorInd[neighID]);
						intersectBoundary(RECTANGULAR_BOX, curSubBound,
							RECTANGULAR_BOX, curNeighBound, boundOverlap);
					} else if (regionArray[neighRegion].spec.shape == SPHERE)
					{
						// Assume that overlap is entire subvolume face
						boundOverlap[0] = 0;
						boundOverlap[1] = h_i;
						boundOverlap[2] = 0;
						boundOverlap[3] = h_i;
						boundOverlap[4] = 0;
						boundOverlap[5] = 0;
					}
					
					for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
					{
						mesoSubArray[curMesoSub].diffRateNeigh[curMolType][curNeighID] =
							2*DIFF_COEF[curRegion][curMolType]/h_i/(h_i + h_j);
						if(fabs(boundOverlap[0] - boundOverlap[1]) > boundAdjError)
							mesoSubArray[curMesoSub].diffRateNeigh[curMolType][curNeighID] *=
								(boundOverlap[1] - boundOverlap[0])/h_i;
						if(fabs(boundOverlap[2] - boundOverlap[3]) > boundAdjError)
							mesoSubArray[curMesoSub].diffRateNeigh[curMolType][curNeighID] *=
								(boundOverlap[3] - boundOverlap[2])/h_i;
						if(fabs(boundOverlap[4] - boundOverlap[5]) > boundAdjError)
							mesoSubArray[curMesoSub].diffRateNeigh[curMolType][curNeighID] *=
								(boundOverlap[5] - boundOverlap[4])/h_i;
						
						if(regionArray[neighRegion].spec.bMicro &&
							DIFF_COEF[curRegion][curMolType] > 0.)
						{ // Include multiplier on propensity
							mesoSubArray[curMesoSub].diffRateNeigh[curMolType][curNeighID]
								*= 2*h_i/sqrt(DIFF_COEF[curRegion][curMolType]
								*PI*regionArray[neighRegion].spec.dt);
						}
					}
				}
				// Correct transition rates with flow parameters
				for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
				{
					if(regionArray[curRegion].spec.bFlow[curMolType])
					{ // Current molecule type can flow. Need to correct propensity if neighboring
					  // subvolume is in flow direction
						bFlow = false;
						flowVal = 0.;
						switch(regionArray[curRegion].spec.flowType[curMolType])
						{
							case FLOW_UNIFORM:
								switch(subNeighDir[curSub][curNeighID])
								{
									case LEFT:
										if(fabs(regionArray[curRegion].spec.flowVector[curMolType][0]) > 0.)
										{
											bFlow = true;
											flowVal = -regionArray[curRegion].spec.flowVector[curMolType][0];
										}
										break;
									case RIGHT:
										if(fabs(regionArray[curRegion].spec.flowVector[curMolType][0]) > 0.)
										{
											bFlow = true;
											flowVal = regionArray[curRegion].spec.flowVector[curMolType][0];
										}
										break;
									case DOWN:
										if(fabs(regionArray[curRegion].spec.flowVector[curMolType][1]) > 0.)
										{
											bFlow = true;
											flowVal = -regionArray[curRegion].spec.flowVector[curMolType][1];
										}
										break;
									case UP:
										if(fabs(regionArray[curRegion].spec.flowVector[curMolType][1]) > 0.)
										{
											bFlow = true;
											flowVal = regionArray[curRegion].spec.flowVector[curMolType][1];
										}
										break;
									case IN:
										if(fabs(regionArray[curRegion].spec.flowVector[curMolType][2]) > 0.)
										{
											bFlow = true;
											flowVal = -regionArray[curRegion].spec.flowVector[curMolType][2];
										}
										break;
									case OUT:
										if(fabs(regionArray[curRegion].spec.flowVector[curMolType][2]) > 0.)
										{
											bFlow = true;
											flowVal = regionArray[curRegion].spec.flowVector[curMolType][2];
										}
										break;
								}
								break;
						}
						if(bFlow)
						{
							// Calculate flow propensity
							if(curRegion == neighRegion)
								flowProp = flowVal/2./h_i;
							else
							{
								relOverlap = 1.;
								if(fabs(boundOverlap[0] - boundOverlap[1]) > boundAdjError)
									relOverlap *= (boundOverlap[1] - boundOverlap[0])/h_i;
								if(fabs(boundOverlap[2] - boundOverlap[3]) > boundAdjError)
									relOverlap *= (boundOverlap[3] - boundOverlap[2])/h_i;
								if(fabs(boundOverlap[4] - boundOverlap[5]) > boundAdjError)
									relOverlap *= (boundOverlap[5] - boundOverlap[4])/h_i;
								
								flowProp = flowVal/(h_i + h_j)*relOverlap;
							}
							
							if(fabs(flowProp) >
								mesoSubArray[curMesoSub].diffRateNeigh[curMolType][curNeighID])
							{ // Flow is "strong" (Peclet number above 2)
								if(flowVal > 0)
									// "Diffusion" rate is only due to flow
									mesoSubArray[curMesoSub].diffRateNeigh[curMolType][curNeighID] =
										2*flowProp;
								else
									// Transitions in this direction will not occur
									mesoSubArray[curMesoSub].diffRateNeigh[curMolType][curNeighID] =	0.;
							} else
							{
								mesoSubArray[curMesoSub].diffRateNeigh[curMolType][curNeighID] +=
									flowProp;
							}								
						}
					}
				}
			}
		}
	}
}

// Reset propensities and reaction times for all subvolumes
void resetMesoSubArray(const uint32_t numMesoSub,
	struct mesoSubvolume3D mesoSubArray[],
	const struct subvolume3D subvolArray[],
	const unsigned short NUM_MOL_TYPES,
	const unsigned short MAX_RXNS,
	const short NUM_REGIONS,
	struct region regionArray[])
{
	uint32_t curMeso, curSub;
	unsigned short curNeigh, curMolType, curRegion, destRegion;
	short curRxn;
	double curDiffRate;
	
	for(curMeso = 0; curMeso < numMesoSub; curMeso++)
	{
		curSub = mesoSubArray[curMeso].subID;
		curRegion = subvolArray[curSub].regionID;
		
		mesoSubArray[curMeso].totalProp = 0.;
		
		// Diffusion reactions
		for(curNeigh = 0; curNeigh < subvolArray[curSub].num_neigh; curNeigh++)
		{
			for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
			{
				destRegion = subvolArray[subvolArray[curSub].neighID[curNeigh]].regionID;
				if(curRegion == destRegion)
					curDiffRate = regionArray[curRegion].diffRate[curMolType];
				else
					curDiffRate = mesoSubArray[curMeso].diffRateNeigh[curMolType][curNeigh];
				mesoSubArray[curMeso].rxnProp[curMolType*subvolArray[curSub].num_neigh+curNeigh] =
					curDiffRate * mesoSubArray[curMeso].num_mol[curMolType];
				mesoSubArray[curMeso].totalProp +=
					mesoSubArray[curMeso].rxnProp[curMolType*subvolArray[curSub].num_neigh+curNeigh];
			}
		}
		
		// Chemical reactions
		for(curRxn = 0; curRxn < regionArray[curRegion].numChemRxn; curRxn++)
		{
			switch(regionArray[curRegion].rxnOrder[curRxn])
			{
				case 0:
					// Reaction is 0th order. Propensity will not change
					// for duration of simulation
					mesoSubArray[curMeso].rxnProp[mesoSubArray[curMeso].firstChemRxn + curRxn] =
						regionArray[curRegion].rxnRate[curRxn];
					break;
				case 1:
					// Reaction is 1st order. Propensity will depend on
					// number of reacting molecules
					mesoSubArray[curMeso].rxnProp[mesoSubArray[curMeso].firstChemRxn + curRxn] =
						regionArray[curRegion].rxnRate[curRxn] *
						mesoSubArray[curMeso].num_mol[regionArray[curRegion].uniReactant[curRxn]];
					break;
				case 2:
					// Reaction is 2nd order. Propensity will depend on
					// number of both types of reacting molecules
					mesoSubArray[curMeso].rxnProp[mesoSubArray[curMeso].firstChemRxn + curRxn] =
						regionArray[curRegion].rxnRate[curRxn] *
						mesoSubArray[curMeso].num_mol[regionArray[curRegion].biReactants[curRxn][0]];
					if(regionArray[curRegion].biReactants[curRxn][0]
						== regionArray[curRegion].biReactants[curRxn][1])
					{ // Reactants are the same type
						if(mesoSubArray[curMeso].num_mol[regionArray[curRegion].biReactants[curRxn][0]] > 0)
						{ // There must be two reactants to have a reaction
							mesoSubArray[curMeso].rxnProp[mesoSubArray[curMeso].firstChemRxn + curRxn] *=
								(mesoSubArray[curMeso].num_mol[regionArray[curRegion].biReactants[curRxn][0]] - 1);
						} else
						{ // Propensity formula not valid; reaction not possible
							mesoSubArray[curMeso].rxnProp[mesoSubArray[curMeso].firstChemRxn + curRxn] = 0;
						}
					} else
					{ // Reactants are different types
						mesoSubArray[curMeso].rxnProp[mesoSubArray[curMeso].firstChemRxn + curRxn] *=
							mesoSubArray[curMeso].num_mol[regionArray[curRegion].biReactants[curRxn][1]];						
					}
					break;
			}
			mesoSubArray[curMeso].totalProp +=
				mesoSubArray[curMeso].rxnProp[mesoSubArray[curMeso].firstChemRxn + curRxn];
		}
		
		// Generate first-reaction time from propensity
		mesoSubArray[curMeso].t_rxn = mesoSubCalcTime(mesoSubArray, curMeso);
	}
}

// Update propensities and next reaction time of subvolume
// NOTE: heapMesoUpdate3D should be called IMMEDIATELY after (i.e., before another
// call to this function), otherwise the heap won't be properly sorted
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
	struct region regionArray[])
{
	double old_time, old_prop; // Old reaction time and propensity
	double delta_prop; // Change in propensity
	double curDiffRate; // Current diffusion rate between meso subvolumes in adjacent regions
	unsigned short i, j;
	unsigned short destRegion[subvolArray[curSub].num_neigh]; // Region(s) of the subvolume's neighbors
	// Find ID of subvolume in list of mesoscopic subvolumes
	uint32_t curMeso = subvolArray[curSub].mesoID;
	unsigned short curRegion = subvolArray[curSub].regionID;
	unsigned short curRxn, curFirstRxn, curSecondRxn;
	unsigned short reactantA, reactantB;
	unsigned short curNumMolTypes;
	uint64_t curMolChange[NUM_MOL_TYPES];
	bool bCurMolAdd[NUM_MOL_TYPES];
	
	double changeMagA, changeMagB; // Will be 0, +1, or -1
	double prop_change;
	
	// Update propensities associated with diffusion and first-order reactions
	old_time = mesoSubArray[curMeso].t_rxn;
	old_prop = mesoSubArray[curMeso].totalProp;
	for(i = 0; i < subvolArray[curSub].num_neigh; i++)
	{
		destRegion[i] = subvolArray[subvolArray[curSub].neighID[i]].regionID;
	}
	
	if(bChemRxn)
	{ // Chemical reaction, so any change in molecules is "possible"		
		curNumMolTypes = NUM_MOL_TYPES;
		for(j = 0; j < NUM_MOL_TYPES; j++)
		{
			curMolChange[j] = numMolChange[j];
			bCurMolAdd[j] = bMolAdd[j];
		}
	} else
	{
		// Change in molecules must be placement of one molecule type or diffusion
		curNumMolTypes = 1;
		for(j = 0; j < NUM_MOL_TYPES; j++)
		{
			if(j == singleMolID)
			{
				curMolChange[j] = numMolChange[0];
				bCurMolAdd[j] = bMolAdd[0];
			} else
			{
				curMolChange[j] = 0;
			}
		}
	}
	
	for(j = 0; j < curNumMolTypes; j++)
	{
		if(!bChemRxn)
			j = singleMolID;
		
		if(curMolChange[j] > 0)
		{	// Number of j molecules changed.
			// Need to update corresponding propensities.
			
			// Updating Diffusion propensities
			for(i = 0; i < subvolArray[curSub].num_neigh; i++)
			{
				if(curRegion == destRegion[i] && !regionArray[curRegion].bFlow)
					curDiffRate = regionArray[curRegion].diffRate[j];
				else
					curDiffRate = mesoSubArray[curMeso].diffRateNeigh[j][i];
				
				mesoSubArray[curMeso].rxnProp[j*subvolArray[curSub].num_neigh+i]
					= curDiffRate*mesoSubArray[curMeso].num_mol[j];
			}
					
			// Update Chemical reaction propensities
			for(curFirstRxn = 0;
				curFirstRxn < regionArray[curRegion].numFirstRxn; curFirstRxn++)
			{
				curRxn = regionArray[curRegion].firstRxn[curFirstRxn];
				if(regionArray[curRegion].uniReactant[curRxn] == j)
				{ // Current reaction is 1st order and jth molecule is reactant
			
					mesoSubArray[curMeso].rxnProp[mesoSubArray[curMeso].firstChemRxn+curRxn]
						= regionArray[curRegion].rxnRate[curRxn]*mesoSubArray[curMeso].num_mol[j];
				}
			}			
			
			mesoSubArray[curMeso].totalProp =
				updateTotalProp(mesoSubArray[curMeso].rxnProp,
				regionArray[curRegion].numChemRxn + mesoSubArray[curMeso].firstChemRxn);
	
			
		}
	}
	
	// Update propensities associated with second order reactions
	for(curSecondRxn = 0; curSecondRxn < regionArray[curRegion].numSecondRxn; curSecondRxn++)
	{
		curRxn = regionArray[curRegion].secondRxn[curSecondRxn];
		
		reactantA = regionArray[curRegion].biReactants[curRxn][0];
		reactantB = regionArray[curRegion].biReactants[curRxn][1];
		prop_change = regionArray[curRegion].rxnRate[curRxn];
		
		if(reactantA == reactantB)
		{ // Both reactants are the same
			if(curMolChange[reactantA] > 0)
			{ // Propensity must be updated
				
				mesoSubArray[curMeso].rxnProp[mesoSubArray[curMeso].firstChemRxn+curRxn] =
					regionArray[curRegion].rxnRate[curRxn] *
					mesoSubArray[curMeso].num_mol[reactantA] *
					(mesoSubArray[curMeso].num_mol[reactantA] - 1);
				mesoSubArray[curMeso].totalProp =
					updateTotalProp(mesoSubArray[curMeso].rxnProp,
					regionArray[curRegion].numChemRxn + mesoSubArray[curMeso].firstChemRxn);
			}
		} else if(curMolChange[reactantA] > 0 || curMolChange[reactantB] > 0)
		{ // Propensity must be updated
			mesoSubArray[curMeso].rxnProp[mesoSubArray[curMeso].firstChemRxn+curRxn] =
				regionArray[curRegion].rxnRate[curRxn] *
				mesoSubArray[curMeso].num_mol[reactantA] *
				mesoSubArray[curMeso].num_mol[reactantB];
			mesoSubArray[curMeso].totalProp =
				updateTotalProp(mesoSubArray[curMeso].rxnProp,
				regionArray[curRegion].numChemRxn + mesoSubArray[curMeso].firstChemRxn);
		}
		
	}
	
	// Update reaction time of subvolume
	if(bAllowCorrectedTime && isfinite(mesoSubArray[curMeso].t_rxn))
	{
		// Subvolume had valid t_rxn. Update appropriately
		// This adjustment is much faster than generating a new time
		mesoSubArray[curMeso].t_rxn = tCur +
			old_prop / mesoSubArray[curMeso].totalProp * (old_time - tCur);
	} else
	{
		mesoSubArray[curMeso].t_rxn = tCur +
			mesoSubCalcTime(mesoSubArray, curMeso);
	}
}

// Update propensities, next reaction times, and heap location of subvolumes
// along micro/meso interface.
// This function is meant to be called AFTER validateMolecules() in
// micro_molecule.c, which populates the bNeedUpdate and numMolFromMicro
// members of the array of regionArray structures
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
	bool heap_childValid[][2])
{
	short curRegion, neighRegion;
	unsigned short curMolType;
	uint32_t curSub;
	
	for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
	{
		if(regionArray[curRegion].spec.bMicro)
			continue;	// Only need to check mesoscopic regions
		
		for(neighRegion = 0; neighRegion < NUM_REGIONS; neighRegion++)
		{
			if(!regionArray[curRegion].isRegionNeigh[neighRegion])
				continue; // Other region isn't a neighbor
			
			if(regionArray[neighRegion].spec.bMicro)
			{	// Other region is a microscopic neighbor
				for(curSub = 0;
					curSub < regionArray[curRegion].numSubRegionNeigh[neighRegion];
					curSub++)
				{ // Check each subvolume along boundary with neighboring region
					if(regionArray[curRegion].bNeedUpdate[neighRegion][curSub])
					{
						// Reset signals for new molecules from micro regime
						for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
						{
							updateMesoSub(mesoSubArray[regionArray[curRegion].neighID[neighRegion][curSub]].subID, false, 
								(uint64_t []){regionArray[curRegion].numMolFromMicro[neighRegion][curSub][curMolType]},
								bTrue, curMolType,
								true, numMesoSub, mesoSubArray, subvolArray,
								tCur, NUM_REGIONS, NUM_MOL_TYPES, regionArray);
							regionArray[curRegion].numMolFromMicro[neighRegion][curSub][curMolType] = 0;
						}
						regionArray[curRegion].bNeedUpdate[neighRegion][curSub] = false;
						heapMesoUpdate(numMesoSub, mesoSubArray, heap_subvolID,
							mesoSubArray[regionArray[curRegion].neighID[neighRegion][curSub]].heapID, heap_childID, heap_childValid);
							
						regionArray[curRegion].bNeedUpdate[neighRegion][curSub] = false;
					}
				}
			}
		}
	}
}

// Sum terms in reaction propensity vector
double updateTotalProp(const double rxnProp[],
	const unsigned short numChemRxn)
{
	double totalProp = 0.;
	for(unsigned short i = 0; i < numChemRxn; i++)
	{
		totalProp += rxnProp[i];
	}
	return totalProp;
}

// Determine next relative subvolume reaction time
double mesoSubCalcTime(struct mesoSubvolume3D mesoSubArray[],
	const uint32_t ID)
{
	return generateExponential(1) / mesoSubArray[ID].totalProp;
}

// Build 2D array listing children of heap elements
void heapMesoFindChildren(const uint32_t numSub,
	uint32_t heap_childID[][2],
	bool heap_childValid[][2])
{
	uint32_t i;
	
	for(i = 0; i < numSub; i++){ // parent == i
		heap_childID[i][0] = (uint32_t) 2*i + 1;
		heap_childID[i][1] = (uint32_t) 2*i + 2;
		
		if (heap_childID[i][1] < numSub){ // Both children valid
			heap_childValid[i][0] = true;
			heap_childValid[i][1] = true;
		} else if (heap_childID[i][0] < numSub){ // First child valid
			heap_childValid[i][0] = true;
			heap_childValid[i][1] = false;
		} else { // No valid children
			heap_childValid[i][0] = false;
			heap_childValid[i][1] = false;
		}			
	}
	
	return;
}

// Free memory of heap arrays
void heapMesoDelete(const uint32_t numElements,
	uint32_t heap_subvolID[],
	uint32_t heap_childID[][2],
	bool heap_childValid[][2])
{	
	if(heap_subvolID != NULL) free(heap_subvolID);
	if(heap_childID != NULL) free(heap_childID);
	if(heap_childValid != NULL) free(heap_childValid);
}

// Build Heap of Next Subvolume Reaction Times
void heapMesoBuild(const uint32_t numSub,
	struct mesoSubvolume3D mesoSubArray[],
	uint32_t heap_subvolID[],
	const unsigned int num_heap_levels,
	uint32_t heap_childID[][2],
	bool heap_childValid[][2])
{
	int i,j;
	uint32_t levelFirst, levelLast; // First and last elements in current level
	uint32_t parent, parentNew, parentOld;
	
	// Initialize heap with ordered (but unsorted) indices
	for(i=0;i<numSub;i++){
		heap_subvolID[i] = i;
		mesoSubArray[i].heapID = i;
	}
	
	// Make heap a "min-heap" using next reaction times
	for(i = num_heap_levels-2; i >= 0; i--){ // For each level in heap
		levelFirst = 1;
		for(j=0; j<i; j++) levelFirst *= 2;
		levelLast = 2*levelFirst - 2; // Last element is 2^(i+1)-2
		levelFirst--; // First element is 2^i - 1
		for(parent = levelFirst; parent <= levelLast; parent++){
			// Compare element with its two children
			parentOld = parent;
			parentNew = heapMesoCompareDown(numSub, mesoSubArray, heap_subvolID,
				parent, heap_childID, heap_childValid);
			while(parentNew != parentOld){
				// Keep comparing with children until this parent is the min
				parentOld = parentNew;
				parentNew = heapMesoCompareDown(numSub, mesoSubArray, heap_subvolID,
					parentOld, heap_childID, heap_childValid);
			}
		}
	}
}

// Update Placement of Single Heap Element Based on Updated Value
uint32_t heapMesoUpdate(const uint32_t numSub,
	struct mesoSubvolume3D mesoSubArray[],
	uint32_t heap_subvolID[],
	const uint32_t heapID,
	uint32_t heap_childID[][2],
	bool heap_childValid[][2])
{
	uint32_t newID, oldID; // Updated current and previous placement in heap
	
	// See if element needs to move "down" (i.e., lower priority)
	oldID = heapID;
	newID = heapMesoCompareDown(numSub, mesoSubArray, heap_subvolID, heapID,
		heap_childID, heap_childValid);
	while(newID != oldID){
		oldID = newID;
		newID = heapMesoCompareDown(numSub, mesoSubArray, heap_subvolID, oldID,
			heap_childID, heap_childValid);
	}
	if (newID == heapID && heapID > 0){
		// Element did not move down and it is not already at the top of the heap.
		// See if it needs to move "up" (i.e., higher priority)
		newID = heapMesoCompareUp(numSub, mesoSubArray, heap_subvolID, heapID);
		while(newID != oldID){
			oldID = newID;
			newID = heapMesoCompareUp(numSub, mesoSubArray, heap_subvolID, oldID);
		}
	}
	return newID;
}

// Compare parent node with its children. Swap if parent does not have smallest value
// Return (possibly new) location of parent node
uint32_t heapMesoCompareDown(const uint32_t numSub,
	struct mesoSubvolume3D mesoSubArray[],
	uint32_t heap_subvolID[],
	const uint32_t parent,
	uint32_t heap_childID[][2],
	bool heap_childValid[][2])
{
		
	uint32_t child1 = heap_childID[parent][0];
	uint32_t child2 = heap_childID[parent][1];
	
	if (heap_childValid[parent][1]){
		if (mesoSubArray[heap_subvolID[parent]].t_rxn <
			mesoSubArray[heap_subvolID[child1]].t_rxn){
			if (mesoSubArray[heap_subvolID[parent]].t_rxn >
				mesoSubArray[heap_subvolID[child2]].t_rxn){
				// child2 is smallest
				heapMesoSwap(mesoSubArray,heap_subvolID,parent,child2);
				return child2;
			} // else parent is smallest (do nothing)
		} else if (mesoSubArray[heap_subvolID[child2]].t_rxn <
			mesoSubArray[heap_subvolID[child1]].t_rxn){
			// child2 is smallest
			heapMesoSwap(mesoSubArray,heap_subvolID,parent,child2);
			return child2;
		} else {
			// child1 is smallest
			heapMesoSwap(mesoSubArray,heap_subvolID,parent,child1);
			return child1;
		}
	} else if (heap_childValid[parent][0] && mesoSubArray[heap_subvolID[parent]].t_rxn >
		mesoSubArray[heap_subvolID[child1]].t_rxn) {
		// child1 is smallest
		heapMesoSwap(mesoSubArray,heap_subvolID,parent,child1);
		return child1;
	}
	return parent;
}

// Compare node with its parent. Swap if parent has a larger value
// Return (possibly new) location of node
uint32_t heapMesoCompareUp(const uint32_t numSub,
	struct mesoSubvolume3D mesoSubArray[],
	uint32_t heap_subvolID[],
	const uint32_t child){
	
	if (child == 0) // Node is at top of heap. No parent exists
		return child;
	
	uint32_t parent = child/2;
	if(child % 2 == 0) // Element is even
		parent--;
	
	if (mesoSubArray[heap_subvolID[parent]].t_rxn >
		mesoSubArray[heap_subvolID[child]].t_rxn){
		// Child is smaller. Swap
		heapMesoSwap(mesoSubArray,heap_subvolID,parent,child);
		return parent;
	} else // Parent is smaller. Do not swap
		return child;
}

// Swap the positions of two elements in the heap
void heapMesoSwap(struct mesoSubvolume3D mesoSubArray[],
	uint32_t heap_subvolID[],
	uint32_t index1,
	uint32_t index2){
	
	uint32_t tempID;
	
	// Update heap IDs associated with the subvolumes whose heap positions are being swapped
	mesoSubArray[heap_subvolID[index1]].heapID = index2;
	mesoSubArray[heap_subvolID[index2]].heapID = index1;
	
	tempID = heap_subvolID[index1];
	heap_subvolID[index1] = heap_subvolID[index2];
	heap_subvolID[index2] = tempID;
}
