/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * subvolume.c - 	structure for storing subvolume properties. Simulation
 *					environment is partitioned into subvolumes
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
 * - added accommodation of spherical (microscopic) subvolumes and actors
 * - moved counting and IDing of neighbor subvolumes in different regions to a separate
 * function to simplify code reuse
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
#include <limits.h> // For SHRT_MAX
#include "subvolume.h" // for "Public" declarations

//
// "Private" Declarations
//

//
// Definitions
//

// Initialize Subvolume Structure
/* A call to this function should have a corresponding call to delete_subvol_array2D,
* usually AFTER build_subvol_array2D has been called */
void allocateSubvolArray(const uint32_t numSub,
	struct subvolume3D ** subvolArray)
{	
	*subvolArray = malloc(numSub*sizeof(struct subvolume3D));
		
	if(*subvolArray == NULL){
		fprintf(stderr, "ERROR: Memory allocation for array of subvolume structures.\n");
		exit(EXIT_FAILURE);
	}
}

// Allocate arrays to help initialization of subvolume and region structure members
void allocateSubvolHelper(const uint32_t numSub,
	uint32_t (** subCoorInd)[3],
	uint32_t ***** subID,
	const short NUM_REGIONS,
	struct region regionArray[])
{	
	short i;
	uint32_t x,y,z;
	*subCoorInd = malloc(numSub*sizeof(uint32_t [3]));
	*subID = malloc(NUM_REGIONS*sizeof(*subID));
		
	if(subCoorInd == NULL || *subID == NULL){
		fprintf(stderr, "ERROR: Memory allocation for temporary subvolume information.\n");
		exit(EXIT_FAILURE);
	}
	
	for(i = 0; i < NUM_REGIONS; i++)
	{
		(*subID)[i] = malloc(regionArray[i].spec.numX*sizeof(*subID));
		if((*subID)[i] == NULL){
			fprintf(stderr, "ERROR: Memory allocation for temporary subvolume information.\n");
			exit(EXIT_FAILURE);
		}
		for(x = 0; x < regionArray[i].spec.numX; x++)
		{
			(*subID)[i][x] = malloc(regionArray[i].spec.numY*sizeof(*subID));
			if((*subID)[i][x] == NULL){
				fprintf(stderr, "ERROR: Memory allocation for temporary subvolume information.\n");
				exit(EXIT_FAILURE);
			}
			for(y = 0; y < regionArray[i].spec.numY; y++)
			{
				(*subID)[i][x][y] = malloc(regionArray[i].spec.numZ*sizeof(*subID));
				if((*subID)[i][x][y] == NULL){
					fprintf(stderr, "ERROR: Memory allocation for temporary subvolume information.\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}
}

// Free Memory of Subvolume Structure
/* This function is needed to free the memory allocated to each
* array member of each structure element before releasing the memory
* allocated to the structure array itself. It should be called after BOTH
* allocate_subvol_array3D AND build_subvol_array3D, but is written to avoid
* any attempts to free memory that was not previously allocated */
void delete_subvol_array3D(const uint32_t numSub,
	struct subvolume3D subvolArray[],
	const unsigned short NUM_MOL_TYPES,
	const short NUM_REGIONS,
	struct region regionArray[])
{
	uint32_t curSub = 0;
	uint32_t curMolType = 0;
	
	if(subvolArray == NULL) // There is no allocated memory to release
		return;
	
	for(curSub; curSub < numSub; curSub++)
	{
		if(subvolArray[curSub].neighID != NULL)	free(subvolArray[curSub].neighID);
		if(!regionArray[subvolArray[curSub].regionID].spec.bMicro)
		{
			if(subvolArray[curSub].num_mol != NULL)	free(subvolArray[curSub].num_mol);
		}
		if(NUM_REGIONS > 1 && subvolArray[curSub].bBoundary
			&& !regionArray[subvolArray[curSub].regionID].spec.bMicro
			&& subvolArray[curSub].diffRateNeigh != NULL)
		{
			for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
			{
				if(subvolArray[curSub].diffRateNeigh[curMolType] != NULL)
					free(subvolArray[curSub].diffRateNeigh[curMolType]);
			}
			free(subvolArray[curSub].diffRateNeigh);
		}
	}
	
	free(subvolArray);
}

void deleteSubvolHelper(uint32_t subCoorInd[][3],
	uint32_t **** subID,
	const short NUM_REGIONS,
	struct region regionArray[])
{	
	short i;
	uint32_t x,y;
	
	if(subCoorInd != NULL) free(subCoorInd);
	
	if(subID != NULL)
	{
		for(i = 0; i < NUM_REGIONS; i++)
		{
			if(subID[i] != NULL)
			{
				for(x = 0; x < regionArray[i].spec.numX; x++)
				{
					if(subID[i][x] != NULL)
					{
						for(y = 0; y < regionArray[i].spec.numY; y++)
						{
							if(subID[i][x][y] != NULL) free(subID[i][x][y]);
						}
						free(subID[i][x]);
					}
				}
				free(subID[i]);
			}
		}
		free(subID);
	}
}

// Construct the array of structures with details of each subvolume
void buildSubvolArray3D(const uint32_t numSub,
	uint32_t * numMesoSub,
	struct subvolume3D subvolArray[],
	const struct spec_region3D subvol_spec[],
	struct region regionArray[],
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const unsigned short MAX_RXNS, 
	const double SUBVOL_BASE_SIZE,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	uint32_t subCoorInd[numSub][3],
	uint32_t **** subID)
{	
	short int i,j, curRegion, neighRegion, sphRegion, rectRegion; // Current Region
	unsigned short curMolType;
	uint32_t curX, curY, curZ; // Coordinates of current subvolume within current region
	bool bRealSub;
	uint32_t curID = 0; // Current Subvolume ID
	uint32_t curBoundID = 0; // Current subvolume in region boundary list
	uint32_t curMesoID = 0; // Current Mesoscopic Subvolume ID
	uint32_t curNeighID = 0; // Current Subvolume neighbour ID
	uint32_t neighID = 0; // Subvolume neighbour ID in master subvolume list
	
	uint32_t sphSub, rectSub; // IDs of subvolumes that are spherical, rectangular
	unsigned short numFaceSph = 0; // Number of faces of subvolume that border a spherical region
	
	// Initialize parameters to efficiently determine neighboring subvolumes in same region
	uint32_t * curSubNeigh = malloc(numSub*sizeof(uint32_t));
	
	bool bMesoToMicro; // Indicate whether we also need to populate region virtual neighbor members
	
	if(curSubNeigh == NULL)
	{
		fprintf(stderr, "ERROR: Memory allocation to temporarily store number of neighbors of each subvolume that are in same region.\n");
		exit(EXIT_FAILURE);
	}
	for(curID = 0; curID < numSub; curID++)
	{
		curSubNeigh[curID] = 0UL;
	}
	
	double curSubBound[6]; // Boundary of current subvolume
	double curNeighBound[6]; // Boundary of prospective neighbor subvolume
	double boundOverlap[6]; // Overlap area of adjacent subvolumes
	for(i=0; i < 6; i++)
	{
		curSubBound[i] = 0.;
		curNeighBound[i] = 0.;
		boundOverlap[i] = 0.;
	}
	unsigned short adjDirection = 0;
	double boundAdjError = SUBVOL_BASE_SIZE * SUB_ADJ_RESOLUTION;
	double subHalfSize[NUM_REGIONS]; // Half of the actual size of subvolumes in a region
	
	double h_i, h_j; // Subvolume sizes (used for finding transition rates)
			
	// Store basic subvolume information based on the specification
	// Populate subvolArray based on subvol_spec
	for(i = 0,curID=0; i < NUM_REGIONS; i++)
	{
		if(regionArray[i].spec.shape == SPHERE) // spheres always have one valid subvolume
		{
			subID[i][0][0][0] = curID;
			subvolArray[curID].regionID = i;
			subCoorInd[curID][0] = regionArray[i].boundary[0];
			subCoorInd[curID][1] = regionArray[i].boundary[1];
			subCoorInd[curID][2] = regionArray[i].boundary[2];
			subvolArray[curID].mesoID = UINT32_MAX; // Spheres must be MICRO
			subvolArray[curID].num_neigh = 0;
			subvolArray[curID].bBoundary = true;
			curID++;
			continue;
		}
		
		for(curZ = 0; curZ < subvol_spec[i].numZ; curZ++)
		{
			for(curY = 0; curY < subvol_spec[i].numY; curY++)
			{
				for(curX = 0; curX < subvol_spec[i].numX; curX++)
				{
					bRealSub = true;
					
					subCoorInd[curID][0] = curX;
					subCoorInd[curID][1] = curY;
					subCoorInd[curID][2] = curZ;
					
					// Check that current location is not in nested region
					for(j = 0; j < regionArray[i].numChildren; j++)
					{
						// Is current location actually within a child?
						if(regionArray[regionArray[i].childrenID[j]].spec.shape == SPHERE)
						{
							// Need real coordinates of subvolume
							findSubvolCoor(curSubBound, regionArray[i], subCoorInd[curID]);
							if(bBoundarySurround(RECTANGULAR_BOX, curSubBound,
								SPHERE, regionArray[regionArray[i].childrenID[j]].boundary, 0.))
							{
								// subvolume is entirely within spherical child
								bRealSub = false;
								subID[i][curX][curY][curZ] = UINT32_MAX;
								continue;
							}
						} else if(curX >= regionArray[i].childrenCoor[j][0]
							&& curX <= regionArray[i].childrenCoor[j][1]
							&& curY >= regionArray[i].childrenCoor[j][2]
							&& curY <= regionArray[i].childrenCoor[j][3]
							&& curZ >= regionArray[i].childrenCoor[j][4]
							&& curZ <= regionArray[i].childrenCoor[j][5])
						{
							bRealSub = false;
							subID[i][curX][curY][curZ] = UINT32_MAX;
							continue;
						}
					}
						
					if(!bRealSub)
					{	// We are currently inside a nested region.
						// Skip subvolume creation
						continue;
					}
					
					subID[i][curX][curY][curZ] = curID;
					
					subvolArray[curID].regionID = i;
					
					// Assign memory to array parameters and confirm success
					
					if(!regionArray[i].spec.bMicro)
					{ // Parameter only needed for mesoscopic regions
						subvolArray[curID].num_mol = malloc(NUM_MOL_TYPES*sizeof(uint64_t));
						if(subvolArray[curID].num_mol == NULL)
						{
							fprintf(stderr, "ERROR: Memory allocation for number of molecules in subvolume %" PRIu32 ".\n", curID);
							exit(EXIT_FAILURE);
						}
						subvolArray[curID].mesoID = curMesoID;
						curMesoID++;
					} else{
						subvolArray[curID].mesoID = UINT32_MAX;				
					}
					
					// Make sure remaining parameters are 0
					subvolArray[curID].num_neigh = 0;
					subvolArray[curID].bBoundary = false;
					
					// Find if an interior subvolume borders a nested region
					for(j = 0; j < regionArray[i].numChildren; j++)
					{
						if(regionArray[regionArray[i].childrenID[j]].spec.shape == SPHERE)
						{
							// Sphere and subvolume are adjacent if they intersect
							// We already determined curSubBound above
							if(bBoundaryIntersect(RECTANGULAR_BOX, curSubBound,
								SPHERE, regionArray[regionArray[i].childrenID[j]].boundary,
								regionArray[i].actualSubSize))
							{
								subvolArray[curID].bBoundary = true;
								continue; // No need to keep comparing with other children
							}							
						} else if(((curX == regionArray[i].childrenCoor[j][0]-1 ||
							curX == regionArray[i].childrenCoor[j][1]+1)
							&& curY >= regionArray[i].childrenCoor[j][2]
							&& curY <= regionArray[i].childrenCoor[j][3]
							&& curZ >= regionArray[i].childrenCoor[j][4]
							&& curZ <= regionArray[i].childrenCoor[j][5])
							||
							((curY == regionArray[i].childrenCoor[j][2]-1 ||
							curY == regionArray[i].childrenCoor[j][3]+1)
							&& curX >= regionArray[i].childrenCoor[j][0]
							&& curX <= regionArray[i].childrenCoor[j][1]
							&& curZ >= regionArray[i].childrenCoor[j][4]
							&& curZ <= regionArray[i].childrenCoor[j][5])
							||
							((curZ == regionArray[i].childrenCoor[j][4]-1 ||
							curZ == regionArray[i].childrenCoor[j][5]+1)
							&& curX >= regionArray[i].childrenCoor[j][0]
							&& curX <= regionArray[i].childrenCoor[j][1]
							&& curY >= regionArray[i].childrenCoor[j][2]
							&& curY <= regionArray[i].childrenCoor[j][3]))
						{
							subvolArray[curID].bBoundary = true;
							continue; // No need to keep comparing with other children
						}
					}
						
					curID++;
				}
			}
		}
	}
		
	// Store basic subvolume information based on the specification
	// Populate subvolArray based on subvol_spec
	for(i = 0,curID=0; i < NUM_REGIONS; i++)
	{
		if(regionArray[i].spec.shape == SPHERE)
		{
			// Already know one subvolume in sphere is valid
			curID++;
			continue;
		}
		
		for(curZ = 0; curZ < subvol_spec[i].numZ; curZ++)
		{
			for(curY = 0; curY < subvol_spec[i].numY; curY++)
			{
				for(curX = 0; curX < subvol_spec[i].numX; curX++)
				{
					if(subID[i][curX][curY][curZ] < UINT32_MAX)
					{	// Current subvolume is valid
					
						if(subCoorInd[curID][0] > 0)
						{ // Neighbor is 1 "x" index "down"
							if(subID[i][curX-1][curY][curZ] < UINT32_MAX)
								subvolArray[curID].num_neigh++;
						} else subvolArray[curID].bBoundary = true;
						if(subCoorInd[curID][0] < subvol_spec[i].numX-1)
						{ // Neighbor is 1 "x" index "up"
							if(subID[i][curX+1][curY][curZ] < UINT32_MAX)
								subvolArray[curID].num_neigh++;
						} else subvolArray[curID].bBoundary = true;
						if(subCoorInd[curID][1] > 0)
						{ // Neighbor is 1 "y" index "down"
							if(subID[i][curX][curY-1][curZ] < UINT32_MAX)
								subvolArray[curID].num_neigh++;
						} else subvolArray[curID].bBoundary = true;
						if(subCoorInd[curID][1] < subvol_spec[i].numY-1)
						{ // Neighbor is 1 "y" index "up"
							if(subID[i][curX][curY+1][curZ] < UINT32_MAX)
								subvolArray[curID].num_neigh++;
						} else subvolArray[curID].bBoundary = true;
						if(subCoorInd[curID][2] > 0)
						{ // Neighbor is 1 "z" index "down"
							if(subID[i][curX][curY][curZ-1] < UINT32_MAX)
								subvolArray[curID].num_neigh++;
						} else subvolArray[curID].bBoundary = true;
						if(subCoorInd[curID][2] < subvol_spec[i].numZ-1)
						{ // Neighbor is 1 "z" index "up"
							if(subID[i][curX][curY][curZ+1] < UINT32_MAX)
								subvolArray[curID].num_neigh++;
						} else subvolArray[curID].bBoundary = true;
						
						curID++;
					}
					
				}
			}
		}
	}
	
	// Assign total number of mesoscopic subvolumes
	*numMesoSub = curMesoID;
	
	// Determine the threshold distances for any two subvolumes to be neighbors
	for(i = 0; i < NUM_REGIONS; i++)
	{
		subHalfSize[i] = regionArray[i].actualSubSize/2;
	}
	
	// Determine the neighbors of each subvolume that are in other regions
	printf("Finding Neighbours of Each Subvolume...\n");
	if(NUM_REGIONS > 1)
	{ // Only need to continue if there is more than one region

		for(curID = 0; curID < numSub; curID++)
		{ // For each subvolume
			if(!subvolArray[curID].bBoundary)
				continue; // This subvolume is not along the boundary
			curRegion = subvolArray[curID].regionID;
			if(curRegion == NUM_REGIONS-1)
				break; // There are no more regions to compare with
			for(curNeighID = regionArray[curRegion+1].firstID;
				curNeighID < numSub; curNeighID++)
			{ // For every remaining subvolume in a different region
				if(!subvolArray[curNeighID].bBoundary)
					continue; // Neighbor is not along region boundary
				
				neighRegion = subvolArray[curNeighID].regionID;
				if(!regionArray[curRegion].isRegionNeigh[neighRegion])
					continue; // Subvolumes are not in neighbouring regions
				
				// Subvolumes are in neighboring regions and each is along its region
				// boundary
				if(checkSubvolNeigh(regionArray, curRegion, neighRegion, &sphRegion,
					&rectRegion, subvolArray, curID, curNeighID, &sphSub, &rectSub,
					numSub, subCoorInd, boundAdjError, &adjDirection, subHalfSize,
					curSubBound, curNeighBound, &numFaceSph))
				{
					// Subvolumes are neighbors. Record IDs
					if (numFaceSph > 0)
					{
						// One subvolume is in a box while the other is in a sphere.
						// The subvolumes can be neighbors along multiple faces
						subvolArray[sphSub].num_neigh++;
						if (regionArray[rectRegion].spec.bMicro)
							subvolArray[rectSub].num_neigh++;
						else
							subvolArray[rectSub].num_neigh += numFaceSph;
						numFaceSph = 0; // Reset value
					} else
					{
						subvolArray[curID].num_neigh++;
						subvolArray[curNeighID].num_neigh++;
					}
				}
				
			}
		}
	}
	
	// Allocate space to store IDs of each subvolumes neighbors
	for(curID = 0; curID < numSub; curID++)
	{
		subvolArray[curID].neighID =
			malloc(subvolArray[curID].num_neigh*sizeof(uint32_t));
		if(subvolArray[curID].neighID == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation for IDs of neighbors of subvolume %" PRIu32 ".\n", curID);
			exit(EXIT_FAILURE);
		}
	}
		
	// Find and store IDs of each subvolume's neighbors
	for(curID = 0; curID < numSub; curID++)
	{
		// Calculate neighbors in the same region
		curRegion = subvolArray[curID].regionID;
				
		if(regionArray[curRegion].spec.shape != SPHERE)
		{
			if(subCoorInd[curID][0] > 0)
			{ // Neighbor is 1 "x" index "down"
				if(subID[curRegion][subCoorInd[curID][0]-1][subCoorInd[curID][1]][subCoorInd[curID][2]] < UINT32_MAX)
				subvolArray[curID].neighID[curSubNeigh[curID]++] =
					subID[curRegion][subCoorInd[curID][0]-1][subCoorInd[curID][1]][subCoorInd[curID][2]];
			}
			if(subCoorInd[curID][0] < subvol_spec[curRegion].numX-1)
			{ // Neighbor is 1 "x" index "up"
				if(subID[curRegion][subCoorInd[curID][0]+1][subCoorInd[curID][1]][subCoorInd[curID][2]] < UINT32_MAX)
				subvolArray[curID].neighID[curSubNeigh[curID]++] =
					subID[curRegion][subCoorInd[curID][0]+1][subCoorInd[curID][1]][subCoorInd[curID][2]];
			}
			if(subCoorInd[curID][1] > 0)
			{ // Neighbor is 1 "y" index "down"
				if(subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]-1][subCoorInd[curID][2]] < UINT32_MAX)
				subvolArray[curID].neighID[curSubNeigh[curID]++] =
					subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]-1][subCoorInd[curID][2]];
			}
			if(subCoorInd[curID][1] < subvol_spec[curRegion].numY-1)
			{ // Neighbor is 1 "y" index "up"
				if(subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]+1][subCoorInd[curID][2]] < UINT32_MAX)
				subvolArray[curID].neighID[curSubNeigh[curID]++] =
					subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]+1][subCoorInd[curID][2]];
			}
			if(subCoorInd[curID][2] > 0)
			{ // Neighbor is 1 "z" index "down"
				if(subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]][subCoorInd[curID][2]-1] < UINT32_MAX)
				subvolArray[curID].neighID[curSubNeigh[curID]++] =
					subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]][subCoorInd[curID][2]-1];
			}
			if(subCoorInd[curID][2] < subvol_spec[curRegion].numZ-1)
			{ // Neighbor is 1 "z" index "up"
				if(subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]][subCoorInd[curID][2]+1] < UINT32_MAX)
				subvolArray[curID].neighID[curSubNeigh[curID]++] =
					subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]][subCoorInd[curID][2]+1];
			}
		}
		
		// Find neighbors in regions with higher indices
		if(curRegion == NUM_REGIONS-1)
			continue; // There are no higher regions to compare with
		if(!subvolArray[curID].bBoundary)
			continue; // This subvolume is not along the boundary
		for(curNeighID = regionArray[curRegion+1].firstID;
			curNeighID < numSub; curNeighID++)
		{ // For every remaining subvolume in a different region
			if(!subvolArray[curNeighID].bBoundary)
				continue; // Neighbor is not along region boundary
			
			neighRegion = subvolArray[curNeighID].regionID;
			if(!regionArray[curRegion].isRegionNeigh[neighRegion])
				continue; // Subvolumes are not in neighbouring regions
			
			// Subvolumes are in neighboring regions and each is along its region
			// boundary
			if(checkSubvolNeigh(regionArray, curRegion, neighRegion, &sphRegion,
				&rectRegion, subvolArray, curID, curNeighID, &sphSub, &rectSub,
				numSub, subCoorInd, boundAdjError, &adjDirection, subHalfSize,
				curSubBound, curNeighBound, &numFaceSph))
			{
				// Subvolumes are neighbors. Record IDs
				if (numFaceSph > 0)
				{
					// One subvolume is in a box while the other is in a sphere.
					// The subvolumes can be neighbors along multiple faces
					subvolArray[sphSub].neighID[curSubNeigh[sphSub]++] = rectSub;
					if (regionArray[rectRegion].spec.bMicro)
						subvolArray[rectSub].neighID[curSubNeigh[rectSub]++] = sphSub;
					else
						for(i = 0; i < numFaceSph; i++)
							subvolArray[rectSub].neighID[curSubNeigh[rectSub]++] = sphSub;
					numFaceSph = 0; // Reset value
				} else
				{
					subvolArray[curID].neighID[curSubNeigh[curID]++] = curNeighID;			
					subvolArray[curNeighID].neighID[curSubNeigh[curNeighID]++] = curID;
				}
			}
		}
	}
		
	if(NUM_REGIONS > 1)
	{ // Only need to enter if there is more than one region

		// Determine transition rates out of mesoscopic subvolumes that are along
		// the boundary of their respective region
		for(curID = 0; curID < numSub; curID++)
		{ // For each subvolume
			if(!subvolArray[curID].bBoundary || regionArray[subvolArray[curID].regionID].spec.bMicro)
				continue; // This subvolume is not meso and along the boundary
			
			// Allocate memory to store transition rate for this subvolume
			subvolArray[curID].diffRateNeigh =
				malloc(NUM_MOL_TYPES*sizeof(subvolArray[curID].diffRateNeigh));
			
			if(subvolArray[curID].diffRateNeigh == NULL){
				fprintf(stderr, "ERROR: Memory allocation for mesoscopic diffusion rates for subvolume %" PRIu32 ", which is alng the boundary or region %u.\n", curID, subvolArray[curID].regionID);
				exit(EXIT_FAILURE);
			}
			
			for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
			{
				subvolArray[curID].diffRateNeigh[curMolType] =
					malloc(subvolArray[curID].num_neigh*sizeof(double));
				if(subvolArray[curID].diffRateNeigh[curMolType] == NULL){
					fprintf(stderr, "ERROR: Memory allocation for mesoscopic diffusion rates for subvolume %" PRIu32 ", which is alng the boundary or region %u.\n", curID, subvolArray[curID].regionID);
					exit(EXIT_FAILURE);
				}
			}
			
			curRegion = subvolArray[curID].regionID;
			h_i = regionArray[curRegion].actualSubSize;
			
			findSubvolCoor(curSubBound, regionArray[curRegion], subCoorInd[curID]);
			
			for(curNeighID = 0; curNeighID < subvolArray[curID].num_neigh;
				curNeighID++)
			{
				// Find actual neighbor index and region
				neighID = subvolArray[curID].neighID[curNeighID];
				neighRegion = subvolArray[neighID].regionID;
				if (curRegion == neighRegion)
				{ // Transition rate only depends on source volume
					for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
					{
						subvolArray[curID].diffRateNeigh[curMolType][curNeighID] =
							DIFF_COEF[curRegion][curMolType]/h_i/h_i;
					}						
				} else
				{
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
						subvolArray[curID].diffRateNeigh[curMolType][curNeighID] =
							2*DIFF_COEF[curRegion][curMolType]/h_i/(h_i + h_j);
						if(fabs(boundOverlap[0] - boundOverlap[1]) > boundAdjError)
							subvolArray[curID].diffRateNeigh[curMolType][curNeighID] *=
								(boundOverlap[1] - boundOverlap[0])/h_i;
						if(fabs(boundOverlap[2] - boundOverlap[3]) > boundAdjError)
							subvolArray[curID].diffRateNeigh[curMolType][curNeighID] *=
								(boundOverlap[3] - boundOverlap[2])/h_i;
						if(fabs(boundOverlap[4] - boundOverlap[5]) > boundAdjError)
							subvolArray[curID].diffRateNeigh[curMolType][curNeighID] *=
								(boundOverlap[5] - boundOverlap[4])/h_i;
					}
				}
			}
		}
		
		// Initialize region knowledge of the subvolumes that are adjacent to it
		initializeRegionSubNeighbor(regionArray, subvol_spec,
			subHalfSize, NUM_REGIONS, NUM_MOL_TYPES, SUBVOL_BASE_SIZE,
			boundAdjError, subvolArray, subCoorInd);
	}
	
	free(curSubNeigh);
	
	return;
}

// Determine whether two subvolumes in neighboring regions are neighbors themselves
// Assert that each subvolume is along its own region's boundary
// If subvolumes are neighbors, then the direction from curID towards neighID
// is written to adjDirection
bool checkSubvolNeigh(struct region regionArray[],
	short int curRegion,
	short int neighRegion,
	short int * sphRegion,
	short int * rectRegion,
	struct subvolume3D subvolArray[],
	uint32_t curID,
	uint32_t curNeighID,
	uint32_t * sphSub,
	uint32_t * rectSub,
	uint32_t numSub,
	uint32_t subCoorInd[numSub][3],
	double boundAdjError,
	unsigned short * adjDirection,
	double subHalfSize[],
	double curSubBound[6],
	double curNeighBound[6],
	unsigned short * numFaceSph)
{
	unsigned short dirArray[6];

	if (regionArray[curRegion].spec.shape == RECTANGULAR_BOX
		&& regionArray[neighRegion].spec.shape == RECTANGULAR_BOX)
	{
		findSubvolCoor(curSubBound, regionArray[curRegion], subCoorInd[curID]);
		findSubvolCoor(curNeighBound, regionArray[neighRegion],
			subCoorInd[curNeighID]);
			
		// Are the two subvolumes close enough to be neighbours?
		if(bBoundaryAdjacent(RECTANGULAR_BOX, curSubBound, RECTANGULAR_BOX,
			curNeighBound, boundAdjError, adjDirection))
		{ // Subvolumes are neighbours. Record
			return true;
		}
	} else if (regionArray[curRegion].spec.shape == SPHERE
		&& regionArray[neighRegion].spec.shape == SPHERE)
	{
		// Two neighboring spherical regions must have neighboring subvolumes
		return true;
	} else if ((regionArray[curRegion].spec.shape == RECTANGULAR_BOX
		&& regionArray[neighRegion].spec.shape == SPHERE)
		|| (regionArray[curRegion].spec.shape == SPHERE
		&& regionArray[neighRegion].spec.shape == RECTANGULAR_BOX))
	{
		// These two regions must be in a parent/child relationship
		if(regionArray[curRegion].spec.shape == RECTANGULAR_BOX)
		{
			*rectSub = curID;
			*rectRegion = curRegion;
			*sphSub = curNeighID;
			*sphRegion = neighRegion;
		} else
		{
			*rectSub = curNeighID;
			*rectRegion = neighRegion;
			*sphSub = curID;
			*sphRegion = curRegion;
		}
		findSubvolCoor(curSubBound, regionArray[*rectRegion],
			subCoorInd[*rectSub]);
			
		
		if(regionArray[*rectRegion].parentID == *sphRegion)
		{
			// Box is inside sphere. If box subvolume in on outer box
			// boundary then the subvolumes might be neighbors
			// SPECIAL CASE: If box is meso, then we need to count one
			// neighbor for each face that is along outer boundary, in order
			// to properly account for all of the virtual subvolumes. The number
			// of faces will be stored to numFaceSph and their directions to dirArray
			return bSubFaceRegion(regionArray, *rectRegion, *sphRegion, curSubBound,
				subHalfSize[*rectRegion], subCoorInd[*rectSub],
				boundAdjError, numFaceSph, dirArray);				
		} else if (regionArray[*sphRegion].parentID == *rectRegion)
		{
			curNeighBound[0] = regionArray[*sphRegion].boundary[0];
			curNeighBound[1] = regionArray[*sphRegion].boundary[1];
			curNeighBound[2] = regionArray[*sphRegion].boundary[2];
			curNeighBound[3] = regionArray[*sphRegion].boundary[3];
			// Sphere is inside box. Need to check intersection of box
			// subvolume with sphere surface
			if(bBoundaryIntersect(RECTANGULAR_BOX, curSubBound,
				SPHERE, curNeighBound, regionArray[*rectRegion].actualSubSize))
			{
				return true;
			}
		} else
		{
			// Something went wrong ...
			fprintf(stderr, "ERROR: Regions %u and %u were determined to be neighbors but this was only valid if one was a child region of the other.\n", curRegion, neighRegion);
			exit(EXIT_FAILURE);
		}
	}
	
	// Subvolumes are not neighbors
	return false;
}

// Calculate cartesian coordinates of rectangular subvolume
void findSubvolCoor(double subBound[6],
	struct region regionSingle,
	uint32_t subCoorInd[3])
{
	subBound[0] = regionSingle.spec.xAnch +
		regionSingle.actualSubSize*subCoorInd[0];
	subBound[1] = subBound[0] + regionSingle.actualSubSize;
	subBound[2] = regionSingle.spec.yAnch +
		regionSingle.actualSubSize*subCoorInd[1];
	subBound[3] = subBound[2] + regionSingle.actualSubSize;
	subBound[4] = regionSingle.spec.zAnch +
		regionSingle.actualSubSize*subCoorInd[2];
	subBound[5] = subBound[4] + regionSingle.actualSubSize;
}