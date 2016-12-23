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
 * Last revised for AcCoRD v1.1 (2016-12-24)
 *
 * Revision history:
 *
 * Revision v1.1 (2016-12-24)
 * - added direction of subvolume neighbors as a standalone 2D array in order
 * to implement fluid flow in the mesoscopic regime
 *
 * Revision v1.0 (2016-10-31)
 * - moved mesoscopic structure fields from subvolume struct to meso subvolume struct
 *
 * Revision v0.7.0.1 (public beta, 2016-08-30)
 * - fixed bug where a molecule with diffusion rate 0 would have an invalid reaction
 * propensity at hybrid interface
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - improved propensity calculation for molecules to leave mesoscopic subvolume and enter
 * virtual microscopic neighbor. From Flegg et al., "Analysis of the two-regime method on
 * square meshes", SIAM Journal of Scientific Computing, vol. 36, no. 3, pp. 561--588, 2014
 *
 * Revision v0.5 (2016-04-15)
 * - corrected memory allocation for subvolume helper arrays
 * - added surface subvolumes
 * - preliminary implementation of interaction between mesoscopic subvolumes and surfaces
 *
 * Revision v0.4.1
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
// A call to this function should have a corresponding call to delete_subvol_array2D,
// usually AFTER build_subvol_array2D has been called
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
	uint32_t (*** subIDSize)[2],
	unsigned short *** subNeighDir,
	const short NUM_REGIONS,
	struct region regionArray[])
{	
	short i;
	uint32_t j,k;
	uint32_t length[3];
	*subCoorInd = malloc(numSub*sizeof(uint32_t [3]));
	*subIDSize = malloc(NUM_REGIONS*sizeof(*subIDSize));
	*subID = malloc(NUM_REGIONS*sizeof(*subID));
	*subNeighDir = malloc(numSub*sizeof(unsigned short *));
	
	if(*subCoorInd == NULL || *subID == NULL ||
		*subIDSize == NULL || *subNeighDir == NULL){
		fprintf(stderr, "ERROR: Memory allocation for temporary subvolume information.\n");
		exit(EXIT_FAILURE);
	}
	
	// Fill in subIDSize for surface regions
	for(i = 0; i < NUM_REGIONS; i++)
	{
		if(regionArray[i].numFace > 0)
		{
			(*subIDSize)[i] = malloc(regionArray[i].numFace*sizeof(uint32_t [2]));
			if((*subIDSize)[i] == NULL){
				fprintf(stderr, "ERROR: Memory allocation for temporary subvolume information.\n");
				exit(EXIT_FAILURE);
			}
			if(regionArray[i].spec.shape == RECTANGULAR_BOX)
			{				
				(*subIDSize)[i][0][0] = regionArray[i].spec.numX;
				(*subIDSize)[i][0][1] = regionArray[i].spec.numY;
				(*subIDSize)[i][1][0] = regionArray[i].spec.numX;
				(*subIDSize)[i][1][1] = regionArray[i].spec.numY;
				(*subIDSize)[i][2][0] = regionArray[i].spec.numX;
				(*subIDSize)[i][2][1] = regionArray[i].spec.numZ;
				(*subIDSize)[i][3][0] = regionArray[i].spec.numX;
				(*subIDSize)[i][3][1] = regionArray[i].spec.numZ;
				(*subIDSize)[i][4][0] = regionArray[i].spec.numY;
				(*subIDSize)[i][4][1] = regionArray[i].spec.numZ;
				(*subIDSize)[i][5][0] = regionArray[i].spec.numY;
				(*subIDSize)[i][5][1] = regionArray[i].spec.numZ;
			} else if (regionArray[i].spec.type == REGION_SURFACE_3D
				|| regionArray[i].spec.type == REGION_NORMAL)
			{ // Region is 2D surface to a 3D shape
				switch(regionArray[i].plane)
				{
					case PLANE_3D: // Spherical surface; will only be 1 subvolume anyway
					case PLANE_XY:
						(*subIDSize)[i][0][0] = regionArray[i].spec.numX;
						(*subIDSize)[i][0][1] = regionArray[i].spec.numY;
						break;
					case PLANE_XZ:
						(*subIDSize)[i][0][0] = regionArray[i].spec.numX;
						(*subIDSize)[i][0][1] = regionArray[i].spec.numZ;
						break;
					case PLANE_YZ:
						(*subIDSize)[i][0][0] = regionArray[i].spec.numY;
						(*subIDSize)[i][0][1] = regionArray[i].spec.numZ;
						break;
				}
			} else
			{ // Region must be a 2D surface
				if(regionArray[i].spec.shape == RECTANGLE)
				{
					for(j = 0; j < 4; j++)
					{
						(*subIDSize)[i][j][1] = 1;
					}
					switch(regionArray[i].plane)
					{
						case PLANE_XY:
							(*subIDSize)[i][0][0] = regionArray[i].spec.numX;
							(*subIDSize)[i][1][0] = regionArray[i].spec.numX;
							(*subIDSize)[i][2][0] = regionArray[i].spec.numY;
							(*subIDSize)[i][3][0] = regionArray[i].spec.numY;
							break;
						case PLANE_XZ:
							(*subIDSize)[i][0][0] = regionArray[i].spec.numX;
							(*subIDSize)[i][1][0] = regionArray[i].spec.numX;
							(*subIDSize)[i][2][0] = regionArray[i].spec.numZ;
							(*subIDSize)[i][3][0] = regionArray[i].spec.numZ;
							break;
						case PLANE_YZ:
							(*subIDSize)[i][0][0] = regionArray[i].spec.numY;
							(*subIDSize)[i][1][0] = regionArray[i].spec.numY;
							(*subIDSize)[i][2][0] = regionArray[i].spec.numZ;
							(*subIDSize)[i][3][0] = regionArray[i].spec.numZ;
							break;
					}
				} else
				{
					fprintf(stderr, "ERROR: Region %u (Label: \"%s\") has an invalid type/shape configuration, which prevents valid subvolume creation.\n", i, regionArray[i].spec.label);
					exit(EXIT_FAILURE);
				}
			}
		}
	}
	
	for(i = 0; i < NUM_REGIONS; i++)
	{
		if(regionArray[i].numFace > 0)
		{ // Region is a surface. First dimension is the number of faces
			length[0] = regionArray[i].numFace;
			if(regionArray[i].numFace == 1)
			{
				length[1] = (*subIDSize)[i][0][0];
				length[2] = (*subIDSize)[i][0][1];
			}
		}
		else
		{ // Region is not a surface. Dimensions are default
			length[0] = regionArray[i].spec.numX;
			length[1] = regionArray[i].spec.numY;
			length[2] = regionArray[i].spec.numZ;
		}
		
		(*subID)[i] = malloc(length[0]*sizeof((*subID)[i]));
		if((*subID)[i] == NULL){
			fprintf(stderr, "ERROR: Memory allocation for temporary subvolume information.\n");
			exit(EXIT_FAILURE);
		}
		
		for(j = 0; j < length[0]; j++)
		{
			if(regionArray[i].numFace > 1)
			{ // If a surface has more than 1 face, then these 2 dimensions vary with each face
				length[1] = (*subIDSize)[i][j][0];
				length[2] = (*subIDSize)[i][j][1];
			}		
			
			(*subID)[i][j] = malloc(length[1]*sizeof((*subID)[i][j]));
			if((*subID)[i][j] == NULL){
				fprintf(stderr, "ERROR: Memory allocation for temporary subvolume information.\n");
				exit(EXIT_FAILURE);
			}
			for(k = 0; k < length[1]; k++)
			{
				(*subID)[i][j][k] = malloc(length[2]*sizeof((*subID)[i][j][k]));
				if((*subID)[i][j][k] == NULL){
					fprintf(stderr, "ERROR: Memory allocation for temporary subvolume information.\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}
}

// Free Memory of Subvolume Structure
// This function is needed to free the memory allocated to each
// array member of each structure element before releasing the memory
// allocated to the structure array itself. It should be called after BOTH
// allocate_subvol_array3D AND build_subvol_array3D, but is written to avoid
// any attempts to free memory that was not previously allocated
void deleteSubvolArray(const uint32_t numSub,
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
		if(!regionArray[subvolArray[curSub].regionID].spec.bMicro)
		{
			if(subvolArray[curSub].neighID != NULL)
				free(subvolArray[curSub].neighID);
		}
	}
	
	free(subvolArray);
}

// Free memory allocated to helper initialization arrays
void deleteSubvolHelper(uint32_t subCoorInd[][3],
	uint32_t **** subID,
	uint32_t (** subIDSize)[2],
	unsigned short ** subNeighDir,
	const short NUM_REGIONS,
	struct region regionArray[],
	const uint32_t numSub)
{	
	short i;
	uint32_t x,y;
	uint32_t length[2];
	
	if(subCoorInd != NULL) free(subCoorInd);
	
	// Free memory used to track subvolume neighbors
	if(subNeighDir != NULL)
	{
		for(x = 0; x < numSub; x++)
		{
			if(subNeighDir[x] != NULL)
				free(subNeighDir[x]);
		}
		free(subNeighDir);
	}
	
	if(subID != NULL && subIDSize != NULL)
	{
		for(i = 0; i < NUM_REGIONS; i++)
		{
			if(regionArray[i].numFace > 0)
			{ // Region is a surface. First dimension is the number of faces
				length[0] = regionArray[i].numFace;
				if(regionArray[i].numFace == 1)
				{
					length[1] = subIDSize[i][0][0];
				}
			}
			else
			{ // Region is not a surface. Dimensions are default
				length[0] = regionArray[i].spec.numX;
				length[1] = regionArray[i].spec.numY;
			}
			
			if(subID[i] != NULL)
			{
				for(x = 0; x < length[0]; x++)
				{
					if(regionArray[i].numFace > 1 && subIDSize[i] != NULL)
					{ // If a surface has more than 1 face, then 2 dimensions vary with each face
						length[1] = subIDSize[i][x][0];
					}
					
					if(subID[i][x] != NULL)
					{
						for(y = 0; y < length[1]; y++)
						{
							if(subID[i][x][y] != NULL) free(subID[i][x][y]);
						}
						free(subID[i][x]);
					}
				}
				free(subID[i]);
				if(regionArray[i].numFace > 0 && subIDSize[i] != NULL) free(subIDSize[i]);
			}
		}
		free(subID);
		free(subIDSize);
	}
}

// Construct the array of structures with details of each subvolume
void buildSubvolArray(const uint32_t numSub,
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
	uint32_t **** subID,
	uint32_t (** subIDSize)[2],
	unsigned short ** subNeighDir)
{	
	short int i,j, curRegion, neighRegion, sphRegion, rectRegion; // Current Region
	unsigned short curMolType;
	uint32_t cur1, cur2, cur3; // Coordinates of current subvolume within current region
	uint32_t length[3]; // Sizes of dimensions used for cur1, cur2, cur3
	uint32_t curCoorInd[3];
	bool bRealSub, bHaveCoor;
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
	unsigned short dirArray[6]; // Array needed when a cubic subvolume neighbors a sphere
	double boundAdjError = SUBVOL_BASE_SIZE * SUB_ADJ_RESOLUTION;
	
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
		
		// Determine lengths for cur1, cur2, cur3 based on number of faces
		if(regionArray[i].numFace > 0)
		{ // Region is a surface. First dimension is the number of faces
			length[0] = regionArray[i].numFace;
			if(regionArray[i].numFace == 1)
			{
				length[1] = subIDSize[i][0][0];
				length[2] = subIDSize[i][0][1];
			}
		}
		else
		{ // Region is not a surface. Dimensions are default
			length[0] = regionArray[i].spec.numX;
			length[1] = regionArray[i].spec.numY;
			length[2] = regionArray[i].spec.numZ;
		}
		
		for(cur1 = 0; cur1 < length[0]; cur1++)
		{
			if(regionArray[i].numFace > 1)
			{ // If a surface has more than 1 face, then these 2 dimensions vary with each face
				length[1] = subIDSize[i][cur1][0];
				length[2] = subIDSize[i][cur1][1];
			}
			
			for(cur2 = 0; cur2 < length[1]; cur2++)
			{
				for(cur3 = 0; cur3 < length[2]; cur3++)				
				{
					bRealSub = true;
					
					// Assign tentative subvolume coordinates
					curCoorInd[0] = cur1;
					curCoorInd[1] = cur2;
					curCoorInd[2] = cur3;
					
					// Check that current location is actually part of region
					bHaveCoor = false;
					
					// Check that current location is not in nested region
					if(regionArray[i].numChildren > 0)
					{						
						for(j = 0; j < regionArray[i].numChildren; j++)
						{
							// Is it possible for child to block subvolumes of parent?
							if((regionArray[i].spec.type != REGION_NORMAL
								&& regionArray[regionArray[i].childrenID[j]].spec.type
								== REGION_NORMAL))
							{
								continue; // This child cannot interfere with the
										  // subvolumes of the parent
							}
														
							if(regionArray[i].numFace > 0)
							{
								if(!bHaveCoor)
								{
									// Need actual subvolume coordinates in order to
									// compare with child regions
									findSubvolCoor(curSubBound, regionArray[i], curCoorInd);
									bHaveCoor = true;
								}
								
								if(bBoundarySurround(regionArray[i].subShape, curSubBound,
									regionArray[regionArray[i].childrenID[j]].spec.shape,
									regionArray[regionArray[i].childrenID[j]].boundary, 0.))
								{
									// subvolume is entirely within child
									bRealSub = false;
									subID[i][cur1][cur2][cur3] = UINT32_MAX;
									break;
								}
							} else if(regionArray[regionArray[i].childrenID[j]].spec.shape == SPHERE)
							{
								// Need real coordinates of subvolume
								if(!bHaveCoor)
								{
									findSubvolCoor(curSubBound, regionArray[i], curCoorInd);
									bHaveCoor = true;
								}
								if(bBoundarySurround(regionArray[i].subShape, curSubBound,
									SPHERE, regionArray[regionArray[i].childrenID[j]].boundary, 0.))
								{
									// subvolume is entirely within spherical child
									bRealSub = false;
									subID[i][cur1][cur2][cur3] = UINT32_MAX;
									break;
								}
							} else if(cur1 >= regionArray[i].childrenCoor[j][0]
								&& cur1 <= regionArray[i].childrenCoor[j][1]
								&& cur2 >= regionArray[i].childrenCoor[j][2]
								&& cur2 <= regionArray[i].childrenCoor[j][3]
								&& cur3 >= regionArray[i].childrenCoor[j][4]
								&& cur3 <= regionArray[i].childrenCoor[j][5])
							{
								bRealSub = false;
								subID[i][cur1][cur2][cur3] = UINT32_MAX;
								break;
							}
						}
					}
						
					if(!bRealSub)
					{	// We are currently inside a nested region.
						// Skip subvolume creation
						continue;
					}
					
					subCoorInd[curID][0] = curCoorInd[0];
					subCoorInd[curID][1] = curCoorInd[1];
					subCoorInd[curID][2] = curCoorInd[2];
					
					subID[i][cur1][cur2][cur3] = curID;
					
					subvolArray[curID].regionID = i;
					
					// Assign memory to array parameters					
					if(!regionArray[i].spec.bMicro)
					{ // Parameter only needed for mesoscopic regions
						subvolArray[curID].mesoID = curMesoID;
						curMesoID++;
					} else{
						subvolArray[curID].mesoID = UINT32_MAX;				
					}
					
					// Make sure remaining parameters are 0
					subvolArray[curID].num_neigh = 0;
					subvolArray[curID].bBoundary = false;
					
					// Find if an interior subvolume borders a nested region
					if(regionArray[i].numFace > 0)
					{ // All 2D subvolumes are boundary subvolumes
						subvolArray[curID].bBoundary = true;
					} else
					{
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
							} else if(((cur1 == regionArray[i].childrenCoor[j][0]-1 ||
								cur1 == regionArray[i].childrenCoor[j][1]+1)
								&& cur2 >= regionArray[i].childrenCoor[j][2]
								&& cur2 <= regionArray[i].childrenCoor[j][3]
								&& cur3 >= regionArray[i].childrenCoor[j][4]
								&& cur3 <= regionArray[i].childrenCoor[j][5])
								||
								((cur2 == regionArray[i].childrenCoor[j][2]-1 ||
								cur2 == regionArray[i].childrenCoor[j][3]+1)
								&& cur1 >= regionArray[i].childrenCoor[j][0]
								&& cur1 <= regionArray[i].childrenCoor[j][1]
								&& cur3 >= regionArray[i].childrenCoor[j][4]
								&& cur3 <= regionArray[i].childrenCoor[j][5])
								||
								((cur3 == regionArray[i].childrenCoor[j][4]-1 ||
								cur3 == regionArray[i].childrenCoor[j][5]+1)
								&& cur1 >= regionArray[i].childrenCoor[j][0]
								&& cur1 <= regionArray[i].childrenCoor[j][1]
								&& cur2 >= regionArray[i].childrenCoor[j][2]
								&& cur2 <= regionArray[i].childrenCoor[j][3]))
							{
								subvolArray[curID].bBoundary = true;
								continue; // No need to keep comparing with other children
							}
						}
					}
					curID++;
				}
			}
		}
	}
		
	// Identify subvolumes on boundary of normal 3D regions and count number
	// of subvolume neighbors in same region
	for(i = 0,curID=0; i < NUM_REGIONS; i++)
	{
		if(regionArray[i].spec.shape == SPHERE)
		{
			// Already know one subvolume in sphere is valid
			curID++;
			continue;
		}
		
		// Determine lengths for cur1, cur2, cur3 based on number of faces
		if(regionArray[i].numFace > 0)
		{ // Region is a surface. First dimension is the number of faces
			length[0] = regionArray[i].numFace;
			if(regionArray[i].numFace == 1)
			{
				length[1] = subIDSize[i][0][0];
				length[2] = subIDSize[i][0][1];
			}
		}
		else
		{ // Region is not a surface. Dimensions are default
			length[0] = regionArray[i].spec.numX;
			length[1] = regionArray[i].spec.numY;
			length[2] = regionArray[i].spec.numZ;
		}
		
		for(cur1 = 0; cur1 < length[0]; cur1++)
		{
			if(regionArray[i].numFace > 1)
			{ // If a surface has more than 1 face, then these 2 dimensions vary with each face
				length[1] = subIDSize[i][cur1][0];
				length[2] = subIDSize[i][cur1][1];
			}
			
			for(cur2 = 0; cur2 < length[1]; cur2++)
			{
				for(cur3 = 0; cur3 < length[2]; cur3++)
				{
					if(subID[i][cur1][cur2][cur3] < UINT32_MAX)
					{	// Current subvolume is valid
					
						if(regionArray[i].numFace == 0)
						{ // Neighbors along cur1 dimension only possible for normal 3D regions
							if(subCoorInd[curID][0] > 0)
							{ // Neighbor is 1 "x" index "down"
								if(subID[i][cur1-1][cur2][cur3] < UINT32_MAX)
									subvolArray[curID].num_neigh++;
							} else subvolArray[curID].bBoundary = true;
							if(subCoorInd[curID][0] < length[0]-1)
							{ // Neighbor is 1 "x" index "up"
								if(subID[i][cur1+1][cur2][cur3] < UINT32_MAX)
									subvolArray[curID].num_neigh++;
							} else subvolArray[curID].bBoundary = true;
						}
						if(subCoorInd[curID][1] > 0)
						{ // Neighbor is 1 "y" index "down"
							if(subID[i][cur1][cur2-1][cur3] < UINT32_MAX)
								subvolArray[curID].num_neigh++;
						} else subvolArray[curID].bBoundary = true;
						if(subCoorInd[curID][1] < length[1]-1)
						{ // Neighbor is 1 "y" index "up"
							if(subID[i][cur1][cur2+1][cur3] < UINT32_MAX)
								subvolArray[curID].num_neigh++;
						} else subvolArray[curID].bBoundary = true;
						if(subCoorInd[curID][2] > 0)
						{ // Neighbor is 1 "z" index "down"
							if(subID[i][cur1][cur2][cur3-1] < UINT32_MAX)
								subvolArray[curID].num_neigh++;
						} else subvolArray[curID].bBoundary = true;
						if(subCoorInd[curID][2] < length[2]-1)
						{ // Neighbor is 1 "z" index "up"
							if(subID[i][cur1][cur2][cur3+1] < UINT32_MAX)
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
				if(checkSubvolNeigh(regionArray, NUM_REGIONS, curRegion, neighRegion,
					&sphRegion,	&rectRegion, curID, curNeighID, &sphSub,
					&rectSub, numSub, subCoorInd, boundAdjError, &adjDirection,
					curSubBound, curNeighBound, &numFaceSph, dirArray))
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
		if(regionArray[subvolArray[curID].regionID].spec.bMicro)
		{
			subNeighDir[curID] = NULL;
			continue; // Don't need to store neighbor IDs for microscopic subvolumes
		}
		
		subvolArray[curID].neighID =
			malloc(subvolArray[curID].num_neigh*sizeof(uint32_t));
		subNeighDir[curID] =
			malloc(subvolArray[curID].num_neigh*sizeof(unsigned short));
		if(subvolArray[curID].neighID == NULL ||
			subNeighDir[curID] == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation for neighbors of subvolume %" PRIu32 ".\n", curID);
			exit(EXIT_FAILURE);
		}
	}
		
	// Find and store ID and direction of each subvolume's neighbors
	for(curID = 0; curID < numSub; curID++)
	{
		// Calculate neighbors in the same region
		curRegion = subvolArray[curID].regionID;
		
		if(regionArray[curRegion].numFace > 0)
		{
			length[0] = regionArray[curRegion].numFace;
			length[1] = subIDSize[curRegion][subCoorInd[curID][0]][0];
			length[2] = subIDSize[curRegion][subCoorInd[curID][0]][1];
		} else
		{
			length[0] = regionArray[curRegion].spec.numX;
			length[1] = regionArray[curRegion].spec.numY;
			length[2] = regionArray[curRegion].spec.numZ;
		}
				
		if(!regionArray[curRegion].spec.bMicro)
		{
			if(regionArray[curRegion].numFace == 0)
			{ // Neighbors along first dimension only possible for normal 3D regions
				if(subCoorInd[curID][0] > 0)
				{ // Neighbor is 1 "x" index "down"
					if(subID[curRegion][subCoorInd[curID][0]-1][subCoorInd[curID][1]][subCoorInd[curID][2]] < UINT32_MAX)
					{
						subNeighDir[curID][curSubNeigh[curID]] = LEFT;
						subvolArray[curID].neighID[curSubNeigh[curID]++] =
							subID[curRegion][subCoorInd[curID][0]-1][subCoorInd[curID][1]][subCoorInd[curID][2]];
					}
				}
				if(subCoorInd[curID][0] < length[0]-1)
				{ // Neighbor is 1 "x" index "up"
					if(subID[curRegion][subCoorInd[curID][0]+1][subCoorInd[curID][1]][subCoorInd[curID][2]] < UINT32_MAX)
					{
						subNeighDir[curID][curSubNeigh[curID]] = RIGHT;
						subvolArray[curID].neighID[curSubNeigh[curID]++] =
							subID[curRegion][subCoorInd[curID][0]+1][subCoorInd[curID][1]][subCoorInd[curID][2]];
					}
				}
			}
			if(subCoorInd[curID][1] > 0)
			{ // Neighbor is 1 "y" index "down"
				if(subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]-1][subCoorInd[curID][2]] < UINT32_MAX)
				{
					subNeighDir[curID][curSubNeigh[curID]] = DOWN;
					subvolArray[curID].neighID[curSubNeigh[curID]++] =
						subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]-1][subCoorInd[curID][2]];
				}
			}
			if(subCoorInd[curID][1] < length[1]-1)
			{ // Neighbor is 1 "y" index "up"
				if(subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]+1][subCoorInd[curID][2]] < UINT32_MAX)
				{
					subNeighDir[curID][curSubNeigh[curID]] = UP;
					subvolArray[curID].neighID[curSubNeigh[curID]++] =
						subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]+1][subCoorInd[curID][2]];
				}
			}
			if(subCoorInd[curID][2] > 0)
			{ // Neighbor is 1 "z" index "down"
				if(subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]][subCoorInd[curID][2]-1] < UINT32_MAX)
				{
					subNeighDir[curID][curSubNeigh[curID]] = IN;
					subvolArray[curID].neighID[curSubNeigh[curID]++] =
						subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]][subCoorInd[curID][2]-1];
				}
			}
			if(subCoorInd[curID][2] < length[2]-1)
			{ // Neighbor is 1 "z" index "up"
				if(subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]][subCoorInd[curID][2]+1] < UINT32_MAX)
				{
					subNeighDir[curID][curSubNeigh[curID]] = OUT;
					subvolArray[curID].neighID[curSubNeigh[curID]++] =
						subID[curRegion][subCoorInd[curID][0]][subCoorInd[curID][1]][subCoorInd[curID][2]+1];
				}
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
			if(checkSubvolNeigh(regionArray, NUM_REGIONS, curRegion, neighRegion, &sphRegion,
				&rectRegion, curID, curNeighID, &sphSub, &rectSub,
				numSub, subCoorInd, boundAdjError, &adjDirection,
				curSubBound, curNeighBound, &numFaceSph, dirArray))
			{
				// Subvolumes are neighbors. Record IDs
				if (numFaceSph > 0)
				{
					// One subvolume is in a box while the other is in a sphere.
					// The subvolumes can be neighbors along multiple faces
					if (!regionArray[rectRegion].spec.bMicro)
						for(i = 0; i < numFaceSph; i++)
						{
							subNeighDir[rectSub][curSubNeigh[rectSub]] = dirArray[i];
							subvolArray[rectSub].neighID[curSubNeigh[rectSub]++] = sphSub;
						}
					numFaceSph = 0; // Reset value
				} else
				{
					if(!regionArray[curRegion].spec.bMicro)
					{
						subNeighDir[curID][curSubNeigh[curID]] = adjDirection;
						subvolArray[curID].neighID[curSubNeigh[curID]++] = curNeighID;
					}
					if(!regionArray[neighRegion].spec.bMicro)
					{
						switch(adjDirection)
						{
							case LEFT:
								subNeighDir[curNeighID][curSubNeigh[curNeighID]] = RIGHT;
									break;
							case RIGHT:
								subNeighDir[curNeighID][curSubNeigh[curNeighID]] = LEFT;
									break;
							case UP:
								subNeighDir[curNeighID][curSubNeigh[curNeighID]] = DOWN;
									break;
							case DOWN:
								subNeighDir[curNeighID][curSubNeigh[curNeighID]] = UP;
									break;
							case IN:
								subNeighDir[curNeighID][curSubNeigh[curNeighID]] = OUT;
									break;
							case OUT:
								subNeighDir[curNeighID][curSubNeigh[curNeighID]] = IN;
									break;
						}
						subvolArray[curNeighID].neighID[curSubNeigh[curNeighID]++] = curID;
					}
				}
			}
		}
	}
		
	if(NUM_REGIONS > 1)
	{ // Only need to enter if there is more than one region

		// Initialize region knowledge of the subvolumes that are adjacent to it
		initializeRegionSubNeighbor(regionArray, subvol_spec,
			NUM_REGIONS, NUM_MOL_TYPES, SUBVOL_BASE_SIZE,
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
	const short NUM_REGIONS,
	short int curRegion,
	short int neighRegion,
	short int * sphRegion,
	short int * rectRegion,
	uint32_t curID,
	uint32_t curNeighID,
	uint32_t * sphSub,
	uint32_t * rectSub,
	uint32_t numSub,
	uint32_t subCoorInd[numSub][3],
	double boundAdjError,
	unsigned short * adjDirection,
	double curSubBound[6],
	double curNeighBound[6],
	unsigned short * numFaceSph,
	unsigned short dirArray[6])
{
	//unsigned short dirArray[6];
	unsigned short surfaceRegion;

	if ((regionArray[curRegion].spec.shape == RECTANGULAR_BOX
		|| regionArray[curRegion].spec.shape == RECTANGLE)
		&& (regionArray[neighRegion].spec.shape == RECTANGULAR_BOX
		|| regionArray[neighRegion].spec.shape == RECTANGLE))
	{
		findSubvolCoor(curSubBound, regionArray[curRegion], subCoorInd[curID]);
		findSubvolCoor(curNeighBound, regionArray[neighRegion],
			subCoorInd[curNeighID]);
			
		// Are the two subvolumes close enough to be neighbours?
		if(bBoundaryAdjacent(regionArray[curRegion].subShape, curSubBound,
			regionArray[neighRegion].subShape,
			curNeighBound, boundAdjError, adjDirection))
		{ // Subvolumes are geometric neighbours			
			if(regionArray[curRegion].spec.type == REGION_SURFACE_3D
				|| regionArray[curRegion].spec.type == REGION_SURFACE_2D
				|| regionArray[neighRegion].spec.type == REGION_SURFACE_3D
				|| regionArray[neighRegion].spec.type == REGION_SURFACE_2D)
			{ // At least one of the region is a surface. Neighboring is valid
				return true;
			} else
			{ // Both regions are normal.
				// Check whether there is a surface region between them				
				return bSurfaceBetweenBoundaries(regionArray, NUM_REGIONS, curRegion,
					neighRegion, curSubBound, curNeighBound, &surfaceRegion);
			}
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
			return bSubFaceRegion(regionArray, *sphRegion, 
				curSubBound, boundAdjError, numFaceSph, dirArray);				
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
	} else
	{
		// Something went wrong ...
		fprintf(stderr, "ERROR: We cannot determine whether the subvolumes of neighboring regions %u and %u (labels: \"%s\" and \"%s\") are also neighbors.\n",
			curRegion, neighRegion, regionArray[curRegion].spec.label,
			regionArray[neighRegion].spec.label);
		exit(EXIT_FAILURE);
	}
	
	// Subvolumes are not neighbors
	return false;
}

// Calculate cartesian coordinates of rectangular subvolume
void findSubvolCoor(double subBound[6],
	const struct region regionSingle,
	uint32_t subCoorInd[3])
{
	if(regionSingle.numFace > 0)
	{ // Region is a surface. First dimension is the number of faces
		if(regionSingle.spec.shape == RECTANGULAR_BOX)
		{ // There are 6 faces
			switch(subCoorInd[0])
			{
				case 0: // Lower XY plane
					subBound[0] = regionSingle.boundary[0] +
						regionSingle.actualSubSize*subCoorInd[1];
					subBound[1] = subBound[0] +	regionSingle.actualSubSize;
					subBound[2] = regionSingle.boundary[2] +
						regionSingle.actualSubSize*subCoorInd[2];
					subBound[3] = subBound[2] +	regionSingle.actualSubSize;
					subBound[4] = regionSingle.boundary[4];
					subBound[5] = regionSingle.boundary[4];
					break;
				case 1: // Upper XY plane
					subBound[0] = regionSingle.boundary[0] +
						regionSingle.actualSubSize*subCoorInd[1];
					subBound[1] = subBound[0] +	regionSingle.actualSubSize;
					subBound[2] = regionSingle.boundary[2] +
						regionSingle.actualSubSize*subCoorInd[2];
					subBound[3] = subBound[2] +	regionSingle.actualSubSize;
					subBound[4] = regionSingle.boundary[5];
					subBound[5] = regionSingle.boundary[5];
					break;
				case 2: // Lower XZ plane
					subBound[0] = regionSingle.boundary[0] +
						regionSingle.actualSubSize*subCoorInd[1];
					subBound[1] = subBound[0] +	regionSingle.actualSubSize;
					subBound[2] = regionSingle.boundary[2];
					subBound[3] = regionSingle.boundary[2];
					subBound[4] = regionSingle.boundary[4] +
						regionSingle.actualSubSize*subCoorInd[2];
					subBound[5] = subBound[4] +	regionSingle.actualSubSize;
					break;
				case 3: // Upper XZ plane
					subBound[0] = regionSingle.boundary[0] +
						regionSingle.actualSubSize*subCoorInd[1];
					subBound[1] = subBound[0] +	regionSingle.actualSubSize;
					subBound[2] = regionSingle.boundary[3];
					subBound[3] = regionSingle.boundary[3];
					subBound[4] = regionSingle.boundary[4] +
						regionSingle.actualSubSize*subCoorInd[2];
					subBound[5] = subBound[4] +	regionSingle.actualSubSize;
					break;
				case 4: // Lower YZ plane
					subBound[0] = regionSingle.boundary[0];
					subBound[1] = regionSingle.boundary[0];
					subBound[2] = regionSingle.boundary[2] +
						regionSingle.actualSubSize*subCoorInd[1];
					subBound[3] = subBound[2] +	regionSingle.actualSubSize;
					subBound[4] = regionSingle.boundary[4] +
						regionSingle.actualSubSize*subCoorInd[2];
					subBound[5] = subBound[4] +	regionSingle.actualSubSize;
					break;
				case 5: // Upper YZ plane
					subBound[0] = regionSingle.boundary[1];
					subBound[1] = regionSingle.boundary[1];
					subBound[2] = regionSingle.boundary[2] +
						regionSingle.actualSubSize*subCoorInd[1];
					subBound[3] = subBound[2] +	regionSingle.actualSubSize;
					subBound[4] = regionSingle.boundary[4] +
						regionSingle.actualSubSize*subCoorInd[2];
					subBound[5] = subBound[4] +	regionSingle.actualSubSize;
					break;
				default:
					// Something went wrong here
					fprintf(stderr, "ERROR: Face %" PRIu32 " is invalid for a rectangular box.\nWe should not have gotten here!\n", subCoorInd[0]);
					exit(EXIT_FAILURE);
			}
			
		} else if(regionSingle.spec.type == REGION_SURFACE_3D
			|| regionSingle.spec.type == REGION_NORMAL)
		{ // Region is 2D surface to a 3D shape
			switch(regionSingle.plane)
			{
				case PLANE_XY:
					subBound[0] = regionSingle.boundary[0] +
						regionSingle.actualSubSize*subCoorInd[1];
					subBound[1] = subBound[0] +	regionSingle.actualSubSize;
					subBound[2] = regionSingle.boundary[2] +
						regionSingle.actualSubSize*subCoorInd[2];
					subBound[3] = subBound[2] +	regionSingle.actualSubSize;
					subBound[4] = regionSingle.boundary[4];
					subBound[5] = regionSingle.boundary[4];
					break;
				case PLANE_XZ:
					subBound[0] = regionSingle.boundary[0] +
						regionSingle.actualSubSize*subCoorInd[1];
					subBound[1] = subBound[0] +	regionSingle.actualSubSize;
					subBound[2] = regionSingle.boundary[2];
					subBound[3] = regionSingle.boundary[2];
					subBound[4] = regionSingle.boundary[4] +
						regionSingle.actualSubSize*subCoorInd[2];
					subBound[5] = subBound[4] +	regionSingle.actualSubSize;
					break;
				case PLANE_YZ:
					subBound[0] = regionSingle.boundary[0];
					subBound[1] = regionSingle.boundary[0];
					subBound[2] = regionSingle.boundary[2] +
						regionSingle.actualSubSize*subCoorInd[1];
					subBound[3] = subBound[2] +	regionSingle.actualSubSize;
					subBound[4] = regionSingle.boundary[4] +
						regionSingle.actualSubSize*subCoorInd[2];
					subBound[5] = subBound[4] +	regionSingle.actualSubSize;
					break;
				default:
					// Something went wrong here
					fprintf(stderr, "ERROR: \"plane\" member of 2D region structure was not properly defined!.\nWe should not have gotten here!\n");
					exit(EXIT_FAILURE);
			}			
		} else if(regionSingle.spec.type == REGION_SURFACE_2D)
		{ // Region is the surface of 2D region
			switch(regionSingle.plane)
			{
				case PLANE_XY:
					switch(subCoorInd[0])
					{
						case 0: // Lower Y
							subBound[0] = regionSingle.boundary[0] +
								regionSingle.actualSubSize*subCoorInd[1];
							subBound[1] = subBound[0] +	regionSingle.actualSubSize;
							subBound[2] = regionSingle.boundary[2];
							subBound[3] = regionSingle.boundary[2];
							break;
						case 1: // Upper Y
							subBound[0] = regionSingle.boundary[0] +
								regionSingle.actualSubSize*subCoorInd[1];
							subBound[1] = subBound[0] +	regionSingle.actualSubSize;
							subBound[2] = regionSingle.boundary[3];
							subBound[3] = regionSingle.boundary[3];
							break;
						case 2: // Lower X
							subBound[0] = regionSingle.boundary[0];
							subBound[1] = regionSingle.boundary[0];
							subBound[2] = regionSingle.boundary[2] +
								regionSingle.actualSubSize*subCoorInd[1];
							subBound[3] = subBound[2] +	regionSingle.actualSubSize;
							break;
						case 3: // Upper X
							subBound[0] = regionSingle.boundary[1];
							subBound[1] = regionSingle.boundary[1];
							subBound[2] = regionSingle.boundary[2] +
								regionSingle.actualSubSize*subCoorInd[1];
							subBound[3] = subBound[2] +	regionSingle.actualSubSize;
							break;
						default:
							// Something went wrong here
							fprintf(stderr, "ERROR: Face %" PRIu32 " is invalid for a rectangle.\nWe should not have gotten here!\n", subCoorInd[0]);
							exit(EXIT_FAILURE);
					}
					subBound[4] = regionSingle.boundary[4];
					subBound[5] = regionSingle.boundary[4];
					break;
				case PLANE_XZ:
					switch(subCoorInd[0])
					{
						case 0: // Lower Z
							subBound[0] = regionSingle.boundary[0] +
								regionSingle.actualSubSize*subCoorInd[1];
							subBound[1] = subBound[0] +	regionSingle.actualSubSize;
							subBound[4] = regionSingle.boundary[4];
							subBound[5] = regionSingle.boundary[5];
							break;
						case 1: // Upper Z
							subBound[0] = regionSingle.boundary[0] +
								regionSingle.actualSubSize*subCoorInd[1];
							subBound[1] = subBound[0] +	regionSingle.actualSubSize;
							subBound[4] = regionSingle.boundary[5];
							subBound[5] = regionSingle.boundary[5];
							break;
						case 2: // Lower X
							subBound[0] = regionSingle.boundary[0];
							subBound[1] = regionSingle.boundary[0];
							subBound[4] = regionSingle.boundary[4] +
								regionSingle.actualSubSize*subCoorInd[1];
							subBound[5] = subBound[4] +	regionSingle.actualSubSize;
							break;
						case 3: // Upper X
							subBound[0] = regionSingle.boundary[1];
							subBound[1] = regionSingle.boundary[1];
							subBound[4] = regionSingle.boundary[4] +
								regionSingle.actualSubSize*subCoorInd[1];
							subBound[5] = subBound[4] +	regionSingle.actualSubSize;
							break;
						default:
							// Something went wrong here
							fprintf(stderr, "ERROR: Face %" PRIu32 " is invalid for a rectangle.\nWe should not have gotten here!\n", subCoorInd[0]);
							exit(EXIT_FAILURE);
					}
					subBound[2] = regionSingle.boundary[2];
					subBound[3] = regionSingle.boundary[2];
					break;
				case PLANE_YZ:
					switch(subCoorInd[0])
					{
						case 0: // Lower Z
							subBound[2] = regionSingle.boundary[2] +
								regionSingle.actualSubSize*subCoorInd[1];
							subBound[3] = subBound[2] +	regionSingle.actualSubSize;
							subBound[4] = regionSingle.boundary[4];
							subBound[5] = regionSingle.boundary[5];
							break;
						case 1: // Upper Z
							subBound[2] = regionSingle.boundary[2] +
								regionSingle.actualSubSize*subCoorInd[1];
							subBound[3] = subBound[2] +	regionSingle.actualSubSize;
							subBound[4] = regionSingle.boundary[5];
							subBound[5] = regionSingle.boundary[5];
							break;
						case 2: // Lower Y
							subBound[2] = regionSingle.boundary[2];
							subBound[3] = regionSingle.boundary[2];
							subBound[4] = regionSingle.boundary[4] +
								regionSingle.actualSubSize*subCoorInd[1];
							subBound[5] = subBound[4] +	regionSingle.actualSubSize;
							break;
						case 3: // Upper Y
							subBound[2] = regionSingle.boundary[3];
							subBound[3] = regionSingle.boundary[3];
							subBound[4] = regionSingle.boundary[4] +
								regionSingle.actualSubSize*subCoorInd[1];
							subBound[5] = subBound[4] +	regionSingle.actualSubSize;
							break;
						default:
							// Something went wrong here
							fprintf(stderr, "ERROR: Face %" PRIu32 " is invalid for a rectangle.\nWe should not have gotten here!\n", subCoorInd[0]);
							exit(EXIT_FAILURE);
					}
					subBound[0] = regionSingle.boundary[0];
					subBound[1] = regionSingle.boundary[0];
					break;
				default:
					// Something went wrong here
					fprintf(stderr, "ERROR: \"plane\" member of 2D region structure was not properly defined!.\nWe should not have gotten here!\n");
					exit(EXIT_FAILURE);
			}
		} else
		{
			// Something went wrong here
			fprintf(stderr, "ERROR: Subvolume coordinates could not be determined for a region with shape %u and type %u.\nWe should not have gotten here!\n",
				regionSingle.spec.shape, regionSingle.spec.type);
			exit(EXIT_FAILURE);
		}
	} else
	{ // Region is not a surface. Dimensions are default
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
}