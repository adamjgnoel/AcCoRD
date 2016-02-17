/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * region.c - 	operations for (microscopic or mesoscopic) regions in
 * 				simulation environment
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
 * - created new structure members to account for meso subvolumes that border
 * the same micro region along more than 1 face
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
#include <string.h> // for strcpy()
#include <limits.h> // For SHRT_MAX
#include <float.h> // For DBL_MAX

#include "region.h" // for "Public" declarations

#include "chem_rxn.h" // For deleting chemical reaction region members
#include "subvolume.h"

//
// "Private" Declarations
//

//
// Definitions
//

// Initialize array of region structs from specification
void initializeRegionArray(struct region regionArray[],
	const struct spec_region3D subvol_spec[],
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const double SUBVOL_BASE_SIZE,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	const unsigned short MAX_RXNS,
	const struct chem_rxn_struct * chem_rxn)
{
	short i,j, curChild; // Current region
	uint32_t curX, curY, curZ; // Coordinates of prospective subvolume
	double curSubBound[6]; // Boundary of prospective subvolume
	uint32_t subCoorInd[3]; // (curX, curY, curZ)
	uint32_t curID; // Subvolume index
	unsigned short curMolType;
	double h_i;
	short curChildArray[NUM_REGIONS];
	
	double boundAdjError = SUBVOL_BASE_SIZE * SUB_ADJ_RESOLUTION;
	
	double boundaryTemp[6]; // Used to store corrected radius for checking valid placement
							// of a round region as child or parent
	
	for(i = 0; i < NUM_REGIONS; i++)
	{
		// Copy main region parameters and outer boundary
		regionArray[i].spec = subvol_spec[i];
		if(subvol_spec[i].shape == RECTANGULAR_BOX ||
			subvol_spec[i].shape == RECTANGLE)
		{
			regionArray[i].actualSubSize = subvol_spec[i].sizeRect*SUBVOL_BASE_SIZE;
			regionArray[i].boundary[0] = subvol_spec[i].xAnch;
			regionArray[i].boundary[1] = regionArray[i].boundary[0] +
				(regionArray[i].actualSubSize*subvol_spec[i].numX);
			regionArray[i].boundary[2] = subvol_spec[i].yAnch;
			regionArray[i].boundary[3] = regionArray[i].boundary[2] +
				(regionArray[i].actualSubSize*subvol_spec[i].numY);
			regionArray[i].boundary[4] = subvol_spec[i].zAnch;
			regionArray[i].boundary[5] = regionArray[i].boundary[4] +
				(regionArray[i].actualSubSize*subvol_spec[i].numZ);
		} else
		{
			regionArray[i].actualSubSize = subvol_spec[i].radius;
			regionArray[i].boundary[0] = subvol_spec[i].xAnch;
			regionArray[i].boundary[1] = subvol_spec[i].yAnch;
			regionArray[i].boundary[2] = subvol_spec[i].zAnch;
			regionArray[i].boundary[3] = subvol_spec[i].radius;
			regionArray[i].boundary[4] = squareDBL(subvol_spec[i].radius);
			regionArray[i].boundary[5] = 0.;
		}		
			
		regionArray[i].numChildren = 0;
		regionArray[i].subResolution = SUBVOL_BASE_SIZE * SUB_ADJ_RESOLUTION;
		curChildArray[i] = 0;
	}
	
	// Determine region nesting
	for(i = 0; i < NUM_REGIONS; i++)
	{	// Assign parents and count number of children in each region
		regionArray[i].bParent = false;
		regionArray[i].parentID = SHRT_MAX;
		if(subvol_spec[i].parent && strlen(subvol_spec[i].parent) > 0)
		{	// Region has a parent. Determine index of parent
			for(j = 0; j < NUM_REGIONS; j++)
			{
				if(!strcmp(subvol_spec[i].parent, subvol_spec[j].label)
					&& i != j)
				{
					// Region j is the parent
					regionArray[i].bParent = true;
					regionArray[i].parentID = j;
					regionArray[j].numChildren++;
					break;
				}
			}
		}
	}	
	for(i = 0; i < NUM_REGIONS; i++)
	{	// Allocate memory to store children indices
		if(regionArray[i].numChildren > 0)
		{
			// Region has children. Allocate children array
			regionArray[i].childrenID = malloc(regionArray[i].numChildren*sizeof(short));
			regionArray[i].childrenCoor = malloc(regionArray[i].numChildren*sizeof(uint32_t [6]));
			if(regionArray[i].childrenID == NULL || regionArray[i].childrenCoor == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation to region %u parameters.\n", i);
				exit(EXIT_FAILURE);
			}
		}
	}	
	for(i = 0; i < NUM_REGIONS; i++)
	{	// Store indices of child regions nested in each parent
		if(regionArray[i].bParent)
		{
			j = regionArray[i].parentID;
			regionArray[j].childrenID[curChildArray[j]++] = i;
		}
	}
	for(i = 0; i < NUM_REGIONS; i++)
	{	// Store relative coordinates of child regions
		for(curChild = 0; curChild < regionArray[i].numChildren; curChild++)
		{
			j = regionArray[i].childrenID[curChild];
			
			if(regionArray[i].spec.shape == RECTANGULAR_BOX &&
				regionArray[j].spec.shape == RECTANGULAR_BOX)
			{
				// Neither parent nor child are spheres
				// Check for flush subvolumes first
				if(fabs((regionArray[j].boundary[0] - regionArray[i].boundary[0])/regionArray[i].actualSubSize - round((regionArray[j].boundary[0] - regionArray[i].boundary[0])/regionArray[i].actualSubSize)) > boundAdjError ||
					fabs((regionArray[j].boundary[1] - regionArray[i].boundary[1])/regionArray[i].actualSubSize - round((regionArray[j].boundary[1] - regionArray[i].boundary[1])/regionArray[i].actualSubSize)) > boundAdjError ||
					fabs((regionArray[j].boundary[2] - regionArray[i].boundary[2])/regionArray[i].actualSubSize - round((regionArray[j].boundary[2] - regionArray[i].boundary[2])/regionArray[i].actualSubSize)) > boundAdjError ||
					fabs((regionArray[j].boundary[3] - regionArray[i].boundary[3])/regionArray[i].actualSubSize - round((regionArray[j].boundary[3] - regionArray[i].boundary[3])/regionArray[i].actualSubSize)) > boundAdjError ||
					fabs((regionArray[j].boundary[4] - regionArray[i].boundary[4])/regionArray[i].actualSubSize - round((regionArray[j].boundary[4] - regionArray[i].boundary[4])/regionArray[i].actualSubSize)) > boundAdjError ||
					fabs((regionArray[j].boundary[5] - regionArray[i].boundary[5])/regionArray[i].actualSubSize - round((regionArray[j].boundary[5] - regionArray[i].boundary[5])/regionArray[i].actualSubSize)) > boundAdjError)
				{
					// Child is not flush with region's subvolumes
					fprintf(stderr, "ERROR: Nested region %u is not properly aligned within parent region %u.\n", j, i);
					fprintf(stderr, "Both regions are rectangular boxes and the outer boundary of the nested region must be flush with subvolumes of the parent region.\n");
					exit(EXIT_FAILURE);
				}
				
				// Confirm that child is inside of parent
				if(!bBoundarySurround(RECTANGULAR_BOX, regionArray[j].boundary,
					RECTANGULAR_BOX, regionArray[i].boundary, -boundAdjError))
				{
					// Child is sticking out of parent
					fprintf(stderr, "ERROR: Outer boundary of nested region %u is not properly surrounded by parent region %u.\n", j, i);
					fprintf(stderr, "Both regions are rectangular boxes.\n");
					exit(EXIT_FAILURE);
				}
				
				regionArray[i].childrenCoor[curChild][0] = (uint32_t)
					round((regionArray[j].boundary[0] - regionArray[i].boundary[0])/regionArray[i].actualSubSize);
				regionArray[i].childrenCoor[curChild][1] = (uint32_t)
					round((regionArray[j].boundary[1] - regionArray[i].boundary[0])/regionArray[i].actualSubSize - 1);
				regionArray[i].childrenCoor[curChild][2] = (uint32_t)
					round((regionArray[j].boundary[2] - regionArray[i].boundary[2])/regionArray[i].actualSubSize);
				regionArray[i].childrenCoor[curChild][3] = (uint32_t)
					round((regionArray[j].boundary[3] - regionArray[i].boundary[2])/regionArray[i].actualSubSize - 1);
				regionArray[i].childrenCoor[curChild][4] = (uint32_t)
					round((regionArray[j].boundary[4] - regionArray[i].boundary[4])/regionArray[i].actualSubSize);
				regionArray[i].childrenCoor[curChild][5] = (uint32_t)
					round((regionArray[j].boundary[5] - regionArray[i].boundary[4])/regionArray[i].actualSubSize - 1);
			} else if(regionArray[i].spec.shape == SPHERE &&
				regionArray[j].spec.shape == RECTANGULAR_BOX)
			{				
				if(!bBoundarySurround(RECTANGULAR_BOX, regionArray[j].boundary, SPHERE, regionArray[i].boundary, regionArray[j].actualSubSize))
				{
					// Child is sticking out of parent
					fprintf(stderr, "ERROR: Outer boundary of nested region %u is not properly surrounded by parent region %u.\n", j, i);
					fprintf(stderr, "Rectangular region %u is inside of spherical region %u and there must be a clearance between the surfaces of at least the subvolume size of region %u which is %.2em.\n", j, i, j, regionArray[j].actualSubSize);
					exit(EXIT_FAILURE);
				}
			} else if(regionArray[i].spec.shape == RECTANGULAR_BOX &&
				regionArray[j].spec.shape == SPHERE)
			{
				if(!regionArray[i].spec.bMicro)
				{
					// A mesoscopic parent cannot have a spherical child
					fprintf(stderr, "ERROR: Mesoscopic region %u has a spherical child region %u.\n", i, j);
					exit(EXIT_FAILURE);					
				}
				
				if(!bBoundarySurround(SPHERE, regionArray[j].boundary, RECTANGULAR_BOX, regionArray[i].boundary, regionArray[i].actualSubSize))
				{
					// Child is sticking out of parent
					fprintf(stderr, "ERROR: Outer boundary of nested region %u is not properly surrounded by parent region %u.\n", j, i);
					fprintf(stderr, "Spherical region %u is inside of rectangular region %u and there must be a clearance between the surfaces of at least the subvolume size of region %u which is %.2em.\n", j, i, i, regionArray[i].actualSubSize);
					exit(EXIT_FAILURE);
				}
			} else if(regionArray[i].spec.shape == SPHERE &&
				regionArray[j].spec.shape == SPHERE)
			{
				if(!bBoundarySurround(SPHERE, regionArray[j].boundary, SPHERE, regionArray[i].boundary, 0.))
				{
					// Child is sticking out of parent
					fprintf(stderr, "ERROR: Outer boundary of spherical nested region %u is not properly surrounded by spherical parent region %u.\n", j, i);
					exit(EXIT_FAILURE);
				}
			} else
			{
				// Combination of child and parent is invalid
				fprintf(stderr, "ERROR: Combination of child/parent for region %u nested in region %u is invalid.\n", j, i);
				exit(EXIT_FAILURE);
			}
		}
	}
	
	curID = 0;
	for(i = 0; i < NUM_REGIONS; i++)
	{
		// Calculate region volume now that know the validity of child regions
		regionArray[i].volume = findRegionVolume(regionArray, i, false);
		
		// Determine number of subvolumes
		regionArray[i].firstID = curID;
		if(regionArray[i].spec.shape == SPHERE)
			regionArray[i].numSub = 1UL; // Spherical regions always have one subvolume
		else
		{
			regionArray[i].numSub = (uint32_t) regionArray[i].spec.numX *
				regionArray[i].spec.numY * regionArray[i].spec.numZ;
			// Subtract subvolumes lost due to each child
			for(curChild = 0; curChild < regionArray[i].numChildren; curChild++)
			{
				j = regionArray[i].childrenID[curChild];
				if(regionArray[j].spec.shape == RECTANGULAR_BOX)
				{
					regionArray[i].numSub -=
						(regionArray[i].childrenCoor[curChild][1] -
						regionArray[i].childrenCoor[curChild][0] + 1) *
						(regionArray[i].childrenCoor[curChild][3] -
						regionArray[i].childrenCoor[curChild][2] + 1) *
						(regionArray[i].childrenCoor[curChild][5] -
						regionArray[i].childrenCoor[curChild][4] + 1);
				} else if (regionArray[j].spec.shape == SPHERE)
				{
					// Need to scan through candidate subvolumes and remove those
					// that would be entirely within the spherical child
					for(curZ = 0; curZ < subvol_spec[i].numZ; curZ++)
					{
						for(curY = 0; curY < subvol_spec[i].numY; curY++)
						{
							for(curX = 0; curX < subvol_spec[i].numX; curX++)
							{
								subCoorInd[0] = curX;
								subCoorInd[1] = curY;
								subCoorInd[2] = curZ;
								findSubvolCoor(curSubBound, regionArray[i], subCoorInd);
								if(bBoundarySurround(RECTANGULAR_BOX, curSubBound,
								SPHERE, regionArray[regionArray[i].childrenID[curChild]].boundary, 0.))
								{
									// This prospective subvolume is entirely within
									// spherical child
									regionArray[i].numSub--;
								}
							}
						}
					}
				}
			}
		}
		
		regionArray[i].isRegionNeigh = malloc(NUM_REGIONS*sizeof(bool));
		regionArray[i].regionNeighDir = malloc(NUM_REGIONS*sizeof(unsigned short));
		if(regionArray[i].isRegionNeigh == NULL || regionArray[i].regionNeighDir == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation for region %u parameters.\n", i);
			exit(EXIT_FAILURE);
		}
		
		// Initialize members used for transitions out of microscopic regions
		if(regionArray[i].spec.bMicro)
		{
			regionArray[i].numRegionNeighFace = malloc(NUM_REGIONS*sizeof(short));
			regionArray[i].boundRegionFaceCoor = 
				malloc(NUM_REGIONS * sizeof(regionArray[i].boundRegionFaceCoor));
			regionArray[i].regionNeighFaceDir = 
				malloc(NUM_REGIONS * sizeof(regionArray[i].regionNeighFaceDir));
				
			if(regionArray[i].numRegionNeighFace == NULL
				|| regionArray[i].boundRegionFaceCoor == NULL
				|| regionArray[i].regionNeighFaceDir == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation for region %u's microscopic parameters.\n", i);
				exit(EXIT_FAILURE);
			}
		}
		
		regionArray[i].numSubRegionNeigh = malloc(NUM_REGIONS*sizeof(uint32_t));
		if(regionArray[i].numSubRegionNeigh == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation for region %u parameters.\n", i);
			exit(EXIT_FAILURE);
		}
		
		for(j = 0; j < NUM_REGIONS; j++)
		{ // Initialize isRegionNeigh to false
			regionArray[i].isRegionNeigh[j] = false;
			regionArray[i].regionNeighDir[j] = UNDEFINED;
		}
		
		regionArray[i].numRegionNeigh = 0; // Initialize # of neighboring regions to 0
		
		// Calculate diffusion rates within region
		if(!regionArray[i].spec.bMicro)
		{
			regionArray[i].diffRate = malloc(NUM_MOL_TYPES*sizeof(double));
			if(regionArray[i].diffRate == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation for region %u diffusion coefficients.\n", i);
				exit(EXIT_FAILURE);
			}
			for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
			{
				h_i = regionArray[i].actualSubSize;
				regionArray[i].diffRate[curMolType] =
					DIFF_COEF[i][curMolType]/h_i/h_i;
			}
		}
		
		curID += regionArray[i].numSub;
	}
	
	// Determine what regions are touching (reduces neighbour search computation)
	findRegionTouch3D(NUM_REGIONS, subvol_spec, regionArray, SUBVOL_BASE_SIZE);
	
	// Define chemical reaction network
	initialize_region_chem_rxn3D(NUM_REGIONS, regionArray,
		NUM_MOL_TYPES, MAX_RXNS, chem_rxn);
}

// Initialize region knowledge of the subvolumes that are adjacent to it
// This function is called by build_subvol_array3D in subvolume.c
void initializeRegionSubNeighbor(struct region regionArray[],
	const struct spec_region3D subvol_spec[],
	const double subHalfSize[],
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const double SUBVOL_BASE_SIZE,
	const double boundAdjError,
	struct subvolume3D * subvolArray,
	uint32_t subCoorInd[][3])
{
	short i, j, k;
	uint32_t curID, curNeighID, neighSubID, curBoundID; // subvolume indices
	unsigned short curMolType;
	unsigned short subDir; // Direction of neighboring subvolume
	unsigned short dirArray[6]; // Direction of neighboring subvolume in case of parent/child
	unsigned short curDir, numDir; // Current neighbor direction in case of parent/child
	double curSubBound[6];
	double neighSubBound[6];
	
	// Determine numSubRegionNeigh - # of subvolumes in each region that are along boundary of each neighboring region
	for(i = 0; i < NUM_REGIONS; i++)
	{
		if(subvol_spec[i].bMicro)
			continue; // Only need to proceed if region is mesoscopic
		
		regionArray[i].neighID = malloc(NUM_REGIONS * sizeof(uint32_t *));
		regionArray[i].boundSubNumFace = malloc(NUM_REGIONS * sizeof(unsigned short *));
		regionArray[i].boundSubCoor =
			malloc(NUM_REGIONS * sizeof(regionArray[i].boundSubCoor));
		regionArray[i].boundVirtualNeighCoor =
			malloc(NUM_REGIONS * sizeof(regionArray[i].boundVirtualNeighCoor));
		regionArray[i].bNeedUpdate =
			malloc(NUM_REGIONS * sizeof(regionArray[i].bNeedUpdate));
		regionArray[i].numMolFromMicro =
			malloc(NUM_REGIONS * sizeof(regionArray[i].numMolFromMicro));
		if(regionArray[i].neighID == NULL
			|| regionArray[i].boundSubNumFace == NULL
			|| regionArray[i].boundSubCoor == NULL
			|| regionArray[i].boundVirtualNeighCoor == NULL
			|| regionArray[i].bNeedUpdate == NULL
			|| regionArray[i].numMolFromMicro == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation for region %u neighbor parameters.\n", i);
			exit(EXIT_FAILURE);
		}
		
		for(j = 0; j < NUM_REGIONS; j++)
		{	// Scan boundary subvolumes to find those bordering this region
			regionArray[i].numSubRegionNeigh[j] = 0;
			
			if(!regionArray[i].isRegionNeigh[j])
				continue; // Regions are the same or they don't border each other
			
			if(!subvol_spec[j].bMicro)
				continue; // Only need to proceed if neighbor is microscopic
			
			/* Need to pass through ith region subvolumes twice. First pass
			* is to find numSubRegionNeigh. Then, once sufficient memory
			* is allocated, a second pass is made to determine the IDs and
			* coordinates of subvolumes that touch region j
			*/
			
			for(curID = regionArray[i].firstID;
				curID < (regionArray[i].firstID + regionArray[i].numSub); curID++)
			{ // For each subvolume in current region
				if(!subvolArray[curID].bBoundary)
					continue; // Subvolume cannot border any regions
				
				findSubvolCoor(curSubBound, regionArray[i], subCoorInd[curID]);
				
				if(bSubFaceRegion(regionArray, i, j, curSubBound,
					subHalfSize[i], subCoorInd[curID], boundAdjError, &numDir, dirArray))
				{
					// This subvolume borders microscopic region j along numDir faces
					regionArray[i].numSubRegionNeigh[j]++;
				}
			}
			
			if(regionArray[i].numSubRegionNeigh[j] < 1)
				continue; // No neighbours found
						  // (failsafe; should not have reached here)
			
			regionArray[i].neighID[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(uint32_t));
			regionArray[i].boundSubNumFace[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(unsigned short));
			regionArray[i].boundSubCoor[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(double [3]));
			regionArray[i].boundVirtualNeighCoor[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(regionArray[i].boundVirtualNeighCoor[j]));
			regionArray[i].bNeedUpdate[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(bool));
			regionArray[i].numMolFromMicro[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(regionArray[i].numMolFromMicro[j]));
			if(regionArray[i].neighID[j] == NULL
				|| regionArray[i].boundSubNumFace[j] == NULL
				|| regionArray[i].boundSubCoor[j] == NULL
				|| regionArray[i].boundVirtualNeighCoor[j] == NULL
				|| regionArray[i].bNeedUpdate[j] == NULL
				|| regionArray[i].numMolFromMicro[j] == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation for region %u's neighbor parameters with region %u.\n", i, j);
				exit(EXIT_FAILURE);
			}
			
			for(curBoundID = 0; curBoundID < regionArray[i].numSubRegionNeigh[j];
				curBoundID++)
			{
				regionArray[i].bNeedUpdate[j][curBoundID] = false;
				regionArray[i].numMolFromMicro[j][curBoundID] =
					malloc(NUM_MOL_TYPES * sizeof(uint64_t));
				if(regionArray[i].numMolFromMicro[j][curBoundID] == NULL)
				{
					fprintf(stderr, "ERROR: Memory allocation for region %u's neighbor parameters with region %u.\n", i, j);
					exit(EXIT_FAILURE);
				}
				// Initialize array values
				for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
				{
					regionArray[i].numMolFromMicro[j][curBoundID][curMolType] = 0ULL;
				}
			}
			
			curBoundID = 0;
			for(curID = regionArray[i].firstID;
				curID < (regionArray[i].firstID + regionArray[i].numSub); curID++)
			{ // For each subvolume in current region
				if(!subvolArray[curID].bBoundary)
					continue; // Subvolume cannot border any regions
				
				findSubvolCoor(curSubBound, regionArray[i], subCoorInd[curID]);
				
				// Determine whether the subvolume faces the neighbor region in any direction
				if(bSubFaceRegion(regionArray, i, j, curSubBound,
					subHalfSize[i], subCoorInd[curID], boundAdjError,
					&numDir, dirArray))
				{
					// This subvolume borders microscopic region j along numDir faces
					regionArray[i].boundSubNumFace[j][curBoundID] = numDir;
					
					regionArray[i].boundVirtualNeighCoor[j][curBoundID] =
						malloc(regionArray[i].boundSubNumFace[j][curBoundID]
						* sizeof(double [3]));
					if(regionArray[i].boundVirtualNeighCoor[j][curBoundID] == NULL)
					{
						fprintf(stderr, "ERROR: Memory allocation for region %u's neighbor parameters with region %u.\n", i, j);
						exit(EXIT_FAILURE);
					}
					
					regionArray[i].neighID[j][curBoundID] = curID;
					regionArray[i].boundSubCoor[j][curBoundID][0] =
						curSubBound[0] + subHalfSize[i];
					regionArray[i].boundSubCoor[j][curBoundID][1] =
						curSubBound[2] + subHalfSize[i];
					regionArray[i].boundSubCoor[j][curBoundID][2] =
						curSubBound[4] + subHalfSize[i];
					
					for(curDir = 0;
						curDir < regionArray[i].boundSubNumFace[j][curBoundID]; curDir++)
					{
						
						switch (dirArray[curDir])
						{
							case LEFT: // New molecule goes to lower x
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][0] =
									curSubBound[0] - regionArray[i].actualSubSize;
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][1] =
									curSubBound[2];
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][2] =
									curSubBound[4];
								break;
							case RIGHT: // New molecule goes to upper x
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][0] =
									curSubBound[0] + regionArray[i].actualSubSize;
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][1] =
									curSubBound[2];
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][2] =
									curSubBound[4];
								break;
							case DOWN: // New molecule goes to lower y
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][0] =
									curSubBound[0];
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][1] =
									curSubBound[2] - regionArray[i].actualSubSize;
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][2] =
									curSubBound[4];
								break;
							case UP: // New molecule goes to upper y
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][0] =
									curSubBound[0];
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][1] =
									curSubBound[2] + regionArray[i].actualSubSize;
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][2] =
									curSubBound[4];
								break;
							case IN: // New molecule goes to lower z
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][0] =
									curSubBound[0];
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][1] =
									curSubBound[2];
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][2] =
									curSubBound[4] - regionArray[i].actualSubSize;
								break;
							case OUT: // New molecule goes to upper z
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][0] =
									curSubBound[0];
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][1] =
									curSubBound[2];
								regionArray[i].boundVirtualNeighCoor[j][curBoundID][curDir][2] =
									curSubBound[4] + regionArray[i].actualSubSize;
								break;
							default:
								fprintf(stderr, "ERROR: Invalid neighbor direction determined for neighboring regions %u and %u.\n", i, j);
						}
					}
					curBoundID++; // Increment the index of the ordered subvolume along region
									// boundary
				}
			}
		}
	}
}

// Free memory of region parameters
void delete_boundary_region_3D(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	struct region regionArray[])
{
	short i,j; // Current region
	uint32_t k; // Current subvolume along boundary
	
	if(regionArray == NULL)
		return;
	
	delete_region_chem_rxn3D(NUM_REGIONS, NUM_MOL_TYPES, regionArray);
	
	for(i = 0; i < NUM_REGIONS; i++)
	{
		if(regionArray[i].numChildren > 0)
		{
			if(regionArray[i].childrenID != NULL) free(regionArray[i].childrenID);
			if(regionArray[i].childrenCoor != NULL) free(regionArray[i].childrenCoor);
		}
		
		if(regionArray[i].regionNeighDir != NULL) free(regionArray[i].regionNeighDir);
		
		if(regionArray[i].spec.bMicro)
		{
			if(regionArray[i].regionNeighID != NULL) free(regionArray[i].regionNeighID);
			for(j = 0; j < NUM_REGIONS; j++)
			{
				if(regionArray[i].numRegionNeighFace[j] > 0)
				{
					if(regionArray[i].boundRegionFaceCoor[j] != NULL)
						free(regionArray[i].boundRegionFaceCoor[j]);
					if(regionArray[i].regionNeighFaceDir[j] != NULL)
						free(regionArray[i].regionNeighFaceDir[j]);
				}
			}
			if(regionArray[i].numRegionNeighFace != NULL)
				free(regionArray[i].numRegionNeighFace);
			if(regionArray[i].boundRegionFaceCoor != NULL)
				free(regionArray[i].boundRegionFaceCoor);
			if(regionArray[i].regionNeighFaceDir != NULL)
				free(regionArray[i].regionNeighFaceDir);
		}
		
		if(!regionArray[i].spec.bMicro)
		{
			if(regionArray[i].diffRate != NULL) free(regionArray[i].diffRate);
		}
		
		if(!regionArray[i].spec.bMicro && NUM_REGIONS > 1)
		{
			for(j = 0; j < NUM_REGIONS; j++)
			{
				if(!regionArray[j].spec.bMicro)
					continue; // Neighbor is mesoscopic; pointers below undefined
				if(!regionArray[i].isRegionNeigh[j])
					continue; // Regions don't border each other
				if(regionArray[i].numSubRegionNeigh[j] < 1)
					continue; // No neighbours found
							  // (failsafe; should not have reached here)
							  
				
				
				for(k = 0; k < regionArray[i].numSubRegionNeigh[j]; k++)
				{
					if(regionArray[i].numMolFromMicro[j][k] != NULL) free(regionArray[i].numMolFromMicro[j][k]);
					if(regionArray[i].boundVirtualNeighCoor[j][k] != NULL)	
						free(regionArray[i].boundVirtualNeighCoor[j][k]);						
				}
				if(regionArray[i].numMolFromMicro[j] != NULL)
					free(regionArray[i].numMolFromMicro[j]);
				
				if(regionArray[i].neighID[j] != NULL) free(regionArray[i].neighID[j]);
				if(regionArray[i].boundSubNumFace[j] != NULL) free(regionArray[i].boundSubNumFace[j]);
				if(regionArray[i].boundSubCoor[j] != NULL)	
					free(regionArray[i].boundSubCoor[j]);
				if(regionArray[i].boundVirtualNeighCoor[j] != NULL)	
					free(regionArray[i].boundVirtualNeighCoor[j]);
				if(regionArray[i].bNeedUpdate[j] != NULL) free(regionArray[i].bNeedUpdate[j]);
			}
			if(regionArray[i].neighID != NULL) free(regionArray[i].neighID);
			if(regionArray[i].boundSubNumFace != NULL) free(regionArray[i].boundSubNumFace);
			if(regionArray[i].boundSubCoor != NULL) free(regionArray[i].boundSubCoor);
			if(regionArray[i].boundVirtualNeighCoor != NULL) free(regionArray[i].boundVirtualNeighCoor);
			if(regionArray[i].bNeedUpdate != NULL) free(regionArray[i].bNeedUpdate);
			if(regionArray[i].numMolFromMicro != NULL) free(regionArray[i].numMolFromMicro);
		}
		
		if(regionArray[i].numSubRegionNeigh != NULL) free(regionArray[i].numSubRegionNeigh);
		if(regionArray[i].isRegionNeigh != NULL) free(regionArray[i].isRegionNeigh);
		
	}
}

// Find volume of specified region (exluding volume of children)
double findRegionVolume(const struct region regionArray[],
	const short curRegion,
	bool bOuter)
{
	short curChild;
	double volume = boundaryArea(regionArray[curRegion].spec.shape, regionArray[curRegion].boundary);
	
	if(!bOuter)
	{
		// Subtract outer volumes of children
		for(curChild = 0; curChild < regionArray[curRegion].numChildren; curChild++)
		{
			volume -= boundaryArea(
				regionArray[regionArray[curRegion].childrenID[curChild]].spec.shape, regionArray[regionArray[curRegion].childrenID[curChild]].boundary);
		}
	}
	
	return volume;
}
	
// Count the cumulative number of subvolumes defined by subvol_spec
uint32_t count_subvol3D(const struct region regionArray[],
	const short NUM_REGIONS)
{
	short i; // Current Region
	uint32_t num_subvolumes = 0;
	for(i = 0; i < NUM_REGIONS; i++)
		num_subvolumes += (uint32_t) regionArray[i].numSub;
	return num_subvolumes;
}

// Determine which regions of subvolumes are adjacent
void findRegionTouch3D(const short NUM_REGIONS,
	const struct spec_region3D subvol_spec[],
	struct region regionArray[],
	const double SUBVOL_BASE_SIZE)
{
	short int i,j,k; // Current Region
	short int curNeigh; // Current neighbor region
	
	unsigned short direction = 0;
	double error = SUBVOL_BASE_SIZE * SUB_ADJ_RESOLUTION;
	
	for(i = 0; i < NUM_REGIONS; i++){
		for(j = 0; j < NUM_REGIONS; j++){
			if(i == j)
			{ // Region is itself and not its neighbor
				continue;
			}
			
			if(regionArray[i].parentID == j)
			{
				// Region j is region i's parent
				regionArray[i].isRegionNeigh[j] = true;
				regionArray[i].regionNeighDir[j] = PARENT;
				regionArray[i].numRegionNeigh++;		
				continue;
			}
			
			if(regionArray[j].parentID == i)
			{
				// Region i is region j's parent
				regionArray[i].isRegionNeigh[j] = true;
				regionArray[i].regionNeighDir[j] = CHILD;
				regionArray[i].numRegionNeigh++;
				continue;
			}
			
			if(regionArray[i].spec.shape == SPHERE ||
				regionArray[j].spec.shape == SPHERE)
			{ // Spherical regions can only touch parent/child regions
				continue;
			}
			
			if(bBoundaryAdjacent(RECTANGULAR_BOX, regionArray[i].boundary, RECTANGULAR_BOX,
				regionArray[j].boundary, error, &direction))
			{	// Regions share a face
				regionArray[i].boundaryRegion[direction] = true;
				regionArray[i].boundaryRegionMicro[direction] = regionArray[j].spec.bMicro;
				regionArray[i].isRegionNeigh[j] = true;
				regionArray[i].regionNeighDir[j] = direction;
				regionArray[i].numRegionNeigh++;
			}
		}
	}
	
	// Make another pass to write overlap coordinates for microscopic regions
	for(i = 0; i < NUM_REGIONS; i++)
	{
		if(!regionArray[i].spec.bMicro)
			continue;
		
		// Assign memory for region neighbor IDs
		regionArray[i].regionNeighID =
			malloc(regionArray[i].numRegionNeigh* sizeof(short));
		if(regionArray[i].regionNeighID == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation for region %u's microscopic parameters.\n", i);
			exit(EXIT_FAILURE);
		}
		curNeigh = 0;
		
		for(j = 0; j < NUM_REGIONS; j++)
		{
			// Initialize number of faces that this region shares with every other region
			regionArray[i].numRegionNeighFace[j] = 0;
			
			if(regionArray[i].isRegionNeigh[j])
			{ // Region i is microscopic and region j is a neighbor
				regionArray[i].regionNeighID[curNeigh++] = j;
				
				switch(regionArray[i].regionNeighDir[j])
				{
					case PARENT:
						// Number of neighboring faces is equal to number of region i faces
						switch(regionArray[i].spec.shape)
						{
							case RECTANGULAR_BOX:
								regionArray[i].numRegionNeighFace[j] = 6;
								break;
							case SPHERE:
								regionArray[i].numRegionNeighFace[j] = 1;
								break;
						}
						break;
					case CHILD:
						// Number of neighboring faces is equal to number of region j faces
						switch(regionArray[j].spec.shape)
						{
							case RECTANGULAR_BOX:
								regionArray[i].numRegionNeighFace[j] = 6;
								break;
							case SPHERE:
								regionArray[i].numRegionNeighFace[j] = 1;
								break;
						}
						break;
					default:
						regionArray[i].numRegionNeighFace[j] = 1;
						break;
				}
				
				regionArray[i].boundRegionFaceCoor[j] = 
					malloc(regionArray[i].numRegionNeighFace[j]* sizeof(double [6]));
				regionArray[i].regionNeighFaceDir[j] = 
					malloc(regionArray[i].numRegionNeighFace[j]* sizeof(short));
				if(regionArray[i].boundRegionFaceCoor[j] == NULL
					|| regionArray[i].regionNeighFaceDir[j] == NULL)
				{
					fprintf(stderr, "ERROR: Memory allocation for region %u's neighbor parameters with region %u.\n", i, j);
					exit(EXIT_FAILURE);
				}
				
				// Determine coordinates of boundary overlap
				switch(regionArray[i].regionNeighDir[j])
				{
					case PARENT:
						// Record all faces of region i as boundary to region j
						// NOTE: This isn't strictly true for cases of adjacent children,
						// but this won't interfere with the correctness of the microscopic
						// simulation
						for(k = 0; k < regionArray[i].numRegionNeighFace[j]; k++)
						{
							recordFace(regionArray[i].spec.shape, regionArray[i].boundary,
								k, regionArray[i].boundRegionFaceCoor[j][k]);
							regionArray[i].regionNeighFaceDir[j][k] = k;
						}
						break;
					case CHILD:
						// Record all faces of region j as boundary to region j
						// NOTE: This isn't strictly true for cases of adjacent children,
						// but this won't interfere with the correctness of the microscopic
						// simulation
						for(k = 0; k < regionArray[i].numRegionNeighFace[j]; k++)
						{
							recordFace(regionArray[j].spec.shape, regionArray[j].boundary,
								k, regionArray[i].boundRegionFaceCoor[j][k]);
							regionArray[i].regionNeighFaceDir[j][k] = k;
						}
						break;
					default:
						// Calculate actual boundary
						intersectBoundary(regionArray[i].spec.shape, regionArray[i].boundary,
							regionArray[j].spec.shape, regionArray[j].boundary,
							regionArray[i].boundRegionFaceCoor[j][0]);
						regionArray[i].regionNeighFaceDir[j][0] =
							regionArray[i].regionNeighDir[j];
						break;
				}
			}
		}
	}
	
	return;
}

// Find index of desired subvolume in list defining region's boundary with another region
uint32_t find_sub_in_bound_list3D(const short curRegion,
	const short destRegion,
	const struct region regionArray[],
	uint32_t curSub)
{
	uint32_t neighID = 0;
	
	while(regionArray[curRegion].neighID[destRegion][neighID] != curSub)
	{
		neighID++;
		if(!(neighID < regionArray[curRegion].numSubRegionNeigh[destRegion]))
			return UINT32_MAX;
	}
	return neighID;
}

// Find closest region for point to be in
unsigned short findNearestValidRegion(const double point[],
	const unsigned short sourceRegion,
	bool * bInRegion,
	const unsigned short NUM_REGIONS,
	const struct region regionArray[])
{
	unsigned short curRegion;
	unsigned short minRegion;
	
	double curDist;
	double minDist = INFINITY;
	
	// First check whether we are inside any region
	for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
	{
		if((regionArray[sourceRegion].isRegionNeigh[curRegion]
			|| curRegion == sourceRegion)
			&& bPointInRegionNotChild(curRegion, regionArray, point))
		{
			// TODO: Check whether reaching this region from source is valid
			// (e.g., is destination reachable from source in one "jump"?)
			// For now, must be neighbor of source region
			*bInRegion = true;
			return curRegion;
		}
	}
	
	// We are not actually inside a region. Find distance to closest region.
	*bInRegion = false;
	for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
	{
		curDist = distanceToBoundary(point, regionArray[curRegion].spec.shape,
			regionArray[curRegion].boundary);
		
		if((regionArray[sourceRegion].isRegionNeigh[curRegion]
			|| curRegion == sourceRegion)
			&& curDist < minDist)
		{
			minDist = curDist;
			minRegion = curRegion;
		}
	}
	return minRegion;
}

// Find the closest subvolume in current region that is along boundary
// of specified neighbor region
uint32_t findNearestSub3D(const short curRegion,
	const struct region regionArray[],
	const short neighRegion,
	double x,
	double y,
	double z)
{
	// Number of subvolumes that border the specified neighbor region
	uint32_t numSub = regionArray[curRegion].numSubRegionNeigh[neighRegion];
	
	uint32_t curSub;
	uint32_t minSub;
	double minDistSq = DBL_MAX;
	double curDistSq;
	for(curSub = 0; curSub < numSub; curSub++)
	{
		curDistSq = (x - regionArray[curRegion].boundSubCoor[neighRegion][curSub][0])
			*(x - regionArray[curRegion].boundSubCoor[neighRegion][curSub][0]);
		if (curDistSq > minDistSq)
			continue; // No need to check y-coordinate
		
		curDistSq += (y - regionArray[curRegion].boundSubCoor[neighRegion][curSub][1])
			*(y - regionArray[curRegion].boundSubCoor[neighRegion][curSub][1]);
		if (curDistSq > minDistSq)
			continue; // No need to check z-coordinate
		
		curDistSq += (z - regionArray[curRegion].boundSubCoor[neighRegion][curSub][2])
			*(z - regionArray[curRegion].boundSubCoor[neighRegion][curSub][2]);
		
		if(curDistSq < minDistSq)
		{ // This subvolume is the closest so far
			minSub = curSub;
			minDistSq = curDistSq;
		}
	}
	return regionArray[curRegion].neighID[neighRegion][minSub];
}

// Volume of intersection within specified region
// Value returned here excludes children
// Intersection must be rectangular. Region can be spherical if it
// fully contains or is contained by the other boundary
// Round children must be entirely inside or outside intersection
double intersectRegionVolume(const short curRegion,
	const struct region regionArray[],
	const int boundary2Type,
	const double boundary2[])
{
	short curChild, childRegion;
	double interBoundary[6];
	double volume, childVolume;
	int intersectType;
	
	// Calculate intersection area including children
	intersectType = intersectBoundary(regionArray[curRegion].spec.shape,
		regionArray[curRegion].boundary,	boundary2Type, boundary2, interBoundary);
	if(intersectType == UNDEFINED_SHAPE)
	{
		// Invalid intersection of boundary with region		
		fprintf(stderr, "ERROR: An intersection with region %u is invalid.\n", curRegion);
		exit(EXIT_FAILURE);			
	}
	volume = boundaryArea(intersectType, interBoundary);
	
	if(volume <= 0)
		return 0.;	// Region does not intersect at all
	
	// Subtract area populated by children
	for(curChild = 0; curChild < regionArray[curRegion].numChildren; curChild++)
	{
		childRegion = regionArray[curRegion].childrenID[curChild];
		intersectType = intersectBoundary(regionArray[childRegion].spec.shape,
			regionArray[childRegion].boundary, boundary2Type, boundary2, interBoundary);
		
		if(intersectType == UNDEFINED_SHAPE)
		{
			// Invalid intersection of boundary with child		
			fprintf(stderr, "ERROR: An intersection with region %u is invalid.\n", curRegion);
			exit(EXIT_FAILURE);			
		} else
		{
			childVolume = boundaryArea(intersectType, interBoundary);
			if(childVolume > 0.)
				volume -= childVolume; // Only subtract positive child volumes
		}
	}
	
	return volume;
}

// Does a line hit the specified region? If so then where and at what distance?
bool bLineHitRegion(const double p1[3],
	const double L[3],
	const double length,
	const short startRegion,
	const short endRegion,
	const struct region regionArray[],
	short * planeID,
	double * d,
	double intersectPoint[3])
{	
	short curPlane;
	double minDist = INFINITY;
	double nearestIntersectPoint[3];
	bool bIntersect = false;
	int boundary1Type;
	bool bInside;
	
	if(regionArray[startRegion].parentID == endRegion)
	{ // Boundary will be shape of the child
		boundary1Type = regionArray[startRegion].spec.shape;
		bInside = true;
	} else
	{ // Boundary will be shape of the destination region
		boundary1Type = regionArray[endRegion].spec.shape;
		bInside = false;
	}
	
	switch(boundary1Type)
	{
		case RECTANGULAR_BOX:
			for(curPlane = 0;
				curPlane < regionArray[startRegion].numRegionNeighFace[endRegion];
				curPlane++)
			{
				if(bLineHitInfinitePlane(p1, L, length, RECTANGULAR_BOX,
					regionArray[startRegion].boundRegionFaceCoor[endRegion][curPlane],
					regionArray[startRegion].regionNeighFaceDir[endRegion][curPlane],
					false, d, intersectPoint)
					&& bPointOnFace(intersectPoint, RECTANGULAR_BOX,
					regionArray[startRegion].boundRegionFaceCoor[endRegion][curPlane], regionArray[startRegion].regionNeighFaceDir[endRegion][curPlane])
					&& *d < minDist)
				{ // Line does intersect this face at a valid distance and it is closest
					bIntersect = true;
					minDist = *d;
					*planeID = regionArray[startRegion].regionNeighFaceDir[endRegion][curPlane];
					nearestIntersectPoint[0] = intersectPoint[0];
					nearestIntersectPoint[1] = intersectPoint[1];
					nearestIntersectPoint[2] = intersectPoint[2];
				}
			}
			if(bIntersect)
			{
				*d = minDist;
				intersectPoint[0] = nearestIntersectPoint[0];
				intersectPoint[1] = nearestIntersectPoint[1];
				intersectPoint[2] = nearestIntersectPoint[2];
				return true;
			} else return false;
		case SPHERE:
			return bLineHitInfinitePlane(p1, L, length, SPHERE,
				regionArray[startRegion].boundRegionFaceCoor[endRegion][0],
				regionArray[startRegion].regionNeighFaceDir[endRegion][0], bInside, d, intersectPoint);
		default:
			fprintf(stderr,"ERROR: Boundary type %d invalid for boundary intersection between regions %u and %u.\n", boundary1Type, startRegion, endRegion);
			return false;		
	}
}

bool bPointInRegionNotChild(const short curRegion,
	const struct region regionArray[],
	const double point[3])
{
	short curChild;
	
	if(bPointInBoundary(point, regionArray[curRegion].spec.shape,
		regionArray[curRegion].boundary))
	{ // Point is within the region's outer boundary
		for(curChild = 0; curChild < regionArray[curRegion].numChildren; curChild++)
		{
			if(bPointInBoundary(point, regionArray[regionArray[curRegion].childrenID[curChild]].spec.shape,
				regionArray[regionArray[curRegion].childrenID[curChild]].boundary))
				return false;
		}
		return true;
	}
	
	return false;
}

// Is point in a region or any of its nested children?
// If true, actualRegion will be the ID of the region with the point
bool bPointInRegionOrChild(const short curRegion,
	const struct region regionArray[],
	const double point[3],
	short * actualRegion)
{
	short curChild;
	
	if(bPointInBoundary(point, regionArray[curRegion].spec.shape,
		regionArray[curRegion].boundary))
	{ // Point is within the region's outer boundary
		for(curChild = 0; curChild < regionArray[curRegion].numChildren; curChild++)
		{
			if(bPointInRegionOrChild(regionArray[curRegion].childrenID[curChild],
				regionArray, point, actualRegion))
			{
				return true;
			}
		}
		// We are actually in curRegion
		* actualRegion = curRegion;
		return true;
	}
	
	return false;
}

// Is a specific face shared between two regions?
bool bSharedBoundary(const short startRegion,
	const short endRegion,
	const struct region regionArray[],
	const short faceID)
{
	if(regionArray[startRegion].spec.shape !=
		regionArray[endRegion].spec.shape)
		return false; // Regions with different shapes cannot share a boundary
	
	switch(regionArray[startRegion].spec.shape)
	{
		case RECTANGULAR_BOX:
			return fabs(regionArray[startRegion].boundary[faceID]
				- regionArray[endRegion].boundary[faceID]) <
				regionArray[startRegion].subResolution;
		case SPHERE:
			// Only one face to compare
			return regionArray[startRegion].boundary[3]
				== regionArray[endRegion].boundary[3];
	}
}

// Lock point coordinate to chosen region face
void lockPointToRegion(double point[3],
	const short startRegion,
	const short endRegion,
	const struct region regionArray[],
	const short faceID)
{
	short boundRegion;
	if (regionArray[startRegion].isRegionNeigh[endRegion]
		&& regionArray[startRegion].regionNeighDir[endRegion] == CHILD)
		boundRegion = endRegion;
	else
		boundRegion = startRegion;
	
	switch(regionArray[boundRegion].spec.shape)
	{
		case RECTANGULAR_BOX:
			if(faceID < 2)
				point[0] = regionArray[boundRegion].boundary[faceID];
			else if (faceID > 3)
				point[2] = regionArray[boundRegion].boundary[faceID];
			else					
				point[1] = regionArray[boundRegion].boundary[faceID];	
			break;
		case SPHERE:
			// Arbitrarily adjust x-coordinate to lock point to surface
			if(point[0] > regionArray[boundRegion].boundary[0])
				point[0] = sqrt(squareDBL(regionArray[boundRegion].boundary[3]) -
					squareDBL(regionArray[boundRegion].boundary[1] - point[1]) -
					squareDBL(regionArray[boundRegion].boundary[2] - point[2]));
			else
				point[0] = -sqrt(squareDBL(regionArray[boundRegion].boundary[3]) -
					squareDBL(regionArray[boundRegion].boundary[1] - point[1]) -
					squareDBL(regionArray[boundRegion].boundary[2] - point[2]));
			point[0] += regionArray[boundRegion].boundary[0];
			break;
		default:
			fprintf(stderr,"ERROR: Point cannot be locked to region %u because it is of type %d.\n", boundRegion, regionArray[boundRegion].spec.shape);
	}
}

// Generate a random cartesian point in the specified region
void generatePointInRegion(const short curRegion,
	const struct region regionArray[],
	double point[3])
{
	bool bNeedPoint = true;
	
	while(bNeedPoint)
	{
		uniformPointVolume(point, regionArray[curRegion].spec.shape,
			regionArray[curRegion].boundary);
		if(bPointInRegionNotChild(curRegion, regionArray, point))
			bNeedPoint = false;
	}
}

// Which region contains given point, excluding children?
short findRegionNotChild(const short NUM_REGIONS,
	const struct region regionArray[],
	double point[3])
{
	short curRegion;
	
	for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
	{
		if(bPointInRegionNotChild(curRegion, regionArray, point))
			return curRegion;
	}
	
	// Point was not found in any region
	return SHRT_MAX;
}

// Does a subvolume face a region? If yes, then along which faces?
// Assert that current subvolume is along its own region boundary, and that
// neighbor region is microscopic
bool bSubFaceRegion(struct region regionArray[],
	const short curRegion,
	const short neighRegion,
	const double curSubBound[6],
	const double subHalfSize,
	const uint32_t subCoorInd[3],
	double boundAdjError,
	unsigned short * numFace,
	unsigned short dirArray[6])
{
	double curNeighBound[3];
	int i;
	
	for(i = 0; i < 3; i++)
	{
		curNeighBound[i] = 0.;
	}
	
	curNeighBound[0] = curSubBound[0] + subHalfSize;
	curNeighBound[1] = curSubBound[2] + subHalfSize;
	curNeighBound[2] = curSubBound[4] + subHalfSize;
	
	*numFace = 0;
	
	// Check to see if point just to "left" is within neighbor
	curNeighBound[0] = curSubBound[0] - boundAdjError;
	if(bPointInRegionNotChild(neighRegion, regionArray,
		curNeighBound))
		dirArray[(*numFace)++] = LEFT;
	curNeighBound[0] = curSubBound[0] + subHalfSize;
	
	// Check to see if point just to "right" is within neighbor
	curNeighBound[0] = curSubBound[1] + boundAdjError;
	if(bPointInRegionNotChild(neighRegion, regionArray,
		curNeighBound))
		dirArray[(*numFace)++] = RIGHT;
	curNeighBound[0] = curSubBound[0] + subHalfSize;
	
	// Check to see if point just to "down" is within neighbor
	curNeighBound[1] = curSubBound[2] - boundAdjError;
	if(bPointInRegionNotChild(neighRegion, regionArray,
		curNeighBound))
		dirArray[(*numFace)++] = DOWN;
	curNeighBound[1] = curSubBound[2] + subHalfSize;
	
	// Check to see if point just to "up" is within neighbor
	curNeighBound[1] = curSubBound[3] + boundAdjError;
	if(bPointInRegionNotChild(neighRegion, regionArray,
		curNeighBound))
		dirArray[(*numFace)++] = UP;
	curNeighBound[1] = curSubBound[2] + subHalfSize;
	
	// Check to see if point just to "in" is within neighbor
	curNeighBound[2] = curSubBound[4] - boundAdjError;
	if(bPointInRegionNotChild(neighRegion, regionArray,
		curNeighBound))
		dirArray[(*numFace)++] = IN;
	curNeighBound[2] = curSubBound[4] + subHalfSize;
	
	// Check to see if point just to "out" is within neighbor
	curNeighBound[2] = curSubBound[5] + boundAdjError;
	if(bPointInRegionNotChild(neighRegion, regionArray,
		curNeighBound))
		dirArray[(*numFace)++] = OUT;
	curNeighBound[2] = curSubBound[4] + subHalfSize;
	
	
	return *numFace > 0;
}