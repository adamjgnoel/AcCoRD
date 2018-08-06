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
 * Last revised for AcCoRD v1.4 (2018-08-06)
 *
 * Revision history:
 *
 * Revision v1.4 (2018-08-06)
 * - added a priori monte carlo (APMC) absorption algorithm as a new surface
 * reaction type. Includes settings for how to define the a priori absorption
 * probability calculation and whether/how to apply a threshold to turn it off
 * - corrected check on molecule transitions between parent and child regions
 * when a molecule is initialized at the parent/child boundary
 *
 * Revision v1.1 (2016-12-24)
 * - added members defining flow parameters for every molecule
 * - added member to denote which molecules can diffuse
 *
 * Revision v1.0 (2016-10-31)
 * - added local diffusion coefficients that can apply to particular region
 * - added specifying diffusion coefficient that applies to specific surface
 * interaction reactions.
 * - moved mesoscopic structure fields from subvolume struct to meso subvolume struct
 *
 * Revision v0.7.0.1 (public beta, 2016-08-30)
 * - corrected calculating region volume when a normal region has a surface child
 * inside
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - added members to track the direction of microscopic subvolumes from mesoscopic
 * subvolumes. Replaced array for storing coordinates of virtual microscopic
 * subvolume with array of boundary coordinates of boundary mesoscopic subvolumes. These
 * changes needed to accommodate improved hybrid transition algorithms
 * - added members to implement bimolecular chemical reactions in microscopic regime.
 * Use a preliminary algorithm where the user must supply the binding and unbinding radii
 * - changed findNearestSub function to return index of subvolume in region's neighID array
 * instead of the global subvolume list. This makes function suitable for more calls.
 * - added case for point boundary when determining intersection with a region
 *
 * Revision v0.5.1 (2016-05-06)
 * - updated function bPointInRegionNotChild to take an extra input to indicate
 * whether to ignore children regions that are surfaces
 * - added bReleaseProduct for surface reactions to indicate whether products are
 * released from the surface
 * - added members to store unique properties needed for absorbing and
 * desorbing reactions
 *
 * Revision v0.5 (2016-04-15)
 * - re-structured region array initialization to nest more code in functions
 * - added more checks on region parameters (including label uniqueness) to verify placement
 * - pushed error exit to end of region initialization so that all errors will be displayed
 * before exiting.
 * - added 2D regions
 * - added type member to spec and plane, dimension, and effectiveDim members to main
 * struct in order to accommodate surface and other 2D regions
 * - added region label when giving errors about region initialization
 * - adjusted clearance between spherical and rectangular regions such that the clearance
 * between them (when one is nested inside the other) is scaled by the subvolume adjacency
 * error
 * - corrected "locking" to spherical region so that the coordinate that is "locked" is
 * the one that is the furthest from the center of the sphere
 *
 * Revision v0.4.1
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
	short i,j; // Current region
	short m;
	unsigned short curMolType;
	double h_i;
	bool bFail = false;	// Fail switch for placement errors (not memory-related)
	
	// Check for unique region labels
	for(i = 0; i < NUM_REGIONS; i++)
	{		
		if(strlen(subvol_spec[i].label) == 0)
			continue; // Current region has no label
		
		for(j = i+1; j < NUM_REGIONS; j++)
		{
			if(strlen(subvol_spec[j].label) == 0)
				continue; // Current region has no label
			
			if(!strcmp(subvol_spec[i].label, subvol_spec[j].label))
			{
				fprintf(stderr, "ERROR: Regions %u and %u have the same label \"%s\".\n",
					i, j, subvol_spec[i].label);
				bFail = true;
			}
		}
	}
	
	// Copy or intialize main region parameters
	// Determine region boundary
	// Determine diffusion rate within region
	for(i = 0; i < NUM_REGIONS; i++)
	{
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
			
			if(subvol_spec[i].shape == RECTANGLE)
			{ // Determine which plane the rectangle is in
				regionArray[i].dimension = DIM_2D;
				if(regionArray[i].actualSubSize*subvol_spec[i].numX == 0)
				{
					regionArray[i].plane = PLANE_YZ;
				} else if (regionArray[i].actualSubSize*subvol_spec[i].numY == 0)
				{
					regionArray[i].plane = PLANE_XZ;
				} else if (regionArray[i].actualSubSize*subvol_spec[i].numZ == 0)
				{
					regionArray[i].plane = PLANE_XY;
				} else
				{ // We shouldn't have gotten here
					fprintf(stderr, "ERROR: Region %u (Label: \"%s\") is not defined properly as a Rectangle.\n", i,	subvol_spec[i].label);
					fprintf(stderr, "There must be 0 subvolumes along one dimension.\n");
					bFail = true;
				}
			} else
			{ // There is no single plane
				regionArray[i].dimension = DIM_3D;
				regionArray[i].plane = PLANE_3D;
			}
		} else
		{
			regionArray[i].actualSubSize = subvol_spec[i].radius;
			regionArray[i].boundary[0] = subvol_spec[i].xAnch;
			regionArray[i].boundary[1] = subvol_spec[i].yAnch;
			regionArray[i].boundary[2] = subvol_spec[i].zAnch;
			regionArray[i].boundary[3] = subvol_spec[i].radius;
			regionArray[i].boundary[4] = squareDBL(subvol_spec[i].radius);
			regionArray[i].boundary[5] = 0.;
			regionArray[i].dimension = DIM_3D;
			regionArray[i].plane = PLANE_3D;
		}
		
		switch(regionArray[i].spec.type)
		{
			case REGION_NORMAL:
				regionArray[i].subShape = regionArray[i].spec.shape;
				regionArray[i].effectiveDim = regionArray[i].dimension;
				if(regionArray[i].plane == PLANE_3D)
				{
					regionArray[i].numFace = 0;
				}
				else
				{
					regionArray[i].numFace = 1;
				}
				break;
			case REGION_SURFACE_3D:
				regionArray[i].effectiveDim = DIM_2D;
				if(subvol_spec[i].shape == RECTANGULAR_BOX)
				{
					regionArray[i].numFace = 6;
					regionArray[i].subShape = RECTANGLE;
				}
				else
				{
					regionArray[i].numFace = 1;
					regionArray[i].subShape = regionArray[i].spec.shape;
				}
				break;
			case REGION_SURFACE_2D:
				regionArray[i].effectiveDim = DIM_1D;
				if(subvol_spec[i].shape == RECTANGLE)
				{
					regionArray[i].numFace = 4;
					regionArray[i].subShape = LINE;
				}
				else
				{
					regionArray[i].numFace = 1;
					regionArray[i].subShape = regionArray[i].spec.shape;
				}
				break;
		}
			
		regionArray[i].numChildren = 0;
		regionArray[i].subResolution = SUBVOL_BASE_SIZE * SUB_ADJ_RESOLUTION;
		regionArray[i].bHasMesoNeigh = false;
		
		// Calculate diffusion rates within region
		if(!regionArray[i].spec.bMicro)
		{
			regionArray[i].diffRate = malloc(NUM_MOL_TYPES*sizeof(double));
			if(regionArray[i].diffRate == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\") diffusion coefficients.\n", i, subvol_spec[i].label);
				exit(EXIT_FAILURE);
			}
			for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
			{
				h_i = regionArray[i].actualSubSize;
				regionArray[i].diffRate[curMolType] =
					DIFF_COEF[i][curMolType]/h_i/h_i;
			}
		}
		
		// Determine which types of molecules cannot diffuse
		regionArray[i].bDiffuse = malloc(NUM_MOL_TYPES*sizeof(bool));
		if(regionArray[i].bDiffuse == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\") diffusion boolean.\n", i, subvol_spec[i].label);
			exit(EXIT_FAILURE);
		}
		for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
		{
			regionArray[i].bDiffuse[curMolType] =
				DIFF_COEF[i][curMolType] > 0.;
		}
	}
	
	// Allocate and initialize flow parameters
	initializeRegionFlow(NUM_REGIONS, NUM_MOL_TYPES, regionArray, subvol_spec);
	
	// Allocate memory for each region's neighbors
	allocateRegionNeighbors(NUM_REGIONS, regionArray);
	
	// Determine region nesting (find each region's parent and children)
	// Includes some (but not exhaustive) checks on spec.type
	initializeRegionNesting(NUM_REGIONS, regionArray, &bFail);
	
	// Determine what regions are touching (reduces neighbour search computation)
	// Includes some (but not exhaustive) checks on spec.type
	findRegionTouch(NUM_REGIONS, subvol_spec, regionArray, SUBVOL_BASE_SIZE);
	
	// Confirm validity of all regions now that neighbors have been identified
	validateRegions(NUM_REGIONS, regionArray, &bFail);
	
	if(bFail)
	{ // There was at least one error in the region configuration
		exit(EXIT_FAILURE);
	}
	
	// Determine number of subvolumes in each region
	findNumRegionSubvolumes(NUM_REGIONS, regionArray);
	
	// Define chemical reaction network
	initializeRegionChemRxn(NUM_REGIONS, regionArray,
		NUM_MOL_TYPES, MAX_RXNS, chem_rxn, DIFF_COEF);
}

// Initialize region knowledge of the subvolumes that are adjacent to it
// This function is called by build_subvol_array3D in subvolume.c
void initializeRegionSubNeighbor(struct region regionArray[],
	const struct spec_region3D subvol_spec[],
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
		regionArray[i].boundSubCenterCoor =
			malloc(NUM_REGIONS * sizeof(regionArray[i].boundSubCenterCoor));
		regionArray[i].boundSubCoor =
			malloc(NUM_REGIONS * sizeof(regionArray[i].boundSubCoor));
		regionArray[i].boundVirtualNeighDir =
			malloc(NUM_REGIONS * sizeof(regionArray[i].boundVirtualNeighDir));
		regionArray[i].bNeedUpdate =
			malloc(NUM_REGIONS * sizeof(regionArray[i].bNeedUpdate));
		regionArray[i].numMolFromMicro =
			malloc(NUM_REGIONS * sizeof(regionArray[i].numMolFromMicro));
		if(regionArray[i].neighID == NULL
			|| regionArray[i].boundSubNumFace == NULL
			|| regionArray[i].boundSubCenterCoor == NULL
			|| regionArray[i].boundSubCoor == NULL
			|| regionArray[i].boundVirtualNeighDir == NULL
			|| regionArray[i].bNeedUpdate == NULL
			|| regionArray[i].numMolFromMicro == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\") neighbor parameters.\n", i, subvol_spec[i].label);
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
				
				if(bSubFaceRegion(regionArray, j, curSubBound,
					boundAdjError, &numDir, dirArray))
				{
					// This subvolume borders microscopic region j along numDir faces
					regionArray[i].numSubRegionNeigh[j]++;
				}
			}
			
			if(regionArray[i].numSubRegionNeigh[j] < 1)
				continue; // No neighbours found
						  // We could reach here if a surface region prevents all possible
						  // subvolumes from being neighbors
			
			regionArray[j].bHasMesoNeigh = true;
			
			regionArray[i].neighID[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(uint32_t));
			regionArray[i].boundSubNumFace[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(unsigned short));
			regionArray[i].boundSubCenterCoor[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(double [3]));
			regionArray[i].boundSubCoor[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(double [6]));
			regionArray[i].boundVirtualNeighDir[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(regionArray[i].boundVirtualNeighDir[j]));
			regionArray[i].bNeedUpdate[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(bool));
			regionArray[i].numMolFromMicro[j] = malloc(regionArray[i].numSubRegionNeigh[j]
				* sizeof(regionArray[i].numMolFromMicro[j]));
			if(regionArray[i].neighID[j] == NULL
				|| regionArray[i].boundSubNumFace[j] == NULL
				|| regionArray[i].boundSubCenterCoor[j] == NULL
				|| regionArray[i].boundSubCoor[j] == NULL
				|| regionArray[i].boundVirtualNeighDir[j] == NULL
				|| regionArray[i].bNeedUpdate[j] == NULL
				|| regionArray[i].numMolFromMicro[j] == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\")'s neighbor parameters with region %u (label: \"%s\").\n", i, subvol_spec[i].label, j, subvol_spec[j].label);
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
					fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\")'s neighbor parameters with region %u (label: \"%s\").\n", i, subvol_spec[i].label, j, subvol_spec[j].label);
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
				if(bSubFaceRegion(regionArray, j, curSubBound,
					boundAdjError, &numDir, dirArray))
				{
					// This subvolume borders microscopic region j along numDir faces
					regionArray[i].boundSubNumFace[j][curBoundID] = numDir;
					
					regionArray[i].boundVirtualNeighDir[j][curBoundID] =
						malloc(regionArray[i].boundSubNumFace[j][curBoundID]
						* sizeof(unsigned short));
					if(regionArray[i].boundVirtualNeighDir[j][curBoundID] == NULL)
					{
						fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\")'s neighbor parameters with region %u (label: \"%s\").\n", i, subvol_spec[i].label, j, subvol_spec[j].label);
						exit(EXIT_FAILURE);
					}
					
					regionArray[i].neighID[j][curBoundID] = subvolArray[curID].mesoID;
					regionArray[i].boundSubCenterCoor[j][curBoundID][0] =
						(curSubBound[0] + curSubBound[1])/2;
					regionArray[i].boundSubCenterCoor[j][curBoundID][1] =
						(curSubBound[2] + curSubBound[3])/2;
					regionArray[i].boundSubCenterCoor[j][curBoundID][2] =
						(curSubBound[4] + curSubBound[5])/2;
					
					for(k = 0; k < 6; k++)
					{
						regionArray[i].boundSubCoor[j][curBoundID][k] =
							curSubBound[k];
					}
					
					for(curDir = 0;
						curDir < regionArray[i].boundSubNumFace[j][curBoundID]; curDir++)
					{
						regionArray[i].boundVirtualNeighDir[j][curBoundID][curDir] =
							dirArray[curDir];
					}
					curBoundID++; // Increment the index of the ordered subvolume along region
									// boundary
				}
			}
		}
	}
}

// Free memory of region parameters
void delete_boundary_region_(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	struct region regionArray[])
{
	short i,j; // Current region
	unsigned short curMolType;
	uint32_t k; // Current subvolume along boundary
	
	if(regionArray == NULL)
		return;
	
	deleteRegionChemRxn(NUM_REGIONS, NUM_MOL_TYPES, regionArray);
	
	for(i = 0; i < NUM_REGIONS; i++)
	{
		if(regionArray[i].spec.bFlow != NULL) free(regionArray[i].spec.bFlow);
		if(regionArray[i].spec.flowType != NULL) free(regionArray[i].spec.flowType);
		if(regionArray[i].spec.flowVector != NULL)
		{
			for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
			{
				if(regionArray[i].spec.flowVector[curMolType] != NULL)
					free(regionArray[i].spec.flowVector[curMolType]);
			}
			free(regionArray[i].spec.flowVector);
		}
		
		if(regionArray[i].numChildren > 0)
		{
			if(regionArray[i].childrenID != NULL) free(regionArray[i].childrenID);
			if(regionArray[i].childrenCoor != NULL) free(regionArray[i].childrenCoor);
		}
		
		if(regionArray[i].regionNeighDir != NULL) free(regionArray[i].regionNeighDir);
		
		if(regionArray[i].spec.bMicro)
		{
			if(regionArray[i].flowConstant != NULL)
			{
				for(j = 0; j < NUM_MOL_TYPES; j++)
				{
					if(regionArray[i].flowConstant[j] != NULL)
						free(regionArray[i].flowConstant[j]);
				}
				free(regionArray[i].flowConstant);
			}
			
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
					if(regionArray[i].numMolFromMicro[j][k] != NULL)
						free(regionArray[i].numMolFromMicro[j][k]);
					if(regionArray[i].boundVirtualNeighDir[j][k] != NULL)	
						free(regionArray[i].boundVirtualNeighDir[j][k]);
				}
				if(regionArray[i].numMolFromMicro[j] != NULL)
					free(regionArray[i].numMolFromMicro[j]);
				
				if(regionArray[i].neighID[j] != NULL) free(regionArray[i].neighID[j]);
				if(regionArray[i].boundSubNumFace[j] != NULL) free(regionArray[i].boundSubNumFace[j]);
				if(regionArray[i].boundSubCenterCoor[j] != NULL)	
					free(regionArray[i].boundSubCenterCoor[j]);
				if(regionArray[i].boundSubCoor[j] != NULL)	
					free(regionArray[i].boundSubCoor[j]);
				if(regionArray[i].boundVirtualNeighDir[j] != NULL)	
					free(regionArray[i].boundVirtualNeighDir[j]);
				if(regionArray[i].bNeedUpdate[j] != NULL) free(regionArray[i].bNeedUpdate[j]);
			}
			if(regionArray[i].neighID != NULL) free(regionArray[i].neighID);
			if(regionArray[i].boundSubNumFace != NULL) free(regionArray[i].boundSubNumFace);
			if(regionArray[i].boundSubCenterCoor != NULL) free(regionArray[i].boundSubCenterCoor);
			if(regionArray[i].boundSubCoor != NULL) free(regionArray[i].boundSubCoor);
			if(regionArray[i].boundVirtualNeighDir != NULL) free(regionArray[i].boundVirtualNeighDir);
			if(regionArray[i].bNeedUpdate != NULL) free(regionArray[i].bNeedUpdate);
			if(regionArray[i].numMolFromMicro != NULL) free(regionArray[i].numMolFromMicro);
		}
		
		if(regionArray[i].numSubRegionNeigh != NULL) free(regionArray[i].numSubRegionNeigh);
		if(regionArray[i].isRegionNeigh != NULL) free(regionArray[i].isRegionNeigh);
		
		if(regionArray[i].bDiffuse != NULL)
			free(regionArray[i].bDiffuse);
	}
}

// Find volume of specified region (exluding volume of children)
double findRegionVolume(const struct region regionArray[],
	const short curRegion,
	bool bOuter)
{
	short curChild, childID;
	double volume;

	switch (regionArray[curRegion].spec.type)
	{
		case REGION_NORMAL:
			// Default. Find volume as usual
			volume = boundaryVolume(regionArray[curRegion].spec.shape, regionArray[curRegion].boundary);
			break;
		case REGION_SURFACE_3D:
			// Region is a 3D surface. If shape is 3D, then find surface area
			if(regionArray[curRegion].plane == PLANE_3D)
				volume = boundarySurfaceArea(regionArray[curRegion].spec.shape, regionArray[curRegion].boundary);
			else
				volume = boundaryVolume(regionArray[curRegion].spec.shape, regionArray[curRegion].boundary);
			break;
		case REGION_SURFACE_2D:
			// Region is a 2D surface. Find perimeter of 2D shape
			volume = boundarySurfaceArea(regionArray[curRegion].spec.shape, regionArray[curRegion].boundary);
			break;
		default:
			// We should not have gotten here
			fprintf(stderr, "ERROR: Region %u (label: \"%s\") has unknown type %d.\n",
				curRegion, regionArray[curRegion].spec.label, regionArray[curRegion].spec.type);
	}
	
	if(!bOuter)
	{
		// Subtract outer volumes of children
		for(curChild = 0; curChild < regionArray[curRegion].numChildren; curChild++)
		{
			childID = regionArray[curRegion].childrenID[curChild];
			switch(regionArray[curRegion].spec.type)
			{
				case REGION_NORMAL:
					// Only find volume if child is of same dimension
					if(regionArray[curRegion].plane == regionArray[childID].plane)
					{
						volume -= boundaryVolume(regionArray[childID].spec.shape,
							regionArray[childID].boundary);
					}
					break;
				case REGION_SURFACE_3D:
					// Subtract surface area if child is 2D
					if(regionArray[childID].plane != PLANE_3D)
					{
						volume -= boundaryVolume(regionArray[childID].spec.shape,
							regionArray[childID].boundary);
					}
					break;
				case REGION_SURFACE_2D:
					// 2D surfaces cannot have children
					break;
			}
		}
	}
	
	return volume;
}
	
// Count the cumulative number of subvolumes defined by subvol_spec
uint32_t countAllSubvolumes(const struct region regionArray[],
	const short NUM_REGIONS)
{
	short i; // Current Region
	uint32_t num_subvolumes = 0;
	for(i = 0; i < NUM_REGIONS; i++)
		num_subvolumes += (uint32_t) regionArray[i].numSub;
	return num_subvolumes;
}

// Determine which regions of subvolumes are adjacent
void findRegionTouch(const short NUM_REGIONS,
	const struct spec_region3D subvol_spec[],
	struct region regionArray[],
	const double SUBVOL_BASE_SIZE)
{
	short int i,j,k; // Current Region
	short int surfRegion, normalRegion;
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
				if(regionArray[j].spec.surfaceType == SURFACE_OUTER
					|| regionArray[i].spec.surfaceType == SURFACE_INNER)
					continue; // Regions are not actually neighbors
				regionArray[i].isRegionNeigh[j] = true;
				regionArray[i].regionNeighDir[j] = PARENT;
				regionArray[i].numRegionNeigh++;		
				continue;
			}
			
			if(regionArray[j].parentID == i)
			{
				// Region i is region j's parent
				if(regionArray[j].spec.surfaceType == SURFACE_INNER
					|| regionArray[i].spec.surfaceType == SURFACE_OUTER)
					continue; // Regions are not actually neighbors
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
			
			if(regionArray[i].plane != regionArray[j].plane)
			{ // Regions can only touch if one is normal 3D and other is 3D surface
				if(!(regionArray[i].plane == PLANE_3D
					&& regionArray[j].spec.type == REGION_SURFACE_3D)
					&& !(regionArray[i].spec.type == REGION_SURFACE_3D
					&& regionArray[j].plane == PLANE_3D))
					continue;
			}
			
			if(bBoundaryAdjacent(regionArray[i].spec.shape, regionArray[i].boundary,
				regionArray[j].spec.shape, regionArray[j].boundary, error, &direction))
			{	// Regions share a geometric face
				// Check for surfaces
				if((regionArray[i].spec.type == REGION_NORMAL
					&& regionArray[j].spec.type != REGION_NORMAL)
					|| (regionArray[i].spec.type != REGION_NORMAL
					&& regionArray[j].spec.type == REGION_NORMAL))
				{
					if(regionArray[i].spec.type == REGION_NORMAL)
					{
						normalRegion = i;
						surfRegion = j;
					} else
					{
						normalRegion = j;
						surfRegion = i;
					}
					switch (regionArray[surfRegion].spec.surfaceType)
					{
						case SURFACE_INNER:
							if((surfRegion == i
								&& (direction == IN || direction == LEFT || direction == DOWN))
								|| (surfRegion == j
								&& (direction == OUT || direction == RIGHT || direction == UP)))
								continue;
							break;
						case SURFACE_OUTER:
							if((surfRegion == j
								&& (direction == IN || direction == LEFT || direction == DOWN))
								|| (surfRegion == i
								&& (direction == OUT || direction == RIGHT || direction == UP)))
								continue;
							break;
					}
				}
				regionArray[i].boundaryRegion[direction] = true;
				regionArray[i].boundaryRegionMicro[direction] = regionArray[j].spec.bMicro;
				regionArray[i].isRegionNeigh[j] = true;
				regionArray[i].regionNeighDir[j] = direction;
				regionArray[i].numRegionNeigh++;
			}
		}
	}

	// Make a pass to check for surface regions being in the same regime as their neighbors
	for(i = 0; i < NUM_REGIONS; i++)
	{
		for(j = i+1; j < NUM_REGIONS; j++)
		{
			if(regionArray[i].isRegionNeigh[j]
				&& (regionArray[i].spec.type != REGION_NORMAL
				|| regionArray[j].spec.type != REGION_NORMAL))
			{ // Regions are neighbors and at least one of them is a surface
				if(regionArray[i].spec.bMicro != regionArray[j].spec.bMicro)
				{
					fprintf(stderr, "ERROR: Regions %u and %u (labels: \"%s\" and \"%s\") are neighbors in different regimes but at least one of them is a surface.\n", i, j, subvol_spec[i].label, subvol_spec[j].label);
					exit(EXIT_FAILURE);
				}
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
			fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\")'s microscopic parameters.\n", i, subvol_spec[i].label);
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
							case RECTANGLE:
								regionArray[i].numRegionNeighFace[j] = 6;
								break;
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
							case RECTANGLE:
								regionArray[i].numRegionNeighFace[j] = 6;
								break;
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
					fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\")'s neighbor parameters with region %u (label: \"%s\").\n", i, subvol_spec[i].label, j, subvol_spec[j].label);
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
uint32_t findSubInBoundaryList(const short curRegion,
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
			&& bPointInRegionNotChild(curRegion, regionArray, point, false))
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
// of specified neighbor region.
// Returned value is index of subvolume in
// regionArray[curRegion].neighID[neighRegion] array
uint32_t findNearestSub(const short curRegion,
	const struct region regionArray[],
	const short neighRegion,
	double x,
	double y,
	double z)
{
	// Number of subvolumes that border the specified neighbor region
	uint32_t numSub = regionArray[curRegion].numSubRegionNeigh[neighRegion];
	
	uint32_t curSub;
	uint32_t minSub = UINT32_MAX;
	double minDistSq = DBL_MAX;
	double curDistSq;
	for(curSub = 0; curSub < numSub; curSub++)
	{
		curDistSq = (x - regionArray[curRegion].boundSubCenterCoor[neighRegion][curSub][0])
			*(x - regionArray[curRegion].boundSubCenterCoor[neighRegion][curSub][0]);
		if (curDistSq > minDistSq)
			continue; // No need to check y-coordinate
		
		curDistSq += (y - regionArray[curRegion].boundSubCenterCoor[neighRegion][curSub][1])
			*(y - regionArray[curRegion].boundSubCenterCoor[neighRegion][curSub][1]);
		if (curDistSq > minDistSq)
			continue; // No need to check z-coordinate
		
		curDistSq += (z - regionArray[curRegion].boundSubCenterCoor[neighRegion][curSub][2])
			*(z - regionArray[curRegion].boundSubCenterCoor[neighRegion][curSub][2]);
		
		if(curDistSq < minDistSq)
		{ // This subvolume is the closest so far
			minSub = curSub;
			minDistSq = curDistSq;
		}
	}
	return minSub;
}

// Determine coordinates of child region as subvolumes of parent
// Applies to Rectangular regions only.
void findChildCoordinates(const short parentRegion,
	const short childRegion,
	const short curChild,
	const struct region regionArray[])
{
	regionArray[parentRegion].childrenCoor[curChild][0] = (uint32_t)
		round((regionArray[childRegion].boundary[0] - regionArray[parentRegion].boundary[0])/regionArray[parentRegion].actualSubSize);
	regionArray[parentRegion].childrenCoor[curChild][1] = (uint32_t)
		round((regionArray[childRegion].boundary[1] - regionArray[parentRegion].boundary[0])/regionArray[parentRegion].actualSubSize - 1);
	regionArray[parentRegion].childrenCoor[curChild][2] = (uint32_t)
		round((regionArray[childRegion].boundary[2] - regionArray[parentRegion].boundary[2])/regionArray[parentRegion].actualSubSize);
	regionArray[parentRegion].childrenCoor[curChild][3] = (uint32_t)
		round((regionArray[childRegion].boundary[3] - regionArray[parentRegion].boundary[2])/regionArray[parentRegion].actualSubSize - 1);
	regionArray[parentRegion].childrenCoor[curChild][4] = (uint32_t)
		round((regionArray[childRegion].boundary[4] - regionArray[parentRegion].boundary[4])/regionArray[parentRegion].actualSubSize);
	regionArray[parentRegion].childrenCoor[curChild][5] = (uint32_t)
		round((regionArray[childRegion].boundary[5] - regionArray[parentRegion].boundary[4])/regionArray[parentRegion].actualSubSize - 1);
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
	
	bool bArea = regionArray[curRegion].effectiveDim != regionArray[curRegion].dimension;
	
	if(bArea && bBoundarySurround(boundary2Type, boundary2, 
		regionArray[curRegion].spec.shape, regionArray[curRegion].boundary, 0.)
		&& !bBoundarySurround(regionArray[curRegion].spec.shape, regionArray[curRegion].boundary,
		boundary2Type, boundary2, 0.))
	{ // Boundary is fully inside region such that there is no intersection with region boundary
		return 0.;		
	}
	
	if(boundary2Type == POINT)
	{
		if(bPointInRegionNotChild(curRegion, regionArray,
			boundary2, false))
			return INFINITY;
		else
			return 0.;
	}
	
	// Calculate overall intersection area (including space occupied by children)
	intersectType = intersectBoundary(regionArray[curRegion].spec.shape,
		regionArray[curRegion].boundary,	boundary2Type, boundary2, interBoundary);
	if(intersectType == UNDEFINED_SHAPE)
	{
		// Invalid intersection of boundary with region		
		fprintf(stderr, "ERROR: An intersection with region %u (label: \"%s\") is invalid.\nShapes are %s and %s",
			curRegion, regionArray[curRegion].spec.label, boundaryString(regionArray[curRegion].spec.shape),
			boundaryString(boundary2Type));
		exit(EXIT_FAILURE);			
	}
	if(regionArray[curRegion].plane != PLANE_3D
		&& intersectType == RECTANGULAR_BOX)
	{
		intersectType = RECTANGLE;
	}
	if(bArea)
	{
		volume = boundarySurfaceArea(intersectType, interBoundary);
	} else
	{
		volume = boundaryVolume(intersectType, interBoundary);
	}
	
	if(volume <= 0)
		return 0.;	// Region does not intersect at all
	
	// Surface regions whose shape is the same dimension that they are a surface
	// for cannot have children that are also surfaces so do not consider their
	// children.
	if(bArea)
	{
		return volume;
	}
	
	// Subtract area populated by children
	for(curChild = 0; curChild < regionArray[curRegion].numChildren; curChild++)
	{
		childRegion = regionArray[curRegion].childrenID[curChild];
		
		// Only consider children that are of the same plane
		// Note that regions whose effective dimension is different from their shape
		// dimension are not considered in this for loop because bArea == true
		if(regionArray[curRegion].plane != regionArray[childRegion].plane)
			continue;
		
		intersectType = intersectBoundary(regionArray[childRegion].spec.shape,
			regionArray[childRegion].boundary, boundary2Type, boundary2, interBoundary);
				
		if(intersectType == UNDEFINED_SHAPE)
		{
			// Invalid intersection of boundary with child		
			fprintf(stderr, "ERROR: An intersection with region %u (label: \"%s\") is invalid.\n", childRegion, regionArray[childRegion].spec.label);	
			fprintf(stderr, "Shapes are %s and %s",
				boundaryString(regionArray[childRegion].spec.shape),
				boundaryString(boundary2Type));
			exit(EXIT_FAILURE);			
		}
		
		if(regionArray[childRegion].plane != PLANE_3D
			&& intersectType == RECTANGULAR_BOX)
		{ // We still want to measure area of 2D region within 3D boundary
			intersectType = regionArray[childRegion].spec.shape;
		}		
		childVolume = boundaryVolume(intersectType, interBoundary);
		if(childVolume > 0.)
			volume -= childVolume; // Only subtract positive child volumes

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
		case RECTANGLE:
		case RECTANGULAR_BOX:
			for(curPlane = 0;
				curPlane < regionArray[startRegion].numRegionNeighFace[endRegion];
				curPlane++)
			{
				if(bLineHitInfinitePlane(p1, L, length, RECTANGULAR_BOX,
					regionArray[startRegion].boundRegionFaceCoor[endRegion][curPlane],
					regionArray[startRegion].regionNeighFaceDir[endRegion][curPlane],
					false, d, intersectPoint, true)
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
				regionArray[startRegion].regionNeighFaceDir[endRegion][0], bInside,
				d, intersectPoint, true);
		default:
			fprintf(stderr,"ERROR: Boundary type %d invalid for boundary intersection between regions %u and %u.\n", boundary1Type, startRegion, endRegion);
			return false;		
	}
}

// Is point in a region and not in any of its children?
// Includes a check to ignore surface regions
bool bPointInRegionNotChild(const short curRegion,
	const struct region regionArray[],
	const double point[3],
	bool bIgnoreSurfaceChildren)
{
	short curChild;
	
	if(bPointInBoundary(point, regionArray[curRegion].spec.shape,
		regionArray[curRegion].boundary))
	{ // Point is within the region's outer boundary
		for(curChild = 0; curChild < regionArray[curRegion].numChildren; curChild++)
		{
			if(bIgnoreSurfaceChildren &&
				regionArray[regionArray[curRegion].childrenID[curChild]].spec.type !=
				REGION_NORMAL)
				continue; // Child is a surface and we are ignoring surfaces
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
	short * actualRegion,
	bool bSurfaceOnly)
{
	short curChild;
	
	if(bPointInBoundary(point, regionArray[curRegion].spec.shape,
		regionArray[curRegion].boundary))
	{ // Point is within the region's outer boundary
		for(curChild = 0; curChild < regionArray[curRegion].numChildren; curChild++)
		{
			if(bSurfaceOnly && regionArray[regionArray[curRegion].childrenID[curChild]].spec.type == REGION_NORMAL)
				continue;
			if(bPointInRegionOrChild(regionArray[curRegion].childrenID[curChild],
				regionArray, point, actualRegion, bSurfaceOnly))
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
		case RECTANGLE:
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

// Is the outer boundary of a child region flush with the subvolumes of its
// parent? Applies to rectangular regions (2D or 3D)
bool bChildRegionNotFlush(const short parentRegion,
	const short childRegion,
	const struct region regionArray[],
	const double boundAdjError)
{
	return fabs((regionArray[childRegion].boundary[0] - regionArray[parentRegion].boundary[0])/regionArray[parentRegion].actualSubSize - round((regionArray[childRegion].boundary[0] - regionArray[parentRegion].boundary[0])/regionArray[parentRegion].actualSubSize)) > boundAdjError ||
		fabs((regionArray[childRegion].boundary[1] - regionArray[parentRegion].boundary[1])/regionArray[parentRegion].actualSubSize - round((regionArray[childRegion].boundary[1] - regionArray[parentRegion].boundary[1])/regionArray[parentRegion].actualSubSize)) > boundAdjError ||
		fabs((regionArray[childRegion].boundary[2] - regionArray[parentRegion].boundary[2])/regionArray[parentRegion].actualSubSize - round((regionArray[childRegion].boundary[2] - regionArray[parentRegion].boundary[2])/regionArray[parentRegion].actualSubSize)) > boundAdjError ||
		fabs((regionArray[childRegion].boundary[3] - regionArray[parentRegion].boundary[3])/regionArray[parentRegion].actualSubSize - round((regionArray[childRegion].boundary[3] - regionArray[parentRegion].boundary[3])/regionArray[parentRegion].actualSubSize)) > boundAdjError ||
		fabs((regionArray[childRegion].boundary[4] - regionArray[parentRegion].boundary[4])/regionArray[parentRegion].actualSubSize - round((regionArray[childRegion].boundary[4] - regionArray[parentRegion].boundary[4])/regionArray[parentRegion].actualSubSize)) > boundAdjError ||
		fabs((regionArray[childRegion].boundary[5] - regionArray[parentRegion].boundary[5])/regionArray[parentRegion].actualSubSize - round((regionArray[childRegion].boundary[5] - regionArray[parentRegion].boundary[5])/regionArray[parentRegion].actualSubSize)) > boundAdjError;
}

// Lock point coordinate to chosen region face
void lockPointToRegion(double point[3],
	const short startRegion,
	const short endRegion,
	const struct region regionArray[],
	const short faceID)
{
	short boundRegion;
	short lockInd;
	double coorSq[3];
	double rSq;
	if (regionArray[startRegion].isRegionNeigh[endRegion]
		&& regionArray[startRegion].regionNeighDir[endRegion] == PARENT)
		boundRegion = startRegion;
	else
		boundRegion = endRegion;
	
	switch(regionArray[boundRegion].spec.shape)
	{
		case RECTANGLE:
		case RECTANGULAR_BOX:
			if(faceID < 2)
				point[0] = regionArray[boundRegion].boundary[faceID];
			else if (faceID > 3)
				point[2] = regionArray[boundRegion].boundary[faceID];
			else					
				point[1] = regionArray[boundRegion].boundary[faceID];	
			break;
		case SPHERE:
			// Arbitrarily adjust "furthest" coordinate to lock point to sphere surface
			coorSq[0] = squareDBL(regionArray[boundRegion].boundary[0] - point[0]);
			coorSq[1] = squareDBL(regionArray[boundRegion].boundary[1] - point[1]);
			coorSq[2] = squareDBL(regionArray[boundRegion].boundary[2] - point[2]);
			rSq = coorSq[0] + coorSq[1] + coorSq[2];
			if(coorSq[0] >= coorSq[1] && coorSq[0] >= coorSq[2])
				lockInd = 0; // x-coordinate is furthest from center
			else if (coorSq[1] >= coorSq[0] && coorSq[1] >= coorSq[2])
				lockInd = 1; // y-coordinate is furthest from center
			else
				lockInd = 2; // z-coordinate is furthest from center
			
			if(point[lockInd] > regionArray[boundRegion].boundary[lockInd])
				point[lockInd] = sqrt(squareDBL(regionArray[boundRegion].boundary[3]) -
					rSq + coorSq[lockInd]);
			else
				point[lockInd] = -sqrt(squareDBL(regionArray[boundRegion].boundary[3]) -
					rSq + coorSq[lockInd]);
			point[lockInd] += regionArray[boundRegion].boundary[lockInd];
			break;
		default:
			fprintf(stderr,"ERROR: Point cannot be locked to region %u (label: \"%s\") because it is of type %d.\n", boundRegion, regionArray[boundRegion].spec.label, regionArray[boundRegion].spec.shape);
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
			regionArray[curRegion].boundary, regionArray[curRegion].dimension ==
			regionArray[curRegion].effectiveDim, regionArray[curRegion].plane);
		if(bPointInRegionNotChild(curRegion, regionArray, point, false))
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
		if(bPointInRegionNotChild(curRegion, regionArray, point, false))
			return curRegion;
	}
	
	// Point was not found in any region
	return SHRT_MAX;
}

// Does a subvolume face a region? If yes, then along which faces?
// Assert that current subvolume is along its own region boundary, and that
// neighbor region is microscopic
bool bSubFaceRegion(struct region regionArray[],
	const short neighRegion,
	const double curSubBound[6],
	double boundAdjError,
	unsigned short * numFace,
	unsigned short dirArray[6])
{
	double neighPoint[3];
	
	*numFace = 0;
	
	// Initialize neighbor point to center of current subvolume
	neighPoint[0] = (curSubBound[0] + curSubBound[1])/2;
	neighPoint[1] = (curSubBound[2] + curSubBound[3])/2;
	neighPoint[2] = (curSubBound[4] + curSubBound[5])/2;
	
	// Check to see if point just to "left" is within neighbor
	neighPoint[0] = curSubBound[0] - boundAdjError;
	if(bPointInRegionNotChild(neighRegion, regionArray, neighPoint, false))
	{
		dirArray[(*numFace)++] = LEFT;
	}
	
	// Check to see if point just to "right" is within neighbor
	neighPoint[0] = curSubBound[1] + boundAdjError;
	if(bPointInRegionNotChild(neighRegion, regionArray, neighPoint, false))
	{
		dirArray[(*numFace)++] = RIGHT;
	}
	neighPoint[0] = (curSubBound[0] + curSubBound[1])/2;
	
	// Check to see if point just to "down" is within neighbor
	neighPoint[1] = curSubBound[2] - boundAdjError;
	if(bPointInRegionNotChild(neighRegion, regionArray, neighPoint, false))
	{
		dirArray[(*numFace)++] = DOWN;
	}
	
	// Check to see if point just to "up" is within neighbor
	neighPoint[1] = curSubBound[3] + boundAdjError;
	if(bPointInRegionNotChild(neighRegion, regionArray, neighPoint, false))
	{
		dirArray[(*numFace)++] = UP;
	}
	neighPoint[1] = (curSubBound[2] + curSubBound[3])/2;
	
	// Check to see if point just to "in" is within neighbor
	neighPoint[2] = curSubBound[4] - boundAdjError;
	if(bPointInRegionNotChild(neighRegion, regionArray, neighPoint, false))
	{
		dirArray[(*numFace)++] = IN;
	}
	
	// Check to see if point just to "out" is within neighbor
	neighPoint[2] = curSubBound[5] + boundAdjError;
	if(bPointInRegionNotChild(neighRegion, regionArray, neighPoint, false))
	{
		dirArray[(*numFace)++] = OUT;
	}
	
	return *numFace > 0;
}

// Initialize the region nesting (i.e., determine each region's parent and
// children regions, if applicable)
void initializeRegionNesting(const short NUM_REGIONS,
	struct region regionArray[],
	bool * bFail)
{
	short i, j, curChild;
	short curChildArray[NUM_REGIONS];	
	bool bFailRectangleChild; // Fail switch for rectangle child on boundary of box
	
	// Determine region nesting (i.e., find each region's parent if it has one)
	for(i = 0; i < NUM_REGIONS; i++)
	{	// Assign parents and count number of children in each region
		regionArray[i].bParent = false;
		regionArray[i].parentID = SHRT_MAX;
		curChildArray[i] = 0;
		if(regionArray[i].spec.parent && strlen(regionArray[i].spec.parent) > 0)
		{	// Region has a parent. Determine index of parent
			for(j = 0; j < NUM_REGIONS; j++)
			{
				if(!strcmp(regionArray[i].spec.parent, regionArray[j].spec.label)
					&& i != j)
				{
					// Region j is the parent
					regionArray[i].bParent = true;
					regionArray[i].parentID = j;
					regionArray[j].numChildren++;
					break;
				}
			}
			if(!regionArray[i].bParent)
			{
				// Parent was not found. Throw error.
				fprintf(stderr, "ERROR: Region %u (Label: \"%s\") has non-existent parent \"%s\".\n",
					i, regionArray[i].spec.label, regionArray[j].spec.parent);
				*bFail = true;
			}
		}
	}
	
	// Allocate memory for each region's children
	for(i = 0; i < NUM_REGIONS; i++)
	{
		if(regionArray[i].numChildren > 0)
		{
			// Region has children. Allocate children array
			regionArray[i].childrenID = malloc(regionArray[i].numChildren*sizeof(short));
			regionArray[i].childrenCoor = malloc(regionArray[i].numChildren*sizeof(uint32_t [6]));
			if(regionArray[i].childrenID == NULL || regionArray[i].childrenCoor == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation to region %u (label: \"%s\") parameters.\n", i, regionArray[i].spec.label);
				exit(EXIT_FAILURE);
			}
		}
	}
	
	// Store IDs of each region's children 
	for(i = 0; i < NUM_REGIONS; i++)
	{
		if(regionArray[i].bParent)
		{
			j = regionArray[i].parentID;
			regionArray[j].childrenID[curChildArray[j]++] = i;
		}
	}
	
	// Store relative coordinates of each region's children 
	for(i = 0; i < NUM_REGIONS; i++)
	{
		for(curChild = 0; curChild < regionArray[i].numChildren; curChild++)
		{
			j = regionArray[i].childrenID[curChild];
			
			if(regionArray[i].spec.shape == RECTANGULAR_BOX &&
				regionArray[j].spec.shape == RECTANGULAR_BOX)
			{				
				// Parent and child are boxes
				// Check for flush subvolumes first
				if(bChildRegionNotFlush(i, j, regionArray, regionArray[i].subResolution))
				{
					// Child is not flush with region's subvolumes
					fprintf(stderr, "ERROR: Nested region %u (label: \"%s\") is not properly aligned within parent region %u (label: \"%s\").\n", j, regionArray[j].spec.label, i, regionArray[i].spec.label);
					fprintf(stderr, "Both regions are rectangular and the outer boundary of the nested region must be flush with subvolumes of the parent region.\n");
					*bFail = true;
					continue;
				}
				
				// If child is an inner surface, then parent must an outer surface that covers
				// it fully
				if(regionArray[j].spec.surfaceType == SURFACE_INNER
					&& (regionArray[i].spec.surfaceType != SURFACE_OUTER
					|| !bBoundarySurround(RECTANGULAR_BOX, regionArray[i].boundary,
					RECTANGULAR_BOX, regionArray[j].boundary, -regionArray[i].subResolution)))
				{
					fprintf(stderr, "ERROR: Nested region %u (label: \"%s\") is an inner-pointing surface with parent region %u (label: \"%s\").\n",
						j, regionArray[j].spec.label, i, regionArray[i].spec.label);
					fprintf(stderr, "A 3D child can only be an inner surface if its parent is an outer surface.\n");
					*bFail = true;
					continue;
				}
				
				// Confirm that child is inside of parent
				if(!bBoundarySurround(RECTANGULAR_BOX, regionArray[j].boundary,
					RECTANGULAR_BOX, regionArray[i].boundary, -regionArray[i].subResolution))
				{
					// Child is sticking out of parent
					fprintf(stderr, "ERROR: Outer boundary of nested region %u (label: \"%s\") is not properly surrounded by parent region %u (label: \"%s\").\n", j, regionArray[j].spec.label, i, regionArray[i].spec.label);
					fprintf(stderr, "Both regions are rectangular.\n");
					*bFail = true;
					continue;
				}
				
				// Determine coordinates of child relative to parent subvolumes
				findChildCoordinates(i, j, curChild, regionArray);
				
			} else if(regionArray[i].spec.shape == RECTANGLE &&
				regionArray[j].spec.shape == RECTANGLE)
			{
				// Check for flush subvolumes				
				if(bChildRegionNotFlush(i, j, regionArray, regionArray[i].subResolution))
				{
					// Child is not flush with region's subvolumes
					fprintf(stderr, "ERROR: Nested region %u (label: \"%s\") is not properly aligned within parent region %u (label: \"%s\").\n", j, regionArray[j].spec.label, i, regionArray[i].spec.label);
					fprintf(stderr, "Both regions are rectangles and the outer boundary of the nested region must be flush with subvolumes of the parent region.\n");
					*bFail = true;
				}
				
				// Confirm that child is inside of parent
				if(regionArray[i].plane != regionArray[j].plane
					|| !bBoundarySurround(RECTANGLE, regionArray[j].boundary,
					RECTANGLE, regionArray[i].boundary, -regionArray[i].subResolution))
				{
					// Child is sticking out of parent
					fprintf(stderr, "ERROR: Outer boundary of nested region %u (label: \"%s\") is not properly surrounded by parent region %u (label: \"%s\").\n", j, regionArray[j].spec.label, i, regionArray[i].spec.label);
					fprintf(stderr, "Both regions are rectangles.\n");
					*bFail = true;
				}
				
				// Determine coordinates of child relative to parent subvolumes
				findChildCoordinates(i, j, curChild, regionArray);
				if(regionArray[i].plane == PLANE_XY)
				{
					regionArray[i].childrenCoor[curChild][4] = 0;
					regionArray[i].childrenCoor[curChild][5] = 0;
				} else if (regionArray[i].plane == PLANE_XZ)
				{
					regionArray[i].childrenCoor[curChild][2] = 0;
					regionArray[i].childrenCoor[curChild][3] = 0;
				} else if (regionArray[i].plane == PLANE_YZ)
				{
					regionArray[i].childrenCoor[curChild][0] = 0;
					regionArray[i].childrenCoor[curChild][1] = 0;
				}
			} else if(regionArray[i].spec.shape == SPHERE &&
				regionArray[j].spec.shape == RECTANGULAR_BOX)
			{
				if(!bBoundarySurround(regionArray[j].spec.shape, regionArray[j].boundary, SPHERE, regionArray[i].boundary, regionArray[j].actualSubSize - regionArray[i].subResolution))
				{
					// Child is sticking out of parent
					fprintf(stderr, "ERROR: Outer boundary of nested region %u (label: \"%s\") is not properly surrounded by parent region %u (label: \"%s\").\n", j, regionArray[j].spec.label, i, regionArray[i].spec.label);
					fprintf(stderr, "Rectangular region %u should be inside of spherical region %u and there must be a clearance between the surfaces of at least the subvolume size of region %u which is %.2em.\n", j, i, j, regionArray[j].actualSubSize);
					*bFail = true;
				}
				
			} else if(regionArray[i].spec.shape == RECTANGULAR_BOX &&
				regionArray[j].spec.shape == SPHERE)
			{
				if(!regionArray[i].spec.bMicro)
				{
					// A mesoscopic parent cannot have a spherical child
					fprintf(stderr, "ERROR: Normal mesoscopic region %u (label: \"%s\") has a spherical child region %u (label: \"%s\").\n", i, regionArray[i].spec.label, j, regionArray[j].spec.label);
					*bFail = true;
				}
				
				if(!bBoundarySurround(SPHERE, regionArray[j].boundary, RECTANGULAR_BOX, regionArray[i].boundary, regionArray[i].actualSubSize - regionArray[i].subResolution))
				{
					// Child is sticking out of parent
					fprintf(stderr, "ERROR: Outer boundary of nested region %u (label: \"%s\") is not properly surrounded by parent region %u (label: \"%s\").\n", j, regionArray[j].spec.label, i, regionArray[i].spec.label);
					fprintf(stderr, "Spherical region %u should be inside of rectangular region %u and there must be a clearance between the surfaces of at least the subvolume size of region %u which is %.2em.\n", j, i, i, regionArray[i].actualSubSize);
					*bFail = true;
				}
			} else if(regionArray[i].spec.shape == SPHERE &&
				regionArray[j].spec.shape == SPHERE)
			{
				if(!bBoundarySurround(SPHERE, regionArray[j].boundary, SPHERE, regionArray[i].boundary, 0.))
				{
					// Child is sticking out of parent
					fprintf(stderr, "ERROR: Outer boundary of spherical nested region %u (label: \"%s\") is not properly surrounded by spherical parent region %u (label: \"%s\").\n", j, regionArray[j].spec.label, i, regionArray[i].spec.label);
					*bFail = true;
				}
			} else if(regionArray[i].spec.shape == RECTANGULAR_BOX &&
				regionArray[j].spec.shape == RECTANGLE)
			{ // A box region can have a rectangle child only if both regions are surfaces
			
				// NOTE: Short-circuiting this option for now!
				// Combination of child and parent is invalid
				fprintf(stderr, "ERROR: Combination of child/parent for region %u (label: \"%s\") nested in region %u (label: \"%s\") is invalid.\n", j, regionArray[j].spec.label, i, regionArray[i].spec.label);
				fprintf(stderr, "Shapes are %s and %s, respectively.\n",
					boundaryString(regionArray[j].spec.shape),
					boundaryString(regionArray[i].spec.shape));
				*bFail = true;
				continue;
			
				if(regionArray[i].spec.type != REGION_SURFACE_3D
					|| regionArray[j].spec.type != REGION_SURFACE_3D)
				{
					// Child is sticking out of parent
					fprintf(stderr, "ERROR: Region %u (label: \"%s\") has improper child region %u (label: \"%s\").\n", i, regionArray[i].spec.label, j, regionArray[j].spec.label);
					fprintf(stderr, "A rectangle can only be the child of a box region if both regions are surfaces.\n");
					*bFail = true;
				}
				
				bFailRectangleChild = false;
				
				// Determine coordinates of child relative to parent subvolumes
				findChildCoordinates(i, j, curChild, regionArray);
				
				switch(regionArray[j].plane)
				{
					case PLANE_XY:
						bFailRectangleChild = !(regionArray[i].childrenCoor[curChild][4]
							|| regionArray[i].childrenCoor[curChild][5]);
						break;
					case PLANE_XZ:
						bFailRectangleChild = !(regionArray[i].childrenCoor[curChild][2]
							|| regionArray[i].childrenCoor[curChild][3]);
						break;
					case PLANE_YZ:
						bFailRectangleChild = !(regionArray[i].childrenCoor[curChild][0]
							|| regionArray[i].childrenCoor[curChild][1]);
						break;
				}
				if(bFailRectangleChild)
				{
					fprintf(stderr, "ERROR: Nested region %u (label: \"%s\") is not on boundary of parent region %u (label: \"%s\").\n", j, regionArray[j].spec.label, i, regionArray[i].spec.label);
					fprintf(stderr, "A rectangular surface region can only have a parent surface if it is directly on the surface of the parent.\n");
					*bFail = true;
					continue;
				}
				
				// Confirm that child is inside of parent
				if(!bBoundarySurround(RECTANGLE, regionArray[j].boundary,
					RECTANGULAR_BOX, regionArray[i].boundary, -regionArray[i].subResolution))
				{
					// Child is sticking out of parent
					fprintf(stderr, "ERROR: Outer boundary of nested region %u (label: \"%s\") is not properly surrounded by parent region %u (label: \"%s\").\n", j, regionArray[j].spec.label, i, regionArray[i].spec.label);
					fprintf(stderr, "Both regions are rectangular.\n");
					*bFail = true;
				}
			} else
			{
				// Combination of child and parent is invalid
				fprintf(stderr, "ERROR: Combination of child/parent for region %u (label: \"%s\") nested in region %u (label: \"%s\") is invalid.\n", j, regionArray[j].spec.label, i, regionArray[i].spec.label);
				fprintf(stderr, "Shapes are %s and %s, respectively.\n",
					boundaryString(regionArray[j].spec.shape),
					boundaryString(regionArray[i].spec.shape));
				*bFail = true;
			}
		}
		
		// Calculate region volume now that we know the validity of child regions
		regionArray[i].volume = findRegionVolume(regionArray, i, false);
	}
}

// Determine number of subvolumes in each region
void findNumRegionSubvolumes(const short NUM_REGIONS,
	struct region regionArray[])
{
	short i, j, curChild;
	uint32_t curID = 0; // Subvolume index
	uint32_t curX, curY, curZ; // Coordinates of prospective subvolume
	uint32_t subCoorInd[3]; // (curX, curY, curZ)
	double curSubBound[6]; // Boundary of prospective subvolume
	
	for(i = 0; i < NUM_REGIONS; i++)
	{		
		// Determine number of subvolumes
		regionArray[i].firstID = curID;
		if(regionArray[i].spec.shape == SPHERE)
			regionArray[i].numSub = 1UL; // Spherical regions always have one subvolume
		else
		{
			if(regionArray[i].plane == PLANE_3D)
			{
				if(regionArray[i].spec.type == REGION_NORMAL)
				{
					regionArray[i].numSub = (uint32_t) regionArray[i].spec.numX *
						regionArray[i].spec.numY * regionArray[i].spec.numZ;
				} else
				{
					regionArray[i].numSub = (uint32_t) 2*regionArray[i].spec.numX *
						regionArray[i].spec.numY + 2*regionArray[i].spec.numY * regionArray[i].spec.numZ + 2*regionArray[i].spec.numX *
						regionArray[i].spec.numZ;
				}
			} else
			{
				switch(regionArray[i].plane)
				{
					case PLANE_XY:
						curX = regionArray[i].spec.numX;
						curY = regionArray[i].spec.numY;
						break;
					case PLANE_XZ:
						curX = regionArray[i].spec.numX;
						curY = regionArray[i].spec.numZ;
						break;
					case PLANE_YZ:
						curX = regionArray[i].spec.numY;
						curY = regionArray[i].spec.numZ;
						break;
				}
				if(regionArray[i].spec.type != REGION_SURFACE_2D)
				{		
					regionArray[i].numSub = (uint32_t) curX * curY;
				} else
				{
					regionArray[i].numSub = (uint32_t) 2*curX + 2*curY;
				}
			}			
			
			// Subtract subvolumes lost due to each child
			for(curChild = 0; curChild < regionArray[i].numChildren; curChild++)
			{
				j = regionArray[i].childrenID[curChild];
				if(regionArray[j].spec.shape == RECTANGULAR_BOX
					&& regionArray[i].spec.type == REGION_NORMAL)
				{	// A 3D child must have a 3D parent, and the subvolumes of
					// the parent are only affected if the parent is a normal region
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
					for(curZ = 0; curZ < regionArray[i].spec.numZ; curZ++)
					{
						for(curY = 0; curY < regionArray[i].spec.numY; curY++)
						{
							for(curX = 0; curX < regionArray[i].spec.numX; curX++)
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
				} else if(regionArray[j].spec.shape == RECTANGLE
					&& (regionArray[i].spec.shape == RECTANGLE
					|| regionArray[i].spec.type != REGION_NORMAL))
				{ // Only subtract area of children if parent is 2D
					switch(regionArray[j].plane)
					{
						case PLANE_XY:
							regionArray[i].numSub -= (uint32_t) regionArray[j].spec.numX *
								regionArray[j].spec.numY;
							break;				
						case PLANE_XZ:
							regionArray[i].numSub -= (uint32_t) regionArray[j].spec.numX *
								regionArray[j].spec.numZ;
							break;				
						case PLANE_YZ:
							regionArray[i].numSub -= (uint32_t) regionArray[j].spec.numY *
								regionArray[j].spec.numZ;
							break;				
						break;
					}
				}
			}
		}		
		curID += regionArray[i].numSub;
	}
}

// Allocate and initialize the flow parameters in each region
void initializeRegionFlow(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	struct region regionArray[],	
	const struct spec_region3D subvol_spec[])
{
	short i, j, m;
	unsigned short curMolType;
	
	for(i = 0; i < NUM_REGIONS; i++)
	{
		regionArray[i].bFlow = false;
		// Copy flow parameters from specification
		// (assignment of subvol_spec copies pointers but not allocated memory!)
		regionArray[i].spec.bFlow = malloc(NUM_MOL_TYPES*sizeof(bool));
		regionArray[i].spec.flowType = malloc(NUM_MOL_TYPES*sizeof(unsigned short));
		regionArray[i].spec.flowVector = malloc(NUM_MOL_TYPES*sizeof(double *));
		if(regionArray[i].spec.bFlow == NULL ||
			regionArray[i].spec.flowType == NULL ||
			regionArray[i].spec.flowVector == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\") flow parameters.\n", i, subvol_spec[i].label);
			exit(EXIT_FAILURE);
		}
		for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
		{
			regionArray[i].spec.bFlow[curMolType] = subvol_spec[i].bFlow[curMolType];
			regionArray[i].spec.flowType[curMolType] = subvol_spec[i].flowType[curMolType];
			switch(regionArray[i].spec.flowType[curMolType])
			{
				case FLOW_NONE:
					j = 0;
					break;
				case FLOW_UNIFORM:
					j = 3;
					break;
			}
			if(j > 0)
			{
				regionArray[i].bFlow = true;
				regionArray[i].spec.flowVector[curMolType] = malloc(j*sizeof(double));
				if(regionArray[i].spec.flowVector[curMolType] == NULL)
				{
					fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\") flow vector for molecule type %d.\n", i, subvol_spec[i].label, curMolType);
					exit(EXIT_FAILURE);
				}				
				switch(regionArray[i].spec.flowType[curMolType])
				{
					case FLOW_UNIFORM:
						for(m = 0; m < 3; m++)
							regionArray[i].spec.flowVector[curMolType][m] =
								subvol_spec[i].flowVector[curMolType][m];
				}
			} else
				regionArray[i].spec.flowVector[curMolType] = NULL;
		}		
		
		// Calculate constant flow parameters for microscopic regions
		if(regionArray[i].spec.bMicro)
		{
			regionArray[i].flowConstant = malloc(NUM_MOL_TYPES*sizeof(double *));
			if(regionArray[i].flowConstant == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\") constant flow parameters.\n", i, subvol_spec[i].label);
				exit(EXIT_FAILURE);
			}
			for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
			{
				switch(regionArray[i].spec.flowType[curMolType])
				{
					case FLOW_NONE:
						j = 0;
						break;
					case FLOW_UNIFORM:
						j = 3;
						break;
				}
				if(j > 0)
				{
					regionArray[i].flowConstant[curMolType] = malloc(j*sizeof(double));
					if(regionArray[i].flowConstant[curMolType] == NULL)
					{
						fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\") constant flow parameters for molecule type %d.\n", i, subvol_spec[i].label, curMolType);
						exit(EXIT_FAILURE);
					}
					switch(regionArray[i].spec.flowType[curMolType])
					{
						case FLOW_UNIFORM:
							for(m = 0; m < 3; m++)
								regionArray[i].flowConstant[curMolType][m] =
									regionArray[i].spec.flowVector[curMolType][m]*regionArray[i].spec.dt;
					}
				} else
				{
					regionArray[i].flowConstant[curMolType] = NULL;
				}
			}			
		} else
			regionArray[i].flowConstant = NULL;
	}
}

// Allocate memory for each region's neighbors
void allocateRegionNeighbors(const short NUM_REGIONS,
	struct region regionArray[])
{
	short i, j;
	
	for(i = 0; i < NUM_REGIONS; i++)
	{		
		regionArray[i].isRegionNeigh = malloc(NUM_REGIONS*sizeof(bool));
		regionArray[i].regionNeighDir = malloc(NUM_REGIONS*sizeof(unsigned short));
		if(regionArray[i].isRegionNeigh == NULL || regionArray[i].regionNeighDir == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\") parameters.\n", i, regionArray[i].spec.label);
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
				fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\")'s microscopic parameters.\n", i, regionArray[i].spec.label);
				exit(EXIT_FAILURE);
			}
		}
		
		regionArray[i].numSubRegionNeigh = malloc(NUM_REGIONS*sizeof(uint32_t));
		if(regionArray[i].numSubRegionNeigh == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation for region %u (label: \"%s\") parameters.\n", i, regionArray[i].spec.label);
			exit(EXIT_FAILURE);
		}
		
		for(j = 0; j < NUM_REGIONS; j++)
		{ // Initialize isRegionNeigh to false
			regionArray[i].isRegionNeigh[j] = false;
			regionArray[i].regionNeighDir[j] = UNDEFINED;
		}
		
		regionArray[i].numRegionNeigh = 0; // Initialize # of neighboring regions to 0
	}
}

// Final check of region overlap and correct adjacency
void validateRegions(const short NUM_REGIONS,
	const struct region regionArray[],
	bool * bFail)
{
	short curRegion, curChild1, curChild2, curChild1ID, curChild2ID;
	short otherRegion, squareRegion, boxRegion;
	short neighID, curNeigh;
	int intersectType;
	double intersection[6];
	double minVol, curVol;
	bool bDoubleSide; // Does a rectangle have neighbors on each side?
	bool bTemp; // Have we found the first 3D neighbor of rectangle yet?
	unsigned short initDir, curDir;
	short dim;
	
	for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
	{		
		// First check overlap of children. Checking overlaps this way reduces
		// the number of comparisons that must be made
		for (curChild1 = 0; curChild1 < regionArray[curRegion].numChildren; curChild1++)
		{ // Compare with every higher-ID child
			curChild1ID = regionArray[curRegion].childrenID[curChild1];
			
			if(regionArray[curChild1ID].spec.shape == RECTANGLE)
				minVol = regionArray[curChild1ID].subResolution *
					regionArray[curChild1ID].subResolution;
			else
				minVol = regionArray[curChild1ID].subResolution *
					regionArray[curChild1ID].subResolution * regionArray[curChild1ID].subResolution;
			
			for(curChild2 = curChild1 + 1;
				curChild2 < regionArray[curRegion].numChildren; curChild2++)
			{
				curChild2ID = regionArray[curRegion].childrenID[curChild2];
				
				intersectType = intersectBoundary(regionArray[curChild1ID].spec.shape,
					regionArray[curChild1ID].boundary, regionArray[curChild2ID].spec.shape,
					regionArray[curChild2ID].boundary, intersection);
				
				if(intersectType == UNDEFINED_SHAPE ||
					boundaryVolume(intersectType, intersection) > minVol)
				{ // Intersection detected					
					fprintf(stderr, "ERROR: Nested region %u (label: \"%s\") overlaps nested region %u (label: \"%s\"). Both regions are children of region %u (label: \"%s\").\n",
						curChild1ID, regionArray[curChild1ID].spec.label, curChild2ID,
						regionArray[curChild2ID].spec.label,
						curRegion, regionArray[curRegion].spec.label);
					*bFail = true;
				}
			}
		}
		
		if(regionArray[curRegion].spec.shape == RECTANGLE)
			minVol = regionArray[curRegion].subResolution *
				regionArray[curRegion].subResolution;
		else
			minVol = regionArray[curRegion].subResolution *
				regionArray[curRegion].subResolution * regionArray[curRegion].subResolution;
		
		// If region is a surface or membrane and it has children,
		// then the children should overlap all of its faces.
		checkSurfaceRegionChildren(curRegion, NUM_REGIONS, regionArray, bFail);
		
		// All children have been checked.
		// If region has no parent, check against remaining regions without a parent
		if(!regionArray[curRegion].bParent)
		{
			for(otherRegion = curRegion + 1; otherRegion < NUM_REGIONS; otherRegion++)
			{
				if(regionArray[otherRegion].bParent)
					continue;
				
				intersectType = intersectBoundary(regionArray[curRegion].spec.shape,
					regionArray[curRegion].boundary, regionArray[otherRegion].spec.shape,
					regionArray[otherRegion].boundary, intersection);
				
				bTemp = false;
				
				// Override intersection if Rectangle and Rectangular Box
				// By default, intersection in this case will be classified as Rectangular Box
				// so that the corresponding volume is 0.
				// Here, we still want to know if corresponding area is non-zero
				if((regionArray[curRegion].spec.shape == RECTANGLE
					&& regionArray[otherRegion].spec.shape == RECTANGULAR_BOX)
					|| (regionArray[otherRegion].spec.shape == RECTANGLE
					&& regionArray[curRegion].spec.shape == RECTANGULAR_BOX))
				{
					if(regionArray[curRegion].spec.shape == RECTANGLE)
					{
						squareRegion = curRegion;
						boxRegion = otherRegion;
					} else
					{
						boxRegion = curRegion;
						squareRegion = otherRegion;
					}
					
					if(boundaryVolume(RECTANGLE, intersection) > 
						squareDBL(regionArray[squareRegion].subResolution))
					{ 	// Intersection is possible.
						// Check whether Rectangle is on outer boundary of box (which is valid)
						switch(regionArray[squareRegion].plane)
						{
							case PLANE_XY:
								dim = 4;
								break;
							case PLANE_XZ:
								dim = 2;
								break;
							case PLANE_YZ:								
								dim = 0;
								break;
						}
						bTemp = !(fabs(regionArray[squareRegion].boundary[dim] -
							regionArray[boxRegion].boundary[dim]) < regionArray[squareRegion].subResolution
							|| fabs(regionArray[squareRegion].boundary[dim] -
							regionArray[boxRegion].boundary[dim+1]) < regionArray[squareRegion].subResolution);
					}
				} else if(intersectType == UNDEFINED_SHAPE ||
					boundaryVolume(intersectType, intersection) > minVol)
				{
					bTemp = true;
				}
				
				if(bTemp)
				{ // Intersection detected
					fprintf(stderr, "ERROR: Region %u (label: \"%s\") overlaps region %u (label: \"%s\"). Both regions have no parent.\n",
						curRegion, regionArray[curRegion].spec.label, otherRegion,
						regionArray[otherRegion].spec.label);
					*bFail = true;
				}
			}
		} else if(regionArray[curRegion].spec.shape == RECTANGLE)
		{ // Region is rectangle and has a parent
			otherRegion = regionArray[curRegion].parentID;
			if(regionArray[otherRegion].plane == PLANE_3D
				&& !regionArray[otherRegion].spec.bMicro)
			{ // Parent is 3D and mesoscopic.
				// Region must lie on boundary of its 3D neighbors
				// We need to know whether there are neighbors on both sides or region
				curVol = 0.;
				bDoubleSide = false;
				bTemp = false;
				for(curNeigh = 0;
					curNeigh < regionArray[curRegion].numRegionNeigh; curNeigh++)
				{
					neighID = regionArray[curRegion].regionNeighID[curNeigh];
					if(neighID == otherRegion
						|| regionArray[neighID].plane != PLANE_3D)
						continue; // Don't check against parent or 2D neighbor
					
					if(!bTemp)
					{
						bTemp = true;
						initDir = regionArray[curRegion].regionNeighDir[curNeigh];
					} else if(initDir != regionArray[curRegion].regionNeighDir[curNeigh])
					{
						bDoubleSide = true;
					}
					
					intersectBoundary(regionArray[curRegion].spec.shape,
					regionArray[curRegion].boundary, regionArray[neighID].spec.shape,
					regionArray[neighID].boundary, intersection);
					
					// Need to force boundary type to rectangle, otherwise
					// area of intersection with 3D neighbor will be 0
					curVol += boundaryVolume(RECTANGLE, intersection);
				}
				
				if((bDoubleSide
					&& fabs(curVol - 2* regionArray[curRegion].volume) > minVol)
					|| fabs(curVol - regionArray[curRegion].volume) > minVol)
				{
					fprintf(stderr, "ERROR: Region %u (label: \"%s\") is a rectangle child of mesoscopic 3D region %u (label: \"%s\") but is not on boundary of region %u's children.\n",
						curRegion, regionArray[curRegion].spec.label, otherRegion,
						regionArray[otherRegion].spec.label, otherRegion);
					*bFail = true;
				}
			}
		}
	}
}

// Check whether the children of a surface region cover all of the
// region's outer boundary
void checkSurfaceRegionChildren(const short curRegion,
	const short NUM_REGIONS,
	const struct region regionArray[],
	bool * bFail)
{	
	double curArea = 0.;
	double surfaceArea = 0.;
	short curChild, curChildID, curFace;
	double faceShared[6];
	
	if(regionArray[curRegion].spec.type == REGION_NORMAL
		|| regionArray[curRegion].numChildren == 0)
		return; // Region is normal or has no children
	
	switch (regionArray[curRegion].spec.shape)
	{
		case SPHERE:
		// Region should have only one child and radii must be the same
		// No need to check centers to be the same because we would have already confirmed that
		// child is valid
			if(regionArray[curRegion].numChildren > 1
				|| regionArray[curRegion].boundary[3] != regionArray[regionArray[curRegion].childrenID[0]].boundary[3]
				|| regionArray[regionArray[curRegion].childrenID[0]].spec.shape != SPHERE)
			{
				fprintf(stderr, "ERROR: Surface spherical region %u (label: \"%s\") has more than 1 child, or its child region %u (label: \"%s\") is not a sphere with the same radius.\n",
					curRegion, regionArray[curRegion].spec.label, regionArray[curRegion].childrenID[0],
					regionArray[regionArray[curRegion].childrenID[0]].spec.label);
				fprintf(stderr, "Radius of region %u is %.4em.\n",
					curRegion, regionArray[curRegion].boundary[3]);
				*bFail = true;
			}
			break;
		case RECTANGLE:
		case RECTANGULAR_BOX:
			// Measure surface area and compare with that shared with its children
			// Only need to check children that are the same shape
			surfaceArea = boundarySurfaceArea(regionArray[curRegion].spec.shape,
				regionArray[curRegion].boundary);
			
			for(curChild = 0; curChild < regionArray[curRegion].numChildren; curChild++)
			{
				curChildID = regionArray[curRegion].childrenID[curChild];
				if(regionArray[curRegion].spec.shape != regionArray[curChildID].spec.shape)
					continue; // Only proceed in this loop if shapes are the same
				
				// Determine if and where child's outer boundary shares outer surface of parent
				for(curFace = 0; curFace < 6; curFace++)
				{
					if(bSharedSurface(regionArray[curRegion].spec.shape,
						regionArray[curRegion].boundary,
						regionArray[curChildID].spec.shape,
						regionArray[curChildID].boundary,
						curFace, faceShared, regionArray[curRegion].subResolution))
					{							
						if(regionArray[curRegion].spec.shape == RECTANGLE)
						{ // Shared face is a line
							curArea += boundaryVolume(LINE, faceShared);
						} else
						{ // Shared face is a rectangle
							curArea += boundaryVolume(RECTANGLE, faceShared);
						}
					}
				}
			}
			
			// Is whole region surface covered by its children?
			if((regionArray[curRegion].spec.shape == RECTANGLE
				&& fabs(surfaceArea - curArea) < squareDBL(regionArray[curRegion].subResolution))
				|| (regionArray[curRegion].spec.shape == RECTANGULAR_BOX
				&& fabs(surfaceArea - curArea) < squareDBL(regionArray[curRegion].subResolution)
				* regionArray[curRegion].subResolution))
			{
				fprintf(stderr, "ERROR: Children of surface region %u (label: \"%s\") do not cover its entire boundary.\n",
					curRegion, regionArray[curRegion].spec.label);
				*bFail = true;
			}
			break;
		default:			
			fprintf(stderr,"ERROR: Surface region %u (label: \"%s\") shape %s invalid.\n", curRegion, regionArray[curRegion].spec.label, boundaryString(regionArray[curRegion].spec.shape));
			*bFail = true;
			return;
	}
}

// Does a one-sided surface exist between 2 rectangular boundaries?
// Such a surface would prevent the boundaries from being neighbors
bool bSurfaceBetweenBoundaries(const struct region regionArray[],
	const short NUM_REGIONS,
	const short region1,
	const short region2,
	const double boundary1[6],
	const double boundary2[6],
	unsigned short * surfaceRegion)
{
	double p1[3];
	double p2[3];
	double trajLine[3];
	double lineLength;
	short planeID;
	double d;
	
	for(*surfaceRegion  = 0; *surfaceRegion < NUM_REGIONS; (*surfaceRegion)++)
	{
		if(regionArray[*surfaceRegion].spec.type != REGION_NORMAL
			&& regionArray[*surfaceRegion].spec.surfaceType != SURFACE_MEMBRANE
			&& ((regionArray[region1].isRegionNeigh[*surfaceRegion]
				&& regionArray[region1].regionNeighDir[*surfaceRegion] != CHILD)
			|| (regionArray[region2].isRegionNeigh[*surfaceRegion]
				&& regionArray[region2].regionNeighDir[*surfaceRegion] != CHILD)))
		{ // surfaceRegion is a "one-sided" surface.
			// Make sure line between these subvolumes doesn't cross this region
			p1[0] = (boundary1[0] + boundary1[1])/2;
			p1[1] = (boundary1[2] + boundary1[3])/2;
			p1[2] = (boundary1[4] + boundary1[5])/2;
			p2[0] = (boundary2[0] + boundary2[1])/2;
			p2[1] = (boundary2[2] + boundary2[3])/2;
			p2[2] = (boundary2[4] + boundary2[5])/2;
			defineLine(p1, p2, trajLine, &lineLength);						
			if(bLineHitBoundary(p1, trajLine, lineLength,
				regionArray[*surfaceRegion].subShape, regionArray[*surfaceRegion].boundary,
				NULL, regionArray[*surfaceRegion].plane, false, &d, p2))
			{ // Subvolumes are not true neighbors
				return false;
			}
		}
	}
	return true;
}

// Does a boundary intersect a specified region?
// Determined by measuring volume of intersection.
// Intersection should be rectangular. Region can be spherical if it
// fully contains or is contained by the other boundary
// Round children must be entirely inside or outside intersection
bool bIntersectRegion(const short curRegion,
	const struct region regionArray[],
	const int boundary2Type,
	const double boundary2[])
{
	double minVol, volume;
	
	if(regionArray[curRegion].effectiveDim == DIM_3D)
	{
		minVol = regionArray[0].subResolution * regionArray[0].subResolution *
			regionArray[0].subResolution;
	} else if (regionArray[curRegion].effectiveDim == DIM_2D)
	{
		minVol = regionArray[0].subResolution * regionArray[0].subResolution;
	} else
	{
		minVol = regionArray[0].subResolution;
	}
	
	volume = intersectRegionVolume(curRegion, regionArray, boundary2Type,
		boundary2);
		
	return volume > minVol;
}