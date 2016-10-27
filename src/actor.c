/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * actor.c - operations on array of actors and its elements
 *
 * Last revised for AcCoRD v1.0 (2016-10-31)
 *
 * Revision history:
 *
 * Revision v1.0 (2016-10-31)
 * - added BURST modulation, which does not modulate binary data but always releases
 * molecules (of all types specified)
 * - moved mesoscopic structure fields from subvolume struct to meso subvolume struct
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - modified random number generation. Now use PCG via a separate interface file.
 * - added active point sources. Can be placed in microscopic or mesoscopic regions.
 * Cannot be on boundary of 2 or more regions or mesoscopic subvolumes
 * - added second order chemical reactions in the microscopic regime via a binding radius
 * and unbinding radius (can also model molecular crowding)
 * - made output of active actor data sequence a user option
 * - added bBits array for user to define a constant active actor bit sequence
 *
 * Revision v0.5.1 (2016-05-06)
 * - updated call to bPointInRegionNotChild to not exclude surface regions
 *
 * Revision v0.5 (2016-04-15)
 * - added ability to define location of actor by a list of regions
 * - added 2D and surface regions. Regions that have an effective dimension different
 * from their actual dimension cannot be intersected by a actor boundary (such regions
 * must be fully inside). Molecules will not be placed by an actor on a 2D region if
 * the actor overlaps at least 1 3D region
 * - tidied up calculations of subvolume coordinates
 *
 * Revision v0.4.1
 * - improved use and format of error messages
 *
 * Revision v0.4
 * - modified use of unsigned long and unsigned long long to uint32_t and uint64_t
 * - removed deprecated debug function
 * - added accommodation of spherical (microscopic) subvolumes and actors
 * - renamed region structure and array of structures to accommodate naming
 * of new boundary structure
 *
 * Revision v0.3.1
 * - header added
*/

#include "actor.h" // for "Public" declarations

//
// "Private" Declarations
//

//
// Definitions
//

/* Initialize array of structures for common actor parameters
* A call to this function should have a corresponding call to deleteActor2D,
* usually AFTER initializeActorCommon2D has been called */
void allocateActorCommonArray(const short NUM_ACTORS,
	struct actorStruct3D ** actorCommonArray)
{
	if(NUM_ACTORS > 0)
	{
		*actorCommonArray = malloc(NUM_ACTORS*sizeof(struct actorStruct3D));
			
		if(*actorCommonArray == NULL){
			fprintf(stderr, "ERROR: Memory allocation for array of actor structures.\n");
			exit(EXIT_FAILURE);
		}		
	} else
	{
		*actorCommonArray = NULL;
	}
	
}

void initializeActorCommon(const short NUM_ACTORS,
	struct actorStruct3D actorCommonArray[],
	const struct actorStructSpec3D actorCommonSpecArray[],
	const struct region regionArray[],
	const short NUM_REGIONS,
	short * NUM_ACTORS_ACTIVE,
	short * numActiveRecord,
	short ** activeRecordID,
	short * NUM_ACTORS_PASSIVE,
	short * numPassiveRecord,
	short ** passiveRecordID,
	const struct subvolume3D subvolArray[],
	uint32_t **** subID,
	uint32_t subCoorInd[][3],
	const double SUBVOL_BASE_SIZE)
{
	short curActor;
	short curRegion, curStr;
	short curInterRegion; // Current intersecting region
	
	int i; // Array index
	
	* NUM_ACTORS_ACTIVE = 0;
	* NUM_ACTORS_PASSIVE = 0;
	
	* numActiveRecord = 0;
	* numPassiveRecord = 0;
	short curPassiveRecord, curActiveRecord;
	
	// Used to find subvolumes inside actor
	uint32_t cur1, cur2, cur3, first1, first2, first3, last1, last2, last3;
	uint32_t curSub;
	uint32_t curInterSub;
	double curSubBound[6];
	
	bool bCurRegionIntersectActor;
	
	for(curActor = 0; curActor < NUM_ACTORS; curActor++)
	{		
		actorCommonArray[curActor].spec = actorCommonSpecArray[curActor];
		
		/*
		* Determine initialization parameters
		*/
		
		if(actorCommonArray[curActor].spec.bDefinedByRegions)
		{
			actorCommonArray[curActor].volume = 0;
		} else if(actorCommonArray[curActor].spec.shape == POINT)
		{ // Actor is a point
			if(!actorCommonArray[curActor].spec.bActive)
			{
				// Only active actors can be points
				fprintf(stderr, "ERROR: Actor %u is passive and defined as a point.\nOnly active actors can be points.\n",
					curActor);
				exit(EXIT_FAILURE);
			}
			actorCommonArray[curActor].volume = 0.;
		} else
		{
			actorCommonArray[curActor].volume =
				boundaryVolume(actorCommonArray[curActor].spec.shape,
					actorCommonArray[curActor].spec.boundary);
		}
		
		actorCommonArray[curActor].numRegion = 0;
		actorCommonArray[curActor].numRegionDim = 0;
		actorCommonArray[curActor].maxDim = 1;
		
		if(actorCommonArray[curActor].spec.bActive)
		{
			actorCommonArray[curActor].activeID = (*NUM_ACTORS_ACTIVE)++;
			if (actorCommonArray[curActor].spec.bWrite)
				(*numActiveRecord)++;
		} else{
			actorCommonArray[curActor].passiveID = (*NUM_ACTORS_PASSIVE)++;
			if (actorCommonArray[curActor].spec.bWrite)
				(*numPassiveRecord)++;
		}
		
		// Find number of regions within actor space
		for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
		{
			bCurRegionIntersectActor = false;
			if(actorCommonArray[curActor].spec.bDefinedByRegions)
			{
				
				for(curStr = 0;
					curStr < actorCommonArray[curActor].spec.numRegion;
					curStr++)
				{
					if(actorCommonArray[curActor].spec.regionLabel[curStr] &&
						strlen(actorCommonArray[curActor].spec.regionLabel[curStr]) > 0 &&
						!strcmp(regionArray[curRegion].spec.label,
						actorCommonArray[curActor].spec.regionLabel[curStr]))
					{ // The actor's location is defined by this region
						bCurRegionIntersectActor = true;
						break;
					}
				}
			} else if(bIntersectRegion(curRegion, regionArray,
				actorCommonArray[curActor].spec.shape,
				actorCommonArray[curActor].spec.boundary))
			{
				// Non-zero volume. Make sure that shape/regime combination is not
				// invalid (i.e., round actor in meso region)
				if(actorCommonArray[curActor].spec.shape == SPHERE
					&& !regionArray[curRegion].spec.bMicro
					&& bBoundarySurround(actorCommonArray[curActor].spec.shape,
					actorCommonArray[curActor].spec.boundary,
					regionArray[curRegion].spec.shape,
					regionArray[curRegion].boundary, 0.))
				{
					// Invalid actor for region
					fprintf(stderr, "ERROR: Round actor %u placed inside mesoscopic region %u.\n",
						curActor, curRegion);
					exit(EXIT_FAILURE);
				}
				
				if(actorCommonArray[curActor].spec.shape == POINT
					&& actorCommonArray[curActor].numRegion > 0)
				{
					fprintf(stderr, "ERROR: Actor %u is a point but is defined on the boundary of multiple regions.\nPoint actors must be defined within a single region.\n",
						curActor);
					exit(EXIT_FAILURE);
				}
				
				// If an actor intersects a 3D surface that is a 3D region, or a
				// 2D surface that is a 2D region, then it must surround the
				// entire region
				if(regionArray[curRegion].dimension != regionArray[curRegion].effectiveDim
					&& !bBoundarySurround(regionArray[curRegion].spec.shape,
					regionArray[curRegion].boundary,
					actorCommonArray[curActor].spec.shape,
					actorCommonArray[curActor].spec.boundary, 0.))
				{
					// Invalid actor for region
					fprintf(stderr, "ERROR: Actor %u intersects surface region %u.\n",
						curActor, curRegion);
					fprintf(stderr, "An actor's surface cannot intersect a 3D surface region that is a 3D shape or a 2D surface region that is a 2D shape.\n");
					exit(EXIT_FAILURE);
				}
				bCurRegionIntersectActor = true;
			}
			
			if(bCurRegionIntersectActor)
			{					
				// Region intersection is valid
				actorCommonArray[curActor].numRegion++;
				if(regionArray[curRegion].effectiveDim == DIM_3D)
				{ // Region is effectively 3D
					if(actorCommonArray[curActor].maxDim < 3)
						actorCommonArray[curActor].maxDim = 3;
				} else if(regionArray[curRegion].effectiveDim == DIM_2D)
				{ // Region is effectively 2D
					if(actorCommonArray[curActor].maxDim < 2)
						actorCommonArray[curActor].maxDim = 2;
				}
			}		
		}
		
		if(actorCommonArray[curActor].numRegion == 0)
		{
			fprintf(stderr,"ERROR: Actor %u placement is completely outside of simulation space.\n", curActor);
			exit(EXIT_FAILURE);
		}
		
		// Allocate structure memory for current actor
		actorCommonArray[curActor].regionID =
			malloc(actorCommonArray[curActor].numRegion*sizeof(unsigned short));
		actorCommonArray[curActor].bRegionInside =
			malloc(actorCommonArray[curActor].numRegion*sizeof(bool));
		actorCommonArray[curActor].numSub =
			malloc(actorCommonArray[curActor].numRegion*sizeof(uint32_t));
		actorCommonArray[curActor].regionInterType =
			malloc(actorCommonArray[curActor].numRegion*sizeof(unsigned short));
		actorCommonArray[curActor].regionInterBound =
			malloc(actorCommonArray[curActor].numRegion*sizeof(double[6]));
		actorCommonArray[curActor].regionInterArea =
			malloc(actorCommonArray[curActor].numRegion*sizeof(double));
		actorCommonArray[curActor].cumFracActorInRegion =
			malloc(actorCommonArray[curActor].numRegion*sizeof(double));
		actorCommonArray[curActor].subID =
			malloc(actorCommonArray[curActor].numRegion*sizeof(uint32_t *));
		if(actorCommonArray[curActor].regionID == NULL
			|| actorCommonArray[curActor].bRegionInside == NULL
			|| actorCommonArray[curActor].numSub == NULL
			|| actorCommonArray[curActor].regionInterType == NULL
			|| actorCommonArray[curActor].regionInterBound == NULL
			|| actorCommonArray[curActor].regionInterArea == NULL
			|| actorCommonArray[curActor].cumFracActorInRegion == NULL
			|| actorCommonArray[curActor].subID == NULL){
			fprintf(stderr,"ERROR: Memory allocation for structure members of actor %u.\n", curActor);
			exit(EXIT_FAILURE);
		}
		
		// If actor is defined by regions, determine actor volume
		if(actorCommonArray[curActor].spec.bDefinedByRegions)
		{
			bCurRegionIntersectActor = false;
			for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
			{
				for(curStr = 0;
					curStr < actorCommonArray[curActor].spec.numRegion;
					curStr++)
				{
					if(actorCommonArray[curActor].spec.regionLabel[curStr] &&
						strlen(actorCommonArray[curActor].spec.regionLabel[curStr]) > 0 &&
						!strcmp(regionArray[curRegion].spec.label,
						actorCommonArray[curActor].spec.regionLabel[curStr]) &&
						actorCommonArray[curActor].maxDim == regionArray[curRegion].effectiveDim)
					{ // The actor's location is defined by this region						
						actorCommonArray[curActor].volume += regionArray[curRegion].volume;
						break;
					}
				}
			}
		}
		
		// Determine IDs of regions within actor space
		curInterRegion = 0;
		for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
		{
			bCurRegionIntersectActor = false;
			if(actorCommonArray[curActor].spec.bDefinedByRegions)
			{				
				for(curStr = 0;
					curStr < actorCommonArray[curActor].spec.numRegion;
					curStr++)
				{
					if(actorCommonArray[curActor].spec.regionLabel[curStr] &&
						strlen(actorCommonArray[curActor].spec.regionLabel[curStr]) > 0 &&
						!strcmp(regionArray[curRegion].spec.label,
						actorCommonArray[curActor].spec.regionLabel[curStr]))
					{ // The actor's location is defined by this region
						actorCommonArray[curActor].regionInterType[curInterRegion] =
							regionArray[curRegion].spec.shape;
						for(i = 0; i < 6; i++)
						{
							actorCommonArray[curActor].regionInterBound[curInterRegion][i] =
								regionArray[curRegion].boundary[i];
						}
						actorCommonArray[curActor].regionInterArea[curInterRegion] =
							regionArray[curRegion].volume;
						
						// Is all of region inside the actor?
						actorCommonArray[curActor].bRegionInside[curInterRegion] = true;
						
						bCurRegionIntersectActor = true;
						break;
					}
				}
			} else if(bIntersectRegion(curRegion, regionArray,
				actorCommonArray[curActor].spec.shape,
				actorCommonArray[curActor].spec.boundary))
			{
				bCurRegionIntersectActor = true;
				
				if(actorCommonArray[curActor].spec.shape == POINT)
				{
					actorCommonArray[curActor].bRegionInside[curInterRegion] = false;
					actorCommonArray[curActor].regionInterArea[curInterRegion] = 0.;
					actorCommonArray[curActor].regionInterType[curInterRegion] = POINT;
					for(i = 0; i < 3; i++)
					{
						actorCommonArray[curActor].regionInterBound[curInterRegion][2*i] =
							actorCommonArray[curActor].spec.boundary[i];
						actorCommonArray[curActor].regionInterBound[curInterRegion][2*i+1] =
							actorCommonArray[curActor].spec.boundary[i];
					}
				} else
				{					
					// Find "outer" boundary of intersection (this includes space of
					// child regions)
					actorCommonArray[curActor].regionInterType[curInterRegion] = 
					intersectBoundary(actorCommonArray[curActor].spec.shape,
						actorCommonArray[curActor].spec.boundary,
						regionArray[curRegion].spec.shape, regionArray[curRegion].boundary,
						actorCommonArray[curActor].regionInterBound[curInterRegion]);
				
					// Find volume of intersection (this does exclude volumes of child regions)
					actorCommonArray[curActor].regionInterArea[curInterRegion] = 
						intersectRegionVolume(curRegion, regionArray,
						actorCommonArray[curActor].spec.shape,
						actorCommonArray[curActor].spec.boundary);
				
					// Is all of region inside the actor?
					actorCommonArray[curActor].bRegionInside[curInterRegion] =
						bBoundarySurround(regionArray[curRegion].spec.shape,
						regionArray[curRegion].boundary,
						actorCommonArray[curActor].spec.shape,
						actorCommonArray[curActor].spec.boundary, 0.);
				}
			}		
			
			// Is current region within actor space?
			if(bCurRegionIntersectActor)
			{
				actorCommonArray[curActor].regionID[curInterRegion] = curRegion;
				actorCommonArray[curActor].numSub[curInterRegion] = 0UL;
				
				
				// Determine (cumulative) fraction of actor in region
				if(curInterRegion > 0)
				{
					actorCommonArray[curActor].cumFracActorInRegion[curInterRegion] =
						actorCommonArray[curActor].cumFracActorInRegion[curInterRegion-1];
				} else
				{
					actorCommonArray[curActor].cumFracActorInRegion[curInterRegion] = 0.;
				}
				if(actorCommonArray[curActor].spec.shape == POINT)
				{ // Point actors must be entirely within their respective region
					actorCommonArray[curActor].cumFracActorInRegion[curInterRegion] = 1.;
				} else if(actorCommonArray[curActor].maxDim ==
					regionArray[curRegion].effectiveDim)
				{ 	// Do not add region volumes that are effectively of a lower dimension
					// than the actor
					actorCommonArray[curActor].cumFracActorInRegion[curInterRegion] +=
						actorCommonArray[curActor].regionInterArea[curInterRegion] /
						actorCommonArray[curActor].volume;
					actorCommonArray[curActor].numRegionDim++;
				}
				
				if(regionArray[curRegion].spec.bMicro)
				{
					// Currently no structure members exclusive to microscopic regions
				} else if(actorCommonArray[curActor].bRegionInside[curInterRegion])
				{
					// All region subvolumes are in the current actor
					actorCommonArray[curActor].numSub[curInterRegion] =
						regionArray[curRegion].numSub;
				} else
				{ // Region is mesoscopic
					
					// Find number of subvolumes in current region
					// that intersect the actor space.
					// An exhaustive search is not necessary, since we can use
					// the intersection boundary to limit the search.
					// Even if actor is spherical, regionInterBound is guaranteed
					// to be coordinates for a rectangular box in this case
					
					findSubSearchRange(regionArray, curRegion, curInterRegion, 
						actorCommonArray, curActor, &first1, &first2,
						&first3, &last1, &last2, &last3,
						false, 0);
					
					for(cur1 = first1; cur1 <= last1; cur1++)
					{
						if(regionArray[curRegion].dimension
							!= regionArray[curRegion].effectiveDim)
						{ // NOTE: We should never actually enter here because 
							// Regions of this type must be fully within the actor
							findSubSearchRange(regionArray, curRegion, curInterRegion,
								actorCommonArray, curActor, &first1, &first2,
								&first3, &last1, &last2, &last3,
								true, cur1);
						}
						
						for(cur2 = first2; cur2 <= last2; cur2++)
						{
							for(cur3 = first3; cur3 <= last3; cur3++)
							{
								// Is subvolume valid?
								if(subID[curRegion][cur1][cur2][cur3] == UINT32_MAX)
									continue; // Subvolume space is within a child
								
								findSubvolCoor(curSubBound, regionArray[curRegion],
									subCoorInd[subID[curRegion][cur1][cur2][cur3]]);
								
								if(bBoundaryIntersect(
									actorCommonArray[curActor].spec.shape,
									actorCommonArray[curActor].spec.boundary,
									regionArray[curRegion].subShape, curSubBound, 0.))
								{ // The subvolume does overlap the actor space
									actorCommonArray[curActor].numSub[curInterRegion]++;
								}
							}
						}
					}
					if(actorCommonArray[curActor].spec.shape == POINT
						&& actorCommonArray[curActor].numSub[curInterRegion] > 1)
					{ // Point actors can only be within one subvolume
						fprintf(stderr, "ERROR: Actor %u is a point in a mesoscopic region %u but is defined on the boundary of multiple subvolumes.\nPoint actors must be defined within a single subvolume if their parent region is mesoscopic.\n",
							curActor, curRegion);
						exit(EXIT_FAILURE);
					}
				}
				
				if(!regionArray[curRegion].spec.bMicro)
				{					
					// Allocate memory for IDs of subvolumes
					actorCommonArray[curActor].subID[curInterRegion] =
						malloc(actorCommonArray[curActor].numSub[curInterRegion]
							*sizeof(uint32_t));
					if(actorCommonArray[curActor].subID[curInterRegion] == NULL){
						fprintf(stderr,"ERROR: Memory allocation for structure members of actor %u.\n", curActor);
						exit(EXIT_FAILURE);
					}
				}
				
				// Memory assigned. Now record IDs of subvolumes within actor
				if(regionArray[curRegion].spec.bMicro)
				{
				} else if(actorCommonArray[curActor].bRegionInside[curInterRegion])
				{
					// All region subvolumes are in the current actor				
					for(curInterSub = 0;
						curInterSub < regionArray[curRegion].numSub; curInterSub++)
					{
						actorCommonArray[curActor].subID[curInterRegion][curInterSub] =
							subvolArray[regionArray[curRegion].firstID + curInterSub].mesoID;
					}
				} else
				{	
					curInterSub = 0;
					for(cur1 = first1; cur1 <= last1; cur1++)
					{
						if(regionArray[curRegion].dimension
							!= regionArray[curRegion].effectiveDim)
						{ // NOTE: We should never actually enter here because 
							// Regions of this type must be fully within the actor
							findSubSearchRange(regionArray, curRegion, curInterRegion,
								actorCommonArray, curActor, &first1, &first2,
								&first3, &last1, &last2, &last3,
								true, cur1);
						}
						
						for(cur2 = first2; cur2 <= last2; cur2++)
						{
							for(cur3 = first3; cur3 <= last3; cur3++)
							{
								// Is subvolume valid?
								if(subID[curRegion][cur1][cur2][cur3] == UINT32_MAX)
									continue; // Subvolume space is within a child
								
								// Confirm that subvolume is within actor space
								findSubvolCoor(curSubBound, regionArray[curRegion],
									subCoorInd[subID[curRegion][cur1][cur2][cur3]]);
								
								if(bBoundaryIntersect(
									actorCommonArray[curActor].spec.shape,
									actorCommonArray[curActor].spec.boundary,
									RECTANGULAR_BOX, curSubBound, 0.))
								{ // The subvolume does intersect the actor space
									curSub = subID[curRegion][cur1][cur2][cur3];
									actorCommonArray[curActor].subID[curInterRegion][curInterSub] = subvolArray[curSub].mesoID;
									curInterSub++;
								}
								
							}
						}
					}
				}
				
				curInterRegion++;
			}
		}
		
		// Check whether an active actor fully covers regions
		if(actorCommonArray[curActor].spec.bActive
			&& actorCommonArray[curActor].cumFracActorInRegion[actorCommonArray[curActor].numRegion-1]
			< 0.9999)
		{ // There is non-negligible actor space in an active actor that is not
			// covering a region.
			fprintf(stderr, "ERROR: Actor %u is active and only %.3f%% of its volume covers regions.\n",
				curActor, 100*actorCommonArray[curActor].cumFracActorInRegion[actorCommonArray[curActor].numRegion-1]);
			exit(EXIT_FAILURE);
		}
	}
	
	// Allocate memory for list of actors that record observations
	*passiveRecordID =
		malloc((*numPassiveRecord)*sizeof(short));
	*activeRecordID =
		malloc((*numActiveRecord)*sizeof(short));
	if((*numPassiveRecord > 0 && *passiveRecordID == NULL)
		|| (*numActiveRecord > 0 && *activeRecordID == NULL)){
		fprintf(stderr, "ERROR: Memory allocation for IDs of actors that will be recorded in the output file.\n");
		exit(EXIT_FAILURE);
	} else{ // There is at least one (passive) actor recording observations
		curPassiveRecord = 0;
		curActiveRecord = 0;
		for(curActor = 0; curActor < NUM_ACTORS; curActor++)
		{
			if(actorCommonArray[curActor].spec.bWrite)
			{
				if(actorCommonArray[curActor].spec.bActive)
					(*activeRecordID)[curActiveRecord++] = curActor;
				else
					(*passiveRecordID)[curPassiveRecord++] = curActor;
			}
		}
	}
}

void allocateActorActivePassiveArray(const short NUM_ACTORS_ACTIVE,
	struct actorActiveStruct3D ** actorActiveArray,
	const short NUM_ACTORS_PASSIVE,
	struct actorPassiveStruct3D ** actorPassiveArray)
{	
	if(NUM_ACTORS_ACTIVE > 0)
	{
		*actorActiveArray =
			malloc(NUM_ACTORS_ACTIVE*sizeof(struct actorActiveStruct3D));
			
		if(*actorActiveArray == NULL){
			fprintf(stderr, "ERROR: Memory allocation for array of active actor structures.\n");
			exit(EXIT_FAILURE);
		}		
	} else
	{
		*actorActiveArray = NULL;
	}
	
	if(NUM_ACTORS_PASSIVE > 0)
	{
		*actorPassiveArray =
			malloc(NUM_ACTORS_PASSIVE*sizeof(struct actorPassiveStruct3D));
			
		if(*actorPassiveArray == NULL){
			fprintf(stderr, "ERROR: Memory allocation for array of passive actor structures.\n");
			exit(EXIT_FAILURE);
		}		
	} else
	{
		*actorPassiveArray = NULL;
	}	
}

void initializeActorActivePassive(const short NUM_ACTORS,
	const struct actorStruct3D actorCommonArray[],
	const unsigned short NUM_MOL_TYPES,
	const struct region regionArray[],
	const short NUM_REGIONS,
	const struct mesoSubvolume3D mesoSubArray[],
	const short NUM_ACTORS_ACTIVE,
	struct actorActiveStruct3D actorActiveArray[],
	const short NUM_ACTORS_PASSIVE,
	struct actorPassiveStruct3D actorPassiveArray[],
	uint32_t subCoorInd[][3])
{
	int i; // loop index
	short curActor, curActive, curPassive;
	short curRegion, curInterRegion;
	short curPassiveRecord;
	uint32_t curSub, curInterSub, gamma;
	double curSubBound[6];
	double curInterSubBound[6];
	unsigned short curMolType, curMolRecord, curMolRecordPos;
	bool bRecordPos;
	
	// Link active and passive actor lists to the common actor list
	for(curActor = 0; curActor < NUM_ACTORS; curActor++)
	{
		if(actorCommonArray[curActor].spec.bActive)
		{
			actorActiveArray[actorCommonArray[curActor].activeID].actorID = curActor;
		} else{
			actorPassiveArray[actorCommonArray[curActor].passiveID].actorID = curActor;
		}
	}
	
	for(curActive = 0; curActive < NUM_ACTORS_ACTIVE; curActive++)
	{
		curActor = actorActiveArray[curActive].actorID;
		
		actorActiveArray[curActive].cumFracActorInSub =
			malloc(actorCommonArray[curActor].numRegion*sizeof(double *));
		actorActiveArray[curActive].molType =
			malloc(NUM_MOL_TYPES*sizeof(double *));
		
		initializeListRelease(&actorActiveArray[curActive].releaseList);
		initializeListData(&actorActiveArray[curActive].binaryData);
		
		if(actorActiveArray[curActive].cumFracActorInSub == NULL ||
			actorActiveArray[curActive].molType == NULL){
			fprintf(stderr,"ERROR: Memory allocation for structure members of active actor %u.\n", curActive);
			exit(EXIT_FAILURE);
		}
		
		actorActiveArray[curActive].alphabetSize = 1;
		for(i = 0; i < actorCommonArray[curActor].spec.modBits; i++)
		{
			actorActiveArray[curActive].alphabetSize *= 2U;
		}
		
		actorActiveArray[curActive].numMolType = 0;
		for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
		{
			if(actorCommonArray[curActor].spec.bReleaseMol[curMolType])
			{
				actorActiveArray[curActive].molType[actorActiveArray[curActive].numMolType++] =
					curMolType;
			}
		}
		
		for(curInterRegion = 0;
			curInterRegion < actorCommonArray[curActor].numRegion;
			curInterRegion++)
		{
			curRegion =
				actorCommonArray[curActor].regionID[curInterRegion];
				
			if(regionArray[curRegion].spec.bMicro
				|| actorCommonArray[curActor].bRegionInside[curInterRegion]
				|| (actorCommonArray[curActor].maxDim != regionArray[curRegion].effectiveDim))
				continue; // Following only applies if region is mesoscopic and not entirely inside actor
						
			actorActiveArray[curActive].cumFracActorInSub[curInterRegion] =
				malloc(actorCommonArray[curActor].numSub[curInterRegion]
				*sizeof(double));
		
			if(actorActiveArray[curActive].cumFracActorInSub[curInterRegion]
				== NULL){
				fprintf(stderr,"ERROR: Memory allocation for structure members of active actor %u.\n", curActive);
				exit(EXIT_FAILURE);
			}
			
			if(actorCommonArray[curActor].spec.shape == POINT)
			{ // Actor is a point. It must be entirely within 1 subvolume
				actorActiveArray[curActive].cumFracActorInSub[curInterRegion][0] = 1.;
				continue;
			}
			
			for(curInterSub = 0;
				curInterSub < actorCommonArray[curActor].numSub[curInterRegion];
				curInterSub++)
			{					
				curSub = mesoSubArray[actorCommonArray[curActor].subID[curInterRegion][curInterSub]].subID;
				
				// Determine boundary of current subvolume
				findSubvolCoor(curSubBound, regionArray[curRegion],
					subCoorInd[curSub]);
				
				// Determine intersection boundary of actor and subvolume
				intersectBoundary(actorCommonArray[curActor].spec.shape,
					actorCommonArray[curActor].spec.boundary,
					regionArray[curRegion].subShape, curSubBound, curInterSubBound);
					
				if (curInterSub > 0)
				{
					actorActiveArray[curActive].cumFracActorInSub[curInterRegion][curInterSub] = 
						actorActiveArray[curActive].cumFracActorInSub[curInterRegion][curInterSub-1];
				} else
				{
					actorActiveArray[curActive].cumFracActorInSub[curInterRegion][0] = 0.;
				}
				actorActiveArray[curActive].cumFracActorInSub[curInterRegion][curInterSub] +=
					boundaryVolume(regionArray[curRegion].subShape, curInterSubBound)
					/ actorCommonArray[curActor].regionInterArea[curInterRegion];
			}
		}
	}
	
	curPassiveRecord = 0;
	for(curPassive = 0; curPassive < NUM_ACTORS_PASSIVE; curPassive++)
	{
		curActor = actorPassiveArray[curPassive].actorID;
		actorPassiveArray[curPassive].numMolRecordID = 0;
		actorPassiveArray[curPassive].numMolRecordPosID = 0;
		actorPassiveArray[curPassive].recordID = SHRT_MAX;
		
		if (actorCommonArray[curActor].spec.bWrite)
		{
			actorPassiveArray[curPassive].recordID = curPassiveRecord++;
			for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
			{
				if (actorCommonArray[curActor].spec.bRecordMol[curMolType])
				{
					actorPassiveArray[curPassive].numMolRecordID++;
					if (actorCommonArray[curActor].spec.bRecordPos[curMolType])
					{
						actorPassiveArray[curPassive].numMolRecordPosID++;				
					}
				}
			}
		}
		
		actorPassiveArray[curPassive].fracSubInActor =
			malloc(actorCommonArray[curActor].numRegion*sizeof(double *));
		actorPassiveArray[curPassive].bRecordMesoAnyPos =
			malloc(actorCommonArray[curActor].numRegion*sizeof(bool));
		actorPassiveArray[curPassive].molRecordID =
			malloc(actorPassiveArray[curPassive].numMolRecordID*sizeof(unsigned short));
		actorPassiveArray[curPassive].molRecordPosID =
			malloc(actorPassiveArray[curPassive].numMolRecordPosID*sizeof(unsigned short));
		actorPassiveArray[curPassive].curMolObs =
			malloc(actorPassiveArray[curPassive].numMolRecordID*sizeof(uint64_t));
		actorPassiveArray[curPassive].subInterBound =
			malloc(actorCommonArray[curActor].numRegion
			*sizeof(double * [6]));
		
		if(actorPassiveArray[curPassive].fracSubInActor == NULL
			|| actorPassiveArray[curPassive].bRecordMesoAnyPos == NULL
			|| actorPassiveArray[curPassive].molRecordID == NULL
			|| actorPassiveArray[curPassive].molRecordPosID == NULL
			|| actorPassiveArray[curPassive].curMolObs == NULL
			|| actorPassiveArray[curPassive].subInterBound == NULL){
			fprintf(stderr,"ERROR: Memory allocation for structure members of active actor %u.\n", curPassive);
			exit(EXIT_FAILURE);
		}
		
		for(curInterRegion = 0;
			curInterRegion < actorCommonArray[actorPassiveArray[curPassive].actorID].numRegion;
			curInterRegion++)
		{
			curRegion =
				actorCommonArray[curActor].regionID[curInterRegion];
				
			if(regionArray[curRegion].spec.bMicro)
				continue; // Following only applies if region is mesoscopic
						
			actorPassiveArray[curPassive].bRecordMesoAnyPos[curInterRegion] = false;
			for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
			{
				if(actorCommonArray[curActor].spec.bWrite
					&& actorCommonArray[curActor].spec.bRecordMol[curMolType]
					&& actorCommonArray[curActor].spec.bRecordPos[curMolType])
				{ // We need to record the boundary of actor and subvolume intersection
					actorPassiveArray[curPassive].bRecordMesoAnyPos[curInterRegion] = true;
					actorPassiveArray[curPassive].subInterBound[curInterRegion] =
						malloc(actorCommonArray[curActor].numSub[curInterRegion]
						*sizeof(double [6]));
					if(actorPassiveArray[curPassive].subInterBound[curInterRegion]
						== NULL){
						fprintf(stderr,"ERROR: Memory allocation for structure members of passive actor %u.\n", curPassive);
						exit(EXIT_FAILURE);
					}
					break;
				}
			}
			
			actorPassiveArray[curPassive].fracSubInActor[curInterRegion] =
				malloc(actorCommonArray[curActor].numSub[curInterRegion]
				*sizeof(double));
		
			if(actorPassiveArray[curPassive].fracSubInActor[curInterRegion]
				== NULL){
				fprintf(stderr,"ERROR: Memory allocation for structure members of passive actor %d.\n", curPassive);
				exit(EXIT_FAILURE);
			}
			
			for(curInterSub = 0;
				curInterSub < actorCommonArray[curActor].numSub[curInterRegion];
				curInterSub++)
			{					
				curSub = mesoSubArray[actorCommonArray[curActor].subID[curInterRegion][curInterSub]].subID;
				
				// Determine boundary of current subvolume
				findSubvolCoor(curSubBound, regionArray[curRegion],
					subCoorInd[curSub]);
				
				// Determine intersection boundary of actor and subvolume
				if(actorCommonArray[curActor].bRegionInside[curInterRegion])
				{ // Entire subvolume must be inside actor
					for(i = 0; i < 6; i++)
						curInterSubBound[i] = curSubBound[i];
				} else
					intersectBoundary(actorCommonArray[curActor].spec.shape,
						actorCommonArray[curActor].spec.boundary,
						regionArray[curRegion].subShape, curSubBound, curInterSubBound);
				
				if(actorPassiveArray[curPassive].bRecordMesoAnyPos[curInterRegion])
				{ // Need to record boundary of intersection of actor and subvolume
					for(i = 0; i < 6; i++)
						actorPassiveArray[curPassive].subInterBound[curInterRegion][curInterSub][i] =
							curInterSubBound[i];
				}
				
				if(actorCommonArray[curActor].bRegionInside[curInterRegion])
					actorPassiveArray[curPassive].fracSubInActor[curInterRegion][curInterSub] = 1.;
				else
				{
					actorPassiveArray[curPassive].fracSubInActor[curInterRegion][curInterSub] =
						boundaryVolume(regionArray[curRegion].subShape, curInterSubBound);
					
					switch(regionArray[curRegion].effectiveDim)
					{ // Scale subvolume volume depending on number of dimensions
						// breaks are deliberately omitted in this control statement
						case DIM_3D:
							actorPassiveArray[curPassive].fracSubInActor[curInterRegion][curInterSub] /= regionArray[curRegion].actualSubSize;
						case DIM_2D:
							actorPassiveArray[curPassive].fracSubInActor[curInterRegion][curInterSub] /= regionArray[curRegion].actualSubSize;
						case DIM_1D:
							actorPassiveArray[curPassive].fracSubInActor[curInterRegion][curInterSub] /= regionArray[curRegion].actualSubSize;
					}
				}
			}
		}
		
		if (actorCommonArray[curActor].spec.bWrite)
		{
			curMolRecord = 0;
			curMolRecordPos = 0;
			for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
			{
				if (actorCommonArray[curActor].spec.bRecordMol[curMolType])
				{
					actorPassiveArray[curPassive].molRecordID[curMolRecord] = curMolType;
					if (actorCommonArray[curActor].spec.bRecordPos[curMolType])
					{
						actorPassiveArray[curPassive].molRecordPosID[curMolRecord] =
							curMolType;
						curMolRecordPos++;
					}
					curMolRecord++;
				}
			}
		}
	}
}

/* Reset actor members that vary with each simulation realization
*
*/
void resetActors(const short NUM_ACTORS,
	struct actorStruct3D actorCommonArray[],
	const unsigned short NUM_MOL_TYPES,
	const struct region regionArray[],
	const short NUM_REGIONS,
	const short NUM_ACTORS_ACTIVE,
	struct actorActiveStruct3D actorActiveArray[],
	const short NUM_ACTORS_PASSIVE,
	struct actorPassiveStruct3D actorPassiveArray[])
{
	short curActor;
	
	for(curActor = 0; curActor < NUM_ACTORS; curActor++)
	{
		actorCommonArray[curActor].nextTime =
			actorCommonArray[curActor].spec.startTime;
		actorCommonArray[curActor].curAction = 0UL;
	}
	
	for(curActor = 0; curActor < NUM_ACTORS_ACTIVE; curActor++)
	{
		if(!isListReleaseEmpty(&actorActiveArray[curActor].releaseList))
		{
			emptyListRelease(&actorActiveArray[curActor].releaseList);
			initializeListRelease(&actorActiveArray[curActor].releaseList);
		}
		if(!isListDataEmpty(&actorActiveArray[curActor].binaryData))
		{
			emptyListData(&actorActiveArray[curActor].binaryData);
			initializeListData(&actorActiveArray[curActor].binaryData);
		}
		
		// First action will be defining a new release and not the actual release of molecules
		actorActiveArray[curActor].curBit = 0;
		actorActiveArray[curActor].nextNewReleaseTime =
			actorCommonArray[actorActiveArray[curActor].actorID].spec.startTime;
		actorActiveArray[curActor].bNextActionNewRelease = true;
		actorActiveArray[curActor].nextEmissionTime = INFINITY;
		actorActiveArray[curActor].nextEmissionIndex = 0;
	}
}

/* Free Memory of Common Actors
* This function is needed to free the memory allocated to each
* array member of each structure element before releasing the memory
* allocated to the structure array itself. It should be called after BOTH
* allocateActorCommonArray3D AND initializeActorCommon3D, but is written to avoid
* any attempts to free memory that was not previously allocated */
void deleteActor(const short NUM_ACTORS,
	struct actorStruct3D actorCommonArray[],
	const struct region regionArray[],
	const short NUM_ACTORS_ACTIVE,
	struct actorActiveStruct3D actorActiveArray[],
	const short NUM_ACTORS_PASSIVE,
	struct actorPassiveStruct3D actorPassiveArray[],
	short passiveRecordID[],
	short activeRecordID[])
{
	short curActor, curActive, curPassive;
	short curInterRegion, curRegion;
	
	if(passiveRecordID != NULL)
		free(passiveRecordID);
	if(activeRecordID != NULL)
		free(activeRecordID);
	
	if(actorActiveArray != NULL)
	{		
		for(curActive = 0; curActive < NUM_ACTORS_ACTIVE; curActive++)
		{
			curActor = actorActiveArray[curActive].actorID;
			for(curInterRegion = 0;
				curInterRegion < actorCommonArray[curActor].numRegion;
				curInterRegion++)
			{
				curRegion =
					actorCommonArray[curActor].regionID[curInterRegion];
							
				if(!regionArray[curRegion].spec.bMicro
					&& !actorCommonArray[curActor].bRegionInside[curInterRegion]
					&& actorCommonArray[curActor].maxDim ==
					regionArray[curRegion].effectiveDim
					&& actorActiveArray[curActive].cumFracActorInSub[curInterRegion] != NULL)
					free(actorActiveArray[curActive].cumFracActorInSub[curInterRegion]);
			}
					
			if(actorActiveArray[curActive].cumFracActorInSub != NULL)
				free(actorActiveArray[curActive].cumFracActorInSub);
			if(actorActiveArray[curActive].molType != NULL)
				free(actorActiveArray[curActive].molType);
			
			emptyListRelease(&actorActiveArray[curActive].releaseList);
			emptyListData(&actorActiveArray[curActive].binaryData);
		}
		free(actorActiveArray);
	}
	
	if(actorPassiveArray != NULL)
	{		
		for(curPassive = 0; curPassive < NUM_ACTORS_PASSIVE; curPassive++)
		{
			curActor = actorPassiveArray[curPassive].actorID;
			for(curInterRegion = 0;
				curInterRegion < actorCommonArray[curActor].numRegion;
				curInterRegion++)
			{
				if(actorCommonArray[curActor].numSub[curInterRegion] > 0)
				{ // Region is mesoscopic and it intersects the actor
					if(actorPassiveArray[curPassive].fracSubInActor[curInterRegion] != NULL)
						free(actorPassiveArray[curPassive].fracSubInActor[curInterRegion]);
					if(actorPassiveArray[curPassive].bRecordMesoAnyPos[curInterRegion]
						&& actorPassiveArray[curPassive].subInterBound[curInterRegion] != NULL)
						free(actorPassiveArray[curPassive].subInterBound[curInterRegion]);
				}
				
			}
					
			if(actorPassiveArray[curPassive].fracSubInActor != NULL)
				free(actorPassiveArray[curPassive].fracSubInActor);
			if(actorPassiveArray[curPassive].bRecordMesoAnyPos != NULL)
				free(actorPassiveArray[curPassive].bRecordMesoAnyPos);
			if(actorPassiveArray[curPassive].molRecordID != NULL)
				free(actorPassiveArray[curPassive].molRecordID);
			if(actorPassiveArray[curPassive].molRecordPosID != NULL)
				free(actorPassiveArray[curPassive].molRecordPosID);
			if(actorPassiveArray[curPassive].curMolObs != NULL)
				free(actorPassiveArray[curPassive].curMolObs);
			if(actorPassiveArray[curPassive].subInterBound != NULL)
				free(actorPassiveArray[curPassive].subInterBound);
		}
		free(actorPassiveArray);
	}
	
	if(actorCommonArray != NULL)
	{		
		for(curActor = 0; curActor < NUM_ACTORS; curActor++)
		{
			for(curInterRegion = 0;
				curInterRegion < actorCommonArray[curActor].numRegion;
				curInterRegion++)
			{
				if(actorCommonArray[curActor].numSub[curInterRegion] > 0
					&& actorCommonArray[curActor].subID[curInterRegion] != NULL)
					free(actorCommonArray[curActor].subID[curInterRegion]);
			}
					
			if(actorCommonArray[curActor].regionID != NULL)
				free(actorCommonArray[curActor].regionID);
			if(actorCommonArray[curActor].bRegionInside != NULL)
				free(actorCommonArray[curActor].bRegionInside);
			if(actorCommonArray[curActor].numSub != NULL)
				free(actorCommonArray[curActor].numSub);
			if(actorCommonArray[curActor].regionInterType != NULL)
				free(actorCommonArray[curActor].regionInterType);
			if(actorCommonArray[curActor].regionInterBound != NULL)
				free(actorCommonArray[curActor].regionInterBound);
			if(actorCommonArray[curActor].regionInterArea != NULL)
				free(actorCommonArray[curActor].regionInterArea);
			if(actorCommonArray[curActor].cumFracActorInRegion != NULL)
				free(actorCommonArray[curActor].cumFracActorInRegion);
			if(actorCommonArray[curActor].subID != NULL)
				free(actorCommonArray[curActor].subID);
		}
		free(actorCommonArray);
	}
}

// Generate a new release and, if necessary, add it to the list
void newRelease(const struct actorStruct3D * actorCommon,
	struct actorActiveStruct3D * actorActive,
	double curTime)
{
	int i; // loop index
	NodeRelease * curRelease;
	
	// Parameters that need to be defined for new release (if there will actually be one)
	double strength, startTime, endTime, frequency;
	unsigned short molType;
	
	// Generate data
	ListData newData;
	initializeListData(&newData);
	
	for(i = 0; i < actorCommon->spec.modBits; i++)
	{
		if (actorCommon->spec.bRandBits)
		{
			if(generateUniform() < actorCommon->spec.probOne)
			{ // Bit is "1"
				if(!addData(&newData, true))
				{
					fprintf(stderr,"ERROR: Memory allocation for new bit.\n");
					exit(EXIT_FAILURE);
				}
			} else
			{ // Bit is "0"
				if(!addData(&newData, false))
				{
					fprintf(stderr,"ERROR: Memory allocation for new bit.\n");
					exit(EXIT_FAILURE);
				}
			}
		} else
		{
			// Bits are not random. Read from pre-determined list
			if(actorCommon->spec.bBits[actorActive->curBit++])
			{ // Bit is "1"
				if(!addData(&newData, true))
				{
					fprintf(stderr,"ERROR: Memory allocation for new bit.\n");
					exit(EXIT_FAILURE);
				}
			} else
			{ // Bit is "0"
				if(!addData(&newData, false))
				{
					fprintf(stderr,"ERROR: Memory allocation for new bit.\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}
	
	// Append newly-generated data to existing data
	transferData(&actorActive->binaryData, &newData);
	
	// Translate new bits into release information, based on modulation scheme
	switch (actorCommon->spec.modScheme)
	{
		case CSK:
			// CSK: Concentration shift keying. Convert data to an unsigned int value
			strength = ((double) binaryToDecimal(&newData, actorActive->alphabetSize))
				* actorCommon->spec.modStrength;
			startTime = 0.;
			endTime = actorCommon->spec.releaseInterval;
			frequency = 0.;
			if(actorCommon->spec.bTimeReleaseRand)
			{ // Emission times are stochastic. First emission time must be generated.
				startTime = generateExponential(1)/strength;
			} else
			{ // Emission times are deterministic. First emission will be at start of interval
				startTime = 0.;
			}
			molType = actorActive->molType[0]; // There is only one kind of molecule in CSK
			break;
		case BURST:
			// BURST: Like CSK but always 1 bit-1.
			// Also, can release multiple types of molecules
			strength = ((double) binaryToDecimal(&newData, actorActive->alphabetSize))
				* actorCommon->spec.modStrength;
			startTime = 0.;
			endTime = actorCommon->spec.releaseInterval;
			frequency = 0.;
			if(actorCommon->spec.bTimeReleaseRand)
			{ // Emission times are stochastic. First emission time must be generated.
				startTime = generateExponential(1)/strength;
			} else
			{ // Emission times are deterministic. First emission will be at start of interval
				startTime = 0.;
			}
			break;
		default:
			fprintf(stderr,
				"ERROR: Modulation scheme %d invalid.\n",
				actorCommon->spec.modScheme);
			exit(EXIT_FAILURE);
	}
	
	// Append release information to list of current releases
	if(strength > 0.)
	{
		switch (actorCommon->spec.modScheme)
		{
			case BURST:
				// This modulation types can release multiple types of molecules
				for(i = 0; i < actorActive->numMolType; i++)
				{
					if(!addRelease(&actorActive->releaseList, strength, actorActive->molType[i],
						curTime + startTime, curTime + endTime, frequency))
					{
						fprintf(stderr,"ERROR: Memory allocation for new active actor release.\n");
						exit(EXIT_FAILURE);
					}
				}
				break;
			default:
				if(!addRelease(&actorActive->releaseList, strength, molType,
					curTime + startTime, curTime + endTime, frequency))
				{
					fprintf(stderr,"ERROR: Memory allocation for new active actor release.\n");
					exit(EXIT_FAILURE);
				}
		}
		
		// Update time of next emission event
		findNextEmission(actorCommon, actorActive);
	}
	
	// Free memory of new data
	emptyListData(&newData);
}

// Find time and index of next release to have an emission
void findNextEmission(const struct actorStruct3D * actorCommon,
	struct actorActiveStruct3D * actorActive)
{
	unsigned int i = 0; // Counter
	NodeRelease * curRelease = actorActive->releaseList;
	actorActive->nextEmissionTime = INFINITY;

	while (curRelease != NULL)
	{
		if (curRelease->item.nextTime < actorActive->nextEmissionTime)
		{
			actorActive->nextEmissionTime = curRelease->item.nextTime;
			actorActive->nextEmissionIndex = i;
		}
		curRelease = curRelease->next;
		i++;
	}
}

// Release molecules for current release at this instant
void fireEmission(const struct actorStruct3D * actorCommon,
	struct actorActiveStruct3D * actorActive,
	struct region region[],
	const short NUM_REGIONS,
	struct subvolume3D subvolArray[],
	struct mesoSubvolume3D * mesoSubArray,
	const uint32_t numMesoSub,
	const unsigned short NUM_MOL_TYPES,
	ListMolRecent3D microMolListRecent[NUM_REGIONS][NUM_MOL_TYPES],
	double tMicro,
	uint32_t * heap_subvolID,
	uint32_t (*heap_childID)[2],
	bool (*b_heap_childValid)[2])
{
	NodeRelease * curRelease = actorActive->releaseList;
	unsigned int curReleaseInd;
	bool bRemoveRelease = false;
	uint64_t numNewMol;
	
	// Traverse release list to the current release
	for(curReleaseInd = 0;
		curReleaseInd < actorActive->nextEmissionIndex;
		curReleaseInd++)
	{
		curRelease = curRelease->next;
	}	
	
	// Fire emission and find time of next emission for this release
	if(actorCommon->spec.bTimeReleaseRand)
	{
		// We only release one molecule right now
		placeMolecules(actorCommon, actorActive, region, NUM_REGIONS,
			subvolArray, mesoSubArray, numMesoSub, (uint64_t) 1,
			curRelease->item.molType, NUM_MOL_TYPES, microMolListRecent,
			curRelease->item.nextTime, tMicro,
			heap_subvolID, heap_childID, b_heap_childValid);
		
		
		// Next emission time must be generated
		curRelease->item.nextTime +=
			generateExponential(1)/curRelease->item.strength;
	} else
	{
		if(actorCommon->spec.bNumReleaseRand)
		{
			numNewMol = generatePoisson(curRelease->item.strength);
		} else
		{
			numNewMol = (uint64_t) curRelease->item.strength;
		}
		// We must place curRelease->item.strength molecules
		placeMolecules(actorCommon, actorActive, region, NUM_REGIONS,
			subvolArray, mesoSubArray, numMesoSub, numNewMol,
			curRelease->item.molType, NUM_MOL_TYPES, microMolListRecent,
			curRelease->item.nextTime, tMicro,
			heap_subvolID, heap_childID, b_heap_childValid);		
		
		// Use slot interval
		if (actorCommon->spec.slotInterval > 0)
		{
			curRelease->item.nextTime += actorCommon->spec.slotInterval;
		} else
		{ // slot interval is 0; no need to check for the end of the slot time
			bRemoveRelease = true;
		}
	}
	
	if(bRemoveRelease
		|| curRelease->item.nextTime > curRelease->item.endTime)
	{ // Release has finished. Remove from list
		deleteRelease(&actorActive->releaseList,
			actorActive->nextEmissionIndex);
	}
	
	// Find time and release of next emission event
	findNextEmission(actorCommon, actorActive);
}

// Place molecules for current emission
void placeMolecules(const struct actorStruct3D * actorCommon,
	const struct actorActiveStruct3D * actorActive,
	struct region region[],
	const short NUM_REGIONS,
	struct subvolume3D subvolArray[],
	struct mesoSubvolume3D * mesoSubArray,
	const uint32_t numMesoSub,
	uint64_t numNewMol,
	const unsigned short curMolType,
	const unsigned short NUM_MOL_TYPES,
	ListMolRecent3D microMolListRecent[NUM_REGIONS][NUM_MOL_TYPES],
	double tCur,
	double tMicro,
	uint32_t * heap_subvolID,
	uint32_t (*heap_childID)[2],
	bool (*b_heap_childValid)[2])
{
	short curRegion, curRegionInter, curRegionDim;
	double uniRV;
	uint64_t curMolecule;
	//double tCur = curRelease->item.nextTime;
	
	if(actorCommon->numRegionDim > 1)
	{ // Molecules can end up in different regions. Place one at a time
		for(curMolecule = 0;
			curMolecule < numNewMol;
			curMolecule++)
		{
			uniRV = generateUniform();
			for(curRegionInter = 0;
				actorCommon->cumFracActorInRegion[curRegionInter] < uniRV;
				curRegionInter++)
			{ // Scanning regions until we find where molecule must be place
				if(curRegionInter >= actorCommon->numRegion)
				{
					fprintf(stderr,"\nWARNING: New molecule does not have a valid region to be placed in.\n");
					break; // end for-loop
				}
			}
			// Molecule will be placed in curRegionInter
			if(curRegionInter < actorCommon->numRegion)
			{
				curRegion = actorCommon->regionID[curRegionInter];
				placeMoleculesInRegion(actorCommon, actorActive, region,
					curRegion, curRegionInter, NUM_REGIONS, subvolArray,
					mesoSubArray, numMesoSub, (uint64_t) 1, curMolType, NUM_MOL_TYPES,
					&microMolListRecent[curRegion][curMolType], tCur, tMicro,
					heap_subvolID, heap_childID, b_heap_childValid);
			}
		}
		
	} else
	{ // All molecules are going in the same region.
		// Find region of same dimension as actor
		for(curRegionInter = 0; curRegionInter < actorCommon->numRegion; curRegionInter++)
		{
			curRegion = actorCommon->regionID[curRegionInter];
			if(region[curRegion].effectiveDim == actorCommon->maxDim)
				break;
		}
		placeMoleculesInRegion(actorCommon, actorActive, region, curRegion,
			curRegionInter, NUM_REGIONS, subvolArray,
			mesoSubArray, numMesoSub, numNewMol, curMolType, NUM_MOL_TYPES,
			&microMolListRecent[curRegion][curMolType], tCur, tMicro,
			heap_subvolID, heap_childID, b_heap_childValid);
	}
}

// Place molecules for current actor emission in specific region
void placeMoleculesInRegion(const struct actorStruct3D * actorCommon,
	const struct actorActiveStruct3D * actorActive,
	struct region regionArray[],
	const short curRegion,
	const short curRegionInter,
	const short NUM_REGIONS,
	struct subvolume3D subvolArray[],
	struct mesoSubvolume3D * mesoSubArray,
	const uint32_t numMesoSub,
	uint64_t numNewMol,
	const unsigned short curMolType,
	const unsigned short NUM_MOL_TYPES,
	ListMolRecent3D * microMolListRecent,
	double tCur,
	double tMicro,
	uint32_t * heap_subvolID,
	uint32_t (*heap_childID)[2],
	bool (*b_heap_childValid)[2])
{
	uint64_t curMolecule;
	double uniRV;
	uint32_t curSub, curSubInter;
	bool bNeedPoint;
	double point[3];
	bool bSurface = regionArray[curRegion].effectiveDim
		!= regionArray[curRegion].dimension;
		
	if(regionArray[curRegion].spec.bMicro)
	{
		// Add molecules uniformly within micro intersect region
		for(curMolecule = 0;
			curMolecule < numNewMol;
			curMolecule++)
		{			
			bNeedPoint = true;
			while(bNeedPoint)
			{
				// Intersection region is defined by actorCommon->regionInterBound
				// Shape of intersection is defined by actorCommon->regionInterType
				uniformPointVolume(point, actorCommon->regionInterType[curRegionInter],
					actorCommon->regionInterBound[curRegionInter], bSurface,
					regionArray[curRegion].plane);
				
				if(actorCommon->spec.shape == POINT
					|| bPointInRegionNotChild(curRegion, regionArray, point, false))
				{
					bNeedPoint = false;
					if(!addMoleculeRecent(microMolListRecent, point[0], point[1],
						point[2], tMicro - tCur))
					{ // Creation of molecule failed
						fprintf(stderr, "ERROR: Memory allocation for new molecule to be placed in region %u.\n", curRegion);
						exit(EXIT_FAILURE);
					}
				}			
			}
		}
	} else
	{
		// Find subvolume in region to add molecule
		if(actorCommon->numSub[curRegionInter] > 1)
		{ // Molecules could end up in different subvolumes. Place one at a time.
			for(curMolecule = 0;
			curMolecule < numNewMol;
			curMolecule++)
			{
				uniRV = generateUniform();
				if(actorCommon->bRegionInside[curRegionInter])
				{ // All subvolumes in the region are equally likely because entire region is in the actor
					curSubInter = (uint32_t) floor(uniRV * actorCommon->numSub[curRegionInter]);
				} else
				{ // Need to consider the individual likelihoods for each subvolume
					for(curSubInter = 0;
					actorActive->cumFracActorInSub[curRegionInter][curSubInter] < uniRV;
					curSubInter++)
					{ // Scanning subvolumes until we find where molecule must be place
						if(curSubInter >= actorCommon->numSub[curRegionInter])
						{
							fprintf(stderr,"\nWARNING: New molecule placed in region %u does not have a valid subvolume to be placed in.\n", curRegion);
							break; // end for-loop
						}
					}
				}
				curSub = actorCommon->subID[curRegionInter][curSubInter];
				// Molecule will be placed in curSubInter
				placeMoleculesInSub(regionArray, subvolArray, mesoSubArray,
					numMesoSub, (uint64_t) 1, curMolType, curSub, NUM_MOL_TYPES,
					tCur, NUM_REGIONS, heap_subvolID, heap_childID, b_heap_childValid);
			}
		} else
		{ // All molecules are going in the same subvolume
			curSub = actorCommon->subID[curRegionInter][0];
			placeMoleculesInSub(regionArray, subvolArray, mesoSubArray,
				numMesoSub, numNewMol, curMolType, curSub, NUM_MOL_TYPES,
				tCur, NUM_REGIONS, heap_subvolID,
				heap_childID, b_heap_childValid);
		}
	}
}

// Add one type of molecule in specific subvolume
void placeMoleculesInSub(struct region regionArray[],
	struct subvolume3D subvolArray[],
	struct mesoSubvolume3D * mesoSubArray,
	const uint32_t numMesoSub,
	uint64_t numNewMol,
	const unsigned short curMolType,
	const uint32_t curMeso,
	const unsigned short NUM_MOL_TYPES,
	double tCur,
	const short NUM_REGIONS,
	uint32_t * heap_subvolID,
	uint32_t (*heap_childID)[2],
	bool (*b_heap_childValid)[2])
{
	uint32_t curSub; // Sub ID in global subvolume list
	
	// Define which molecules are added where
	mesoSubArray[curMeso].num_mol[curMolType] += numNewMol;
	curSub = mesoSubArray[curMeso].subID;
	
	// Place molecule(s) and update heap for subvolumes
	updateMesoSub(curSub, false, (uint64_t []){numNewMol},
		(bool []){true}, curMolType,
		true, numMesoSub, mesoSubArray, subvolArray, tCur,
		NUM_REGIONS, NUM_MOL_TYPES, regionArray);
	heapMesoUpdate(numMesoSub, mesoSubArray, heap_subvolID,
		mesoSubArray[curMeso].heapID, heap_childID, b_heap_childValid);
}

// Find range of subvolumes to search over for intersection with an actor
void findSubSearchRange(const struct region regionArray[],
	const short curRegion,
	const short curInterRegion,
	const struct actorStruct3D actorCommonArray[],
	const short curActor,
	uint32_t * first1,
	uint32_t * first2,
	uint32_t * first3,
	uint32_t * last1,
	uint32_t * last2,
	uint32_t * last3,
	bool bHaveCur1,
	uint32_t cur1)
{
	uint32_t maxSize[3];
	short dim[6];
	double yAnch, zAnch;
	
	if(bHaveCur1)
	{
		// Region has more than one face and we have current face
		// Applies to 3D surfaces that are 3D regions and 2D surfaces that are 2D regions
		// These kinds of regions are always surrounded by actor, so we have to search
		// all of the subvolumes
		*first2 = 0;
		*first3 = 0;
		if(regionArray[curRegion].spec.shape == RECTANGULAR_BOX)
		{
			switch(cur1)
			{
				case 0:
				case 1:
					*last2 = regionArray[curRegion].spec.numX-1;
					*last3 = regionArray[curRegion].spec.numY-1;
					break;
				case 2:
				case 3:
					*last2 = regionArray[curRegion].spec.numX-1;
					*last3 = regionArray[curRegion].spec.numZ-1;
					break;
				case 4:
				case 5:
					*last2 = regionArray[curRegion].spec.numY-1;
					*last3 = regionArray[curRegion].spec.numZ-1;
					break;
				default:
					// Something went wrong
					fprintf(stderr, "ERROR: Intersection between actor %u and meso subvolumes of region %u could not be determined.\n",
						curActor, curRegion);
					exit(EXIT_FAILURE);
			}
		} else if (regionArray[curRegion].spec.shape == RECTANGLE)
		{
			// Third dimension will always have just length 1
			*last3 = 1;
			switch(regionArray[curRegion].plane)
			{
				case PLANE_XY:
					switch(cur1)
					{
						case 0:
						case 1:
							*last2 = regionArray[curRegion].spec.numX-1;
						case 2:
						case 3:
							*last2 = regionArray[curRegion].spec.numY-1;
						default:
							// Something went wrong
							fprintf(stderr, "ERROR: Intersection between actor %u and meso subvolumes of region %u could not be determined.\n",
								curActor, curRegion);
							exit(EXIT_FAILURE);
					}
					break;
				case PLANE_XZ:
					switch(cur1)
					{
						case 0:
						case 1:
							*last2 = regionArray[curRegion].spec.numX-1;
						case 2:
						case 3:
							*last2 = regionArray[curRegion].spec.numZ-1;
						default:
							// Something went wrong
							fprintf(stderr, "ERROR: Intersection between actor %u and meso subvolumes of region %u could not be determined.\n",
								curActor, curRegion);
							exit(EXIT_FAILURE);
					}
					break;
				case PLANE_YZ:
					switch(cur1)
					{
						case 0:
						case 1:
							*last2 = regionArray[curRegion].spec.numY-1;
						case 2:
						case 3:
							*last2 = regionArray[curRegion].spec.numZ-1;
						default:
							// Something went wrong
							fprintf(stderr, "ERROR: Intersection between actor %u and meso subvolumes of region %u could not be determined.\n",
								curActor, curRegion);
							exit(EXIT_FAILURE);
					}
					break;
				default:
					// Something went wrong
					fprintf(stderr, "ERROR: Intersection between actor %u and meso subvolumes of region %u could not be determined.\n",
						curActor, curRegion);
					exit(EXIT_FAILURE);
			}
		}
		return;
	} else if(regionArray[curRegion].numFace > 0)
	{
		// Indices along first dimension will be all of the faces
		*first1 = 0;
		*last1 = regionArray[curRegion].numFace-1;
		maxSize[0] = *last1;
		if(regionArray[curRegion].numFace == 1)
		{
			// Meso subvolumes with 1 face must be in a rectangle region
			switch(regionArray[curRegion].plane)
			{
				case PLANE_XY:
					dim[2] = 0;
					dim[3] = 1;
					dim[4] = 2;
					dim[5] = 3;
					yAnch = regionArray[curRegion].spec.xAnch;
					zAnch = regionArray[curRegion].spec.yAnch;
					maxSize[1] = regionArray[curRegion].spec.numX-1;
					maxSize[2] = regionArray[curRegion].spec.numY-1;
					break;
				case PLANE_XZ:
					dim[2] = 0;
					dim[3] = 1;
					dim[4] = 4;
					dim[5] = 5;
					yAnch = regionArray[curRegion].spec.xAnch;
					zAnch = regionArray[curRegion].spec.zAnch;
					maxSize[1] = regionArray[curRegion].spec.numX-1;
					maxSize[2] = regionArray[curRegion].spec.numZ-1;
					break;
				case PLANE_YZ:
					dim[2] = 2;
					dim[3] = 3;
					dim[4] = 4;
					dim[5] = 5;
					yAnch = regionArray[curRegion].spec.yAnch;
					zAnch = regionArray[curRegion].spec.zAnch;
					maxSize[1] = regionArray[curRegion].spec.numY-1;
					maxSize[2] = regionArray[curRegion].spec.numZ-1;
					break;
				default:
					// Something went wrong
					fprintf(stderr, "ERROR: Intersection between actor %u and meso subvolumes of region %u could not be determined.\n",
						curActor, curRegion);
					exit(EXIT_FAILURE);
			}
		} else
		{ // Cannot yet determine the remaining indices
			return;
		}		
	}
	else
	{ // Region is not a surface. Dimensions are default
		*first1 = (uint32_t) floor(
			(actorCommonArray[curActor].regionInterBound[curInterRegion][0]
			- regionArray[curRegion].spec.xAnch) /regionArray[curRegion].actualSubSize);
		*last1 = (uint32_t) ceil(
			(actorCommonArray[curActor].regionInterBound[curInterRegion][1]
			- regionArray[curRegion].spec.xAnch - regionArray[curRegion].actualSubSize)
			/regionArray[curRegion].actualSubSize);
			
		dim[2] = 2;
		dim[3] = 3;
		dim[4] = 4;
		dim[5] = 5;	
		yAnch = regionArray[curRegion].spec.yAnch;
		zAnch = regionArray[curRegion].spec.zAnch;
		maxSize[0] = regionArray[curRegion].spec.numX-1;
		maxSize[1] = regionArray[curRegion].spec.numY-1;
		maxSize[2] = regionArray[curRegion].spec.numZ-1;
	}
	
	*first2 = (uint32_t) floor(
		(actorCommonArray[curActor].regionInterBound[curInterRegion][dim[2]] - yAnch)
		/regionArray[curRegion].actualSubSize);
	*first3 = (uint32_t) floor(
		(actorCommonArray[curActor].regionInterBound[curInterRegion][dim[4]] - zAnch)
		/regionArray[curRegion].actualSubSize);
	*last2 = (uint32_t) ceil(
		(actorCommonArray[curActor].regionInterBound[curInterRegion][dim[3]] - yAnch - regionArray[curRegion].actualSubSize)
		/regionArray[curRegion].actualSubSize);
	*last3 = (uint32_t) ceil(
		(actorCommonArray[curActor].regionInterBound[curInterRegion][dim[5]] - zAnch - regionArray[curRegion].actualSubSize)
		/regionArray[curRegion].actualSubSize);
	
	// Correct for rounding that takes us beyond number of subvolumes in any dimension
	if(*last1 > maxSize[0]) *last1 = maxSize[0];
	if(*last2 > maxSize[1]) *last2 = maxSize[1];
	if(*last3 > maxSize[2]) *last3 = maxSize[2];
}