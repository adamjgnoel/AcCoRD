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
 * Last revised for AcCoRD v0.4.1
 *
 * Revision history:
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
void allocateActorCommonArray3D(const short NUM_ACTORS,
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

void initializeActorCommon3D(const short NUM_ACTORS,
	struct actorStruct3D actorCommonArray[],
	const struct actorStructSpec3D actorCommonSpecArray[],
	const struct region regionArray[],
	const short NUM_REGIONS,
	short * NUM_ACTORS_ACTIVE,
	short * NUM_ACTORS_PASSIVE,
	short * numActorRecord,
	short ** actorRecordID,
	uint32_t **** subID,
	const double SUBVOL_BASE_SIZE)
{
	short curActor;
	short curRegion;
	short curInterRegion; // Current intersecting region
	
	int i; // Array index
	
	* NUM_ACTORS_ACTIVE = 0;
	* NUM_ACTORS_PASSIVE = 0;
	
	* numActorRecord = 0;
	short curActorRecord;
	
	// Used to find subvolumes inside actor
	uint32_t curX, curY, curZ, firstX, firstY, firstZ, lastX, lastY, lastZ, firstSub;
	uint32_t curSub;
	uint32_t curInterSub;
	double curSubBound[6];
	double minVolume = regionArray[0].subResolution * regionArray[0].subResolution
		* regionArray[0].subResolution;
	
	for(curActor = 0; curActor < NUM_ACTORS; curActor++)
	{		
		actorCommonArray[curActor].spec = actorCommonSpecArray[curActor];
		
		/*
		* Determine initialization parameters
		*/
		
		actorCommonArray[curActor].volume =
			boundaryVolume(actorCommonArray[curActor].spec.shape,
				actorCommonArray[curActor].spec.boundary);
		
		actorCommonArray[curActor].numRegion = 0;
		
		if(actorCommonArray[curActor].spec.bActive)
		{
			actorCommonArray[curActor].activeID = (*NUM_ACTORS_ACTIVE)++;
		} else{
			actorCommonArray[curActor].passiveID = (*NUM_ACTORS_PASSIVE)++;
			if (actorCommonArray[curActor].spec.bWrite)
				(*numActorRecord)++;
		}
		
		// Find number of regions within actor space
		for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
		{
			// Do the region and actor overlap?
			if(intersectRegionVolume(curRegion, regionArray,
				actorCommonArray[curActor].spec.shape,
				actorCommonArray[curActor].spec.boundary) > minVolume)
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
				actorCommonArray[curActor].numRegion++;
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
		
		// Determine IDs of regions within actor space
		curInterRegion = 0;
		for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
		{
			// Is current region within actor space?
			if(intersectRegionVolume(curRegion, regionArray,
				actorCommonArray[curActor].spec.shape,
				actorCommonArray[curActor].spec.boundary) > minVolume)
			{
				actorCommonArray[curActor].regionID[curInterRegion] = curRegion;
				actorCommonArray[curActor].numSub[curInterRegion] = 0UL;
				
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
				actorCommonArray[curActor].bRegionInside[curInterRegion] = fabs(actorCommonArray[curActor].regionInterArea[curInterRegion] -
					regionArray[curRegion].volume) < SUBVOL_BASE_SIZE*
					SUBVOL_BASE_SIZE*SUBVOL_BASE_SIZE*SUB_ADJ_RESOLUTION;
				
				// Determine (cumulative) fraction of actor in region
				if(curInterRegion > 0)
				{
					actorCommonArray[curActor].cumFracActorInRegion[curInterRegion] =
						actorCommonArray[curActor].cumFracActorInRegion[curInterRegion-1];
				} else
				{
					actorCommonArray[curActor].cumFracActorInRegion[curInterRegion] = 0.;
				}
				actorCommonArray[curActor].cumFracActorInRegion[curInterRegion] +=
					actorCommonArray[curActor].regionInterArea[curInterRegion] /
					actorCommonArray[curActor].volume;
				
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
					firstX = (uint32_t) floor(
						(actorCommonArray[curActor].regionInterBound[curInterRegion][0] - regionArray[curRegion].spec.xAnch)
						/regionArray[curRegion].actualSubSize);
					firstY = (uint32_t) floor(
						(actorCommonArray[curActor].regionInterBound[curInterRegion][2] - regionArray[curRegion].spec.yAnch)
						/regionArray[curRegion].actualSubSize);
					firstZ = (uint32_t) floor(
						(actorCommonArray[curActor].regionInterBound[curInterRegion][4] - regionArray[curRegion].spec.zAnch)
						/regionArray[curRegion].actualSubSize);
					lastX = (uint32_t) ceil(
						(actorCommonArray[curActor].regionInterBound[curInterRegion][1] - regionArray[curRegion].spec.xAnch - regionArray[curRegion].actualSubSize)
						/regionArray[curRegion].actualSubSize);
					lastY = (uint32_t) ceil(
						(actorCommonArray[curActor].regionInterBound[curInterRegion][3] - regionArray[curRegion].spec.yAnch - regionArray[curRegion].actualSubSize)
						/regionArray[curRegion].actualSubSize);
					lastZ = (uint32_t) ceil(
						(actorCommonArray[curActor].regionInterBound[curInterRegion][5] - regionArray[curRegion].spec.zAnch - regionArray[curRegion].actualSubSize)
						/regionArray[curRegion].actualSubSize);
					
					// Find "Anchor" subvolume in intersection
					//firstSub = region2D[curRegion].firstID
					//	+ firstY*region2D[curRegion].spec.numX + firstX;
					
					for(curZ = firstZ; curZ <= lastZ; curZ++)
					{
						for(curY = firstY; curY <= lastY; curY++)
						{
							for(curX = firstX; curX <= lastX; curX++)
							{
								// Is subvolume valid?
								if(subID[curRegion][curX][curY][curZ] == UINT32_MAX)
									continue; // Subvolume space is within a child
								
								// Confirm that subvolume is within actor space
								curSubBound[0] = regionArray[curRegion].spec.xAnch
									+ curX*regionArray[curRegion].actualSubSize;
								curSubBound[1] = curSubBound[0]
									+ regionArray[curRegion].actualSubSize;
								curSubBound[2] = regionArray[curRegion].spec.yAnch
									+ curY*regionArray[curRegion].actualSubSize;
								curSubBound[3] = curSubBound[2]
									+ regionArray[curRegion].actualSubSize;
								curSubBound[4] = regionArray[curRegion].spec.zAnch
									+ curZ*regionArray[curRegion].actualSubSize;
								curSubBound[5] = curSubBound[4]
									+ regionArray[curRegion].actualSubSize;
								
								if(bBoundaryIntersect(
									actorCommonArray[curActor].spec.shape,
									actorCommonArray[curActor].spec.boundary,
									RECTANGULAR_BOX, curSubBound, 0.))
								{ // The subvolume does overlap the actor space
									actorCommonArray[curActor].numSub[curInterRegion]++;
								}
							}
						}
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
							regionArray[curRegion].firstID + curInterSub;
					}
				} else
				{	
					curInterSub = 0;
					for(curZ = firstZ; curZ <= lastZ; curZ++)
					{
						for(curY = firstY; curY <= lastY; curY++)
						{
							for(curX = firstX; curX <= lastX; curX++)
							{
								// Is subvolume valid?
								if(subID[curRegion][curX][curY][curZ] == UINT32_MAX)
									continue; // Subvolume space is within a child
								
								// Confirm that subvolume is within actor space
								curSubBound[0] = regionArray[curRegion].spec.xAnch
									+ curX*regionArray[curRegion].actualSubSize;
								curSubBound[1] = curSubBound[0]
									+ regionArray[curRegion].actualSubSize;
								curSubBound[2] = regionArray[curRegion].spec.yAnch
									+ curY*regionArray[curRegion].actualSubSize;
								curSubBound[3] = curSubBound[2]
									+ regionArray[curRegion].actualSubSize;
								curSubBound[4] = regionArray[curRegion].spec.zAnch
									+ curZ*regionArray[curRegion].actualSubSize;
								curSubBound[5] = curSubBound[4]
									+ regionArray[curRegion].actualSubSize;
								
								if(bBoundaryIntersect(
									actorCommonArray[curActor].spec.shape,
									actorCommonArray[curActor].spec.boundary,
									RECTANGULAR_BOX, curSubBound, 0.))
								{ // The subvolume does intersect the actor space
									curSub = subID[curRegion][curX][curY][curZ];
									actorCommonArray[curActor].subID[curInterRegion][curInterSub] = curSub;
									curInterSub++;
								}
								
							}
						}
					}
				}
				
				curInterRegion++;
			}
		}
		
		/*
		* Allocate simulation parameters
		* These parameters must be reset for each repeat
		*/
	}
	
	// Allocate memory for list of actors that record observations
	*actorRecordID =
		malloc((*numActorRecord)*sizeof(short));
	if(*numActorRecord > 0 && *actorRecordID == NULL){
		fprintf(stderr, "ERROR: Memory allocation for IDs of actors that will be recorded in the output file.\n");
		exit(EXIT_FAILURE);
	} else{ // There is at least one (passive) actor recording observations
		curActorRecord = 0;
		for(curActor = 0; curActor < NUM_ACTORS; curActor++)
		{
			if(!actorCommonArray[curActor].spec.bActive
				&& actorCommonArray[curActor].spec.bWrite)
				(*actorRecordID)[curActorRecord++] = curActor;
		}
	}
}

void allocateActorActivePassiveArray3D(const short NUM_ACTORS_ACTIVE,
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

void initializeActorActivePassive3D(const short NUM_ACTORS,
	const struct actorStruct3D actorCommonArray[],
	const unsigned short NUM_MOL_TYPES,
	const struct region regionArray[],
	const short NUM_REGIONS,
	const short NUM_ACTORS_ACTIVE,
	struct actorActiveStruct3D actorActiveArray[],
	const short NUM_ACTORS_PASSIVE,
	struct actorPassiveStruct3D actorPassiveArray[],
	uint32_t subCoorInd[][3])
{
	int i; // loop index
	short curActor, curActive, curPassive;
	short curRegion, curInterRegion;
	short curActorRecord;
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
		actorActiveArray[curActive].cumFracActorInSub = malloc(actorCommonArray[curActor].numRegion*sizeof(double *));
		
		initializeListRelease(&actorActiveArray[curActive].releaseList);
		initializeListData(&actorActiveArray[curActive].binaryData);
		
		if(actorActiveArray[curActive].cumFracActorInSub == NULL){
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
				|| actorCommonArray[curActor].bRegionInside[curInterRegion])
				continue; // Following only applies if region is mesoscopic and not entirely inside actor
						
			actorActiveArray[curActive].cumFracActorInSub[curInterRegion] =
				malloc(actorCommonArray[curActor].numSub[curInterRegion]
				*sizeof(double));
		
			if(actorActiveArray[curActive].cumFracActorInSub[curInterRegion]
				== NULL){
				fprintf(stderr,"ERROR: Memory allocation for structure members of active actor %u.\n", curActive);
				exit(EXIT_FAILURE);
			}
			
			for(curInterSub = 0;
				curInterSub < actorCommonArray[curActor].numSub[curInterRegion];
				curInterSub++)
			{					
				curSub = actorCommonArray[curActor].subID[curInterRegion][curInterSub];
				
				// Determine boundary of current subvolume
				curSubBound[0] = regionArray[curRegion].spec.xAnch
					+ regionArray[curRegion].actualSubSize*subCoorInd[curSub][0];
				curSubBound[1] = curSubBound[0]
					+ regionArray[curRegion].actualSubSize;
				curSubBound[2] = regionArray[curRegion].spec.yAnch
					+ regionArray[curRegion].actualSubSize*subCoorInd[curSub][1];
				curSubBound[3] = curSubBound[2]
					+ regionArray[curRegion].actualSubSize;
				curSubBound[4] = regionArray[curRegion].spec.zAnch
					+ regionArray[curRegion].actualSubSize*subCoorInd[curSub][2];
				curSubBound[5] = curSubBound[4]
					+ regionArray[curRegion].actualSubSize;
					
				// Determine intersection boundary of actor and subvolume
				intersectBoundary(actorCommonArray[curActor].spec.shape,
					actorCommonArray[curActor].spec.boundary,
					RECTANGULAR_BOX, curSubBound, curInterSubBound);
					
				if (curInterSub > 0)
				{
					actorActiveArray[curActive].cumFracActorInSub[curInterRegion][curInterSub] = 
						actorActiveArray[curActive].cumFracActorInSub[curInterRegion][curInterSub-1];
				} else
				{
					actorActiveArray[curActive].cumFracActorInSub[curInterRegion][0] = 0.;
				}
				actorActiveArray[curActive].cumFracActorInSub[curInterRegion][curInterSub] +=
					boundaryVolume(RECTANGULAR_BOX, curInterSubBound)
					/ actorCommonArray[curActor].regionInterArea[curInterRegion];
			}
		}
	}
	
	curActorRecord = 0;
	for(curPassive = 0; curPassive < NUM_ACTORS_PASSIVE; curPassive++)
	{
		curActor = actorPassiveArray[curPassive].actorID;
		actorPassiveArray[curPassive].numMolRecordID = 0;
		actorPassiveArray[curPassive].numMolRecordPosID = 0;
		actorPassiveArray[curPassive].recordID = SHRT_MAX;
		
		if (actorCommonArray[curActor].spec.bWrite)
		{
			actorPassiveArray[curPassive].recordID = curActorRecord++;
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
			
			// TODO: Could also skip the following if entire region is
			// inside of actor, if indicated by an additional member
			
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
				curSub = actorCommonArray[curActor].subID[curInterRegion][curInterSub];
				
				// Determine boundary of current subvolume
				curSubBound[0] = regionArray[curRegion].spec.xAnch
					+ regionArray[curRegion].actualSubSize*subCoorInd[curSub][0];
				curSubBound[1] = curSubBound[0]
					+ regionArray[curRegion].actualSubSize;
				curSubBound[2] = regionArray[curRegion].spec.yAnch
					+ regionArray[curRegion].actualSubSize*subCoorInd[curSub][1];
				curSubBound[3] = curSubBound[2]
					+ regionArray[curRegion].actualSubSize;
				curSubBound[4] = regionArray[curRegion].spec.zAnch
					+ regionArray[curRegion].actualSubSize*subCoorInd[curSub][2];
				curSubBound[5] = curSubBound[4]
					+ regionArray[curRegion].actualSubSize;
				
				// Determine intersection boundary of actor and subvolume
				if(actorCommonArray[curActor].bRegionInside[curInterRegion])
				{ // Entire subvolume must be inside actor
					for(i = 0; i < 6; i++)
						curInterSubBound[i] = curSubBound[i];
				} else
					intersectBoundary(actorCommonArray[curActor].spec.shape,
						actorCommonArray[curActor].spec.boundary,
						RECTANGULAR_BOX, curSubBound, curInterSubBound);
				
				if(actorPassiveArray[curPassive].bRecordMesoAnyPos[curInterRegion])
				{ // Need to record boundary of intersection of actor and subvolume
					for(i = 0; i < 6; i++)
						actorPassiveArray[curPassive].subInterBound[curInterRegion][curInterSub][i] =
							curInterSubBound[i];
				}
				
				if(actorCommonArray[curActor].bRegionInside[curInterRegion])
					actorPassiveArray[curPassive].fracSubInActor[curInterRegion][curInterSub] = 1.;
				else
					actorPassiveArray[curPassive].fracSubInActor[curInterRegion][curInterSub] =
						boundaryVolume(RECTANGULAR_BOX, curInterSubBound)
						/ regionArray[curRegion].actualSubSize
						/ regionArray[curRegion].actualSubSize
						/ regionArray[curRegion].actualSubSize;
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
void resetActors3D(const short NUM_ACTORS,
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
void deleteActor3D(const short NUM_ACTORS,
	struct actorStruct3D actorCommonArray[],
	const struct region regionArray[],
	const short NUM_ACTORS_ACTIVE,
	struct actorActiveStruct3D actorActiveArray[],
	const short NUM_ACTORS_PASSIVE,
	struct actorPassiveStruct3D actorPassiveArray[],
	short actorRecordID[])
{
	short curActor, curActive, curPassive;
	short curInterRegion, curRegion;
	
	free(actorRecordID);
	
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
					&& actorActiveArray[curActive].cumFracActorInSub[curInterRegion] != NULL)
					free(actorActiveArray[curActive].cumFracActorInSub[curInterRegion]);
			}
					
			if(actorActiveArray[curActive].cumFracActorInSub != NULL)
				free(actorActiveArray[curActive].cumFracActorInSub);
			
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
void newRelease3D(const struct actorStruct3D * actorCommon,
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
			if(mt_drand() < actorCommon->spec.probOne)
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
			// Bits are not random. Read from pre-determined list. TODO
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
				startTime = -log(mt_drand())/strength;
			} else
			{ // Emission times are deterministic. First emission will be at start of interval
				startTime = 0.;
			}
			molType = actorActive->molType[0]; // There is only one kind of molecule in CSK
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
		if(!addRelease(&actorActive->releaseList, strength, molType,
			curTime + startTime, curTime + endTime, frequency))
		{
			fprintf(stderr,"ERROR: Memory allocation for new active actor release.\n");
			exit(EXIT_FAILURE);
		}
		
		// Update time of next emission event
		findNextEmission3D(actorCommon, actorActive);
	}
	
	// Free memory of new data
	emptyListData(&newData);
}

// Find time and index of next release to have an emission
void findNextEmission3D(const struct actorStruct3D * actorCommon,
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

void fireEmission3D(const struct actorStruct3D * actorCommon,
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
		placeMolecules3D(actorCommon, actorActive, region, NUM_REGIONS,
			subvolArray, mesoSubArray, numMesoSub, 1ULL,
			curRelease->item.molType, NUM_MOL_TYPES, microMolListRecent,
			curRelease->item.nextTime, tMicro,
			heap_subvolID, heap_childID, b_heap_childValid);
		
		
		// Next emission time must be generated
		curRelease->item.nextTime +=
			-log(mt_drand())/curRelease->item.strength;
	} else
	{
		if(actorCommon->spec.bNumReleaseRand)
		{
			numNewMol = rd_poisson(curRelease->item.strength);
		} else
		{
			numNewMol = (uint64_t) curRelease->item.strength;
		}
		// We must place curRelease->item.strength molecules
		placeMolecules3D(actorCommon, actorActive, region, NUM_REGIONS,
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
	findNextEmission3D(actorCommon, actorActive);
}

// Place molecules for current emission
void placeMolecules3D(const struct actorStruct3D * actorCommon,
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
	short curRegion, curRegionInter;
	double uniRV;
	uint64_t curMolecule;
	//double tCur = curRelease->item.nextTime;
	
	if(actorCommon->numRegion > 1)
	{ // Molecules can end up in different regions. Place one at a time
		for(curMolecule = 0;
			curMolecule < numNewMol;
			curMolecule++)
		{
			uniRV = mt_drand();
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
				placeMoleculesInRegion3D(actorCommon, actorActive, region,
					curRegion, curRegionInter, NUM_REGIONS, subvolArray,
					mesoSubArray, numMesoSub, 1ULL, curMolType, NUM_MOL_TYPES,
					&microMolListRecent[curRegion][curMolType], tCur, tMicro,
					heap_subvolID, heap_childID, b_heap_childValid);
			}
		}
		
	} else
	{ // All molecules are going in the same region.
		curRegion = actorCommon->regionID[0];
		placeMoleculesInRegion3D(actorCommon, actorActive, region, curRegion,
			0, NUM_REGIONS, subvolArray,
			mesoSubArray, numMesoSub, numNewMol, curMolType, NUM_MOL_TYPES,
			&microMolListRecent[curRegion][curMolType], tCur, tMicro,
			heap_subvolID, heap_childID, b_heap_childValid);
	}
}

// Place molecules for current actor emission in specific region
void placeMoleculesInRegion3D(const struct actorStruct3D * actorCommon,
	const struct actorActiveStruct3D * actorActive,
	struct region region[],
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
		
	if(region[curRegion].spec.bMicro)
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
					actorCommon->regionInterBound[curRegionInter]);
				
				if(bPointInRegionNotChild(curRegion, region, point))
				{
					bNeedPoint = false;
					if(!addMoleculeRecent3D(microMolListRecent, point[0], point[1],
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
				uniRV = mt_drand();
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
				placeMoleculesInSub3D(region, subvolArray, mesoSubArray,
					numMesoSub, 1ULL, curMolType, curSub, NUM_MOL_TYPES,
					tCur, NUM_REGIONS, heap_subvolID, heap_childID, b_heap_childValid);
			}
		} else
		{ // All molecules are going in the same subvolume
			curSub = actorCommon->subID[curRegionInter][0];
			placeMoleculesInSub3D(region, subvolArray, mesoSubArray,
				numMesoSub, numNewMol, curMolType, curSub, NUM_MOL_TYPES,
				tCur, NUM_REGIONS, heap_subvolID,
				heap_childID, b_heap_childValid);
		}
	}
}

// Add one type of molecule in specific subvolume
void placeMoleculesInSub3D(struct region region[],
	struct subvolume3D subvolArray[],
	struct mesoSubvolume3D * mesoSubArray,
	const uint32_t numMesoSub,
	uint64_t numNewMol,
	const unsigned short curMolType,
	const uint32_t curSub,
	const unsigned short NUM_MOL_TYPES,
	double tCur,
	const short NUM_REGIONS,
	uint32_t * heap_subvolID,
	uint32_t (*heap_childID)[2],
	bool (*b_heap_childValid)[2])
{
	uint32_t curMeso;
	
	// Define which molecules are added where
	subvolArray[curSub].num_mol[curMolType] += numNewMol;
	curMeso = subvolArray[curSub].mesoID;
	
	// Place molecule(s) and update heap for subvolumes
	updateMesoSub3D(curSub, false, (uint64_t []){numNewMol},
		(bool []){true}, curMolType,
		true, numMesoSub, mesoSubArray, subvolArray, tCur,
		NUM_REGIONS, NUM_MOL_TYPES, region);
	heapMesoUpdate3D(numMesoSub, mesoSubArray, heap_subvolID,
		mesoSubArray[curMeso].heapID, heap_childID, b_heap_childValid);
}