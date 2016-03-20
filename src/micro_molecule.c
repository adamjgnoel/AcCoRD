/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * micro_molecule.c - 	linked list of individual molecules in same
 * 						microscopic region
 *
 * Last revised for AcCoRD LATEST_RELEASE
 *
 * Revision history:
 *
 * Revision LATEST_RELEASE
 * - added surface reactions, including membrane transitions
 * - corrected distance to end point when a molecule is "pushed" into a neighboring
 * region
 * - added fail check to while loop when a molecule is "pushed" into a 
 * neighboring region. Error will display if we did not end up in specified
 * region or one of its children.
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

#include "micro_molecule.h"

// Local Function Prototypes

static void traverse(const ListMol3D * p_list, void (* p_fun)(ItemMol3D item));

static void traverseRecent(const ListMolRecent3D * p_list, void (* p_fun)(ItemMolRecent3D item));

static bool addItem(ItemMol3D item, ListMol3D * p_list);

static bool addItemRecent(ItemMolRecent3D item, ListMolRecent3D * p_list);

static void removeItem(NodeMol3D * prevNode, NodeMol3D * curNode);

static void removeItemRecent(NodeMolRecent3D * prevNode, NodeMolRecent3D * curNode);

static void copyToNode(ItemMol3D item, NodeMol3D * p_node);

static void copyToNodeRecent(ItemMolRecent3D item, NodeMolRecent3D * p_node);

// Specific Definitions

// Create new molecule at specified coordinates
bool addMolecule(ListMol3D * p_list, double x, double y, double z)
{
	ItemMol3D new_molecule = {x,y,z, true};
	return addItem(new_molecule, p_list);
}

// Create new molecule at specified coordinates
bool addMoleculeRecent(ListMolRecent3D * p_list, double x, double y, double z, double dt_partial)
{
	ItemMolRecent3D new_molecule = {x,y,z, dt_partial};
	return addItemRecent(new_molecule, p_list);
}

// Move one molecule to the specified coordinates
void moveMolecule(ItemMol3D * molecule, double x, double y, double z)
{
	molecule->x = x;
	molecule->y = y;
	molecule->z = z;
}

// Move one molecule to the specified coordinates
void moveMoleculeRecent(ItemMolRecent3D * molecule, double x, double y, double z)
{
	molecule->x = x;
	molecule->y = y;
	molecule->z = z;
}

// Move ALL molecules in the list by the same standard deviation
void diffuseMolecules(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	ListMol3D p_list[NUM_REGIONS][NUM_MOL_TYPES],
	ListMolRecent3D p_listRecent[NUM_REGIONS][NUM_MOL_TYPES],
	const struct region regionArray[],
	struct mesoSubvolume3D mesoSubArray[],
	struct subvolume3D subvolArray[],
	double sigma[NUM_REGIONS][NUM_MOL_TYPES],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES])
{
	NodeMol3D * curNode, * prevNode, * nextNode;
	NodeMolRecent3D * curNodeR;
	
	double oldPoint[3];
	double newPoint[3];
	short curRegion, curType, newRegion, newType, transRegion;
	uint32_t newSub, curBoundSub;
	bool bInRegion;
	int curMol = 0;
	
	bool bReaction;
	unsigned short curRxn, curProd;
	
	// Indicate that every microscopic molecule in a "normal" list
	// needs to be moved.
	// We do this to avoid moving a molecule more than once if it is moved
	// to a different region.
	for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
	{
		for(curType = 0; curType < NUM_MOL_TYPES; curType++)
		{
			if(isListMol3DEmpty(&p_list[curRegion][curType])
				|| sigma[curRegion][curType] == 0.)
				continue; // No need to validate an empty list of molecules or ones that can't move
			
			curNode = p_list[curRegion][curType];
			
			while(curNode != NULL)
			{
				curNode->item.bNeedUpdate = true;
				curNode = curNode->next;
			}
		}
	}
	
	// Diffuse molecule in "regular" molecule lists
	for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
	{
		for(curType = 0; curType < NUM_MOL_TYPES; curType++)
		{
			if(isListMol3DEmpty(&p_list[curRegion][curType])
				|| sigma[curRegion][curType] == 0.)
				continue; // No need to validate an empty list of molecules
			
			curNode = p_list[curRegion][curType];
			prevNode = NULL;
			
			while(curNode != NULL)
			{
				nextNode = curNode->next;
				
				if(curNode->item.bNeedUpdate)
				{ // Molecule needs to be moved
					curNode->item.bNeedUpdate = false;
				
					oldPoint[0] = curNode->item.x;
					oldPoint[1] = curNode->item.y;
					oldPoint[2] = curNode->item.z;
					
					// Diffuse molecule
					diffuseOneMolecule(&curNode->item, sigma[curRegion][curType]);
					
					newPoint[0] = curNode->item.x;
					newPoint[1] = curNode->item.y;
					newPoint[2] = curNode->item.z;
								
					bReaction = false;
					if(validateMolecule(newPoint, oldPoint, NUM_REGIONS, curRegion,
						&newRegion, &transRegion, regionArray, curType, &bReaction, &curRxn))
					{
						// Molecule is still within region and no further action is required
						prevNode = curNode;
					} else
					{
						// newPoint tells us where to place a molecule in newRegion
						if(newRegion == curRegion)
						{ //just update molecule location
							moveMolecule(&curNode->item, newPoint[0], newPoint[1], newPoint[2]);
							prevNode = curNode;
						} else
						{
							// Molecule is now in a different region
							if(regionArray[newRegion].spec.bMicro)
							{ // New region is microscopic. Move to appropriate list
								
								if(bReaction)
								{ // We need to fire the corresponding reaction curRxn
									if(regionArray[newRegion].numRxnProducts[curRxn] > 0)
									{
										for(curProd = 0;
											curProd < regionArray[newRegion].numRxnProducts[curRxn];
											curProd++)
										{
											// Add the (curProd)th product to the corresponding molecule list
											if(!addMolecule(
												&p_list[newRegion][regionArray[newRegion].productID[curRxn][curProd]],
												newPoint[0], newPoint[1], newPoint[2]))
											{ // Creation of molecule failed
												fprintf(stderr, "ERROR: Memory allocation to create molecule of type %u from reaction %u.\n",
												regionArray[newRegion].productID[curRxn][curProd], curRxn);
												exit(EXIT_FAILURE);						
											}
											// Indicate that product molecule doesn't need to be
											// moved again
											p_list[newRegion][regionArray[newRegion].productID[curRxn][curProd]]->item.bNeedUpdate = false;
										}
									}									
								} else if(!addMolecule(&p_list[newRegion][curType],
									newPoint[0], newPoint[1], newPoint[2]))
								{
									fprintf(stderr, "ERROR: Memory allocation to move molecule between microscopic regions %u and %u.\n", curRegion, newRegion);
									exit(EXIT_FAILURE);
								} else
								{
									// Indicate that molecule doesn't need to be moved again
									p_list[newRegion][curType]->item.bNeedUpdate = false;
								}
							} else
							{ // New region is mesoscopic. Find nearest subvolume to new point
								newSub = findNearestSub(newRegion, regionArray,
									transRegion, newPoint[0], newPoint[1], newPoint[2]);
								subvolArray[newSub].num_mol[curType]++;
								 // TODO: Keep track of subvolumes that need updated propensities
								 // mesoID is subvolArray[newSub].mesoID
								 //
								curBoundSub = 0;
								while(regionArray[newRegion].neighID[transRegion][curBoundSub]
									!= newSub)
								{
									curBoundSub++;
								}
								regionArray[newRegion].bNeedUpdate[transRegion][curBoundSub] = true;
								regionArray[newRegion].numMolFromMicro[transRegion][curBoundSub][curType]++;
							}
							
							// Remove molecule from current region molecule list
							removeItem(prevNode,curNode);
							
							if(prevNode == NULL)
							{ // We removed first molecule in list.
							  // nextNode is now the start of the list
							  // (i.e., we must update pointer to list)
							  p_list[curRegion][curType] = nextNode;
							}
						}
					}					
				} else
				{
					prevNode = curNode;
				}
				curNode = nextNode;
			}
		}
	}
	
	// Diffuse molecule in "recent" molecule lists
	// These molecules all get added to the "normal" lists (if region is microscopic)
	for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
	{
		for(curType = 0; curType < NUM_MOL_TYPES; curType++)
		{
			if(isListMol3DRecentEmpty(&p_listRecent[curRegion][curType]))
				continue; // No need to validate an empty list of molecules
			
			curNodeR = p_listRecent[curRegion][curType];
			curMol = 0;
			while(curNodeR != NULL)
			{				
				oldPoint[0] = curNodeR->item.x;
				oldPoint[1] = curNodeR->item.y;
				oldPoint[2] = curNodeR->item.z;
				
				// Diffuse molecule
				diffuseOneMoleculeRecent(&curNodeR->item, DIFF_COEF[curRegion][curType]);
				
				newPoint[0] = curNodeR->item.x;
				newPoint[1] = curNodeR->item.y;
				newPoint[2] = curNodeR->item.z;
				
				// Once molecule is validated, we can proceed directly to transferring
				// it to the relevant "normal" list and remove it from this list				
				bReaction = false;
				validateMolecule(newPoint, oldPoint, NUM_REGIONS, curRegion,
					&newRegion, &transRegion, regionArray, curType, &bReaction, &curRxn);
					
				if(regionArray[newRegion].spec.bMicro)
				{ // Region is microscopic. Move to appropriate list
					if(bReaction)
					{ // We need to fire the corresponding reaction curRxn
						if(regionArray[newRegion].numRxnProducts[curRxn] > 0)
						{
							for(curProd = 0;
								curProd < regionArray[newRegion].numRxnProducts[curRxn];
								curProd++)
							{
								// Add the (curProd)th product to the corresponding molecule list
								if(!addMolecule(
									&p_list[newRegion][regionArray[newRegion].productID[curRxn][curProd]],
									newPoint[0], newPoint[1], newPoint[2]))
								{ // Creation of molecule failed
									fprintf(stderr, "ERROR: Memory allocation to create molecule of type %u from reaction %u.\n",
									regionArray[newRegion].productID[curRxn][curProd], curRxn);
									exit(EXIT_FAILURE);						
								}
								// Indicate that product molecule doesn't need to be
								// moved again
								p_list[newRegion][regionArray[newRegion].productID[curRxn][curProd]]->item.bNeedUpdate = false;
							}
						}									
					} else if(!addMolecule(&p_list[newRegion][curType],
						newPoint[0], newPoint[1], newPoint[2]))
					{
						fprintf(stderr, "ERROR: Memory allocation to move molecule between recent molecule list of region %u and list of region %u.\n", curRegion, newRegion);
						exit(EXIT_FAILURE);
					}
				} else
				{ // New region is mesoscopic. Find nearest subvolume to new point
					newSub = findNearestSub(newRegion, regionArray,
						transRegion, newPoint[0], newPoint[1], newPoint[2]);
					subvolArray[newSub].num_mol[curType]++;
					 // TODO: Keep track of subvolumes that need updated propensities
					 // mesoID is subvolArray[newSub].mesoID
					 //
					curBoundSub = 0;
					while(regionArray[newRegion].neighID[transRegion][curBoundSub]
						!= newSub)
					{
						curBoundSub++;
					}
					regionArray[newRegion].bNeedUpdate[transRegion][curBoundSub] = true;
					regionArray[newRegion].numMolFromMicro[transRegion][curBoundSub][curType]++;
				}
			
				curNodeR = curNodeR->next;
				curMol++;
			}
		}
	}
	
	// Empty recent lists
	for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
	{
		for(curType = 0; curType < NUM_MOL_TYPES; curType++)
		{
			if(isListMol3DRecentEmpty(&p_listRecent[curRegion][curType]))
				continue; // No need to validate an empty list of molecules
			
			emptyListMol3DRecent(&p_listRecent[curRegion][curType]);
			initializeListMolRecent(&p_listRecent[curRegion][curType]);
		}
	}
}

// Move one molecule by some standard deviation
void diffuseOneMolecule(ItemMol3D * molecule, double sigma)
{
	molecule->x = rd_normal(molecule->x, sigma);
	molecule->y = rd_normal(molecule->y, sigma);
	molecule->z = rd_normal(molecule->z, sigma);
}

// Move one molecule by some standard deviation
void diffuseOneMoleculeRecent(ItemMolRecent3D * molecule, double DIFF_COEF)
{
	double sigma = sqrt(2*molecule->dt_partial*DIFF_COEF);
	molecule->x = rd_normal(molecule->x, sigma);
	molecule->y = rd_normal(molecule->y, sigma);
	molecule->z = rd_normal(molecule->z, sigma);
}

// Check first order reactions for all molecules in list
void rxnFirstOrder(ListMol3D * p_list,
	const struct region regionArray,
	unsigned short curMolType,
	const unsigned short NUM_MOL_TYPES,
	ListMolRecent3D pRecentList[NUM_MOL_TYPES])
{
	NodeMol3D * curNode = *p_list;
	NodeMol3D * prevNode = NULL;
	NodeMol3D * nextNode;
	double curRand;
	double curTime;
	int i; // Reaction loop index
	unsigned short curProd;
	unsigned short curRxn;
	bool bRemove;
	
	while(curNode != NULL)
	{
		nextNode = curNode->next;
		bRemove = false;
		
		// Generate RV to see which/whether first order reaction occurred
		curRand = mt_drand();
		
		for(i = 0; i < regionArray.numFirstCurReactant[curMolType]; i++)
		{
			if(curRand < regionArray.uniCumProb[curMolType][i])
			{
				// Reaction i took place
				curRxn = regionArray.firstRxnID[curMolType][i];
				
				// Generate products (if necessary)
				// TODO: May need to consider associating a separation distance with
				// reactions that have more than one product.
				if(regionArray.numRxnProducts[curRxn] > 0)
				{
					// If there are any products, then we need to generate a
					// corresponding reaction time so we know when the products
					// were generated
					curTime =
						-log((1.-regionArray.minRxnTimeRV[curMolType])*mt_drand()
						+ regionArray.minRxnTimeRV[curMolType])
						/ regionArray.uniSumRate[curMolType];
					for(curProd = 0;
						curProd < regionArray.numRxnProducts[curRxn];
						curProd++)
					{
						// Add the (curProd)th product to the corresponding molecule list
						if(!addMoleculeRecent(
							&pRecentList[regionArray.productID[curRxn][curProd]],
							curNode->item.x, curNode->item.y, curNode->item.z,
							regionArray.spec.dt-curTime))
						{ // Creation of molecule failed
							fprintf(stderr, "ERROR: Memory allocation to create molecule of type %u from reaction %u.\n", regionArray.productID[curRxn][curProd], curRxn);
							exit(EXIT_FAILURE);						
						}
					}
				}
				
				// Remove current molecule from list
				removeItem(prevNode, curNode);
				bRemove = true;
				break; // exit for loop
			}
		}
		
		if(prevNode == NULL && bRemove)
		{	// prevNode does not change, but we removed first molecule in list.
			// nextNode is now the start of the list
			// (i.e., we must update pointer to list)
			*p_list = nextNode;
		} else if(!bRemove){
			prevNode = curNode;
		}
		curNode = nextNode;
	}
}

// Check first order reactions for all molecules in recent molecule list
// The key differences between this function and rxnFirstOrder are as follows:
// - counts the number of product molecules
// - IF this is called with bCheckCount == true, then only the specified number
// of molecules are checked. This is used when this function is called repeatedly
// to deal with molecules that were created in previous calls in the same time step
void rxnFirstOrderRecent(const unsigned short NUM_MOL_TYPES,
	ListMolRecent3D pRecentList[NUM_MOL_TYPES],
	const struct region regionArray,
	unsigned short curMolType,
	bool bCheckCount,
	uint32_t numMolCheck[NUM_MOL_TYPES])
{
	NodeMolRecent3D * curNode = pRecentList[curMolType];
	NodeMolRecent3D * prevNode = NULL;
	NodeMolRecent3D * nextNode;
	double curRand;
	double curTime;
	int i; // Reaction loop index
	unsigned short curProd;
	unsigned short curRxn;
	bool bRemove;
	double uniCumProb[regionArray.numFirstCurReactant[curMolType]];
	double minRxnTimeRV;
	uint32_t numMolOrig = numMolCheck[curMolType];
	uint32_t curMol = 0;
	unsigned short prodID;
	
	while(curNode != NULL)
	{
		if(bCheckCount && ++curMol > numMolOrig)
		{ // We only need to check the first numMolOrig molecules in this list
		  // and we have already done so.
			break;
		}
		
		nextNode = curNode->next;
		bRemove = false;
		
		// Generate RV to see which/whether first order reaction occurred
		curRand = mt_drand();
		
		// Need to generate probabilities for current molecule based on the partial
		// time step size
		
		for(i = 0; i < regionArray.numFirstCurReactant[curMolType]; i++)
		{
			curRxn = regionArray.firstRxnID[curMolType][i];
			
			if(i == 0)
			{
				uniCumProb[i] = regionArray.uniRelativeRate[curMolType][i]
					*(1.-exp(-curNode->item.dt_partial*regionArray.uniSumRate[curMolType]));
			} else
			{
				uniCumProb[i] = uniCumProb[i-1]
					+ regionArray.uniRelativeRate[curMolType][i]
					*(1.-exp(-curNode->item.dt_partial*regionArray.uniSumRate[curMolType]));
			}
			
			if(curRand < uniCumProb[i])
			{
				// Generate products (if necessary)
				// TODO: May need to consider associating a separation distance with
				// reactions that have more than one product.
				if(regionArray.numRxnProducts[curRxn] > 0)
				{
					// minRxnTimeRV = exp(-curNode->item.dt_partial*regionArray.rxnRate[curRxn]);
					minRxnTimeRV = exp(-curNode->item.dt_partial*regionArray.uniSumRate[curMolType]);
					curTime =
						-log((1.-minRxnTimeRV)*mt_drand()
						+ minRxnTimeRV)
						/ regionArray.uniSumRate[curMolType];
					for(curProd = 0;
						curProd < regionArray.numRxnProducts[curRxn];
						curProd++)
					{
						prodID = regionArray.productID[curRxn][curProd];
						// Add the (curProd)th product to the corresponding molecule list
						if(!addMoleculeRecent(
							&pRecentList[prodID],
							curNode->item.x, curNode->item.y, curNode->item.z,
							curNode->item.dt_partial-curTime))
						{ // Creation of molecule failed
							fprintf(stderr, "ERROR: Memory allocation to create molecule of type %u from reaction %u.\n", regionArray.productID[curRxn][curProd], curRxn);;
							exit(EXIT_FAILURE);						
						}
						if(prevNode == NULL
							&& prodID == curMolType)
						{ // We've just added a molecule to the start of the current
							// molecule list. Update prevNode
							prevNode = pRecentList[curMolType];
						}
						// Update the number of molecules that need to be checked
						// in the next round
						// NOTE: We don't need to check more reactions if this product
						// isn't the reactant for any first-order reactions, but it's
						// somewhat inefficient to check here
						if(bCheckCount || prodID <= curMolType)
						{
							numMolCheck[prodID]++;
						}
					}
				}
				
				// Remove current molecule from list
				removeItemRecent(prevNode, curNode);
				bRemove = true;
				break; // exit for loop
			}
		}
		
		if(prevNode == NULL && bRemove)
		{	// prevNode does not change, but we removed first molecule in list.
			// nextNode is now the start of the list
			// (i.e., we must update pointer to list)
			pRecentList[curMolType] = nextNode;
		} else if(!bRemove){
			prevNode = curNode;
		}
		curNode = nextNode;
	}
	
	if(bCheckCount)
	{ // Correct the number of molecules to check in next round
		numMolCheck[curMolType] -= numMolOrig;
	}
}

// Empty recent list and add corresponding molecules to "normal" list
void transferMolecules(ListMolRecent3D * molListRecent, ListMol3D * molList)
{
	NodeMolRecent3D * p_node = *molListRecent;
	
	// Copy list of molecules to normal list
	while(p_node != NULL)
	{
		if(!addMolecule(molList, p_node->item.x, p_node->item.y, p_node->item.z))
		{
			// Creation of molecule failed
			fprintf(stderr, "ERROR: Memory allocation to create molecule when transferring from recent list to regular list.\n");
			exit(EXIT_FAILURE);
		}
		p_node = p_node->next;
	}
	
	// Empty the recent list and re-initialize
	emptyListMol3DRecent(molListRecent);
	initializeListMolRecent(molListRecent);
}

bool validateMolecule(double newPoint[3],
	double oldPoint[3],
	const short NUM_REGIONS,
	const short curRegion,
	short * newRegion,
	short * transRegion,
	const struct region regionArray[],
	short molType,
	bool * bReaction,
	unsigned short * curRxn)
{
	double trajLine[3];
	double lineLength;
	
	// Default for newRegion is curRegion
	*newRegion = curRegion;
	
	short trajRegion = curRegion; // Regions passed through by the molecule's trajectory
	
	if(bPointInRegionNotChild(curRegion, regionArray, newPoint)
		&& regionArray[curRegion].numChildren < 1)
	{ // This is simplest case. No region boundary interactions
		return true;
	} else
	{ // Molecule may have left region's outer boundary or went through a child region
		
		// Define trajectory vector
		defineLine(oldPoint, newPoint, trajLine, &lineLength);
		
		return followMolecule(oldPoint, newPoint, trajLine,
			lineLength, curRegion,
			newRegion, transRegion, regionArray, molType, bReaction, curRxn, 0);
	}
}

// Recursively follow a molecule's path through region boundaries from its diffusion
// start and end points
// Return false if molecule path had to be changed
bool followMolecule(const double startPoint[3],
	double endPoint[3],
	double lineVector[3],
	double lineLength,
	const short startRegion,
	short * endRegion,
	short * transRegion,
	const struct region regionArray[],
	short molType,
	bool * bReaction,
	unsigned short * curRxn,
	unsigned int depth)
{
	short curNeigh, curRegion, closestNormal;
	short curFace, nearestFace;
	double minDist, curDist;
	double curIntersectPoint[3];
	double nearestIntersectPoint[3];
	double newEndPoint[3];
	double newLineVector[3];
	bool bCurIntersect = false;
	
	double curRand;
	int i;
	unsigned short curProd;
	* bReaction = false;
	bool bContinue = false;
	bool bReflect = false;
	
	double pushFrac = SUB_ADJ_RESOLUTION;
	
	// First check all neighbor regions to see which (if any) are intersected first
	minDist = INFINITY;
	* endRegion = SHRT_MAX;
	closestNormal = SHRT_MAX;
	for(curNeigh = 0; curNeigh < regionArray[startRegion].numRegionNeigh; curNeigh++)
	{
		curRegion = regionArray[startRegion].regionNeighID[curNeigh];
		
		switch(regionArray[startRegion].regionNeighDir[curRegion])
		{
			case PARENT:
				// curRegion is startRegion's parent
				// Need to check against all faces of startRegion
				bCurIntersect = bLineHitRegion(startPoint, lineVector, lineLength,
					startRegion, curRegion, regionArray,
					&curFace, &curDist, curIntersectPoint)
					&& !bSharedBoundary(startRegion, curRegion, regionArray, curFace);
				break;
			case CHILD:
				// startRegion is curRegion's parent
				// Need to check against all faces of curRegion
				bCurIntersect = bLineHitRegion(startPoint, lineVector, lineLength,
					startRegion, curRegion, regionArray,
					&curFace, &curDist, curIntersectPoint);
				break;
			default:
				// Regions are adjacent. Only check adjacent plane
				curFace = regionArray[startRegion].regionNeighDir[curRegion];
				bCurIntersect = bLineHitInfinitePlane(startPoint, lineVector,
					lineLength, RECTANGULAR_BOX,
					regionArray[startRegion].boundRegionFaceCoor[curRegion][0],
					curFace, false, &curDist, curIntersectPoint)
					&& bPointOnFace(curIntersectPoint, RECTANGULAR_BOX,
					regionArray[startRegion].boundRegionFaceCoor[curRegion][0],
					curFace);
				break;
		}
		
		if(bCurIntersect && curDist < minDist+regionArray[startRegion].subResolution)
		{ // There is intersection with this face and it is ~closest so far
			if(regionArray[curRegion].spec.type == REGION_NORMAL)
			{
				closestNormal = curRegion;
			}
			// Check region priority
			if(* endRegion < SHRT_MAX
				&& regionArray[*endRegion].spec.type != REGION_NORMAL
				&& minDist-curDist < regionArray[startRegion].subResolution)
			{ // curRegion is not effectively any closer than a surface region, so do not override
				continue;
			}
			minDist = curDist;
			* endRegion = curRegion;
			nearestFace = curFace;
			nearestIntersectPoint[0] = curIntersectPoint[0];
			nearestIntersectPoint[1] = curIntersectPoint[1];
			nearestIntersectPoint[2] = curIntersectPoint[2];
		}
	}
	
	if(* endRegion < SHRT_MAX)
	{
		// Molecule has entered another region
		// "Lock" location at exact boundary
		lockPointToRegion(nearestIntersectPoint, startRegion, *endRegion, regionArray,
				nearestFace);		
		
		if(regionArray[*endRegion].spec.type != REGION_NORMAL)
		{ // Closest region is some kind of surface
			// Check whether we are actually in child of surface
			if(!bPointInRegionOrChild(*endRegion, regionArray, nearestIntersectPoint, endRegion))
			{ // Something went wrong here because we are not actually
				// in the endRegion region (or its children)
				fprintf(stderr, "ERROR: Invalid transition from region %u (Label: \"%s\") to region %u (Label: \"%s\"). Leaving molecule in \"%s\"\n",
					startRegion, regionArray[startRegion].spec.label,
					*endRegion, regionArray[*endRegion].spec.label,
					regionArray[startRegion].spec.label);
				endPoint[0] = startPoint[0];
				endPoint[1] = startPoint[1];
				endPoint[2] = startPoint[2];
				*endRegion = startRegion;
				return false;
			}
			
			if(regionArray[*endRegion].numFirstRxn > 0)
			{ // There is at least 1 possible reaction
				curRand = mt_drand();
				for(i = 0; i < regionArray[*endRegion].numFirstCurReactant[molType]; i++)
				{
					if(curRand < regionArray[*endRegion].uniCumProb[molType][i])
					{
						// Reaction i took place
						*curRxn = regionArray[*endRegion].firstRxnID[molType][i];
						*bReaction = true;
						continue;
					}
				}
			}
			if(*bReaction)
			{
				switch(regionArray[*endRegion].spec.surfaceType)
				{
					case SURFACE_INNER:
					case SURFACE_OUTER:
						// Molecule needs to stick
						bContinue = false;
						break;
					case SURFACE_MEMBRANE:
						// Molecule can continue into the closest normal region
						bContinue = true;
						*endRegion = closestNormal;
						break;
				}
			} else
			{
				// Need to reflect
				bReflect = true;
			}
		} else
		{
			bContinue = true;
		}
		
		if(bReflect)
		{
			// Reflect the point back into its starting region
			*endRegion = startRegion;
		} else
		{
			// Assume that transition is valid and proceed to follow point into neighbor
			
			if(bContinue)
			{
				lineLength -= minDist;
			
				// "Push" slightly into region to confirm which region we are really in
				// Could be child of *endRegion
				pushPoint(nearestIntersectPoint, curIntersectPoint, lineLength*pushFrac, lineVector);		
				while(!bPointInBoundary(curIntersectPoint, regionArray[*endRegion].spec.shape,
					regionArray[*endRegion].boundary) && pushFrac > 0.)
				{
					pushFrac *= 0.1;
					pushPoint(nearestIntersectPoint, curIntersectPoint, lineLength*pushFrac, lineVector);
				}
				
				// Determine region that we should actually be in (could be child of *endRegion)
				if(!bPointInRegionOrChild(*endRegion, regionArray, curIntersectPoint, endRegion))
				{ // Something went wrong here because we are not actually
					// in the endRegion region (or its children)
					fprintf(stderr, "ERROR: Invalid transition from region %u (Label: \"%s\") to region %u (Label: \"%s\"). Leaving molecule in \"%s\"\n",
						startRegion, regionArray[startRegion].spec.label,
						*endRegion, regionArray[*endRegion].spec.label,
						regionArray[startRegion].spec.label);
					endPoint[0] = startPoint[0];
					endPoint[1] = startPoint[1];
					endPoint[2] = startPoint[2];
					*endRegion = startRegion;
					return false;
				}
					
				if(!regionArray[*endRegion].spec.bMicro)
				{
					// Molecule is entering mesoscopic region. Stop following here.
					endPoint[0] = curIntersectPoint[0];
					endPoint[1] = curIntersectPoint[1];
					endPoint[2] = curIntersectPoint[2];
					* transRegion = startRegion; // Indicate from which micro region we came from
					return false;
				}			
				lineLength -= lineLength*pushFrac; // Correct line length for having been pushed	
				return followMolecule(curIntersectPoint, endPoint, lineVector,
					lineLength, *endRegion, endRegion, transRegion, regionArray, molType,
					bReaction, curRxn, depth+1)
					&& startRegion == *endRegion;
			} else
			{
				endPoint[0] = nearestIntersectPoint[0];
				endPoint[1] = nearestIntersectPoint[1];
				endPoint[2] = nearestIntersectPoint[2];
				return false;
			}
		}
	} else if(bPointInRegionNotChild(startRegion, regionArray, endPoint))
	{
		// Point is still in current region so diffusion is valid with no further modification
		*endRegion = startRegion;
		return true;
	}
	
	// Point needs to be reflected off of its own region boundary	
	if(!reflectPoint(startPoint, lineVector, lineLength, endPoint, newEndPoint,
		nearestIntersectPoint, &nearestFace, regionArray[startRegion].spec.shape,
		regionArray[startRegion].boundary, true))
	{ // Reflection failed. Leave point at nearestIntersectPoint
		*endRegion = startRegion;
		endPoint[0] = nearestIntersectPoint[0];
		endPoint[1] = nearestIntersectPoint[1];
		endPoint[2] = nearestIntersectPoint[2];
		return false;		
	}
	// Lock point to actual boundary only if actual reflection occurred
	lockPointToRegion(nearestIntersectPoint, startRegion, startRegion, regionArray,
		nearestFace);
	// Follow point past reflection. Define new line unit vector
	defineLine(nearestIntersectPoint, newEndPoint, lineVector, &lineLength);
	
	followMolecule(nearestIntersectPoint, newEndPoint, lineVector,
			lineLength, startRegion, endRegion, transRegion, regionArray,
			molType, bReaction, curRxn, depth+1);
	
	endPoint[0] = newEndPoint[0];
	endPoint[1] = newEndPoint[1];
	endPoint[2] = newEndPoint[2];
	return false;
}

// Count number of molecules inside observer
uint64_t countMolecules(ListMol3D * p_list,
	int obsType,
	double boundary[])
{
	NodeMol3D * p_node = *p_list;
	uint64_t curCount = 0ULL;
	
	while(p_node != NULL)
	{
		if (isMoleculeObserved(&p_node->item, obsType, boundary))
			curCount++;
		p_node = p_node->next;
	}
	return curCount;
}

// Count number of molecules inside observer
uint64_t countMoleculesRecent(ListMolRecent3D * p_list,
	int obsType,
	double boundary[])
{
	NodeMolRecent3D * p_node = *p_list;
	uint64_t curCount = 0ULL;
	
	while(p_node != NULL)
	{
		if (isMoleculeObservedRecent(&p_node->item, obsType, boundary))
			curCount++;
		p_node = p_node->next;
	}
	return curCount;
}

// Count molecules within boundary and add 
uint64_t recordMolecules(ListMol3D * p_list,
	ListMol3D * recordList,
	int obsType,
	double boundary[],
	bool bRecordPos)
{
	NodeMol3D * p_node = *p_list;
	uint64_t curCount = 0ULL;
	
	while(p_node != NULL)
	{
		if (isMoleculeObserved(&p_node->item, obsType, boundary))
		{
			curCount++;
			if(bRecordPos && !addItem(p_node->item, recordList))
			{
				fprintf(stderr,"\nERROR: Memory allocation for recording molecule positions.\n");
				exit(EXIT_FAILURE);
			}
		}
		p_node = p_node->next;
	}
	return curCount;
}

// Count molecules within boundary and add 
uint64_t recordMoleculesRecent(ListMolRecent3D * p_list,
	ListMol3D * recordList,
	int obsType,
	double boundary[],
	bool bRecordPos)
{
	NodeMolRecent3D * p_node = *p_list;
	uint64_t curCount = 0ULL;
	
	while(p_node != NULL)
	{
		if (isMoleculeObservedRecent(&p_node->item, obsType, boundary))
		{
			curCount++;
			if(bRecordPos && !addMolecule(recordList, p_node->item.x, p_node->item.y, p_node->item.z))
			{
				fprintf(stderr,"\nERROR: Memory allocation for recording molecule positions.\n");
				exit(EXIT_FAILURE);
			}
		}
		p_node = p_node->next;
	}
	return curCount;
}

// Is molecule inside the specified boundary?
bool isMoleculeObserved(ItemMol3D * molecule, int obsType, double boundary[])
{
	switch (obsType)
	{
		case CIRCLE: // Observer is circular
			// boundary has format {r^2, x, y}
			return false; // TODO: Complete this case
		case RECTANGLE:
		case RECTANGULAR_BOX: // Observer is a box
			// boundary has format {xMin, xMax, yMin, yMax, zMin, zMax}
			return (molecule->x >= boundary[0] && molecule->x <= boundary[1] &&
				molecule->y >= boundary[2] && molecule->y <= boundary[3] &&
				molecule->z >= boundary[4] && molecule->z <= boundary[5]);
		case SPHERE:
			// boundary has format {x, y, z, r, r^2}
			return squareDBL(molecule->x - boundary[0])
				+ squareDBL(molecule->y - boundary[1])
				+ squareDBL(molecule->z - boundary[2])
				<= boundary[4];
		default:
			return false;
	}
}

// Is molecule inside the specified boundary?
bool isMoleculeObservedRecent(ItemMolRecent3D * molecule, int obsType, double boundary[])
{
	switch (obsType)
	{
		case CIRCLE: // Observer is circular
			// boundary has format {r^2, x, y}
			return false; // TODO: Complete this case
		case RECTANGLE:
		case RECTANGULAR_BOX: // Observer is a box
			// boundary has format {xMin, xMax, yMin, yMax, zMin, zMax}
			return (molecule->x >= boundary[0] && molecule->x <= boundary[1] &&
				molecule->y >= boundary[2] && molecule->y <= boundary[3] &&
				molecule->z >= boundary[4] && molecule->z <= boundary[5]);
		case SPHERE:
			// boundary has format {x, y, z, r, r^2}
			return squareDBL(molecule->x - boundary[0])
				+ squareDBL(molecule->y - boundary[1])
				+ squareDBL(molecule->z - boundary[2])
				<= boundary[4];
		default:
			return false;
	}
}

// General Definitions

// Initialize list
void initializeListMol(ListMol3D * p_list)
{
	* p_list = NULL;
}

// Initialize Recent list
void initializeListMolRecent(ListMolRecent3D * p_list)
{
	* p_list = NULL;
}

// Is the list empty?
bool isListMol3DEmpty(const ListMol3D * p_list)
{
	if(*p_list == NULL) return true;
	return false;
}

// Is the list empty?
bool isListMol3DRecentEmpty(const ListMolRecent3D * p_list)
{
	if(*p_list == NULL) return true;
	return false;
}

// Visit each node and execute the function pointed to by p_fun
// p_fun must be a void function that takes an item as input
void traverse(const ListMol3D * p_list, void (* p_fun)(ItemMol3D item))
{
	NodeMol3D * p_node = *p_list;
	
	while(p_node != NULL)
	{
		(*p_fun)(p_node->item); // Apply function
		p_node = p_node->next;  // Advance to next item
	}
}

// Visit each node and execute the function pointed to by p_fun
// p_fun must be a void function that takes an item as input
void traverseRecent(const ListMolRecent3D * p_list, void (* p_fun)(ItemMolRecent3D item))
{
	NodeMolRecent3D * p_node = *p_list;
	
	while(p_node != NULL)
	{
		(*p_fun)(p_node->item); // Apply function
		p_node = p_node->next;  // Advance to next item
	}
}

// De-allocate memory of all nodes in the list
void emptyListMol(ListMol3D * p_list)
{
	NodeMol3D * p_save;
	
	while(*p_list != NULL)
	{
		p_save = (*p_list)->next;	// Save address of next node
		free(*p_list);				// Free memory of current node
		*p_list = p_save;			// Advance to next node
	}
}

// De-allocate memory of all nodes in the list
void emptyListMol3DRecent(ListMolRecent3D * p_list)
{
	NodeMolRecent3D * p_save;
	
	while(*p_list != NULL)
	{
		p_save = (*p_list)->next;	// Save address of next node
		free(*p_list);				// Free memory of current node
		*p_list = p_save;			// Advance to next node
	}
}

// Create new node to hold item and add it to the start of the list
bool addItem(ItemMol3D item, ListMol3D * p_list)
{
	NodeMol3D * p_new;
	
	p_new = malloc(sizeof(NodeMol3D));
	if (p_new == NULL)
		return false;	// Quit on failure of malloc
		
	copyToNode(item, p_new);
	p_new->next = * p_list;
	* p_list = p_new; // Point start of list to new node
	return true;
}

// Create new node to hold item and add it to the start of the list
bool addItemRecent(ItemMolRecent3D item, ListMolRecent3D * p_list)
{
	NodeMolRecent3D * p_new;
	
	p_new = malloc(sizeof(NodeMolRecent3D));
	if (p_new == NULL)
		return false;	// Quit on failure of malloc
		
	copyToNodeRecent(item, p_new);
	p_new->next = * p_list;
	* p_list = p_new; // Point start of list to new node
	return true;
}

// Remove node from list. In order to do this efficiently, need to pass Nodes
// and not the start of the list.
void removeItem(NodeMol3D * prevNode, NodeMol3D * curNode)
{
	if (curNode != NULL)
	{
		if (prevNode != NULL)
		{
			// Update pointer of previous node to point to following node
			prevNode->next = curNode->next;		
		}
		free(curNode); // Free memory of node
	}
}

// Remove node from list. In order to do this efficiently, need to pass Nodes
// and not the start of the list.
void removeItemRecent(NodeMolRecent3D * prevNode, NodeMolRecent3D * curNode)
{
	if (curNode != NULL)
	{
		if (prevNode != NULL)
		{
			// Update pointer of previous node to point to following node
			prevNode->next = curNode->next;		
		}
		free(curNode); // Free memory of node
	}
}

// Copy an item to a node
static void copyToNode(ItemMol3D item, NodeMol3D * p_node)
{
	p_node->item = item; // Structure copy
}

// Copy an item to a node
static void copyToNodeRecent(ItemMolRecent3D item, NodeMolRecent3D * p_node)
{
	p_node->item = item; // Structure copy
}
