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
 * Last revised for AcCoRD v0.5.1 (2016-05-06)
 *
 * Revision history:
 *
 * Revision v0.5.1 (2016-05-06)
 * - updated first order reaction functions to account for surface reactions that
 * release products from the surface. This is done in a common function (for both old
 * and recent molecules). Placement of products depends on user configuration
 * - fixed membrane transition reactions. They can be reversible and should not have
 * product molecules explicitly defined.
 * - added calls to new functions to determine adsorption/desorption probabilities
 * for recent molecules
 * - corrected how molecules are locked to region boundary when they cross regions
 * - updating reaction probabilities for surface reactions so that user has
 * choices for what calculation to use.
 *
 * Revision v0.5 (2016-04-15)
 * - added surface reactions, including membrane transitions
 * - added switch to record all molecules in a region instead of just those
 * within some specified boundary
 * - corrected distance to end point when a molecule is "pushed" into a neighboring
 * region
 * - added fail check to while loop when a molecule is "pushed" into a 
 * neighboring region. Error will display if we did not end up in specified
 * region or one of its children.
 * - corrected molecule diffusion validation algorithm to reflect off of the correct
 * surface boundary when a reflection is needed
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
	const struct chem_rxn_struct chem_rxn[],
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
					if(validateMolecule(newPoint, oldPoint, NUM_REGIONS, NUM_MOL_TYPES, curRegion,
						&newRegion, &transRegion, regionArray, curType, &bReaction,
						false, regionArray[curRegion].spec.dt, chem_rxn, DIFF_COEF, &curRxn))
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
				validateMolecule(newPoint, oldPoint, NUM_REGIONS, NUM_MOL_TYPES, curRegion,
					&newRegion, &transRegion, regionArray, curType, &bReaction,
					true, curNodeR->item.dt_partial, chem_rxn, DIFF_COEF, &curRxn);
					
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
void rxnFirstOrder(const unsigned short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	unsigned short curRegion,
	ListMol3D p_list[NUM_REGIONS][NUM_MOL_TYPES],
	const struct region regionArray[],
	unsigned short curMolType,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	ListMolRecent3D pRecentList[NUM_REGIONS][NUM_MOL_TYPES])
{
	NodeMol3D * curNode = p_list[curRegion][curMolType];
	NodeMolRecent3D * curNodeRecent = pRecentList[curRegion][curMolType]; // Just used as a placeholder here
	NodeMol3D * prevNode = NULL;
	NodeMol3D * nextNode;
	double curRand;
	double curTime;
	int i; // Reaction loop index
	unsigned short curProd, curProdID;
	unsigned short curRxn;
	unsigned short destRegion;
	double point[3];
	bool bRemove;
	
	if(regionArray[curRegion].bUseRxnOutProb[curMolType])
	{ // This molecule has an independently-calculated desorption-type reaction
		curRxn = regionArray[curRegion].rxnOutID[curMolType];
		while(curNode != NULL)
		{
			nextNode = curNode->next;
			bRemove = false;
			
			// Generate RV to see which/whether first order reaction occurred
			curRand = mt_drand();
			if(curRand < regionArray[curRegion].surfRxnOutProb[curMolType])
			{				
				// Generate products (if necessary)
				// TODO: May need to consider associating a separation distance with
				// reactions that have more than one product.
				if(regionArray[curRegion].numRxnProducts[curRxn] > 0)
				{
					rxnFirstOrderProductPlacement(curNode, curNodeRecent,
						curRxn, NUM_REGIONS, NUM_MOL_TYPES, curRegion,
						p_list, pRecentList, regionArray, curMolType, DIFF_COEF, false);
				}
				
				// Remove current molecule from list
				removeItem(prevNode, curNode);
				bRemove = true;
			}
		
			if(prevNode == NULL && bRemove)
			{	// prevNode does not change, but we removed first molecule in list.
				// nextNode is now the start of the list
				// (i.e., we must update pointer to list)
				p_list[curRegion][curMolType] = nextNode;
			} else if(!bRemove){
				prevNode = curNode;
			}
			curNode = nextNode;
		}
		
		// Reset molecule list parameters if there are other 1st order reactions
		if(regionArray[curRegion].numFirstCurReactant[curMolType] > 1)
		{
			curNode = p_list[curRegion][curMolType];
			prevNode = NULL;
		} else
			return; // No need to scan again
	}
		
	while(curNode != NULL)
	{
		nextNode = curNode->next;
		bRemove = false;
		
		// Generate RV to see which/whether first order reaction occurred
		curRand = mt_drand();
		
		for(i = 0; i < regionArray[curRegion].numFirstCurReactant[curMolType]; i++)
		{			
			if(curRand < regionArray[curRegion].uniCumProb[curMolType][i])
			{
				// Reaction i took place
				curRxn = regionArray[curRegion].firstRxnID[curMolType][i];
				
				// Generate products (if necessary)
				// TODO: May need to consider associating a separation distance with
				// reactions that have more than one product.
				if(regionArray[curRegion].numRxnProducts[curRxn] > 0)
				{					
					rxnFirstOrderProductPlacement(curNode, curNodeRecent,
						curRxn, NUM_REGIONS, NUM_MOL_TYPES, curRegion,
						p_list, pRecentList, regionArray, curMolType, DIFF_COEF, false);
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
			p_list[curRegion][curMolType] = nextNode;
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
void rxnFirstOrderRecent(const unsigned short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	unsigned short curRegion,
	ListMolRecent3D pRecentList[NUM_REGIONS][NUM_MOL_TYPES],
	ListMol3D p_list[NUM_REGIONS][NUM_MOL_TYPES],
	const struct region regionArray[],
	unsigned short curMolType,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	bool bCheckCount,
	uint32_t numMolCheck[NUM_REGIONS][NUM_MOL_TYPES])
{
	NodeMolRecent3D * curNode = pRecentList[curRegion][curMolType];
	NodeMol3D * curNodeOld = p_list[curRegion][curMolType];
	NodeMolRecent3D * prevNode = NULL;
	NodeMolRecent3D * nextNode;
	double curRand, curProb;
	double curTime;
	int i; // Reaction loop index
	unsigned short curProd;
	unsigned short curRxn;
	bool bRemove;
	double uniCumProb[regionArray[curRegion].numFirstCurReactant[curMolType]];
	double minRxnTimeRV;
	uint32_t numMolOrig = numMolCheck[curRegion][curMolType];
	// Number of molecules to check in second while loop
	uint32_t numMolReCheck = numMolCheck[curRegion][curMolType];
	uint32_t curMol = 0;
	unsigned short prodID;
	unsigned short destRegion;
	double point[3];
	
	if(regionArray[curRegion].bUseRxnOutProb[curMolType])
	{ // This molecule has an independently-calculated desorption-type reaction
		curRxn = regionArray[curRegion].rxnOutID[curMolType];
		while(curNode != NULL)
		{			
			if(bCheckCount && ++curMol > numMolOrig)
			{ // We only need to check the first numMolOrig molecules in this list
			  // and we have already done so.
				break;
			}
			
			nextNode = curNode->next;
			bRemove = false;
			
			// Generate RV and reaction probability to see whether
			// desorption reaction occurred
			calculateDesorptionProb(&curProb, curRegion, curMolType,
				regionArray[curRegion].rxnOutID[curMolType],
				curNode->item.dt_partial, NUM_REGIONS,
				regionArray, NUM_MOL_TYPES, DIFF_COEF);
			curRand = mt_drand();
			if(curRand < curProb)
			{				
				// Generate products (if necessary)
				// TODO: May need to consider associating a separation distance with
				// reactions that have more than one product.
				if(regionArray[curRegion].numRxnProducts[curRxn] > 0)
				{
					rxnFirstOrderProductPlacement(curNodeOld, curNode,
						curRxn, NUM_REGIONS, NUM_MOL_TYPES, curRegion,
						p_list, pRecentList, regionArray, curMolType, DIFF_COEF, true);
				}
				
				// Remove current molecule from list
				numMolReCheck--; // One less molecule to check for regular reactions
				removeItemRecent(prevNode, curNode);
				bRemove = true;
			}
		
			if(prevNode == NULL && bRemove)
			{	// prevNode does not change, but we removed first molecule in list.
				// nextNode is now the start of the list
				// (i.e., we must update pointer to list)
				pRecentList[curRegion][curMolType] = nextNode;
			} else if(!bRemove){
				prevNode = curNode;
			}
			curNode = nextNode;
		}
		
		// Reset molecule list parameters if there are other 1st order reactions
		if(regionArray[curRegion].numFirstCurReactant[curMolType] > 1)
		{
			curNode = pRecentList[curRegion][curMolType];
			prevNode = NULL;
		} else
			return; // No need to scan again
	}
	
	curMol = 0;
	while(curNode != NULL)
	{
		if(bCheckCount && ++curMol > numMolReCheck)
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
		
		for(i = 0; i < regionArray[curRegion].numFirstCurReactant[curMolType]; i++)
		{
			curRxn = regionArray[curRegion].firstRxnID[curMolType][i];
			
			if(i == 0)
			{
				uniCumProb[i] = regionArray[curRegion].uniRelativeRate[curMolType][i]
					*(1.-exp(-curNode->item.dt_partial*regionArray[curRegion].uniSumRate[curMolType]));
			} else
			{
				uniCumProb[i] = uniCumProb[i-1]
					+ regionArray[curRegion].uniRelativeRate[curMolType][i]
					*(1.-exp(-curNode->item.dt_partial*regionArray[curRegion].uniSumRate[curMolType]));
			}
			
			if(curRand < uniCumProb[i])
			{
				// Generate products (if necessary)
				// TODO: May need to consider associating a separation distance with
				// reactions that have more than one product.
				if(regionArray[curRegion].numRxnProducts[curRxn] > 0)
				{					
					rxnFirstOrderProductPlacement(curNodeOld, curNode,
						curRxn, NUM_REGIONS, NUM_MOL_TYPES, curRegion,
						p_list, pRecentList, regionArray, curMolType, DIFF_COEF, true);
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
			pRecentList[curRegion][curMolType] = nextNode;
		} else if(!bRemove){
			prevNode = curNode;
		}
		curNode = nextNode;
	}
	
	if(bCheckCount)
	{ // Correct the number of molecules to check in next round
		numMolCheck[curRegion][curMolType] -= numMolOrig;
	}
}

// Place products of 1st order reaction
void rxnFirstOrderProductPlacement(const NodeMol3D * curMol,
	const NodeMolRecent3D * curMolRecent,
	const unsigned short curRxn,
	const unsigned short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	unsigned short curRegion,
	ListMol3D p_list[NUM_REGIONS][NUM_MOL_TYPES],
	ListMolRecent3D pRecentList[NUM_REGIONS][NUM_MOL_TYPES],
	const struct region regionArray[],
	unsigned short curMolType,	
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	const bool bRecent)
{
	double curTime, timeLeft; // Time of reaction and time remaining in time step
	double curRand, dist;
	double minRxnTimeRV;
	unsigned short curProd, curProdID;
	unsigned short destRegion;
	double point[3];
	double newPoint[3];
	double lineVector[3];
	double lineLength;
	int reflectFace, reflectDim;
	
	if(bRecent)
	{
		point[0] = curMolRecent->item.x;
		point[1] = curMolRecent->item.y;
		point[2] = curMolRecent->item.z;
	} else
	{
		point[0] = curMol->item.x;
		point[1] = curMol->item.y;
		point[2] = curMol->item.z;
	}
	
	if(regionArray[curRegion].bUseRxnOutProb[curMolType])
	{
		timeLeft = 0;
		if(bRecent)
			curTime = curMolRecent->item.dt_partial;
		else
			curTime = regionArray[curRegion].spec.dt;
	} else
	{
		if(bRecent)
		{
			minRxnTimeRV = exp(-curMolRecent->item.dt_partial*regionArray[curRegion].uniSumRate[curMolType]);
			curTime =
				-log((1.-minRxnTimeRV)*mt_drand() + minRxnTimeRV)
				/ regionArray[curRegion].uniSumRate[curMolType];
			timeLeft = curMolRecent->item.dt_partial - curTime;
		} else
		{
			curTime =
				-log((1.-regionArray[curRegion].minRxnTimeRV[curMolType])*mt_drand()
				+ regionArray[curRegion].minRxnTimeRV[curMolType])
				/ regionArray[curRegion].uniSumRate[curMolType];
			timeLeft = regionArray[curRegion].spec.dt - curTime;
		}
	}
	
	for(curProd = 0; curProd < regionArray[curRegion].numRxnProducts[curRxn];
		curProd++)
	{
		// Find actual ID of current product molecule
		curProdID = regionArray[curRegion].productID[curRxn][curProd];
		
		// Reset location of new product
		newPoint[0] = point[0];
		newPoint[1] = point[1];
		newPoint[2] = point[2];
		
		if(regionArray[curRegion].bReleaseProduct[curRxn][curProd])
		{ // Region is a surface and molecule must be released
			switch(regionArray[curRegion].releaseType[curRxn])
			{
				case PROD_PLACEMENT_LEAVE:
					// Molecule doesn't move. Can diffuse normally for timeLeft
					dist = 0.;
					break;
				case PROD_PLACEMENT_FORCE:
					// Force diffusion of time timeLeft away from surface
					dist = 2*fabs(rd_normal(0,
						sqrt(2*(timeLeft)*DIFF_COEF[curRegion][curProdID])));
					break;
				case PROD_PLACEMENT_STEADY_STATE:
					// Force diffusion of time timeLeft assuming steady state
					// with reverse reaction
					curRand = mt_drand();
					if(regionArray[curRegion].bReversible[curRxn])
					{
						dist = sqrt(2*DIFF_COEF[curRegion][curProdID]*(curTime))*
							(0.729614*curRand - 0.70252*curRand*curRand)/
							(1 - 1.47494*curRand + 0.484371*curRand*curRand);
					} else
					{
						dist = sqrt(2*DIFF_COEF[curRegion][curProdID]*(curTime))*
							(0.571825*curRand - 0.552246*curRand*curRand)/
							(1 - 1.53908*curRand + 0.546424*curRand*curRand);
					}
					break;
			}
			
			if(dist > 0.)
			{ // Need to determine new point in direction away from surface
				switch(regionArray[curRegion].spec.shape)
				{
					case RECTANGLE:
						switch(regionArray[curRegion].plane)
						{
							case PLANE_XY:
								reflectDim = 2; // Push off along z-axis
							case PLANE_XZ:
								reflectDim = 1; // Push off along y-axis
							case PLANE_YZ:
								reflectDim = 0; // Push off along x-axis
						}
						if(regionArray[curRegion].spec.surfaceType
							== SURFACE_INNER)
						{ // Move molecule in positive direction
							newPoint[reflectDim] += dist;
						} else
							newPoint[reflectDim] -= dist;
						break;
					case RECTANGULAR_BOX:
						// We need to know exactly what face we're on
						reflectFace = closestFace(point, RECTANGULAR_BOX,
							regionArray[curRegion].boundary);
							
						// Record dimension of face
						if(reflectFace == LEFT || reflectFace == RIGHT)
							reflectDim = 0;
						else if(reflectFace == DOWN || reflectFace == UP)
							reflectDim = 1;
						else if(reflectFace == IN || reflectFace == OUT)
							reflectDim = 2;
						
						// Push away from face in correct direction
						if(reflectFace == LEFT || reflectFace == DOWN
							|| reflectFace == IN)
						{
							if(regionArray[curRegion].spec.surfaceType
								== SURFACE_INNER)
							{ // Move molecule in positive direction
								newPoint[reflectDim] += dist;
							} else
								newPoint[reflectDim] -= dist;							
						} else
						{
							if(regionArray[curRegion].spec.surfaceType
								== SURFACE_INNER)
							{ // Move molecule in positive direction
								newPoint[reflectDim] -= dist;
							} else
								newPoint[reflectDim] += dist;	
						}
						break;
					case SPHERE:
						// Define line from center of sphere to current point
						if(regionArray[curRegion].spec.surfaceType
							== SURFACE_INNER)
						{ // "Starting point" for vector is on sphere
							defineLine(point, regionArray[curRegion].boundary,
								lineVector, &lineLength);
						} else
						{ // "Starting point" for vector is center of sphere
							defineLine(regionArray[curRegion].boundary, point,
								lineVector, &lineLength);
						}
						
						// Now follow line by specified distance
						pushPoint(point, newPoint, dist, lineVector);
						break;
					default:
						break;
				}
			}	
			destRegion = findDestRegion(newPoint, curRegion, regionArray);
		} else
		{ // Leave product molecule in current region
			destRegion = curRegion;
		}
		if(!addMoleculeRecent(&pRecentList[destRegion][curProdID],
			newPoint[0], newPoint[1], newPoint[2], timeLeft))
		{ // Creation of molecule failed
			fprintf(stderr, "ERROR: Memory allocation to create molecule of type %u from reaction %u.\n",
				regionArray[curRegion].productID[curRxn][curProd], curRxn);
			exit(EXIT_FAILURE);						
		}
	}
}

// If a molecule is the product of a surface reaction and it is supposed
// to be released from the surface, find the destination normal region
unsigned short findDestRegion(const double point[3],
	const unsigned short curRegion,
	const struct region regionArray[])
{
	unsigned short curNeigh, curNeighRegion;
	unsigned short minRegion;
	double curDist;
	double minDist = INFINITY;
	
	if(regionArray[curRegion].spec.type == REGION_NORMAL)
	{ // Something went wrong
		fprintf(stderr, "ERROR: Molecule in region %u (Label: \"%s\") is supposed to be the product of a surface region reaction, but this region is not a surface.\n", curRegion,
			regionArray[curRegion].spec.label);
		fprintf(stderr,"Leaving molecule in current region.\n");
		return curRegion;
	}
	
	// First check whether point is already inside a valid normal region
	for(curNeigh = 0; curNeigh < regionArray[curRegion].numRegionNeigh; curNeigh++)
	{
		curNeighRegion = regionArray[curRegion].regionNeighID[curNeigh];
		
		if(regionArray[curNeighRegion].spec.type == REGION_NORMAL)
		{ // Current neighbor region is not a surface
			if(bPointInRegionNotChild(curNeighRegion, regionArray, point, true))
			{ // Point is inside this region. Place here
				return curNeighRegion;
			}
		}
	}
	
	// Point is not already inside a valid region. Find closest normal region
	// that is neighbor of the surface
	for(curNeigh = 0; curNeigh < regionArray[curRegion].numRegionNeigh; curNeigh++)
	{
		curNeighRegion = regionArray[curRegion].regionNeighID[curNeigh];
		
		if(regionArray[curNeighRegion].spec.type == REGION_NORMAL)
		{ // Current neighbor region is not a surface
			curDist = distanceToBoundary(point, regionArray[curNeighRegion].spec.shape,
				regionArray[curNeighRegion].boundary);
			
			if(curDist < minDist)
			{
				minDist = curDist;
				minRegion = curNeighRegion;
			}
		}
	}
	
	if(isfinite(minDist))
		return minRegion;
	
	// We could not find a destination region
	fprintf(stderr, "ERROR: Molecule in region %u (Label: \"%s\") is supposed to be the product of a surface region reaction, but a destination region could not be found.\n", curRegion,
		regionArray[curRegion].spec.label);
	fprintf(stderr,"Leaving molecule in current region.\n");
	return curRegion;
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
	const unsigned short NUM_MOL_TYPES,
	const short curRegion,
	short * newRegion,
	short * transRegion,
	const struct region regionArray[],
	short molType,
	bool * bReaction,
	bool bRecent,
	double dt,
	const struct chem_rxn_struct chem_rxn[],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	unsigned short * curRxn)
{
	double trajLine[3];
	double lineLength;
	
	// Default for newRegion is curRegion
	*newRegion = curRegion;
	
	short trajRegion = curRegion; // Regions passed through by the molecule's trajectory
	
	if(regionArray[curRegion].numChildren < 1
		&& bPointInRegionNotChild(curRegion, regionArray, newPoint, false))
	{ // This is simplest case. No region boundary interactions
		return true;
	} else
	{ // Molecule may have left region's outer boundary or went through a child region
		
		// Define trajectory vector
		defineLine(oldPoint, newPoint, trajLine, &lineLength);
		
		return followMolecule(oldPoint, newPoint, trajLine,
			lineLength, curRegion,
			newRegion, transRegion, NUM_REGIONS, NUM_MOL_TYPES,
			regionArray, molType, bReaction, curRxn,
			bRecent, dt, chem_rxn, DIFF_COEF, 0);
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
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const struct region regionArray[],
	short molType,
	bool * bReaction,
	unsigned short * curRxn,
	bool bRecent,
	double dt,
	const struct chem_rxn_struct chem_rxn[],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	unsigned int depth)
{
	short curNeigh, curRegion, closestNormal;
	short curFace, nearestFace;
	double minDist, curDist, minNormalDist;
	double curIntersectPoint[3];
	double nearestIntersectPoint[3];
	double newEndPoint[3];
	double newLineVector[3];
	bool bCurIntersect = false;
	
	double curRand;
	int i;
	unsigned short curProd, rxnID;
	* bReaction = false;
	bool bContinue = false;
	bool bReflect = false;
	double rxnProb, kPrime, kminus1Prime;
	complex double c1, c2;
	
	double pushFrac = SUB_ADJ_RESOLUTION;
	
	bool bReflectInside;
	short reflectRegion;
	
	// First check all neighbor regions to see which (if any) are intersected first
	minDist = INFINITY;
	minNormalDist = INFINITY;
	* endRegion = SHRT_MAX;
	closestNormal = SHRT_MAX;
	for(curNeigh = 0; curNeigh < regionArray[startRegion].numRegionNeigh; curNeigh++)
	{
		curRegion = regionArray[startRegion].regionNeighID[curNeigh];
		
		switch(regionArray[startRegion].regionNeighDir[curRegion])
		{
			case PARENT:
			case CHILD:
				// one region is parent of other
				// Need to check against all faces of child
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
					curFace, false, &curDist, curIntersectPoint, true)
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
				minNormalDist = curDist;
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
		
		// Check if we started on region boundary
		if(minDist == 0)
		{ // We started on neighbor region boundary. Check whether we diffused "back" into startRegion
			pushPoint(nearestIntersectPoint, curIntersectPoint, lineLength*pushFrac, lineVector);
			if(bPointInRegionNotChild(startRegion, regionArray, curIntersectPoint, false))
			{
				// We didn't actually "cross" region boundary
				// Keep following molecule within own region
				lineLength -= lineLength*pushFrac; // Correct line length for having been pushed
				return followMolecule(curIntersectPoint, endPoint, lineVector,
					lineLength, startRegion, endRegion, transRegion,
					NUM_REGIONS, NUM_MOL_TYPES, regionArray, molType,
					bReaction, curRxn, bRecent, dt, chem_rxn, DIFF_COEF, depth+1);
			}
		}
		
		// "Lock" location at exact boundary
		if(regionArray[startRegion].regionNeighDir[*endRegion] != PARENT
			&& regionArray[startRegion].regionNeighDir[*endRegion] != CHILD)
		{ // Correct nearestFace to refer to nearest face on destination region
			nearestFace =  regionArray[*endRegion].regionNeighDir[startRegion];
		}
		lockPointToRegion(nearestIntersectPoint, startRegion, *endRegion, regionArray,
			nearestFace);
		
		if(regionArray[*endRegion].spec.type != REGION_NORMAL)
		{ // Closest region is some kind of surface
			// Check whether we are actually in child of surface
			if(regionArray[*endRegion].spec.shape != SPHERE
				&& !bPointInRegionOrChild(*endRegion, regionArray, nearestIntersectPoint,
				endRegion, true))
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
			
			switch(regionArray[*endRegion].spec.surfaceType)
			{
				case SURFACE_INNER:
				case SURFACE_OUTER:
					// Need to check for surface absorption
					if(regionArray[*endRegion].bSurfRxnIn[molType])
					{ // Absorption is possible
						*curRxn = regionArray[*endRegion].rxnInID[molType];
						if(bRecent)
						{
							// Need to calculate absorption probability
							rxnProb = calculateAbsorptionProb(*endRegion,
								molType, *curRxn,
								dt, NUM_REGIONS, regionArray, NUM_MOL_TYPES,
								DIFF_COEF);
						} else
						{ // Use pre-calculated probability
							rxnProb = regionArray[*endRegion].surfRxnInProb[molType];
						}
					} else
						rxnProb = 0.;
					break;
				case SURFACE_MEMBRANE:
					// Need to check relative direction so that correct probability
					// is used
					switch(regionArray[startRegion].regionNeighDir[* endRegion])
					{
						case LEFT:
						case DOWN:
						case IN:
						case PARENT:
							// Use "Inner" membrane transition probability
							if(regionArray[*endRegion].bSurfRxnIn[molType])
							{
								*curRxn = regionArray[*endRegion].rxnInID[molType];
								if(bRecent)
									rxnProb = calculateMembraneProb(*endRegion,
										molType, *curRxn,
										dt, NUM_REGIONS, regionArray, NUM_MOL_TYPES,
										DIFF_COEF);
								else
									rxnProb = regionArray[*endRegion].surfRxnInProb[molType];
							} else
								rxnProb = 0.;
							break;
						default:
							// Use "Outer" membrane transition probability
							if(regionArray[*endRegion].bSurfRxnOut[molType])
							{
								*curRxn = regionArray[*endRegion].rxnOutID[molType];
								if(bRecent)
									rxnProb = calculateMembraneProb(*endRegion,
										molType, *curRxn,
										dt, NUM_REGIONS, regionArray, NUM_MOL_TYPES,
										DIFF_COEF);
								else
									rxnProb = regionArray[*endRegion].surfRxnOutProb[molType];
							} else
								rxnProb = 0.;
							break;
					}
					break;
			}
			
			curRand = mt_drand();
			if(curRand < rxnProb)
			{
				// Reaction curRxn took place
				*bReaction = true;
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
						if(minNormalDist > minDist + regionArray[startRegion].subResolution)
						{ // Normal region at same distance as membrane was not detected
							// Check neighbors of membrane for closest normal region
							// Scan normal neighbors of membrane that are children or parent
							for(curNeigh = 0; curNeigh < regionArray[*endRegion].numRegionNeigh; curNeigh++)
							{
								curRegion = regionArray[*endRegion].regionNeighID[curNeigh];
								
								if(curRegion == startRegion
									|| regionArray[curRegion].spec.type != REGION_NORMAL)
									continue;
								
								switch(regionArray[*endRegion].regionNeighDir[curRegion])
								{
									case PARENT:
									case CHILD:
										// One region is other's parent
										// Need to check against all faces of child
										bCurIntersect = bLineHitRegion(startPoint, lineVector, lineLength,
											*endRegion, curRegion, regionArray,
											&curFace, &curDist, curIntersectPoint);
										break;
									default:
										// Regions are adjacent. Only check adjacent plane
										curFace = regionArray[*endRegion].regionNeighDir[curRegion];
										bCurIntersect = bLineHitInfinitePlane(startPoint, lineVector,
											lineLength, RECTANGULAR_BOX,
											regionArray[*endRegion].boundRegionFaceCoor[curRegion][0],
											curFace, false, &curDist, curIntersectPoint, true)
											&& bPointOnFace(curIntersectPoint, RECTANGULAR_BOX,
											regionArray[*endRegion].boundRegionFaceCoor[curRegion][0],
											curFace);
										break;
								}
								
								if(bCurIntersect && curDist < minNormalDist)
								{ // There is intersection with this face and it is closest normal region so far
									minNormalDist = curDist;
									closestNormal = curRegion;
									nearestFace = curFace;
									nearestIntersectPoint[0] = curIntersectPoint[0];
									nearestIntersectPoint[1] = curIntersectPoint[1];
									nearestIntersectPoint[2] = curIntersectPoint[2];
								}
							}
							if(minNormalDist > minDist + regionArray[startRegion].subResolution)
							{ // Could not find region in membrane neighbor list
								// Region is probably a grandparent or grandchild to membrane
								fprintf(stderr, "ERROR: Could not evaluate transition from region %u (Label: \"%s\") through membrane region %u (Label: \"%s\").\n",
								startRegion, regionArray[startRegion].spec.label,
								*endRegion, regionArray[*endRegion].spec.label);
								exit(EXIT_FAILURE);
							}
						}
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
			reflectRegion = *endRegion;
			*endRegion = startRegion;
			bReflectInside = false;
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
				if(!bPointInRegionOrChild(*endRegion, regionArray, curIntersectPoint, endRegion, false))
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
					lineLength, *endRegion, endRegion, transRegion,
					NUM_REGIONS, NUM_MOL_TYPES, regionArray, molType,
					bReaction, curRxn, bRecent, dt, chem_rxn, DIFF_COEF, depth+1)
					&& startRegion == *endRegion;
			} else
			{
				endPoint[0] = nearestIntersectPoint[0];
				endPoint[1] = nearestIntersectPoint[1];
				endPoint[2] = nearestIntersectPoint[2];
				return false;
			}
		}
	} else if(bPointInRegionNotChild(startRegion, regionArray, endPoint, false))
	{
		// Point is still in current region so diffusion is valid with no further modification
		*endRegion = startRegion;
		return true;
	} else
	{
		bReflectInside = true;
		reflectRegion = startRegion;
	}
	
	// Point needs to be reflected off of its own region boundary	
	if(!reflectPoint(startPoint, lineVector, lineLength, endPoint, newEndPoint,
		nearestIntersectPoint, &nearestFace, regionArray[reflectRegion].spec.shape,
		regionArray[reflectRegion].boundary, bReflectInside))
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
			lineLength, startRegion, endRegion, transRegion,
			NUM_REGIONS, NUM_MOL_TYPES, regionArray,
			molType, bReaction, curRxn, bRecent, dt, chem_rxn, DIFF_COEF, depth+1);
	
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
	bool bRecordPos,
	bool bRecordAll)
{
	NodeMol3D * p_node = *p_list;
	uint64_t curCount = 0ULL;
	
	while(p_node != NULL)
	{
		if (bRecordAll || isMoleculeObserved(&p_node->item, obsType, boundary))
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
	bool bRecordPos,
	bool bRecordAll)
{
	NodeMolRecent3D * p_node = *p_list;
	uint64_t curCount = 0ULL;
	
	while(p_node != NULL)
	{
		if (bRecordAll || isMoleculeObservedRecent(&p_node->item, obsType, boundary))
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
