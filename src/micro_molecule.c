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
 * Last revised for AcCoRD v1.4 (2018-08-06)
 *
 * Revision history:
 *
 * Revision v1.4 (2018-08-06)
 * - added a priori monte carlo (APMC) absorption algorithm as a new surface
 * reaction type. Includes settings for how to define the a priori absorption
 * probability calculation and whether/how to apply a threshold to turn it off
 * - corrected missing assignment for first order reactions that don't have any
 * product molecules to place. Led to memory error when first molecule in a list
 * has to be removed
 *
 * Revision v1.2 (2018-05-30)
 * - fixed implementation of replication reactions, where a first order reactant produces
 * at least one copy of itself. If such a reactant reacted again within the same
 * microscopic time step, then the new molecule(s) previously went missing.
 *
 * Revision v1.1 (2016-12-24)
 * - simplified detection of whether molecules flow or diffuse in each region
 * - added uniform flow to the diffusion algorithm
 * - modified meso-to-micro hybrid transition algorithm when a molecule is placed
 * in the microscopic regime. Now, the trajectory of the molecule will be tracked
 * to make sure that it can reach its intended destination. Reflections are added
 * as necessary. Molecule is assumed to start from the middle of the mesoscopic
 * subvolume.
 *
 * Revision v1.0 (2016-10-31)
 * - added specifying diffusion coefficient that applies to specific surface
 * interaction reactions.
 * - moved mesoscopic structure fields from subvolume struct to meso subvolume struct
 *
 * Revision v0.7 (public beta, 2016-07-09)
 * - set microscopic partial time step to 0 when creating new molecule from meso
 * diffusion.
 * - corrected what region is indexed when new microscopic molecule diffuses away
 * from hybrid interface
 * - use molecule's partial time step when determining whether it entered and exited
 * mesoscopic regime within a time step
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - added check for molecules entering mesoscopic regime "during" a microscopic time step,
 * i.e., when molecule is in micro at start and end of diffusion step.
 * - added function for placing molecules in microscopic regime when they come from the
 * mesoscopic region
 * - modified random number generation. Now use PCG via a separate interface file.
 * - added bimolecular chemical reactions in microscopic regime. Transition of reactants
 * to reaction site is tested if reaction site is not in same region. Transition of
 * product molecules from reaction site is also tested if unbinding radius is used.
 * Molecules are compared with all potential reactants in current region and neighboring
 * regions where the reaction is also valid. A molecule cannot participate in more than
 * one bimolecular reaction is a single time step (as a reactant or product). Reaction
 * site and displacement of products depends on relative diffusion coefficients
 * - added check on total number of chemical reactions in a region before checking for
 * surface reactions during diffusion validation, since many surface reaction structure
 * members are not initialized unless there is at least one reaction at the surface
 * - corrected bug where a molecule is not correctly reflected off of a child region
 * surface
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

// Move ALL molecules in the list via diffusion and flow
void diffuseMolecules(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	ListMol3D p_list[NUM_REGIONS][NUM_MOL_TYPES],
	ListMolRecent3D p_listRecent[NUM_REGIONS][NUM_MOL_TYPES],
	const struct region regionArray[],
	struct mesoSubvolume3D mesoSubArray[],
	double sigma[NUM_REGIONS][NUM_MOL_TYPES],
	const struct chem_rxn_struct chem_rxn[],
	const double HYBRID_DIST_MAX,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES])
{
	NodeMol3D * curNode, * prevNode, * nextNode;
	NodeMolRecent3D * curNodeR;
	
	double oldPoint[3];
	double newPoint[3];
	unsigned short curType, newType;
	short curRegion, newRegion, transRegion;
	uint32_t newSub, curBoundSub;
	bool bInRegion, bPointChange;
	int curMol = 0;
	
	bool bReaction;
	unsigned short curRxn, curProd;

	bool bRemove, bValidDiffusion;
	uint32_t minSub;
	
	// A Priori surface reaction parameters
	bool bApmc; // Is there an A Priori surface reaction that we need to consider?
	bool bApmcCur; // Did current molecule have (or attempt) an A Priori surface reaction?
	unsigned short curGlobalRxn;
	
	// Indicate that every microscopic molecule in a "normal" list
	// needs to be moved.
	// We do this to avoid moving a molecule more than once if it is moved
	// to a different region.
	for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
	{
		for(curType = 0; curType < NUM_MOL_TYPES; curType++)
		{
			if(isListMol3DEmpty(&p_list[curRegion][curType])
				|| (!regionArray[curRegion].bDiffuse[curType]
				&& !regionArray[curRegion].spec.bFlow[curType]))
				continue; // No need to validate an empty list of molecules
						  // or molecules that cannot move
			
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
				|| (!regionArray[curRegion].bDiffuse[curType]
				&& !regionArray[curRegion].spec.bFlow[curType]))
				continue; // No need to validate an empty list of molecules
						  // or molecules that cannot move
			
			curNode = p_list[curRegion][curType];
			prevNode = NULL;
			
			// Check whether we have a possible A Priori surface reaction
			if(regionArray[curRegion].spec.surfaceType == NO_SURFACE
				&& regionArray[curRegion].numApmcRxn[curType] > 0)
				bApmc = true;
			else
			{
				bApmc = false;
				bApmcCur = false;
			}
			
			while(curNode != NULL)
			{
				nextNode = curNode->next;
				
				if(curNode->item.bNeedUpdate)
				{ // Molecule needs to be moved
					curNode->item.bNeedUpdate = false;
				
					oldPoint[0] = curNode->item.x;
					oldPoint[1] = curNode->item.y;
					oldPoint[2] = curNode->item.z;
					
					if (bApmc)
					{// Apply A Priori surface reaction test
						bApmcCur = testApmcRxn(oldPoint, newPoint, curRegion, &newRegion, &curRxn,
							curType, regionArray[curRegion].spec.dt, NUM_REGIONS,
							regionArray, NUM_MOL_TYPES, chem_rxn, DIFF_COEF, &curGlobalRxn);
					}
					
					if(!bApmc || !bApmcCur)
					{ // Diffusion for this molecule can proceed
						while(true)
						{
							// Diffuse molecule
							if(regionArray[curRegion].bDiffuse[curType])
								diffuseOneMolecule(&curNode->item, sigma[curRegion][curType]);
							
							// Move molecule via flow
							if(regionArray[curRegion].spec.bFlow[curType])
								flowTransportOneMolecule(&curNode->item,
									regionArray[curRegion].spec.flowType[curType],
									regionArray[curRegion].flowConstant[curType]);
							
							newPoint[0] = curNode->item.x;
							newPoint[1] = curNode->item.y;
							newPoint[2] = curNode->item.z;
										
							bReaction = false;
							bValidDiffusion = validateMolecule(newPoint, oldPoint, NUM_REGIONS,
								NUM_MOL_TYPES, curRegion, &newRegion, &transRegion, &bPointChange,
								regionArray, curType, &bReaction, &bApmcCur,
								false, regionArray[curRegion].spec.dt, chem_rxn, DIFF_COEF, &curRxn);
							
							if(bApmc && bApmcCur)
							{ // Molecule hit region excluded by A Priori test
								// Need to revert diffusion step and re-attempt
								moveMolecule(&curNode->item, oldPoint[0], oldPoint[1], oldPoint[2]);
							}
							else
								break; // We can break out of this while loop
						}
					}
					
					if(bApmc && bApmcCur)
					{ // Need to fire corresponding surface reaction
						bRemove = true;
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
					} else if(regionArray[newRegion].spec.bMicro)
					{
						// Check for entering meso region within time step, even though
						// we weren't in meso region at end of time step
						if(regionArray[newRegion].bHasMesoNeigh
							&& bEnterMesoIndirect(NUM_REGIONS, NUM_MOL_TYPES,
							regionArray,	 curType, newRegion, &transRegion, oldPoint, newPoint,
							&minSub, HYBRID_DIST_MAX, regionArray[curRegion].spec.dt, DIFF_COEF))
						{
							newSub = regionArray[transRegion].neighID[newRegion][minSub];
							mesoSubArray[newSub].num_mol[curType]++;							
							regionArray[transRegion].bNeedUpdate[newRegion][minSub] = true;
							regionArray[transRegion].numMolFromMicro[newRegion][minSub][curType]++;
							
							bRemove = true; // Remove molecule from curRegion list
						} else if(bValidDiffusion)
						{ // Molecule is still within region and point did not change
							bRemove = false;
							prevNode = curNode;
						} else if(newRegion == curRegion)
						{
							// Molecule is still within region. Just need to update position
							moveMolecule(&curNode->item, newPoint[0], newPoint[1], newPoint[2]);
							bRemove = false;
							prevNode = curNode;
						} else
						{
							// Molecule is in a new microscopic region. Move to appropriate list
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
							bRemove = true;
						}
					} else
					{
						// New region is mesoscopic. Find nearest subvolume to new point
						newSub = regionArray[newRegion].neighID[transRegion]
							[findNearestSub(newRegion, regionArray,
							transRegion, newPoint[0], newPoint[1], newPoint[2])];
						mesoSubArray[newSub].num_mol[curType]++;
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
						
						bRemove = true; // Remove molecule from curRegion list
					}
					
					if(bRemove)
					{
						// Remove molecule from current region molecule list
						removeItem(prevNode,curNode);
						
						if(prevNode == NULL)
						{ // We removed first molecule in list.
						  // nextNode is now the start of the list
						  // (i.e., we must update pointer to list)
						  p_list[curRegion][curType] = nextNode;
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
			
			// Check whether we have a possible A Priori surface reaction
			if(regionArray[curRegion].spec.surfaceType == NO_SURFACE
				&& regionArray[curRegion].numApmcRxn[curType] > 0)
				bApmc = true;
			else
			{
				bApmc = false;
				bApmcCur = false;
			}
			
			while(curNodeR != NULL)
			{				
				oldPoint[0] = curNodeR->item.x;
				oldPoint[1] = curNodeR->item.y;
				oldPoint[2] = curNodeR->item.z;
					
				if (bApmc)
				{// Apply A Priori surface reaction test
					bApmcCur = testApmcRxn(oldPoint, newPoint, curRegion, &newRegion, &curRxn,
						curType, curNodeR->item.dt_partial, NUM_REGIONS,
						regionArray, NUM_MOL_TYPES, chem_rxn, DIFF_COEF, &curGlobalRxn);
				}
				
				
				if(!bApmc || !bApmcCur)
				{ // Diffusion for this molecule can proceed
					while(true)
					{
						
						// Diffuse molecule
						if(regionArray[curRegion].bDiffuse[curType])
							diffuseOneMoleculeRecent(&curNodeR->item, DIFF_COEF[curRegion][curType]);
							
						// Move molecule via flow
						if(regionArray[curRegion].spec.bFlow[curType])
							flowTransportOneMoleculeRecent(&curNodeR->item,
							regionArray[curRegion].spec.flowType[curType],
							regionArray[curRegion].spec.flowVector[curType]);
						
						newPoint[0] = curNodeR->item.x;
						newPoint[1] = curNodeR->item.y;
						newPoint[2] = curNodeR->item.z;
						
						// Once molecule is validated, we can proceed directly to transferring
						// it to the relevant "normal" list and remove it from this list				
						bReaction = false;
						validateMolecule(newPoint, oldPoint, NUM_REGIONS, NUM_MOL_TYPES, curRegion,
							&newRegion, &transRegion, &bPointChange,
							regionArray, curType, &bReaction, &bApmcCur,
							true, curNodeR->item.dt_partial, chem_rxn, DIFF_COEF, &curRxn);
							
						if(bApmc && bApmcCur)
						{ // Molecule hit region excluded by A Priori test
							// Need to revert diffusion step and re-attempt
							moveMoleculeRecent(&curNodeR->item, oldPoint[0], oldPoint[1], oldPoint[2]);
						}
						else
							break; // We can break out of this while loop
					}
				}
				
				
				if(bApmc && bApmcCur)
				{ // Need to fire corresponding surface reaction
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
						}
					}						
				} else if(regionArray[newRegion].spec.bMicro)
				{ // Region is microscopic.
					
					if(regionArray[newRegion].bHasMesoNeigh
						&& bEnterMesoIndirect(NUM_REGIONS, NUM_MOL_TYPES,
						regionArray,	 curType, newRegion, &transRegion, oldPoint, newPoint,
						&minSub, HYBRID_DIST_MAX, curNodeR->item.dt_partial, DIFF_COEF))
					{
						newSub = regionArray[transRegion].neighID[newRegion][minSub];
						mesoSubArray[newSub].num_mol[curType]++;							
						regionArray[transRegion].bNeedUpdate[newRegion][minSub] = true;
						regionArray[transRegion].numMolFromMicro[newRegion][minSub][curType]++;
						
					} else if(bReaction)
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
							}
						}									
					} else
					{
						if(!addMolecule(&p_list[newRegion][curType],
							newPoint[0], newPoint[1], newPoint[2]))
						{
							fprintf(stderr, "ERROR: Memory allocation to move molecule between recent molecule list of region %u and list of region %u.\n", curRegion, newRegion);
							exit(EXIT_FAILURE);
						}
					}
				} else
				{ // New region is mesoscopic. Find nearest subvolume to new point
					newSub = regionArray[newRegion].neighID[transRegion]
						[findNearestSub(newRegion, regionArray,
						transRegion, newPoint[0], newPoint[1], newPoint[2])];
					mesoSubArray[newSub].num_mol[curType]++;
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
	molecule->x = generateNormal(molecule->x, sigma);
	molecule->y = generateNormal(molecule->y, sigma);
	molecule->z = generateNormal(molecule->z, sigma);
}

// Move one molecule according to a flow vector
void flowTransportOneMolecule(ItemMol3D * molecule,
	const unsigned short flowType,
	double *flowConstant)
{
	switch(flowType)
	{
		case FLOW_UNIFORM:
			molecule->x += flowConstant[0];
			molecule->y += flowConstant[1];
			molecule->z += flowConstant[2];
			break;
	}
}

// Move one molecule by some standard deviation
void diffuseOneMoleculeRecent(ItemMolRecent3D * molecule, double DIFF_COEF)
{
	double sigma = sqrt(2*molecule->dt_partial*DIFF_COEF);
	molecule->x = generateNormal(molecule->x, sigma);
	molecule->y = generateNormal(molecule->y, sigma);
	molecule->z = generateNormal(molecule->z, sigma);
}

// Move one molecule according to a flow vector
void flowTransportOneMoleculeRecent(ItemMolRecent3D * molecule,
	const unsigned short flowType,
	double *flowVector)
{
	switch(flowType)
	{
		case FLOW_UNIFORM:
			molecule->x += flowVector[0]*molecule->dt_partial;
			molecule->y += flowVector[1]*molecule->dt_partial;
			molecule->z += flowVector[2]*molecule->dt_partial;
			break;
	}
}

// Did molecule enter mesoscopic region while diffusing?
bool bEnterMesoIndirect(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const struct region regionArray[],
	const short curType,
	const short curRegion,
	short *mesoRegion,
	const double oldPoint[3],
	const double newPoint[3],
	uint32_t * newSub,
	const double HYBRID_DIST_MAX,
	const double tLeft,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES])
{
	double minDistSq, curDistSq, initDist, finalDist, minDist, curProb;
	short curNeigh, newRegion;
	unsigned short curFace, minFace;
	uint32_t curSub, minSub, numSub;
	uint32_t nearestSub;
	
	if(!(HYBRID_DIST_MAX > 0))
		return false;
	
	for(curNeigh = 0; curNeigh < regionArray[curRegion].numRegionNeigh; curNeigh++)
	{
		newRegion = regionArray[curRegion].regionNeighID[curNeigh];
		if(regionArray[newRegion].spec.bMicro)
			continue; // Only need to check mesoscopic neighbors
				
		// Find distance from this neighbor to final point
		minSub = findNearestSub(newRegion, regionArray, curRegion,
			newPoint[0], newPoint[1], newPoint[2]);
		finalDist = distanceToBoundary(newPoint, RECTANGULAR_BOX,
			regionArray[newRegion].boundSubCoor[curRegion][minSub]);
		
		if(finalDist > HYBRID_DIST_MAX)
			continue;
										
		// Find distance from this neighbor to starting point
		minSub = findNearestSub(newRegion, regionArray, curRegion,
			oldPoint[0], oldPoint[1], oldPoint[2]);
		initDist = distanceToBoundary(oldPoint, RECTANGULAR_BOX,
			regionArray[newRegion].boundSubCoor[curRegion][minSub]);
		
		if(initDist > HYBRID_DIST_MAX)
			continue;
		
		curProb = exp(-finalDist * initDist
			/DIFF_COEF[curRegion][curType]/tLeft);
		if(generateUniform()<curProb)
		{ // Molecule enters Meso region
			*mesoRegion = newRegion;
			*newSub = minSub;
			return true;
		}
	}
	
	return false;
}

// Place a molecule entering microscopic region from a mesoscopic subvolume
bool placeInMicroFromMeso(const unsigned short curRegion,
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const unsigned short destRegion,
	uint32_t * newSub,
	const struct region regionArray[],
	const uint32_t curBoundSub,
	const bool bSmallSub,
	const unsigned short curMolType,
	ListMolRecent3D pRecentList[NUM_REGIONS][NUM_MOL_TYPES],
	const struct chem_rxn_struct chem_rxn[],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES])
{
	double randCoor[3];
	double newPoint[3];
	double length; // Length of line to the new point
	double lengthToInterface;
	double L[3]; // Line vector between subvolume center and new point
	double intersectPoint[3];
	short newRegion, transRegion;
	bool bPointChange, bReaction, bApmcRevert;
	unsigned short curRxn;
	
	unsigned short faceDir;
	double curRand;
	
	// Determine direction that new molecule should be placed
	if(regionArray[curRegion].boundSubNumFace[destRegion][curBoundSub] > 0)
	{
		// More than one face of this subvolume faces this micro region
		faceDir = (unsigned short) floor(generateUniform()*
			regionArray[curRegion].boundSubNumFace[destRegion][curBoundSub]);
	} else
		faceDir = 0;
	
	curRand = generateUniform();
	randCoor[0] = sqrt(2*DIFF_COEF[curRegion][curMolType]*regionArray[destRegion].spec.dt)*
		(0.729614*curRand - 0.70252*curRand*curRand)/
		(1 - 1.47494*curRand + 0.484371*curRand*curRand);
	if(bSmallSub)
	{ // Assume size of subvolume is on order of microscopic diffusion step
		randCoor[1] = regionArray[curRegion].actualSubSize
			*(generateTriangular()-0.5);
		randCoor[2] = regionArray[curRegion].actualSubSize
			*(generateTriangular()-0.5);
	} else
	{ // Assume size of subvolume is much larger than microscopic diffusion step
		randCoor[1] = regionArray[curRegion].actualSubSize
			*generateUniform();
		randCoor[2] = regionArray[curRegion].actualSubSize
			*generateUniform();
	}
	
	// Determine coordinates of new molecule in microscopic region
	switch(regionArray[curRegion].boundVirtualNeighDir[destRegion][curBoundSub][faceDir])
	{
		case LEFT:
			newPoint[0] =
				regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][0]
				- randCoor[0];
			newPoint[1] = regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][2]
				+ randCoor[1];
			newPoint[2] = regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][4]
				+ randCoor[2];
			break;
		case RIGHT:
			newPoint[0] =
				regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][1]
				+ randCoor[0];
			newPoint[1] = regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][2]
				+ randCoor[1];
			newPoint[2] = regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][4]
				+ randCoor[2];
			break;
		case DOWN:
			newPoint[1] =
				regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][2]
				- randCoor[0];
			newPoint[0] = regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][0]
				+ randCoor[1];
			newPoint[2] = regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][4]
				+ randCoor[2];
			break;
		case UP:
			newPoint[1] =
				regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][3]
				+ randCoor[0];
			newPoint[0] = regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][0]
				+ randCoor[1];
			newPoint[2] = regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][4]
				+ randCoor[2];
			break;
		case IN:
			newPoint[2] =
				regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][4]
				- randCoor[0];
			newPoint[0] = regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][0]
				+ randCoor[1];
			newPoint[1] = regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][2]
				+ randCoor[2];
			break;
		case OUT:
			newPoint[2] =
				regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][5]
				+ randCoor[0];
			newPoint[0] = regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][0]
				+ randCoor[1];
			newPoint[1] = regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][2]
				+ randCoor[2];
			break;
		default:
			fprintf(stderr,"ERROR: Invalid direction defined for subvolume in region %u having virtual subvolume neighbor in region %u.\n",
			curRegion, destRegion);
			exit(EXIT_FAILURE);
	}
	
	defineLine(regionArray[curRegion].boundSubCenterCoor[destRegion][curBoundSub],
		newPoint, L, &length);
	
	// Check trajectory of molecule from center of subvolume to new point
	if(bLineHitInfinitePlane(
		regionArray[curRegion].boundSubCenterCoor[destRegion][curBoundSub],
		L, length, RECTANGULAR_BOX,
		regionArray[curRegion].boundSubCoor[destRegion][curBoundSub],
		regionArray[curRegion].boundVirtualNeighDir[destRegion][curBoundSub][faceDir],
		true, &lengthToInterface, intersectPoint, false))
	{
		// Check whether the intersection point is on the subvolume face
		// If not, then correct
		if(intersectPoint[0] <
			regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][0])
			intersectPoint[0] =
				regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][0];
		else if(intersectPoint[0] >
			regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][1])
			intersectPoint[0] =
				regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][1];
				
		if(intersectPoint[1] <
			regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][2])
			intersectPoint[1] =
				regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][2];
		else if(intersectPoint[1] >
			regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][3])
			intersectPoint[1] =
				regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][3];
				
		if(intersectPoint[2] <
			regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][4])
			intersectPoint[2] =
				regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][4];
		else if(intersectPoint[2] >
			regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][5])
			intersectPoint[2] =
				regionArray[curRegion].boundSubCoor[destRegion][curBoundSub][5];
		
		// "Push" slightly into microscopic region
		pushPoint(intersectPoint, intersectPoint, 0.01*(length-lengthToInterface), L);
		
		// Follow trajectory to the new molecule point
		validateMolecule(newPoint, intersectPoint, NUM_REGIONS, NUM_MOL_TYPES,
			destRegion, &newRegion, &transRegion, &bPointChange,
			regionArray, curMolType, &bReaction, &bApmcRevert, true, 0, chem_rxn,
			DIFF_COEF, &curRxn);
		
		if(regionArray[newRegion].spec.bMicro)
		{
			if(!addMoleculeRecent(&pRecentList[newRegion][curMolType],
				newPoint[0], newPoint[1], newPoint[2], 0))
			{ // Creation of molecule failed
				fprintf(stderr,"ERROR: Memory allocation to create molecule of type %u transitioning from region %u to region %u.\n",
					curMolType, curRegion, destRegion);
				exit(EXIT_FAILURE);
			}
			return true;
		} else
		{ // Molecule ended up back in mesoscopic regime
			*newSub = regionArray[newRegion].neighID[transRegion]
				[findNearestSub(newRegion, regionArray,
				transRegion, newPoint[0], newPoint[1], newPoint[2])];
			return false;
		}
	}
	
	
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
			curRand = generateUniform();
			if(curRand < regionArray[curRegion].surfRxnOutProb[curMolType])
			{				
				// Generate products (if necessary)
				// TODO: May need to consider associating a separation distance with
				// reactions that have more than one product.
				if(regionArray[curRegion].numRxnProducts[curRxn] > 0)
				{
					rxnFirstOrderProductPlacement(curNode, curNodeRecent,
						curRxn, NUM_REGIONS, NUM_MOL_TYPES, curRegion,
						p_list, pRecentList, regionArray, curMolType, DIFF_COEF, false, NULL);
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
		if(regionArray[curRegion].numFirstRxnWithReactant[curMolType] > 1)
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
		curRand = generateUniform();
		
		for(i = 0; i < regionArray[curRegion].numFirstRxnWithReactant[curMolType]; i++)
		{			
			if(curRand < regionArray[curRegion].uniCumProb[curMolType][i])
			{
				// Reaction i took place
				curRxn = regionArray[curRegion].firstRxnWithReactantID[curMolType][i];
				
				// Generate products (if necessary)
				// TODO: May need to consider associating a separation distance with
				// reactions that have more than one product.
				if(regionArray[curRegion].numRxnProducts[curRxn] > 0)
				{					
					rxnFirstOrderProductPlacement(curNode, curNodeRecent,
						curRxn, NUM_REGIONS, NUM_MOL_TYPES, curRegion,
						p_list, pRecentList, regionArray, curMolType, DIFF_COEF, false, NULL);
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
	bool bProductIsReactant; // bool to catch cases where first molecule
		// in list reacts but its product should go at the start of the list
	double uniCumProb[regionArray[curRegion].numFirstRxnWithReactant[curMolType]];
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
				regionArray, NUM_MOL_TYPES);
			curRand = generateUniform();
			if(curRand < curProb)
			{				
				// Generate products (if necessary)
				// TODO: May need to consider associating a separation distance with
				// reactions that have more than one product.
				if(regionArray[curRegion].numRxnProducts[curRxn] > 0)
				{
					rxnFirstOrderProductPlacement(curNodeOld, curNode,
						curRxn, NUM_REGIONS, NUM_MOL_TYPES, curRegion,
						p_list, pRecentList, regionArray, curMolType, DIFF_COEF, true, &bProductIsReactant);
				} else
					bProductIsReactant = false;
				
				// Remove current molecule from list
				numMolReCheck--; // One less molecule to check for regular reactions
				removeItemRecent(prevNode, curNode);
				bRemove = true;
			}
		
			if(prevNode == NULL && bRemove)
			{	// prevNode does not change, but we removed first molecule in list.
				if(bProductIsReactant)
				{ 	// A new molecule is at the start of the list. prevNode must be update to point to previous node
					prevNode = pRecentList[curRegion][curMolType];
					while(prevNode->next != curNode)
						prevNode = prevNode->next;
					prevNode->next = nextNode;
					bProductIsReactant = false;
				} else
				{ 	// nextNode is now the start of the list
					// (i.e., we must update pointer to list)
					pRecentList[curRegion][curMolType] = nextNode;
				}			
			} else if(!bRemove){
				prevNode = curNode;
			}
			curNode = nextNode;
		}
		
		// Reset molecule list parameters if there are other 1st order reactions
		if(regionArray[curRegion].numFirstRxnWithReactant[curMolType] > 1)
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
		curRand = generateUniform();
		
		// Need to generate probabilities for current molecule based on the partial
		// time step size
		
		for(i = 0; i < regionArray[curRegion].numFirstRxnWithReactant[curMolType]; i++)
		{
			curRxn = regionArray[curRegion].firstRxnWithReactantID[curMolType][i];
			
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
						p_list, pRecentList, regionArray, curMolType, DIFF_COEF, true, &bProductIsReactant);
				} else
					bProductIsReactant = false;
				
				// Remove current molecule from list
				removeItemRecent(prevNode, curNode);
				bRemove = true;
				break; // exit for loop
			}
		}
		
		if(prevNode == NULL && bRemove)
		{	// prevNode does not change, but we removed first molecule in list.
			if(bProductIsReactant)
			{ 	// A new molecule is at the start of the list. prevNode must be update to point to previous node
				prevNode = pRecentList[curRegion][curMolType];
				while(prevNode->next != curNode)
					prevNode = prevNode->next;
				prevNode->next = nextNode;
				bProductIsReactant = false;
			} else
			{ 	// nextNode is now the start of the list
				// (i.e., we must update pointer to list)
				pRecentList[curRegion][curMolType] = nextNode;
			}
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
	const bool bRecent,
	bool * bProductIsReactant)
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
		*bProductIsReactant = false; // bProductIsReactant only exists if bRecent is true
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
				-log((1.-minRxnTimeRV)*generateUniform() + minRxnTimeRV)
				/ regionArray[curRegion].uniSumRate[curMolType];
			timeLeft = curMolRecent->item.dt_partial - curTime;
		} else
		{
			curTime =
				-log((1.-regionArray[curRegion].minRxnTimeRV[curMolType])*generateUniform()
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
		
		if(bRecent && curProdID == curMolType)
			*bProductIsReactant = true;
		
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
					dist = 2*fabs(generateNormal(0,
						sqrt(2*(timeLeft)*DIFF_COEF[curRegion][curProdID])));
					break;
				case PROD_PLACEMENT_STEADY_STATE:
					// Force diffusion of time timeLeft assuming steady state
					// with reverse reaction
					curRand = generateUniform();
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

// Check all second order reactions for microscopic regions
void rxnSecondOrder(const unsigned short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	ListMol3D p_list[NUM_REGIONS][NUM_MOL_TYPES],
	const struct region regionArray[],
	struct mesoSubvolume3D mesoSubArray[],
	const struct chem_rxn_struct chem_rxn[],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES])
{
	short curRegion, neighRegion, rxnRegion, transRegion, destRegion;
	bool bSameRegion;		// Are current and neighboring molecules in same region?
	short curMolType, secondMolType;
	short curRxnSecond, curRxnNeighSecond, curRxnRegion, curRxnNeighRegion, curRxn;
	unsigned short diffRxn; 		// ID of reaction while a reactant or product is diffusing
						// to/from site
	short curProd, curProdRxn;
	NodeMol3D * curNode, * prevNode, * nextNode;
	NodeMol3D * curNeighNode, * prevNeighNode, * nextNeighNode;
	double rxnCoor[3];		// Location where reaction "occurs"
	double sphDef[] = {0., 0., 0., 1.}; 	// Use sphere when we need a random direction
							// Needed when we need to place more than 2 reaction products
	double rxnProdCoor[3];	// Location where product molecule is placed
	double prodDirVec[4]; 	// Unit vector in direction of where product molecule is to be placed
							// Needed when there is an unbinding radius
	double prodDist; 		// Distance of product molecule from
							// "reaction site"
	double relDiff; 			// Ratio of 1st reactant's diffusion coefficient to mutual diffusion
	double sumProdDiff;		// Sum of product diffusion coefficients
	double * relDiffProd;	// Each product's diffusion coefficient scaled by mutual diffusion
							// and unbinding radius
	bool bNeedUnbind;		// Need to apply unbinding radius
	bool bNewRegion;			// We placed a product molecule but destination is not in a valid region
	bool bAddProdMicro;		// Product molecule is going in a microscopic region
	uint32_t newSub, curBoundSub; // Index of subvolume of product that ends up in meso region
	bool bRxn;				// Current molecule reacted
	bool bApmcRxn; 			// Current molecule went to a failed A Priori surface reaction
	short numReactant1Add; 	// # of products same type as reactant 1 in same region
	short numReactant2Add; 	// # of products same type as reactant 2 in same region
	
	// Variables needed if product molecule leaves reactant regions
	double reactantPoint[3];
	bool bDiffRxn; // reactant or product undergoes surface reaction while diffusing
	bool bPointChange; // trajectory of diffusing reactant or product changed
	unsigned short curDiffRxnProd;
	
	// Initialize all potential reactants as capable of reacting
	for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
	{
		if(!regionArray[curRegion].spec.bMicro || regionArray[curRegion].numSecondRxn == 0)
			continue;
		
		for(curMolType = 0; curMolType < NUM_MOL_TYPES; curMolType++)
		{
			// Could this molecule react and are there any in region?
			if(regionArray[curRegion].numSecondRxnWithReactant[curMolType] == 0
				|| isListMol3DEmpty(&p_list[curRegion][curMolType]))
				continue;
			
			// Indicate that molecules in list are potential reactants
			curNode = p_list[curRegion][curMolType];
			while(curNode != NULL)
			{
				curNode->item.bNeedUpdate = true;
				curNode = curNode->next;
			}
		}
	}
	
	// Find and execute second order reactions where valid
	for(curRegion = 0; curRegion < NUM_REGIONS; curRegion++)
	{
		if(!regionArray[curRegion].spec.bMicro || regionArray[curRegion].numSecondRxn == 0)
			continue;
		
		for(curRxnSecond = 0;
			curRxnSecond < regionArray[curRegion].numSecondRxn; curRxnSecond++)
		{
			curRxnRegion = regionArray[curRegion].secondRxn[curRxnSecond];
			curRxn = regionArray[curRegion].globalRxnID[curRxnRegion];
			curMolType = regionArray[curRegion].biReactants[curRxnRegion][0];
			secondMolType = regionArray[curRegion].biReactants[curRxnRegion][1];
			
			if(isListMol3DEmpty(&p_list[curRegion][curMolType]))
				continue;
			
			relDiff = DIFF_COEF[curRegion][curMolType]
				/(DIFF_COEF[curRegion][curMolType] + DIFF_COEF[curRegion][secondMolType]);
			
			// Determine relative unbinding distance if there is more than 1 product and
			// a nonzero unbinding radius
			if(regionArray[curRegion].numRxnProducts[curRxnRegion] > 1
				&& regionArray[curRegion].rUnbind[curRxnRegion] > 0.)
			{
				bNeedUnbind = true;
				relDiffProd = malloc(regionArray[curRegion].numRxnProducts[curRxnRegion] * sizeof(double));
				if(relDiffProd == NULL)
				{
					fprintf(stderr, "ERROR: Memory allocation for relative diffusion of products of reaction %u failed.\n", curRxn);
					exit(EXIT_FAILURE);
				}
				
				sumProdDiff = 0.;
				for(curProdRxn = 0;
				curProdRxn < regionArray[curRegion].numRxnProducts[curRxnRegion];
				curProdRxn++)
				{
					curProd = regionArray[curRegion].productID[curRxnRegion][curProdRxn];
					sumProdDiff += DIFF_COEF[curRegion][curProd];
				}
				
				for(curProdRxn = 0;
				curProdRxn < regionArray[curRegion].numRxnProducts[curRxnRegion];
				curProdRxn++)
				{
					curProd = regionArray[curRegion].productID[curRxnRegion][curProdRxn];
					relDiffProd[curProdRxn] =  regionArray[curRegion].rUnbind[curRxnRegion]
						*DIFF_COEF[curRegion][curProd]/sumProdDiff;
				}
			} else
				bNeedUnbind = false;
			
			// Scan all "first reactant" molecules against "second reactants" in current
			// and remaining regions
			for(neighRegion = curRegion; neighRegion < NUM_REGIONS; neighRegion++)
			{
				bSameRegion = neighRegion == curRegion;
				
				// Is region a neighbor?
				if(!bSameRegion
					&& !regionArray[curRegion].isRegionNeigh[neighRegion])
					continue;
				
				// Is other region microscopic, contain the second molecule type,
				// and allow the reaction to occur?
				if(!regionArray[neighRegion].spec.bMicro
					|| isListMol3DEmpty(&p_list[neighRegion][secondMolType])
					|| regionArray[neighRegion].numChemRxn == 0
					|| !regionArray[neighRegion].bGlobalRxnID[curRxn])
					continue;
				
				// Reaction is possible. Scan molecules in both regions
				curNode = p_list[curRegion][curMolType];
				prevNode = NULL;
				while(curNode != NULL)
				{
					nextNode = curNode->next;
					
					bRxn = false;					
					if(curNode->item.bNeedUpdate)
					{ // Molecule has not yet already reacted
						if(bSameRegion
							&& curMolType == secondMolType)
						{ // We are scanning molecules in the same list
							curNeighNode = nextNode;
							prevNeighNode = curNode;
						}
						else
						{
							curNeighNode = p_list[neighRegion][secondMolType];
							prevNeighNode = NULL;
						}

						while(curNeighNode != NULL)
						{
							nextNeighNode = curNeighNode->next;
							
							if(curNeighNode->item.bNeedUpdate)
							{ // Neighbor molecule has not yet already reacted
								
								// Check distance between molecules with binding radius
								if(!moleculeSeparation(&curNode->item,
									&curNeighNode->item,
									regionArray[curRegion].rBindSq[curRxnRegion]))
								{ // Molecule separation is less than binding radius
								  // Execute reaction
									numReactant1Add = 0;
									numReactant2Add = 0;
									
									// Find ''central point'' of reaction
									rxnCoor[0] = curNode->item.x + relDiff*(curNeighNode->item.x
										- curNode->item.x);
									rxnCoor[1] = curNode->item.y + relDiff*(curNeighNode->item.y
										- curNode->item.y);
									rxnCoor[2] = curNode->item.z + relDiff*(curNeighNode->item.z
										- curNode->item.z);
									
									// Can 1st reactant reach reaction site?
									reactantPoint[0] = curNode->item.x;
									reactantPoint[1] = curNode->item.y;
									reactantPoint[2] = curNode->item.z;
									validateMolecule(rxnCoor, reactantPoint, NUM_REGIONS,
										NUM_MOL_TYPES, curRegion, &rxnRegion, &transRegion,
										&bPointChange, regionArray, curMolType,
										&bDiffRxn, &bApmcRxn, true, 0., chem_rxn, DIFF_COEF,
										&diffRxn);
									
									if(bPointChange)
									{
										// 1st reactant can't reach reaction site
										prevNeighNode = curNeighNode;
										curNeighNode = nextNeighNode;
										continue;
									}
									
									// Can 2nd reactant reach reaction site?
									reactantPoint[0] = curNeighNode->item.x;
									reactantPoint[1] = curNeighNode->item.y;
									reactantPoint[2] = curNeighNode->item.z;
									validateMolecule(rxnCoor, reactantPoint, NUM_REGIONS,
										NUM_MOL_TYPES, neighRegion, &rxnRegion, &transRegion,
										&bPointChange, regionArray, secondMolType,
										&bDiffRxn, &bApmcRxn, true, 0., chem_rxn, DIFF_COEF,
										&diffRxn);
									
									if(bPointChange
										|| !regionArray[rxnRegion].bGlobalRxnID[curRxn])
									{
										// 2nd reactant can't reach reaction site
										// Or reaction takes place in an invalid region
										prevNeighNode = curNeighNode;
										curNeighNode = nextNeighNode;
										continue;
									}
									
									if(regionArray[curRegion].numRxnProducts[curRxnRegion] > 0)
									{										
										for(curProdRxn = 0;
										curProdRxn < regionArray[curRegion].numRxnProducts[curRxnRegion];
										curProdRxn++)
										{
											curProd = regionArray[curRegion].productID[curRxnRegion][curProdRxn];
											
											// Determine location of product molecule
											if(bNeedUnbind)
											{ // We must determine directions of unbinding
												if(regionArray[curRegion].numRxnProducts[curRxnRegion] == 2)
												{
													// Since there are 2 products and 2 reactants, use
													// direction of reactants as direction of products
													if(curProdRxn == 0)
													{
														// Placing first product. Need vector of line towards
														// first reactant
														defineLine2(rxnCoor, curNode->item.x, curNode->item.y,
															curNode->item.z, prodDirVec);
														rxnProdCoor[0] = rxnCoor[0] +
															relDiffProd[curProdRxn]*prodDirVec[0];
														rxnProdCoor[1] = rxnCoor[1] +
															relDiffProd[curProdRxn]*prodDirVec[1];
														rxnProdCoor[2] = rxnCoor[2] +
															relDiffProd[curProdRxn]*prodDirVec[2];
													} else
													{
														// Second product goes in opposite direction from first
														rxnProdCoor[0] = rxnCoor[0] -
															relDiffProd[curProdRxn]*prodDirVec[0];
														rxnProdCoor[1] = rxnCoor[1] -
															relDiffProd[curProdRxn]*prodDirVec[1];
														rxnProdCoor[2] = rxnCoor[2] -
															relDiffProd[curProdRxn]*prodDirVec[2];
													}
												} else
												{
													// There are more than 2 products. Choose a random direction
													uniformPointVolume(prodDirVec, SPHERE, sphDef, true, 0);
													rxnProdCoor[0] = rxnCoor[0] +
														relDiffProd[curProdRxn]*prodDirVec[0];
													rxnProdCoor[1] = rxnCoor[1] +
														relDiffProd[curProdRxn]*prodDirVec[1];
													rxnProdCoor[2] = rxnCoor[2] +
														relDiffProd[curProdRxn]*prodDirVec[2];
												}
												
												// Determine what region product ends up in
												if(bPointInRegionNotChild(rxnRegion, regionArray,
													rxnProdCoor, false))
													destRegion = rxnRegion;
												else
												{
													bDiffRxn = false;
													validateMolecule(rxnProdCoor, rxnCoor, NUM_REGIONS,
														NUM_MOL_TYPES, rxnRegion, &destRegion, &transRegion,
														&bPointChange, regionArray, curProd,
														&bDiffRxn, &bApmcRxn, true, 0., chem_rxn, DIFF_COEF,
														&diffRxn);
													
													if(bDiffRxn &&
													regionArray[destRegion].numRxnProducts[diffRxn] > 0)
													{ // Product molecule went on to react with a surface
														for(curDiffRxnProd = 0;
															curDiffRxnProd < regionArray[destRegion].numRxnProducts[diffRxn];
															curDiffRxnProd++)
														{
															// Add the (curDiffRxnProd)th product to the corresponding molecule list
															if(!addMolecule(
																&p_list[destRegion][regionArray[destRegion].productID[diffRxn][curDiffRxnProd]],
																rxnProdCoor[0], rxnProdCoor[1], rxnProdCoor[2]))
															{ // Creation of molecule failed
																fprintf(stderr, "ERROR: Memory allocation to create molecule of type %u from reaction %u.\n",
																regionArray[destRegion].productID[diffRxn][curDiffRxnProd], diffRxn);
																exit(EXIT_FAILURE);						
															}
															// Indicate that product molecule can't react again
															p_list[destRegion][regionArray[destRegion].productID[diffRxn][curDiffRxnProd]]->item.bNeedUpdate = false;
														}
														continue; // Don't try to place bimolecular product
																	// since it went on to react
													}
												}
											} else
											{ // Place product where reaction occured
												rxnProdCoor[0] = rxnCoor[0];
												rxnProdCoor[1] = rxnCoor[1];
												rxnProdCoor[2] = rxnCoor[2];
												destRegion = rxnRegion;
											}
											
											// Place new product molecule
											if(regionArray[destRegion].spec.bMicro)
											{
												bAddProdMicro = true;
											} else
											{ // Add product to subvolume
												// Confirm meso region is neighbor of curRegion
												if(!regionArray[curRegion].isRegionNeigh[destRegion])
												{ // Destination mesoscopic region is not a neighbor of
													// curRegion. Force product to be placed in curRegion
													bAddProdMicro = true;
													destRegion = curRegion;
												} else
												{
													bAddProdMicro = false;
												
													newSub = regionArray[destRegion].neighID[curRegion]
														[findNearestSub(destRegion, regionArray,
														curRegion, rxnProdCoor[0], rxnProdCoor[1], rxnProdCoor[2])];
													mesoSubArray[newSub].num_mol[curProd]++;
													curBoundSub = 0;
													while(regionArray[destRegion].neighID[curRegion][curBoundSub]
														!= newSub)
													{
														curBoundSub++;
													}
													regionArray[destRegion].bNeedUpdate[curRegion][curBoundSub] = true;
													regionArray[destRegion].numMolFromMicro[curRegion][curBoundSub][curProd]++;
												}
											}
											
											if(curProd == curMolType
												&& curRegion == destRegion)
												numReactant1Add++;
											
											if(curProd == secondMolType
												&& neighRegion == destRegion)
												numReactant2Add++;
											
											if(bAddProdMicro)
											{ // Add current product to corresponding molecule list
												if(!addMolecule(
													&p_list[destRegion][curProd],
													rxnProdCoor[0], rxnProdCoor[1], rxnProdCoor[2]))
												{ // Creation of molecule failed
													fprintf(stderr, "ERROR: Memory allocation to create molecule of type %u from reaction %u.\n",
													curProd, curRxn);
													exit(EXIT_FAILURE);						
												}
												// Indicate that product molecule can't react again
												p_list[destRegion][curProd]->item.bNeedUpdate = false;
											}
										}
									}
									
									// Update start of list or previous node if applicable
									if(prevNode == NULL)
									{
										switch (numReactant1Add)
										{
											case 0:
												p_list[curRegion][curMolType] = nextNode;
												break;
											case 1:
												prevNode = p_list[curRegion][curMolType];
												break;
											case 2:
												prevNode =
													p_list[curRegion][curMolType]->next;
												break;
										}
									}
									if(prevNeighNode == NULL)
									{
										switch (numReactant2Add)
										{
											case 0:
												p_list[neighRegion][secondMolType] = nextNeighNode;
												break;
											case 1:
												prevNeighNode = p_list[neighRegion][secondMolType];
												break;
											case 2:
												prevNeighNode =
													p_list[neighRegion][secondMolType]->next;
												break;
										}
									}
									
									// Remove reactants from their respective lists
									if(curNode->next == curNeighNode)
									{
										// Reactants are adjacent molecules in same molecule list
										// We will be removing both curNode and nextNode
										curNode->next = nextNeighNode;
										nextNode = nextNeighNode;
										prevNeighNode = prevNode;
										if(prevNode == NULL)
										{ // Need to correct start of list again
											p_list[curRegion][curMolType] = nextNode;
										}
									}
									removeItem(prevNode, curNode);
									removeItem(prevNeighNode, curNeighNode);
									
									// Break out of inner while loop
									bRxn = true;
									break;
								}
							}
							
							prevNeighNode = curNeighNode;
							curNeighNode = nextNeighNode;
						}
					}
					if(!bRxn)
						prevNode = curNode;
					curNode = nextNode;
				}
			}
			
			if(bNeedUnbind
				&& relDiffProd != NULL)
			{
				free(relDiffProd);
			}
		}
	}
}

// Compare distance between 2 molecules with threshold
// Return true if distance is greater than threshold
bool moleculeSeparation(ItemMol3D * molecule1, ItemMol3D * molecule2, double threshSq)
{
	double distSq;
	
	distSq = squareDBL(molecule1->x - molecule2->x);
	if(distSq > threshSq)
		return true;
	
	distSq += squareDBL(molecule1->y - molecule2->y);
	if(distSq > threshSq)
		return true;
	
	distSq += squareDBL(molecule1->z - molecule2->z);
	if(distSq > threshSq)
		return true;
	
	return false;
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
	bool * bPointChange,
	const struct region regionArray[],
	unsigned short molType,
	bool * bReaction,
	bool * bApmcRevert,
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
		*bPointChange = false;
		return true;
	} else
	{ // Molecule may have left region's outer boundary or went through a child region
		
		// Define trajectory vector
		defineLine(oldPoint, newPoint, trajLine, &lineLength);
		
		return followMolecule(oldPoint, newPoint, trajLine,
			lineLength, curRegion,
			newRegion, transRegion, bPointChange, NUM_REGIONS, NUM_MOL_TYPES,
			regionArray, molType, bReaction, curRxn, bApmcRevert,
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
	bool * bPointChange,
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const struct region regionArray[],
	unsigned short molType,
	bool * bReaction,
	unsigned short * curRxn,
	bool * bApmcRevert,
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
	* bApmcRevert = false;
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
	* bPointChange = false;
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
					lineLength, startRegion, endRegion, transRegion, bPointChange,
					NUM_REGIONS, NUM_MOL_TYPES, regionArray, molType,
					bReaction, curRxn, bApmcRevert, bRecent, dt,
					chem_rxn, DIFF_COEF, depth+1);
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
				*bPointChange = true;
				return false;
			}
			
			if(regionArray[*endRegion].numChemRxn > 0)
			{				
				switch(regionArray[*endRegion].spec.surfaceType)
				{
					case SURFACE_INNER:
					case SURFACE_OUTER:
						// Need to check for surface absorption
						if(regionArray[*endRegion].bSurfRxnIn[molType])
						{ // Absorption is possible
							*curRxn = regionArray[*endRegion].rxnInID[molType];
							// Check for A Priori surface reaction
							if(chem_rxn[regionArray[*endRegion].globalRxnID[*curRxn]].surfRxnType == RXN_A_PRIORI_ABSORBING)
							{ // Reaction is A Priori. Need to invalidate diffusion
								// step and re-try
								endPoint[0] = startPoint[0];
								endPoint[1] = startPoint[1];
								endPoint[2] = startPoint[2];
								*bPointChange = true;
								*bApmcRevert = true;
								return false;
							}
							if(bRecent)
							{
								// Need to calculate absorption probability
								rxnProb = calculateAbsorptionProb(*endRegion,
									molType, *curRxn,
									dt, NUM_REGIONS, regionArray, NUM_MOL_TYPES);
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
											dt, NUM_REGIONS, regionArray, NUM_MOL_TYPES);
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
											dt, NUM_REGIONS, regionArray, NUM_MOL_TYPES);
									else
										rxnProb = regionArray[*endRegion].surfRxnOutProb[molType];
								} else
									rxnProb = 0.;
								break;
						}
						break;
				}
			
				curRand = generateUniform();
				if(curRand < rxnProb)
				{
					// Reaction curRxn took place
					*bReaction = true;
				}
			} else
				*bReaction = false; // Surface region has no chemical reactions; must reflect
			
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
					*bPointChange = true;
					return false;
				}
					
				if(!regionArray[*endRegion].spec.bMicro)
				{
					// Molecule is entering mesoscopic region. Stop following here.
					endPoint[0] = curIntersectPoint[0];
					endPoint[1] = curIntersectPoint[1];
					endPoint[2] = curIntersectPoint[2];
					* transRegion = startRegion; // Indicate from which micro region we came from
					* bPointChange = true;
					return false;
				}
				lineLength -= lineLength*pushFrac; // Correct line length for having been pushed	
				return followMolecule(curIntersectPoint, endPoint, lineVector,
					lineLength, *endRegion, endRegion, transRegion, bPointChange,
					NUM_REGIONS, NUM_MOL_TYPES, regionArray, molType,
					bReaction, curRxn, bApmcRevert, bRecent, dt,
					chem_rxn, DIFF_COEF, depth+1)
					&& startRegion == *endRegion;
			} else
			{
				endPoint[0] = nearestIntersectPoint[0];
				endPoint[1] = nearestIntersectPoint[1];
				endPoint[2] = nearestIntersectPoint[2];
				* bPointChange = true;
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
		*bPointChange = true;
		return false;		
	}
	// Lock point to actual boundary only if actual reflection occurred
	lockPointToRegion(nearestIntersectPoint, startRegion, reflectRegion, regionArray,
		nearestFace);
	// Follow point past reflection. Define new line unit vector
	defineLine(nearestIntersectPoint, newEndPoint, lineVector, &lineLength);
	
	followMolecule(nearestIntersectPoint, newEndPoint, lineVector,
			lineLength, startRegion, endRegion, transRegion, bPointChange,
			NUM_REGIONS, NUM_MOL_TYPES, regionArray,
			molType, bReaction, curRxn, bApmcRevert, bRecent,
			dt, chem_rxn, DIFF_COEF, depth+1);
	
	endPoint[0] = newEndPoint[0];
	endPoint[1] = newEndPoint[1];
	endPoint[2] = newEndPoint[2];
	*bPointChange = true;
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
