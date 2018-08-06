/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * chem_rxn.c - structure for storing chemical reaction properties
 *
 * Last revised for AcCoRD v1.4 (2018-08-06)
 *
 * Revision history:
 *
 * Revision v1.4 (2018-08-06)
 * - added a priori monte carlo (APMC) absorption algorithm as a new surface
 * reaction type. Includes settings for how to define the a priori absorption
 * probability calculation and whether/how to apply a threshold to turn it off
 *
 * Revision v1.0 (2016-10-31)
 * - enabled local diffusion coefficients. Chemical reactions involving surface
 * interactions can specify the diffusion coefficient to use in transition
 * probabilities as a reaction parameter (default is the molecule's default
 * diffusion coefficient)
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - preliminary implementation of bimolecular reactions in microscopic regime
 * (based on binding and unbinding radii). Can also model molecular crowding
 * - added new members to region array structure to facilitate microscopic
 * bimolecular reactions
 * - added indicator for whether reactions in the global reaction lists can occur
 * in given region, which is needed to asses bimolecular reactions where the
 * reactants are in different regions
 *
 * Revision v0.5.1 (2016-05-06)
 * - added label, bReversible, and labelCoupled so that reactions can be named
 * and coupled together
 * - added bReleaseProduct for surface reactions to indicate which products are
 * released from the surface
 * - updated reaction probabilities for surface reactions so that user has
 * choices for what calculation to use. Added types to store user choices
 * - adsorption, desorption, and membrane probability calculations are mostly based on
 * S.S. Andrews, "Accurate particle-based simulation of adsorption, desorption
 * and partial transmission" Physical Biology, vol. 6, p.046015, 2009
 * - constrained absorbing and desorbing reactions to one per type of molecule
 * at a given region. In many cases these reactions are now treated separately
 * from other types of 1st order reactions
 * - constrained membrane reactions to one inner and one outer reaction per type of
 * molecule at a given membrane region.
 *
 * Revision v0.5 (2016-04-15)
 * - removed limit on number of molecule types
 * - removed limit on number of products in a reaction
 * - added region label when giving errors about region initialization
 * - added bSurface, surfRxnType members to reaction structure to add specialized reactions
 * - implementations of some surface reactions
 *
 * Revision v0.4.1
 * - improved use and format of error messages
 *
 * Revision v0.4
 * - modified use of unsigned long and unsigned long long to uint32_t and uint64_t
 * - modified propensity updates to do a full re-calculation in order to avoid
 * numerical underflow problems
 * - removed deprecated debug function
 * - added restriction of chemical reactions to specific regions
 * - renamed region structure and array of structures to accommodate naming
 * of new boundary structure
 *
 * Revision v0.3.1
 * - header added
*/

#include "chem_rxn.h" // for "Public" declarations

#include <stdio.h>
#include <stdlib.h> // for exit(), malloc
#include <stdbool.h> // for C++ bool conventions
#include <string.h> // for strlen(), strcmp()


//
// "Private" Declarations
//

//
// Definitions
//

/* Allocate space for the regionArray members needed for chemical
 * reactions and initialize their values based on the chem_rxn
 * array of chem_rxn_struct
*/
void initializeRegionChemRxn(const short NUM_REGIONS,
	struct region regionArray[],
	const unsigned short NUM_MOL_TYPES,
	const unsigned short MAX_RXNS,
	const struct chem_rxn_struct chem_rxn[],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES])
{
	short i; // Current region
	unsigned short j; // Current reaction
	unsigned short curRxn; // Current reaction in chem_rxn array (when we iterate over region reactions)
	unsigned short k; // Current molecule type OR reaction exception
	unsigned short curProd; // Index of current product (if a reaction has multiple
							// products of the same molecule type)
	unsigned short n; // Additional counting index
	uint32_t num_reactants; // Number of reactants in a reaction
	uint32_t numTotalProd; // Total number of products in a reaction
	
	bool bFoundReactant; // Have we found the first reactant yet?
	
	unsigned short numInfRxn; // Number of 1st order reactions for current reactant that have infinite
								// reaction rates
	
	bool (* bRxnInRegion)[NUM_REGIONS]; // Which reactions can occur in which regions?
	unsigned short (* rxnInRegionID)[NUM_REGIONS]; // IDs of reactions that can occur in each region
	
	double kPrime; // Normalized reaction rate for surface reactions
	double kminus1Prime;
	double complex c1, c2;
	bool bDesorptionSS = false; // Is a desorption probability based on steay state?
	
	// Build arrays to indicate where each reaction can take place
	bRxnInRegion = malloc(MAX_RXNS * sizeof(bool [NUM_REGIONS]));
	rxnInRegionID = malloc(MAX_RXNS * sizeof(unsigned short [NUM_REGIONS]));
	if(bRxnInRegion == NULL || rxnInRegionID == NULL)
	{
		fprintf(stderr, "ERROR: Memory allocation for region reaction matrix.\n");
		exit(EXIT_FAILURE);
	}
	for(i = 0; i < NUM_REGIONS; i++)
	{
		regionArray[i].numChemRxn = 0;
		for(j = 0; j < MAX_RXNS; j++)
		{
			if(chem_rxn[j].bEverywhere &&
				((chem_rxn[j].bSurface && regionArray[i].spec.type != REGION_NORMAL)
				|| (!chem_rxn[j].bSurface && regionArray[i].spec.type == REGION_NORMAL)))
				bRxnInRegion[j][i] = true;
			else
				bRxnInRegion[j][i] = false;
			if(regionArray[i].spec.label && strlen(regionArray[i].spec.label) > 0)
			{ // Need to check for reaction exceptions since region has a label
				for(k = 0; k < chem_rxn[j].numRegionExceptions; k++)
				{
					if(chem_rxn[j].regionExceptionLabel[k] &&
						strlen(chem_rxn[j].regionExceptionLabel[k]) > 0 &&
						!strcmp(regionArray[i].spec.label, chem_rxn[j].regionExceptionLabel[k]))
					{ // This region is listed as an exception
						bRxnInRegion[j][i] = !bRxnInRegion[j][i];
					}
				}
			}
			// Count the reaction if it can take place in region
			if(bRxnInRegion[j][i])
			{
				rxnInRegionID[regionArray[i].numChemRxn++][i] = j;
			}
		}
	}
		
	// Allocate memory
	for(i = 0; i < NUM_REGIONS; i++)
	{	// TODO: Separate members that are only needed for meso or micro regions

		regionArray[i].numZerothRxn = 0;
		regionArray[i].numFirstRxn = 0;
		regionArray[i].numSecondRxn = 0;
				
		if(regionArray[i].numChemRxn == 0)
			continue; // No reactions in this region; no need to allocate memory
		
		regionArray[i].globalRxnID =
			malloc(regionArray[i].numChemRxn*sizeof(unsigned short));
		regionArray[i].bGlobalRxnID =
			malloc(MAX_RXNS*sizeof(bool));
		regionArray[i].bReversible =
			malloc(regionArray[i].numChemRxn*sizeof(bool));
		regionArray[i].reverseRxnID =
			malloc(regionArray[i].numChemRxn*sizeof(unsigned short));
		regionArray[i].numMolChange =
			malloc(regionArray[i].numChemRxn * sizeof(uint64_t *));
		regionArray[i].bMolAdd =
			malloc(regionArray[i].numChemRxn * sizeof(bool *));
		regionArray[i].numRxnProducts =
			malloc(regionArray[i].numChemRxn*sizeof(uint32_t));
		regionArray[i].productID =
			malloc(regionArray[i].numChemRxn * sizeof(unsigned short *));
		regionArray[i].bReleaseProduct =
			malloc(regionArray[i].numChemRxn * sizeof(bool *));
		regionArray[i].releaseType =
			malloc(regionArray[i].numChemRxn * sizeof(short));
		regionArray[i].bUpdateProp =
			malloc(regionArray[i].numChemRxn * sizeof(bool *));
		regionArray[i].rxnOrder =
			malloc(regionArray[i].numChemRxn*sizeof(uint32_t));
		regionArray[i].rxnRate =
			malloc(regionArray[i].numChemRxn*sizeof(double));
		regionArray[i].zerothRxn =
			malloc(regionArray[i].numChemRxn*sizeof(unsigned short));
		regionArray[i].firstRxn =
			malloc(regionArray[i].numChemRxn*sizeof(unsigned short));
		regionArray[i].secondRxn =
			malloc(regionArray[i].numChemRxn*sizeof(unsigned short));
		regionArray[i].tZeroth =
			malloc(regionArray[i].numChemRxn*sizeof(double));
		regionArray[i].rxnRateZerothMicro =
			malloc(regionArray[i].numChemRxn*sizeof(double));
		regionArray[i].uniReactant =
			malloc(regionArray[i].numChemRxn*sizeof(unsigned short));
		regionArray[i].numFirstRxnWithReactant =
			malloc(NUM_MOL_TYPES*sizeof(unsigned short));
		regionArray[i].firstRxnWithReactantID =
			malloc(NUM_MOL_TYPES*sizeof(unsigned short *));
		regionArray[i].numSecondRxnWithReactant =
			malloc(NUM_MOL_TYPES*sizeof(unsigned short));
		regionArray[i].uniSumRate =
			malloc(NUM_MOL_TYPES*sizeof(double));
		regionArray[i].uniCumProb =
			malloc(NUM_MOL_TYPES*sizeof(double *));
		regionArray[i].uniRelativeRate =
			malloc(NUM_MOL_TYPES*sizeof(double *));
		regionArray[i].minRxnTimeRV =
			malloc(NUM_MOL_TYPES*sizeof(double));
		regionArray[i].rBind =
			malloc(regionArray[i].numChemRxn*sizeof(double));
		regionArray[i].rBindSq =
			malloc(regionArray[i].numChemRxn*sizeof(double));
		regionArray[i].rUnbind =
			malloc(regionArray[i].numChemRxn*sizeof(double));
		regionArray[i].biReactants =
			malloc(regionArray[i].numChemRxn*sizeof(unsigned short [2]));
		regionArray[i].rxnProbType =
			malloc(regionArray[i].numChemRxn*sizeof(short));
		regionArray[i].rxnDiffCoef =
			malloc(regionArray[i].numChemRxn*sizeof(double));
		regionArray[i].bSurfRxnIn =
			malloc(NUM_MOL_TYPES*sizeof(bool));
		regionArray[i].rxnInID =
			malloc(NUM_MOL_TYPES*sizeof(unsigned short));
		regionArray[i].surfRxnInProb =
			malloc(NUM_MOL_TYPES*sizeof(double));
		regionArray[i].bSurfRxnOut =
			malloc(NUM_MOL_TYPES*sizeof(bool));
		regionArray[i].rxnOutID =
			malloc(NUM_MOL_TYPES*sizeof(unsigned short));
		regionArray[i].surfRxnOutProb =
			malloc(NUM_MOL_TYPES*sizeof(double));
		regionArray[i].bUseRxnOutProb =
			malloc(NUM_MOL_TYPES*sizeof(bool));
		
		if(regionArray[i].globalRxnID == NULL
			|| regionArray[i].bGlobalRxnID == NULL
			|| regionArray[i].bReversible == NULL
			|| regionArray[i].reverseRxnID == NULL
			|| regionArray[i].numMolChange == NULL
			|| regionArray[i].bMolAdd == NULL
			|| regionArray[i].numRxnProducts == NULL
			|| regionArray[i].productID == NULL
			|| regionArray[i].bReleaseProduct == NULL
			|| regionArray[i].releaseType == NULL
			|| regionArray[i].bUpdateProp == NULL
			|| regionArray[i].rxnOrder == NULL
			|| regionArray[i].rxnRate == NULL
			|| regionArray[i].zerothRxn == NULL
			|| regionArray[i].firstRxn == NULL
			|| regionArray[i].secondRxn == NULL
			|| regionArray[i].tZeroth == NULL
			|| regionArray[i].rxnRateZerothMicro == NULL
			|| regionArray[i].uniReactant == NULL
			|| regionArray[i].numFirstRxnWithReactant == NULL
			|| regionArray[i].firstRxnWithReactantID == NULL
			|| regionArray[i].numSecondRxnWithReactant == NULL
			|| regionArray[i].uniSumRate == NULL
			|| regionArray[i].uniCumProb == NULL
			|| regionArray[i].uniRelativeRate == NULL
			|| regionArray[i].minRxnTimeRV == NULL
			|| regionArray[i].rBind == NULL
			|| regionArray[i].rBindSq == NULL
			|| regionArray[i].rUnbind == NULL
			|| regionArray[i].biReactants == NULL
			|| regionArray[i].rxnProbType == NULL
			|| regionArray[i].rxnDiffCoef == NULL
			|| regionArray[i].bSurfRxnIn == NULL
			|| regionArray[i].rxnInID == NULL
			|| regionArray[i].surfRxnInProb == NULL
			|| regionArray[i].bSurfRxnOut == NULL
			|| regionArray[i].rxnOutID == NULL
			|| regionArray[i].surfRxnOutProb == NULL
			|| regionArray[i].bUseRxnOutProb == NULL)
		{
			fprintf(stderr, "ERROR: Memory allocation for chemical reactions in region %u (label: \"%s\").\n", i, regionArray[i].spec.label);
			exit(EXIT_FAILURE);
		}
		
		// Initialize region's global reaction bool array
		for(j = 0; j < MAX_RXNS; j++)
		{
			regionArray[i].bGlobalRxnID[j] = false;
		}
		
		for(j = 0; j < regionArray[i].numChemRxn; j++)
		{
			numTotalProd = 0;
			curRxn = rxnInRegionID[j][i];
			regionArray[i].globalRxnID[j] = curRxn;
			regionArray[i].bGlobalRxnID[curRxn] = true;
			
			// Determine total number of products for this reaction
			// Then we know how much memory to assign to productID
			for(k = 0; k < NUM_MOL_TYPES; k++)
			{
				numTotalProd += chem_rxn[curRxn].products[k];
			}
			
			regionArray[i].numMolChange[j] =
				malloc(NUM_MOL_TYPES * sizeof(uint64_t));
			regionArray[i].bMolAdd[j] =
				malloc(NUM_MOL_TYPES * sizeof(bool));
			regionArray[i].productID[j] =
				malloc(numTotalProd * sizeof(unsigned short));
			regionArray[i].bReleaseProduct[j] =
				malloc(numTotalProd * sizeof(bool));
			regionArray[i].bUpdateProp[j] =
				malloc(NUM_MOL_TYPES * sizeof(bool));
			if(regionArray[i].numMolChange[j] == NULL
				|| regionArray[i].bMolAdd[j] == NULL
				|| regionArray[i].productID[j] == NULL
				|| regionArray[i].bReleaseProduct[j] == NULL
				|| regionArray[i].bUpdateProp[j] == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation for chemical reaction %u in region %u (label: \"%s\").\n", j, i, regionArray[i].spec.label);
				exit(EXIT_FAILURE);
			}
		}		
		
		for(j = 0; j < NUM_MOL_TYPES; j++)
		{
			regionArray[i].firstRxnWithReactantID[j] =
				malloc(regionArray[i].numChemRxn * sizeof(unsigned short));
			regionArray[i].uniCumProb[j] =
				malloc(regionArray[i].numChemRxn * sizeof(double));
			regionArray[i].uniRelativeRate[j] =
				malloc(regionArray[i].numChemRxn * sizeof(double));
			if(regionArray[i].firstRxnWithReactantID[j] == NULL
				|| regionArray[i].uniCumProb[j] == NULL
				|| regionArray[i].uniRelativeRate[j] == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation for chemical reactions in region %u (label: \"%s\").\n", i, regionArray[i].spec.label);
				exit(EXIT_FAILURE);
			}
			
			regionArray[i].bSurfRxnIn[j] = false;
			regionArray[i].rxnInID[j] = USHRT_MAX;
			regionArray[i].surfRxnInProb[j] = 0.;
			regionArray[i].bSurfRxnOut[j] = false;
			regionArray[i].rxnOutID[j] = USHRT_MAX;
			regionArray[i].surfRxnOutProb[j] = 0.;
			regionArray[i].bUseRxnOutProb[j] = false;
		}
	}
	
	// Initialize the allocated elements based on chem_rxn array
	for(i = 0; i < NUM_REGIONS; i++)
	{
		regionArray[i].rBindMax = 0;
		if(regionArray[i].numChemRxn == 0)
			continue; // No reactions in this region; no need to proceed
		
		// Initialize elements that are indexed by reaction
		for(j = 0; j < regionArray[i].numChemRxn; j++)
		{
			num_reactants = 0;
			bFoundReactant = false;
			regionArray[i].numRxnProducts[j] = 0;
			curRxn = rxnInRegionID[j][i]; // Current reaction in chem_rxn array
			
			regionArray[i].rBind[j] = chem_rxn[curRxn].rBind;
			regionArray[i].rBindSq[j] = squareDBL(regionArray[i].rBind[j]);
			regionArray[i].rUnbind[j] = chem_rxn[curRxn].rUnbind;
			if(regionArray[i].rBind[j] > regionArray[i].rBindMax)
				regionArray[i].rBindMax = regionArray[i].rBind[j]; // Found new maximum binding radius
			
			regionArray[i].rxnProbType[j] = chem_rxn[curRxn].rxnProbType;
			regionArray[i].rxnDiffCoef[j] = chem_rxn[curRxn].diffusion;
			regionArray[i].releaseType[j] = chem_rxn[curRxn].releaseType;
			
			if((chem_rxn[curRxn].surfRxnType == RXN_MEMBRANE_IN
				|| chem_rxn[curRxn].surfRxnType == RXN_MEMBRANE_OUT)
				&& regionArray[i].spec.surfaceType != SURFACE_MEMBRANE)
			{
				fprintf(stderr, "ERROR: Chemical reaction %u is a membrane reaction but is defined for region %u (Label: \"%s\"), which is not a membrane region.\n",
					j, i, regionArray[i].spec.label);
				exit(EXIT_FAILURE);
			} else if(regionArray[i].spec.surfaceType == SURFACE_MEMBRANE
				&& !(chem_rxn[curRxn].surfRxnType == RXN_MEMBRANE_IN
				|| chem_rxn[curRxn].surfRxnType == RXN_MEMBRANE_OUT))
			{
				fprintf(stderr, "ERROR: Chemical reaction %u is not a membrane reaction but is defined for membrane region %u (Label: \"%s\").\n",
					j, i, regionArray[i].spec.label);
				exit(EXIT_FAILURE);
			}
			
			for(k = 0; k < NUM_MOL_TYPES; k++)
			{
				num_reactants += chem_rxn[curRxn].reactants[k];
								
				switch(chem_rxn[curRxn].reactants[k])
				{
					case 1:
						regionArray[i].uniReactant[j] = k;
						if(bFoundReactant)
							regionArray[i].biReactants[j][1] = k;
						else
							regionArray[i].biReactants[j][0] = k;
						bFoundReactant = true;
						break;
					case 2:
						regionArray[i].biReactants[j][0] = k;
						regionArray[i].biReactants[j][1] = k;
						break;
				}
				
				// Is there a net gain or loss of molecule k?
				if(chem_rxn[curRxn].reactants[k] < chem_rxn[curRxn].products[k])
				{	// Net gain of molecule type k
					regionArray[i].numMolChange[j][k] =
						chem_rxn[curRxn].products[k] - chem_rxn[curRxn].reactants[k];
					regionArray[i].bMolAdd[j][k] = true;
				} else
				{	// Net loss of molecule type k
					regionArray[i].numMolChange[j][k] =
						chem_rxn[curRxn].reactants[k] - chem_rxn[curRxn].products[k];
					regionArray[i].bMolAdd[j][k] = false;					
				}
				
				// Is molecule k a reactant?
				// If it is, then we must update propensity if number
				// of k molecules changes
				if(chem_rxn[curRxn].reactants[k] > 0)
					regionArray[i].bUpdateProp[j][k] = true;
				else
					regionArray[i].bUpdateProp[j][k] = false;
				
				// Is molecule k a product?
				// If it is, then we must list it in the listing of reaction products
				// for the current reaction
				// We also check for whether product must be released from surface regions
				for(curProd = 0; curProd < chem_rxn[curRxn].products[k]; curProd++)
				{
					regionArray[i].bReleaseProduct[j][regionArray[i].numRxnProducts[j]] =
						chem_rxn[curRxn].bReleaseProduct[k];
					regionArray[i].productID[j][regionArray[i].numRxnProducts[j]++] = k;
				}
			}
			
			regionArray[i].rxnOrder[j] = num_reactants;
			
			switch(num_reactants)
			{
				case 0:
					if(chem_rxn[curRxn].bSurface && chem_rxn[curRxn].surfRxnType != RXN_NORMAL)
					{
						fprintf(stderr, "ERROR: Chemical reaction %u is 0th order and must be defined as a normal surface reaction.\n", j);
						exit(EXIT_FAILURE);
					}
				
					// meso rxnRate must be calculated for one subvolume
					if (regionArray[i].plane == PLANE_3D
						&& regionArray[i].spec.type == REGION_NORMAL)
						regionArray[i].rxnRate[j] = chem_rxn[curRxn].k
							*regionArray[i].actualSubSize * regionArray[i].actualSubSize * regionArray[i].actualSubSize;
					else if (regionArray[i].spec.type == REGION_NORMAL
						|| regionArray[i].spec.type == REGION_SURFACE_3D) // Region is 2D so need area and not volume
						regionArray[i].rxnRate[j] = chem_rxn[curRxn].k
							*regionArray[i].actualSubSize * regionArray[i].actualSubSize;
					else
						regionArray[i].rxnRate[j] = chem_rxn[curRxn].k
							*regionArray[i].actualSubSize; // Region is 1D so need length
					// microscopic 0th order rate depends on total region volume
					regionArray[i].rxnRateZerothMicro[regionArray[i].numZerothRxn] =
						chem_rxn[curRxn].k * regionArray[i].volume;
					regionArray[i].zerothRxn[regionArray[i].numZerothRxn++] = j;
					break;
				case 1:
					regionArray[i].rxnRate[j] = chem_rxn[curRxn].k;
					regionArray[i].firstRxn[regionArray[i].numFirstRxn++] = j;
					break;
				case 2:
					if(chem_rxn[curRxn].surfRxnType != RXN_NORMAL)
					{
						fprintf(stderr, "ERROR: Chemical reaction %u is a 2nd order reaction and so it must be defined as a normal reaction.\n", j);
						exit(EXIT_FAILURE);
					}
					if(regionArray[i].spec.bMicro)
					{ // Microscopic bimolecular reaction
						// TODO: For now just copy reaction rate, even though we rely
						// on binding radius alone. Warn if binding radius is 0.
						regionArray[i].rxnRate[j] = chem_rxn[curRxn].k;
						if(chem_rxn[curRxn].rBind == 0.)
						{
							fprintf(stderr, "WARNING: Chemical reaction %u is 2nd order and can occur in the microscopic regime but has a binding radius of 0.\n", j);
						}
					} else
					{ // Mesoscopic bimolecular reaction
						if (regionArray[i].plane == PLANE_3D
							&& regionArray[i].spec.type == REGION_NORMAL)
							regionArray[i].rxnRate[j] = chem_rxn[curRxn].k
								/regionArray[i].actualSubSize / regionArray[i].actualSubSize / regionArray[i].actualSubSize;
						else if (regionArray[i].spec.type == REGION_NORMAL
							|| regionArray[i].spec.type == REGION_SURFACE_3D) // Region is 2D so need area and not volume
							regionArray[i].rxnRate[j] = chem_rxn[curRxn].k
								/regionArray[i].actualSubSize / regionArray[i].actualSubSize;
						else
							regionArray[i].rxnRate[j] = chem_rxn[curRxn].k
								/regionArray[i].actualSubSize; // Region is 1D so need length
					}
					regionArray[i].secondRxn[regionArray[i].numSecondRxn++] = j;
					break;
				default:
					fprintf(stderr, "ERROR: Chemical reaction %u has too many reactants.\n", j);
					exit(EXIT_FAILURE);
			}
		}
		
		// Identify reversible reactions where needed
		for(j = 0; j < regionArray[i].numChemRxn; j++)
		{
			// Check for reaction being reversible and if so then find reverse reaction
			curRxn = rxnInRegionID[j][i]; // Current reaction in chem_rxn array
			regionArray[i].bReversible[j] = chem_rxn[curRxn].bReversible;
			if(regionArray[i].bReversible[j])
			{
				regionArray[i].reverseRxnID[j] = j;
				for(k = 0; k < regionArray[i].numChemRxn; k++)
				{
					if(!strcmp(chem_rxn[curRxn].labelCoupled, chem_rxn[rxnInRegionID[k][i]].label)
					&& j != k)
					{
						// Reaction rxnInRegionID[k][i] is the reverse
						regionArray[i].reverseRxnID[j] = k;
						break;
					}
				}
				if(regionArray[i].reverseRxnID[j] == j)
				{ // Reverse reaction was not found
					fprintf(stderr, "ERROR: Reaction %u (Label: \"%s\") in region %u (Label: \"%s\") has non-existent reverse reaction \"%s\".\n",
					j, chem_rxn[curRxn].label, i, regionArray[i].spec.label,
					chem_rxn[curRxn].labelCoupled);
					exit(EXIT_FAILURE);
				}
			}
		}
		
		// Initialize elements that are sorted by reactant
		for(j = 0; j < NUM_MOL_TYPES; j++)
		{
			regionArray[i].numFirstRxnWithReactant[j] = 0;
			regionArray[i].numSecondRxnWithReactant[j] = 0;
			regionArray[i].uniSumRate[j] = 0.;
			regionArray[i].uniCumProb[j][0] = 0.;
			numInfRxn = 0;
		
			// Scan reactions to see which have current molecule as a reactant
			for(k = 0; k < regionArray[i].numChemRxn; k++)
			{
				curRxn = rxnInRegionID[k][i]; // Reaction ID in chem_rxn	
				switch(regionArray[i].rxnOrder[k])
				{
					case 1:
						if(chem_rxn[curRxn].reactants[j] > 0)
						{
							
							regionArray[i].firstRxnWithReactantID[j][regionArray[i].numFirstRxnWithReactant[j]] = k;
							
							regionArray[i].numFirstRxnWithReactant[j]++;
							
							switch(chem_rxn[curRxn].surfRxnType)
							{
								case RXN_ABSORBING:
								case RXN_A_PRIORI_ABSORBING:
									if(regionArray[i].bSurfRxnIn[j])
									{ // Molecule already has an absorbing reaction
										fprintf(stderr, "ERROR: Molecule type %u in region %u (Label: \"%s\") has more than one absorbing reaction specified.\n",
											j, i, regionArray[i].spec.label);
										exit(EXIT_FAILURE);
									}
									regionArray[i].bSurfRxnIn[j] = true;
									regionArray[i].rxnInID[j] = k;
									
									if(chem_rxn[curRxn].surfRxnType == RXN_ABSORBING)
									{
										// Determine reaction probability given full time step
										regionArray[i].surfRxnInProb[j] =
											calculateAbsorptionProb(i, j, k,
											regionArray[i].spec.dt, NUM_REGIONS,
											regionArray, NUM_MOL_TYPES);								
									}
									break;
								case RXN_MEMBRANE_IN:
									if(regionArray[i].bSurfRxnIn[j])
									{ // Molecule already has a membrane reaction
										fprintf(stderr, "ERROR: Molecule type %u in region %u (Label: \"%s\") has more than one inner membrane reaction specified.\n",
											j, i, regionArray[i].spec.label);
										exit(EXIT_FAILURE);
									}
									regionArray[i].bSurfRxnIn[j] = true;
									regionArray[i].rxnInID[j] = k;
									
									// Determine reaction probability given full time step
									regionArray[i].surfRxnInProb[j] =
										calculateMembraneProb(i, j, k,
										regionArray[i].spec.dt, NUM_REGIONS,
										regionArray, NUM_MOL_TYPES);
									break;
								case RXN_MEMBRANE_OUT:
									if(regionArray[i].bSurfRxnOut[j])
									{ // Molecule already has a membrane reaction
										fprintf(stderr, "ERROR: Molecule type %u in region %u (Label: \"%s\") has more than one outer membrane reaction specified.\n",
											j, i, regionArray[i].spec.label);
										exit(EXIT_FAILURE);
									}
									regionArray[i].bSurfRxnOut[j] = true;
									regionArray[i].rxnOutID[j] = k;
									
									// Determine reaction probability given full time step
									regionArray[i].surfRxnOutProb[j] =
										calculateMembraneProb(i, j, k,
										regionArray[i].spec.dt, NUM_REGIONS,
										regionArray, NUM_MOL_TYPES);
									break;
								case RXN_DESORBING:
									if(regionArray[i].bSurfRxnOut[j])
									{ // Molecule already has a desorbing reaction
										fprintf(stderr, "ERROR: Molecule type %u in region %u (Label: \"%s\") has more than one desorbing reaction specified.\n",
											j, i, regionArray[i].spec.label);
										exit(EXIT_FAILURE);
									}
									regionArray[i].bSurfRxnOut[j] = true;
									regionArray[i].rxnOutID[j] = k;
									
									// Determine reaction probability given full time step
									regionArray[i].bUseRxnOutProb[j] =
										calculateDesorptionProb(&regionArray[i].surfRxnOutProb[j],
										i, j, k, regionArray[i].spec.dt,
										NUM_REGIONS, regionArray, NUM_MOL_TYPES);
									if(regionArray[i].bUseRxnOutProb[j])
										continue;
									// No break here because most desorbing cases have
									// reaction probability treated as a regular 1st order reaction
								case RXN_RECEPTOR:
								case RXN_NORMAL:
									regionArray[i].uniSumRate[j] += chem_rxn[curRxn].k;
									
									if(chem_rxn[curRxn].k == INFINITY)
										numInfRxn++;
									break;
								default:
									fprintf(stderr, "ERROR: Chemical reaction %u has invalid 1st order reaction type %u.\n", curRxn, chem_rxn[curRxn].surfRxnType);
									exit(EXIT_FAILURE);
							}
							
						}
						break;
					case 2:
						if(regionArray[i].biReactants[k][0] == j)
							regionArray[i].numSecondRxnWithReactant[j]++;
						if(regionArray[i].biReactants[k][1] == j)
							regionArray[i].numSecondRxnWithReactant[j]++;
						break;
				}
			}
			
			// Scan reactions with current molecule as reactant to determine the
			// cumulative reaction probabilities
			for(k = 0; k < regionArray[i].numFirstRxnWithReactant[j]; k++)
			{
				curRxn = rxnInRegionID[regionArray[i].firstRxnWithReactantID[j][k]][i];
				switch(chem_rxn[curRxn].surfRxnType)
				{
					case RXN_RECEPTOR:
					case RXN_NORMAL:
					case RXN_DESORBING:
						if(k > 0)
							regionArray[i].uniCumProb[j][k] = regionArray[i].uniCumProb[j][k-1];
								
						// Relative rate will be used when we need to recalculate the probabilities
						// due to smaller time step size
						if(regionArray[i].rxnRate[regionArray[i].firstRxnWithReactantID[j][k]] == INFINITY)
						{
							regionArray[i].uniRelativeRate[j][k] = 1/numInfRxn;
							regionArray[i].uniCumProb[j][k] += 1/numInfRxn;
						}
						else
						{
							
							if(chem_rxn[curRxn].surfRxnType == RXN_DESORBING
								&& regionArray[i].rxnProbType[regionArray[i].firstRxnWithReactantID[j][k]] == RXN_PROB_STEADY_STATE
								&& regionArray[i].bReversible[regionArray[i].firstRxnWithReactantID[j][k]])
							{ // This reaction is reversible steady state desorption
							 // It was excluding from calculation of uniSumRate
							 // Ignore its reaction probability
								regionArray[i].uniRelativeRate[j][k] = 0.;
							} else
							{
								regionArray[i].uniRelativeRate[j][k] =
								regionArray[i].rxnRate[regionArray[i].firstRxnWithReactantID[j][k]]
								/ regionArray[i].uniSumRate[j];
							}							
								
							// The following cumulative probability calculation is valid for when
							// a molecule existed from the end of the last region time step
							regionArray[i].uniCumProb[j][k] +=
								regionArray[i].uniRelativeRate[j][k]
								* (1 - exp(-regionArray[i].spec.dt*regionArray[i].uniSumRate[j]));
						}
						
						break;
					case RXN_ABSORBING:
					case RXN_A_PRIORI_ABSORBING:
					case RXN_MEMBRANE_IN:
					case RXN_MEMBRANE_OUT:
						// Absorbing reactions aren't included in calculation of cumulative
						// reaction probabilities						
						break;
					default:
						fprintf(stderr, "ERROR: Chemical reaction %u has invalid 1st order reaction type %u.\n", j, chem_rxn[curRxn].surfRxnType);
						exit(EXIT_FAILURE);
				}	
			}
			
			regionArray[i].minRxnTimeRV[j] =
				exp(-regionArray[i].spec.dt*regionArray[i].uniSumRate[j]);
		}
		
	}
	
	// Allocate member for A Priori surface reactions
	// A Priori reactions have to be treated separately because they do not
	// occur in the region(s) for which they are defined
	for(i = 0; i < NUM_REGIONS; i++)
	{		
		if(regionArray[i].spec.surfaceType == NO_SURFACE)
		{
			regionArray[i].numApmcRxn =
				malloc(NUM_MOL_TYPES*sizeof(unsigned short));
			regionArray[i].apmcRxnRegion =
				malloc(NUM_MOL_TYPES*sizeof(unsigned short *));
			regionArray[i].apmcRxnID =
				malloc(NUM_MOL_TYPES*sizeof(unsigned short *));
			regionArray[i].apmcGlobalRxnID =
				malloc(NUM_MOL_TYPES*sizeof(unsigned short *));
			regionArray[i].uniCumProbApmc =
				malloc(NUM_MOL_TYPES*sizeof(double *));
			regionArray[i].apmcAlpha =
				malloc(NUM_MOL_TYPES*sizeof(double));
			
			if(regionArray[i].numApmcRxn == NULL
				|| regionArray[i].apmcRxnRegion == NULL
				|| regionArray[i].apmcRxnID == NULL
				|| regionArray[i].apmcGlobalRxnID == NULL
				|| regionArray[i].uniCumProbApmc == NULL
				|| regionArray[i].apmcAlpha == NULL)
			{
				fprintf(stderr, "ERROR: Memory allocation for a priori chemical reactions in region %u (label: \"%s\").\n", i, regionArray[i].spec.label);
				exit(EXIT_FAILURE);
			}
			
			for(j = 0; j < NUM_MOL_TYPES; j++)
				regionArray[i].numApmcRxn[j] = 0;
				
			// Count number of A Priori reactions that could apply for each molecule type
			// Each surface region with some a priori reaction is counted separately
			for(j = 0; j < NUM_REGIONS; j++)
			{
				if(regionArray[j].spec.surfaceType == NO_SURFACE
					|| regionArray[j].numChemRxn == 0)
					continue; // A Priori reactions will only be associated with surface regions
					
				for(k = 0; k < regionArray[j].numChemRxn; k++)
				{
					curRxn = rxnInRegionID[k][j]; // Current reaction in chem_rxn array
					if(chem_rxn[curRxn].bSurface
						&& chem_rxn[curRxn].surfRxnType == RXN_A_PRIORI_ABSORBING)
					{ // We have an a priori type reaction. Check whether it applies
						if(chem_rxn[curRxn].bRxnThreshold
							&& chem_rxn[curRxn].rxnThresholdType == RXN_THRESHOLD_REGION
							&& !regionArray[j].isRegionNeigh[i])
							continue; // A Priori reaction is controlled by region proximity
										// but these two regions aren't neighbours
										
						// Current A Priori reaction applies. Add to count for corresponding
						// molecule type
						regionArray[i].numApmcRxn[regionArray[j].uniReactant[k]]++;
					}
				}
			}
			
			// Allocate memory for the A Priori reactions
			for(j = 0; j < NUM_MOL_TYPES; j++)
			{
				if(regionArray[i].numApmcRxn[j] > 0)
				{
					regionArray[i].apmcRxnRegion[j] =
						malloc(regionArray[i].numApmcRxn[j]*sizeof(unsigned short));
					regionArray[i].apmcRxnID[j] =
						malloc(regionArray[i].numApmcRxn[j]*sizeof(unsigned short));
					regionArray[i].apmcGlobalRxnID[j] =
						malloc(regionArray[i].numApmcRxn[j]*sizeof(unsigned short));
					regionArray[i].uniCumProbApmc[j] =
						malloc(regionArray[i].numApmcRxn[j]*sizeof(double));
						
					if(regionArray[i].apmcRxnRegion[j] == NULL
						|| regionArray[i].apmcRxnID[j] == NULL
						|| regionArray[i].apmcGlobalRxnID[j] == NULL
						|| regionArray[i].uniCumProbApmc[j] == NULL)
					{
						fprintf(stderr, "ERROR: Memory allocation for a priori chemical reactions in region %u (label: \"%s\").\n", i, regionArray[i].spec.label);
						exit(EXIT_FAILURE);
					}
					
					for(k = 0; k < regionArray[i].numApmcRxn[j]; k++)
					{
						regionArray[i].apmcRxnRegion[j][k] = USHRT_MAX;
					}
				}
			}
			
			// Assign A Priori reactions
			for(j = 0; j < NUM_REGIONS; j++)
			{
				if(regionArray[j].spec.surfaceType == NO_SURFACE
					|| regionArray[j].numChemRxn == 0)
					continue; // A Priori reactions will only be associated with surface regions
					
				for(k = 0; k < regionArray[j].numChemRxn; k++)
				{
					curRxn = rxnInRegionID[k][j]; // Current reaction in chem_rxn array
					if(chem_rxn[curRxn].bSurface
						&& chem_rxn[curRxn].surfRxnType == RXN_A_PRIORI_ABSORBING)
					{ // We have an a priori type reaction. Check whether it applies
						if(chem_rxn[curRxn].bRxnThreshold
							&& chem_rxn[curRxn].rxnThresholdType == RXN_THRESHOLD_REGION
							&& !regionArray[j].isRegionNeigh[i])
							continue; // A Priori reaction is controlled by region proximity
										// but these two regions aren't neighbours
										
						// Current A Priori reaction applies. Assign details to current region
						n = 0;
						while(regionArray[i].apmcRxnRegion[regionArray[j].uniReactant[k]][n] < USHRT_MAX)
							n++;
						regionArray[i].apmcRxnRegion[regionArray[j].uniReactant[k]][n] = j;
						regionArray[i].apmcRxnID[regionArray[j].uniReactant[k]][n] = k;
						regionArray[i].apmcGlobalRxnID[regionArray[j].uniReactant[k]][n] = curRxn;
						regionArray[i].uniCumProbApmc[regionArray[j].uniReactant[k]][n] = 0;
						regionArray[i].apmcAlpha[regionArray[j].uniReactant[k]] =
							(chem_rxn[curRxn].k*regionArray[j].boundary[3] + DIFF_COEF[i][regionArray[j].uniReactant[k]])/DIFF_COEF[i][regionArray[j].uniReactant[k]]/regionArray[j].boundary[3];
					}
				}
			}
		} else
		{
			regionArray[i].numApmcRxn = NULL;
			regionArray[i].apmcRxnRegion = NULL;
			regionArray[i].apmcRxnID = NULL;
			regionArray[i].apmcGlobalRxnID = NULL;
			regionArray[i].uniCumProbApmc = NULL;
			regionArray[i].apmcAlpha = NULL;
		}
	}
	
	// Free memory of temporary parameters
	if(bRxnInRegion != NULL)
		free(bRxnInRegion);
	if(rxnInRegionID != NULL)
		free(rxnInRegionID);
}

// Free memory of region parameters
void deleteRegionChemRxn(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	struct region regionArray[])
{
	short i; // Current region
	unsigned short j; // Current reaction
	
	if(regionArray == NULL)
		return;
	
	for(i = 0; i < NUM_REGIONS; i++)
	{
		
		for(j = 0; j < regionArray[i].numChemRxn; j++)
		{
			if(regionArray[i].numMolChange[j] != NULL)
				free(regionArray[i].numMolChange[j]);
			if(regionArray[i].bMolAdd[j] != NULL)
				free(regionArray[i].bMolAdd[j]);
			if(regionArray[i].productID[j] != NULL)
				free(regionArray[i].productID[j]);
			if(regionArray[i].bReleaseProduct[j] != NULL)
				free(regionArray[i].bReleaseProduct[j]);
			if(regionArray[i].bUpdateProp[j] != NULL)
				free(regionArray[i].bUpdateProp[j]);
		}
		
		for(j = 0; j < NUM_MOL_TYPES; j++)
		{
			if(regionArray[i].numChemRxn > 0)
			{
				if(regionArray[i].firstRxnWithReactantID[j] != NULL)
					free(regionArray[i].firstRxnWithReactantID[j]);
				if(regionArray[i].uniCumProb[j] != NULL)
					free(regionArray[i].uniCumProb[j]);
				if(regionArray[i].uniRelativeRate[j] != NULL)
					free(regionArray[i].uniRelativeRate[j]);
			}
			
			// A Priori surface reactions
			if(regionArray[i].spec.surfaceType == NO_SURFACE
				&& regionArray[i].numApmcRxn[j] > 0)
			{
				if(regionArray[i].apmcRxnRegion != NULL && regionArray[i].apmcRxnRegion[j] != NULL)
					free(regionArray[i].apmcRxnRegion[j]);
				if(regionArray[i].apmcRxnID != NULL && regionArray[i].apmcRxnID[j] != NULL)
					free(regionArray[i].apmcRxnID[j]);
				if(regionArray[i].apmcGlobalRxnID != NULL && regionArray[i].apmcGlobalRxnID[j] != NULL)
					free(regionArray[i].apmcGlobalRxnID[j]);
				if(regionArray[i].uniCumProbApmc != NULL && regionArray[i].uniCumProbApmc[j] != NULL)
					free(regionArray[i].uniCumProbApmc[j]);
			}
		}
		
		if(regionArray[i].numChemRxn > 0)
		{
			if(regionArray[i].globalRxnID != NULL) free(regionArray[i].globalRxnID);
			if(regionArray[i].bGlobalRxnID != NULL) free(regionArray[i].bGlobalRxnID);
			if(regionArray[i].bReversible != NULL) free(regionArray[i].bReversible);
			if(regionArray[i].reverseRxnID != NULL) free(regionArray[i].reverseRxnID);
			if(regionArray[i].numMolChange != NULL) free(regionArray[i].numMolChange);
			if(regionArray[i].bMolAdd != NULL) free(regionArray[i].bMolAdd);
			if(regionArray[i].numRxnProducts != NULL) free(regionArray[i].numRxnProducts);
			if(regionArray[i].productID != NULL) free(regionArray[i].productID);
			if(regionArray[i].bReleaseProduct != NULL) free(regionArray[i].bReleaseProduct);
			if(regionArray[i].releaseType != NULL) free(regionArray[i].releaseType);
			if(regionArray[i].bUpdateProp != NULL) free(regionArray[i].bUpdateProp);
			if(regionArray[i].rxnOrder != NULL) free(regionArray[i].rxnOrder);
			if(regionArray[i].rxnRate != NULL) free(regionArray[i].rxnRate);
			if(regionArray[i].zerothRxn != NULL) free(regionArray[i].zerothRxn);
			if(regionArray[i].firstRxn != NULL) free(regionArray[i].firstRxn);
			if(regionArray[i].secondRxn != NULL) free(regionArray[i].secondRxn);
			if(regionArray[i].tZeroth != NULL) free(regionArray[i].tZeroth);
			if(regionArray[i].rxnRateZerothMicro != NULL) free(regionArray[i].rxnRateZerothMicro);
			if(regionArray[i].uniReactant != NULL) free(regionArray[i].uniReactant);
			if(regionArray[i].numFirstRxnWithReactant != NULL)
				free(regionArray[i].numFirstRxnWithReactant);
			if(regionArray[i].firstRxnWithReactantID != NULL)
				free(regionArray[i].firstRxnWithReactantID);
			if(regionArray[i].numSecondRxnWithReactant != NULL)
				free(regionArray[i].numSecondRxnWithReactant);
			if(regionArray[i].uniSumRate != NULL) free(regionArray[i].uniSumRate);
			if(regionArray[i].uniCumProb != NULL) free(regionArray[i].uniCumProb);
			if(regionArray[i].uniRelativeRate != NULL) free(regionArray[i].uniRelativeRate);
			if(regionArray[i].minRxnTimeRV != NULL) free(regionArray[i].minRxnTimeRV);
			if(regionArray[i].rBind != NULL) free(regionArray[i].rBind);
			if(regionArray[i].rBindSq != NULL) free(regionArray[i].rBindSq);
			if(regionArray[i].rUnbind != NULL) free(regionArray[i].rUnbind);
			if(regionArray[i].biReactants != NULL) free(regionArray[i].biReactants);
			if(regionArray[i].rxnProbType != NULL) free(regionArray[i].rxnProbType);
			if(regionArray[i].rxnDiffCoef != NULL) free(regionArray[i].rxnDiffCoef);
			if(regionArray[i].bSurfRxnIn != NULL) free(regionArray[i].bSurfRxnIn);
			if(regionArray[i].rxnInID != NULL) free(regionArray[i].rxnInID);
			if(regionArray[i].surfRxnInProb != NULL) free(regionArray[i].surfRxnInProb);
			if(regionArray[i].bSurfRxnOut != NULL) free(regionArray[i].bSurfRxnOut);
			if(regionArray[i].rxnOutID != NULL) free(regionArray[i].rxnOutID);
			if(regionArray[i].surfRxnOutProb != NULL) free(regionArray[i].surfRxnOutProb);
			if(regionArray[i].bUseRxnOutProb != NULL) free(regionArray[i].bUseRxnOutProb);
		}
		
		// A Priori surface reactions
		if(regionArray[i].spec.surfaceType == NO_SURFACE)
		{
			if(regionArray[i].numApmcRxn != NULL) free(regionArray[i].numApmcRxn);
			if(regionArray[i].apmcRxnRegion != NULL) free(regionArray[i].apmcRxnRegion);
			if(regionArray[i].apmcRxnID != NULL) free(regionArray[i].apmcRxnID);
			if(regionArray[i].apmcGlobalRxnID != NULL) free(regionArray[i].apmcGlobalRxnID);
			if(regionArray[i].uniCumProbApmc != NULL) free(regionArray[i].uniCumProbApmc);
			if(regionArray[i].apmcAlpha != NULL) free(regionArray[i].apmcAlpha);
		}
	}
}

// Calculate probability of absorption reaction for specified time step
double calculateAbsorptionProb(const short curRegion,
	const unsigned short curMolType,
	const unsigned short curRegionRxn,
	const double dt,
	const short NUM_REGIONS,
	const struct region regionArray[],
	const unsigned short NUM_MOL_TYPES)
{
	double kPrime, kminus1Prime;
	complex double c1, c2;
	double rxnProb;
	
	switch(regionArray[curRegion].rxnProbType[curRegionRxn])
	{
		case RXN_PROB_NORMAL:
			// Use independent 1st order reaction rate probability
			// Not accurate for surface transitions
			rxnProb = (1 - exp(-dt
				*regionArray[curRegion].rxnRate[curRegionRxn]));
			break;
		case RXN_PROB_MIXED:
			// Use well-mixed reaction probability
			// S.S. Andrews Physical Biology 2009 Eq. 1
			rxnProb = regionArray[curRegion].rxnRate[curRegionRxn] *
				sqrt(PI*dt/regionArray[curRegion].rxnDiffCoef[curRegionRxn]);
			break;
		case RXN_PROB_STEADY_STATE:
			// Use steady state reaction probabilities
			kPrime = regionArray[curRegion].rxnRate[curRegionRxn]*
				sqrt(dt/regionArray[curRegion].rxnDiffCoef[curRegionRxn]/2);
			if(regionArray[curRegion].bReversible[curRegionRxn])
			{ // S.S. Andrews Physical Biology 2009 Eq. 37
				kminus1Prime =
					regionArray[curRegion].
					rxnRate[regionArray[curRegion].reverseRxnID[curRegionRxn]]
					*dt;
				c1 = (kPrime - csqrt(C(kPrime*kPrime -2*kminus1Prime,0)))/sqrt(2);
				c2 = (kPrime + csqrt(C(kPrime*kPrime -2*kminus1Prime,0)))/sqrt(2);
				rxnProb = cabs(kPrime*sqrt(2*PI)*
					(c2-c1 - c2*cerfcx(c1) + c1*cerfcx(c2))/(c1*c2*(c2-c1)));
			} else
			{ // S.S. Andrews Physical Biology 2009 Eq. 21
				rxnProb =
					kPrime*sqrt(2*PI) - 3.33321*kPrime*kPrime
					+ 3.35669*kPrime*kPrime*kPrime
					- 1.52092*kPrime*kPrime*kPrime*kPrime;
			}
			break;
		default:
			fprintf(stderr, "ERROR: Chemical reaction %u of region %u (Label: \"%s\") has invalid absorbing reaction probability type %u.\n",
				curRegionRxn, curRegion, regionArray[curRegion].spec.label,
				regionArray[curRegion].rxnProbType[curRegionRxn]);
			exit(EXIT_FAILURE);
	}
	
	return rxnProb;
}

// Calculate probability of desorption reaction for specified time step
bool calculateDesorptionProb(double * rxnProb,
	const short curRegion,
	const unsigned short curMolType,
	const unsigned short curRegionRxn,
	const double dt,
	const short NUM_REGIONS,
	const struct region regionArray[],
	const unsigned short NUM_MOL_TYPES)
{
	unsigned short curProd;
	double kPrime, kminus1Prime;
	complex double c1, c2;
	
	switch(regionArray[curRegion].rxnProbType[curRegionRxn])
	{
		case RXN_PROB_NORMAL:
		case RXN_PROB_MIXED:
			// No need to calculate this probability separately
			* rxnProb = 0.;
			return false;
		case RXN_PROB_STEADY_STATE:
			// Use steady state reaction probabilities
			if(regionArray[curRegion].bReversible[curRegionRxn])
			{ // S.S. Andrews Physical Biology 2009 Eq. 37
				curProd = regionArray[curRegion].productID[curRegionRxn][0];
				
				* rxnProb = calculateAbsorptionProb(curRegion,
					curProd, regionArray[curRegion].reverseRxnID[curRegionRxn],
					dt, NUM_REGIONS, regionArray, NUM_MOL_TYPES)
					* regionArray[curRegion].rxnRate[curRegionRxn]
					/ regionArray[curRegion].rxnRate[regionArray[curRegion].reverseRxnID[curRegionRxn]]
					*sqrt(regionArray[curRegion].rxnDiffCoef[curRegionRxn]*dt/PI);
				return true;
			} else
			{ // No need to calculate this probability separately
				* rxnProb = 0.;
				return false;
			}
		default:
			fprintf(stderr, "ERROR: Chemical reaction %u of region %u (Label: \"%s\") has invalid desorbing reaction probability type %u.\n",
				curRegionRxn, curRegion, regionArray[curRegion].spec.label,
				regionArray[curRegion].rxnProbType[curRegionRxn]);
			exit(EXIT_FAILURE);
	}
}

// Calculate probability of membrane transition for specified time step
double calculateMembraneProb(const short curRegion,
	const unsigned short curMolType,
	const unsigned short curRegionRxn,
	const double dt,
	const short NUM_REGIONS,
	const struct region regionArray[],
	const unsigned short NUM_MOL_TYPES)
{
	double kFor, kBack, kSum;
	double rxnProb;
	
	switch(regionArray[curRegion].rxnProbType[curRegionRxn])
	{
		case RXN_PROB_NORMAL:
			// Use independent 1st order reaction rate probability
			// Not accurate for surface transitions
			rxnProb = (1 - exp(-dt
				*regionArray[curRegion].rxnRate[curRegionRxn]));
			break;
		case RXN_PROB_MIXED:
			// Use well-mixed reaction probability
			// S.S. Andrews Physical Biology 2009 Eq. 1
			rxnProb = regionArray[curRegion].rxnRate[curRegionRxn] *
				sqrt(PI*dt/regionArray[curRegion].rxnDiffCoef[curRegionRxn]);
			break;
		case RXN_PROB_STEADY_STATE:
			// Use steady state reaction probabilities
			kFor = regionArray[curRegion].rxnRate[curRegionRxn]*
				sqrt(dt/regionArray[curRegion].rxnDiffCoef[curRegionRxn]/2);
			if(regionArray[curRegion].bReversible[curRegionRxn])
			{ // S.S. Andrews Physical Biology 2009 Eq. 47
				kBack =
					regionArray[curRegion].rxnRate[regionArray[curRegion].reverseRxnID[curRegionRxn]]*
					sqrt(dt/regionArray[curRegion].rxnDiffCoef[curRegionRxn]/2);
				kSum = kFor + kBack;
				rxnProb = kFor/kSum/kSum*(2*kSum - sqrt(PI/2)
					+ sqrt(PI/2)*erfcx(sqrt(2)*kSum));
			} else
			{ // S.S. Andrews Physical Biology 2009 Eq. 21
				rxnProb =
					kFor*sqrt(2*PI) - 3.33321*kFor*kFor
					+ 3.35669*kFor*kFor*kFor
					- 1.52092*kFor*kFor*kFor*kFor;
			}
			break;
		default:
			fprintf(stderr, "ERROR: Chemical reaction %u of region %u (Label: \"%s\") has invalid membrane reaction probability type %u.\n",
				curRegionRxn, curRegion, regionArray[curRegion].spec.label,
				regionArray[curRegion].rxnProbType[curRegionRxn]);
			exit(EXIT_FAILURE);
	}
	
	return rxnProb;
}

// Test A Priori surface reaction(s) for given molecule
bool testApmcRxn(const double oldPoint[3],
	double newPoint[3],
	const short curRegion,
	short * newRegion,
	unsigned short * newRegionRxn,
	const unsigned short curMolType,
	const double dt,
	const short NUM_REGIONS,
	const struct region regionArray[],
	const unsigned short NUM_MOL_TYPES,
	const struct chem_rxn_struct chem_rxn[],
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	unsigned short * curGlobalRxn)
{
	unsigned short curApmcRxn = 0; // Current reaction indexed in apmcRxnID
	//unsigned short curGlobalRxn = 0; // Current reaction indexed in global reaction list
	double curProb;
	unsigned short curSurfRegion;
	double dist; // Distance from point to current surface
	double curRand;
	
	// Calculate probabilities associated with each prospective reaction
	regionArray[curRegion].uniCumProbApmc[curMolType][0] = 0;
	for(curApmcRxn = 0; curApmcRxn < regionArray[curRegion].numApmcRxn[curMolType];
		curApmcRxn++)
	{
		if (curApmcRxn > 0)
			regionArray[curRegion].uniCumProbApmc[curMolType][curApmcRxn] =
				regionArray[curRegion].uniCumProbApmc[curMolType][curApmcRxn-1];
		
		// Surface region with current A Priori surface reaction
		curSurfRegion = regionArray[curRegion].apmcRxnRegion[curMolType][curApmcRxn];
		*curGlobalRxn = regionArray[curRegion].apmcGlobalRxnID[curMolType][curApmcRxn];
		
		// Calculate distance to surface
		dist = distanceToBoundary(oldPoint, regionArray[curSurfRegion].spec.shape,
			regionArray[curSurfRegion].boundary);
		
		if(chem_rxn[*curGlobalRxn].rxnThresholdType == RXN_THRESHOLD_DISTANCE
			&& dist > chem_rxn[*curGlobalRxn].rxnThreshold)
			continue; // This reaction is beyond threshold distance
		
		// Calculate (ideal and independent) absorption probability
		switch(chem_rxn[*curGlobalRxn].rxnProbType)
		{
			case RXN_PROB_A_PRIORI_SPHERE:
				// TODO: Add version with flow
				// TODO: Add reference to equation source
				if (chem_rxn[*curGlobalRxn].k < Inf)
				{ // Schulten and Kosztin, "Lectures in Theoretical Biophysics", Eq. (3.114)					
					curProb = 0; // Placeholder (TODO: Need to resolve numerical issues here)
				} else
				{ // Schulten and Kosztin, "Lectures in Theoretical Biophysics", Eq. (3.116)
					curProb = regionArray[curSurfRegion].boundary[3] /
						(dist + regionArray[curSurfRegion].boundary[3]) *
						erfc(dist/sqrt(4*DIFF_COEF[curRegion][curMolType]*dt));
				}
				break;
			case RXN_PROB_A_PRIORI_INFINITE_PLANE:
				// TODO: Find flow along projection from point to boundary
				// TODO: Add reference to equation source
				curProb = erfc(dist/sqrt(4*DIFF_COEF[curRegion][curMolType]*dt));
				break;
			default:
				fprintf(stderr, "ERROR: Chemical reaction %u has invalid A Priori reaction probability type %u.\n",
					*curGlobalRxn, chem_rxn[*curGlobalRxn].rxnProbType);
				exit(EXIT_FAILURE);			
		}
		
		if(chem_rxn[*curGlobalRxn].rxnThresholdType == RXN_THRESHOLD_PROB
			&& curProb < chem_rxn[*curGlobalRxn].rxnThreshold)
			continue; // This reaction is less likely than threshold probability
			
		regionArray[curRegion].uniCumProbApmc[curMolType][curApmcRxn] += curProb;
	}
	
	// Correct cumulative probabilities if the total sum is greater than 1
	if (regionArray[curRegion].uniCumProbApmc[curMolType][curApmcRxn-1] > 1.)
	{
		for(curApmcRxn = 0; curApmcRxn < regionArray[curRegion].numApmcRxn[curMolType];
			curApmcRxn++)
		{
			regionArray[curRegion].uniCumProbApmc[curMolType][curApmcRxn] /=
				regionArray[curRegion].uniCumProbApmc[curMolType][regionArray[curRegion].numApmcRxn[curMolType]-1];
		}
	}
	
	// Test whether a reaction occurred
	curRand = generateUniform();
	for(curApmcRxn = 0; curApmcRxn < regionArray[curRegion].numApmcRxn[curMolType];
		curApmcRxn++)
	{
		if(curRand < regionArray[curRegion].uniCumProbApmc[curMolType][curApmcRxn])
		{ // Current reaction took place
			*curGlobalRxn = regionArray[curRegion].apmcGlobalRxnID[curMolType][curApmcRxn];
			*newRegionRxn = regionArray[curRegion].apmcRxnID[curMolType][curApmcRxn];
			*newRegion = regionArray[curRegion].apmcRxnRegion[curMolType][curApmcRxn];
			
			// Place reacting molecule at point on reacting surface that is closest
			// to current molecule location
			closestPoint(oldPoint, newPoint, regionArray[*newRegion].spec.shape,
				regionArray[*newRegion].boundary);			
			return true;
		}
	}
	
	return false;
}