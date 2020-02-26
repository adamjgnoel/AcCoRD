/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * accord.c - main file
 *
 * Last revised for AcCoRD v1.2 (2018-05-30)
 *
 * Revision history:
 *
 * Revision v1.2 (2018-05-30)
 * - corrected limit on the number of actions made by passive actors
 * - fixed typo in comment describing implementation of reactions for "new"
 * microscopic molecules (created in current time step)
 *
 * Revision v1.1 (2016-12-24)
 * - added direction of subvolume neighbors as a standalone 2D array in order
 * to implement fluid flow in the mesoscopic regime
 * - modified meso-to-micro transition algorithm to allow a molecule to immediately
 * reflect back into the mesoscopic regime if it hits a boundary while being placed
 * in the microscopic regime
 *
 * Revision v1.0 (2016-10-31)
 * - enabled local diffusion coefficients
 * - moved mesoscopic structure fields from subvolume struct to meso subvolume struct
 *
 * Revision v0.7.0.1 (public beta, 2016-08-30)
 * - added measurement of simulation runtime to be written to simulation output
 * - corrected passive actor indexing
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - improved placement of molecules from mesoscopic regime into microscopic regime. User
 * can choose between small time step / large subvolume or large time step / small subvolume
 * algorithms via config. Actual placement moved to function in micro_molecule source file
 * - added bimolecular chemical reactions in microscopic regime (via binding and unbinding
 * radii)
 * - modified random number generation. Now use PCG via a separate interface file.
 * - made output of active actor data sequence a user option
 *
 * Revision v0.5.1 (2016-05-06)
 * - updated placement of molecules created by 0th order reactions when they are
 * created at a surface reaction
 * - re-structured how microscopic 1st order reactions are iteratively performed to
 * account for release of products from microscopic regions
 *
 * Revision v0.5 (2016-04-15)
 * - corrected display of simulation end time
 * - added display of initialization start time
 *
 * Revision v0.4.1
 * - improved use and format of error messages
 * - added more simulation parameters to display while initializing a simulation
 * - added copyright notice and code repository to terminal display
 * - added simulation start time and end time to terminal display
 *
 * Revision v0.4
 * - modified use of unsigned long and unsigned long long to uint32_t and uint64_t
 * - corrected placement of molecule from mesoscopic regime to the microscopic regime when
 * the index of the molecule is greater than 0
 * - corrected calculation for number of levels in heaps (both mesoscopic subvolumes and
 * simulation timers). Previous calculation was incorrect if number of nodes was equal to a
 * power of 2.
 * - renamed region structure and array of structures to accommodate naming
 * of new boundary structure
 *
 * Revision v0.3.1
 * - header added
*/

#include <stdio.h>
#include <stdlib.h> // for exit(), malloc
#include <inttypes.h> // for extended integer type macros
#include <stdbool.h> // for C++ bool naming, requires C99
#include <time.h> // For time record keeping
#include <limits.h> // For SHRT_MAX
#include <math.h> // For ceil(), isfin()
#include "rand_accord.h" // For PRNGs
#include "region.h" // for subvolume definitions, operations
#include "subvolume.h" // for subvolume definitions, operations
#include "meso.h" // for subvolume definitions, operations
#include "micro_molecule.h" // for individual molecule definitions, operations
#include "chem_rxn.h" // for chemical reaction definitions and operations
#include "actor.h" // for actor definitions and operations
#include "mol_release.h" // for actor emission structure (linked list)
#include "actor_data.h" // for active actor data structure (linked list)
#include "observations.h" // for observation structure (linked list)
#include "timer_accord.h" // for timer creation and sorting
#include "global_param.h" // for common global parameters
#include "file_io.h" // For I/O with config and output files

const char CONFIG_NAME[] = "accord_config_sample.txt"; // TEMP - will be loaded from input

int main(int argc, char *argv[])
{
	int i, j; // generic indices
	uint64_t k;
	uint32_t numSub; 		// Total number of subvolumes in system
	uint32_t curSub;		// Index of subvolume where next reaction occurs
	uint32_t curSubID; 	// Current subvolume in actor list
	uint32_t destSub;		// Index of destination subvolume (for diffusion)
	uint32_t destMeso;  	// Index of destination in mesoscopic list
	unsigned short curRegion;	// Region of current subvolume
	unsigned short destRegion;	// Region of destination subvolume (for diffusion)
	unsigned short curRegionID; // Current region in actor list
	unsigned short faceDir; 		// Index in direction array to place molecule
								// (from micro to meso)
	unsigned int curRepeat; 	// Current simulation realization
	double tCur; 				// Current overall simulation time
	double tMeso, tMicro; 		// MESO and MICRO regime simulation times
	double point[3];				// Coordinates of new micro molecules created by 0th order rxn
	double curRand;			// Uniform RV generated for MESO to MICRO transition
	double randCoor[3];			// Random relative coordinates for MESO to MICRO transition
	bool bNeedPoint;			// Need to keep looking for a valid micro location
	
	// Timer and progress variables
	time_t timer;
	char timeBuffer[26]; 	// Array for clock time
	struct tm* timeInfo; 	// Stucture of timer info
	unsigned int updateFreq;
	double fracComplete;
	
	printf("AcCoRD (Actor-based Communication via Reaction-Diffusion)\n");
	printf("Version v1.4.2 (2020-02-12)\n");
	printf("Copyright 2016-2020 Adam Noel. All rights reserved.\n");
	printf("Source code at https://github.com/adamjgnoel/AcCoRD\n");
	printf("User Documentation at https://warwick.ac.uk/fac/sci/eng/staff/ajgn/software/accord/\n");
	
	
	time(&timer);
	timeInfo = localtime(&timer);
	strftime(timeBuffer, 26, "%Y-%m-%d %H:%M:%S", timeInfo);
	printf("Starting initialization at %s.\n", timeBuffer);
	
	//
	// STEP 1: Load Configuration
	//
	// Consider using an interactive environment to specify input file, set message verbosity, and check configuration
	printf("Loading configuration parameters ");
	struct simSpec3D spec;
	if (argc > 2)
	{
		printf("from file \"%s\" using seed offset %d\n", argv[1], atoi(argv[2]));
		loadConfig(argv[1], atoi(argv[2]), &spec);
	}
	else if (argc > 1)
	{
		printf("from file \"%s\" using seed offset defined in that file\n", argv[1]);
		printf("NOTE: To specify a different seed offset, call AcCoRD with the offset as the 2nd argument.\n");
		loadConfig(argv[1], 0, &spec);
	}
	else
	{
		printf("from default configuration file \"%s\"\n", CONFIG_NAME);
		printf("NOTE: To specify a different configuration file, call AcCoRD from command line in format ACCORD_EXE PATH_TO_CONFIG\n");
		printf("NOTE: To specify a different seed offset, call AcCoRD with the offset as the 2nd argument.\n");
		loadConfig(CONFIG_NAME, 0, &spec);
	}
	
	uint64_t numMolChange[spec.NUM_MOL_TYPES]; // Number of molecules changed by mesoscopic event
	for(i = 0; i < spec.NUM_MOL_TYPES; i++)
		numMolChange[i] = 0; // Initialize elements to 0
	bool bTrue[spec.NUM_MOL_TYPES]; // Array of all true values
	bool bFalse[spec.NUM_MOL_TYPES];
	for(i = 0; i < spec.NUM_MOL_TYPES; i++)
	{
		bTrue[i] = true;
		bFalse[i] = false;
	}
	
	// MESO parameters
	uint32_t numMesoSub;
	uint32_t curMeso, curDestMeso;
	double new_time, old_prop, old_time;
	double prop_sum, rand_prop;
	unsigned short curRxn;		// Index of current reaction to fire
	uint64_t numMesoSteps;
	unsigned short curMolType;  // Current type of molecule diffusing
	
	// MICRO parameters
	ListMol3D microMolList[spec.NUM_REGIONS][spec.NUM_MOL_TYPES];
	ListMolRecent3D microMolListRecent[spec.NUM_REGIONS][spec.NUM_MOL_TYPES]; // For molecules created at
															   // at an arbitrary time
	for(i = 0; i < spec.NUM_REGIONS; i++)
	{
		for(j = 0; j < spec.NUM_MOL_TYPES; j++)
		{
			initializeListMol(&microMolList[i][j]);
			initializeListMolRecent(&microMolListRecent[i][j]);
		}
	}
	double micro_sigma[spec.NUM_REGIONS][spec.NUM_MOL_TYPES];
	uint64_t curMol;
	uint64_t numMol;
	unsigned short curZerothRxn, curFirstRxn, curSecondRxn;
	bool bCheckCount;
	uint32_t numMicroMolCheck[spec.NUM_REGIONS][spec.NUM_MOL_TYPES];
	uint32_t sumMicroMolCheck;
	
	// Hybrid parameters
	uint32_t curBoundSub;  // Index of subvolume in region boundary list
	
	// Actor parameters
	short NUM_ACTORS_ACTIVE, NUM_ACTORS_PASSIVE;
	short curActor, curPassive, curActive;
	unsigned short curMolPassive;
	short numPassiveRecord, curPassiveRecord;
	short numActiveRecord, curActiveRecord;
	short * passiveRecordID; // Array of IDs of passive actors whose observations are recorded
	short * activeRecordID; // Array of IDs of active actors whose observations are recorded
						   // TODO: Should these IDs be from actor list or passive list?
	
	clock_t startTime, endTime; // Time record keeping
	double runTime; 				// Runtime for realizations
	
	double DIFF_COEF [spec.NUM_REGIONS][spec.NUM_MOL_TYPES];
	for(i = 0; i < spec.NUM_REGIONS; i++)
	{
		if(spec.subvol_spec[i].bLocalDiffusion)
		{
			for(j = 0; j < spec.NUM_MOL_TYPES; j++)
			{
				DIFF_COEF[i][j] = spec.subvol_spec[i].diffusion[j];
			}
		} else
		{
			for(j = 0; j < spec.NUM_MOL_TYPES; j++)
			{
				DIFF_COEF[i][j] = spec.DIFF_COEF[j];
			}
		}
	}
		
	//
	// STEP 2: Validate Configuration
	//
	
	//
	// STEP 3: Initialize Simulation Environment
	//
	
	//
	// 3-A Initialize Mesoscopic Environment
	//
	
	// Initialize array of region information
	printf("Initializing region parameters.\n");
	printf("Number of regions: %u\n", spec.NUM_REGIONS);
	struct region regionArray[spec.NUM_REGIONS];
	initializeRegionArray(regionArray,	spec.subvol_spec, spec.NUM_REGIONS,
		spec.NUM_MOL_TYPES,	spec.SUBVOL_BASE_SIZE, DIFF_COEF,
		spec.MAX_RXNS, spec.chem_rxn);
	
	
	// Define subvolume array
	printf("Initializing microscopic and mesoscopic subvolume parameters.\n");
	numSub = countAllSubvolumes(regionArray, spec.NUM_REGIONS);
	printf("Number of subvolumes: %" PRIu32 "\n", numSub);
	struct subvolume3D * subvolArray;
	allocateSubvolArray(numSub,&subvolArray);
	
	// Allocate temporary arrays for managing subvolume validity and placement
	uint32_t (* subCoorInd)[3]; // List (within region) index coordinates for each subvolume
	uint32_t **** subID; // Subvolume indices in master list
	uint32_t (** subIDSize)[2]; // Reader array for subID	
	unsigned short ** subNeighDir; // Direction of each of a subvolume's neighbors
	allocateSubvolHelper(numSub, &subCoorInd, &subID, &subIDSize, &subNeighDir,
		spec.NUM_REGIONS, regionArray);
	
	// Delete temporary arrays for managing subvolume validity and placement
	//deleteSubvolHelper(subCoorInd, subID, subIDSize, spec.NUM_REGIONS, regionArray);
	//allocateSubvolHelper(numSub, &subCoorInd, &subID, &subIDSize, spec.NUM_REGIONS, regionArray);
	
	// Build subvolume array from subvolume specifications
	buildSubvolArray(numSub, &numMesoSub, subvolArray,
		spec.subvol_spec, regionArray,
		spec.NUM_REGIONS, spec.NUM_MOL_TYPES,
		spec.MAX_RXNS, spec.SUBVOL_BASE_SIZE,
		DIFF_COEF, subCoorInd, subID, subIDSize, subNeighDir);
		
	// Determine rates associated with chemical reaction events
	// Record dependencies (i.e., if a given rxn fires, what propensities need to be updated?)
	// 0th order - k*vol (indp of number of molecules)
	// 1st order - k*A (indp of volume)
	// 2nd order A+B - k*A*B/vol
	// 2nd order A+A - k*A*(A-1)/vol
	
	// Build array of mesoscopic subvolumes (with parameters only needed f
	printf("Number of mesoscopic subvolumes: %" PRIu32 "\n", numMesoSub);
	printf("Initializing mesoscopic subvolumes and meso reaction heap...\n");
	struct mesoSubvolume3D * mesoSubArray;
	allocateMesoSubArray(numMesoSub,&mesoSubArray);
	initializeMesoSubArray(numMesoSub, numSub, mesoSubArray, subvolArray,
		spec.SUBVOL_BASE_SIZE, spec.NUM_MOL_TYPES, spec.MAX_RXNS, regionArray,
		spec.NUM_REGIONS, subCoorInd, subNeighDir, DIFF_COEF);
	
	// Build heap for mesoscopic subvolumes and associated arrays
	uint32_t * heap_subvolID;
	uint32_t (*heap_childID)[2];
	bool (*b_heap_childValid)[2];
	allocateMesoHeapArray(numMesoSub, &heap_subvolID, &heap_childID, &b_heap_childValid);
	
	heapMesoFindChildren(numMesoSub, heap_childID, b_heap_childValid);
	unsigned int num_heap_levels = (unsigned int) ceil(log2(numMesoSub+1));
	
	// Build actor array
	printf("Initializing simulation actors...\n");
	printf("Number of actors: %u\n", spec.NUM_ACTORS);
	struct actorStruct3D * actorCommonArray;
	struct actorActiveStruct3D * actorActiveArray = NULL;
	struct actorPassiveStruct3D * actorPassiveArray = NULL;
	allocateActorCommonArray(spec.NUM_ACTORS, &actorCommonArray);
	initializeActorCommon(spec.NUM_ACTORS, actorCommonArray,
		spec.actorSpec,
		regionArray, spec.NUM_REGIONS,
		&NUM_ACTORS_ACTIVE, &numActiveRecord, &activeRecordID,
		&NUM_ACTORS_PASSIVE, &numPassiveRecord, &passiveRecordID,
		subvolArray, subID, subCoorInd, spec.SUBVOL_BASE_SIZE);
	allocateActorActivePassiveArray(NUM_ACTORS_ACTIVE, &actorActiveArray,
		NUM_ACTORS_PASSIVE, &actorPassiveArray);
	initializeActorActivePassive(spec.NUM_ACTORS, actorCommonArray,
		spec.NUM_MOL_TYPES,
		regionArray, spec.NUM_REGIONS, mesoSubArray, NUM_ACTORS_ACTIVE, actorActiveArray,
		NUM_ACTORS_PASSIVE, actorPassiveArray, subCoorInd);
	printf("Number of active actors: %u\n", NUM_ACTORS_ACTIVE);
	printf("Number of passive actors: %u\n", NUM_ACTORS_PASSIVE);
	
	// Delete temporary arrays for managing subvolume validity and placement
	deleteSubvolHelper(subCoorInd, subID, subIDSize, subNeighDir,
		spec.NUM_REGIONS, regionArray, numSub);
	
	ListMol3D molListPassive3D[spec.NUM_MOL_TYPES];
	for(j = 0; j < spec.NUM_MOL_TYPES; j++)
	{
		initializeListMol(&molListPassive3D[j]);
	}
	
	// Create arrays to store the maximum number of bits of each
	// active actor and each recorded passive actor.
	// Will be appended to output file to assist importing into Matlab
	uint32_t maxActiveBits[numActiveRecord];
	uint32_t maxPassiveObs[numPassiveRecord];
	
	// Create array of linked lists for recording actor observations
	ListObs3D observationArray[numPassiveRecord];
	for(curActor = 0; curActor < numPassiveRecord; curActor++)
	{
		maxPassiveObs[curActor] = 0;
		initializeListObs(&observationArray[curActor],
			actorPassiveArray[actorCommonArray[passiveRecordID[curActor]].passiveID].numMolRecordID);
	}
	for(curActor = 0; curActor < numActiveRecord; curActor++)
	{
		maxActiveBits[curActor] = 0;
	}
	
	// Create timer heap	
	short NUM_TIMERS =
		spec.NUM_ACTORS + 1 + 1; 	// NUM_ACTORS + (ANY MESO?) + (ANY MICRO?)
	short * heapTimer;           		   	// Heap of timer IDs
	struct timerStruct * timerArray;		// timer values for simulation
	short MESO_TIMER_ID = spec.NUM_ACTORS; 	// ID of meso timer in timer array
	short MICRO_TIMER_ID =
		spec.NUM_ACTORS + 1; 	// ID of (first) micro timer in timer array
	short curTimer; 					   	// ID of current timer
	bool bRecordPos;
	short (*heapTimerChildID)[2];
	bool (*b_heapTimerChildValid)[2];
	allocateTimerArray(NUM_TIMERS, &timerArray);
	initializeTimerArray(NUM_TIMERS, timerArray);
	allocateTimerHeapArray(NUM_TIMERS, &heapTimer,
		&heapTimerChildID, &b_heapTimerChildValid);	
	heapTimerFindChildren(NUM_TIMERS, heapTimerChildID, b_heapTimerChildValid);
	unsigned int numHeapTimerLevels = (unsigned int) ceil(log2(NUM_TIMERS+1));
	
	// Open output text file	
	FILE * out, * outSummary;

	if (argc > 1)
		initializeOutput(&out, &outSummary, argv[1], spec);
	else
		initializeOutput(&out, &outSummary, CONFIG_NAME, spec);
	
	//
	// 3-B Initialize Microscopic Environment
	//
	
	// Create a list for every type of molecule in each microscopic region
	for(i = 0; i < spec.NUM_REGIONS; i++)
	{
		for(j = 0; j < spec.NUM_MOL_TYPES; j++)
		{
			micro_sigma[i][j] = sqrt(2.*spec.DT_MICRO*DIFF_COEF[i][j]);
		}
	}
		
	// Initialize random number generation
	printf("Starting up random number generator with seed offset: %u\n", spec.SEED);
	rngInitialize(spec.SEED);
	
	//
	// STEP 4: Run Simulation
	//
	
	// Initialize variables to track simulation progress
	updateFreq = (unsigned int)	ceil((double) spec.NUM_REPEAT / spec.MAX_UPDATES);
	time(&timer);
	timeInfo = localtime(&timer);
	strftime(timeBuffer, 26, "%Y-%m-%d %H:%M:%S", timeInfo);
	
	printf("Starting simulation at %s.\n", timeBuffer);
	startTime = clock();
	for(curRepeat = 0; curRepeat < spec.NUM_REPEAT; curRepeat++)
	{
		
		// Initialize current realization
		
		//
		// Mesoscopic Initialization
		//
		
		// Perform initial placement of molecules
		// FILL IN mesoscopic subvolumes with 0 molecules
		for(i = 0; i < numMesoSub; i++)
		{
			for(j = 0; j < spec.NUM_MOL_TYPES; j++)
				mesoSubArray[i].num_mol[j] = 0ULL;
		}
		
		// Determine Reaction propensities based on subvolArray[].num_mol
		resetMesoSubArray(numMesoSub, mesoSubArray, subvolArray,
			spec.NUM_MOL_TYPES, spec.MAX_RXNS, spec.NUM_REGIONS, regionArray);
		// Initialize heap for next subvolume method
		heapMesoBuild(numMesoSub, mesoSubArray, heap_subvolID, num_heap_levels,
			heap_childID, b_heap_childValid);
		
		//
		// Microscopic Initialization
		//
		for(i = 0; i < spec.NUM_REGIONS; i++)
		{
			for(j = 0; j < spec.NUM_MOL_TYPES; j++)
			{
				if(!isListMol3DEmpty(&microMolList[i][j]))
				{ // Molecule list is not empty
					// Free memory and re-allocate empty list
					emptyListMol(&microMolList[i][j]);
					initializeListMol(&microMolList[i][j]);
				}
				if(!isListMol3DRecentEmpty(&microMolListRecent[i][j]))
				{ // Molecule list is not empty
					// Free memory and re-allocate empty list
					emptyListMol3DRecent(&microMolListRecent[i][j]);
					initializeListMolRecent(&microMolListRecent[i][j]);
				}
			}
		}
		// Reset start times for zeroth order reactions in micro regions
		for(curRegion = 0; curRegion < spec.NUM_REGIONS; curRegion++)
		{
			if(!regionArray[curRegion].spec.bMicro)
				continue;

			for(curZerothRxn = 0; curZerothRxn < regionArray[curRegion].numZerothRxn; curZerothRxn++)
				regionArray[curRegion].tZeroth[curZerothRxn] =
					generateExponential(1)/regionArray[curRegion].rxnRateZerothMicro[curZerothRxn];
		}
		
		//
		// Actor initialization
		//
		
		// Reset actor simulation parameters
		resetActors(spec.NUM_ACTORS, actorCommonArray,
			spec.NUM_MOL_TYPES, regionArray,
			spec.NUM_REGIONS, NUM_ACTORS_ACTIVE, actorActiveArray, NUM_ACTORS_PASSIVE,
			actorPassiveArray);
			
		// Reset timer array
		if(numMesoSub > 0)
			tMeso = mesoSubArray[heap_subvolID[0]].t_rxn; // Time of next event in MESO regime
		else
			tMeso = INFINITY;
		resetTimerArray(NUM_TIMERS, timerArray, spec.NUM_ACTORS, actorCommonArray,
			MESO_TIMER_ID, tMeso, MICRO_TIMER_ID, spec.DT_MICRO);
		// Initialize heap for timers
		heapTimerBuild(NUM_TIMERS, timerArray, heapTimer, numHeapTimerLevels,
			heapTimerChildID, b_heapTimerChildValid);
	
		// Reset observation lists
		for(curActor = 0; curActor < numPassiveRecord; curActor++)
		{
			if(!isListEmptyObs(&observationArray[curActor]))
			{
				emptyListObs(&observationArray[curActor]);
				initializeListObs(&observationArray[curActor],
					actorPassiveArray[actorCommonArray[passiveRecordID[curActor]].passiveID].numMolRecordID);				
			}
		}
		
		// Initialize runtime parameters
		tCur = 0.;			// TODO: Allow definition by user. Should be min(startTime)
		tMicro = spec.DT_MICRO; // Time of next event in MICRO regime
		
		numMesoSteps = 0ULL;
		
		while(timerArray[heapTimer[0]].nextTime <= spec.TIME_FINAL)
		{	
			if (numMesoSub > 0)
			{
				// Update meso timer in timer heap
				heapTimerUpdate(NUM_TIMERS, timerArray, heapTimer, timerArray[MESO_TIMER_ID].heapID,
					heapTimerChildID, b_heapTimerChildValid);
			}
			
			// Determine the next type of step in the simulation
			if (heapTimer[0] < spec.NUM_ACTORS)
			{	// Next step is by an Actor
		
				tCur = timerArray[heapTimer[0]].nextTime;
		
				if(actorCommonArray[heapTimer[0]].spec.bActive)
				{	// Actor is active. Place molecules or create a release object (which
					// will place molecules) as specified
					
					curActive = actorCommonArray[heapTimer[0]].activeID;
					
					// Is next action the start of a new release?
					if (actorActiveArray[curActive].bNextActionNewRelease)
					{ // Determine the parameters of the new release
						actorCommonArray[heapTimer[0]].curAction++;
						
						newRelease(&actorCommonArray[heapTimer[0]],
							&actorActiveArray[curActive], timerArray[heapTimer[0]].nextTime);
						if (actorCommonArray[heapTimer[0]].spec.bIndependent)
						{
							if(!actorCommonArray[heapTimer[0]].spec.bMaxAction
								|| actorCommonArray[heapTimer[0]].curAction < actorCommonArray[heapTimer[0]].spec.numMaxAction)
								actorActiveArray[curActive].nextNewReleaseTime +=
									actorCommonArray[heapTimer[0]].spec.actionInterval;
							else
								actorActiveArray[curActive].nextNewReleaseTime = INFINITY;
						} else
						{
							actorActiveArray[curActive].nextNewReleaseTime = INFINITY;
						}
					} else
					{ // Next action is the release of molecules from a current release
						fireEmission(&actorCommonArray[heapTimer[0]],
							&actorActiveArray[curActive], regionArray,
							spec.NUM_REGIONS, subvolArray, mesoSubArray,
							numMesoSub, spec.NUM_MOL_TYPES, 
							microMolListRecent, tMicro,
							heap_subvolID, heap_childID,
							b_heap_childValid);
						// Updated mesoscopic time
						if (numMesoSub > 0)
						{
						timerArray[MESO_TIMER_ID].nextTime =
							mesoSubArray[heap_subvolID[0]].t_rxn;
						tMeso = mesoSubArray[heap_subvolID[0]].t_rxn;
						}
					}
					
					// Update next time associated with actor
					if(actorActiveArray[curActive].nextNewReleaseTime <
						actorActiveArray[curActive].nextEmissionTime)
					{ // Next event will be a new release
						actorActiveArray[curActive].bNextActionNewRelease = true;
						timerArray[heapTimer[0]].nextTime =
							actorActiveArray[curActive].nextNewReleaseTime;
					} else
					{ // Next event will be an emission from a current release
						actorActiveArray[curActive].bNextActionNewRelease = false;
						timerArray[heapTimer[0]].nextTime =
							actorActiveArray[curActive].nextEmissionTime;
					}
					actorCommonArray[heapTimer[0]].nextTime =
						timerArray[heapTimer[0]].nextTime;
					
				} else
				{	// Actor is passive. Make required observations as specified
					
					actorCommonArray[heapTimer[0]].curAction++;
					
					// Initialize molecule list for coordinates
					for(j = 0; j < spec.NUM_MOL_TYPES; j++)
					{
						initializeListMol(&molListPassive3D[j]);
					}
					
					curPassive = actorCommonArray[heapTimer[0]].passiveID;
					curPassiveRecord = actorPassiveArray[curPassive].recordID;
					for(curMolPassive = 0;
						curMolPassive < actorPassiveArray[curPassive].numMolRecordID;
						curMolPassive++)
					{	// Determine type of molecule being observed
						curMolType =
							actorPassiveArray[curPassive].molRecordID[curMolPassive];
						// Will molecule coordinates be recorded
						bRecordPos =
							actorCommonArray[heapTimer[0]].spec.bRecordPos[curMolType];
												
						actorPassiveArray[curPassive].curMolObs[curMolPassive] = 0ULL;
						// Search for molecules in each region in actor
						for(curRegionID = 0;
							curRegionID < actorCommonArray[heapTimer[0]].numRegion;
							curRegionID++)
						{
							curRegion =
								actorCommonArray[heapTimer[0]].regionID[curRegionID];
							if(regionArray[curRegion].spec.bMicro)
							{	// Search through region's molecule list
								actorPassiveArray[curPassive].curMolObs[curMolPassive] +=
									recordMolecules(&microMolList[curRegion][curMolType], &molListPassive3D[curMolPassive], actorCommonArray[heapTimer[0]].regionInterType[curRegionID], actorCommonArray[heapTimer[0]].regionInterBound[curRegionID], bRecordPos,
									actorCommonArray[heapTimer[0]].bRegionInside[curRegionID]) +
									recordMoleculesRecent(&microMolListRecent[curRegion][curMolType], &molListPassive3D[curMolPassive], actorCommonArray[heapTimer[0]].regionInterType[curRegionID], actorCommonArray[heapTimer[0]].regionInterBound[curRegionID], bRecordPos,
									actorCommonArray[heapTimer[0]].bRegionInside[curRegionID]);
							} else
							{	// Search through subvolumes inside actor
								for(curSubID = 0;
									curSubID < actorCommonArray[heapTimer[0]].numSub[curRegionID];
									curSubID++)
								{
									curSub = actorCommonArray[heapTimer[0]].subID[curRegionID][curSubID];
									numMol = mesoSubArray[curSub].num_mol[curMolType];
									if(actorPassiveArray[curPassive].fracSubInActor[curRegionID][curSubID] < 1.)
									{ // Roll die to see whether each molecule is in actor
										for(curMol = 0;
											curMol < numMol;
											curMol++)
										{
											// Roll die to see whether molecule is within actor
											if(generateUniform() < actorPassiveArray[curPassive].fracSubInActor[curRegionID][curSubID])
											{
												actorPassiveArray[curPassive].curMolObs[curMolPassive]++;
												
												// Add molecule position
												if(bRecordPos)
												{
													uniformPointVolume(point, regionArray[curRegion].subShape,
													actorPassiveArray[curPassive].subInterBound[curRegionID][curSubID], false, 0);
													if(!addMolecule(&molListPassive3D[curMolPassive],
														point[0], point[1], point[2]))
													{ // Creation of molecule failed
														fprintf(stderr,"ERROR: Memory allocation to create molecule %"
														PRIu64 " of %" PRIu64" of type %u being observed by actor %u in region %u.\n",
														curMol, numMol, curMolType, heapTimer[0], curRegion);
														exit(EXIT_FAILURE);
													}
												}
											}
										}
									} else
									{ // All molecules in subvolume are within actor
										actorPassiveArray[curPassive].curMolObs[curMolPassive] += numMol;
										
										// Add molecule position
										if(bRecordPos)
										{
											for(curMol = 0;
												curMol < numMol;
												curMol++)
											{
												uniformPointVolume(point, regionArray[curRegion].subShape,
													actorPassiveArray[curPassive].subInterBound[curRegionID][curSubID], false, 0);
												if(!addMolecule(&molListPassive3D[curMolPassive],
													point[0], point[1], point[2]))
												{ // Creation of molecule failed
													fprintf(stderr,"ERROR: Memory allocation to create molecule %"
														PRIu64 " of %" PRIu64" of type %u being created by actor %u in region %u.\n",
														curMol, numMol, curMolType, heapTimer[0], curRegion);
													exit(EXIT_FAILURE);
												}
											}
										}
										
									}
								} // Search over subvolumes in region
							}
								
						} // Search over regions in actor
					}
					
					if(curPassiveRecord < SHRT_MAX)
					{ // Add observation data to the observation list
						addObservation(&observationArray[curPassiveRecord],
							((actorCommonArray[heapTimer[0]].spec.bRecordTime)? 1 : 0),
							actorPassiveArray[curPassive].numMolRecordID,
							&actorCommonArray[heapTimer[0]].nextTime,
							actorPassiveArray[curPassive].curMolObs,
							molListPassive3D);
					}
					// Empty molecule list for coordinates
					for(j = 0; j < spec.NUM_MOL_TYPES; j++)
					{
						emptyListMol(&molListPassive3D[j]);
					}
					
					// Update timer structure array
					if (actorCommonArray[heapTimer[0]].spec.bIndependent)
					{	// Actor is independent. Next action time is known
						if(!actorCommonArray[heapTimer[0]].spec.bMaxAction
							|| actorCommonArray[heapTimer[0]].curAction < actorCommonArray[heapTimer[0]].spec.numMaxAction)
						{
							actorCommonArray[heapTimer[0]].nextTime +=
								actorCommonArray[heapTimer[0]].spec.actionInterval;
							timerArray[heapTimer[0]].nextTime +=
								actorCommonArray[heapTimer[0]].spec.actionInterval;
						}
						else
						{
							actorCommonArray[heapTimer[0]].nextTime = INFINITY;
							timerArray[heapTimer[0]].nextTime = INFINITY;
						}
					} else
					{	// Actor is dependent. Next action time is unknown
						actorCommonArray[heapTimer[0]].nextTime = INFINITY;
						timerArray[heapTimer[0]].nextTime = INFINITY;
					}
				}
						
			} else if (heapTimer[0] > spec.NUM_ACTORS)
			{ 	// Next step is in Micro regime
				
				// Update Overall Time
				tCur = tMicro;
				
				// Execute zeroth order reactions to generate new molecules
				// and add to list of recently-created molecules
				// TODO: Move into a chem_rxn.c or micro_molecule.c function
				for(i = 0; i < spec.NUM_REGIONS; i++)
				{
					if(!regionArray[i].spec.bMicro)
						continue;
						
					for(curZerothRxn = 0; curZerothRxn < regionArray[i].numZerothRxn; curZerothRxn++)
					{ // For current 0th order reaction in this region
						while(regionArray[i].tZeroth[curZerothRxn] < tCur)
						{ // Next reaction time was before end of last Micro time step
							
							// Determine reaction index in master region reaction list
							curRxn = regionArray[i].zerothRxn[curZerothRxn];
							
							// Create molecules as specified by
							// regionArray[i].numMolChange[curRxn]
							// (In zeroth order, molecules can only be created)
							for(j = 0; j < spec.NUM_MOL_TYPES; j++)
							{
								for(k = 0; k < regionArray[i].numMolChange[curRxn][j]; k++)
								{
									generatePointInRegion(i, regionArray, point);
									if(regionArray[i].bReleaseProduct[curRxn][j])
									{ // Region is a surface and molecule must be released
										destRegion = findDestRegion(point, i, regionArray);
									} else
									{ // Place product molecule in current region
										destRegion = i;
									}
									if(!addMoleculeRecent(&microMolListRecent[destRegion][j], point[0],
										point[1], point[2],
										tMicro - regionArray[i].tZeroth[curZerothRxn]))
									{ // Creation of molecule failed
										fprintf(stderr,"ERROR: Memory allocation to create molecule %"
											PRIu64 " of %" PRIu64" of type %u being created by reaction %u in region %u and destination region %u.\n",
											k, regionArray[i].numMolChange[curRxn][j], j, curRxn, i, destRegion);
									}
								}									
							}
							
							regionArray[i].tZeroth[curZerothRxn] +=
								generateExponential(1)/regionArray[i].rxnRateZerothMicro[curZerothRxn];
						}
					}
				}
				
				// Execute first order reactions on all molecules
				// (both recently-created and "old")
				for(i = 0; i < spec.NUM_REGIONS; i++)
				{
					if(!regionArray[i].spec.bMicro || regionArray[i].numFirstRxn < 1)
						continue;
					
					for(j = 0; j < spec.NUM_MOL_TYPES; j++)
					{ // For current molecule type in this region
					
						if(!isListMol3DEmpty(&microMolList[i][j])
							&& regionArray[i].numFirstRxnWithReactant[j] > 0)
						{ // Check 1st order reactions of "old" molecules
							rxnFirstOrder(spec.NUM_REGIONS, spec.NUM_MOL_TYPES, i,
								microMolList, regionArray, j,
								DIFF_COEF, microMolListRecent);
						}
						numMicroMolCheck[i][j] = 0;
					}
				}
					
				sumMicroMolCheck = 1;
				bCheckCount = false;
				while(sumMicroMolCheck > 0)
				{
					for(i = 0; i < spec.NUM_REGIONS; i++)
					{
						if(!regionArray[i].spec.bMicro || regionArray[i].numFirstRxn < 1)
							continue;
						
						for(j = 0; j < spec.NUM_MOL_TYPES; j++)
						{ // For current molecule type in this region
						
							if(!isListMol3DRecentEmpty(&microMolListRecent[i][j])
								&& regionArray[i].numFirstRxnWithReactant[j] > 0)
							{ // Check 1st order reactions of "new" molecules
								rxnFirstOrderRecent(spec.NUM_REGIONS,
									spec.NUM_MOL_TYPES, i, microMolListRecent,
									microMolList, regionArray, j, DIFF_COEF,
									bCheckCount, numMicroMolCheck);
							}
						}
						
						for(j = 0; j < spec.NUM_MOL_TYPES; j++)
						{ // For current molecule type in this region
						
							if(regionArray[i].numFirstRxnWithReactant[j] < 1)
							{ 
								numMicroMolCheck[i][j] = 0;
							}
						}
					}
						
					// Find total number of molecules to re-check.
					// We do this in a new for-loop because numMicroMolCheck
					// of any element could change within rxnFirstOrderRecent
					sumMicroMolCheck = 0;
					for(i = 0; i < spec.NUM_REGIONS; i++)
					{
						if(!regionArray[i].spec.bMicro || regionArray[i].numFirstRxn < 1)
							continue;
						for(j = 0; j < spec.NUM_MOL_TYPES; j++)
						{
							sumMicroMolCheck += numMicroMolCheck[i][j];
						}
					}
					bCheckCount = true;
				}
				
				// Diffuse all microscopic molecules to valid locations and merge
				// 2 sets of molecule lists into 1. Also imposes flow.
				diffuseMolecules(spec.NUM_REGIONS, spec.NUM_MOL_TYPES, microMolList,
					microMolListRecent, regionArray, mesoSubArray,
					micro_sigma, spec.chem_rxn, spec.MAX_HYBRID_DIST, DIFF_COEF);
				
				// Execute 2nd Order Reactions
				rxnSecondOrder(spec.NUM_REGIONS, spec.NUM_MOL_TYPES, microMolList,
					regionArray, mesoSubArray, spec.chem_rxn, DIFF_COEF);
				
				if(numMesoSub > 0)
				{
					// Check whether any subvolumes must be updated due to added molecules
					// from microscopic regime
					updateMesoSubBoundary(numSub, numMesoSub, mesoSubArray,
						subvolArray, regionArray, bTrue, tCur, spec.NUM_REGIONS,
						spec.NUM_MOL_TYPES, heap_subvolID, heap_childID, b_heap_childValid);
					
					// Update timer structure array
					timerArray[MESO_TIMER_ID].nextTime =
						mesoSubArray[heap_subvolID[0]].t_rxn;
					tMeso = mesoSubArray[heap_subvolID[0]].t_rxn;
				}
				
				// Update Time of next MICRO event
				tMicro += spec.DT_MICRO;
				
				// Update timer structure array
				timerArray[MICRO_TIMER_ID].nextTime += spec.DT_MICRO;
			} else
			{	// Next step is in Meso regime
				numMesoSteps++;
				// Update Overall Time
				tCur = tMeso;
				
				// Find subvolume where next reaction occurs
				curMeso = heap_subvolID[0];
				curSub = mesoSubArray[curMeso].subID;
				// Determine region of next reaction
				// (needed for diffusion to another subvolume or to
				// check for a valid chemical reaction)
				curRegion = subvolArray[curSub].regionID;
				// Determine ID of next reaction
				curRxn = 0;
				prop_sum = mesoSubArray[curMeso].rxnProp[curRxn];
				rand_prop = generateUniform()*mesoSubArray[curMeso].totalProp;
				while (prop_sum < rand_prop)
				{
					prop_sum += mesoSubArray[curMeso].rxnProp[++curRxn];
				}
				
				// "Fire" reaction event				
				if(curRxn < mesoSubArray[curMeso].firstChemRxn)
				{ // Event is diffusion to a neighbor
				  // Need to determine type of molecule and destination
					if(spec.NUM_MOL_TYPES > 1)
					{ // Need to use division operators to determine destination
						curMolType = curRxn / subvolArray[curSub].num_neigh;
						destSub =
							subvolArray[curSub].neighID[curRxn % subvolArray[curSub].num_neigh];
					} else
					{ // Only one type of molecule
						curMolType = 0;
						destSub = subvolArray[curSub].neighID[curRxn];
					}
					
					// Remove molecule from current subvolume
					mesoSubArray[curMeso].num_mol[curMolType] -= 1ULL;
					numMolChange[0] = 1ULL;
					
					// Update next reaction time and heap structure location
					updateMesoSub(curSub, false, numMolChange,
						bFalse, curMolType, false, numMesoSub,
						mesoSubArray, subvolArray, tCur,
						spec.NUM_REGIONS, spec.NUM_MOL_TYPES,
						regionArray);
					heapMesoUpdate(numMesoSub, mesoSubArray, heap_subvolID, 0UL,
						heap_childID, b_heap_childValid); // ALWAYS head node
					
					// Add molecule to destination
					if(regionArray[subvolArray[destSub].regionID].spec.bMicro)
					{ // Destination is in a microscopic region
						
						// Find ID of current subvolume in region boundary list
						destRegion = subvolArray[destSub].regionID;
						curBoundSub = findSubInBoundaryList(curRegion, destRegion, regionArray,
							curMeso);
						
						if(curBoundSub < UINT32_MAX)
						{ // Add new molecule to random location next to source subvolume
					
							if(!placeInMicroFromMeso(curRegion, spec.NUM_REGIONS,
								spec.NUM_MOL_TYPES, destRegion, &destSub,
								regionArray,
								curBoundSub, spec.B_HYBRID_SMALL_SUB, curMolType,
								microMolListRecent, spec.chem_rxn, DIFF_COEF))
							{ // Molecule "bounced back" to meso
							  // destSub is now the mesoID
								curDestMeso = destSub;
								destSub = mesoSubArray[curDestMeso].subID;
								mesoSubArray[curDestMeso].num_mol[curMolType] += 1ULL;
								numMolChange[0] = 1ULL;
								
								// Update propensities of destination subvolume
								// and its location in the heap
								updateMesoSub(destSub, false, numMolChange, bTrue, curMolType, true,
									numMesoSub,	mesoSubArray, subvolArray, tCur, spec.NUM_REGIONS,
									spec.NUM_MOL_TYPES, regionArray);
								heapMesoUpdate(numMesoSub, mesoSubArray, heap_subvolID,
									mesoSubArray[subvolArray[destSub].mesoID].heapID, heap_childID, b_heap_childValid);
							}
						} else
						{
							// Error
							fprintf(stderr,"ERROR: Molecule was supposed to transition out of subvolume %" PRIu32 " in region %u and into region %u but this subvolume could not be found in the list of subvolumes along the boundary between the two regions.\n", curSub, curRegion, destRegion);
							exit(EXIT_FAILURE);
						}
					} else
					{ // Destination is also a mesoscopic subvolume
						curDestMeso = subvolArray[destSub].mesoID;
						mesoSubArray[curDestMeso].num_mol[curMolType] += 1ULL;
						numMolChange[0] = 1ULL;
						
						// Update propensities of destination subvolume
						// and its location in the heap
						updateMesoSub(destSub, false, numMolChange, bTrue, curMolType, true,
							numMesoSub,	mesoSubArray, subvolArray, tCur, spec.NUM_REGIONS,
							spec.NUM_MOL_TYPES, regionArray);
						heapMesoUpdate(numMesoSub, mesoSubArray, heap_subvolID,
							mesoSubArray[subvolArray[destSub].mesoID].heapID, heap_childID, b_heap_childValid);
							
					}
					
					// Reset numMolChange array for next diffusion event
					numMolChange[0] = 0ULL;
					
				} else if(curRxn < mesoSubArray[curMeso].firstChemRxn
					+ regionArray[curRegion].numChemRxn)
				{ // Event is a chemical reaction
					// Shift reaction ID by number of diffusion events
					curRxn -= mesoSubArray[curMeso].firstChemRxn;
					// TODO: Simplify by storing reactants and products of each reaction
					for(curMolType = 0; curMolType < spec.NUM_MOL_TYPES; curMolType++)
					{
						if(regionArray[curRegion].bMolAdd[curRxn][curMolType])
							mesoSubArray[curMeso].num_mol[curMolType] +=
								regionArray[curRegion].numMolChange[curRxn][curMolType];
						else
							mesoSubArray[curMeso].num_mol[curMolType] -=
								regionArray[curRegion].numMolChange[curRxn][curMolType];
					}
										
					// Update next reaction time and heap structure location
					updateMesoSub(curSub, true, regionArray[curRegion].numMolChange[curRxn],
						regionArray[curRegion].bMolAdd[curRxn], 0, false, numMesoSub,
						mesoSubArray, subvolArray, tCur, spec.NUM_REGIONS,
						spec.NUM_MOL_TYPES, regionArray);
					heapMesoUpdate(numMesoSub, mesoSubArray, heap_subvolID, 0UL,
						heap_childID, b_heap_childValid); // ALWAYS head node
						
				} else
				{ // ID of reaction is beyond valid range. Error					
					fprintf(stderr,"ERROR: Current mesoscopic event %u in subvolume %" PRIu32 " is invalid.\nSubvolume is number %" PRIu32 " in the mesoscopic list\n", curRxn, curSub, curMeso);
					exit(EXIT_FAILURE);
				}
				// Update timer structure array
				timerArray[MESO_TIMER_ID].nextTime =
					mesoSubArray[heap_subvolID[0]].t_rxn;
				tMeso = mesoSubArray[heap_subvolID[0]].t_rxn;
			}
			
			// Update timer heap
			heapTimerUpdate(NUM_TIMERS, timerArray, heapTimer, 0,
				heapTimerChildID, b_heapTimerChildValid);
			
		}
		
		// Write realization observations to output file
		printOneTextRealization(out, spec, curRepeat, observationArray,
			numPassiveRecord, passiveRecordID, numActiveRecord, activeRecordID,
			actorCommonArray, actorActiveArray, actorPassiveArray,
			maxActiveBits, maxPassiveObs);
		
		if((curRepeat+1) % updateFreq == 0U)
		{
			fracComplete = (double) (curRepeat+1)/spec.NUM_REPEAT;
			printf("Simulation %.1f%% complete (%u of %u repeats). Est. time left: %.f sec.\n",
				fracComplete*100,curRepeat+1,spec.NUM_REPEAT,
				(double) (clock() - startTime)*(1/fracComplete - 1)/CLOCKS_PER_SEC);
		}
	}
	time(&timer);
	timeInfo = localtime(&timer);
	strftime(timeBuffer, 26, "%Y-%m-%d %H:%M:%S", timeInfo);
	printf("Ending simulation at %s.\n", timeBuffer);
	endTime = clock();
	runTime = (double) (endTime-startTime)/CLOCKS_PER_SEC;
	printf("Simulation ran in %f seconds\n", runTime);
	
	//
	// STEP 5: Save Output Summary
	//
	
	printf("Writing simulation summary file ...\n");
	
	// Print end time and info used to help Matlab importing
	printTextEnd(outSummary, numActiveRecord, numPassiveRecord, actorCommonArray,
		actorActiveArray, actorPassiveArray,
		passiveRecordID, activeRecordID, maxActiveBits, maxPassiveObs, runTime);
		
	//
	// STEP 6: Free Memory
	//
	printf("Memory cleanup ...\n");
	
	if (fclose(out) != 0)
		fprintf(stderr,"ERROR: Could not close output file \"%s.txt\".\n",spec.OUTPUT_NAME);
	if (fclose(outSummary) != 0)
		fprintf(stderr,"ERROR: Could not close output summary file \"%s_summary.txt\".\n", spec.OUTPUT_NAME);
	
	for(curActor = 0; curActor < numPassiveRecord; curActor++)
	{
		if(!isListEmptyObs(&observationArray[curActor]))
		{
			emptyListObs(&observationArray[curActor]);			
		}
	}
	for(i = 0; i < spec.NUM_REGIONS; i++)
	{
		for(j = 0; j < spec.NUM_MOL_TYPES; j++)
		{
			emptyListMol(&microMolList[i][j]);
			emptyListMol3DRecent(&microMolListRecent[i][j]);
		}
	}
	deleteActor(spec.NUM_ACTORS, actorCommonArray, regionArray,
		NUM_ACTORS_ACTIVE, actorActiveArray, NUM_ACTORS_PASSIVE, actorPassiveArray,
		passiveRecordID, activeRecordID);
	heapMesoDelete(numSub, heap_subvolID, heap_childID, b_heap_childValid);
	heapTimerDelete(heapTimer, heapTimerChildID, b_heapTimerChildValid);
	deleteTimerHeapArray(timerArray);
	deleteMesoSubArray(numMesoSub, mesoSubArray, subvolArray, spec.NUM_MOL_TYPES, spec.NUM_REGIONS);
	deleteSubvolArray(numSub, subvolArray, spec.NUM_MOL_TYPES, spec.NUM_REGIONS, regionArray);
	delete_boundary_region_(spec.NUM_REGIONS,
		spec.NUM_MOL_TYPES, regionArray);
	
	deleteConfig(spec);
	
	printf("Done!\n\n");
	return 0;
}
