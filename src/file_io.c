/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * file_io.c - interface with JSON configuration files
 *
 * Last revised for AcCoRD v1.4.2 (2020-02-12)
 *
 * Revision history:
 *
 * Revision v1.4.2 (2020-02-12)
 * - updated file output file directory creation to accommodate Mac OS
 *
 * Revision v1.4 (2018-08-06)
 * - added a priori monte carlo (APMC) absorption algorithm as a new surface
 * reaction type. Includes settings for how to define the a priori absorption
 * probability calculation and whether/how to apply a threshold to turn it off
 *
 * Revision v1.3 (2018-07-31)
 * - changed checks for the "Random Molecule Release Times?" "Slot Interval" parameters
 * for active actors, since whether they are needed depends on the values of other parameters.
 *
 * Revision v1.1 (2016-12-24)
 * - added flow parameters (based on either no flow or uniform flow).
 * Global settings apply to any subset of molecule types but can be modified for any
 * individual molecule type in any region.
 * - corrected default subvolume base size if is is not defined correctly
 * - corrected detection of missing surface reaction probability for surface reactions
 *
 * Revision v1.0 (2016-10-31)
 * - added warnings for defining passive actor parameters for active actors or
 * vice versa.
 * - added "BURST" modulation, which does not rely on bits and always releases molecules
 * in every action interval. Also is able to release multiple types of molecules.
 * - added simpler methods for defining region anchor coordinates and the number of
 * subvolumes along each dimension of rectangular regions.
 * - added local diffusion coefficients as region parameters. Any region can
 * over-ride the default diffusion coefficients defined for all molecules
 * - added specifying diffusion coefficient to use for surface reaction transition
 * probabilities. By default, the reactant's default diffusion coefficient is used
 * for absorption and membrane reactions, and the first product's default diffusion
 * coefficient is used for desorption reactions. Warning appears if coefficient is
 * defined for reaction types that cannot use it.
 *
 * Revision v0.7.0.1 (public beta, 2016-08-30)
 * - added measurement of simulation runtime to be written to simulation output
 * - fixed bug in writing index of active actor in simulation summary
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - added option for user to select small subvolume or big subvolume assumption
 * at hybrid interfaces. Only needed if there is at least 1 microscopic and 1
 * mesoscopic region defined.
 * - added option for user to choose maximum distance beyond which a micro to mesoscopic
 * substep transition will not be considered (either initial or final point should be beyond
 * this distance)
 * - added binding and unbinding radii as bimolecular reaction parameters (needed in
 * microscopic regime)
 * - made output of active actor data sequence a user option
 * - added option for user to define a constant active actor bit sequence
 * - added point active actors defined by their 3D coordinate.
 * - added check for positive (and not just non-negative) radii for spherical regions and
 * actors (i.e., 0 is not a valid radius)
 * - added warnings for unnecessary active actor parameters depending on values
 * of other active actor parameters
 *
 * Revision v0.5.1 (2016-05-06)
 * - added bReleaseProduct to chemical reaction. Applies to surface reactions
 * - added chemical reaction properties to define coupled reversible reactions
 * - added chemical reaction properties to configure absorbing, desorbing, and membrane
 * reactions, including how transition probabilities are calculated and how
 * desorbed molecules are placed. Membrane reactions must have no product molecules
 * specified (product is always the same as the reactant)
 * - shortened string used to indicate how actor location is defined. This avoids
 * a MATLAB warning when loaded simulation configuration files
 *
 * Revision v0.5 (2016-04-15)
 * - added ability to define location of actor by a list of regions
 * - modified check on number of subvolumes along each dimension of a rectangular region
 * - added type and surfaceType properties to region. Default values are REGION_NORMAL and
 * NO_SURFACE, respectively
 * - added bSurface and surfRxnType properties to chemical reaction. Default values are false and RXN_NORMAL, respectively.
 * - added stringWrite function to nest some calls to stringAllocate, strlen,
 * and strcpy
 * - removed NUM_DIM parameter from simulation spec
 * - removed upper limit on number of molecule types
 *
 * Revision v0.4.1
 * - added search for configuration file. First checks current directory, then
 * config subdirectory, and then config sibling directory
 * - added search for results folder in current directory and then as subdirectory
 * of parent. If it cannot be found, results folder is created in current directory
 * and output is placed there.
 * - improved use and format of error messages
 *
 * Revision v0.4
 * - modified use of unsigned long and unsigned long long to uint32_t and uint64_t
 * - added options to accommodate spherical regions and actors
 * - added restriction of chemical reactions to specific regions
 * - renamed region structure and array of structures to accommodate naming
 * of new boundary structure
 *
 * Revision v0.3.1.1
 * - corrected warning message when loading a configuration file to display the correct
 *   keys to quit or continue
 *
 * Revision v0.3.1
 * - check all JSON members for existence and correct format before reading them in
 * - header added
*/

#include "file_io.h"

// Load configuration file
void loadConfig(const char * CONFIG_NAME,
	uint32_t customSEED,
	struct simSpec3D * curSpec)
{
	FILE * configFile;
	long fileLength;
	cJSON * configJSON;
	cJSON *simControl, *environment, *regionSpec,
		*curObj, *curObjInner, *curArray,
		*actorSpec, *actorShape, *actorModScheme,
		*diffCoef, *flowVec, *flowLocal,
		*chemSpec, *rxnSpec;
	int arrayLen;
	char * tempString;
	int curArrayItem;
	unsigned short curMolType, i;
	unsigned short numFlowExcept, curFlowExcept;
	char * configContent;
	char * configNameFull;
	unsigned int dirLength, nameLength;
	bool bWarn = false; // Was there at least one warning?
	bool bWarnOverride; // Did the user request to override pause for warnings?
	int numWarn = 0; // Number of warnings found in configuration file.
	int ch;
	size_t temp; // Garbage variable for discarded file content length
	int minSubDim = 0; // Minimum # of subvolumes along each dimension for a rectangular region
	bool bNeedReleaseType; // For surface reaction, does user need to specify what happens to products?
	bool bHasMicro = false; // There is at least one microscopic region
	bool bHasMeso = false; // There is at least one mesoscopic region
	
	// Construct full name of configuration file
	nameLength = strlen(CONFIG_NAME);
	configNameFull = malloc(strlen("../config/") + nameLength + 1);
	if(configNameFull == NULL)
	{
		fprintf(stderr,"ERROR: Memory could not be allocated to build the configuration file name.\n");
		exit(EXIT_FAILURE);
	}
	
	// Open configuration file.
	// First check current directory, then check "config" folder,
	// then check "config" folder in parent directory
	configFile = fopen(CONFIG_NAME, "r");
	if(configFile == NULL)
	{
		strcpy(configNameFull,"config/");
		strcat(configNameFull,CONFIG_NAME);
		configFile = fopen( configNameFull, "r");
		if(configFile == NULL)
		{
			strcpy(configNameFull,"../config/");
			strcat(configNameFull,CONFIG_NAME);
			configFile = fopen( configNameFull, "r");
			if(configFile == NULL)
			{
				fprintf(stderr,"ERROR: Configuration file \"%s\" not found.\n",CONFIG_NAME);
				fprintf(stderr,"AcCoRD searches 1) in current directory, 2) in \"config\" subdirectory, and then 3) in \"..\\config\\\" directory.\n");
				exit(EXIT_FAILURE);
			}
		}
	} else
	{
		strcpy(configNameFull,CONFIG_NAME);
	}
	printf("Successfully opened configuration file at \"%s\".\n", configNameFull);
	
	// Read in contents of configuration file
	fseek(configFile, 0, SEEK_END);
	fileLength = ftell(configFile);
	fseek(configFile,0,SEEK_SET);
	
	configContent = malloc(fileLength + 1);
	if(configContent == NULL)
	{
		fprintf(stderr,"ERROR: Memory could not be allocated to store the configuration information.\n");
		exit(EXIT_FAILURE);
	}
	temp = fread(configContent,1,fileLength,configFile);
	fclose(configFile);
	
	// Convert file contents into a JSON object
	configJSON = cJSON_Parse(configContent);
	if (!configJSON)
	{
		fprintf(stderr, "ERROR: Invalid configuration file formatting in area of: [%s]\n",cJSON_GetErrorPtr());
		fprintf(stderr, "Could not convert file contents into a valid JSON object.\n");
		fprintf(stderr, "Please see AcCoRD documentation on how to write a configuration file.\n");
		exit(EXIT_FAILURE);
	}
	
	// Check Existence of Primary Structure Objects
	if(!cJSON_bItemValid(configJSON,"Simulation Control", cJSON_Object))
	{
		fprintf(stderr,"ERROR: Configuration file is missing \"Simulation Control\" object.\n");
		exit(EXIT_FAILURE);
	}
	if(!cJSON_bItemValid(configJSON,"Environment", cJSON_Object))
	{
		fprintf(stderr,"ERROR: Configuration file is missing \"Environment\" object.\n");
		exit(EXIT_FAILURE);
	}
	if(!cJSON_bItemValid(configJSON,"Chemical Properties", cJSON_Object))
	{
		fprintf(stderr,"ERROR: Configuration file is missing \"Chemical Properties\" object.\n");
		exit(EXIT_FAILURE);
	}
	
	// Check for warning override
	if(!cJSON_bItemValid(configJSON,"Warning Override", cJSON_True))
	{
		printf("WARNING %d: Configuration file is missing \"Warning Override\" boolean. Simulation will require user confirmation to execute.\n", numWarn++);
		bWarnOverride = false;
		bWarn = true;
	} else
	{
		bWarnOverride = cJSON_GetObjectItem(configJSON,"Warning Override")->valueint;
		if(bWarnOverride)
		{
			printf("NOTE: Warning override enabled. Simulation will execute automatically, even if warnings appear in the configuration file.\n");
		} else
		{
			printf("NOTE: Warning override disabled. Simulation will require user confirmation to execute if warnings appear in the configuration file.\n");
		}
	}
	
	//
	// Transfer JSON content to Simulation Structure	
	//
	
	// Load Simulation Control Object
	simControl = cJSON_GetObjectItem(configJSON,"Simulation Control");
	if(customSEED > 0)
	{ // User specified a seed when simulation was called
		curSpec->SEED = customSEED;
	} else
	{ // Use seed listed in the configuration, if it exists
		if(cJSON_bItemValid(simControl,"Random Number Seed", cJSON_Number))
		{
			curSpec->SEED =
				cJSON_GetObjectItem(simControl,"Random Number Seed")->valueint;
		} else
		{
			bWarn = true;
			printf("WARNING %d: \"Random Number Seed\" not defined and no custom seed specified. Assigning default value of \"0\".\n", numWarn++);
			curSpec->SEED = 0;
		}		
	}
	
	if(!cJSON_bItemValid(configJSON,"Output Filename", cJSON_String) ||
		strlen(cJSON_GetObjectItem(configJSON,"Output Filename")->valuestring) < 1)
	{ // Config file does not list a valid Output Filename
		bWarn = true;
		printf("WARNING %d: \"Output Filename\" not defined or has length zero. Assigning default value of \"test\".\n", numWarn++);
		curSpec->OUTPUT_NAME = stringAllocate(strlen("test") + 15);
		sprintf(curSpec->OUTPUT_NAME, "%s_SEED%d", "test", curSpec->SEED);
	} else{
		curSpec->OUTPUT_NAME = stringAllocate(strlen(cJSON_GetObjectItem(configJSON,
			"Output Filename")->valuestring) + 15);
		sprintf(curSpec->OUTPUT_NAME, "%s_SEED%d", cJSON_GetObjectItem(configJSON,
			"Output Filename")->valuestring, curSpec->SEED);
	}
	
	if(!cJSON_bItemValid(simControl,"Number of Repeats", cJSON_Number) ||
		cJSON_GetObjectItem(simControl,"Number of Repeats")->valueint < 0)
	{ // Config file does not list a valid Number of Repeats
		bWarn = true;
		printf("WARNING %d: \"Number of Repeats\" not defined or has invalid value. Assigning default value of \"1\" realization.\n", numWarn++);
		curSpec->NUM_REPEAT = 1;
	} else{
		curSpec->NUM_REPEAT =
			cJSON_GetObjectItem(simControl,"Number of Repeats")->valueint;
	}
		
	if(!cJSON_bItemValid(simControl,"Final Simulation Time", cJSON_Number) ||
		cJSON_GetObjectItem(simControl,"Final Simulation Time")->valuedouble < 0)
	{ // Config file does not list a valid Final Simulation Time
		bWarn = true;
		printf("WARNING %d: \"Final Simulation Time\" not defined or has invalid value. Assigning default value of \"0\" seconds.\n", numWarn++);
		curSpec->TIME_FINAL = 0.;
	} else{
		curSpec->TIME_FINAL =
			cJSON_GetObjectItem(simControl,"Final Simulation Time")->valuedouble;
	}
		
	if(!cJSON_bItemValid(simControl,"Global Microscopic Time Step", cJSON_Number) ||
		cJSON_GetObjectItem(simControl,"Global Microscopic Time Step")->valuedouble < 0)
	{ // Config file does not list a valid Global Microscopic Time Step
		bWarn = true;
		printf("WARNING %d: \"Global Microscopic Time Step\" not defined or has invalid value. Assigning default value of \"0\" seconds.\n", numWarn++);
		curSpec->DT_MICRO = 0.;
	} else{
		curSpec->DT_MICRO =
			cJSON_GetObjectItem(simControl,"Global Microscopic Time Step")->valuedouble;
	}
	
	if(!cJSON_bItemValid(simControl,"Max Number of Progress Updates", cJSON_Number) ||
		cJSON_GetObjectItem(simControl,"Max Number of Progress Updates")->valueint < 0)
	{ // Config file does not list a valid Number of Repeats
		bWarn = true;
		printf("WARNING %d: \"Max Number of Progress Updates\" not defined or has invalid value. Assigning default value of \"10\" updates.\n", numWarn++);
		curSpec->MAX_UPDATES = 10;
	} else{
		curSpec->MAX_UPDATES =
			cJSON_GetObjectItem(simControl,"Max Number of Progress Updates")->valueint;
	}
	
	// Load Chemical Properties Object
	chemSpec = cJSON_GetObjectItem(configJSON,"Chemical Properties");
		
	if(!cJSON_bItemValid(chemSpec,"Number of Molecule Types", cJSON_Number) ||
		cJSON_GetObjectItem(chemSpec,"Number of Molecule Types")->valueint < 1)
	{ // Config file does not list a valid Number of Molecule Types
		bWarn = true;
		printf("WARNING %d: \"Number of Molecule Types\" not defined or has invalid value. Assigning default value of \"1\" type.\n", numWarn++);
		curSpec->NUM_MOL_TYPES = 1;
	} else{
		curSpec->NUM_MOL_TYPES =
			cJSON_GetObjectItem(chemSpec,"Number of Molecule Types")->valueint;
	}
	
	curSpec->DIFF_COEF = malloc(curSpec->NUM_MOL_TYPES * sizeof(double));
	if(curSpec->DIFF_COEF == NULL)
	{
		fprintf(stderr,"ERROR: Memory could not be allocated to store diffusion coefficients\n");
		exit(EXIT_FAILURE);
	}
	
	if(!cJSON_bItemValid(chemSpec,"Diffusion Coefficients", cJSON_Array) ||
		cJSON_GetArraySize(cJSON_GetObjectItem(chemSpec,"Diffusion Coefficients")) != curSpec->NUM_MOL_TYPES)
	{ // Config file does not list a valid Diffusion Coefficients array
		bWarn = true;
		printf("WARNING %d: \"Diffusion Coefficients\" not defined or not of correct length. Assigning default value of \"0\" to each molecule type.\n", numWarn++);
		for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
		{
			curSpec->DIFF_COEF[curMolType] = 0.;
		}
	} else{
		diffCoef = cJSON_GetObjectItem(chemSpec, "Diffusion Coefficients");
		for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES;
			curMolType++)
		{
			if(!cJSON_bArrayItemValid(diffCoef,curMolType, cJSON_Number) ||
				cJSON_GetArrayItem(diffCoef,curMolType)->valuedouble < 0)
			{
				bWarn = true;
				printf("WARNING %d: \"Diffusion Coefficients\" item %d not defined or has an invalid value. Assigning default value of \"0\".\n", numWarn++, curMolType);
				curSpec->DIFF_COEF[curMolType] = 0;
			} else
			{
				curSpec->DIFF_COEF[curMolType] =
					cJSON_GetArrayItem(diffCoef, curMolType)->valuedouble;
			}
		}
	}
	
	// Load global flow parameters
	curSpec->B_GLOBAL_MOL_FLOW = malloc(curSpec->NUM_MOL_TYPES * sizeof(bool));
	if(curSpec->B_GLOBAL_MOL_FLOW == NULL)
	{
		fprintf(stderr,"ERROR: Memory could not be allocated to store boolean vector indicating which molecules can flow\n");
		exit(EXIT_FAILURE);
	}
	curSpec->GLOBAL_FLOW_VECTOR = NULL;
	if(!cJSON_bItemValid(chemSpec,"Global Flow Type", cJSON_String))
	{ // Config file does not list a valid Global Flow Type
		bWarn = true;
		printf("WARNING %d: \"Global Flow Type\" not defined or has invalid value. Assigning default value \"None\".\n", numWarn++);
		curSpec->GLOBAL_FLOW_TYPE = FLOW_NONE;
		
		arrayLen = 0;
	} else{
		tempString = stringWrite(cJSON_GetObjectItem(chemSpec,"Global Flow Type")->valuestring);
		if(strcmp(tempString,"None") == 0)
		{ // No molecule flow
			curSpec->GLOBAL_FLOW_TYPE = FLOW_NONE;
			arrayLen = 0;
		} else if(strcmp(tempString,"Uniform") == 0)
		{ // Steady uniform flow
			curSpec->GLOBAL_FLOW_TYPE = FLOW_UNIFORM;
			arrayLen = 3;
		} else
		{
			bWarn = true;
			printf("WARNING %d: \"Global Flow Type\" has invalid value. Assigning default value \"None\".\n", numWarn++);
			curSpec->GLOBAL_FLOW_TYPE = FLOW_NONE;
			arrayLen = 0;
		}
		free(tempString);
	}
	
	// Load flow parameters
	if(arrayLen > 0)
	{
		curSpec->GLOBAL_FLOW_VECTOR = malloc(arrayLen * sizeof(double));
		if(curSpec->GLOBAL_FLOW_VECTOR == NULL)
		{
			fprintf(stderr,"ERROR: Memory could not be allocated to store global flow vector.\n");
			exit(EXIT_FAILURE);
		}
		
		if(!cJSON_bItemValid(chemSpec,"Global Flow Vector", cJSON_Array) ||
			cJSON_GetArraySize(cJSON_GetObjectItem(chemSpec,"Global Flow Vector"))
			!= arrayLen)
		{
			bWarn = true;
			printf("WARNING %d: \"Global Flow Vector\" was not defined or has the wrong length. Assigning default value of \"0\" to each element.\n", numWarn++);
			for(i = 0; i < arrayLen; i++)
			{
				curSpec->GLOBAL_FLOW_VECTOR[i] = 0.;						
			}
		} else
		{
			flowVec = cJSON_GetObjectItem(chemSpec, "Global Flow Vector");
			for(i = 0; i < arrayLen; i++)
			{
				if(!cJSON_bArrayItemValid(flowVec,i, cJSON_Number))
				{
					bWarn = true;
					printf("WARNING %d: \"Global Flow Vector\" item %d not defined or has an invalid value. Assigning default value of \"0\".\n", numWarn++, i);
					curSpec->GLOBAL_FLOW_VECTOR[i] = 0.;
				} else
				{
					curSpec->GLOBAL_FLOW_VECTOR[i] =
						cJSON_GetArrayItem(flowVec, i)->valuedouble;
				}
			}
		}
		
		// Determine which molecule types the global flow applies to
		if(cJSON_bItemValid(chemSpec,"Does Molecule Type Flow?", cJSON_Array))
		{
			if(cJSON_GetArraySize(cJSON_GetObjectItem(chemSpec,"Does Molecule Type Flow?")) != curSpec->NUM_MOL_TYPES)
			{
				bWarn = true;
				printf("WARNING %d: \"Does Molecule Type Flow?\" does not have the correct length. Applying global flow to each type of molecule.\n", numWarn++);
				for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
				{
					curSpec->B_GLOBAL_MOL_FLOW[curMolType] = true;
				}
			} else
			{
				curObj = cJSON_GetObjectItem(chemSpec,"Does Molecule Type Flow?");
				for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
				{
					if(!cJSON_bArrayItemValid(curObj,curMolType, cJSON_True))
					{
						bWarn = true;
						printf("WARNING %d: Molecule type %d does not have a valid global \"Does Molecule Type Flow?\" value. Setting to true.\n", numWarn++, curMolType);
						curSpec->B_GLOBAL_MOL_FLOW[curMolType] = true;
					} else
					{
						curSpec->B_GLOBAL_MOL_FLOW[curMolType] = 
							cJSON_GetArrayItem(curObj, curMolType)->valueint;
					}
				}
			}
		} else
		{
			// By default, global flow applies to every type of molecule
			printf("Note: \"Does Molecule Type Flow?\" array was not defined to specify which molecule types experience global flow. Applying to all.\n");
			for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
			{
				curSpec->B_GLOBAL_MOL_FLOW[curMolType] = true;
			}
		}
	} else
	{
		for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
		{
			curSpec->B_GLOBAL_MOL_FLOW[curMolType] = false;
		}
		
		if(cJSON_bItemValid(chemSpec,"Global Flow Vector", cJSON_Array) ||
			cJSON_bItemValid(chemSpec,"Does Molecule Type Flow?", cJSON_Array))
		{
			bWarn = true;
			printf("WARNING %d: Global flow parameters were defined but are not needed because there is no global flow. Ignoring.\n", numWarn++);
		}		
	}			
	
	if(!cJSON_bItemValid(chemSpec,"Chemical Reaction Specification", cJSON_Array))
	{
		bWarn = true;
		printf("WARNING %d: Configuration file is missing \"Chemical Reaction Specification\" array. Assuming that no chemical reactions are possible.", numWarn++);
		curSpec->MAX_RXNS = 0;
	} else
	{
		rxnSpec = cJSON_GetObjectItem(chemSpec, "Chemical Reaction Specification");
		curSpec->MAX_RXNS =
			cJSON_GetArraySize(rxnSpec);
		curSpec->chem_rxn =
			malloc(curSpec->MAX_RXNS *
			sizeof(struct chem_rxn_struct));
		if(curSpec->chem_rxn == NULL)
		{
			fprintf(stderr,"ERROR: Memory could not be allocated to store chemical reaction details\n");
			exit(EXIT_FAILURE);
		}
		for(curArrayItem = 0;
			curArrayItem < curSpec->MAX_RXNS; curArrayItem++)
		{
			curSpec->chem_rxn[curArrayItem].reactants =
				malloc(curSpec->NUM_MOL_TYPES * sizeof(uint32_t));
			curSpec->chem_rxn[curArrayItem].products =
				malloc(curSpec->NUM_MOL_TYPES * sizeof(uint32_t));
			curSpec->chem_rxn[curArrayItem].bReleaseProduct =
				malloc(curSpec->NUM_MOL_TYPES * sizeof(uint32_t));
			if(curSpec->chem_rxn[curArrayItem].reactants == NULL
				|| curSpec->chem_rxn[curArrayItem].products == NULL
				|| curSpec->chem_rxn[curArrayItem].bReleaseProduct == NULL)
			{
				fprintf(stderr,"ERROR: Memory could not be allocated to store chemical reaction reactants and products\n");
				exit(EXIT_FAILURE);
			}
			
			if(!cJSON_bArrayItemValid(rxnSpec, curArrayItem, cJSON_Object))
			{
				bWarn = true;
				printf("WARNING %d: \"Chemical Reaction Specification\" item %d is not a JSON object. Creating empty reaction.\n", numWarn++, curArrayItem);
				curSpec->chem_rxn[curArrayItem].k = 0.;
				curSpec->chem_rxn[curArrayItem].bSurface = false;
				curSpec->chem_rxn[curArrayItem].bEverywhere = false;
				curSpec->chem_rxn[curArrayItem].numRegionExceptions = 0;
				curSpec->chem_rxn[curArrayItem].regionExceptionLabel = NULL;
				for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES;
				curMolType++)
				{
					curSpec->chem_rxn[curArrayItem].reactants[curMolType] = 0;
					curSpec->chem_rxn[curArrayItem].products[curMolType] = 0;
					curSpec->chem_rxn[curArrayItem].bReleaseProduct[curMolType] = false;
				}
			} else
			{
				
				curObj = cJSON_GetArrayItem(rxnSpec, curArrayItem);
				if(!cJSON_bItemValid(curObj,"Reaction Rate", cJSON_Number) ||
					!cJSON_bItemValid(curObj,"Reactants", cJSON_Array) ||
					!cJSON_bItemValid(curObj,"Products", cJSON_Array) ||
					cJSON_GetArraySize(cJSON_GetObjectItem(curObj,"Reactants")) != curSpec->NUM_MOL_TYPES ||
					cJSON_GetArraySize(cJSON_GetObjectItem(curObj,"Products")) != curSpec->NUM_MOL_TYPES ||
					cJSON_GetObjectItem(curObj,"Reaction Rate")->valuedouble < 0.)
				{
					bWarn = true;
					printf("WARNING %d: \"Chemical Reaction Specification\" item %d has missing parameters, an invalid reaction rate, or an incorrect number of molecule types. Creating empty reaction.\n", numWarn++, curArrayItem);
					curSpec->chem_rxn[curArrayItem].k = 0.;
					curSpec->chem_rxn[curArrayItem].bSurface = false;
					curSpec->chem_rxn[curArrayItem].bReleaseProduct = NULL;
					curSpec->chem_rxn[curArrayItem].bEverywhere = false;
					curSpec->chem_rxn[curArrayItem].numRegionExceptions = 0;
					curSpec->chem_rxn[curArrayItem].regionExceptionLabel = NULL;
					for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES;
					curMolType++)
					{
						curSpec->chem_rxn[curArrayItem].reactants[curMolType] = 0;
						curSpec->chem_rxn[curArrayItem].products[curMolType] = 0;
						curSpec->chem_rxn[curArrayItem].bReleaseProduct[curMolType] = false;
					}
				} else
				{		
					// Reaction label
					if(!cJSON_bItemValid(curObj,"Label", cJSON_String))
					{ // Reaction does not have a defined Label
						bWarn = true;
						printf("WARNING %d: Chemical reaction %d does not have a defined \"Label\". Assigning empy string.\n", numWarn++, curArrayItem);
						curSpec->chem_rxn[curArrayItem].label = '\0';
					} else{
						curSpec->chem_rxn[curArrayItem].label =
							stringWrite(cJSON_GetObjectItem(curObj,"Label")->valuestring);
					}
					
					// Is reaction reversible?
					if(!cJSON_bItemValid(curObj,"Is Reaction Reversible?", cJSON_True))
					{ // Reaction does not have a valid Is Reaction Reversible?
						bWarn = true;
						printf("WARNING %d: Chemical reaction %d does not have a valid \"Is Reaction Reversible?\". Assigning default value \"false\".\n", numWarn++, curArrayItem);
						curSpec->chem_rxn[curArrayItem].bReversible = false;
					} else
					{
						curSpec->chem_rxn[curArrayItem].bReversible = 
							cJSON_GetObjectItem(curObj, "Is Reaction Reversible?")->valueint;
					
						if(curSpec->chem_rxn[curArrayItem].bReversible)
						{ // Reaction is reversible. Store name of reverse reaction
							if(!cJSON_bItemValid(curObj,"Reverse Reaction Label", cJSON_String))
							{ // Reaction does not have a defined Reverse Reaction Label
								bWarn = true;
								printf("WARNING %d: Chemical reaction %d is reversible but does not have a defined \"Reverse Reaction Label\". Assigning empy string.\n", numWarn++, curArrayItem);
								curSpec->chem_rxn[curArrayItem].labelCoupled = '\0';
							} else{
								curSpec->chem_rxn[curArrayItem].labelCoupled =
									stringWrite(cJSON_GetObjectItem(curObj,"Reverse Reaction Label")->valuestring);
							}
						} else if(cJSON_bItemValid(curObj,"Reverse Reaction Label", cJSON_String))
						{
							bWarn = true;
							printf("WARNING %d: Chemical reaction %d does not need \"Reverse Reaction Label\" defined. Ignoring.\n", numWarn++, curArrayItem);
							curSpec->chem_rxn[curArrayItem].labelCoupled = '\0';
						}
					}
					
					
					if(!cJSON_bItemValid(curObj,"Surface Reaction?", cJSON_True))
					{ // Reaction does not have a valid Surface Reaction?
						bWarn = true;
						printf("WARNING %d: Chemical reaction %d does not have a valid \"Surface Reaction?\". Assigning default value \"false\".\n", numWarn++, curArrayItem);
						curSpec->chem_rxn[curArrayItem].bSurface = false;
						curSpec->chem_rxn[curArrayItem].surfRxnType = RXN_NORMAL;
					} else
					{
						curSpec->chem_rxn[curArrayItem].bSurface = 
							cJSON_GetObjectItem(curObj, "Surface Reaction?")->valueint;
					}
					
					if(curSpec->chem_rxn[curArrayItem].bSurface)
					{ // We have a surface reaction. Determine what type
						if(!cJSON_bItemValid(curObj,"Surface Reaction Type", cJSON_String))
						{ // Reaction does not have a defined Surface Reaction Type
							bWarn = true;
							printf("WARNING %d: Chemical reaction %d does not have a defined \"Surface Reaction Type\". Setting to default value \"Normal\".\n", numWarn++, curArrayItem);
							tempString =
								stringWrite("Normal");
						} else{
							tempString =
								stringWrite(cJSON_GetObjectItem(curObj,
								"Surface Reaction Type")->valuestring);
						}
						
						if(strcmp(tempString,"Normal") == 0)
							curSpec->chem_rxn[curArrayItem].surfRxnType = RXN_NORMAL;
						else if(strcmp(tempString,"Absorbing") == 0)
							curSpec->chem_rxn[curArrayItem].surfRxnType = RXN_ABSORBING;
						else if(strcmp(tempString,"A Priori Absorbing") == 0)
							curSpec->chem_rxn[curArrayItem].surfRxnType = RXN_A_PRIORI_ABSORBING;
						else if(strcmp(tempString,"Desorbing") == 0)
							curSpec->chem_rxn[curArrayItem].surfRxnType = RXN_DESORBING;
						else if(strcmp(tempString,"Receptor Binding") == 0)
							curSpec->chem_rxn[curArrayItem].surfRxnType = RXN_RECEPTOR;
						else if(strcmp(tempString,"Membrane Inner") == 0)
							curSpec->chem_rxn[curArrayItem].surfRxnType = RXN_MEMBRANE_IN;
						else if(strcmp(tempString,"Membrane Outer") == 0)
							curSpec->chem_rxn[curArrayItem].surfRxnType = RXN_MEMBRANE_OUT;
						else
						{
							bWarn = true;
							printf("WARNING %d: Chemical reaction %d has an invalid \"Surface Reaction Type\". Setting to default value \"Normal\".\n", numWarn++, curArrayItem);
							curSpec->chem_rxn[curArrayItem].surfRxnType = RXN_NORMAL;
						}
						free(tempString);
						
						switch(curSpec->chem_rxn[curArrayItem].surfRxnType)
						{
							case RXN_ABSORBING:
							case RXN_A_PRIORI_ABSORBING:
							case RXN_DESORBING:
							case RXN_MEMBRANE_IN:
							case RXN_MEMBRANE_OUT:
								// These surface reactions have choices in how reaction
								// probability is calculated
								if(!cJSON_bItemValid(curObj,"Surface Transition Probability", cJSON_String))
								{ // Reaction does not have a defined Surface Transition Probability
									bWarn = true;
									printf("WARNING %d: Chemical reaction %d does not have a defined \"Surface Transition Probability\". Setting to default value \"Normal\".\n", numWarn++, curArrayItem);
									tempString =
										stringWrite("Normal");
								} else{
									tempString =
										stringWrite(cJSON_GetObjectItem(curObj,
										"Surface Transition Probability")->valuestring);
								}
								
								if(strcmp(tempString,"Normal") == 0)
									curSpec->chem_rxn[curArrayItem].rxnProbType = RXN_PROB_NORMAL;
								else if(strcmp(tempString,"Mixed") == 0)
									curSpec->chem_rxn[curArrayItem].rxnProbType = RXN_PROB_MIXED;
								else if(strcmp(tempString,"Steady State") == 0)
									curSpec->chem_rxn[curArrayItem].rxnProbType = RXN_PROB_STEADY_STATE;
								else if(strcmp(tempString,"A Priori Sphere") == 0)
									curSpec->chem_rxn[curArrayItem].rxnProbType = RXN_PROB_A_PRIORI_SPHERE;
								else if(strcmp(tempString,"A Priori Plane") == 0)
									curSpec->chem_rxn[curArrayItem].rxnProbType = RXN_PROB_A_PRIORI_INFINITE_PLANE;
								else
								{
									bWarn = true;
									printf("WARNING %d: Chemical reaction %d has an invalid \"Surface Transition Probability\". Setting to default value \"Normal\".\n", numWarn++, curArrayItem);
									curSpec->chem_rxn[curArrayItem].rxnProbType = RXN_PROB_NORMAL;
								}
								free(tempString);
								break;
							default:
								if(cJSON_bItemValid(curObj,"Surface Transition Probability", cJSON_String))
								{ // Reaction does not have a defined Surface Transition Probability
									bWarn = true;
									printf("WARNING %d: Chemical reaction %d does not need a \"Surface Transition Probability\". Ignoring.\n", numWarn++, curArrayItem);
								}
						}
						
						switch(curSpec->chem_rxn[curArrayItem].surfRxnType)
						{
							case RXN_ABSORBING:
							case RXN_A_PRIORI_ABSORBING:						
							case RXN_DESORBING:
								// Reaction needs bReleaseProduct array defined
								bNeedReleaseType = false; // By default, we don't need to define release type
								if(!cJSON_bItemValid(curObj,"Products Released?", cJSON_Array) ||
									cJSON_GetArraySize(cJSON_GetObjectItem(curObj,
									"Products Released?")) != curSpec->NUM_MOL_TYPES)
								{ // Reaction does not have a valid Products Released? array
									bWarn = true;
									printf("WARNING %d: Chemical reaction %d does not have a valid \"Products Released?\" array. Setting all values to \"false\".\n", numWarn++, curArrayItem);	
									for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES;
									curMolType++)
									{
										curSpec->chem_rxn[curArrayItem].bReleaseProduct[curMolType] = false;
									}
									curSpec->chem_rxn[curArrayItem].releaseType =
										PROD_PLACEMENT_LEAVE;
								} else{
									curObjInner = cJSON_GetObjectItem(curObj,"Products Released?");
									for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES;
									curMolType++)
									{
										if(!cJSON_bArrayItemValid(curObjInner,curMolType, cJSON_True))
										{
											bWarn = true;
											printf("WARNING %d: Molecule type %d does not have a valid \"Products Released?\" value for chemical reaction %d. Setting to default value \"false\".\n", numWarn++, curMolType, curArrayItem);
											curSpec->chem_rxn[curArrayItem].bReleaseProduct[curMolType] = false;
										} else
										{
											curSpec->chem_rxn[curArrayItem].bReleaseProduct[curMolType] =
												cJSON_GetArrayItem(curObjInner, curMolType)->valueint;
											if(cJSON_GetArrayItem(curObjInner, curMolType)->valueint)
											{
												bNeedReleaseType = true;
											}
										}
									}
									
									// How are released molecules placed?
									if(bNeedReleaseType)
									{
										if(!cJSON_bItemValid(curObj,"Release Placement Type", cJSON_String))
										{ // Reaction does not have a defined Release Placement Type
											bWarn = true;
											printf("WARNING %d: Chemical reaction %d does not have a defined \"Release Placement Type\". Setting to default value \"Leave\".\n", numWarn++, curArrayItem);
											tempString =
												stringWrite("Leave");
										} else{
											tempString =
												stringWrite(cJSON_GetObjectItem(curObj,
												"Release Placement Type")->valuestring);
										}
										
										if(strcmp(tempString,"Leave") == 0)
											curSpec->chem_rxn[curArrayItem].releaseType =
												PROD_PLACEMENT_LEAVE;
										else if(strcmp(tempString,"Full Diffusion") == 0)
											curSpec->chem_rxn[curArrayItem].releaseType = PROD_PLACEMENT_FORCE;
										else if(strcmp(tempString,"Steady State Diffusion") == 0)
											curSpec->chem_rxn[curArrayItem].releaseType = PROD_PLACEMENT_STEADY_STATE;
										else
										{
											bWarn = true;
											printf("WARNING %d: Chemical reaction %d has an invalid \"Release Placement Type\". Setting to default value \"Leave\".\n", numWarn++, curArrayItem);
											curSpec->chem_rxn[curArrayItem].releaseType = PROD_PLACEMENT_LEAVE;
										}
										free(tempString);
									} else
									{
										if(cJSON_bItemValid(curObj,"Release Placement Type", cJSON_String))
										{ // Reaction did not need defined Release Placement Type
											bWarn = true;
											printf("WARNING %d: Chemical reaction %d has \"Release Placement Type\" defined but there are no products released. Ignoring.\n", numWarn++, curArrayItem);
										}
									}
								}
								break;
							default:
								// Reaction should not have bReleaseProduct array defined
								if(cJSON_bItemValid(curObj,"Products Released?", cJSON_Array))
								{
									bWarn = true;
									printf("WARNING %d: Reaction %d is a surface reaction but does not need \"Products Released?\" defined. Ignoring.\n", numWarn++, curArrayItem);
								}
						}
						
						// Check for the additional parameters needed for A Priori absorption
						if(curSpec->chem_rxn[curArrayItem].surfRxnType == RXN_A_PRIORI_ABSORBING)
						{
							if(!cJSON_bItemValid(curObj,"Surface Reaction Threshold Type", cJSON_String))
							{ // Reaction does not have a valid Surface Reaction Threshold Type
								bWarn = true;
								printf("WARNING %d: Chemical reaction %d does not have a valid \"Surface Reaction Threshold Type\". Assigning default value \"None\".\n", numWarn++, curArrayItem);
								tempString = stringWrite("None");								
							} else
							{
								tempString =
									stringWrite(cJSON_GetObjectItem(curObj,
									"Surface Reaction Threshold Type")->valuestring);
							}
							
							if(strcmp(tempString,"None") == 0)
							{
								curSpec->chem_rxn[curArrayItem].bRxnThreshold = false;
								curSpec->chem_rxn[curArrayItem].rxnThresholdType = RXN_THRESHOLD_DISTANCE;
								curSpec->chem_rxn[curArrayItem].rxnThreshold = Inf;
							} else
							{
								curSpec->chem_rxn[curArrayItem].bRxnThreshold = true;
								
								if(!cJSON_bItemValid(curObj,"Surface Reaction Threshold Value", cJSON_Number)
									|| cJSON_GetObjectItem(curObj,"Surface Reaction Threshold Value")->valuedouble < 0.)
								{ // Reaction does not have a valid Surface Reaction Threshold Value
									bWarn = true;
									printf("WARNING %d: Chemical reaction %d does not have a \"Surface Reaction Threshold Value\". Disabling threshold.\n", numWarn++, curArrayItem);
									curSpec->chem_rxn[curArrayItem].bRxnThreshold = false;
									curSpec->chem_rxn[curArrayItem].rxnThreshold = Inf;
								}
								else
								{
									curSpec->chem_rxn[curArrayItem].rxnThreshold =
										cJSON_GetObjectItem(curObj,"Surface Reaction Threshold Value")->valuedouble;
										
									if(strcmp(tempString,"Distance") == 0)
										curSpec->chem_rxn[curArrayItem].rxnThresholdType = RXN_THRESHOLD_DISTANCE;
									else if(strcmp(tempString,"Probability") == 0)
										curSpec->chem_rxn[curArrayItem].rxnThresholdType = RXN_THRESHOLD_PROB;
									else if(strcmp(tempString,"Current Region") == 0)
										curSpec->chem_rxn[curArrayItem].rxnThresholdType = RXN_THRESHOLD_REGION;
									else
									{
										bWarn = true;
										printf("WARNING %d: Chemical reaction %d has an invalid \"Surface Reaction Threshold Type\". Setting to default value \"None\".\n", numWarn++, curArrayItem);
										curSpec->chem_rxn[curArrayItem].bRxnThreshold = false;
										curSpec->chem_rxn[curArrayItem].rxnThresholdType = RXN_THRESHOLD_DISTANCE;
										curSpec->chem_rxn[curArrayItem].rxnThreshold = Inf;
									}	
								}								
							}
							free(tempString);								
								
							// Reaction must have an APMC-type transition probability defined
							if (curSpec->chem_rxn[curArrayItem].rxnProbType != RXN_PROB_A_PRIORI_SPHERE
								&& curSpec->chem_rxn[curArrayItem].rxnProbType != RXN_PROB_A_PRIORI_INFINITE_PLANE)
							{
								bWarn = true;
								printf("WARNING %d: Chemical reaction %d is an A Priori type reaction but has a non-A Priori type \"Surface Transition Probability\". Setting to default value \"A Priori Plane\".\n", numWarn++, curArrayItem);
								curSpec->chem_rxn[curArrayItem].rxnProbType = RXN_PROB_A_PRIORI_INFINITE_PLANE;
							}
						} else
						{
							// Reaction should not have the surface reaction threshold parameters defined
							if(cJSON_bItemValid(curObj,"Surface Reaction Threshold Type", cJSON_String)
								|| cJSON_bItemValid(curObj,"Surface Reaction Threshold Value", cJSON_String))
							{
								bWarn = true;
								printf("WARNING %d: Reaction %d is a surface reaction that does not need either \"Surface Reaction Threshold Type\" or \"Surface Reaction Threshold Value\" defined. Ignoring.\n", numWarn++, curArrayItem);
							}
							
							// Reaction should not have APMC-type transition probabilities
							if (curSpec->chem_rxn[curArrayItem].rxnProbType == RXN_PROB_A_PRIORI_SPHERE
								|| curSpec->chem_rxn[curArrayItem].rxnProbType == RXN_PROB_A_PRIORI_INFINITE_PLANE)
							{
								bWarn = true;
								printf("WARNING %d: Chemical reaction %d is not an A Priori type reaction but has an A Priori type \"Surface Transition Probability\". Setting to default value \"Normal\".\n", numWarn++, curArrayItem);
								curSpec->chem_rxn[curArrayItem].rxnProbType = RXN_PROB_NORMAL;
							}
						}
						
					} else
					{
						curSpec->chem_rxn[curArrayItem].surfRxnType = RXN_NORMAL;
						for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES;
						curMolType++)
						{
							curSpec->chem_rxn[curArrayItem].bReleaseProduct[curMolType] = false;
						}
						// Check for existence of unnecessary parameters and display warnings if they are defined
						if(cJSON_bItemValid(curObj,"Surface Reaction Type", cJSON_String))
						{
							bWarn = true;
							printf("WARNING %d: Reaction %d is not a surface reaction and so does not need \"Surface Reaction Type\" defined. Ignoring.\n", numWarn++, curArrayItem);
						}
						if(cJSON_bItemValid(curObj,"Products Released?", cJSON_Array))
						{
							bWarn = true;
							printf("WARNING %d: Reaction %d is not a surface reaction and so does not need \"Products Released?\" defined. Ignoring.\n", numWarn++, curArrayItem);
						}						
					}
					
					if(!cJSON_bItemValid(curObj,"Default Everywhere?", cJSON_True))
					{ // Reaction does not have a valid Default Everywhere?
						bWarn = true;
						printf("WARNING %d: Chemical reaction %d does not have a valid \"Default Everywhere?\". Assigning default value \"true\".\n", numWarn++, curArrayItem);
						curSpec->chem_rxn[curArrayItem].bEverywhere = true;
					} else
					{
						curSpec->chem_rxn[curArrayItem].bEverywhere = 
							cJSON_GetObjectItem(curObj, "Default Everywhere?")->valueint;
					}
						
					// Record exceptions to the default reaction location
					if(!cJSON_bItemValid(curObj,"Exception Regions", cJSON_Array))
					{ // Chemical reaction does not have an Exception Regions array
						bWarn = true;
						printf("WARNING %d: Chemical reaction %d has a missing or invalid \"Default Everywhere?\". Assigning default value of \"0\" exceptions.\n", numWarn++, curArrayItem);
						curSpec->chem_rxn[curArrayItem].numRegionExceptions = 0;
					} else
					{
						// Read number of exceptions
						curSpec->chem_rxn[curArrayItem].numRegionExceptions =
							cJSON_GetArraySize(cJSON_GetObjectItem(curObj,"Exception Regions"));
						
						curSpec->chem_rxn[curArrayItem].regionExceptionLabel =
							malloc(curSpec->chem_rxn[curArrayItem].numRegionExceptions * sizeof(char *));
						if(curSpec->chem_rxn[curArrayItem].regionExceptionLabel == NULL)
						{
							fprintf(stderr,"ERROR: Memory could not be allocated to store chemical reaction region exceptions\n");
							exit(EXIT_FAILURE);
						}
						
						// Read in names of exception regions							
						curObjInner = cJSON_GetObjectItem(curObj,"Exception Regions");
						for(i = 0; i < curSpec->chem_rxn[curArrayItem].numRegionExceptions; i++)
						{
							if(!cJSON_bArrayItemValid(curObjInner,i, cJSON_String))
							{ // Exception region is not a valid string. Ignore
								bWarn = true;
								printf("WARNING %d: Chemical reaction %d exception region %d is not a valid string. Assigning empty string.\n", numWarn++, curArrayItem, i);
								curSpec->chem_rxn[curArrayItem].regionExceptionLabel[i] = '\0';
							} else
							{
								curSpec->chem_rxn[curArrayItem].regionExceptionLabel[i] =
									stringWrite(cJSON_GetArrayItem(curObjInner,i)->valuestring);
							}					
						}
					}
					
					curSpec->chem_rxn[curArrayItem].k = 
						cJSON_GetObjectItem(curObj, "Reaction Rate")->valuedouble;
					
					// Check for binding and unbinding radii
					if(cJSON_bItemValid(curObj,"Binding Radius", cJSON_Number))
					{
						if(cJSON_GetObjectItem(curObj, "Binding Radius")->valuedouble >= 0.)
						{
							curSpec->chem_rxn[curArrayItem].rBind = 
								cJSON_GetObjectItem(curObj, "Binding Radius")->valuedouble;
						} else
						{
							bWarn = true;
							printf("WARNING %d: Reaction %d has an invalid binding radius. Setting to default value of \"0\".\n", numWarn++, curArrayItem);
							curSpec->chem_rxn[curArrayItem].rBind = 0.;
						}
					} else
						curSpec->chem_rxn[curArrayItem].rBind = 0.;
					
					if(cJSON_bItemValid(curObj,"Unbinding Radius", cJSON_Number))
					{						
						if(cJSON_GetObjectItem(curObj, "Unbinding Radius")->valuedouble >= 0.)
						{
							curSpec->chem_rxn[curArrayItem].rUnbind = 
								cJSON_GetObjectItem(curObj, "Unbinding Radius")->valuedouble;
						} else
						{
							bWarn = true;
							printf("WARNING %d: Reaction %d has an invalid unbinding radius. Setting to default value of \"0\".\n", numWarn++, curArrayItem);
							curSpec->chem_rxn[curArrayItem].rUnbind = 0.;
						}
					} else
						curSpec->chem_rxn[curArrayItem].rUnbind = 0.;
					
					// Read array of reactant stoichiometry
					curObjInner = cJSON_GetObjectItem(curObj,"Reactants");
					for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES;
					curMolType++)
					{
						if(!cJSON_bArrayItemValid(curObjInner,curMolType, cJSON_Number) ||
							cJSON_GetArrayItem(curObjInner,curMolType)->valueint < 0)
						{
							bWarn = true;
							printf("WARNING %d: Molecule type %d has an incorrect number of reactants in reaction %d. Setting to default value of \"0\".\n", numWarn++, curMolType, curArrayItem);
							curSpec->chem_rxn[curArrayItem].reactants[curMolType] = 0;
						} else
						{
							curSpec->chem_rxn[curArrayItem].reactants[curMolType] =
								cJSON_GetArrayItem(curObjInner, curMolType)->valueint;
						}
					}
					
					// Read array of product stoichiometry
					curObjInner = cJSON_GetObjectItem(curObj,"Products");
					for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES;
					curMolType++)
					{
						if(!cJSON_bArrayItemValid(curObjInner,curMolType, cJSON_Number) ||
							cJSON_GetArrayItem(curObjInner,curMolType)->valueint < 0)
						{
							bWarn = true;
							printf("WARNING %d: Molecule type %d has an incorrect number of products in reaction %d. Setting to default value of \"0\".\n", numWarn++, curMolType, curArrayItem);
							curSpec->chem_rxn[curArrayItem].products[curMolType] = 0;
						} else if((curSpec->chem_rxn[curArrayItem].surfRxnType == RXN_MEMBRANE_IN
							|| curSpec->chem_rxn[curArrayItem].surfRxnType == RXN_MEMBRANE_OUT)
							&& cJSON_GetArrayItem(curObjInner,curMolType)->valueint > 0)
						{
							bWarn = true;
							printf("WARNING %d: Molecule type %d in reaction %d has a non-zero number of product molecules defined. Membrane reactions should have no product molecules defined. Setting to default value of \"0\".\n", numWarn++, curMolType, curArrayItem);
							curSpec->chem_rxn[curArrayItem].products[curMolType] = 0;
						} else
						{
							curSpec->chem_rxn[curArrayItem].products[curMolType] =
								cJSON_GetArrayItem(curObjInner, curMolType)->valueint;
						}
					}
					
					// Determine diffusion coefficient to use for surface reactions
					if(curSpec->chem_rxn[curArrayItem].bSurface)
					{
						if(cJSON_bItemValid(curObj,"Surface Reaction Diffusion Coefficient", cJSON_Number) &&
						cJSON_GetObjectItem(curObj,"Surface Reaction Diffusion Coefficient")->valuedouble > 0.)
						{ // Read specified diffusion coefficient
							curSpec->chem_rxn[curArrayItem].diffusion =
								cJSON_GetObjectItem(curObj,"Surface Reaction Diffusion Coefficient")->valuedouble;
						} else
						{ // Determine diffusion coefficient from reactant or product
							switch(curSpec->chem_rxn[curArrayItem].surfRxnType)
							{
								case RXN_DESORBING:
									// Find product
									for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES;
										curMolType++)
									{
										if(curSpec->chem_rxn[curArrayItem].products[curMolType] > 0)
										{
											curSpec->chem_rxn[curArrayItem].diffusion =
												curSpec->DIFF_COEF[curMolType];
											break;
										}
									}
									break;
								default:
									// Find reactant
									for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES;
										curMolType++)
									{
										if(curSpec->chem_rxn[curArrayItem].reactants[curMolType] > 0)
										{
											curSpec->chem_rxn[curArrayItem].diffusion =
												curSpec->DIFF_COEF[curMolType];
											break;
										}
									}
									break;
							}
						}
					} else if(cJSON_bItemValid(curObj,"Surface Reaction Diffusion Coefficient", cJSON_Number))
					{
						bWarn = true;
						printf("WARNING %d: Reaction %d is not a surface reaction and so does not need \"Surface Reaction Diffusion Coefficient\" defined. Ignoring.\n", numWarn++, curArrayItem);
					}
				}
			}
		}
	}
	
	// Load Environment Object
	environment = cJSON_GetObjectItem(configJSON,"Environment");
	if(!cJSON_bItemValid(environment,"Region Specification", cJSON_Array) ||
		cJSON_GetArraySize(cJSON_GetObjectItem(environment,"Region Specification")) < 1)
	{
		fprintf(stderr,"ERROR: Configuration file is missing \"Region Specification\" array in \"Environment\" object or it has a length less than 1.\n");
		exit(EXIT_FAILURE);
	}
	if(!cJSON_bItemValid(environment,"Actor Specification", cJSON_Array) ||
		cJSON_GetArraySize(cJSON_GetObjectItem(environment,"Actor Specification")) < 1)
	{
		fprintf(stderr,"ERROR: Configuration file is missing \"Actor Specification\" array in \"Environment\" object or it has a length less than 1.\n");
		exit(EXIT_FAILURE);
	}
	
	if(!cJSON_bItemValid(environment,"Subvolume Base Size", cJSON_Number) ||
		cJSON_GetObjectItem(environment,"Subvolume Base Size")->valuedouble <= 0)
	{
		bWarn = true;
		printf("WARNING %d: \"Subvolume Base Size\" not defined or is invalid. Setting to default value of \"1\".\n", numWarn++);
		curSpec->SUBVOL_BASE_SIZE = 1.;
	} else
	{
		curSpec->SUBVOL_BASE_SIZE =
			cJSON_GetObjectItem(environment,"Subvolume Base Size")->valuedouble;
	}	
	
	regionSpec = cJSON_GetObjectItem(environment,"Region Specification");
	actorSpec = cJSON_GetObjectItem(environment,"Actor Specification");
	curSpec->NUM_REGIONS =
		cJSON_GetArraySize(regionSpec);
	curSpec->NUM_ACTORS =
		cJSON_GetArraySize(actorSpec);

	curSpec->subvol_spec =
		malloc(curSpec->NUM_REGIONS *
		sizeof(struct spec_region3D));
	curSpec->actorSpec =
		malloc(curSpec->NUM_ACTORS *
		sizeof(struct actorStructSpec3D));
	if(curSpec->subvol_spec == NULL
		|| curSpec->actorSpec == NULL)
	{
		fprintf(stderr,"ERROR: Memory could not be allocated to load region or actor details\n");
		exit(EXIT_FAILURE);
	}
	for(curArrayItem = 0;
		curArrayItem < curSpec->NUM_REGIONS; curArrayItem++)
	{
		if(!cJSON_bArrayItemValid(regionSpec, curArrayItem, cJSON_Object))
		{
			fprintf(stderr, "ERROR: Region %d is not described by a JSON object.\n", curArrayItem);
			exit(EXIT_FAILURE);
		}
		
		curObj = cJSON_GetArrayItem(regionSpec, curArrayItem);
		
		// Region label
		if(!cJSON_bItemValid(curObj,"Label", cJSON_String))
		{ // Region does not have a defined Label
			bWarn = true;
			printf("WARNING %d: Region %d does not have a defined \"Label\". Assigning empy string.\n", numWarn++, curArrayItem);
			curSpec->subvol_spec[curArrayItem].label = '\0';
		} else{
			curSpec->subvol_spec[curArrayItem].label =
				stringWrite(cJSON_GetObjectItem(curObj,"Label")->valuestring);
		}
		
		// Local diffusion coefficients
		curSpec->subvol_spec[curArrayItem].bLocalDiffusion = false;
		if(cJSON_bItemValid(curObj,"Local Diffusion Coefficients", cJSON_Array))
		{ // Region has custom local diffusion coefficients
			if(cJSON_GetArraySize(cJSON_GetObjectItem(curObj,"Local Diffusion Coefficients")) != curSpec->NUM_MOL_TYPES)
			{
				bWarn = true;
				printf("WARNING %d: Region %d has defined \"Local Diffusion Coefficients\" but not of the correct length. Ignoring.\n", numWarn++, curArrayItem);
			} else{
				curSpec->subvol_spec[curArrayItem].bLocalDiffusion = true;
				diffCoef = cJSON_GetObjectItem(curObj, "Local Diffusion Coefficients");
				curSpec->subvol_spec[curArrayItem].diffusion =
					malloc(curSpec->NUM_MOL_TYPES * sizeof(double));
				for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES;
				curMolType++)
				{
					if(!cJSON_bArrayItemValid(diffCoef,curMolType, cJSON_Number) ||
						cJSON_GetArrayItem(diffCoef,curMolType)->valuedouble < 0)
					{
						bWarn = true;
						printf("WARNING %d: \"Diffusion Coefficients\" item %d not defined or has an invalid value. Assigning default value of \"0\".\n", numWarn++, curMolType);
						curSpec->subvol_spec[curArrayItem].diffusion[curMolType] = 0;
					} else
					{
						curSpec->subvol_spec[curArrayItem].diffusion[curMolType] =
							cJSON_GetArrayItem(diffCoef, curMolType)->valuedouble;
					}
				}
			}
		}
		
		// Local flow parameters
		curSpec->subvol_spec[curArrayItem].bFlow =
			malloc(curSpec->NUM_MOL_TYPES * sizeof(bool));
		curSpec->subvol_spec[curArrayItem].bFlowLocal =
			malloc(curSpec->NUM_MOL_TYPES * sizeof(bool));
		curSpec->subvol_spec[curArrayItem].flowType =
			malloc(curSpec->NUM_MOL_TYPES * sizeof(unsigned short));
		curSpec->subvol_spec[curArrayItem].flowVector =
			malloc(curSpec->NUM_MOL_TYPES * sizeof(double *));
		if(curSpec->subvol_spec[curArrayItem].bFlow == NULL ||
			curSpec->subvol_spec[curArrayItem].bFlowLocal == NULL ||
			curSpec->subvol_spec[curArrayItem].flowType == NULL ||
			curSpec->subvol_spec[curArrayItem].flowVector == NULL)
		{
			fprintf(stderr,"ERROR: Memory could not be allocated to store region flow parameters.\n");
			exit(EXIT_FAILURE);
		}
		for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
		{ // Assign default (global) flow parameters
			curSpec->subvol_spec[curArrayItem].bFlowLocal[curMolType] = false;
			curSpec->subvol_spec[curArrayItem].flowVector[curMolType] = NULL;
			if(curSpec->B_GLOBAL_MOL_FLOW[curMolType])
			{
				curSpec->subvol_spec[curArrayItem].bFlow[curMolType] = true;
				curSpec->subvol_spec[curArrayItem].flowType[curMolType] =
					curSpec->GLOBAL_FLOW_TYPE;
			} else
			{
				curSpec->subvol_spec[curArrayItem].bFlow[curMolType] = false;
				curSpec->subvol_spec[curArrayItem].flowType[curMolType] =
					FLOW_NONE;
			}
		}
		// Check for region exceptions to global flow parameters
		if(cJSON_bItemValid(curObj,"Local Flow", cJSON_Array))
		{ // Region has custom local flow parameters
			flowLocal = cJSON_GetObjectItem(curObj, "Local Flow");
			numFlowExcept = cJSON_GetArraySize(flowLocal);
			for(curFlowExcept = 0; curFlowExcept < numFlowExcept; curFlowExcept++)
			{
				// Read current flow exception
				if(!cJSON_bArrayItemValid(flowLocal, curFlowExcept, cJSON_Object))
				{
					bWarn = true;
					printf("WARNING %d: \"Local Flow\" item %d in region %d is not a JSON object. Ignoring.\n", numWarn++, curFlowExcept, curArrayItem);
					continue;
				}
				
				// Check for valid exception
				curObjInner = cJSON_GetArrayItem(flowLocal, curFlowExcept);
				if(!cJSON_bItemValid(curObjInner,"Is Molecule Type Affected?", cJSON_Array) ||
					cJSON_GetArraySize(cJSON_GetObjectItem(curObjInner,"Is Molecule Type Affected?")) != curSpec->NUM_MOL_TYPES ||
					!cJSON_bItemValid(curObjInner,"Flow Type", cJSON_String))
				{
					bWarn = true;
					printf("WARNING %d: \"Local Flow\" item %d in region %d has missing parameters. Ignoring\n", numWarn++, curFlowExcept, curArrayItem);
					continue;
				}
				
				curArray = cJSON_GetObjectItem(curObjInner,"Is Molecule Type Affected?");
				// Apply exception to valid molecule types where applicable
				for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
				{
					if(!cJSON_bArrayItemValid(curArray,curMolType, cJSON_True))
					{
						bWarn = true;
						printf("WARNING %d: Molecule type %d does not have a valid \"Is Molecule Type Affected?\" value for \"Local Flow Parameters\" item %d in region %d. Ignoring.\n", numWarn++, curMolType, curFlowExcept, curArrayItem);
						continue;
					}
					if(cJSON_GetArrayItem(curArray, curMolType)->valueint)
					{ // Exception is supposed to apply to this type of molecule
						if(curSpec->subvol_spec[curArrayItem].bFlowLocal[curMolType])
						{ // Exception was already defined for this type of molecule!
							bWarn = true;
							printf("WARNING %d: Molecule type %d already had local flow parameters defined for region %d, but they are being defined again by \"Local Flow Parameters\" item %d. Ignoring.\n", numWarn++, curMolType, curArrayItem, curFlowExcept);
							continue;
						}
						curSpec->subvol_spec[curArrayItem].bFlowLocal[curMolType] = true;
						// Load flow type for current molecule in this region
						if(!cJSON_bItemValid(curObjInner,"Flow Type", cJSON_String))
						{ // Config file does not list a valid Flow Type
							bWarn = true;
							printf("WARNING %d: \"Flow Type\" in local flow exception %d of region %d not defined or has invalid value. Assigning default value \"None\" to molecule type %d in this region.\n",
								numWarn++, curFlowExcept, curArrayItem, curMolType);
							curSpec->subvol_spec[curArrayItem].flowType[curMolType] = FLOW_NONE;
							continue;
						}
						tempString =
							stringWrite(cJSON_GetObjectItem(curObjInner,"Flow Type")->valuestring);
						if(strcmp(tempString,"None") == 0)
						{ // No molecule flow
							curSpec->subvol_spec[curArrayItem].flowType[curMolType] = FLOW_NONE;
							curSpec->subvol_spec[curArrayItem].bFlow[curMolType] = false;
							arrayLen = 0;
						} else if(strcmp(tempString,"Uniform") == 0)
						{ // Steady uniform flow
							curSpec->subvol_spec[curArrayItem].flowType[curMolType] = FLOW_UNIFORM;
							curSpec->subvol_spec[curArrayItem].bFlow[curMolType] = true;
							arrayLen = 3;
						} else
						{
							bWarn = true;
							printf("WARNING %d: \"Flow Type\" in local flow exception %d of region %d not defined or has invalid value. Assigning default value \"None\" to molecule type %d in this region.\n",
								numWarn++, curFlowExcept, curArrayItem, curMolType);
							curSpec->subvol_spec[curArrayItem].flowType[curMolType] = FLOW_NONE;
							curSpec->subvol_spec[curArrayItem].bFlow[curMolType] = false;
							arrayLen = 0;
						}
						free(tempString);
						
						// Load flow vector for this type of molecule
						if(arrayLen > 0)
						{
							curSpec->subvol_spec[curArrayItem].flowVector[curMolType] =
								malloc(arrayLen * sizeof(double));
							if(curSpec->subvol_spec[curArrayItem].flowVector[curMolType] == NULL)
							{
								fprintf(stderr,"ERROR: Memory could not be allocated to store local flow vector.\n");
								exit(EXIT_FAILURE);
							}
							
							if(!cJSON_bItemValid(curObjInner,"Flow Vector", cJSON_Array) ||
								cJSON_GetArraySize(cJSON_GetObjectItem(curObjInner,"Flow Vector"))
								!= arrayLen)
							{
								bWarn = true;
								printf("WARNING %d: \"Flow Vector\" in local flow exception %d of region %d was not defined or has the wrong length. Assigning default value of \"0\" to each element.\n", numWarn++, curFlowExcept, curArrayItem);
								for(i = 0; i < arrayLen; i++)
								{
									curSpec->subvol_spec[curArrayItem].flowVector[curMolType][i] = 0.;						
								}
								continue;
							}
											
							flowVec = cJSON_GetObjectItem(curObjInner, "Flow Vector");
							for(i = 0; i < arrayLen; i++)
							{
								if(!cJSON_bArrayItemValid(flowVec,i, cJSON_Number))
								{
									bWarn = true;
									printf("WARNING %d: \"Flow Vector\" item %d in local flow exception %d of region %d not defined or has an invalid value. Assigning default value of \"0\".\n", numWarn++, i, curFlowExcept, curArrayItem);
									curSpec->subvol_spec[curArrayItem].flowVector[curMolType][i] = 0.;
								} else
								{
									curSpec->subvol_spec[curArrayItem].flowVector[curMolType][i] =
										cJSON_GetArrayItem(flowVec, i)->valuedouble;
								}
							}
						} else if(cJSON_bItemValid(curObjInner,"Flow Vector", cJSON_Array))
						{
							bWarn = true;
							printf("WARNING %d: \"Flow Vector\" was defined in local flow exception %d of region %d but is not needed. Ignoring.\n", numWarn++, curFlowExcept, curArrayItem);
						}
					}
				}
			}
		}
		// Apply global flow settings to remaining molecules
		switch(curSpec->GLOBAL_FLOW_TYPE)
		{
			case FLOW_NONE:
				arrayLen = 0;
				break;
			case FLOW_UNIFORM:
				arrayLen = 3;
				break;
		}
		if(arrayLen > 0)
		{			
			for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
			{
				if(curSpec->subvol_spec[curArrayItem].bFlowLocal[curMolType])
					continue; // An exception was already recorded for this type
				curSpec->subvol_spec[curArrayItem].flowVector[curMolType] =
					malloc(arrayLen * sizeof(double));
				if(curSpec->subvol_spec[curArrayItem].flowVector[curMolType] == NULL)
				{
					fprintf(stderr,"ERROR: Memory could not be allocated to store local flow vector for molecule type % in region %d.\n", curMolType, curArrayItem);
					exit(EXIT_FAILURE);
				}
				for(i = 0; i < arrayLen; i++)
				{
					curSpec->subvol_spec[curArrayItem].flowVector[curMolType][i] =
						curSpec->GLOBAL_FLOW_VECTOR[i];
				}
			}
		}
		
		// Region Parent
		if(!cJSON_bItemValid(curObj,"Parent Label", cJSON_String))
		{ // Region does not have a defined Parent Label
			bWarn = true;
			printf("WARNING %d: Region %d does not have a defined \"Parent Label\". Assigning empy string.\n", numWarn++, curArrayItem);
			curSpec->subvol_spec[curArrayItem].parent = '\0';
		} else{
			curSpec->subvol_spec[curArrayItem].parent =
				stringWrite(cJSON_GetObjectItem(curObj,"Parent Label")->valuestring);
		}
		
		// Region Shape
		if(!cJSON_bItemValid(curObj,"Shape", cJSON_String))
		{ // Region does not have a defined Shape
			bWarn = true;
			printf("WARNING %d: Region %d does not have a defined \"Shape\". Setting to default value \"Rectangular Box\".\n", numWarn++, curArrayItem);
			tempString =
				stringWrite("Rectangular Box");
		} else{
			tempString =
				stringWrite(cJSON_GetObjectItem(curObj,
				"Shape")->valuestring);
		}
		
		if(strcmp(tempString,"Rectangle") == 0)
			curSpec->subvol_spec[curArrayItem].shape = RECTANGLE;
		else if(strcmp(tempString,"Circle") == 0)
			curSpec->subvol_spec[curArrayItem].shape = CIRCLE;
		else if(strcmp(tempString,"Rectangular Box") == 0)
			curSpec->subvol_spec[curArrayItem].shape = RECTANGULAR_BOX;
		else if(strcmp(tempString,"Sphere") == 0)
			curSpec->subvol_spec[curArrayItem].shape = SPHERE;
		else
		{
			bWarn = true;
			printf("WARNING %d: Region %d has an invalid \"Shape\". Setting to default value \"Rectangular Box\".\n", numWarn++, curArrayItem);
			curSpec->subvol_spec[curArrayItem].shape = RECTANGULAR_BOX;
		}
		free(tempString);
		
		// Region Type
		if(!cJSON_bItemValid(curObj,"Type", cJSON_String))
		{ // Region does not have a defined Type
			bWarn = true;
			printf("WARNING %d: Region %d does not have a defined \"Type\". Setting to default value \"Normal\".\n", numWarn++, curArrayItem);
			tempString =
				stringWrite("Normal");
		} else{
			tempString =
				stringWrite(cJSON_GetObjectItem(curObj,
				"Type")->valuestring);
		}
		
		if(strcmp(tempString,"Normal") == 0)
		{
			curSpec->subvol_spec[curArrayItem].type = REGION_NORMAL;
			curSpec->subvol_spec[curArrayItem].surfaceType = NO_SURFACE;
			if(cJSON_bItemValid(curObj,"Surface Type", cJSON_String))
			{
				bWarn = true;
				printf("WARNING %d: Region %d does not need \"Surface Type\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
		}
		else if(strcmp(tempString,"3D Surface") == 0)
			curSpec->subvol_spec[curArrayItem].type = REGION_SURFACE_3D;
		else if(strcmp(tempString,"2D Surface") == 0)
		{
			if(curSpec->subvol_spec[curArrayItem].shape == RECTANGLE)
				curSpec->subvol_spec[curArrayItem].type = REGION_SURFACE_2D;
			else
			{
				bWarn = true;
				printf("WARNING %d: Region %d is a 3D shape but was classified as a 2D surface. Changing to \"3D surface\".\n", numWarn++, curArrayItem);
				curSpec->subvol_spec[curArrayItem].type = REGION_SURFACE_3D;
			}
		}
		else
		{
			bWarn = true;
			printf("WARNING %d: Region %d has an invalid \"Type\". Setting to default value \"Normal\".\n", numWarn++, curArrayItem);
			curSpec->subvol_spec[curArrayItem].type = REGION_NORMAL;
			curSpec->subvol_spec[curArrayItem].surfaceType = NO_SURFACE;
		}
		free(tempString);
		
		if(curSpec->subvol_spec[curArrayItem].type != REGION_NORMAL)
		{
			if(!cJSON_bItemValid(curObj,"Surface Type", cJSON_String))
			{
				bWarn = true;
				printf("WARNING %d: Region %d does not have a valid \"Surface Type\". Assigning default value \"Membrane\".\n", numWarn++, curArrayItem);
				tempString = stringWrite("Membrane");
				curSpec->subvol_spec[curArrayItem].surfaceType = SURFACE_MEMBRANE;
			} else
			{
				tempString = stringWrite(cJSON_GetObjectItem(curObj,
				"Surface Type")->valuestring);
			}
			if(strcmp(tempString,"Membrane") == 0)
				curSpec->subvol_spec[curArrayItem].surfaceType = SURFACE_MEMBRANE;
			else if(strcmp(tempString,"Inner") == 0)
				curSpec->subvol_spec[curArrayItem].surfaceType = SURFACE_INNER;
			else if(strcmp(tempString,"Outer") == 0)
				curSpec->subvol_spec[curArrayItem].surfaceType = SURFACE_OUTER;
			else
			{
				bWarn = true;
				printf("WARNING %d: Region %d has an invalid \"Surface Type\". Setting to default value \"Membrane\".\n", numWarn++, curArrayItem);
				curSpec->subvol_spec[curArrayItem].surfaceType = SURFACE_MEMBRANE;
			}
			free(tempString);
		}
		
		// Region Position
		if(!cJSON_bItemValid(curObj,"Anchor Coordinate", cJSON_Array) ||
			cJSON_GetArraySize(cJSON_GetObjectItem(curObj,"Anchor Coordinate")) != 3)
		{ // Region is not using ``newer'' anchor definition format
			printf("Region %d does not correctly use the new format for \"Anchor Coordinate\". Checking for old format.\n", curArrayItem);
			
			if(!cJSON_bItemValid(curObj,"Anchor X Coordinate", cJSON_Number))
			{ // Region does not have a valid Anchor X Coordinate
				bWarn = true;
				printf("WARNING %d: Region %d does not have a valid \"Anchor X Coordinate\". Assigning default value \"0\".\n", numWarn++, curArrayItem);
				curSpec->subvol_spec[curArrayItem].xAnch = 0;
			} else
			{
				curSpec->subvol_spec[curArrayItem].xAnch = 
					cJSON_GetObjectItem(curObj, "Anchor X Coordinate")->valuedouble;
			}
			
			if(!cJSON_bItemValid(curObj,"Anchor Y Coordinate", cJSON_Number))
			{ // Region does not have a valid Anchor Y Coordinate
				bWarn = true;
				printf("WARNING %d: Region %d does not have a valid \"Anchor Y Coordinate\". Assigning default value \"0\".\n", numWarn++, curArrayItem);
				curSpec->subvol_spec[curArrayItem].yAnch = 0;
			} else
			{
				curSpec->subvol_spec[curArrayItem].yAnch = 
					cJSON_GetObjectItem(curObj, "Anchor Y Coordinate")->valuedouble;
			}
			
			if(!cJSON_bItemValid(curObj,"Anchor Z Coordinate", cJSON_Number))
			{ // Region does not have a valid Anchor Z Coordinate
				bWarn = true;
				printf("WARNING %d: Region %d does not have a valid \"Anchor Z Coordinate\". Assigning default value \"0\".\n", numWarn++, curArrayItem);
				curSpec->subvol_spec[curArrayItem].zAnch = 0;
			} else
			{
				curSpec->subvol_spec[curArrayItem].zAnch = 
					cJSON_GetObjectItem(curObj, "Anchor Z Coordinate")->valuedouble;
			}
			
		} else
		{
			curObjInner = cJSON_GetObjectItem(curObj,"Anchor Coordinate");
			
			if(!cJSON_bArrayItemValid(curObjInner,0, cJSON_Number))
			{ // Region does not have a valid Anchor X Coordinate
				bWarn = true;
				printf("WARNING %d: Region %d does not have a valid X coordinate in \"Anchor Coordinate\" array. Assigning default value \"0\".\n", numWarn++, curArrayItem);
				curSpec->subvol_spec[curArrayItem].xAnch = 0;
			} else
			{
				curSpec->subvol_spec[curArrayItem].xAnch = 
					cJSON_GetArrayItem(curObjInner, 0)->valuedouble;
			}
			
			if(!cJSON_bArrayItemValid(curObjInner,1, cJSON_Number))
			{ // Region does not have a valid Anchor Y Coordinate
				bWarn = true;
				printf("WARNING %d: Region %d does not have a valid Y coordinate in \"Anchor Coordinate\" array. Assigning default value \"0\".\n", numWarn++, curArrayItem);
				curSpec->subvol_spec[curArrayItem].yAnch = 0;
			} else
			{
				curSpec->subvol_spec[curArrayItem].yAnch = 
					cJSON_GetArrayItem(curObjInner, 1)->valuedouble;
			}
			
			if(!cJSON_bArrayItemValid(curObjInner,2, cJSON_Number))
			{ // Region does not have a valid Anchor Z Coordinate
				bWarn = true;
				printf("WARNING %d: Region %d does not have a valid Z coordinate in \"Anchor Coordinate\" array. Assigning default value \"0\".\n", numWarn++, curArrayItem);
				curSpec->subvol_spec[curArrayItem].zAnch = 0;
			} else
			{
				curSpec->subvol_spec[curArrayItem].zAnch = 
					cJSON_GetArrayItem(curObjInner, 2)->valuedouble;
			}
		}
		
		
		
		if(cJSON_bItemValid(curObj,"Time Step", cJSON_Number))
		{
			bWarn = true;
			printf("WARNING %d: Region %d does not need \"Time Step\" defined. This will be implemented in a future version. Ignoring.\n", numWarn++, curArrayItem);
		}
		
		// Load remaining parameters depending on region shape
		if(curSpec->subvol_spec[curArrayItem].shape == RECTANGULAR_BOX ||
			curSpec->subvol_spec[curArrayItem].shape == RECTANGLE)
		{
			curSpec->subvol_spec[curArrayItem].radius = 0;
			// Check for existence of unnecessary parameters and display
			// warnings if they are defined.
			if(cJSON_bItemValid(curObj,"Radius", cJSON_Number))
			{
				bWarn = true;
				printf("WARNING %d: Region %d does not need \"Radius\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
		
			// Width of subvolumes in region (multiple of SUBVOL_BASE_SIZE)
			if(!cJSON_bItemValid(curObj,"Integer Subvolume Size", cJSON_Number) ||
				cJSON_GetObjectItem(curObj,"Integer Subvolume Size")->valueint < 1)
			{ // Region does not have a valid Integer Subvolume Size
				bWarn = true;
				printf("WARNING %d: Region %d does not have a valid \"Integer Subvolume Size\". Assigning default value \"1\".\n", numWarn++, curArrayItem);
				curSpec->subvol_spec[curArrayItem].sizeRect = 1;
			} else
			{
				curSpec->subvol_spec[curArrayItem].sizeRect = 
					cJSON_GetObjectItem(curObj, "Integer Subvolume Size")->valueint;
			}
			
			// Is region microscopic or mesoscopic?
			if(!cJSON_bItemValid(curObj,"Is Region Microscopic?", cJSON_True))
			{ // Region does not have a valid Is Region Microscopic?
				bWarn = true;
				printf("WARNING %d: Region %d does not have a valid \"Is Region Microscopic?\". Assigning default value \"false\".\n", numWarn++, curArrayItem);
				curSpec->subvol_spec[curArrayItem].bMicro = false;
			} else
			{
				curSpec->subvol_spec[curArrayItem].bMicro = 
					cJSON_GetObjectItem(curObj, "Is Region Microscopic?")->valueint;
			}
			
			// Track whether we need might need to consider hybrid interfaces
			if(curSpec->subvol_spec[curArrayItem].bMicro)
				bHasMicro = true;
			else
				bHasMeso = true;
			
			if(curSpec->subvol_spec[curArrayItem].shape == RECTANGLE)
				minSubDim = 0;
			else
				minSubDim = 1;
			
			// Number of subvolumes along each dimension
			if(!cJSON_bItemValid(curObj,"Number of Subvolumes Per Dimension", cJSON_Array) ||
				cJSON_GetArraySize(cJSON_GetObjectItem(curObj,"Number of Subvolumes Per Dimension")) != 3)
			{ // Region is not using ``newer'' subvolume definition format
				printf("Region %d does not correctly use the new format for \"Number of Subvolumes Per Dimension\". Checking for old format.\n", curArrayItem);
				
				if(!cJSON_bItemValid(curObj,"Number of Subvolumes Along X", cJSON_Number) ||
					cJSON_GetObjectItem(curObj,"Number of Subvolumes Along X")->valueint < minSubDim)
				{ // Region does not have a valid Number of Subvolumes Along X
					bWarn = true;
					printf("WARNING %d: Region %d does not have a valid \"Number of Subvolumes Along X\". Assigning default value \"1\".\n", numWarn++, curArrayItem);
					curSpec->subvol_spec[curArrayItem].numX = 1;
				} else
				{
					curSpec->subvol_spec[curArrayItem].numX = 
						cJSON_GetObjectItem(curObj, "Number of Subvolumes Along X")->valueint;
				}
				
				if(!cJSON_bItemValid(curObj,"Number of Subvolumes Along Y", cJSON_Number) ||
					cJSON_GetObjectItem(curObj,"Number of Subvolumes Along Y")->valueint < minSubDim)
				{ // Region does not have a valid Number of Subvolumes Along Y
					bWarn = true;
					printf("WARNING %d: Region %d does not have a valid \"Number of Subvolumes Along Y\". Assigning default value \"1\".\n", numWarn++, curArrayItem);
					curSpec->subvol_spec[curArrayItem].numY = 1;
				} else
				{
					curSpec->subvol_spec[curArrayItem].numY = 
						cJSON_GetObjectItem(curObj, "Number of Subvolumes Along Y")->valueint;
				}
				
				if(!cJSON_bItemValid(curObj,"Number of Subvolumes Along Z", cJSON_Number) ||
					cJSON_GetObjectItem(curObj,"Number of Subvolumes Along Z")->valueint < minSubDim)
				{ // Region does not have a valid Number of Subvolumes Along Z
					bWarn = true;
					printf("WARNING %d: Region %d does not have a valid \"Number of Subvolumes Along Z\". Assigning default value \"1\".\n", numWarn++, curArrayItem);
					curSpec->subvol_spec[curArrayItem].numZ = 1;
				} else
				{
					curSpec->subvol_spec[curArrayItem].numZ = 
						cJSON_GetObjectItem(curObj, "Number of Subvolumes Along Z")->valueint;
				}
				
			} else
			{
				curObjInner =
					cJSON_GetObjectItem(curObj,"Number of Subvolumes Per Dimension");
				
				if(!cJSON_bArrayItemValid(curObjInner, 0, cJSON_Number) ||
					cJSON_GetArrayItem(curObjInner, 0)->valueint < minSubDim)
				{ // Region does not have a valid Number of Subvolumes Along X
					bWarn = true;
					printf("WARNING %d: Region %d does not have a valid X-value in \"Number of Subvolumes Per Dimension\". Assigning default value \"1\".\n", numWarn++, curArrayItem);
					curSpec->subvol_spec[curArrayItem].numX = 1;
				} else
				{
					curSpec->subvol_spec[curArrayItem].numX = 
						cJSON_GetArrayItem(curObjInner, 0)->valueint;
				}
				
				if(!cJSON_bArrayItemValid(curObjInner, 1, cJSON_Number) ||
					cJSON_GetArrayItem(curObjInner, 1)->valueint < minSubDim)
				{ // Region does not have a valid Number of Subvolumes Along Y
					bWarn = true;
					printf("WARNING %d: Region %d does not have a valid Y-value in \"Number of Subvolumes Per Dimension\". Assigning default value \"1\".\n", numWarn++, curArrayItem);
					curSpec->subvol_spec[curArrayItem].numY = 1;
				} else
				{
					curSpec->subvol_spec[curArrayItem].numY = 
						cJSON_GetArrayItem(curObjInner, 1)->valueint;
				}
				
				if(!cJSON_bArrayItemValid(curObjInner, 2, cJSON_Number) ||
					cJSON_GetArrayItem(curObjInner, 2)->valueint < minSubDim)
				{ // Region does not have a valid Number of Subvolumes Along Z
					bWarn = true;
					printf("WARNING %d: Region %d does not have a valid Z-value in \"Number of Subvolumes Per Dimension\". Assigning default value \"1\".\n", numWarn++, curArrayItem);
					curSpec->subvol_spec[curArrayItem].numZ = 1;
				} else
				{
					curSpec->subvol_spec[curArrayItem].numZ = 
						cJSON_GetArrayItem(curObjInner, 2)->valueint;
				}
			}
			
			// Confirm that a rectangle region is actually 2D
			if(curSpec->subvol_spec[curArrayItem].shape == RECTANGLE)
			{
				if((curSpec->subvol_spec[curArrayItem].numX == 0
					&& (curSpec->subvol_spec[curArrayItem].numY < 1
					|| curSpec->subvol_spec[curArrayItem].numZ < 1))
					|| (curSpec->subvol_spec[curArrayItem].numY == 0
					&& (curSpec->subvol_spec[curArrayItem].numX < 1
					|| curSpec->subvol_spec[curArrayItem].numZ < 1))
					|| (curSpec->subvol_spec[curArrayItem].numZ == 0
					&& (curSpec->subvol_spec[curArrayItem].numY < 1
					|| curSpec->subvol_spec[curArrayItem].numX < 1))
					|| (curSpec->subvol_spec[curArrayItem].numX > 0
					&& curSpec->subvol_spec[curArrayItem].numY > 0
					&& curSpec->subvol_spec[curArrayItem].numZ > 0))
				{
					bWarn = true;
					printf("WARNING %d: Region %d is not properly defined as a Rectangle. Defining along XY plane with 1 subvolume along X and Y.\n", numWarn++, curArrayItem);
					curSpec->subvol_spec[curArrayItem].numX = 1;
					curSpec->subvol_spec[curArrayItem].numY = 1;
					curSpec->subvol_spec[curArrayItem].numZ = 0;
				}
			}
		} else // Region is round
		{
			curSpec->subvol_spec[curArrayItem].sizeRect = 0;
			curSpec->subvol_spec[curArrayItem].bMicro = true;
			bHasMicro = true;
			curSpec->subvol_spec[curArrayItem].numX = 1;
			curSpec->subvol_spec[curArrayItem].numY = 1;
			curSpec->subvol_spec[curArrayItem].numZ = 1;
			// Check for existence of unnecessary parameters and display
			// warnings if they are defined.
			if(cJSON_bItemValid(curObj,"Integer Subvolume Size", cJSON_Number))
			{
				bWarn = true;
				printf("WARNING %d: Region %d does not need \"Integer Subvolume Size\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Is Region Microscopic?", cJSON_True))
			{
				bWarn = true;
				printf("WARNING %d: Region %d does not need \"Is Region Microscopic?\" defined. This region must be microscopic. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Number of Subvolumes Along X", cJSON_Number))
			{
				bWarn = true;
				printf("WARNING %d: Region %d does not need \"Number of Subvolumes Along X\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Number of Subvolumes Along Y", cJSON_Number))
			{
				bWarn = true;
				printf("WARNING %d: Region %d does not need \"Number of Subvolumes Along Y\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Number of Subvolumes Along Z", cJSON_Number))
			{
				bWarn = true;
				printf("WARNING %d: Region %d does not need \"Number of Subvolumes Along Z\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Number of Subvolumes Per Dimension", cJSON_Array))
			{
				bWarn = true;
				printf("WARNING %d: Region %d does not need \"Number of Subvolumes Per Dimension.\n", numWarn++, curArrayItem);
			}
		
			// Region radius
			if(!cJSON_bItemValid(curObj,"Radius", cJSON_Number) ||
				cJSON_GetObjectItem(curObj,"Radius")->valuedouble <= 0.)
			{ // Region does not have a valid Radius
				bWarn = true;
				printf("WARNING %d: Region %d does not have a valid \"Radius\". Assigning value of \"Subvolume Base Size\".\n", numWarn++, curArrayItem);
				curSpec->subvol_spec[curArrayItem].radius = curSpec->SUBVOL_BASE_SIZE;
			} else
			{
				curSpec->subvol_spec[curArrayItem].radius = 
					cJSON_GetObjectItem(curObj, "Radius")->valuedouble;
			}
			
		}
		
		// Override region time step with global one
		curSpec->subvol_spec[curArrayItem].dt = curSpec->DT_MICRO;
		//curSpec->subvol_spec[curArrayItem].dt = 
		//	cJSON_GetObjectItem(curObj,
		//	"Time Step")->valuedouble;
	}
	
	// Check whether we might need a hybrid interface
	if(bHasMeso && bHasMicro)
	{
		if(cJSON_bItemValid(simControl,"Small Subvolumes at Hybrid Interface?", cJSON_True))
		{
			curSpec->B_HYBRID_SMALL_SUB =
				cJSON_GetObjectItem(simControl, "Small Subvolumes at Hybrid Interface?")->valueint;
		} else
		{ // Simulation does not have a valid Small Subvolumes at Hybrid Interface?
			bWarn = true;
			printf("WARNING %d: Simulation has at least 1 microscopic region and 1 mesoscopic region, but does not have a valid \"Small Subvolumes at Hybrid Interface?\". Assigning default value \"true\".\n", numWarn++);
			curSpec->B_HYBRID_SMALL_SUB = true;
		}
		
		if(cJSON_bItemValid(simControl,"Max Intrastep Micro to Meso Distance", cJSON_Number) &&
			cJSON_GetObjectItem(simControl, "Max Intrastep Micro to Meso Distance")->valuedouble >= 0.)
		{
			curSpec->MAX_HYBRID_DIST =
				cJSON_GetObjectItem(simControl, "Max Intrastep Micro to Meso Distance")->valuedouble;
		} else
		{ // Simulation does not have a valid Min Intrastep Micro to Meso Probability
			bWarn = true;
			printf("WARNING %d: Simulation has at least 1 microscopic region and 1 mesoscopic region, but does not have a valid \"Max Intrastep Micro to Meso Distance\". Assigning default value \"INFINITY\".\n", numWarn++);
			curSpec->MAX_HYBRID_DIST = INFINITY;
		}
	} else
	{
		// Warn if any hybrid properties are defined
		if(cJSON_bItemValid(simControl,
			"Small Subvolumes at Hybrid Interface?", cJSON_True))
		{
			bWarn = true;
			printf("WARNING %d: Simulation cannot have a hybrid interface, but \"Small Subvolumes at Hybrid Interface?\" was defined. Ignoring.\n", numWarn++);
		}
		if(cJSON_bItemValid(simControl,"Max Intrastep Micro to Meso Distance", cJSON_Number))
		{
			bWarn = true;
			printf("WARNING %d: Simulation cannot have a hybrid interface, but \"Max Intrastep Micro to Meso Distance\" was defined. Ignoring.\n", numWarn++);
		}
	}
	
	// Load actor details
	for(curArrayItem = 0;
		curArrayItem < curSpec->NUM_ACTORS; curArrayItem++)
	{			
		if(!cJSON_bArrayItemValid(actorSpec, curArrayItem, cJSON_Object))
		{
			fprintf(stderr, "ERROR: Actor %d is not described by a JSON object.\n", curArrayItem);
			exit(EXIT_FAILURE);
		}
		
		curObj = cJSON_GetArrayItem(actorSpec, curArrayItem);		
		
		if(!cJSON_bItemValid(curObj,"Is Location Defined by Regions?", cJSON_True))
		{ // Actor does not have a valid Is Location Defined by Regions?
			bWarn = true;
			printf("WARNING %d: Actor %d does not have a valid \"Is Location Defined by Regions?\". Assigning default value \"false\".\n", numWarn++, curArrayItem);
			curSpec->actorSpec[curArrayItem].bDefinedByRegions = false;
		} else
		{
			curSpec->actorSpec[curArrayItem].bDefinedByRegions = 
				cJSON_GetObjectItem(curObj, "Is Location Defined by Regions?")->valueint;
		}
		
		if(curSpec->actorSpec[curArrayItem].bDefinedByRegions)
		{
			// Set actor parameters that are not needed and see if config file still
			// defines them
			curSpec->actorSpec[curArrayItem].shape = UNDEFINED_SHAPE;
						
			if(cJSON_bItemValid(curObj,"Shape", cJSON_String))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Shape\" defined because its location is defined by regions. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Outer Boundary", cJSON_Array))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Outer Boundary\" defined because its location is defined by regions. Ignoring.\n", numWarn++, curArrayItem);
			}
			
			// Read regions that define location of actor
			if(!cJSON_bItemValid(curObj,"List of Regions Defining Location", cJSON_Array))
			{ // Actor does not have a List of Regions Defining Location array
				bWarn = true;
				printf("WARNING %d: Actor %d has a missing or invalid \"List of Regions Defining Location\". Assigning default value of \"0\" regions.\n", numWarn++, curArrayItem);
				curSpec->actorSpec[curArrayItem].numRegion = 0;
			} else
			{
				// Read number of regions
				curSpec->actorSpec[curArrayItem].numRegion =
					cJSON_GetArraySize(cJSON_GetObjectItem(curObj,"List of Regions Defining Location"));
				
				curSpec->actorSpec[curArrayItem].regionLabel =
					malloc(curSpec->actorSpec[curArrayItem].numRegion * sizeof(char *));
				if(curSpec->actorSpec[curArrayItem].regionLabel == NULL)
				{
					fprintf(stderr,"ERROR: Memory could not be allocated to store region list to define actor %d\n", curArrayItem);
					exit(EXIT_FAILURE);
				}
				
				// Read in names of regions							
				curObjInner = cJSON_GetObjectItem(curObj,"List of Regions Defining Location");
				for(i = 0; i < curSpec->actorSpec[curArrayItem].numRegion; i++)
				{
					if(!cJSON_bArrayItemValid(curObjInner,i, cJSON_String))
					{ // Exception region is not a valid string. Ignore
						bWarn = true;
						printf("WARNING %d: Actor %d region %d is not a valid string. Assigning empty string.\n", numWarn++, curArrayItem, i);
						curSpec->actorSpec[curArrayItem].regionLabel[i] = '\0';
					} else
					{
						curSpec->actorSpec[curArrayItem].regionLabel[i] =
							stringWrite(cJSON_GetArrayItem(curObjInner,i)->valuestring);
					}					
				}
			}
		} else
		{
			// Set actor parameters that are not needed and see if config file still
			// defines them
			curSpec->actorSpec[curArrayItem].numRegion = 0;
			curSpec->actorSpec[curArrayItem].regionLabel = NULL;
			
			if(cJSON_bItemValid(curObj,"List of Regions Defining Location", cJSON_Array))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"List of Regions Defining Location\" defined because its location is defined by an explicit shape. Ignoring.\n", numWarn++, curArrayItem);
			}
			
			if(!cJSON_bItemValid(curObj,"Shape", cJSON_String))
			{ // Actor does not have a defined Shape
				bWarn = true;
				printf("WARNING %d: Actor %d does not have a defined \"Shape\". Setting to default value \"Rectangular Box\".\n", numWarn++, curArrayItem);
				tempString =
					stringWrite("Rectangular Box");
			} else{
				tempString =
					stringWrite(cJSON_GetObjectItem(curObj,
					"Shape")->valuestring);
			}
			
			if(strcmp(tempString,"Rectangle") == 0)
			{
				curSpec->actorSpec[curArrayItem].shape = RECTANGLE;
				arrayLen = 6;
			}
			else if(strcmp(tempString,"Circle") == 0)
			{
				curSpec->actorSpec[curArrayItem].shape = CIRCLE;
				arrayLen = 4;
			}
			else if(strcmp(tempString,"Rectangular Box") == 0)
			{
				curSpec->actorSpec[curArrayItem].shape = RECTANGULAR_BOX;
				arrayLen = 6;
			}
			else if(strcmp(tempString,"Sphere") == 0)
			{
				curSpec->actorSpec[curArrayItem].shape = SPHERE;
				arrayLen = 4;
			}
			else if(strcmp(tempString,"Point") == 0)
			{
				curSpec->actorSpec[curArrayItem].shape = POINT;
				arrayLen = 3;
			}
			else
			{
				bWarn = true;
				printf("WARNING %d: Actor %d has an invalid \"Shape\". Setting to default value \"Rectangular Box\".\n", numWarn++, curArrayItem);
				curSpec->actorSpec[curArrayItem].shape = RECTANGULAR_BOX;
				arrayLen = 6;
			}
			
			if(!cJSON_bItemValid(curObj,"Outer Boundary", cJSON_Array) ||
				cJSON_GetArraySize(cJSON_GetObjectItem(curObj,"Outer Boundary")) != arrayLen)
			{
				bWarn = true;
				printf("WARNING %d: Actor %d has a missing or invalid \"Outer Boundary\". Setting to default value all \"0\"s.\n", numWarn++, curArrayItem);
				for(i = 0; i < arrayLen; i++)
				{
					curSpec->actorSpec[curArrayItem].boundary[i] = 0;
				}
			} else
			{
				curObjInner = cJSON_GetObjectItem(curObj,"Outer Boundary");
				for(i = 0; i < arrayLen; i++)
				{
					if(!cJSON_bArrayItemValid(curObjInner,i, cJSON_Number))
					{
						bWarn = true;
						printf("WARNING %d: Actor %d has an invalid \"Outer Boundary\" parameter %d. Setting to default value \"0\".\n", numWarn++, curArrayItem, i);
						curSpec->actorSpec[curArrayItem].boundary[i] = 0;
					} else
					{
						curSpec->actorSpec[curArrayItem].boundary[i] =
							cJSON_GetArrayItem(curObjInner, i)->valuedouble;
					}					
				}
			}
			
			// For spherical boundaries, confirm radius validity and add r^2 term
			if(strcmp(tempString,"Sphere") == 0)
			{
				if(curSpec->actorSpec[curArrayItem].boundary[3] <= 0.)
				{
					bWarn = true;
					printf("WARNING %d: Actor %d is a sphere but has invalid radius %f. Setting to default value \"1e-6\".\n",
						numWarn++, curArrayItem,
						curSpec->actorSpec[curArrayItem].boundary[3]);
					curSpec->actorSpec[curArrayItem].boundary[3] = 1e-6;
				}
				curSpec->actorSpec[curArrayItem].boundary[4] =
					squareDBL(curSpec->actorSpec[curArrayItem].boundary[3]);
			}
			free(tempString);
		}
		
		if(!cJSON_bItemValid(curObj,"Is Actor Active?", cJSON_True))
		{ // Actor does not have a valid Is Actor Active?
			bWarn = true;
			printf("WARNING %d: Actor %d does not have a valid \"Is Actor Active?\". Assigning default value \"false\".\n", numWarn++, curArrayItem);
			curSpec->actorSpec[curArrayItem].bActive = false;
		} else
		{
			curSpec->actorSpec[curArrayItem].bActive = 
				cJSON_GetObjectItem(curObj, "Is Actor Active?")->valueint;
		}
		
		if(!cJSON_bItemValid(curObj,"Start Time", cJSON_Number))
		{ // Actor does not have a valid Start Time
			bWarn = true;
			printf("WARNING %d: Actor %d does not have a valid \"Start Time\". Assigning default value \"0\".\n", numWarn++, curArrayItem);
			curSpec->actorSpec[curArrayItem].startTime = 0;
		} else
		{
			curSpec->actorSpec[curArrayItem].startTime = 
				cJSON_GetObjectItem(curObj, "Start Time")->valuedouble;
		}
		
		if(!cJSON_bItemValid(curObj,"Is There Max Number of Actions?", cJSON_True))
		{ // Actor does not have a valid Is There Max Number of Actions?
			bWarn = true;
			printf("WARNING %d: Actor %d does not have a valid \"Is There Max Number of Actions?\". Assigning default value \"false\".\n", numWarn++, curArrayItem);
			curSpec->actorSpec[curArrayItem].bMaxAction = false;
		} else
		{
			curSpec->actorSpec[curArrayItem].bMaxAction = 
				cJSON_GetObjectItem(curObj, "Is There Max Number of Actions?")->valueint;
		}
		
		if(curSpec->actorSpec[curArrayItem].bMaxAction)
		{
			if(!cJSON_bItemValid(curObj,"Max Number of Actions", cJSON_Number) ||
				cJSON_GetObjectItem(curObj,"Max Number of Actions")->valueint < 1)
			{ // Region does not have a valid Max Number of Actions
				bWarn = true;
				printf("WARNING %d: Region %d does not have a valid \"Max Number of Actions\". Assigning default value \"1\".\n", numWarn++, curArrayItem);
				curSpec->actorSpec[curArrayItem].numMaxAction = 1;
			} else
			{
				curSpec->actorSpec[curArrayItem].numMaxAction = 
					cJSON_GetObjectItem(curObj, "Max Number of Actions")->valueint;
			}
		} else
		{
			curSpec->actorSpec[curArrayItem].numMaxAction = 0;
			if(cJSON_bItemValid(curObj,"Max Number of Actions", cJSON_Number))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Max Number of Actions\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
		}
		
		if(!cJSON_bItemValid(curObj,"Is Actor Independent?", cJSON_True))
		{ // Actor does not have a valid Is Actor Independent?
			bWarn = true;
			printf("WARNING %d: Actor %d does not have a valid \"Is Actor Independent?\". Assigning default value \"true\".\n", numWarn++, curArrayItem);
			curSpec->actorSpec[curArrayItem].bIndependent = true;
		} else
		{
			curSpec->actorSpec[curArrayItem].bIndependent = 
				cJSON_GetObjectItem(curObj, "Is Actor Independent?")->valueint;
		}
		
		if(!cJSON_bItemValid(curObj,"Action Interval", cJSON_Number) ||
			cJSON_GetObjectItem(curObj,"Action Interval")->valuedouble <= 0)
		{ // Actor does not have a valid Action Interval
			bWarn = true;
			printf("WARNING %d: Actor %d does not have a valid \"Action Interval\". Assigning default value \"1\".\n", numWarn++, curArrayItem);
			curSpec->actorSpec[curArrayItem].actionInterval = 1.;
		} else
		{
			curSpec->actorSpec[curArrayItem].actionInterval = 
				cJSON_GetObjectItem(curObj, "Action Interval")->valuedouble;
		}
		
		if(!cJSON_bItemValid(curObj,"Is Actor Activity Recorded?", cJSON_True))
		{ // Actor does not have a valid Is Actor Activity Recorded?
			bWarn = true;
			printf("WARNING %d: Actor %d does not have a valid \"Is Actor Activity Recorded?\". Assigning default value \"true\".\n", numWarn++, curArrayItem);
			curSpec->actorSpec[curArrayItem].bWrite = true;
		} else
		{
			curSpec->actorSpec[curArrayItem].bWrite = 
				cJSON_GetObjectItem(curObj, "Is Actor Activity Recorded?")->valueint;
		}
		
		if(curSpec->actorSpec[curArrayItem].bActive)
		{ // Actor is active. Check for all active parameters
	
			// Check for existence of unnecessary parameters and display
			// warnings if they are defined.
			if(cJSON_bItemValid(curObj,"Is Time Recorded with Activity?", cJSON_True))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Is Time Recorded with Activity?\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Is Molecule Type Observed?", cJSON_Array))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Is Molecule Type Observed?\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Is Molecule Position Observed?", cJSON_Array))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Is Molecule Position Observed?\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			
			if(!cJSON_bItemValid(curObj,"Random Number of Molecules?", cJSON_True))
			{ // Actor does not have a valid value for Random Number of Molecules?
				bWarn = true;
				printf("WARNING %d: Actor %d does not have a valid value for \"Random Number of Molecules?\". Assigning default value \"false\".\n", numWarn++, curArrayItem);
				curSpec->actorSpec[curArrayItem].bNumReleaseRand = false;
			} else
			{
				curSpec->actorSpec[curArrayItem].bNumReleaseRand = 
					cJSON_GetObjectItem(curObj, "Random Number of Molecules?")->valueint;
			}
			
			// Random release times is only needed if "Random Number of Molecules?" is true
			if(curSpec->actorSpec[curArrayItem].bNumReleaseRand == true)
			{
				if(!cJSON_bItemValid(curObj,"Random Molecule Release Times?", cJSON_True))
				{ // Actor does not have a valid value for Random Molecule Release Times?
					bWarn = true;
					printf("WARNING %d: Actor %d does not have a valid value for \"Random Molecule Release Times?\". Assigning default value \"false\".\n", numWarn++, curArrayItem);
					curSpec->actorSpec[curArrayItem].bTimeReleaseRand = false;
				} else
				{
					curSpec->actorSpec[curArrayItem].bTimeReleaseRand = 
						cJSON_GetObjectItem(curObj, "Random Molecule Release Times?")->valueint;
				}
			} else
			{
				curSpec->actorSpec[curArrayItem].bTimeReleaseRand = false;
				if(cJSON_bItemValid(curObj,"Random Molecule Release Times?", cJSON_Number))
				{
					bWarn = true;
					printf("WARNING %d: Actor %d does not need \"Random Molecule Release Times?\" defined. Ignoring.\n", numWarn++, curArrayItem);
				}
			}
			
		
			if(!cJSON_bItemValid(curObj,"Release Interval", cJSON_Number) ||
				cJSON_GetObjectItem(curObj,"Release Interval")->valuedouble < 0)
			{ // Actor does not have a valid Action Interval
				bWarn = true;
				printf("WARNING %d: Actor %d does not have a valid \"Release Interval\". Assigning default value \"0\" seconds.\n", numWarn++, curArrayItem);
				curSpec->actorSpec[curArrayItem].releaseInterval = 0.;
			} else
			{
				curSpec->actorSpec[curArrayItem].releaseInterval = 
					cJSON_GetObjectItem(curObj, "Release Interval")->valuedouble;
			}
		
			// Slot interval is only needed if "Random Molecule Release Times?" is false
			if(curSpec->actorSpec[curArrayItem].bTimeReleaseRand == true)
			{ // Release times within release interval are random; no slot interval needed
				if(cJSON_bItemValid(curObj,"Slot Interval", cJSON_Number))
				{
					bWarn = true;
					printf("WARNING %d: Actor %d does not need \"Slot Interval\" defined. Ignoring.\n", numWarn++, curArrayItem);
				}
			} else
			{
				if(!cJSON_bItemValid(curObj,"Slot Interval", cJSON_Number) ||
					cJSON_GetObjectItem(curObj,"Slot Interval")->valuedouble < 0)
				{ // Actor does not have a valid Action Interval
					bWarn = true;
					printf("WARNING %d: Actor %d does not have a valid \"Slot Interval\". Assigning default value \"0\" seconds.\n", numWarn++, curArrayItem);
					curSpec->actorSpec[curArrayItem].slotInterval = 0.;
				} else
				{
					curSpec->actorSpec[curArrayItem].slotInterval = 
						cJSON_GetObjectItem(curObj, "Slot Interval")->valuedouble;
				}
			}
		
			if(!cJSON_bItemValid(curObj,"Modulation Scheme", cJSON_String))
			{ // Actor does not have a defined Modulation Scheme
				bWarn = true;
				printf("WARNING %d: Actor %d does not have a defined \"Modulation Scheme\". Setting to default value \"CSK\".\n", numWarn++, curArrayItem);
				tempString =
					stringWrite("CSK");
			} else{
				tempString =
					stringWrite(cJSON_GetObjectItem(curObj,
					"Modulation Scheme")->valuestring);
			}			
			
			if(strcmp(tempString,"CSK") == 0)
				curSpec->actorSpec[curArrayItem].modScheme = CSK;
			else if(strcmp(tempString,"Burst") == 0)
				curSpec->actorSpec[curArrayItem].modScheme = BURST;
			else
			{
				bWarn = true;
				printf("WARNING %d: Actor %d has an invalid \"Modulation Scheme\". Setting to default value \"CSK\".\n", numWarn++, curArrayItem);
				curSpec->actorSpec[curArrayItem].modScheme = CSK;
			}
			free(tempString);
		
			if(curSpec->actorSpec[curArrayItem].modScheme == BURST)
			{ // Modulation type does not use random bits
				
				curSpec->actorSpec[curArrayItem].bRandBits = true;
				curSpec->actorSpec[curArrayItem].probOne = 1.;
				curSpec->actorSpec[curArrayItem].modBits = 1;
				if(cJSON_bItemValid(curObj,"Bits Random?", cJSON_True))
				{
					bWarn = true;
					printf("WARNING %d: Actor %d does not need \"Bits Random?\" defined. Ignoring.\n", numWarn++, curArrayItem);
				}
				if(cJSON_bItemValid(curObj,"Bit Sequence", cJSON_Array))
				{
					bWarn = true;
					printf("WARNING %d: Actor %d does not need \"Bit Sequence\" defined. Ignoring.\n", numWarn++, curArrayItem);
				}
				if(cJSON_bItemValid(curObj,"Probability of Bit 1", cJSON_Number))
				{
					bWarn = true;
					printf("WARNING %d: Actor %d does not need \"Probability of Bit 1\" defined. Ignoring.\n", numWarn++, curArrayItem);
				}
				if(cJSON_bItemValid(curObj,"Modulation Bits", cJSON_Number))
				{
					bWarn = true;
					printf("WARNING %d: Actor %d does not need \"Modulation Bits\" defined. Ignoring.\n", numWarn++, curArrayItem);
				}
				
			} else
			{ // Modulation type uses random bits
			
				if(!cJSON_bItemValid(curObj,"Bits Random?", cJSON_True))
				{ // Actor does not have a valid value for Bits Random?
					bWarn = true;
					printf("WARNING %d: Actor %d does not have a valid value for \"Bits Random?\". Assigning default value \"true\".\n", numWarn++, curArrayItem);
					curSpec->actorSpec[curArrayItem].bRandBits = true;
				} else
				{
					curSpec->actorSpec[curArrayItem].bRandBits = 
						cJSON_GetObjectItem(curObj, "Bits Random?")->valueint;
				}
				
				if(curSpec->actorSpec[curArrayItem].bRandBits)
				{ // Bits are random. Need probability of bit being 1.
			
					if(!cJSON_bItemValid(curObj,"Probability of Bit 1", cJSON_Number) ||
						cJSON_GetObjectItem(curObj,"Probability of Bit 1")->valuedouble < 0. ||
						cJSON_GetObjectItem(curObj,"Probability of Bit 1")->valuedouble > 1.)
					{ // Actor does not have a valid Probability of Bit 1
						bWarn = true;
						printf("WARNING %d: Actor %d does not have a valid \"Probability of Bit 1\". Assigning default value \"0.5\".\n", numWarn++, curArrayItem);
						curSpec->actorSpec[curArrayItem].probOne = 0.5;
					} else
					{
						curSpec->actorSpec[curArrayItem].probOne = 
							cJSON_GetObjectItem(curObj, "Probability of Bit 1")->valuedouble;
					}
					
					// Actor should not have a defined bit sequence
					if(cJSON_bItemValid(curObj,"Bit Sequence", cJSON_Array))
					{
						bWarn = true;
						printf("WARNING %d: Actor %d does not need \"Bit Sequence\" defined. Ignoring.\n", numWarn++, curArrayItem);
					}
				} else
				{ // Need sequence of bits to use
					if(!cJSON_bItemValid(curObj,"Bit Sequence", cJSON_Array))
					{ // Actor does not have a valid Bit Sequence
						bWarn = true;
						printf("WARNING %d: Actor %d does not have a valid \"Bit Sequence\". Setting Max Number of Actions to \"0\".\n", numWarn++, curArrayItem);
						curSpec->actorSpec[curArrayItem].bBits = NULL;
						curSpec->actorSpec[curArrayItem].numMaxAction = 0;
						curSpec->actorSpec[curArrayItem].bMaxAction = true;
						arrayLen = 0;
					} else
					{ // Compare length of bit sequence with maximum number of actions
						arrayLen =
							cJSON_GetArraySize(cJSON_GetObjectItem(curObj,"Bit Sequence"));
						curSpec->actorSpec[curArrayItem].bBits = 
							malloc(arrayLen * sizeof(bool));
						
						if(curSpec->actorSpec[curArrayItem].bBits == NULL)
						{
							fprintf(stderr,"ERROR: Memory could not be allocated to store bit sequence of actor %d\n", curArrayItem);
							exit(EXIT_FAILURE);
						}
						
						curObjInner = cJSON_GetObjectItem(curObj,"Bit Sequence");
						for(i = 0; i < arrayLen; i++)
						{
							if(!(cJSON_bArrayItemValid(curObjInner,i, cJSON_Number) ||
								cJSON_bArrayItemValid(curObjInner,i, cJSON_True)) ||
								(cJSON_GetArrayItem(curObjInner,i)->valueint != 0 &&
								cJSON_GetArrayItem(curObjInner,i)->valueint != 1))
							{
								bWarn = true;
								printf("WARNING %d: Bit %d in the sequence of actor %d has an invalid value %d. Setting to default value of \"0\".\n", numWarn++, i, curArrayItem, cJSON_GetArrayItem(curObjInner,i)->valueint);
								curSpec->actorSpec[curArrayItem].bBits[i] = false;
							} else
							{
								curSpec->actorSpec[curArrayItem].bBits[i] =
									cJSON_GetArrayItem(curObjInner,i)->valueint;
							}
						}
					}
					
					// Actor should not have a defined bit probability
					if(cJSON_bItemValid(curObj,"Probability of Bit 1", cJSON_Number))
					{
						bWarn = true;
						printf("WARNING %d: Actor %d does not need \"Probability of Bit 1\" defined. Ignoring.\n", numWarn++, curArrayItem);
					}
				}		
				if(!cJSON_bItemValid(curObj,"Modulation Bits", cJSON_Number) ||
					cJSON_GetObjectItem(curObj,"Modulation Bits")->valueint < 1)
				{ // Region does not have a valid Modulation Bits
					bWarn = true;
					printf("WARNING %d: Region %d does not have a valid \"Modulation Bits\". Assigning default value \"1\".\n", numWarn++, curArrayItem);
					curSpec->actorSpec[curArrayItem].modBits = 1;
				} else
				{
					curSpec->actorSpec[curArrayItem].modBits = 
						cJSON_GetObjectItem(curObj, "Modulation Bits")->valueint;
				}
				
				if(!curSpec->actorSpec[curArrayItem].bRandBits)
				{ // Compare defined sequence length with number of actions and bits
					if(arrayLen % curSpec->actorSpec[curArrayItem].modBits != 0)
					{ // Sequence length does not accommodate number of modulation bits
						bWarn = true;
						printf("WARNING %d: Actor %d has a \"Bit Sequence\" with length %d but number of modulation bits is %d. Sequence has an insufficient number of bits for the last symbol, which will not be released.\n",
							numWarn++, curArrayItem, arrayLen, curSpec->actorSpec[curArrayItem].modBits);
					}
					if(curSpec->actorSpec[curArrayItem].bMaxAction
						&& curSpec->actorSpec[curArrayItem].numMaxAction
						!= arrayLen / curSpec->actorSpec[curArrayItem].modBits)
					{ // Sequence length does not correspond to maximum number of actions
						bWarn = true;
						printf("WARNING %d: Actor %d has a \"Bit Sequence\" with length %d but specified maximum number of actions is %d. Overriding maximum number of actions to \"%u\".\n",
							numWarn++, curArrayItem, arrayLen, curSpec->actorSpec[curArrayItem].numMaxAction, arrayLen / curSpec->actorSpec[curArrayItem].modBits);
					} else
					{
						curSpec->actorSpec[curArrayItem].bMaxAction = true;
					}
					curSpec->actorSpec[curArrayItem].numMaxAction =
						arrayLen / curSpec->actorSpec[curArrayItem].modBits;
				}
			}
		
			if(!cJSON_bItemValid(curObj,"Modulation Strength", cJSON_Number) ||
				cJSON_GetObjectItem(curObj,"Modulation Strength")->valuedouble <= 0.)
			{ // Actor does not have a valid Modulation Strength
				bWarn = true;
				printf("WARNING %d: Actor %d does not have a valid \"Modulation Strength\". Assigning default value \"1\".\n", numWarn++, curArrayItem);
				curSpec->actorSpec[curArrayItem].modStrength = 1.;
			} else
			{
				curSpec->actorSpec[curArrayItem].modStrength = 
					cJSON_GetObjectItem(curObj, "Modulation Strength")->valuedouble;
			}
			
			curSpec->actorSpec[curArrayItem].bReleaseMol =
				malloc(curSpec->NUM_MOL_TYPES * sizeof(bool));
			if(curSpec->actorSpec[curArrayItem].bReleaseMol == NULL)
			{
				fprintf(stderr,"ERROR: Memory could not be allocated to store array of booleans for actor %d, which is active.\n", curArrayItem);
				exit(EXIT_FAILURE);
			}
			
			if(!cJSON_bItemValid(curObj,"Is Molecule Type Released?", cJSON_Array) ||
				cJSON_GetArraySize(cJSON_GetObjectItem(curObj,"Is Molecule Type Released?")) != curSpec->NUM_MOL_TYPES)
			{ // Config file does not list a valid Is Molecule Type Released? array
				bWarn = true;
				printf("WARNING %d: Actor %d does not have a valid \"Is Molecule Type Released?\" array or not of correct length. Assigning default value \"true\" to first molecule type.\n", numWarn++, curArrayItem);
				curSpec->actorSpec[curArrayItem].bReleaseMol[0] = true;
				for(curMolType = 1; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
				{
					curSpec->actorSpec[curArrayItem].bReleaseMol[curMolType] = false;
				}
			} else{
				curObjInner = cJSON_GetObjectItem(curObj, "Is Molecule Type Released?");
				for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
				{
					if(!cJSON_bArrayItemValid(curObjInner,curMolType, cJSON_True))
					{
						bWarn = true;
						printf("WARNING %d: \"Is Molecule Type Released?\" %d of Actor %d not defined or has an invalid value. Assigning default value of \"false\".\n", numWarn++, curMolType, curArrayItem);
						curSpec->actorSpec[curArrayItem].bReleaseMol[curMolType] = false;
					} else
					{
						curSpec->actorSpec[curArrayItem].bReleaseMol[curMolType] =
							cJSON_GetArrayItem(curObjInner, curMolType)->valueint;
					}
				}
			}
		} else
		{ // Actor is passive. Check for all passive parameters
				
			// Check for existence of unnecessary parameters and display
			// warnings if they are defined.
			if(cJSON_bItemValid(curObj,"Random Number of Molecules?", cJSON_True))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Random Number of Molecules?\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Random Molecule Release Times?", cJSON_True))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Random Molecule Release Times?\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Release Interval", cJSON_Number))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Release Interval\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Slot Interval", cJSON_Number))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Slot Interval\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Bits Random?", cJSON_True))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Bits Random?\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Bit Sequence", cJSON_Array))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Bit Sequence\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Probability of Bit 1", cJSON_Number))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Probability of Bit 1\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Modulation Scheme", cJSON_String))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Modulation Scheme\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Modulation Bits", cJSON_Number))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Modulation Bits\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Modulation Strength", cJSON_Number))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Modulation Strength\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			if(cJSON_bItemValid(curObj,"Is Molecule Type Released?", cJSON_Array))
			{
				bWarn = true;
				printf("WARNING %d: Actor %d does not need \"Is Molecule Type Released?\" defined. Ignoring.\n", numWarn++, curArrayItem);
			}
			
			
			if(!cJSON_bItemValid(curObj,"Is Time Recorded with Activity?", cJSON_True))
			{ // Actor does not have a valid Is Time Recorded with Activity?
				bWarn = true;
				printf("WARNING %d: Actor %d does not have a valid \"Is Time Recorded with Activity?\". Assigning default value \"false\".\n", numWarn++, curArrayItem);
				curSpec->actorSpec[curArrayItem].bRecordTime = false;
			} else
			{
				curSpec->actorSpec[curArrayItem].bRecordTime = 
					cJSON_GetObjectItem(curObj, "Is Time Recorded with Activity?")->valueint;
			}
			
			curSpec->actorSpec[curArrayItem].bRecordMol =
				malloc(curSpec->NUM_MOL_TYPES * sizeof(bool));
			curSpec->actorSpec[curArrayItem].bRecordPos =
				malloc(curSpec->NUM_MOL_TYPES * sizeof(bool));
			if(curSpec->actorSpec[curArrayItem].bRecordMol == NULL ||
				curSpec->actorSpec[curArrayItem].bRecordPos == NULL)
			{
				fprintf(stderr,"ERROR: Memory could not be allocated to store array of booleans for actor %d, which is passive.\n", curArrayItem);
				exit(EXIT_FAILURE);
			}
			
			
			if(!cJSON_bItemValid(curObj,"Is Molecule Type Observed?", cJSON_Array) ||
				cJSON_GetArraySize(cJSON_GetObjectItem(curObj,"Is Molecule Type Observed?")) != curSpec->NUM_MOL_TYPES)
			{ // Config file does not list a valid Is Molecule Type Observed? array
				bWarn = true;
				printf("WARNING %d: Actor %d does not have a valid \"Is Molecule Type Observed?\" array or not of correct length. Assigning default value \"true\" to each molecule type.\n", numWarn++, curArrayItem);
				for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
				{
					curSpec->actorSpec[curArrayItem].bRecordMol[curMolType] = true;
				}
			} else{
				curObjInner = cJSON_GetObjectItem(curObj, "Is Molecule Type Observed?");
				for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
				{
					if(!cJSON_bArrayItemValid(curObjInner,curMolType, cJSON_True))
					{
						bWarn = true;
						printf("WARNING %d: \"Is Molecule Type Observed?\" %d of Actor %d not defined or has an invalid value. Assigning default value of \"true\".\n", numWarn++, curMolType, curArrayItem);
						curSpec->actorSpec[curArrayItem].bRecordMol[curMolType] = true;
					} else
					{
						curSpec->actorSpec[curArrayItem].bRecordMol[curMolType] =
							cJSON_GetArrayItem(curObjInner, curMolType)->valueint;
					}
				}
			}
			
			if(!cJSON_bItemValid(curObj,"Is Molecule Position Observed?", cJSON_Array) ||
				cJSON_GetArraySize(cJSON_GetObjectItem(curObj,"Is Molecule Position Observed?")) != curSpec->NUM_MOL_TYPES)
			{ // Config file does not list a valid Is Molecule Position Observed? array
				bWarn = true;
				printf("WARNING %d: Actor %d does not have a valid \"Is Molecule Position Observed?\" array or not of correct length. Assigning default value \"false\" to each molecule type.\n", numWarn++, curArrayItem);
				for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
				{
					curSpec->actorSpec[curArrayItem].bRecordPos[curMolType] = false;
				}
			} else{
				curObjInner = cJSON_GetObjectItem(curObj, "Is Molecule Position Observed?");
				for(curMolType = 0; curMolType < curSpec->NUM_MOL_TYPES; curMolType++)
				{
					if(!cJSON_bArrayItemValid(curObjInner,curMolType, cJSON_True))
					{
						bWarn = true;
						printf("WARNING %d: \"Is Molecule Position Observed?\" %d of Actor %d not defined or has an invalid value. Assigning default value of \"false\".\n", numWarn++, curMolType, curArrayItem);
						curSpec->actorSpec[curArrayItem].bRecordPos[curMolType] = false;
					} else
					{
						curSpec->actorSpec[curArrayItem].bRecordPos[curMolType] =
							cJSON_GetArrayItem(curObjInner, curMolType)->valueint;
					}
				}
			}
		}
	}
	
	// Cleanup
	cJSON_Delete(configJSON);
	free(configContent);
	free(configNameFull);
	
	// Pause for warnings if needed
	printf("Configuration file has %d warning(s). ", numWarn);
	if(bWarn  && !bWarnOverride)
	{
		printf("Press \'Enter\' to continue the simulation or \'q\'+\'Enter\' to quit.\n");
		ch = getchar();
		if(ch == 'q')
			exit(EXIT_FAILURE);
	} else printf("\n");
}

// Release memory allocated to configuration settings
void deleteConfig(struct simSpec3D curSpec)
{
	unsigned short curRegion, curActor, curMolType;
	unsigned short curRxn;
	
	if(curSpec.DIFF_COEF != NULL)
		free(curSpec.DIFF_COEF);
	if(curSpec.GLOBAL_FLOW_VECTOR != NULL)
		free(curSpec.GLOBAL_FLOW_VECTOR);
	if(curSpec.OUTPUT_NAME != NULL)
		free(curSpec.OUTPUT_NAME);
	if(curSpec.B_GLOBAL_MOL_FLOW != NULL)
		free(curSpec.B_GLOBAL_MOL_FLOW);
	
	if(curSpec.chem_rxn != NULL)
	{
		for(curRxn = 0; curRxn < curSpec.MAX_RXNS; curRxn++)
		{			
			if(curSpec.chem_rxn[curRxn].label != NULL)
				free(curSpec.chem_rxn[curRxn].label);
			if(curSpec.chem_rxn[curRxn].bReversible
				&& curSpec.chem_rxn[curRxn].labelCoupled != NULL)
				free(curSpec.chem_rxn[curRxn].labelCoupled);
			
			if(curSpec.chem_rxn[curRxn].reactants != NULL)
				free(curSpec.chem_rxn[curRxn].reactants);
			if(curSpec.chem_rxn[curRxn].products != NULL)
				free(curSpec.chem_rxn[curRxn].products);
			if(curSpec.chem_rxn[curRxn].bReleaseProduct != NULL)
				free(curSpec.chem_rxn[curRxn].bReleaseProduct);
			
			if(curSpec.chem_rxn[curRxn].regionExceptionLabel != NULL)
			{
				for(curRegion = 0; curRegion < curSpec.chem_rxn[curRxn].numRegionExceptions; curRegion++)
				{
					if(curSpec.chem_rxn[curRxn].regionExceptionLabel[curRegion] != NULL)
						free(curSpec.chem_rxn[curRxn].regionExceptionLabel[curRegion]);
				}
				free(curSpec.chem_rxn[curRxn].regionExceptionLabel);
			}
		}
		free(curSpec.chem_rxn);
	}
	
	if(curSpec.subvol_spec != NULL)
	{
		for(curRegion = 0; curRegion < curSpec.NUM_REGIONS; curRegion++)
		{
			if(curSpec.subvol_spec[curRegion].label != NULL)
				free(curSpec.subvol_spec[curRegion].label);
			if(curSpec.subvol_spec[curRegion].parent != NULL)
				free(curSpec.subvol_spec[curRegion].parent);
			if(curSpec.subvol_spec[curRegion].bLocalDiffusion &&
				curSpec.subvol_spec[curRegion].diffusion != NULL)
				free(curSpec.subvol_spec[curRegion].diffusion);
			for(curMolType = 0; curMolType < curSpec.NUM_MOL_TYPES; curMolType++)
			{
				if(curSpec.subvol_spec[curRegion].flowVector[curMolType] != NULL)
					free(curSpec.subvol_spec[curRegion].flowVector[curMolType]);
			}
			if(curSpec.subvol_spec[curRegion].bFlowLocal != NULL)
				free(curSpec.subvol_spec[curRegion].bFlowLocal);
			if(curSpec.subvol_spec[curRegion].bFlow != NULL)
				free(curSpec.subvol_spec[curRegion].bFlow);
			if(curSpec.subvol_spec[curRegion].flowType != NULL)
				free(curSpec.subvol_spec[curRegion].flowType);
			if(curSpec.subvol_spec[curRegion].flowVector != NULL)
				free(curSpec.subvol_spec[curRegion].flowVector);
		}
		free(curSpec.subvol_spec);
	}
	
	if(curSpec.actorSpec != NULL)
	{
		for(curActor = 0; curActor < curSpec.NUM_ACTORS; curActor++)
		{
			if(curSpec.actorSpec[curActor].bDefinedByRegions)
			{
				if(curSpec.actorSpec[curActor].regionLabel != NULL)
				{
					for(curRegion = 0; curRegion < curSpec.actorSpec[curActor].numRegion; curRegion++)
					{
						if(curSpec.actorSpec[curActor].regionLabel[curRegion] != NULL)
							free(curSpec.actorSpec[curActor].regionLabel[curRegion]);
					}
					free(curSpec.actorSpec[curActor].regionLabel);
				}
			}
			
			if(curSpec.actorSpec[curActor].bActive)
			{
				if(curSpec.actorSpec[curActor].bReleaseMol != NULL)
					free(curSpec.actorSpec[curActor].bReleaseMol);
				if(!curSpec.actorSpec[curActor].bRandBits
					&& curSpec.actorSpec[curActor].bBits != NULL)
					free(curSpec.actorSpec[curActor].bBits);
			} else
			{
				if(curSpec.actorSpec[curActor].bRecordMol != NULL)
					free(curSpec.actorSpec[curActor].bRecordMol);
				if(curSpec.actorSpec[curActor].bRecordPos != NULL)
					free(curSpec.actorSpec[curActor].bRecordPos);
			}
		}
		free(curSpec.actorSpec);
	}
}

// Initialize the simulation output file
void initializeOutput(FILE ** out,
	FILE ** outSummary,
	const char * CONFIG_NAME,
	const struct simSpec3D curSpec)
{
	time_t timer;
	char timeBuffer[26];
	struct tm* timeInfo;
	cJSON * root;
	char * outText;
	char * outDir2 = "../results";
	char * outDir1 = "results";
	bool bUseDir1;
	int mkdirOutput;
	char * outputNameFull;
	char * outputSummaryNameFull;
	unsigned int dirLength, nameLength;
	struct stat sb;
	
	time(&timer);
	timeInfo = localtime(&timer);
	
	strftime(timeBuffer, 26, "%Y-%m-%d %H:%M:%S", timeInfo);
	
	// Check existence of results folder and create it if it does not exist
	if (stat(outDir1, &sb) == 0 && S_ISDIR(sb.st_mode))
    { // Directory "results/" exists. Use it for output
		bUseDir1 = true;
    } else if (stat(outDir2, &sb) == 0 && S_ISDIR(sb.st_mode))
    { // Directory "../results/" exists. Use it for output
		bUseDir1 = false;
    } else
	{ // Create directory "results/" and use for output
		printf("NOTE: \"results\" directory could not be found. Trying to create.\n");
		#if defined(__linux__) || defined(__APPLE__)
			mkdirOutput = mkdir(outDir1, S_IRWXU);
		#else
			mkdirOutput = _mkdir(outDir1);
		#endif
		if(mkdirOutput == -1)
		{ // Directory could not be created
			fprintf(stderr,"ERROR: \"results\" directory could not be created.\n");
			exit(EXIT_FAILURE);
		}
		bUseDir1 = true;
	}
	
	// Construct full name of output file
	if(bUseDir1)
	{
		dirLength = strlen(outDir1)+1;
	} else
	{
		dirLength = strlen(outDir2)+1;
	}
	
	nameLength = strlen(curSpec.OUTPUT_NAME);
	outputNameFull = malloc(dirLength + nameLength + 5);
	outputSummaryNameFull = malloc(dirLength + nameLength + 23);
	if(outputNameFull == NULL
		|| outputSummaryNameFull == NULL)
	{
		fprintf(stderr,"ERROR: Memory could not be allocated to store the configuration file name\n");
		exit(EXIT_FAILURE);
	}
	outputNameFull[0] = '\0'; // Initialize full name to the empty string
	outputSummaryNameFull[0] = '\0';
	if(bUseDir1)
	{
		strcat(outputNameFull,outDir1);
		strcat(outputNameFull,"/");
	} else
	{
		strcat(outputNameFull,outDir2);
		strcat(outputNameFull,"/");
	}
	strcat(outputNameFull,curSpec.OUTPUT_NAME);
	strcat(outputSummaryNameFull,outputNameFull);
	strcat(outputNameFull,".txt");
	strcat(outputSummaryNameFull,"_summary.txt");
	
	printf("Simulation output will be written to \"%s\".\n", outputNameFull);
	printf("Simulation summary will be written to \"%s\".\n", outputSummaryNameFull);
	
	if((*out = fopen(outputNameFull, "w")) == NULL)
	{
		fprintf(stderr,"ERROR: Cannot create output file \"%s\".\n",outputNameFull);
		exit(EXIT_FAILURE);
	}
	if((*outSummary = fopen(outputSummaryNameFull, "w")) == NULL)
	{
		fprintf(stderr,"ERROR: Cannot create output summary file \"%s\".\n",
			outputSummaryNameFull);
		exit(EXIT_FAILURE);
	}
	
	root = cJSON_CreateObject();
	cJSON_AddStringToObject(root, "ConfigFile", CONFIG_NAME);
	cJSON_AddNumberToObject(root, "SEED", curSpec.SEED);
	cJSON_AddNumberToObject(root, "NumRepeat", curSpec.NUM_REPEAT);
	cJSON_AddStringToObject(root, "StartTime", timeBuffer);
	
	outText = cJSON_Print(root);
	fprintf(*outSummary, "%s", outText);
	fprintf(*outSummary, "\n");
	cJSON_Delete(root);
	free(outText);
	free(outputNameFull);
	free(outputSummaryNameFull);	
}

// Copy string (with memory allocation)
char * stringWrite(char * src)
{
	char * string = stringAllocate(strlen(src));
	strcpy(string, src);
	
	return string;
}

// Allocate memory for a string
char * stringAllocate(long stringLength)
{
	char * string = malloc(stringLength+1);
	if(string == NULL)
	{
		fprintf(stderr,"ERROR: Memory could not be allocated for string copy.\n");
		exit(EXIT_FAILURE);
	}
	
	return string;
}

// Print simulation output from one realization
void printOneTextRealization(FILE * out,
	const struct simSpec3D curSpec,
	unsigned int curRepeat,
	ListObs3D observationArray[],
	short numPassiveRecord,
	short * passiveRecordID,
	short numActiveRecord,
	short * activeRecordID,
	const struct actorStruct3D actorCommonArray[],
	const struct actorActiveStruct3D actorActiveArray[],
	const struct actorPassiveStruct3D actorPassiveArray[],	
	uint32_t maxActiveBits[],
	uint32_t maxPassiveObs[])
{
	short curActor, curActorPassive, curActorRecord, curActorActive;
	char * outText;
	NodeData * curData;
	NodeObs3D * curObs;
	unsigned short curMolInd, curMolType;
	ListMol3D * curMolList;
	NodeMol3D * curMolNode;
	uint32_t curActiveBits, curPassiveObs;
	
	fprintf(out, "Realization %u:\n", curRepeat);
	
	// Record active actor binary data
	for(curActorRecord = 0; curActorRecord < numActiveRecord;
		curActorRecord++)
	{
		curActor = activeRecordID[curActorRecord];
		curActorActive = actorCommonArray[curActor].activeID;
		curData =  actorActiveArray[curActorActive].binaryData.head;
		curActiveBits = 0;
		fprintf(out, "\tActiveActor %u:\n\t\t", curActor);
		while(curData != NULL)
		{
			fprintf(out, "%u ", curData->item.bit);
			curData = curData->next;
			curActiveBits++;
		}		
		fprintf(out, "\n");
		
		if (curActiveBits > maxActiveBits[curActorActive])
			maxActiveBits[curActorActive] = curActiveBits;
	}
	
	// Record observations by passive actors that are being recorded
	for(curActorRecord = 0; curActorRecord < numPassiveRecord; curActorRecord++)
	{
		// Actor in common actor list is passiveRecordID[curActorRecord]
		curActor = passiveRecordID[curActorRecord];
		curObs = (&observationArray[curActorRecord])->head;
		fprintf(out, "\tPassiveActor %u:\n", curActor);
		
		// Preliminary scan to find number of observations and compare
		// with largest number of observations made thus far in any realization
		curPassiveObs = 0;
		while(curObs != NULL)
		{
			curPassiveObs++;
			curObs = curObs->next;
		}		
		if (curPassiveObs > maxPassiveObs[curActorRecord])
			maxPassiveObs[curActorRecord] = curPassiveObs;
		
		// Record actor observation times (if being recorded)
		curObs = (&observationArray[curActorRecord])->head;
		if (actorCommonArray[curActor].spec.bRecordTime)
		{
			fprintf(out, "\t\tTime:\n\t\t\t");
			while(curObs != NULL)
			{
				fprintf(out, "%.4e ", curObs->item.paramDouble[0]);
				curObs = curObs->next;
			}
			fprintf(out, "\n");
		}
		
		// Record observations associated with each type of molecule being recorded
		curActorPassive = actorCommonArray[curActor].passiveID;
		for(curMolInd = 0;
			curMolInd < actorPassiveArray[curActorPassive].numMolRecordID;
			curMolInd++)
		{
			curMolType = actorPassiveArray[curActorPassive].molRecordID[curMolInd];
			fprintf(out, "\t\tMolID %u:\n\t\t\tCount:\n\t\t\t\t", curMolType);
			
			// Record molecule counts made by observer
			curObs = (&observationArray[curActorRecord])->head;
			while(curObs != NULL)
			{
				fprintf(out, "%" PRIu64 " ", curObs->item.paramUllong[curMolInd]);
				curObs = curObs->next;
			}
			fprintf(out, "\n");
			
			// Record molecule coordinates if specified
			if (actorCommonArray[curActor].spec.bRecordPos[curMolType])
			{
				curObs = (&observationArray[curActorRecord])->head;
				fprintf(out, "\t\t\tPosition:");
				while(curObs != NULL)
				{		
					fprintf(out, "\n\t\t\t\t");
					// Each observation will have the positions of some number of molecules
					fprintf(out, "(");
					curMolList = curObs->item.molPos[curMolInd];
					if(!isListMol3DEmpty(curMolList))
					{
						curMolNode = *curMolList;
						while(curMolNode != NULL)
						{
							fprintf(out, "(%e, %e, %e) ",
								curMolNode->item.x, curMolNode->item.y, curMolNode->item.z);
							curMolNode = curMolNode->next;
						}
					}						
					fprintf(out, ")");		
					curObs = curObs->next;
				}
				fprintf(out, "\n");
			}
		}
	}
	fprintf(out, "\n");
}

// Print end of simulation data
void printTextEnd(FILE * out,	
	short numActiveRecord,
	short numPassiveRecord,
	const struct actorStruct3D actorCommonArray[],
	const struct actorActiveStruct3D actorActiveArray[],
	const struct actorPassiveStruct3D actorPassiveArray[],
	short * passiveRecordID,
	short * activeRecordID,
	uint32_t maxActiveBits[],
	uint32_t maxPassiveObs[],
	double runTime)
{
	time_t timer;
	char timeBuffer[26];
	struct tm* timeInfo;
	short curActor, curPassive, curActorRecord;
	cJSON * root;
	cJSON * curArray, *curItem, *newItem, *newActor, *innerArray;
	char * outText;
	unsigned short curMolInd;
	
	time(&timer);
	timeInfo = localtime(&timer);
	
	strftime(timeBuffer, 26, "%Y-%m-%d %H:%M:%S", timeInfo);
	
	root = cJSON_CreateObject();
	
	// Store information about the active actors
	cJSON_AddNumberToObject(root, "NumberActiveRecord", numActiveRecord);
	cJSON_AddItemToObject(root, "ActiveInfo", curArray=cJSON_CreateArray());
	for(curActorRecord = 0; curActorRecord < numActiveRecord; curActorRecord++)
	{
		curActor = activeRecordID[curActorRecord];
		newActor = cJSON_CreateObject();
		cJSON_AddNumberToObject(newActor, "ID", curActor);
		cJSON_AddNumberToObject(newActor, "MaxBitLength",
			maxActiveBits[curActorRecord]);
		cJSON_AddItemToArray(curArray, newActor);
	}
	
	// Store information about the passive actors that were recorded
	cJSON_AddNumberToObject(root, "NumberPassiveRecord", numPassiveRecord);	
	cJSON_AddItemToObject(root, "RecordInfo", curArray=cJSON_CreateArray());
	for(curActorRecord = 0; curActorRecord < numPassiveRecord; curActorRecord++)
	{
		curActor = passiveRecordID[curActorRecord];
		curPassive = actorCommonArray[curActor].passiveID;
		newActor = cJSON_CreateObject();
		// Record Passive Actor IDs that are being recorded
		cJSON_AddNumberToObject(newActor, "ID", curActor);
		cJSON_AddNumberToObject(newActor, "bRecordTime",
			actorCommonArray[curActor].spec.bRecordTime);
		// Record maximum number of observations made by each recorded actor
		cJSON_AddNumberToObject(newActor, "MaxCountLength",
			maxPassiveObs[curActorRecord]);
		cJSON_AddNumberToObject(newActor, "NumMolTypeObs",
			actorPassiveArray[curPassive].numMolRecordID);
		cJSON_AddItemToObject(newActor, "MolObsID", innerArray=cJSON_CreateArray());
		for(curMolInd = 0; curMolInd < actorPassiveArray[curPassive].numMolRecordID;
			curMolInd++)
		{
			cJSON_AddItemToArray(innerArray,
				cJSON_CreateNumber(
				actorPassiveArray[curPassive].molRecordID[curMolInd]));
		}
		cJSON_AddItemToObject(newActor, "bRecordPos", innerArray=cJSON_CreateArray());
		for(curMolInd = 0; curMolInd < actorPassiveArray[curPassive].numMolRecordID;
			curMolInd++)
		{
			cJSON_AddItemToArray(innerArray,
				cJSON_CreateNumber(
				actorCommonArray[curActor].spec.bRecordPos[actorPassiveArray[curPassive].molRecordID[curMolInd]]));
		}		
		cJSON_AddItemToArray(curArray, newActor);
	}
	
	cJSON_AddStringToObject(root, "EndTime", timeBuffer);
	
	cJSON_AddNumberToObject(root, "RunTime", runTime);
	
	outText = cJSON_Print(root);
	fprintf(out, "%s", outText);
	cJSON_Delete(root);
	free(outText);
}

