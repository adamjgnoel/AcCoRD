/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * actor.h - operations on array of actors and its elements
 * Last revised for AcCoRD v0.4
 *
 * Revision history:
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

#ifndef ACTOR_H
#define ACTOR_H

#include <stdio.h>
#include <stdlib.h> // for exit(), malloc
#include <stdbool.h> // for C++ bool naming, requires C99
#include <limits.h> // For SHRT_MAX
#include "region.h"
#include "base.h" // For region adjacency
#include "micro_molecule.h" // For molecule coordinate linked list
#include "global_param.h" // For region adjacency
#include "mol_release.h" // For each active actor's linked list of active emissions
#include "actor_data.h" // For each active actor's linked list of binary data
#include "subvolume.h"  // For subvolume structure

/*
* Data Type Declarations
*/

/* The actorStructSpec3D structure contains the user-defined actor parameters
* for a single 3D actor
*/
struct actorStructSpec3D { // Configuration parameters
	
	/*
	*  Parameters Common to Any Actor
	*/
	
	// Actor shape type
	// Values defined in global_param.h
	int shape;
	
	// Outer boundary of this actor, defined in absolute space
	// TODO: Account for shapes besides rectangles (at least circles)
	double boundary[6];
	
	// Actors are "Active" or "Passive"
	// Active -> change a state of the simulation environment
	// Passive -> observe a state of the simulation environment
	bool bActive;
	
	// Actor's "Start time", i.e., when it is first supposed to act
	double startTime;
	
	// Does actor have a maximum number of times to act?
	bool bMaxAction;
	
	// Maximum number of actions (if bMaxAction == true)
	uint32_t numMaxAction;
	
	// Is actor behavior independent of the system state?
	bool bIndependent;
	
	// Frequency of independent behavior (if bIndependent == true)
	double actionInterval;
	
	/*
	*  Parameters for Active Actors
	*/
	
	// Is the number of molecules released stochastic?
	bool bNumReleaseRand;
	
	// Are the molecule release times stochastic?
	bool bTimeReleaseRand;
	
	// Period for release activity associated with one actionInterval
	// Could (in general) be less than or greater than actionInterval
	double releaseInterval;
	
	// Slot period (for "slotted" release, as opposed to kinetic)
	double slotInterval;
	
	// Is the actor activity defined by independent random bits?
	// TODO: Need array for defined bits if bRandBits == FALSE
	bool bRandBits;
	
	// If bRandBits, what is the probability of one bit having value 1
	double probOne;
	
	// Modulation scheme. Values are defined in global_param.h
	unsigned short modScheme;
	
	// Modulation bits (i.e., modulation order is 2^modBits)
	unsigned short modBits;
	
	// Modulation base strength. Exact meaning depends on the
	// modulation scheme
	double modStrength;
	
	// Which types of molecules are released?
	// The number of types with value 1 should be consistent with the modulation scheme
	bool bReleaseMol[MAX_MOL_TYPES];
	
	// TODO: Consider parameters for recording the actions of this actor (i.e., the
	// randomly-modulated bits, random signal strength, action times, etc)
		
	/*
	*  Parameters for Passive Actors
	*/
	
	// Record observations of actor to output file
	bool bWrite;
	
	// Log observation time (if bWrite == true)
	bool bRecordTime;
	
	// Which molecule types are observed? (if bWrite == true)
	bool bRecordMol[MAX_MOL_TYPES];
	
	// Which molecule types have positions recorded?
	// (if bWrite == true AND bRecordMol[ID] == true)
	bool bRecordPos[MAX_MOL_TYPES];
};

/* The actorStruct3D structure contains all parameters specific to any 3D
* single simulation actor.
*/
struct actorStruct3D { // Common actor parameters

	// The user-defined common actor parameters
	struct actorStructSpec3D spec;

	//
	// "Initialization" parameters (Determined from configuration parameters at runtime)
	//
	
	// IDs of actor in active/passive list
	short activeID;
	short passiveID;
	
	// Size of actor
	// TODO: May not need to be a permanent structure member; initialization only
	double volume;
	
	// The number of regions that are (at least partially) within actor space
	unsigned short numRegion;
	
	// Array of IDs of regions that are (at least partially) within actor space
	// Length is numRegion
	unsigned short * regionID;
	
	// Array of bools indicating whether the actor is fully inside an
	// overlapping region or not
	// Only applies if the actor is spherical, and if false then the region is
	// fully inside the actor
	// Length is numRegion
	bool * bRegionInside;
	
	// The number of subvolumes in each region that are (at least partially)
	// within actor space. Only applies if region is mesoscopic.
	// Length is numRegion
	uint32_t * numSub;
	
	// Shape of intersection of actor and region
	// Length is numRegion
	unsigned short * regionInterType;
	
	// Outer boundary of intersection of actor and region
	// TODO: May not need to be a permanent COMMON structure member; initialization
	// only (unless region is MICRO)
	// Length is numRegion x 6
	double (* regionInterBound)[6];
	
	// Area of intersection of actor and region
	// Length is numRegion
	// TODO: May not need to be a permanent structure member; initialization only
	double * regionInterArea;
	
	// Fraction of actor in intersection of actor and region
	// Length is numRegion
	double * cumFracActorInRegion;
	
	// The IDs of subvolumes in each region that are (at least partially)
	// within actor space. Only applies if region is mesoscopic.
	// Length is numRegion x numSub
	uint32_t ** subID;
	
	//
	// "Simulation" parameters (Determined at simulation time)
	//
	
	// NOTE: These parameters would need to be reset for each repeat

	// (Current) next "Action" time of actor
	double nextTime;
	
	// Index of current action (needed if bMaxAction == true)
	uint32_t curAction;
	
	// FUTURE MEMBERS (POTENTIAL)
	// Indicator for how next time is determined
	// Controls for dependent actors
};

/* The actorActiveStruct3D structure contains all parameters specific to an active 3D
* simulation actor.
*/
struct actorActiveStruct3D { // Active actor parameters

	//
	// "Initialization" parameters (Determined from configuration parameters at runtime)
	//
	
	// ID of actor in common (active and passive) list
	short actorID;
	
	// Size of symbol alphabet
	unsigned int alphabetSize;
	
	// Number of types of molecules that actor could release
	unsigned short numMolType;
	
	// Indices of the molecule types that actor could release
	// NOTE: The molecule types need to be defined in the order corresponding to the bit
	unsigned short molType[MAX_MOL_TYPES];
	
	// The cumulative fraction of actor that is within a given mesoscopic subvolume, given that the
	// subvolume is within actor space. Only applies if region is mesoscopic.
	// Length is numRegion x numSub
	double ** cumFracActorInSub;
	
	//
	// "Simulation" parameters (Determined at simulation time)
	//
	
	// Time for start of next NEW release by the actor (i.e., start of next symbol interval)
	double nextNewReleaseTime;
	
	// Time for next emission by a current release of the actor
	double nextEmissionTime;
	
	// Is next action time for a new symbol?
	bool bNextActionNewRelease;

	// Index of next emission by a current release of the actor
	unsigned int nextEmissionIndex;
	
	// Linked list for currently active release intervals
	ListRelease releaseList;
	
	// Random data associated with the actor
	ListData binaryData;
};

/* The actorPassiveStruct3D structure contains all parameters specific to a passive 3D
* simulation actor.
*/
struct actorPassiveStruct3D { // Passive actor parameters

	//
	// "Initialization" parameters (Determined from configuration parameters at runtime)
	//
	
	// ID of actor in common (active and passive) list
	short actorID;
	
	// The fraction of subvolume that is within actor space, given that a mesoscopic
	// subvolume is within actor space. Only applies if region is mesoscopic.
	// Length is numRegion x numSub
	double ** fracSubInActor;
	
	// Index of actor in list of passive actors that record observations
	short recordID;
	
	// Are we writing the position of at LEAST one molecule type in a mesoscopic region?
	// Indicates whether we need to define subInterBound
	// Length is numRegion
	bool * bRecordMesoAnyPos;
	
	// Outer boundary of intersection of actor and subvolume
	// Only allocated if region is mesoscopic and molecule positions are being recorded
	// Length is numRegion x numSub x 6
	double (** subInterBound)[6];
	
	// Number of molecule types that this actor records (if bWrite == true)
	unsigned short numMolRecordID;
	
	// IDs of molecules that this actor records counts of (if bWrite == true)
	// Length is numMolRecordID
	unsigned short * molRecordID;
	
	// Number of molecule types that this actor records positions
	// (if bWrite == true AND bRecordMol == true)
	unsigned short numMolRecordPosID;
	
	// IDs of molecules that this actor records positions of
	// (if bWrite == true AND bRecordMol == true)
	// Length is numMolRecordPosID
	unsigned short * molRecordPosID;
	
	//
	// "Simulation" parameters (Determined at simulation time)
	//
	
	// Number of molecules observed in the current observation
	// Length is numMolRecordID
	uint64_t * curMolObs;
};

/*
* Function Declarations
*/

void allocateActorCommonArray3D(const short NUM_ACTORS,
	struct actorStruct3D ** actorCommonArray);

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
	const double SUBVOL_BASE_SIZE);

void allocateActorActivePassiveArray3D(const short NUM_ACTORS_ACTIVE,
	struct actorActiveStruct3D ** actorActiveArray,
	const short NUM_ACTORS_PASSIVE,
	struct actorPassiveStruct3D ** actorPassiveArray);

void initializeActorActivePassive3D(const short NUM_ACTORS,
	const struct actorStruct3D actorCommonArray[],
	const unsigned short NUM_MOL_TYPES,
	const struct region regionArray[],
	const short NUM_REGIONS,
	const short NUM_ACTORS_ACTIVE,
	struct actorActiveStruct3D actorActiveArray[],
	const short NUM_ACTORS_PASSIVE,
	struct actorPassiveStruct3D actorPassiveArray[],
	uint32_t subCoorInd[][3]);

void resetActors3D(const short NUM_ACTORS,
	struct actorStruct3D actorCommonArray[],
	const unsigned short NUM_MOL_TYPES,
	const struct region regionArray[],
	const short NUM_REGIONS,
	const short NUM_ACTORS_ACTIVE,
	struct actorActiveStruct3D actorActiveArray[],
	const short NUM_ACTORS_PASSIVE,
	struct actorPassiveStruct3D actorPassiveArray[]);

void deleteActor3D(const short NUM_ACTORS,
	struct actorStruct3D actorCommonArray[],
	const struct region regionArray[],
	const short NUM_ACTORS_ACTIVE,
	struct actorActiveStruct3D actorActiveArray[],
	const short NUM_ACTORS_PASSIVE,
	struct actorPassiveStruct3D actorPassiveArray[],
	short actorRecordID[]);

// TODO: Move these release functions into a "higher" level file
	
void newRelease3D(const struct actorStruct3D * actorCommon,
	struct actorActiveStruct3D * actorActive,
	double curTime);

void findNextEmission3D(const struct actorStruct3D * actorCommon,
	struct actorActiveStruct3D * actorActive);

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
	bool (*b_heap_childValid)[2]);
	
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
	bool (*b_heap_childValid)[2]);

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
	bool (*b_heap_childValid)[2]);
	
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
	bool (*b_heap_childValid)[2]);

#endif // ACTOR_H
