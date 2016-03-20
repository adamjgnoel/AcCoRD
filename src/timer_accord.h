/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * timer_accord.h - 	manage passage of time in simulation, including
 *						heap to store the next action time by all actors
 * 						and the micro and meso simulation regimes
 *
 * Last revised for AcCoRD v0.4.1
 *
 * Revision history:
 *
 * Revision v0.4.1
 * - improved use and format of error messages
 *
 * Revision v0.4
 * - removed deprecated debug functions
 * - renamed region structure and array of structures to accommodate naming
 * of new boundary structure
 *
 * Revision v0.3.1
 * - header added
 *
 * Created 2015-03-02
*/
#ifndef TIMER_ACCORD_H
#define TIMER_ACCORD_H

#include <stdio.h> // for printf, scanf
#include <stdlib.h> // for exit(), malloc
#include <stdbool.h> // for C++ bool conventions
#include <math.h> // for INFINITY
#include "actor.h"  // for actor start times
#include "region.h" // for region parameters (each micro region has an associated time)

//
// Constant definitions
//


//
// Subvolume data type declarations
//


/* The timerStruct structure lists properties for a single timer,
* The entire simulation environment has only one array of these structures.
*/
struct timerStruct{
	
	short heapID; // Index of timer in timer heap
	
	double nextTime; // Timer value (time of next action)
	
	// FUTURE MEMBERS
};

//
// Function Declarations
//

// Allocate Array of pointers to Timer values
void allocateTimerArray(const short NUM_TIMERS,
	struct timerStruct ** timerArray);

// Allocate Heap (sorted timer IDs) and associated arrays
void allocateTimerHeapArray(const short NUM_TIMERS,
	short ** heapTimer,
	short (**heapTimerChildID)[2],
	bool (**b_heapTimerChildValid)[2]);

// Free memory allocated to array of Timer values
void deleteTimerHeapArray(struct timerStruct timerArray[]);

// Initialize Array of pointers to Timer values
void initializeTimerArray(const short NUM_TIMERS,
	struct timerStruct timerArray[]);
	
// Reset initial times for all timers
void resetTimerArray(const short NUM_TIMERS,
	struct timerStruct timerArray[],
	const short NUM_ACTORS,
	const struct actorStruct3D actorCommonArray[],
	const short MESO_TIMER_ID,
	const double tMeso,
	const short MICRO_TIMER_ID,
	const double DT_MICRO);

// Update next time of timer
void updateTimer(const short curTimer,
	const short NUM_TIMERS,
	struct timerStruct timerArray[],
	double tCur);

// Heap Operations
void heapTimerFindChildren(const short NUM_TIMERS,
	short heapTimerChildID[][2],
	bool b_heapTimerChildValid[][2]);
	
void heapTimerDelete(short heapTimer[],
	short heapTimerChildID[][2],
	bool b_heapTimerChildValid[][2]);
	
void heapTimerBuild(const short NUM_TIMERS,
	struct timerStruct timerArray[],
	short heapTimer[],
	const unsigned int numHeapLevels,
	short heapTimerChildID[][2],
	bool b_heapTimerChildValid[][2]);
	
short heapTimerUpdate(const short NUM_TIMERS,
	struct timerStruct timerArray[],
	short heapTimer[],
	const short heapID,
	short heapTimerChildID[][2],
	bool b_heapTimerChildValid[][2]);
	
short heapTimerCompareDown(const short NUM_TIMERS,
	struct timerStruct timerArray[],
	short heapTimer[],
	const short parent,
	short heapTimerChildID[][2],
	bool b_heapTimerChildValid[][2]);
	
short heapTimerCompareUp(const short NUM_TIMERS,
	struct timerStruct timerArray[],
	short heapTimer[],
	const short child);
	
void heapTimerSwap(struct timerStruct timerArray[],
	short heapTimer[],
	short index1, short index2);
	
#endif // TIMER_ACCORD_H
