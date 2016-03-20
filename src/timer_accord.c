/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * timer_accord.c - 	manage passage of time in simulation, including
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
#include "timer_accord.h"

//
// "Private" Declarations
//

//
// Definitions
//

/* Allocate space for an array of timer values
*/
void allocateTimerArray(const short NUM_TIMERS,
	struct timerStruct ** timerArray)
{
	
	*timerArray = malloc(NUM_TIMERS*sizeof(struct timerStruct));
		
	if(*timerArray == NULL){
		fprintf(stderr, "ERROR: Memory allocation for array of timer structures.\n");
		exit(EXIT_FAILURE);
	}
}

/* Allocate space for the heap of IDs pointing to the array of timers
 * timer structures, a 2D array listing the child IDs of each heap element,
 * and a 2D array of bools indicating whether each child is valid
*/
void allocateTimerHeapArray(const short NUM_TIMERS,
	short ** heapTimer,
	short (**heapTimerChildID)[2],
	bool (**b_heapTimerChildValid)[2])
{
	
	*heapTimer = malloc(NUM_TIMERS*sizeof(short));
	*heapTimerChildID = malloc(NUM_TIMERS*sizeof(short [2]));
	*b_heapTimerChildValid = malloc(NUM_TIMERS*sizeof(bool [2]));
		
	if(*heapTimer == NULL || *heapTimerChildID == NULL || *b_heapTimerChildValid == NULL)
	{
		fprintf(stderr, "ERROR: Memory allocation for heap of timer structures.\n");
		exit(EXIT_FAILURE);
	}
}

void deleteTimerHeapArray(struct timerStruct timerArray[])
{
	if(timerArray != NULL) free(timerArray);
}

// Initialize Array of pointers to Timer values
void initializeTimerArray(const short NUM_TIMERS,
	struct timerStruct timerArray[])
{
	short curTimer;
	
	for(curTimer=0; curTimer < NUM_TIMERS; curTimer++)
	{
		timerArray[curTimer].nextTime = 0.;
		timerArray[curTimer].heapID = 0;
	}
}

// Reset initial times for all timers
void resetTimerArray(const short NUM_TIMERS,
	struct timerStruct timerArray[],
	const short NUM_ACTORS,
	const struct actorStruct3D actorCommonArray[],
	const short MESO_TIMER_ID,
	const double tMeso,
	const short MICRO_TIMER_ID,
	const double DT_MICRO)
{
	short curTimer;
	
	for(curTimer = 0; curTimer < NUM_ACTORS; curTimer++)
	{
		timerArray[curTimer].nextTime = actorCommonArray[curTimer].spec.startTime;
	}
	
	timerArray[MESO_TIMER_ID].nextTime = tMeso;
	timerArray[MICRO_TIMER_ID].nextTime = DT_MICRO;
}

// Update next time of timer
// NOTE: heapTimerUpdate should be called IMMEDIATELY after (i.e., before another
// call to this function), otherwise the heap won't be properly sorted
void updateTimer(const short curTimer,
	const short NUM_TIMERS,
	struct timerStruct timerArray[],
	double tCur)
{
	return;
}

// Build 2D array listing children of heap elements
void heapTimerFindChildren(const short NUM_TIMERS,
	short heapTimerChildID[][2],
	bool b_heapTimerChildValid[][2])
{
	short i;
	
	for(i = 0; i < NUM_TIMERS; i++){ // parent == i
		heapTimerChildID[i][0] = (short) 2*i + 1;
		heapTimerChildID[i][1] = (short) 2*i + 2;
		
		if (heapTimerChildID[i][1] < NUM_TIMERS){ // Both children valid
			b_heapTimerChildValid[i][0] = true;
			b_heapTimerChildValid[i][1] = true;
		} else if (heapTimerChildID[i][0] < NUM_TIMERS){ // First child valid
			b_heapTimerChildValid[i][0] = true;
			b_heapTimerChildValid[i][1] = false;
		} else { // No valid children
			b_heapTimerChildValid[i][0] = false;
			b_heapTimerChildValid[i][1] = false;
		}			
	}
	
	return;
}

// Free memory of heap arrays
void heapTimerDelete(short heapTimer[],
	short heapTimerChildID[][2],
	bool b_heapTimerChildValid[][2])
{	
	if(heapTimer != NULL) free(heapTimer);
	if(heapTimerChildID != NULL) free(heapTimerChildID);
	if(b_heapTimerChildValid != NULL) free(b_heapTimerChildValid);
}

// Build Heap of Timers
void heapTimerBuild(const short NUM_TIMERS,
	struct timerStruct timerArray[],
	short heapTimer[],
	const unsigned int numHeapLevels,
	short heapTimerChildID[][2],
	bool b_heapTimerChildValid[][2])
{
	int i,j;
	short levelFirst, levelLast; // First and last elements in current level
	short parent, parentNew, parentOld;
	
	// Initialize heap with ordered (but unsorted) indices
	for(i=0;i<NUM_TIMERS;i++){
		heapTimer[i] = i;
		timerArray[i].heapID = i;
	}
	
	// Make heap a "min-heap" using next reaction times
	for(i = numHeapLevels-2; i >= 0; i--){ // For each level in heap
		levelFirst = 1;
		for(j=0; j<i; j++) levelFirst *= 2;
		levelLast = 2*levelFirst - 2; // Last element is 2^(i+1)-2
		levelFirst--; // First element is 2^i - 1
		for(parent = levelFirst; parent <= levelLast; parent++){
			// Compare element with its two children
			parentOld = parent;
			parentNew = heapTimerCompareDown(NUM_TIMERS, timerArray, heapTimer,
				parent, heapTimerChildID, b_heapTimerChildValid);
			while(parentNew != parentOld){
				// Keep comparing with children until this parent is the min
				parentOld = parentNew;
				parentNew = heapTimerCompareDown(NUM_TIMERS, timerArray, heapTimer,
					parentOld, heapTimerChildID, b_heapTimerChildValid);
			}
		}
	}
}

// Update Placement of Single Heap Element Based on Updated Value
short heapTimerUpdate(const short NUM_TIMERS,
	struct timerStruct timerArray[],
	short heapTimer[],
	const short heapID,
	short heapTimerChildID[][2],
	bool b_heapTimerChildValid[][2])
{
	short newID, oldID; // Updated current and previous placement in heap
	
	// See if element needs to move "down" (i.e., lower priority)
	oldID = heapID;
	newID = heapTimerCompareDown(NUM_TIMERS, timerArray, heapTimer, heapID,
		heapTimerChildID, b_heapTimerChildValid);
	while(newID != oldID){
		oldID = newID;
		newID = heapTimerCompareDown(NUM_TIMERS, timerArray, heapTimer, oldID,
			heapTimerChildID, b_heapTimerChildValid);
	}
	if (newID == heapID && heapID > 0){
		// Element did not move down and it is not already at the top of the heap.
		// See if it needs to move "up" (i.e., higher priority)
		newID = heapTimerCompareUp(NUM_TIMERS, timerArray, heapTimer, heapID);
		while(newID != oldID){
			oldID = newID;
			newID = heapTimerCompareUp(NUM_TIMERS, timerArray, heapTimer, oldID);
		}
	}
	return newID;
}

// Compare parent node with its children. Swap if parent does not have smallest value
// Return (possibly new) location of parent node
short heapTimerCompareDown(const short NUM_TIMERS,
	struct timerStruct timerArray[],
	short heapTimer[],
	const short parent,
	short heapTimerChildID[][2],
	bool b_heapTimerChildValid[][2])
{
		
	short child1 = heapTimerChildID[parent][0];
	short child2 = heapTimerChildID[parent][1];
	
	if (b_heapTimerChildValid[parent][1]){
		if (timerArray[heapTimer[parent]].nextTime <
			timerArray[heapTimer[child1]].nextTime){
			if (timerArray[heapTimer[parent]].nextTime >
				timerArray[heapTimer[child2]].nextTime){
				// child2 is smallest
				heapTimerSwap(timerArray,heapTimer,parent,child2);
				return child2;
			} // else parent is smallest (do nothing)
		} else if (timerArray[heapTimer[child2]].nextTime <
			timerArray[heapTimer[child1]].nextTime){
			// child2 is smallest
			heapTimerSwap(timerArray,heapTimer,parent,child2);
			return child2;
		} else {
			// child1 is smallest
			heapTimerSwap(timerArray,heapTimer,parent,child1);
			return child1;
		}
	} else if (b_heapTimerChildValid[parent][0]
		&& timerArray[heapTimer[parent]].nextTime >
		timerArray[heapTimer[child1]].nextTime) {
		// child1 is smallest
		heapTimerSwap(timerArray,heapTimer,parent,child1);
		return child1;
	}
	return parent;
}

// Compare node with its parent. Swap if parent has a larger value
// Return (possibly new) location of node
short heapTimerCompareUp(const short NUM_TIMERS,
	struct timerStruct timerArray[],
	short heapTimer[],
	const short child){
	
	if (child == 0) // Node is at top of heap. No parent exists
		return child;
	
	short parent = child/2;
	if(child % 2 == 0) // Element is even
		parent--;
	
	if (timerArray[heapTimer[parent]].nextTime >
		timerArray[heapTimer[child]].nextTime){
		// Child is smaller. Swap
		heapTimerSwap(timerArray,heapTimer,parent,child);
		return parent;
	} else // Parent is smaller. Do not swap
		return child;
}

// Swap the positions of two elements in the heap
void heapTimerSwap(struct timerStruct timerArray[],
	short heapTimer[],
	short index1,
	short index2){
	
	short tempID;
	
	// Update heap IDs associated with the timers whose heap positions are being swapped
	timerArray[heapTimer[index1]].heapID = index2;
	timerArray[heapTimer[index2]].heapID = index1;
	
	tempID = heapTimer[index1];
	heapTimer[index1] = heapTimer[index2];
	heapTimer[index2] = tempID;
}
