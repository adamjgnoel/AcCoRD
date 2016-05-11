/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2015 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * mol_release.h - 	linked list of "current" molecule releases associated
 * 					with an active actor
 *
 * Last revised for AcCoRD v0.3.1
 *
 * Revision history:
 *
 * Revision v0.3.1
 * - header added
 *
 * Created 2015-04-08
*/
#ifndef MOL_RELEASE_H
#define MOL_RELEASE_H

#include <stdio.h> // to create and edit files
#include <stdlib.h> // for exit(), malloc, free, NULL
#include <stdbool.h> // for C++ bool naming, requires C99
#include "global_param.h" // For modulation scheme #define

// release specific declarations

struct singleRelease {
	double strength; // Number of molecules to release or release rate
	unsigned short molType; // Index of molecule being released
	double nextTime; // Global next time for a molecule release corresponding to this emission
	double endTime; // Global end time associated with this emission
	double frequency; // Release frequency (if modulation has an associated frequency)
};

// General (linked list) type declarations

typedef struct singleRelease ItemRelease;

typedef struct node_Release{
	ItemRelease item;
	struct node_Release * next;
} NodeRelease;

typedef NodeRelease * ListRelease;

// observations specific Prototypes

bool addRelease(ListRelease * list,
	const double strength,
	const unsigned short molType,
	double startTime,
	double endTime,
	double frequency);

void deleteRelease(ListRelease * list,
	const unsigned int releaseInd);

// General Prototypes

void initializeListRelease(ListRelease * list);

bool isListReleaseEmpty(const ListRelease * list);

void emptyListRelease(ListRelease * list);


#endif // MOL_RELEASE_H
