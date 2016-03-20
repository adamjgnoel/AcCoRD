/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2015 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * observations.h - 	linked list of observations made by a passive actor
 *
 * Last revised for AcCoRD v0.4.1
 *
 * Revision history:
 *
 * Revision v0.4.1
 * - improved use and format of error messages
 *
 * Revision v0.3.1
 * - header added
 *
 * Created 2015-02-23
*/
#ifndef OBSERVATIONS_H
#define OBSERVATIONS_H

#include <stdio.h> // to create and edit files
#include <stdlib.h> // for exit(), malloc, free, NULL
#include <stdbool.h> // for C++ bool naming, requires C99
#include "micro_molecule.h" // for list of molecule positions

// observations specific declarations

// paramUllong = {numMol1 ... numMolN}
// paramDouble = {recordTime, mol1Pos1X, mol1Pos1Y, mol1Pos1Z,
// 		mol1Pos2X, mol1Pos2Y, mol1Pos2Z, ... molNPosMX, molNPosMY, molNPosMZ}
struct observation_list3D {
	unsigned short numDouble; // number of double parameters
	unsigned short numUllong; // number of uint64_t parameters
	double * paramDouble; // array of double parameters
	uint64_t * paramUllong; // array of uint64_t parameters
	ListMol3D ** molPos; // Linked list of molecule coordinates
};

// General (linked list) type declarations

typedef struct observation_list3D ItemObs3D;

typedef struct node_Obs3D{
	ItemObs3D item;
	struct node_Obs3D * next;
} NodeObs3D;

typedef struct list_Obs3D{
	unsigned short numMolTypeObs; // Number of types of molecules being observed
	NodeObs3D * head; // Pointer to first observation
	NodeObs3D * tail; // Pointer to most recent observation
} ListObs3D;

// observations specific Prototypes

bool addObservation(ListObs3D * list,
	const unsigned short numDouble,
	const unsigned short numUllong,
	double paramDouble[],
	uint64_t paramUllong[],
	ListMol3D molPos[]);

// General Prototypes

void initializeListObs(ListObs3D * list,
	unsigned short numMolTypeObs);

bool isListEmptyObs(const ListObs3D * list);

void emptyListObs(ListObs3D * list);


#endif // OBSERVATIONS_H
