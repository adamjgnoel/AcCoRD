/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2015 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * actor_data.h - linked list of binary data associated with active actor
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
 * Created 2015-04-08
*/
#ifndef ACTOR_DATA_H
#define ACTOR_DATA_H

#include <stdio.h> // to create and edit files
#include <stdlib.h> // for exit(), malloc, free, NULL
#include <stdbool.h> // for C++ bool naming, requires C99

// data specific declarations
struct curData {
	bool bit;
};

// General (linked list) type declarations

typedef struct curData ItemData;

typedef struct node_Data{
	ItemData item;
	struct node_Data * next;
} NodeData;

typedef struct list_data{
	NodeData * head; // Pointer to first observation
	NodeData * tail; // Pointer to most recent observation
} ListData;

// data specific Prototypes

bool addData(ListData * list,
	const bool bit);

void transferData(ListData * oldData,
	ListData * newData);
	
unsigned int binaryToDecimal(const ListData * list,
	unsigned int alphabetSize);

unsigned int binarySubListToDecimal(const NodeData * firstNode,
	unsigned int * maxValue);

// General Prototypes

void initializeListData(ListData * list);

bool isListDataEmpty(const ListData * list);

void emptyListData(ListData * list);


#endif // ACTOR_DATA_H
