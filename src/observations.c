/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2015 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * observations.c - 	linked list of observations made by a passive actor
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

#include "observations.h"

// Local Function Prototypes

static bool addItem(ItemObs3D item, ListObs3D * list);

static void removeItem(NodeObs3D * prevNode, NodeObs3D * curNode);

static void traverse(const ListObs3D * list, void (* p_fun)(ItemObs3D item));

static void copyToNode(ItemObs3D item, NodeObs3D * p_node);

// Specific Definitions

// Create new observation with the corresponding information
bool addObservation(ListObs3D * list,
	const unsigned short numDouble,
	const unsigned short numUllong,
	double * paramDouble,
	uint64_t * paramUllong,
	ListMol3D molPos [])
{
	NodeMol3D * molPosList;
	unsigned short curMolInd;
	
	// Allocate memory
	unsigned short curData;
	double * paramDoubleNew = malloc(numDouble * sizeof(double));
	uint64_t * paramUllongNew =
		malloc(numUllong * sizeof(uint64_t));
	ListMol3D ** molPosListNew = malloc(list->numMolTypeObs * sizeof(ListMol3D *));
	if(paramDoubleNew == NULL || paramDoubleNew == NULL || molPosListNew == NULL)
		return false;
	
	for(curMolInd = 0; curMolInd < list->numMolTypeObs; curMolInd++)
	{
		molPosListNew[curMolInd] = malloc(sizeof(ListMol3D));
		if(molPosListNew[curMolInd] == NULL)
			return false;
		initializeListMol(molPosListNew[curMolInd]);
	}
	
	
	// Copy array data
	for(curData = 0; curData < numDouble; curData++)
	{
		paramDoubleNew[curData] = paramDouble[curData];
	}
	for(curData = 0; curData < numUllong; curData++)
	{
		paramUllongNew[curData] = paramUllong[curData];
	}
	
	// Copy molecule coordinate lists
	for(curMolInd = 0; curMolInd < list->numMolTypeObs; curMolInd++)
	{
		molPosList = molPos[curMolInd];
		if(!isListMol3DEmpty(&molPos[curMolInd]))
		{
			while(molPosList != NULL)
			{
				if(!addMolecule(molPosListNew[curMolInd],
					molPosList->item.x, molPosList->item.y, molPosList->item.z))
				{
					// Creation of molecule failed
					fprintf(stderr, "ERROR: Memory allocation to record molecule positions.\n");
					exit(EXIT_FAILURE);
				}
				molPosList = molPosList->next;
			}
		} else
		{
			*(molPosListNew[curMolInd]) = NULL;
		}
	}
	
	ItemObs3D newObs3D = {numDouble,numUllong,paramDoubleNew,paramUllongNew,molPosListNew};
	return addItem(newObs3D, list);
}

// General Definitions

// Initialize list
void initializeListObs(ListObs3D * list,
	unsigned short numMolTypeObs)
{
	list->numMolTypeObs = numMolTypeObs;
	list->head = NULL;
	list->tail = NULL;
}

// Is the list empty?
bool isListEmptyObs(const ListObs3D * list)
{
	if(list->head == NULL) return true;
	return false;
}

// Create new node to hold item and add it to the end of the list
bool addItem(ItemObs3D item, ListObs3D * list)
{
	NodeObs3D * p_new;
	
	p_new = malloc(sizeof(NodeObs3D));
	if (p_new == NULL)
		return false;	// Quit on failure of malloc
		
	copyToNode(item, p_new);
	p_new->next = NULL;
	if(list->tail == NULL)
	{ // List was empty. This item is first in list
		list->head = p_new;
		list->tail = p_new;
	} else
	{
		list->tail->next = p_new; // Point end of list to new node
		list->tail = p_new;
	}
	return true;
}

// Remove node from list. In order to do this efficiently, need to pass Nodes
// and not the start of the list.
// TODO: Update this function for current list type (it does not currently have the
// correct free calls but the function is not called so it is not a problem)
void removeItem(NodeObs3D * prevNode, NodeObs3D * curNode)
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

// Visit each node and execute the function pointed to by p_fun
// p_fun must be a void function that takes an item as input
void traverse(const ListObs3D * list, void (* p_fun)(ItemObs3D item))
{
	NodeObs3D * p_node = list->head;
	
	while(p_node != NULL)
	{
		(*p_fun)(p_node->item); // Apply function
		p_node = p_node->next;  // Advance to next item
	}
}

// De-allocate memory of all nodes in the list
void emptyListObs(ListObs3D * list)
{
	NodeObs3D * p_save;
	NodeObs3D * p_cur;
	unsigned short curMolInd;
	
	while(list->head != NULL)
	{
		p_save = list->head->next;	// Save address of next node
		free(list->head->item.paramDouble);
		free(list->head->item.paramUllong);
		for(curMolInd = 0; curMolInd < list->numMolTypeObs; curMolInd++)
		{
			if(!isListMol3DEmpty(list->head->item.molPos[curMolInd]))
				emptyListMol(list->head->item.molPos[curMolInd]);
			free(list->head->item.molPos[curMolInd]);
		}
		free(list->head->item.molPos);
		free(list->head);				// Free memory of current node
		list->head = p_save;			// Advance to next node
	}
}

// Copy an item to a node
static void copyToNode(ItemObs3D item, NodeObs3D * p_node)
{
	p_node->item = item; // Structure copy
}