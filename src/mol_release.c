/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2015 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * mol_release.c - 	linked list of "current" molecule releases associated
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

#include "mol_release.h"

// Local Function Prototypes

static bool addItem(ItemRelease item, ListRelease * list);

static void removeItem(NodeRelease * prevNode, NodeRelease * curNode);

static void traverse(const ListRelease * list, void (* p_fun)(ItemRelease item));

static void copyToNode(ItemRelease item, NodeRelease * p_node);

// Specific Definitions

// Create new molecule emission with the corresponding information
bool addRelease(ListRelease * list,
	const double strength,
	const unsigned short molType,
	double startTime,
	double endTime,
	double frequency)
{	
	ItemRelease newRelease = {strength, molType, startTime, endTime, frequency};
	return addItem(newRelease, list);
}

void deleteRelease(ListRelease * list,
	const unsigned int releaseInd)
{
	NodeRelease * curRelease = *list;
	NodeRelease * prevRelease = NULL;
	NodeRelease * nextRelease;
	unsigned int curReleaseInd;
	
	// Traverse release list to the current release
	for(curReleaseInd = 0;
		curReleaseInd < releaseInd;
		curReleaseInd++)
	{
		prevRelease = curRelease;
		curRelease = curRelease->next;
	}
	nextRelease = curRelease->next;
	
	removeItem(prevRelease, curRelease);
	if(prevRelease == NULL)
	{ // We removed the first release. nextRelease is now the start of the list
		* list = nextRelease;
	}
}

// General Definitions

// Initialize list
void initializeListRelease(ListRelease * list)
{
	* list = NULL;
}

// Is the list empty?
bool isListReleaseEmpty(const ListRelease * list)
{
	if(list == NULL) return true;
	return false;
}

// Create new node to hold item and add it to the start of the list
bool addItem(ItemRelease item, ListRelease * list)
{
	NodeRelease * p_new;
	
	p_new = malloc(sizeof(NodeRelease));
	if (p_new == NULL)
		return false;	// Quit on failure of malloc
		
	copyToNode(item, p_new);
	p_new->next = * list;
	* list = p_new; // Point start of list to new node
	return true;
}

// Remove node from list. In order to do this efficiently, need to pass Nodes
// and not the start of the list.
void removeItem(NodeRelease * prevNode, NodeRelease * curNode)
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
void traverse(const ListRelease * list, void (* p_fun)(ItemRelease item))
{
	NodeRelease * p_node = *list;
	
	while(p_node != NULL)
	{
		(*p_fun)(p_node->item); // Apply function
		p_node = p_node->next;  // Advance to next item
	}
}

// De-allocate memory of all nodes in the list
void emptyListRelease(ListRelease * list)
{
	NodeRelease * p_save;
	
	while(*list != NULL)
	{
		p_save = (*list)->next;	// Save address of next node
		free(*list);				// Free memory of current node
		*list = p_save;			// Advance to next node
	}
}

// Copy an item to a node
static void copyToNode(ItemRelease item, NodeRelease * p_node)
{
	p_node->item = item; // Structure copy
}
