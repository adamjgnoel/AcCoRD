/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2015 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * actor_data.c - linked list of binary data associated with active actor
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

#include "actor_data.h"

// Local Function Prototypes

static bool addItem(ItemData item, ListData * list);

static void removeItem(NodeData * prevNode, NodeData * curNode);

static void traverse(const ListData * list, void (* p_fun)(ItemData item));

static void copyToNode(ItemData item, NodeData * p_node);

// Specific Definitions

// Create new data bit
bool addData(ListData * list,
	const bool bit)
{
	
	ItemData newData = {bit};
	return addItem(newData, list);
}

// Copy new list of data bits to end of old list
void transferData(ListData * oldData,
	ListData * newData)
{
	NodeData * p_node = newData->head;
	
	while(p_node != NULL)
	{
		if(!addData(oldData, p_node->item.bit))
		{
			// Creation of bit failed
			fprintf(stderr, "ERROR: Memory could not be allocated to add bit to data sequence.\n");
			exit(EXIT_FAILURE);
		}
		p_node = p_node->next;
	}
}

// Convert data in binary list to decimal
unsigned int binaryToDecimal(const ListData * list,
	unsigned int alphabetSize)
{
	unsigned int curValue = 0U;
	alphabetSize /= 2U;
	NodeData * curNode = list->head;
	
	while (curNode != NULL)
	{
		curValue += ((unsigned int) curNode->item.bit) * alphabetSize;
		alphabetSize /= 2U;
		curNode = curNode->next;
	}
		
	return curValue;
}

// Recursively convert partial data list from binary to decimal. TODO: May not need
unsigned int binarySubListToDecimal(const NodeData * firstNode,
	unsigned int * maxValue)
{
	unsigned int curMaxValue, curValue;
	if (firstNode == NULL)
	{
		*maxValue = 1U;
		return 0;
	}
	
	curValue = binarySubListToDecimal(firstNode->next, &curMaxValue);
	curValue += ((unsigned int) firstNode->item.bit) * curMaxValue;
	*maxValue = 2*curMaxValue;
	return curValue;
}
	
// General Definitions

// Initialize list
void initializeListData(ListData * list)
{
	list->head = NULL;
	list->tail = NULL;
}

// Is the list empty?
bool isListDataEmpty(const ListData * list)
{
	if(list->head == NULL) return true;
	return false;
}

// Create new node to hold item and add it to the end of the list
bool addItem(ItemData item, ListData * list)
{
	NodeData * p_new;
	
	p_new = malloc(sizeof(NodeData));
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
void removeItem(NodeData * prevNode, NodeData * curNode)
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
void traverse(const ListData * list, void (* p_fun)(ItemData item))
{
	NodeData * p_node = list->head;
	
	while(p_node != NULL)
	{
		(*p_fun)(p_node->item); // Apply function
		p_node = p_node->next;  // Advance to next item
	}
}

// De-allocate memory of all nodes in the list
void emptyListData(ListData * list)
{
	NodeData * p_save;
	NodeData * p_cur;
	
	while(list->head != NULL)
	{
		p_save = list->head->next;	// Save address of next node
		free(list->head);				// Free memory of current node
		list->head = p_save;			// Advance to next node
	}
}

// Copy an item to a node
static void copyToNode(ItemData item, NodeData * p_node)
{
	p_node->item = item; // Structure copy
}