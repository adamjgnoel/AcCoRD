/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2015 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * global_param.h - global parameters that are independent of a specific
 * 					simulation
 * Last revised for AcCoRD v0.3.1
 *
 * Revision history:
 *
 * Revision v0.3.1
 * - header added
*/
#ifndef GLOBAL_PARAM_H
#define GLOBAL_PARAM_H

// Define global parameters that are independent of a specific
// simulation.

#define PI 3.14159265359

// Maximum number of types of molecules in a single simulation.
// Used to initialize the structure defining chemical reactions.
#define MAX_MOL_TYPES 20

// Maximum number of products that a single chemical reaction can have
#define MAX_RXN_PRODUCTS 4

// Shape Indicators
// NOTE: Changes to names must be reflected in file_io.c
#define RECTANGLE 0
#define CIRCLE 1
#define RECTANGULAR_BOX 2
#define SPHERE 3
#define UNDEFINED_SHAPE 4

// Resolution of adjacency
// (i.e., what fraction of base subvolume size do the edges of shapes need to be
// in order to be declared adjacent?)
#define SUB_ADJ_RESOLUTION 0.01

// Directions
#define LEFT 0
#define RIGHT 1
#define DOWN 2
#define UP 3
#define IN 4
#define OUT 5
#define PARENT 6 // (i.e., outer)
#define CHILD 7 // (i.e., inner)
#define UNDEFINED 8

// Modulation schemes
// NOTE: Changes to names must be reflected in file_io.c
#define CSK 0

#endif // GLOBAL_PARAM_H
