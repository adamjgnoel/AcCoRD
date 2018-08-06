/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * global_param.h - global parameters that are independent of a specific
 * 					simulation
 *
 * Last revised for AcCoRD v1.4 (2018-08-06)
 *
 * Revision history:
 *
 * Revision v1.4 (2018-08-06)
 * - added parameters to implement a priori surface reactions, including chemical
 * reaction type RXN_A_PRIORI_ABSORBING, surface reaction probability types
 * RXN_PROB_A_PRIORI_SPHERE and RXN_PROB_A_PRIORI_INFINITE_PLANE, and types of
 * a priori surface reaction thresholds (RXN_THRESHOLD_DISTANCE and
 * RXN_THRESHOLD_PROB)
 *
 * Revision v1.1 (2016-12-24)
 * - added flow types FLOW_NONE and FLOW_UNIFORM
 *
 * Revision v1.0 (2016-10-31)
 * - added BURST modulation
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - added POINT as a shape
 *
 * Revision v0.5.1 (2016-05-06)
 * - added types of chemical reaction probability calculations
 * - removed MAX_MOL_TYPES, MAX_RXN_PRODUCTS
 *
 * Revision v0.5 (2016-04-15)
 * - removed use of MAX_MOL_TYPES
 * - removed use of MAX_RXN_PRODUCTS
 * - added types of regions and surfaces
 *
 * Revision v0.3.1
 * - header added
*/
#ifndef GLOBAL_PARAM_H
#define GLOBAL_PARAM_H

// Define global parameters that are independent of a specific
// simulation.

#define PI 3.14159265359

// Dimensions
#define DIM_X 0
#define DIM_Y 1
#define DIM_Z 2

// Shape Indicators
// NOTE: Changes to list of names must be reflected in file_io.c
// if regions or actors can have the updated shapes
#define RECTANGLE 0
#define CIRCLE 1
#define RECTANGULAR_BOX 2
#define SPHERE 3
#define LINE 4
#define POINT 5
#define UNDEFINED_SHAPE 6

// Types of regions
// NOTE: Changes to list of names must be reflected in file_io.c
#define REGION_NORMAL 0
#define REGION_SURFACE_3D 1
#define REGION_SURFACE_2D 2

// Types of surfaces
// NOTE: Changes to list of names must be reflected in file_io.c
#define NO_SURFACE 0
#define SURFACE_INNER 1
#define SURFACE_OUTER 2
#define SURFACE_MEMBRANE 3

// Planes for 2d shapes
// NOTE: Changes to list of names must be reflected in region.c
#define PLANE_3D 0
#define PLANE_XY 1
#define PLANE_XZ 2
#define PLANE_YZ 3

// Effective dimension for volume/area
// Depends on shape and type
#define DIM_1D 1
#define DIM_2D 2
#define DIM_3D 3

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
// NOTE: Changes to list of names must be reflected in file_io.c
#define CSK 0
#define BURST 1

// 1st Order Surface reaction types
// NOTE: Changes to list of names must be reflected in file_io.c
#define RXN_NORMAL 0
#define RXN_ABSORBING 1
#define RXN_DESORBING 2
#define RXN_RECEPTOR 3
#define RXN_MEMBRANE_IN 4
#define RXN_MEMBRANE_OUT 5
#define RXN_A_PRIORI_ABSORBING 6

// Types of reaction probability calculations
// NOTE: Changes to list of names must be reflected in file_io.c
#define RXN_PROB_NORMAL 0
#define RXN_PROB_MIXED 1
#define RXN_PROB_STEADY_STATE 2
#define RXN_PROB_A_PRIORI_SPHERE 3
#define RXN_PROB_A_PRIORI_INFINITE_PLANE 4

// Types of thresholding for A Priori surface reactions
#define RXN_THRESHOLD_DISTANCE 0
#define RXN_THRESHOLD_PROB 1
#define RXN_THRESHOLD_REGION 2

// Types of molecule placement strategies when leaving surface
// NOTE: Changes to list of names must be reflected in file_io.c
#define PROD_PLACEMENT_LEAVE 0 // Leave molecule next to surface
#define PROD_PLACEMENT_FORCE 1 // Force diffusion away from surface
#define PROD_PLACEMENT_STEADY_STATE 2 // Force diffuse based on steady-state values

// Flow types
// NOTE: Changes to list of names must be reflected in file_io.c
#define FLOW_NONE 0
#define FLOW_UNIFORM 1

#endif // GLOBAL_PARAM_H
