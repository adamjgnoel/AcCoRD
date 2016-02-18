/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * base.h - general utility functions that can apply to different simulation data
 * 			structures
 * Last revised for AcCoRD v0.4.1
 *
 * Revision history:
 *
 * Revision v0.4.1
 * - improved use and format of error messages
 *
 * Revision v0.4
 * - filled in cases for spherical boundaries to accommodate spherical regions and actors
 * - added clearance to bBoundaryIntersect to account for cases where we need to know
 * whether two shapes are within a specified distance of each other and not just
 * overlapping
 * - added functions for
 * 		- distance between 2 points
 *		- determining whether one shape is completely inside another shape
 * 		- squaring values
 *		- distance to nearest point on boundary
 * 		- determining intersection of line and boundary
 *
 * Revision v0.3.1
 * - header added
 *
 * Created 2015-02-19
*/

#ifndef BASE_H
#define BASE_H

#include <stdio.h>
#include <stdbool.h> // for C++ bool naming, requires C99
#include <stdlib.h> // for exit(), malloc
#include <math.h> // for fabs()
#include "randistrs.h" // For PRNGs
#include "global_param.h"

/*
* Data Type Declarations
*/

/*
* Function Declarations
*/

// Is point inside of boundary?
bool bPointInBoundary(const double point[3],
	const int boundary1Type,
	const double boundary1[]);
	
/* Do two sets of boundaries intersect?
*  Boundary Formats:
*  RECTANGLE 	-> [x_lower,x_upper,y_lower,y_upper]
*  CIRCLE    	-> [x_center,y_center,radius,radius^2]
*/
bool bBoundaryIntersect(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	const double clearance);

/* Are two sets of boundaries adjacent?
*  Both boundaries must be rectangular (either 2D or 3D)
*/
bool bBoundaryAdjacent(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	const double error,
	unsigned short * direction);

/* Is first boundary entirely inside the second?
*/
bool bBoundarySurround(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	const double clearance);

// Does a point lie within box created by two other points?
bool bPointBetween(const double p1[3],
	const double p2[3],
	const double newPoint[3]);

// Does a line segment intersect some boundary face? If so then which one and where?
// Returns the closest intersecting face from point p1 in positive direction along
// unit vector L
bool bLineHitBoundary(const double p1[3],
	const double L[3],
	const double length,
	const int boundary1Type,
	const double boundary1[],
	short * planeID,
	const bool bInside,
	double * d,
	double intersectPoint[3]);

// Does a line segment hit an infinite plane? If so then where?
bool bLineHitInfinitePlane(const double p1[3],
	const double L[3],
	const double length,
	const int boundary1Type,
	const double boundary1[],
	const short planeID,
	const bool bInside,
	double * d,
	double intersectPoint[3]);

// Is point that is in infinite plane also on boundary face?
// Assert that point is already on corresponding plane
bool bPointOnFace(const double p1[3],
	const int boundary1Type,
	const double boundary1[],
	const short planeID);

// Record specified face of boundary
void recordFace(const int boundary1Type,
	const double boundary1[],
	const short faceID,
	double boundaryFace[]);

// What is the value of the plane equation?
double planeEquation(const double point[3],
	const double plane[4]);

// Reflect point against a boundary.
// oldPoint is used to determine direction of reflection if needed
// bReflectInside indicates whether point should be reflected into boundary
// Returns false if point did not intersect boundary
bool reflectPoint(const double oldPoint[3],
	const double L[3],
	const double length,
	const double curPoint[3],
	double newPoint[3],
	double intersectPoint[3],
	short * planeID,
	const int boundary1Type,
	const double boundary1[],
	bool bReflectInside);

// "Push" a point along a line
void pushPoint(double p1[3],
	double p2[3],
	const double dist,
	const double L[3]);
	
// Determine distance from point to a boundary
double distanceToBoundary(const double point[3],
	const int boundary1Type,
	const double boundary1[]);

/* Determine boundary of intersection of two boundaries
*  Only valid for rectangular boundaries (rectangles or boxes) or spherical intersections.
*  Boundary Formats:
*  RECTANGLE 	-> [x_lower,x_upper,y_lower,y_upper]
*/
int intersectBoundary(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	double intersection[6]);

// Define unit vector pointing from one point to another
void defineLine(const double p1[3],
	const double p2[3],
	double L[3],
	double * length);

/* Determine area of boundary
*/
double boundaryArea(const int boundary1Type,
	const double boundary1[]);

/* Find a random coordinate within the specified range
*/
double uniformPoint(double rangeMin,
	double rangeMax);

// Find a random coordinate within the specified boundary
void uniformPointVolume(double point[3],
	const int boundaryType,
	const double boundary1[]);

/* Find distance between 2 3D points
*/
double pointDistance3D(const double point1[3],
	const double point2[3]);

/* Square a double value
*/
double squareDBL(double v);
	
#endif // BASE_H