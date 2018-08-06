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
 *
 * Last revised for AcCoRD v1.4 (2018-08-06)
 *
 * Revision history:
 *
 * Revision v1.4 (2018-08-06)
 * - added function to find closest point on a boundary from some other point
 *
 * Revision v1.1 (2016-12-24)
 * - renamed some switch cases to use direction macros instead of explicit integers
 *
 * Revision v0.7.0.1 (public beta, 2016-08-30)
 * - corrected calculation for spherical volume (was partially integer)
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - modified random number generation. Now use PCG via a separate interface file.
 * - added implementation of point shapes. Implementation is not comprehensive, but enough
 * to account for active point actors.
 * - added new version of function for the distance between 2 points where the coordinates
 * of one point are defined as separate variables
 * - added sphere as trivial case in function for finding a shape's closest face
 *
 * Revision v0.5.1 (2016-05-06)
 * - added 2D rectangle case to point reflection. Actually only works for surface cases,
 * since definition of faces are for the 3D case
 * - added closestFace and distanceToFace functions to find closest boundary face
 * to point (in direction of face normal only)
 * - updated check on a line hitting an infinite plane where acceptance of the distance = 0
 * case is passed as an argument
 *
 * Revision v0.5 (2016-04-15)
 * - filling in cases for 2D Rectangles
 * - added function to calculate boundary surface area. Renamed boundaryArea
 * function to boundaryVolume to avoid name confusion
 * - added function to return string of boundary name and integrated with error
 * messages
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
#include "rand_accord.h" // For PRNGs
#include "global_param.h"

//
// Data Type Declarations
//

//
// Function Declarations
//

// Is point inside of boundary?
bool bPointInBoundary(const double point[3],
	const int boundary1Type,
	const double boundary1[]);
	
// Do two sets of boundaries intersect?
bool bBoundaryIntersect(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	const double clearance);

// Are two sets of boundaries adjacent? Intersections will not be detected.
// Both boundaries must be rectangular (either 2D or 3D)
bool bBoundaryAdjacent(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	const double error,
	unsigned short * direction);

// Is first boundary entirely inside the second?
//
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
	const short planeIDConst,
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
	double intersectPoint[3],
	bool bZeroDist);

// Is point that is in infinite plane also on boundary face?
// Assert that point is already on corresponding plane
bool bPointOnFace(const double p1[3],
	const int boundary1Type,
	const double boundary1[],
	const short planeID);

// Do 2 boundaries share the same given surface?
// If so, faceShared specifies where they overlap
// This function is distinct from bBoundaryAdjacent because the shared
// face must be the same on both boundaries (e.g., lower x)
bool bSharedSurface(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	const short faceID,
	double faceShared[],
	const double error);

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

// Determine closest point on boundary to current point
// Unlike closestFace, uses actual straightline distance
void closestPoint(const double point[3],
	double newPoint[3],
	const int boundary1Type,
	const double boundary1[]);

// Determine closest boundary face from point
// Distance is checked along face normals only (i.e., we assume that we're already
// at one of the faces)
int closestFace(const double point[3],
	const int boundary1Type,
	const double boundary1[]);

// Determine distance from point to a boundary face
// Distance is along the face normal direction only
double distanceToFace(const double point[3],
	const int boundary1Type,
	const double boundary1[],
	const int faceID);

// Determine boundary of intersection of two boundaries
//  Only valid for rectangular boundaries (rectangles or boxes) or spherical intersections.
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

// Define unit vector pointing from one point to another
// This version defines 2nd point as individual arguments
// and does not return the length of the line
void defineLine2(const double p1[3],
	const double p2x,
	const double p2y,
	const double p2z,
	double L[3]);

// Determine volume of boundary
//
double boundaryVolume(const int boundary1Type,
	const double boundary1[]);

// Determine boundary surface Area
double boundarySurfaceArea(const int boundary1Type,
	const double boundary1[]);

// Find a random coordinate within the specified range
//
double uniformPoint(double rangeMin,
	double rangeMax);

// Find a random coordinate within the specified boundary
void uniformPointVolume(double point[3],
	const int boundaryType,
	const double boundary1[],
	bool bSurface,
	const short planeID);

// Find distance between 2 3D points
//
double pointDistance(const double point1[3],
	const double point2[3]);

// Square a double value
//
double squareDBL(double v);

// Return string with name of boundary
// Uses static memory strings in case output is not assigned
// to allocated memory
const char * boundaryString(const int boundaryType);
	
#endif // BASE_H