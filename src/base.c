/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * base.c - general utility functions that can apply to different simulation data
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
 * Revision v0.3.1.1
 * - added check in boundaryArea to identify invalid rectangle and box definitions. Will
 * 	 now return a zero value if boundary is invalid.
 *
 * Revision v0.3.1
 * - header added
 *
 * Created 2015-02-19
*/

#include "base.h" // for "Public" declarations

//
// "Private" Declarations
//

//
// Definitions
//

// Is point inside of boundary?
bool bPointInBoundary(const double point[3],
	const int boundary1Type,
	const double boundary1[])
{
	switch (boundary1Type)
	{
		case POINT:
			return point[0] == boundary1[0]
				&& point[1] == boundary1[1]
				&& point[2] == boundary1[2];
		case RECTANGLE:
		case RECTANGULAR_BOX:
			return (point[0] >= boundary1[0]
				&& point[0] <= boundary1[1]
				&& point[1] >= boundary1[2]
				&& point[1] <= boundary1[3]
				&& point[2] >= boundary1[4]
				&& point[2] <= boundary1[5]);
		case SPHERE:
			return (pointDistance(point, boundary1) <= boundary1[3]);
		default:
			fprintf(stderr,"ERROR: Cannot find point in shape type %s.\n", boundaryString(boundary1Type));
			return false;
	}
}

// Do two sets of boundaries overlap?
bool bBoundaryIntersect(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	const double clearance)
{
	double d;
	switch (boundary1Type)
	{
		case POINT:
			// Just need to check whether point is within boundary2
			return bPointInBoundary(boundary1, boundary2Type, boundary2);
		case RECTANGLE:
		case RECTANGULAR_BOX:
			switch (boundary2Type)
			{
				case RECTANGULAR_BOX:
					return (boundary1[2] < boundary2[3]
						&& boundary1[3] > boundary2[2]
						&& boundary1[0] < boundary2[1]
						&& boundary1[1] > boundary2[0]
						&& boundary1[4] < boundary2[5]
						&& boundary1[5] > boundary2[4]);
				case SPHERE:
					d = 0;
					if(boundary2[0] < boundary1[0])
						d += squareDBL(boundary2[0] - boundary1[0]);
					else if (boundary2[0] > boundary1[1])
						d += squareDBL(boundary2[0] - boundary1[1]);
					if(boundary2[1] < boundary1[2])
						d += squareDBL(boundary2[1] - boundary1[2]);
					else if (boundary2[1] > boundary1[3])
						d += squareDBL(boundary2[1] - boundary1[3]);
					if(boundary2[2] < boundary1[4])
						d += squareDBL(boundary2[2] - boundary1[4]);
					else if (boundary2[2] > boundary1[5])
						d += squareDBL(boundary2[2] - boundary1[5]);
					
					return (d < squareDBL(boundary2[3] + clearance) &&
						!bBoundarySurround(RECTANGULAR_BOX, boundary1,
						SPHERE, boundary2, 0.) &&
						!bBoundarySurround(SPHERE, boundary2,
						RECTANGULAR_BOX, boundary1, 0.) );
				default:
					fprintf(stderr,
						"ERROR: Cannot determine the intersection of a %s and a %s.\n",
						boundaryString(boundary2Type), boundaryString(boundary1Type));
					return false;
			}
		case SPHERE:
			switch (boundary2Type)
			{
				case SPHERE:
					d = pointDistance(boundary1, boundary2);
					return (d < boundary1[3] + boundary2[3] + clearance
						&& d > fabs(boundary1[3] - boundary2[3]));
				case RECTANGLE:
				case RECTANGULAR_BOX:
					d = 0;
					if(boundary1[0] < boundary2[0])
						d += squareDBL(boundary2[0] - boundary1[0]);
					else if (boundary1[0] > boundary2[1])
						d += squareDBL(boundary1[0] - boundary2[1]);
					if(boundary1[1] < boundary2[2])
						d += squareDBL(boundary1[1] - boundary2[2]);
					else if (boundary1[1] > boundary2[3])
						d += squareDBL(boundary1[1] - boundary2[3]);
					if(boundary1[2] < boundary2[4])
						d += squareDBL(boundary1[2] - boundary2[4]);
					else if (boundary1[2] > boundary2[5])
						d += squareDBL(boundary1[2] - boundary2[5]);
					
					return (d < squareDBL(boundary1[3] + clearance) &&
						!bBoundarySurround(RECTANGULAR_BOX, boundary2,
						SPHERE, boundary1, 0.) &&
						!bBoundarySurround(SPHERE, boundary1,
						RECTANGULAR_BOX, boundary2, 0.) );
				default:
					fprintf(stderr,
						"ERROR: Cannot determine the intersection of a %s and a %s.\n",
						boundaryString(boundary2Type), boundaryString(boundary1Type));
					return false;
			}
		default:
			fprintf(stderr,"ERROR: Cannot find intersection with shape %s.\n", boundaryString(boundary1Type));
			return false;
	}
}

// Are two sets of boundaries adjacent? Intersections will not be detected.
// Both boundaries must be rectangular (either 2D or 3D)
bool bBoundaryAdjacent(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	const double distError,
	unsigned short * direction)
{	

	if((boundary1Type == RECTANGULAR_BOX
		&& boundary2Type == RECTANGULAR_BOX) || 
		(boundary1Type == RECTANGLE
		&& boundary2Type == RECTANGULAR_BOX) || 
		(boundary1Type == RECTANGULAR_BOX
		&& boundary2Type == RECTANGLE))
	{
		if( // Do boxes share face along xy-plane?
			(boundary1[1] > boundary2[0] + distError) && (boundary2[1] > boundary1[0] + distError)
				&& (boundary1[3] > boundary2[2] + distError) && (boundary2[3] > boundary1[2] + distError))
		{
			if(fabs(boundary1[4] - boundary2[5]) < distError)
			{ // Boundary 2 is adjacent to boundary 1 along 1's lower z
				*direction = IN;
				return true;
			} else if (fabs(boundary2[4] - boundary1[5]) < distError)
			{ // Boundary 2 is adjacent to boundary 1 along 1's upper z
				*direction = OUT;
				return true;
			}
		} else if( // Do boxes share face along zy-plane?
			(boundary1[3] > boundary2[2] + distError) && (boundary2[3] > boundary1[2] + distError)
				&& (boundary1[5] > boundary2[4] + distError) && (boundary2[5] > boundary1[4] + distError))
		{
			if(fabs(boundary1[0] - boundary2[1]) < distError)
			{ // Boundary 2 is adjacent to boundary 1 along 1's lower x
				*direction = LEFT;
				return true;
			} else if (fabs(boundary2[0] - boundary1[1]) < distError)
			{ // Boundary 2 is adjacent to boundary 1 along 1's upper x
				*direction = RIGHT;
				return true;
			}
		} else if( // Do boxes share face along zx-plane?
			(boundary1[1] > boundary2[0] + distError) && (boundary2[1] > boundary1[0] + distError)
				&& (boundary1[5] > boundary2[4] + distError) && (boundary2[5] > boundary1[4] + distError))
		{
			if(fabs(boundary1[2] - boundary2[3]) < distError)
			{ // Boundary 2 is adjacent to boundary 1 along 1's lower y
				*direction = DOWN;
				return true;
			} else if (fabs(boundary2[2] - boundary1[3]) < distError)
			{ // Boundary 2 is adjacent to boundary 1 along 1's upper y
				*direction = UP;
				return true;
			}
		}
	} else if(boundary1Type == RECTANGLE
		&& boundary2Type == RECTANGLE)
	{ // Boundaries are both rectangles. They must lie in same plane to have adjacency
		if(boundary1[0] == boundary1[1]
			&& fabs(boundary1[0] - boundary2[0]) < distError
			&& fabs(boundary1[0] - boundary2[1]) < distError)
		{ // boundaries are both in YZ plane
			if ((boundary1[3] > boundary2[2] + distError)
				&& (boundary2[3] > boundary1[2] + distError))
			{ // There is overlap along Y
				if(fabs(boundary1[4] - boundary2[5]) < distError)
				{ // Boundary 2 is adjacent to boundary 1 along 1's lower z
					*direction = IN;
					return true;
				} else if (fabs(boundary2[4] - boundary1[5]) < distError)
				{ // Boundary 2 is adjacent to boundary 1 along 1's upper z
					*direction = OUT;
					return true;
				}
			} else if((boundary1[5] > boundary2[4] + distError)
				&& (boundary2[5] > boundary1[4] + distError))
			{ // There is overlap along Z
				if(fabs(boundary1[2] - boundary2[3]) < distError)
				{ // Boundary 2 is adjacent to boundary 1 along 1's lower y
					*direction = DOWN;
					return true;
				} else if (fabs(boundary2[2] - boundary1[3]) < distError)
				{ // Boundary 2 is adjacent to boundary 1 along 1's upper y
					*direction = UP;
					return true;
				}
			}
		} else if(boundary1[2] == boundary1[3]
			&& fabs(boundary1[2] - boundary2[2]) < distError
			&& fabs(boundary1[2] - boundary2[3]) < distError)
		{ // boundaries are both in XZ plane
			if((boundary1[1] > boundary2[0] + distError)
				&& (boundary2[1] > boundary1[0] + distError))
			{ // There is overlap along X
				if(fabs(boundary1[4] - boundary2[5]) < distError)
				{ // Boundary 2 is adjacent to boundary 1 along 1's lower z
					*direction = IN;
					return true;
				} else if (fabs(boundary2[4] - boundary1[5]) < distError)
				{ // Boundary 2 is adjacent to boundary 1 along 1's upper z
					*direction = OUT;
					return true;
				}
			} else if((boundary1[5] > boundary2[4] + distError)
				&& (boundary2[5] > boundary1[4] + distError))
			{ // There is overlap along Z
				if(fabs(boundary1[0] - boundary2[1]) < distError)
				{ // Boundary 2 is adjacent to boundary 1 along 1's lower x
					*direction = LEFT;
					return true;
				} else if (fabs(boundary2[0] - boundary1[1]) < distError)
				{ // Boundary 2 is adjacent to boundary 1 along 1's upper x
					*direction = RIGHT;
					return true;
				}
			}
		} else if(boundary1[4] == boundary1[5]
			&& fabs(boundary1[4] - boundary2[4]) < distError
			&& fabs(boundary1[4] - boundary2[5]) < distError)
		{ // boundaries are both in XY plane
			if((boundary1[1] > boundary2[0] + distError)
				&& (boundary2[1] > boundary1[0] + distError))
			{ // There is overlap along X
				if(fabs(boundary1[2] - boundary2[3]) < distError)
				{ // Boundary 2 is adjacent to boundary 1 along 1's lower y
					*direction = DOWN;
					return true;
				} else if (fabs(boundary2[2] - boundary1[3]) < distError)
				{ // Boundary 2 is adjacent to boundary 1 along 1's upper y
					*direction = UP;
					return true;
				}
			} else if ((boundary1[3] > boundary2[2] + distError)
				&& (boundary2[3] > boundary1[2] + distError))
			{ // There is overlap along Y
				if(fabs(boundary1[0] - boundary2[1]) < distError)
				{ // Boundary 2 is adjacent to boundary 1 along 1's lower x
					*direction = LEFT;
					return true;
				} else if (fabs(boundary2[0] - boundary1[1]) < distError)
				{ // Boundary 2 is adjacent to boundary 1 along 1's upper x
					*direction = RIGHT;
					return true;
				}
			}
		}
	}
	
	return false;
}

// Is first boundary entirely inside the second?
bool bBoundarySurround(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	const double clearance)
{
	double p1[3];
	double p2[3];
	
	switch (boundary1Type)
	{ // Is boundary1 inside of boundary2?
		case POINT:
			return bPointInBoundary(boundary1, boundary2Type, boundary2);
		case RECTANGLE:
		case RECTANGULAR_BOX:
			switch (boundary2Type)
			{
				case POINT:
					return false; // No object can be inside a point
				case RECTANGLE:
				case RECTANGULAR_BOX:
					return (boundary1[0] >= boundary2[0] + clearance
						&& boundary1[1] <= boundary2[1] - clearance
						&& boundary1[2] >= boundary2[2] + clearance
						&& boundary1[3] <= boundary2[3] - clearance
						&& boundary1[4] >= boundary2[4] + clearance
						&& boundary1[5] <= boundary2[5] - clearance);
				case SPHERE:				
					p1[0] = boundary1[0];
					p1[1] = boundary1[2];
					p1[2] = boundary1[4];
					if(boundary2[3] < pointDistance(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[0];
					p1[1] = boundary1[2];
					p1[2] = boundary1[5];
					if(boundary2[3] < pointDistance(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[0];
					p1[1] = boundary1[3];
					p1[2] = boundary1[4];
					if(boundary2[3] < pointDistance(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[0];
					p1[1] = boundary1[3];
					p1[2] = boundary1[5];
					if(boundary2[3] < pointDistance(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[1];
					p1[1] = boundary1[2];
					p1[2] = boundary1[4];
					if(boundary2[3] < pointDistance(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[1];
					p1[1] = boundary1[2];
					p1[2] = boundary1[5];
					if(boundary2[3] < pointDistance(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[1];
					p1[1] = boundary1[3];
					p1[2] = boundary1[4];
					if(boundary2[3] < pointDistance(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[1];
					p1[1] = boundary1[3];
					p1[2] = boundary1[5];
					if(boundary2[3] < pointDistance(p1, boundary2) + clearance)
						return false;
					// All fail cases have been tried
					return true;
				default:
					fprintf(stderr,
						"ERROR: Cannot determine whether a %s is inside of a %s.\n",
						boundaryString(boundary2Type), boundaryString(boundary1Type));
					return false;
			}
		case SPHERE:
			switch (boundary2Type)
			{
				case POINT:
					return false; // No object can be inside a point
				case RECTANGLE:
					return false; // A 3D object cannot be inside of a 2D object
				case RECTANGULAR_BOX:
					return(boundary1[3] <= (boundary1[0] - boundary2[0] - clearance) &&
						boundary1[3] <= (boundary2[1] - boundary1[0] - clearance) &&
						boundary1[3] <= (boundary1[1] - boundary2[2] - clearance) &&
						boundary1[3] <= (boundary2[3] - boundary1[1] - clearance) &&
						boundary1[3] <= (boundary1[2] - boundary2[4] - clearance) &&
						boundary1[3] <= (boundary2[5] - boundary1[2] - clearance));
				case SPHERE:
					return(boundary2[3] >= (boundary1[3] +
						pointDistance(boundary1, boundary2) + clearance));
				default:
					fprintf(stderr,
						"ERROR: Cannot determine whether a %s is inside of a %s.\n",
						boundaryString(boundary2Type), boundaryString(boundary1Type));
					return false;
			}
		default:
			fprintf(stderr,"ERROR: Cannot determine whether shape %s is inside another boundary.\n", boundaryString(boundary1Type));
			return false;
	}
}

// Does a point lie within box created by two other points?
bool bPointBetween(const double p1[3],
	const double p2[3],
	const double newPoint[3])
{
	int i;
	
	for(i = 0; i < 3; i++)
	{
		if(p1[i] > p2[i])
		{
			if(newPoint[i] < p2[i] || newPoint[i] > p1[i])
				return false;
		} else
		{
			if(newPoint[i] > p2[i] || newPoint[i] < p1[i])
				return false;
		}
	}
	return true;
}

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
	double intersectPoint[3])
{
	short curPlane;
	double minDist = INFINITY;
	double nearestIntersectPoint[3];
	bool bIntersect = false;
	
	switch(boundary1Type)
	{
		case RECTANGLE:
			if(bLineHitInfinitePlane(p1, L, length, RECTANGLE, boundary1,
				planeIDConst, false, d, intersectPoint, false)
				&& bPointOnFace(intersectPoint, RECTANGLE, boundary1, planeIDConst)
				&& *d < minDist)
			{
				return true;
			}
			return false;
		case RECTANGULAR_BOX:
			for(curPlane = 0; curPlane < 6; curPlane++)
			{
				if(bLineHitInfinitePlane(p1, L, length, RECTANGULAR_BOX, boundary1,
					curPlane, false, d, intersectPoint, false)
					&& bPointOnFace(intersectPoint, RECTANGULAR_BOX, boundary1, curPlane)
					&& *d < minDist)
				{ // Line does intersect this face at a valid distance and it is closest
					bIntersect = true;
					*planeID = curPlane;
					minDist = *d;
					nearestIntersectPoint[0] = intersectPoint[0];
					nearestIntersectPoint[1] = intersectPoint[1];
					nearestIntersectPoint[2] = intersectPoint[2];
				}
			}
			if(bIntersect)
			{
				*d = minDist;
				intersectPoint[0] = nearestIntersectPoint[0];
				intersectPoint[1] = nearestIntersectPoint[1];
				intersectPoint[2] = nearestIntersectPoint[2];
				return true;
			}
			return false;
		case SPHERE:
			return bLineHitInfinitePlane(p1, L, length, SPHERE, boundary1,
					curPlane, bInside, d, intersectPoint, false);
		default:
			fprintf(stderr,"ERROR: Cannot determine whether shape %s intersects another shape.\n", boundaryString(boundary1Type));
			return false;		
	}
}

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
	bool bZeroDist)
{	
	double centerToP1[3];
	double LDotCenterToP1;
	
	switch(boundary1Type)
	{
		case RECTANGLE:
		/*	switch(planeID)
			{
				case PLANE_XY:
					*d = (boundary1[4]-p1[2])/L[2];
					break;
				case PLANE_XZ:
					*d = (boundary1[2]-p1[1])/L[1];
					break;
				case PLANE_YZ:
					*d = (boundary1[0]-p1[0])/L[0];
					break;
				default:
					fprintf(stderr,"ERROR: Plane ID %d invalid for rectangle.\n", planeID);
					return false;	
			}
			intersectPoint[0] = (*d)*L[0] + p1[0];
			intersectPoint[1] = (*d)*L[1] + p1[1];
			intersectPoint[2] = (*d)*L[2] + p1[2];
			break;*/
		case RECTANGULAR_BOX:
			switch(planeID)
			{
				case LEFT:
					*d = (boundary1[0]-p1[0])/L[0];
					break;
				case RIGHT:
					*d = (boundary1[1]-p1[0])/L[0];
					break;
				case DOWN:
					*d = (boundary1[2]-p1[1])/L[1];
					break;
				case UP:
					*d = (boundary1[3]-p1[1])/L[1];
					break;
				case IN:
					*d = (boundary1[4]-p1[2])/L[2];
					break;
				case OUT:
					*d = (boundary1[5]-p1[2])/L[2];
					break;
				default:
					fprintf(stderr,"ERROR: Plane ID %d invalid for rectangular box.\n", planeID);
					return false;					
			}
			intersectPoint[0] = (*d)*L[0] + p1[0];
			intersectPoint[1] = (*d)*L[1] + p1[1];
			intersectPoint[2] = (*d)*L[2] + p1[2];
			break;
		case SPHERE:
			
			centerToP1[0] = p1[0] - boundary1[0];
			centerToP1[1] = p1[1] - boundary1[1];
			centerToP1[2] = p1[2] - boundary1[2];
			
			LDotCenterToP1 = L[0]*centerToP1[0] + L[1]*centerToP1[1] +
				L[2]*centerToP1[2];
				
			*d = sqrt(squareDBL(LDotCenterToP1) + squareDBL(boundary1[3])
				- squareDBL(centerToP1[0]) - squareDBL(centerToP1[1])
				- squareDBL(centerToP1[2]));
			
			if(bInside)
				*d = -LDotCenterToP1 + *d;
			else
				*d = -LDotCenterToP1 - *d;
			
			intersectPoint[0] = p1[0] + L[0]*(*d);
			intersectPoint[1] = p1[1] + L[1]*(*d);
			intersectPoint[2] = p1[2] + L[2]*(*d);
			break;		
		default:
			fprintf(stderr,"ERROR: Cannot determine whether a line hits the plane of a %s.\n", boundaryString(boundary1Type));
			*d = 0;
			return false;	
	}
	
	if(bZeroDist)
		return *d >= 0. && *d <= length;
	else
		return *d > 0. && *d <= length;
}

// Is point that is in infinite plane also on boundary face?
// Assert that point is already on corresponding plane
bool bPointOnFace(const double p1[3],
	const int boundary1Type,
	const double boundary1[],
	const short planeID)
{
	switch(boundary1Type)
	{
		case RECTANGLE:
		/*	switch(planeID)
			{
				case PLANE_XY:
					return (p1[1] >= boundary1[2]
						&& p1[1] <= boundary1[3]
						&& p1[0] >= boundary1[0]
						&& p1[0] <= boundary1[1]);
				case PLANE_XZ:
					return (p1[0] >= boundary1[0]
						&& p1[0] <= boundary1[1]
						&& p1[2] >= boundary1[4]
						&& p1[2] <= boundary1[5]);
				case PLANE_YZ:
					return (p1[1] >= boundary1[2]
						&& p1[1] <= boundary1[3]
						&& p1[2] >= boundary1[4]
						&& p1[2] <= boundary1[5]);					
			}*/
		case RECTANGULAR_BOX:
			switch(planeID)
			{
				case LEFT: // yz plane
				case RIGHT:
					return (p1[1] >= boundary1[2]
						&& p1[1] <= boundary1[3]
						&& p1[2] >= boundary1[4]
						&& p1[2] <= boundary1[5]);
				case DOWN: // xz plane
				case UP:
					return (p1[0] >= boundary1[0]
						&& p1[0] <= boundary1[1]
						&& p1[2] >= boundary1[4]
						&& p1[2] <= boundary1[5]);
				case IN: // xy plane
				case OUT:
					return (p1[1] >= boundary1[2]
						&& p1[1] <= boundary1[3]
						&& p1[0] >= boundary1[0]
						&& p1[0] <= boundary1[1]);
			}
		case SPHERE:
			// Trivially true
			return true;
		default:
			fprintf(stderr,"ERROR: Cannot determine whether point is on the face of a %s.\n", boundaryString(boundary1Type));
			return false;	
	}
}

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
	const double error)
{
	short i;
	short dim[2];
	switch (boundary1Type)
	{
		case RECTANGLE:
			switch (boundary2Type)
			{
				case RECTANGLE:
					// dim[0] will define the plane that the rectangles are on
					// dim[1] will define the shared varying coordinate
					
					// What plane are we on?
					if(boundary1[0] == boundary1[1])
						dim[0] = 0;
					else if(boundary1[2] == boundary1[3])
						dim[0] = 2;
					else if(boundary1[4] == boundary1[5])
						dim[0] = 4;
				
					// Is specified face normal to rectangle perimeter?
					// Are rectangles defined on the same plane?
					if(faceID == dim[0] || faceID == (dim[0] + 1) 
						|| boundary2[dim[0]] != boundary2[dim[0]+1])
						return false; // Shared face not possible
						
					switch(faceID)
					{
						case 0:
						case 1:
							if(dim[0] == 2) dim[1] = 4;
							if(dim[0] == 4) dim[1] = 2;
							break;
						case 2:
						case 3:
							if(dim[0] == 0) dim[1] = 4;
							if(dim[0] == 4) dim[1] = 0;
							break;
						case 4:
						case 5:
							if(dim[0] == 0) dim[1] = 2;
							if(dim[0] == 2) dim[1] = 0;
							break;
						default:
							fprintf(stderr,
								"ERROR: Specified face %d invalid for 2 Rectangles.\n", faceID);
							return false;							
					}
					
					// Is the line actually shared?
					if(fabs(boundary1[faceID] - boundary2[faceID]) > error)
						return false; // Lines are different
					
					if(boundary1[dim[1]] >= boundary2[dim[1]+1]
						|| boundary1[dim[1]+1] <= boundary2[dim[1]])
						return false; // The segments do not overlap
					
					// We have overlap. Determine shared segment
					for(i = 0; i < 6; i++)
						faceShared[i] = boundary1[i];
					
					if(boundary1[dim[1]] < boundary2[dim[1]])
						faceShared[dim[1]] = boundary2[dim[1]];
					else
						faceShared[dim[1]] = boundary1[dim[1]];
					
					if(boundary1[dim[1]+1] < boundary2[dim[1]+1])
						faceShared[dim[1]+1] = boundary1[dim[1]+1];
					else
						faceShared[dim[1]+1] = boundary2[dim[1]+1];
					
					return true;
				default:
					fprintf(stderr,
						"ERROR: Boundary types %s and %s are not allowed to share a surface.\n",
						boundaryString(boundary1Type), boundaryString(boundary2Type));
					return false;
			}
		case RECTANGULAR_BOX:
			switch (boundary2Type)
			{
				case RECTANGULAR_BOX:
					// dim[0] and dim[1] will define the plane that the shared surface
					// would be on (if it exists)
					switch(faceID)
					{
						case 0:
						case 1:
							dim[0] = 2;
							dim[1] = 4;
							break;
						case 2:
						case 3:
							dim[0] = 0;
							dim[1] = 4;
							break;
						case 4:
						case 5:
							dim[0] = 0;
							dim[1] = 2;
							break;
						default:
							fprintf(stderr,
								"ERROR: Specified face %d invalid for 2 Rectanglular Boxes.\n", faceID);
							return false;							
					}
					
					// Are the 2 faces on the same plane?
					if(fabs(boundary1[faceID] - boundary2[faceID]) > error)
						return false; // Planes are different
					
					// Do the 2 faces overlap?
					if(boundary1[dim[0]] >= boundary2[dim[0]+1]
						|| boundary1[dim[0]+1] <= boundary2[dim[0]]
						|| boundary1[dim[1]] >= boundary2[dim[1]+1]
						|| boundary1[dim[1]+1] <= boundary2[dim[1]])
						return false; // The faces do not overlap
					
					// We have overlap. Determine shared rectangle
					for(i = 0; i < 6; i++)
						faceShared[i] = boundary1[i];
					
					for(i = 0; i < 2; i++)
					{
						if(boundary1[dim[i]] < boundary2[dim[i]])
							faceShared[dim[i]] = boundary2[dim[i]];
						else
							faceShared[dim[i]] = boundary1[dim[i]];
						
						if(boundary1[dim[i]+1] < boundary2[dim[i]+1])
							faceShared[dim[i]+1] = boundary1[dim[i]+1];
						else
							faceShared[dim[i]+1] = boundary2[dim[i]+1];
					}
					
					return true;
				default:
					fprintf(stderr,
						"ERROR: Boundary types %s and %s are not allowed to share a surface.\n",
						boundaryString(boundary1Type), boundaryString(boundary2Type));
					return false;
			}
		case SPHERE: // Only one face on a sphere; no need to check faceID
			switch (boundary2Type)
			{
				case SPHERE:
					for(short i = 0; i < 3; i++)
					{
						if(boundary1[i] == boundary2[i])
							faceShared[i] = boundary1[i];
						else
							return false;
					}
					return true;
				default:
					fprintf(stderr,
						"ERROR: Boundary types %s and %s are not allowed to share a surface.\n",
						boundaryString(boundary1Type), boundaryString(boundary2Type));
					return false;
			}
		default:
			fprintf(stderr,"ERROR: Boundary type %s invalid to share a surface.\n", boundaryString(boundary1Type));
			return false;
	}
}

// Record specified face of boundary
void recordFace(const int boundary1Type,
	const double boundary1[],
	const short faceID,
	double boundaryFace[])
{
	
	switch(boundary1Type)
	{
		case RECTANGULAR_BOX:
		case RECTANGLE:
			switch(faceID)
			{
				case 0: // lower yz plane
					boundaryFace[0] = boundary1[0];
					boundaryFace[1] = boundary1[0];
					boundaryFace[2] = boundary1[2];
					boundaryFace[3] = boundary1[3];
					boundaryFace[4] = boundary1[4];
					boundaryFace[5] = boundary1[5];
					return;
				case 1: // upper yz plane
					boundaryFace[0] = boundary1[1];
					boundaryFace[1] = boundary1[1];
					boundaryFace[2] = boundary1[2];
					boundaryFace[3] = boundary1[3];
					boundaryFace[4] = boundary1[4];
					boundaryFace[5] = boundary1[5];
					return;
				case 2: // lower xz plane
					boundaryFace[0] = boundary1[0];
					boundaryFace[1] = boundary1[1];
					boundaryFace[2] = boundary1[2];
					boundaryFace[3] = boundary1[2];
					boundaryFace[4] = boundary1[4];
					boundaryFace[5] = boundary1[5];
					return;
				case 3: // upper xz plane
					boundaryFace[0] = boundary1[0];
					boundaryFace[1] = boundary1[1];
					boundaryFace[2] = boundary1[3];
					boundaryFace[3] = boundary1[3];
					boundaryFace[4] = boundary1[4];
					boundaryFace[5] = boundary1[5];
					return;
				case 4: // lower xy plane
					boundaryFace[0] = boundary1[0];
					boundaryFace[1] = boundary1[1];
					boundaryFace[2] = boundary1[2];
					boundaryFace[3] = boundary1[3];
					boundaryFace[4] = boundary1[4];
					boundaryFace[5] = boundary1[4];
					return;
				case 5: // lower xy plane
					boundaryFace[0] = boundary1[0];
					boundaryFace[1] = boundary1[1];
					boundaryFace[2] = boundary1[2];
					boundaryFace[3] = boundary1[3];
					boundaryFace[4] = boundary1[5];
					boundaryFace[5] = boundary1[5];
					return;
				default:
					fprintf(stderr,"ERROR: Face ID %d invalid for a %s.\n", faceID, boundaryString(boundary1Type));
					return;
			}
			return;
		case SPHERE:
			boundaryFace[0] = boundary1[0];
			boundaryFace[1] = boundary1[1];
			boundaryFace[2] = boundary1[2];
			boundaryFace[3] = boundary1[3];
			return;
		default:
			fprintf(stderr,"ERROR: Cannot record the face boundary of shape %s.\n", boundaryString(boundary1Type));
			return;	
	}
}


// What is the value of the plane equation for a given point?
double planeEquation(const double point[3],
	const double plane[4])
{
	return point[0]*plane[0] + point[1]*plane[1] + point[2]*plane[2] + plane[3];
}

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
	bool bReflectInside)
{
	double dist;			// Distance from oldPoint to boundary along line to curPoint
	double LDotCenterToOld;	// Dot product of L and centerToOld
	double pIntMinusC[3];
	double dDistance;	// Distance from intersectPoint to newPoint
	
	newPoint[0] = curPoint[0];
	newPoint[1] = curPoint[1];
	newPoint[2] = curPoint[2];
	
	if(!bLineHitBoundary(oldPoint, L, length, boundary1Type, boundary1,
		planeID, *planeID, bReflectInside, &dist, intersectPoint))
	{ // Line did not hit the boundary that it needs to reflect off of
		// We should just lock boundary closest to endPoint
		if(!bLineHitBoundary(oldPoint, L, INFINITY, boundary1Type, boundary1,
			planeID, *planeID, bReflectInside, &dist, intersectPoint))
		{ // Assume that point is already at boundary we want to reflect off of
		  // Just keep point at start
			intersectPoint[0] = oldPoint[0];
			intersectPoint[1] = oldPoint[1];
			intersectPoint[2] = oldPoint[2];
		} // Else line does eventually hit boundary. Just place point at that intersection
		newPoint[0] = intersectPoint[0];
		newPoint[1] = intersectPoint[1];
		newPoint[2] = intersectPoint[2];
		return false;
	}
	
	switch(boundary1Type)
	{
		case RECTANGLE:
		case RECTANGULAR_BOX:		
				
			switch(*planeID)
			{
				case LEFT:
					// Reflect off of lower x
					newPoint[0] = boundary1[0] + boundary1[0] - curPoint[0];
					return true;
				case RIGHT:
					// Reflect off of upper x
					newPoint[0] = boundary1[1] + boundary1[1] - curPoint[0];
					return true;
				case DOWN:
					// Reflect off of lower y
					newPoint[1] = boundary1[2] + boundary1[2] - curPoint[1];
					return true;
				case UP:
					// Reflect off of upper y
					newPoint[1] = boundary1[3] + boundary1[3] - curPoint[1];
					return true;
				case IN:
					// Reflect off of lower z
					newPoint[2] = boundary1[4] + boundary1[4] - curPoint[2];
					return true;
				case OUT:
					// Reflect off of upper z
					newPoint[2] = boundary1[5] + boundary1[5] - curPoint[2];
					return true;
				default:
					fprintf(stderr,"WARNING: Plane intersection ID %d invalid for a %s.\n", *planeID, boundaryString(boundary1Type));
					return false;					
			}
		case SPHERE:
			*planeID = 0; // There's only one surface on a sphere
			
			pIntMinusC[0] = intersectPoint[0] - boundary1[0];
			pIntMinusC[1] = intersectPoint[1] - boundary1[1];
			pIntMinusC[2] = intersectPoint[2] - boundary1[2];
			
			dDistance = 2*((curPoint[0]-intersectPoint[0])*pIntMinusC[0]
				+ (curPoint[1]-intersectPoint[1])*pIntMinusC[1]
				+ (curPoint[2]-intersectPoint[2])*pIntMinusC[2])
				/(squareDBL(pIntMinusC[0]) + squareDBL(pIntMinusC[1])
				+ squareDBL(pIntMinusC[2]));
			
			newPoint[0] -= dDistance*pIntMinusC[0];
			newPoint[1] -= dDistance*pIntMinusC[1];
			newPoint[2] -= dDistance*pIntMinusC[2];
			
			return true;
		default:
			fprintf(stderr,"ERROR: Cannot reflect a point off of a %s.\n", boundaryString(boundary1Type));
			return false;
	}
}

// "Push" a point along a line
void pushPoint(double p1[3],
	double p2[3],
	const double dist,
	const double L[3])
{
	p2[0] = p1[0] + dist*L[0];
	p2[1] = p1[1] + dist*L[1];
	p2[2] = p1[2] + dist*L[2];
}

// Determine distance from point to a boundary
double distanceToBoundary(const double point[3],
	const int boundary1Type,
	const double boundary1[])
{
	double dist = 0.;
	double dist2;
	
	switch (boundary1Type)
	{
		case RECTANGULAR_BOX:
			if(bPointInBoundary(point, boundary1Type, boundary1))
			{
				// Point is inside box find closest face
				dist = point[0] - boundary1[0];
				dist2 = boundary1[1] - point[0];
				if(dist2 < dist)
					dist = dist2;
				dist2 = point[1] - boundary1[2];
				if(dist2 < dist)
					dist = dist2;
				dist2 = boundary1[3] - point[1];
				if(dist2 < dist)
					dist = dist2;
				dist2 = point[2] - boundary1[4];
				if(dist2 < dist)
					dist = dist2;
				dist2 = boundary1[5] - point[2];
				if(dist2 < dist)
					dist = dist2;
				return dist;
			} else
			{ // Point is outside box
				if(point[0] < boundary1[0])
					dist += squareDBL(boundary1[0] - point[0]);
				else if(point[0] > boundary1[1])
					dist += squareDBL(boundary1[1] - point[0]);
				if(point[1] < boundary1[2])
					dist += squareDBL(boundary1[2] - point[1]);
				else if(point[1] > boundary1[3])
					dist += squareDBL(boundary1[3] - point[1]);
				if(point[2] < boundary1[4])
					dist += squareDBL(boundary1[4] - point[2]);
				else if(point[2] > boundary1[5])
					dist += squareDBL(boundary1[5] - point[2]);
				return sqrt(dist);
			}
			return 0.;
		case SPHERE:
			dist = pointDistance(point, boundary1) - boundary1[3];
			if(dist < 0)
				dist = -dist;
			return dist;
		default:
			fprintf(stderr,"ERROR: Cannot determine the distance from a point to a %s.\n", boundaryString(boundary1Type));
			return 0.;
	}
}

// Determine closest point on boundary to current point
// Unlike closestFace, uses actual straightline distance
void closestPoint(const double point[3],
	double newPoint[3],
	const int boundary1Type,
	const double boundary1[])
{
	double dist = 0.;
	double dist2;
	unsigned int face;
	double L[3];
	
	switch (boundary1Type)
	{
		case RECTANGULAR_BOX:
			if(bPointInBoundary(point, boundary1Type, boundary1))
			{
				// Point is inside box find closest face
				dist = point[0] - boundary1[0];
				dist2 = boundary1[1] - point[0];
				face = 0;
				newPoint[0] = point[0];
				newPoint[1] = point[1];
				newPoint[2] = point[2];
				if(dist2 < dist)
				{
					dist = dist2;
					face = 1;
				}
				dist2 = point[1] - boundary1[2];
				if(dist2 < dist)
				{
					dist = dist2;
					face = 2;
				}
				dist2 = boundary1[3] - point[1];
				if(dist2 < dist)
				{
					dist = dist2;
					face = 3;
				}
				dist2 = point[2] - boundary1[4];
				if(dist2 < dist)
				{
					dist = dist2;
					face = 4;
				}
				dist2 = boundary1[5] - point[2];
				if(dist2 < dist)
				{
					dist = dist2;
					face = 5;
				}
				switch(face)
				{
					case 0:
						newPoint[0] = boundary1[0];
						break;
					case 1:
						newPoint[0] = boundary1[1];
						break;
					case 2:
						newPoint[1] = boundary1[2];
						break;
					case 3:
						newPoint[1] = boundary1[3];
						break;
					case 4:
						newPoint[2] = boundary1[4];
						break;
					case 5:
						newPoint[2] = boundary1[5];
						break;
				}
			} else
			{ // Point is outside box
				if(point[0] < boundary1[0])
					newPoint[0] = boundary1[0];
				else if(point[0] > boundary1[1])
					newPoint[0] = boundary1[1];
				else
					newPoint[0] = point[0];
				if(point[1] < boundary1[2])
					newPoint[1] = boundary1[2];
				else if(point[1] > boundary1[3])
					newPoint[1] = boundary1[3];
				else
					newPoint[1] = point[1];					
				if(point[2] < boundary1[4])
					newPoint[2] = boundary1[4];
				else if(point[2] > boundary1[5])
					newPoint[2] = boundary1[5];
				else
					newPoint[2] = point[2];
			}
			return;
		case SPHERE:
			// Define line from point to center of sphere
			defineLine(point, boundary1, L, &dist);
			if(dist > boundary1[3])
			{ // Point is outside sphere
				bLineHitInfinitePlane(point, L, dist, SPHERE, boundary1, 0, false, &dist,
					newPoint, false);
			} else
			{ // Point is inside sphere
				bLineHitInfinitePlane(point, L, 2*boundary1[3], SPHERE, boundary1, 0, true, &dist,
					newPoint, false);
			}			
			return;
		default:
			fprintf(stderr,"ERROR: Cannot determine the distance from a point to a %s.\n", boundaryString(boundary1Type));
			return;
	}
}

// Determine closest boundary face from point
// Distance is checked along face normals only (i.e., we assume that we're already
// at one of the faces)
int closestFace(const double point[3],
	const int boundary1Type,
	const double boundary1[])
{
	int curFace;
	int minFace = UNDEFINED;
	double curDist;
	double minDist = INFINITY;
	
	switch(boundary1Type)
	{
		case RECTANGLE:
		case RECTANGULAR_BOX:
			for(curFace = 0; curFace < 6; curFace++)
			{
				curDist = distanceToFace(point, boundary1Type, boundary1, curFace);
				if(curDist < minDist)
				{
					minDist = curDist;
					minFace = curFace;
				}
			}
			return minFace;
		case SPHERE:
			// Sphere has only 1 face so it must be the closest
			return 0;
		default:
			fprintf(stderr,"ERROR: Cannot determine the distance from a %s.\n",
				boundaryString(boundary1Type));
			return UNDEFINED;
	}
}

// Determine distance from point to a boundary face
// Distance is along the face normal direction only
double distanceToFace(const double point[3],
	const int boundary1Type,
	const double boundary1[],
	const int faceID)
{
	double dist;
	
	switch(boundary1Type)
	{
		case RECTANGULAR_BOX:
			switch(faceID)
			{
				case LEFT:
					return fabs(point[0] - boundary1[0]);
				case RIGHT:
					return fabs(point[0] - boundary1[1]);
				case DOWN:
					return fabs(point[1] - boundary1[2]);
				case UP:
					return fabs(point[1] - boundary1[3]);
				case IN:
					return fabs(point[2] - boundary1[4]);
				case OUT:
					return fabs(point[2] - boundary1[5]);
				default:
				fprintf(stderr,"ERROR: Cannot determine the distance to face %d of a %s.\n",
					faceID, boundaryString(boundary1Type));
				return INFINITY;
			}
			break;
		case SPHERE:
			dist = pointDistance(point, boundary1);
			if(dist > boundary1[3])
				return dist - boundary1[3];
			else
				return dist;
		default:
			fprintf(stderr,"ERROR: Cannot determine the distance from a %s.\n",
				boundaryString(boundary1Type));
			return INFINITY;
	}
}

// Determine boundary of intersection of two boundaries
//  Only valid for rectangular boundaries (rectangles or boxes) or spherical intersections.
int intersectBoundary(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	double intersection[6])
{
	
	if((boundary1Type == RECTANGULAR_BOX || boundary1Type == RECTANGLE)
		&&	(boundary2Type == RECTANGULAR_BOX || boundary2Type == RECTANGLE))
	{
		intersection[0] = (boundary1[0] > boundary2[0])? boundary1[0] : boundary2[0];
		intersection[1] = (boundary1[1] < boundary2[1])? boundary1[1] : boundary2[1];
		intersection[2] = (boundary1[2] > boundary2[2])? boundary1[2] : boundary2[2];
		intersection[3] = (boundary1[3] < boundary2[3])? boundary1[3] : boundary2[3];
		intersection[4] = (boundary1[4] > boundary2[4])? boundary1[4] : boundary2[4];
		intersection[5] = (boundary1[5] < boundary2[5])? boundary1[5] : boundary2[5];
		if(boundary1Type == RECTANGLE && boundary2Type == RECTANGLE)
			return RECTANGLE;
		else
			return RECTANGULAR_BOX;
	} else if(boundary1Type == SPHERE || boundary2Type == SPHERE)
	{
		// At least one of the boundaries is a sphere. One boundary must be
		// contained fully within the other boundary
		if(bBoundarySurround(boundary1Type, boundary1, boundary2Type, boundary2, 0.))
		{
			// boundary 1 is within boundary 2
			intersection[0] = boundary1[0];
			intersection[1] = boundary1[1];
			intersection[2] = boundary1[2];
			intersection[3] = boundary1[3];
			intersection[4] = boundary1[4];
			intersection[5] = boundary1[5];
			return boundary1Type;
		} else if(bBoundarySurround(boundary2Type, boundary2, boundary1Type, boundary1, 0.))
		{
			// boundary 2 is within boundary 1
			intersection[0] = boundary2[0];
			intersection[1] = boundary2[1];
			intersection[2] = boundary2[2];
			intersection[3] = boundary2[3];
			intersection[4] = boundary2[4];
			intersection[5] = boundary2[5];
			return boundary2Type;
		} else if(!bBoundaryIntersect(boundary2Type, boundary2, boundary1Type, boundary1, 0.))
		{
			// Boundaries do not intersect at all
			intersection[0] = 0.;			
			intersection[1] = 0.;
			intersection[2] = 0.;
			intersection[3] = 0.;
			intersection[4] = 0.;
			intersection[5] = 0.;
			return RECTANGULAR_BOX;
		} else
		{
			// Intersection is invalid
			fprintf(stderr,"ERROR: Intersection of two boundaries is invalid. At least one boundary is spherical and hits the other boundary.\n");
			return UNDEFINED_SHAPE;
		}
	} else
	{	// Intersection for combination of boundary types is unknown
		fprintf(stderr,"ERROR: Cannot determine the intersection of a %s and a %s.\n", boundaryString(boundary2Type), boundaryString(boundary1Type));
		return UNDEFINED_SHAPE;
	}
	
}

// Define unit vector pointing from one point to another
void defineLine(const double p1[3],
	const double p2[3],
	double L[3],
	double * length)
{
	*length = sqrt(squareDBL(p2[0]-p1[0]) + squareDBL(p2[1]-p1[1]) + squareDBL(p2[2]-p1[2]));
	
	if (*length > 0.)
	{
		L[0] = (p2[0]-p1[0])/(*length);
		L[1] = (p2[1]-p1[1])/(*length);
		L[2] = (p2[2]-p1[2])/(*length);
	} else
	{
		L[0] = 0.;
		L[1] = 0.;
		L[2] = 0.;
		*length = 0.;
	}
}

// Define unit vector pointing from one point to another
// This version defines 2nd point as individual arguments
// and does not return the length of the line
void defineLine2(const double p1[3],
	const double p2x,
	const double p2y,
	const double p2z,
	double L[3])
{
	double length = sqrt(squareDBL(p2x-p1[0]) + squareDBL(p2y-p1[1]) + squareDBL(p2z-p1[2]));
	
	if (length > 0.)
	{
		L[0] = (p2x-p1[0])/length;
		L[1] = (p2y-p1[1])/length;
		L[2] = (p2z-p1[2])/length;
	} else
	{
		L[0] = 0.;
		L[1] = 0.;
		L[2] = 0.;
	}
}

// Determine volume of boundary
double boundaryVolume(const int boundary1Type,
	const double boundary1[])
{
	switch (boundary1Type)
	{
		case POINT:
			return 0.;
		case RECTANGLE:
			if(boundary1[1] < boundary1[0] || boundary1[3] < boundary1[2]
				 || boundary1[5] < boundary1[4])
				return 0.;
			else
				if(boundary1[0] == boundary1[1])
					return (boundary1[5]-boundary1[4])*(boundary1[3]-boundary1[2]);
				if(boundary1[2] == boundary1[3])
					return (boundary1[1]-boundary1[0])*(boundary1[5]-boundary1[4]);
				if(boundary1[4] == boundary1[5])
					return (boundary1[1]-boundary1[0])*(boundary1[3]-boundary1[2]);
		case RECTANGULAR_BOX:
			if(boundary1[1] < boundary1[0]
				|| boundary1[3] < boundary1[2]
				|| boundary1[5] < boundary1[4])
				return 0.;
			else
				return (boundary1[1]-boundary1[0])*(boundary1[3]-boundary1[2])
					*(boundary1[5]-boundary1[4]);
		case CIRCLE:
			return PI*squareDBL(boundary1[3]);
		case SPHERE:
			return 4./3.*PI*boundary1[3]*boundary1[3]*boundary1[3];
		case LINE:
			return sqrt(squareDBL(boundary1[1]-boundary1[0])
				+ squareDBL(boundary1[3]-boundary1[2])
				+ squareDBL(boundary1[5]-boundary1[4]));
		default:
			fprintf(stderr,"ERROR: Cannot determine the volume of a %s.\n", boundaryString(boundary1Type));
			return 0;
	}
}

// Determine boundary surface Area
double boundarySurfaceArea(const int boundary1Type,
	const double boundary1[])
{
	double area = 0.;
	
	switch (boundary1Type)
	{
		case RECTANGLE:
			if(boundary1[1] < boundary1[0] || boundary1[3] < boundary1[2]
				 || boundary1[5] < boundary1[4])
				return 0.;
			
			area += 2*(boundary1[1]-boundary1[0]);
			area += 2*(boundary1[3]-boundary1[2]);
			area += 2*(boundary1[5]-boundary1[4]);
			return area;
			
		case RECTANGULAR_BOX:
			if(boundary1[1] < boundary1[0]
				|| boundary1[3] < boundary1[2]
				|| boundary1[5] < boundary1[4])
				return 0.;
			
			area += 2* (boundary1[1]-boundary1[0]) * (boundary1[3]-boundary1[2]);
			area += 2* (boundary1[1]-boundary1[0]) * (boundary1[5]-boundary1[4]);
			area += 2* (boundary1[5]-boundary1[4]) * (boundary1[3]-boundary1[2]);
			return area;
			
		case CIRCLE:
			return 2*PI*boundary1[3];
		case SPHERE:
			return 4*PI*boundary1[3]*boundary1[3];
		default:
			fprintf(stderr,"ERROR: Boundary type %s invalid.\n", boundaryString(boundary1Type));
			return 0;
	}
}

// Find a random coordinate within the specified range
double uniformPoint(double rangeMin,
	double rangeMax)
{
	return (rangeMin + (rangeMax-rangeMin)*generateUniform());
}

// Find a random coordinate within the specified boundary
void uniformPointVolume(double point[3],
	const int boundaryType,
	const double boundary1[],
	bool bSurface,
	const short planeID)
{
	bool bNeedPoint;
	short curFace;
	double r, rSq;
	
	switch (boundaryType)
	{
		case POINT:
			// Simplest case. No need to generate a random point.
			// Assume that boundary still defined in "Rectangular" format
			point[0] = boundary1[0];
			point[1] = boundary1[2];
			point[2] = boundary1[4];
			return;
		case RECTANGLE:
			if(bSurface)
			{
				curFace = (short) floor(4*generateUniform());
				switch(planeID)
				{
					case PLANE_XY:
						point[2] = boundary1[4];
						switch(curFace)
						{
							case 0:
							case 1:
								point[0] = boundary1[curFace];
								point[1] = uniformPoint(boundary1[2], boundary1[3]);
								break;
							case 2:
							case 3:
								point[0] = uniformPoint(boundary1[0], boundary1[1]);
								point[1] = boundary1[curFace];
								break;
						}
						break;
					case PLANE_XZ:
						point[1] = boundary1[2];
						switch(curFace)
						{
							case 0:
							case 1:
								point[0] = boundary1[curFace];
								point[2] = uniformPoint(boundary1[4], boundary1[5]);
								break;
							case 2:
							case 3:
								point[0] = uniformPoint(boundary1[0], boundary1[1]);
								point[1] = boundary1[curFace+2];
								break;
						}
						break;
					case PLANE_YZ:
						point[0] = boundary1[0];
						switch(curFace)
						{
							case 0:
							case 1:
								point[2] = boundary1[curFace+4];
								point[1] = uniformPoint(boundary1[2], boundary1[3]);
								break;
							case 2:
							case 3:
								point[2] = uniformPoint(boundary1[4], boundary1[5]);
								point[1] = boundary1[curFace];
								break;
						}
						break;
					default:
						// Something went wrong
						fprintf(stderr,"ERROR: Cannot generate a uniform random point on plane %u of a rectangle.\n", planeID);
						return;
				}
				return;
			}
			switch(planeID)
			{
				case PLANE_XY:
					point[0] = uniformPoint(boundary1[0], boundary1[1]);
					point[1] = uniformPoint(boundary1[2], boundary1[3]);
					point[2] = boundary1[4];
					break;
				case PLANE_XZ:
					point[0] = uniformPoint(boundary1[0], boundary1[1]);
					point[1] = boundary1[2];	
					point[2] = uniformPoint(boundary1[4], boundary1[5]);
					break;
				case PLANE_YZ:
					point[0] = boundary1[0];
					point[1] = uniformPoint(boundary1[2], boundary1[3]);
					point[2] = uniformPoint(boundary1[4], boundary1[5]);
					break;
				default:
					// Something went wrong
					fprintf(stderr,"ERROR: Cannot generate a uniform random point on plane %u of a rectangle.\n", planeID);
					return;
			}
		case RECTANGULAR_BOX:
			if(bSurface)
			{
				curFace = (short) floor(6*generateUniform());
				switch(curFace)
				{
					case 0:
					case 1:
						point[0] = boundary1[curFace];
						point[1] = uniformPoint(boundary1[2], boundary1[3]);
						point[2] = uniformPoint(boundary1[4], boundary1[5]);
						break;
					case 2:
					case 3:
						point[0] = uniformPoint(boundary1[0], boundary1[1]);
						point[1] = boundary1[curFace];
						point[2] = uniformPoint(boundary1[4], boundary1[5]);
						break;
					case 4:
					case 5:
						point[0] = uniformPoint(boundary1[0], boundary1[1]);
						point[1] = uniformPoint(boundary1[2], boundary1[3]);
						point[2] = boundary1[curFace];
						break;
				}
				return;
			}
			point[0] = uniformPoint(boundary1[0], boundary1[1]);
			point[1] = uniformPoint(boundary1[2], boundary1[3]);
			point[2] = uniformPoint(boundary1[4], boundary1[5]);
			return;
		case CIRCLE:
			return;
		case SPHERE:			
			// Use rejection method to create point in sphere
			bNeedPoint = true;
			while(bNeedPoint)
			{
				point[0] = generateUniform();
				point[1] = generateUniform();
				point[2] = generateUniform();
				
				rSq = squareDBL(point[0]) + squareDBL(point[1])
					+ squareDBL(point[2]);
				
				if (rSq < 1.)
				{
					// Found valid point. Scale as needed and randomize sign
					if(generateUniform() > 0.5)
						point[0] = -point[0];
					if(generateUniform() > 0.5)
						point[1] = -point[1];
					if(generateUniform() > 0.5)
						point[2] = -point[2];
					
					if(bSurface)
					{
						r = sqrt(rSq);
						point[0] /= r;
						point[1] /= r;
						point[2] /= r;
					}					
					point[0] = boundary1[0] + point[0]*boundary1[3];
					point[1] = boundary1[1] + point[1]*boundary1[3];
					point[2] = boundary1[2] + point[2]*boundary1[3];
					
					bNeedPoint = false;
				}
			}
			return;
		default:
			fprintf(stderr,"ERROR: Cannot generate a uniform random point in a %s.\n", boundaryString(boundaryType));
			return;
	}
}

// Find distance between 2 3D points
double pointDistance(const double point1[3],
	const double point2[3])
{
	return sqrt(squareDBL(point2[0] - point1[0]) +
		squareDBL(point2[1] - point1[1]) +
		squareDBL(point2[2] - point1[2]));
}

// Square a double value
double squareDBL(double v)
{
	return v*v;
}

// Return string with name of boundary
// Uses static memory strings in case output is not assigned
// to allocated memory
const char * boundaryString(const int boundaryType)
{
	static char pointString[] = "Point";
	static char rectString[] = "Rectangle";
	static char boxString[] = "Rectangular Box";
	static char circleString[] = "Circle";
	static char sphereString[] = "Sphere";
	static char emptyString[] = "";
	
	switch (boundaryType)
	{
		case POINT:
			return pointString;
		case RECTANGLE:
			return rectString;
		case RECTANGULAR_BOX:
			return boxString;
		case CIRCLE:
			return circleString;
		case SPHERE:
			return sphereString;
		default:
			fprintf(stderr,"ERROR: Shape type %d does not have an associated name.\n", boundaryType);
			return emptyString;
	}
}