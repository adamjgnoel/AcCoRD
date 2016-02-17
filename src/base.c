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
 * Last revised for AcCoRD v0.5
 *
 * Revision history:
 *
 * Revision v0.5
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
		case RECTANGLE:
			return (point[0] >= boundary1[0]
				&& point[0] <= boundary1[1]
				&& point[1] >= boundary1[2]
				&& point[1] <= boundary1[3]
				&& point[2] >= boundary1[4]
				&& point[2] <= boundary1[5]);
		case RECTANGULAR_BOX:
			return (point[0] >= boundary1[0]
				&& point[0] <= boundary1[1]
				&& point[1] >= boundary1[2]
				&& point[1] <= boundary1[3]
				&& point[2] >= boundary1[4]
				&& point[2] <= boundary1[5]);
		case SPHERE:
			return (pointDistance3D(point, boundary1) < boundary1[3]);
		default:
			fprintf(stderr,"ERROR: Boundary type %d invalid.\n", boundary1Type);
			return false;
	}
}

//Do two sets of boundaries overlap?
//  Boundary Formats:
//  RECTANGLE 	-> [x_lower,x_upper,y_lower,y_upper]
//  CIRCLE    	-> [x_center,y_center,radius,radius^2]
//
bool bBoundaryIntersect(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	const double clearance)
{
	double d;
	switch (boundary1Type)
	{
		case RECTANGLE:
			switch (boundary2Type)
			{
				case RECTANGLE:
					return (boundary1[2] < boundary2[3]
						&& boundary1[3] > boundary2[2]
						&& boundary1[0] < boundary2[1]
						&& boundary1[1] > boundary2[0]);
				default:
					fprintf(stderr,
						"ERROR: Boundary type combination %d and %d invalid.\n",
						boundary1Type, boundary2Type);
					return false;
			}
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
						"ERROR: Boundary type combination %d and %d invalid.\n",
						boundary1Type, boundary2Type);
					return false;
			}
		case SPHERE:
			switch (boundary2Type)
			{
				case SPHERE:
					d = pointDistance3D(boundary1, boundary2);
					return (d < boundary1[3] + boundary2[3] + clearance
						&& d > fabs(boundary1[3] - boundary2[3]));
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
						"ERROR: Boundary type combination %d and %d invalid.\n",
						boundary1Type, boundary2Type);
					return false;
			}
		default:
			fprintf(stderr,"ERROR: Boundary type %d invalid.\n", boundary1Type);
			return false;
	}
}

// Are two sets of boundaries adjacent?
//  Both boundaries must be rectangular (either 2D or 3D)
bool bBoundaryAdjacent(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	const double distError,
	unsigned short * direction)
{	
	if(boundary1Type == RECTANGLE && boundary2Type == RECTANGLE)
	{
		if( // Do rectangles share face along y-axis?
			(boundary1[3] > boundary2[2] + distError) && (boundary2[3] > boundary1[2] + distError)
		)
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
		} else if( // Do rectangles share face along x-axis?
			(boundary1[1] > boundary2[0] + distError) && (boundary2[1] > boundary1[0] + distError)
		)
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
	} else if(boundary1Type == RECTANGULAR_BOX && boundary2Type == RECTANGULAR_BOX)
	{
		if( // Do rectangles share face along xy-plane?
			(boundary1[1] > boundary2[0] + distError) && (boundary2[1] > boundary1[0] + distError
				&& boundary1[3] > boundary2[2] + distError) && (boundary2[3] > boundary1[2] + distError)
		)
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
		} else if( // Do rectangles share face along zy-plane?
			(boundary1[3] > boundary2[2] + distError) && (boundary2[3] > boundary1[2] + distError
				&& boundary1[5] > boundary2[4] + distError) && (boundary2[5] > boundary1[4] + distError)
		)
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
		} else if( // Do rectangles share face along zx-plane?
			(boundary1[1] > boundary2[0] + distError) && (boundary2[1] > boundary1[0] + distError
				&& boundary1[5] > boundary2[4] + distError) && (boundary2[5] > boundary1[4] + distError)
		)
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
		case RECTANGLE:
			switch (boundary2Type)
			{
				case RECTANGLE:
					// TODO: Clearance doesn't apply to 1 dimension of 2D shapes
					return (boundary1[0] >= boundary2[0] + clearance
						&& boundary1[1] <= boundary2[1] - clearance
						&& boundary1[2] >= boundary2[2] + clearance
						&& boundary1[3] <= boundary2[3] - clearance
						&& boundary1[4] >= boundary2[4] + clearance
						&& boundary1[5] <= boundary2[5] - clearance);
				default:
					fprintf(stderr,
						"ERROR: Boundary type combination %d and %d invalid.\n",
						boundary1Type, boundary2Type);
					return false;
			}
		case RECTANGULAR_BOX:
			switch (boundary2Type)
			{
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
					if(boundary2[3] < pointDistance3D(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[0];
					p1[1] = boundary1[2];
					p1[2] = boundary1[5];
					if(boundary2[3] < pointDistance3D(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[0];
					p1[1] = boundary1[3];
					p1[2] = boundary1[4];
					if(boundary2[3] < pointDistance3D(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[0];
					p1[1] = boundary1[3];
					p1[2] = boundary1[5];
					if(boundary2[3] < pointDistance3D(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[1];
					p1[1] = boundary1[2];
					p1[2] = boundary1[4];
					if(boundary2[3] < pointDistance3D(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[1];
					p1[1] = boundary1[2];
					p1[2] = boundary1[5];
					if(boundary2[3] < pointDistance3D(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[1];
					p1[1] = boundary1[3];
					p1[2] = boundary1[4];
					if(boundary2[3] < pointDistance3D(p1, boundary2) + clearance)
						return false;
					p1[0] = boundary1[1];
					p1[1] = boundary1[3];
					p1[2] = boundary1[5];
					if(boundary2[3] < pointDistance3D(p1, boundary2) + clearance)
						return false;
					// All fail cases have been tried
					return true;
				default:
					fprintf(stderr,
						"ERROR: Boundary type combination %d and %d invalid.\n",
						boundary1Type, boundary2Type);
					return false;
			}
		case SPHERE:
			switch (boundary2Type)
			{
				case RECTANGULAR_BOX:
					return(boundary1[3] <= (boundary1[0] - boundary2[0] - clearance) &&
						boundary1[3] <= (boundary2[1] - boundary1[0] - clearance) &&
						boundary1[3] <= (boundary1[1] - boundary2[2] - clearance) &&
						boundary1[3] <= (boundary2[3] - boundary1[1] - clearance) &&
						boundary1[3] <= (boundary1[2] - boundary2[4] - clearance) &&
						boundary1[3] <= (boundary2[5] - boundary1[2] - clearance));
				case SPHERE:
					return(boundary2[3] >= (boundary1[3] +
						pointDistance3D(boundary1, boundary2) + clearance));
				default:
					fprintf(stderr,
						"ERROR: Boundary type combination %d and %d invalid.\n",
						boundary1Type, boundary2Type);
					return false;
			}
		default:
			fprintf(stderr,"ERROR: Boundary type %d invalid.\n", boundary1Type);
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
		case RECTANGULAR_BOX:
			for(curPlane = 0; curPlane < 6; curPlane++)
			{
				if(bLineHitInfinitePlane(p1, L, length, RECTANGULAR_BOX, boundary1,
					curPlane, false, d, intersectPoint)
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
					curPlane, bInside, d, intersectPoint);
		default:
			fprintf(stderr,"ERROR: Boundary type %d invalid for boundary intersection.\n", boundary1Type);
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
	double intersectPoint[3])
{	
	double centerToP1[3];
	double LDotCenterToP1;
	
	switch(boundary1Type)
	{
		case RECTANGULAR_BOX:
			switch(planeID)
			{
				case 0:
					*d = (boundary1[0]-p1[0])/L[0];
					break;
				case 1:
					*d = (boundary1[1]-p1[0])/L[0];
					break;
				case 2:
					*d = (boundary1[2]-p1[1])/L[1];
					break;
				case 3:
					*d = (boundary1[3]-p1[1])/L[1];
					break;
				case 4:
					*d = (boundary1[4]-p1[2])/L[2];
					break;
				case 5:
					*d = (boundary1[5]-p1[2])/L[2];
					break;
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
			fprintf(stderr,"ERROR: Boundary type %d invalid for single plane intersection.\n", boundary1Type);
			*d = 0;
			return false;	
	}
			
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
		case RECTANGULAR_BOX:
			switch(planeID)
			{
				case 0: // yz plane
				case 1:
					return (p1[1] >= boundary1[2]
						&& p1[1] <= boundary1[3]
						&& p1[2] >= boundary1[4]
						&& p1[2] <= boundary1[5]);
				case 2: // xz plane
				case 3:
					return (p1[0] >= boundary1[0]
						&& p1[0] <= boundary1[1]
						&& p1[2] >= boundary1[4]
						&& p1[2] <= boundary1[5]);
				case 4: // xy plane
				case 5:
					return (p1[1] >= boundary1[2]
						&& p1[1] <= boundary1[3]
						&& p1[0] >= boundary1[0]
						&& p1[0] <= boundary1[1]);
			}
		case SPHERE:
			// Trivially true
			return true;
		default:
			fprintf(stderr,"ERROR: Boundary type %d invalid.\n", boundary1Type);
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
					fprintf(stderr,"ERROR: Face ID %d invalid.\n", faceID);
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
			fprintf(stderr,"ERROR: Boundary type %d invalid.\n", boundary1Type);
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
		planeID, bReflectInside, &dist, intersectPoint))
	{ // Line did not hit the boundary that it needs to reflect off of
		// We should just lock boundary closest to endPoint
		if(!bLineHitBoundary(oldPoint, L, INFINITY, boundary1Type, boundary1,
			planeID, bReflectInside, &dist, intersectPoint))
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
		case RECTANGULAR_BOX:		
				
			switch(*planeID)
			{
				case 0:
					// Reflect off of lower x
					newPoint[0] = boundary1[0] + boundary1[0] - curPoint[0];
					return true;
				case 1:
					// Reflect off of upper x
					newPoint[0] = boundary1[1] + boundary1[1] - curPoint[0];
					return true;
				case 2:
					// Reflect off of lower y
					newPoint[1] = boundary1[2] + boundary1[2] - curPoint[1];
					return true;
				case 3:
					// Reflect off of upper y
					newPoint[1] = boundary1[3] + boundary1[3] - curPoint[1];
					return true;
				case 4:
					// Reflect off of lower z
					newPoint[2] = boundary1[4] + boundary1[4] - curPoint[2];
					return true;
				case 5:
					// Reflect off of upper z
					newPoint[2] = boundary1[5] + boundary1[5] - curPoint[2];
					return true;
				default:
					fprintf(stderr,"WARNING: Plane intersection ID %d invalid.\n", *planeID);
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
			fprintf(stderr,"ERROR: Boundary type %d invalid.\n", boundary1Type);
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
			dist = pointDistance3D(point, boundary1) - boundary1[3];
			if(dist < 0)
				dist = -dist;
			return dist;
		default:
			fprintf(stderr,"ERROR: Boundary type %d invalid.\n", boundary1Type);
			return 0.;
	}
}

// Determine boundary of intersection of two boundaries
//  Only valid for rectangular boundaries (rectangles or boxes) or spherical intersections.
//  Boundary Formats:
//  RECTANGLE 	-> [x_lower,x_upper,y_lower,y_upper]
//
int intersectBoundary(const int boundary1Type,
	const double boundary1[],
	const int boundary2Type,
	const double boundary2[],
	double intersection[6])
{
	
	if(boundary1Type == RECTANGLE && boundary2Type == RECTANGLE)
	{
		intersection[0] = (boundary1[0] > boundary2[0])? boundary1[0] : boundary2[0];
		intersection[1] = (boundary1[1] < boundary2[1])? boundary1[1] : boundary2[1];
		intersection[2] = (boundary1[2] > boundary2[2])? boundary1[2] : boundary2[2];
		intersection[3] = (boundary1[3] < boundary2[3])? boundary1[3] : boundary2[3];
		intersection[4] = 0.;
		intersection[5] = 0.;
		return RECTANGLE;
	} else if(boundary1Type == RECTANGULAR_BOX && boundary2Type == RECTANGULAR_BOX)
	{
		intersection[0] = (boundary1[0] > boundary2[0])? boundary1[0] : boundary2[0];
		intersection[1] = (boundary1[1] < boundary2[1])? boundary1[1] : boundary2[1];
		intersection[2] = (boundary1[2] > boundary2[2])? boundary1[2] : boundary2[2];
		intersection[3] = (boundary1[3] < boundary2[3])? boundary1[3] : boundary2[3];
		intersection[4] = (boundary1[4] > boundary2[4])? boundary1[4] : boundary2[4];
		intersection[5] = (boundary1[5] < boundary2[5])? boundary1[5] : boundary2[5];
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
			fprintf(stderr,"ERROR: Error: Intersection of two boundaries is invalid. At least one boundary is spherical and hits the other boundary.\n");
			return UNDEFINED_SHAPE;
		}
	} else
	{	// Intersection for combination of boundary types is unknown
		fprintf(stderr,"ERROR: Error: Intersection between Boundary type %d and %d unknown.\n", boundary1Type, boundary2Type);
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

// Determine area of boundary
double boundaryArea(const int boundary1Type,
	const double boundary1[])
{
	switch (boundary1Type)
	{
		case RECTANGLE:
			if(boundary1[1] < boundary1[0] || boundary1[3] < boundary1[2])
				return 0.;
			else
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
			return 4/3*PI*boundary1[3]*boundary1[3]*boundary1[3];
		default:
			fprintf(stderr,"ERROR: Boundary type %d invalid.\n", boundary1Type);
			return 0;
	}
}

// Find a random coordinate within the specified range
double uniformPoint(double rangeMin,
	double rangeMax)
{
	return (rangeMin + (rangeMax-rangeMin)*mt_drand());
}

// Find a random coordinate within the specified boundary
void uniformPointVolume(double point[3],
	const int boundaryType,
	const double boundary1[])
{
	bool bNeedPoint;
	
	switch (boundaryType)
	{
		case RECTANGLE:
			return;
		case RECTANGULAR_BOX:
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
				point[0] = mt_drand();
				point[1] = mt_drand();
				point[2] = mt_drand();
				if (squareDBL(point[0]) + squareDBL(point[1])
					+ squareDBL(point[2]) < 1.)
				{
					// Found valid point. Scale as needed and randomize sign
					if(mt_drand() > 0.5)
						point[0] = -point[0];
					point[0] = boundary1[0] + point[0]*boundary1[3];
					if(mt_drand() > 0.5)
						point[1] = -point[1];
					point[1] = boundary1[1] + point[1]*boundary1[3];
					if(mt_drand() > 0.5)
						point[2] = -point[2];
					point[2] = boundary1[2] + point[2]*boundary1[3];
					bNeedPoint = false;
				}
			}
			return;
		default:
			fprintf(stderr,"ERROR: Boundary type %d invalid.\n", boundaryType);
			return;
	}
}

// Find distance between 2 3D points
double pointDistance3D(const double point1[3],
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