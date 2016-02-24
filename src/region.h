/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * region.h - 	operations for (microscopic or mesoscopic) regions in
 * 				simulation environment
 *
 * Last revised for AcCoRD LATEST_RELEASE
 *
 * Revision history:
 *
 * Revision LATEST_RELEASE
 * - added 2D regions
 * - added type member to spec and plane member to main struct in order to accommodate
 * surface and other 2D regions
 * - added region label when giving errors about region initialization
 * - adjusted clearance between spherical and rectangular regions such that the clearance
 * between them (when one is nested inside the other) is adjusted by the subvolume adjacency
 * error
 *
 * Revision v0.4.1
 * - improved use and format of error messages
 *
 * Revision v0.4
 * - modified use of unsigned long and unsigned long long to uint32_t and uint64_t
 * - added accommodation of spherical (microscopic) subvolumes and actors
 * - created new structure members to account for meso subvolumes that border
 * the same micro region along more than 1 face
 * - renamed region structure and array of structures to accommodate naming
 * of new boundary structure
 *
 * Revision v0.3.1
 * - header added
*/
#ifndef REGION_H
#define REGION_H

//#include "subvolume.h"
#include "base.h" // For region adjacency
#include "global_param.h" // For region adjacency

/*
* Data Type Declarations
*/
struct subvolume3D; // Forward declaration (declared in subvolume.h)
					// This declaration avoids header coupling (since subvolume.h
					// includes region.h), but means we can only refer to subvolume3D
					// as a pointer in this file
struct chem_rxn_struct; // (declared in chem_rxn.h)

/* The spec_region3D structure includes the parameters for a region that are
* defined by the user input configuration
*/
struct spec_region3D { // Used to define a region of subvolumes
	// The "anchor" point is the "lower" coordinates that define the
	// location of the region, or the centre of a round region
	double xAnch, yAnch, zAnch;
	
	// label is an optional string to name the region.
	// parent is an optional string that indicates the label of the region
	// that this region is nested in.
	char * label;
	char * parent;
	
	// Is the region actually a surface?
	// Surfaces are hollow regions with a number of unique properties
	// The primary characteristic is that they control transitions
	// between regions that are on each side of the surface.
	// Surfaces can have parents if they are embedded within a region(s)
	// Surfaces can have children if their shape is 3D
	bool bSurface;
	
	// Is the region microscopic?
	bool bMicro;
	
	// Region shape. There are specific restrictions on shape, depending on
	// whether it is micro, whether its parent is micro, and whether it is
	// a surface. Generally, RECTANGULAR_BOX is fine, but round shapes are
	// restricted to microscopic regions that do not have mesoscopic parents
	int shape;
	
	// Region type. Default is REGION_NORMAL, which means that the region
	// occupies the entire volume defined and does not impede diffusion across
	// its boundary.
	// REGION_SURFACE means that the region is "hollow" and impedes (though
	// may not necessarily prevent) diffusion across its boundary
	// REGION_MEMBRANE is similar to REGION_SURFACE but is specifically a
	// selective diffusion barrier (may not need to keep)
	int type;
	
	// The simulation configuration has a defined base subvolume size,
	// SUBVOL_BASE_SIZE. Square/cube subvolumes have a length that is a multiple
	// of SUBVOL_BASE_SIZE.
	// NOTE: Does NOT apply to round regions
	unsigned int sizeRect;
	
	// Radius of round regions. Independent of SUBVOL_BASE_SIZE, since we do not need
	// to worry about perfect alignment of round regions, and round regions always
	// have one subvolume
	double radius;
	
	// The number of subvolumes along each dimension of the region. If the
	// region has no exclusions, then the product of these parameters is
	// the number of subvolumes in the region.
	// NOTE: Does NOT apply to round regions, which always have one subvolume
	unsigned int numX, numY, numZ;
	
	// Region time step
	double dt;
	
	// FUTURE MEMBERS (POTENTIAL)
	// Indicator for presence of system boundary
	// Details of region-specific reactions
	// Exclusion zones to more easily surround smaller regions
};

/* The region structure contains all parameters specific to a single
* region, including the user-defined parameters defined in spec_region3D. The
* structure members that describe the region's location relative to other regions
* are determined at runtime.
*/
struct region { // Region boundary parameters
	// The user-defined region parameters
	struct spec_region3D spec;
	
	// What plane is the region in?
	// Applies particularly to 2D regions.
	// Values are defined in global_param.h
	short plane;
	
	// Is this region nested inside another region
	bool bParent;
	
	// ID of region that this region is directly nested in
	short parentID;
	
	// # of regions nested inside this region
	short numChildren;
	
	// IDs of regions nested inside this region
	// Length is numChildren
	short * childrenID;
	
	// Relative outer boundary of children, written as indexing of the
	// regions subvolumes, i.e., lists range of subvolume coordinates that
	// excluded due to each child.
	// Size is numChildrenx6
	uint32_t (* childrenCoor)[6];
	
	// Volume (or area in 2D)
	double volume;
	
	// Resolution size
	// Used to determine how close subvolume faces must be for them to be adjacent
	double subResolution;
	
	// Index of the first subvolume in this region in the overall subvolume list
	uint32_t firstID;
	
	// Number of subvolumes in this region
	uint32_t numSub;
	
	// "True" subvolume size (i.e., spec.sizeRect scaled by SUBVOL_BASE_SIZE if rectangular)
	double actualSubSize;
	
	// Coordinates of "Outer" Region boundary
	double boundary[6];
	
	// Is there at least one boundary region along each outer boundary face?
	// TODO: Following Members may only be needed for boundary reactions
	bool boundaryRegion[6];
	bool boundaryRegionMicro[6]; // Are outer boundary regions microscopic?
	
	// Is each other region a neighbor of this region? The element in this array
	// for the same region is assigned false
	// Length is simulation-defined NUM_REGIONS
	bool * isRegionNeigh; // Are other regions neighbours of this region?
	
	// Number of other regions that are a neighbor of this region
	short numRegionNeigh;
	
	// IDs of regions that are neighbors. Length is value of numRegionNeigh
	// TODO: Might not need
	short * regionNeighID;
	
	// If region is a neighbour, in what direction is it? If region is not
	// a neighbor, corresponding element is undefined. Length is NUM_REGIONS
	// Values indicate direction as follows:
	// 0 - in -ve x-direction
	// 1 - in +ve x-direction
	// 2 - in -ve y-direction
	// 3 - in +ve y-direction
	// 4 - in -ve z-direction
	// 5 - in +ve z-direction
	unsigned short * regionNeighDir;
	
	// Number of faces of each neighbor that are boundaries to this region
	// Only needed if this region is microscopic
	// Length is NUM_REGIONS
	// TODO: Might not need since only have multiple faces in case of parent/child
	// and this is easy to check
	short * numRegionNeighFace;
	
	// 2D array of boundary coordinates for faces in boundary regions that overlap
	// this region.
	// Only needed if this region is microscopic
	// Length is NUM_REGIONS x numRegionNeighFace x 6
	// TODO: Might not need as a 2D array since we only have multiple faces in the case
	// of parent/child and those are easy to check
	double (** boundRegionFaceCoor)[6];
	
	// 2D array of boundary directions for faces in boundary regions that overlap
	// this region.
	// This is distinct from regionNeighDir because it also handles directions of
	// child and parent faces
	// Only needed if this region is microscopic
	// Length is NUM_REGIONS x numRegionNeighFace
	// TODO: Might not need as a 2D array since we only have multiple faces in the case
	// of parent/child and those are easy to check
	short ** regionNeighFaceDir;
	
	// Number of subvolumes in current region that are along boundary shared
	// with each neighboring region. Length is NUM_REGIONS
	// Only needed if region is mesoscopic and neighbor is microscopic
	uint32_t * numSubRegionNeigh;
	
	// 2D array of IDs of subvolumes in current region that border each neighbor.
	// Only needed if region is mesoscopic and neighbor is microscopic.
	// ID is from GLOBAL subvolume list and NOT mesoscopic list
	// Size is NUM_REGIONS x numSubRegionNeigh
	// For details of memory allocation to 2D array with a different number of
	// columns for each row: http://c-faq.com/aryptr/dynmuldimary.html
	uint32_t ** neighID;
	
	// 2D array indicating how many faces of a given mesoscopic subvolume
	// are adjacent to a given microscopic neighbor region.
	// Only needed if region is mesoscopic and neighbor is microscopic.
	// Size is NUM_REGIONS x numSubRegionNeigh
	unsigned short ** boundSubNumFace;
	
	// 2D array of subvolume coordinates for subvolumes in current region that
	// border each neighbour. Only needed if region is mesoscopic and neighbor
	// is microscopic. Size is NUM_REGIONS x numSubRegionNeigh x 3
	double (** boundSubCoor)[3];
	
	// 3D array of virtual subvolume anchor coordinates associated with each
	// subvolume in current region that borders each neighbor.	
	// Only needed if region is mesoscopic and neighbor
	// is microscopic.
	// Size is NUM_REGIONS x numSubRegionNeigh x boundSubNumFace x 3
	double (*** boundVirtualNeighCoor)[3];
	
	// 2D array of booleans to indicate whether a boundary subvolume needs its
	// propensities updated due to molecules entering from the microscopic regime.
	// Only needed if region is mesoscopic and neighbor is microscopic
	// Size is NUM_REGIONS x numSubRegionNeigh
	bool ** bNeedUpdate;
	
	// 3D array of number of molecules entering boundary subvolumes from
	// microscopic neighbor.
	// Only needed if region is mesoscopic and neighbor is microscopic
	// Size is NUM_REGIONS x numSubRegionNeigh x NUM_MOL_TYPES
	uint64_t *** numMolFromMicro;
	
	// Array of diffusion rates of molecules between subvolumes within region
	// Only needed if region is mesoscopic
	// Length is NUM_MOL_TYPES
	double * diffRate;
	
	//
	// Chemical Reaction Parameters
	// TODO: Some parameters are specific to meso or micro regions only. They could
	// re-located into other structures
	//
	
	// Number of chemical reactions possible in this region
	unsigned short numChemRxn;
	unsigned short numZerothRxn;
	unsigned short numFirstRxn;
	unsigned short numSecondRxn;
	
	// Magnitude of stoichiometry changes associated with each possible chemical reaction
	// Size is numChemRxn x NUM_MOL_TYPES
	uint64_t ** numMolChange;
	
	// Sign of stoichiometry changes associated with each possible chemical reaction
	// Size is numChemRxn x NUM_MOL_TYPES
	bool ** bMolAdd;
	
	// Number of products associated with each chemical reaction
	// Size is numChemRxn; maximum value is MAX_RXN_PRODUCTS (defined in global_param.h)
	uint32_t * numRxnProducts;
	
	// Molecule ID of each product associated with each chemical reaction
	// Size is numChemRxn x MAX_RXN_PRODUCTS (defined in global_param.h)
	unsigned short ** productID;
	
	// Indicator of what molecule number changes mean that a reaction propensity
	// must be updated.
	// This must be specified because the reactants cannot always be inferred from
	// stoichiometry changes alone (e.g., if a reactant is also a product)
	// Size is numChemRxn x NUM_MOL_TYPES
	bool ** bUpdateProp;
	
	// Order of each chemical reaction. This must also be specified because the
	// reactants cannot always be inferred from the stoichiometry changes.
	// Size is numChemRxn
	uint32_t * rxnOrder;
	
	// Base rate of each chemical reaction. This must be defined at the region level
	// since 0th and 2nd order reactions have rates that depend on the subvolume size
	// Size is numChemRxn
	double * rxnRate;
	
	// Indices of reactions that are 0th order.
	// Size is numChemRxn; elements numZerothRxn and greater are undefined
	unsigned short * zerothRxn;
	
	// Indices of reactions that are 1st order.
	// Size is numChemRxn; elements numFirstRxn and greater are undefined
	unsigned short * firstRxn;
	
	// Indices of reactions that are 2nd order.
	// Size is numChemRxn; elements numSecondRxn and greater are undefined
	unsigned short * secondRxn;
	
	// Reaction time and rate of zeroth order reactions
	// Used in micro regions only
	// Size is numChemRxn; elements numZerothRxn and greater are undefined
	double * tZeroth;
	double * rxnRateZerothMicro;
	
	// For unimolecular reactions, what is the index of the reactant?
	// Undefined for reactions that are not unimolecular.
	// Size is numChemRxn
	unsigned short * uniReactant;
	
	// For unimolecular reactions, what are the cumulative reaction
	// probabilities for each molecule type?
	
	// How many first order reactions have the current molecule as the reactant?
	// TODO: Improve member name
	unsigned short * numFirstCurReactant; // Length NUM_MOL_TYPES
	
	// What are the IDs of the chemical reactions that have the current molecule
	// as the reactant.
	// TODO: Improve member name
	unsigned short ** firstRxnID; // Length NUM_MOL_TYPES x numChemRxn
	
	// Used in micro regions only	
	double * uniSumRate; // Length NUM_MOL_TYPES
	double ** uniCumProb; // Length NUM_MOL_TYPES x numChemRxn
	double ** uniRelativeRate; // Length NUM_MOL_TYPES x numChemRxn
	double * minRxnTimeRV; // uniform RV must be generated in the range [minRxnTimeRV,1]
						 // in order to correspond to a reaction time that is less than
						 // the region's micro time step.
						 // Length NUM_MOL_TYPES
	
	// For bimolecular reactions, what are the indices of the two reactants?
	// Undefined for reactions that are not bimolecular.
	// Size is numChemRxn x 2
	unsigned short (* biReactants)[2];
};

/*
* Function Declarations
*/

// Initialize array of region structs from specification
void initializeRegionArray(struct region regionArray[],
	const struct spec_region3D subvol_spec[],
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const double SUBVOL_BASE_SIZE,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	const unsigned short MAX_RXNS,
	const struct chem_rxn_struct * chem_rxn);

// Initialize region knowledge of the subvolumes that are adjacent to it
void initializeRegionSubNeighbor(struct region regionArray[],
	const struct spec_region3D subvol_spec[],
	const double subHalfSize[],
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const double SUBVOL_BASE_SIZE,
	const double boundAdjError,
	struct subvolume3D * subvolArray,
	uint32_t subCoorInd[][3]);

// Free memory of region parameters
void delete_boundary_region_3D(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	struct region regionArray[]);

// Find volume of specified region (exluding volume of children)
double findRegionVolume(const struct region regionArray[],
	const short curRegion,
	bool bOuter);
	
// Count the cumulative number of subvolumes defined by subvol_spec
uint32_t count_subvol3D(const struct region regionArray[],
	const short NUM_REGIONS);

// Determine which regions of subvolumes are adjacent
void findRegionTouch3D(const short NUM_REGIONS,
	const struct spec_region3D subvol_spec[],
	struct region regionArray[],
	const double SUBVOL_BASE_SIZE);

// Find index of desired subvolume in list defining region's boundary with another region
uint32_t find_sub_in_bound_list3D(const short curRegion,
	const short destRegion,
	const struct region regionArray[],
	uint32_t curSub);

// Find closest region for point to be in
unsigned short findNearestValidRegion(const double point[],
	const unsigned short sourceRegion,
	bool * bInRegion,
	const unsigned short NUM_REGIONS,
	const struct region regionArray[]);

// Find the closest subvolume in current region that is along boundary
// of specified neighbor region
uint32_t findNearestSub3D(const short curRegion,
	const struct region regionArray[],
	const short neighRegion,
	double x,
	double y,
	double z);

// Determine coordinates of child region as subvolumes of parent
// Applies to Rectangular regions only.
void findChildCoordinates(const short parentRegion,
	const short childRegion,
	const short curChild,
	const struct region regionArray[]);

// Volume of intersection within specified region
// Value returned here excludes children
// Intersection must be rectangular
// Round children must be entirely inside or outside intersection
double intersectRegionVolume(const short curRegion,
	const struct region regionArray[],
	const int boundary2Type,
	const double boundary2[]);

// Does a line hit the specified region? If so then where and at what distance?
bool bLineHitRegion(const double p1[3],
	const double L[3],
	const double length,
	const short startRegion,
	const short endRegion,
	const struct region regionArray[],
	short * planeID,
	double * d,
	double intersectPoint[3]);

// Is point inside region and not one of its children?
bool bPointInRegionNotChild(const short curRegion,
	const struct region regionArray[],
	const double point[3]);

// Is point in a region or any of its nested children?
// If true, actualRegion will be the ID of the region with the point
bool bPointInRegionOrChild(const short curRegion,
	const struct region regionArray[],
	const double point[3],
	short * actualRegion);

// Is a specific face shared between two regions?
bool bSharedBoundary(const short startRegion,
	const short endRegion,
	const struct region regionArray[],
	const short faceID);

// Is the outer boundary of a child region flush with the subvolumes of its
// parent? Applies to rectangular regions (2D or 3D)
bool bChildRegionNotFlush(const short parentRegion,
	const short childRegion,
	const struct region regionArray[],
	const double boundAdjError);

// Lock point coordinate to chosen region face
void lockPointToRegion(double point[3],
	const short startRegion,
	const short endRegion,
	const struct region regionArray[],
	const short faceID);

// Generate a random cartesian point in the specified region
void generatePointInRegion(const short curRegion,
	const struct region regionArray[],
	double point[3]);

// Which region contains given point, excluding children?
short findRegionNotChild(const short NUM_REGIONS,
	const struct region regionArray[],
	double point[3]);

// Does a subvolume face a region? If yes, then along which faces?
// Assert that current subvolume is mesoscopic, along its own region boundary, and that
// neighbor region is microscopic
bool bSubFaceRegion(struct region regionArray[],
	const short curRegion,
	const short neighRegion,
	const double curSubBound[6],
	const double subHalfSize,
	const uint32_t subCoorInd[3],
	double boundAdjError,
	unsigned short * numFace,
	unsigned short dirArray[6]);

// Final check of region overlap and correct adjacency
void validateRegions(const short NUM_REGIONS,
	const struct region regionArray[]);

#endif // REGION_H
