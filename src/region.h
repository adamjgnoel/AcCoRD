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
 * Last revised for AcCoRD v1.4 (2018-08-06)
 *
 * Revision history:
 *
 * Revision v1.4 (2018-08-06)
 * - added a priori monte carlo (APMC) absorption algorithm as a new surface
 * reaction type. Includes settings for how to define the a priori absorption
 * probability calculation and whether/how to apply a threshold to turn it off
 * - corrected check on molecule transitions between parent and child regions
 * when a molecule is initialized at the parent/child boundary
 *
 * Revision v1.1 (2016-12-24)
 * - added members defining flow parameters for every molecule
 * - added member to denote which molecules can diffuse
 *
 * Revision v1.0 (2016-10-31)
 * - added local diffusion coefficients that can apply to particular region
 * - added specifying diffusion coefficient that applies to specific surface
 * interaction reactions.
 * - moved mesoscopic structure fields from subvolume struct to meso subvolume struct
 *
 * Revision v0.7.0.1 (public beta, 2016-08-30)
 * - corrected calculating region volume when a normal region has a surface child
 * inside
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - added members to track the direction of microscopic subvolumes from mesoscopic
 * subvolumes. Replaced array for storing coordinates of virtual microscopic
 * subvolume with array of boundary coordinates of boundary mesoscopic subvolumes. These
 * changes needed to accommodate improved hybrid transition algorithms
 * - added members to implement bimolecular chemical reactions in microscopic regime.
 * Use a preliminary algorithm where the user must supply the binding and unbinding radii
 * - changed findNearestSub function to return index of subvolume in region's neighID array
 * instead of the global subvolume list. This makes function suitable for more calls.
 * - added case for point boundary when determining intersection with a region
 *
 * Revision v0.5.1 (2016-05-06)
 * - updated function bPointInRegionNotChild to take an extra input to indicate
 * whether to ignore children regions that are surfaces
 * - added bReleaseProduct for surface reactions to indicate whether products are
 * released from the surface
 * - added members to store unique properties needed for absorbing and
 * desorbing reactions
 *
 * Revision v0.5 (2016-04-15)
 * - re-structured region array initialization to nest more code in functions
 * - added more checks on region parameters (including label uniqueness) to verify placement
 * - pushed error exit to end of region initialization so that all errors will be displayed
 * before exiting.
 * - added 2D regions
 * - added type member to spec and plane, dimension, and effectiveDim members to main
 * struct in order to accommodate surface and other 2D regions
 * - added region label when giving errors about region initialization
 * - adjusted clearance between spherical and rectangular regions such that the clearance
 * between them (when one is nested inside the other) is scaled by the subvolume adjacency
 * error
 * - corrected "locking" to spherical region so that the coordinate that is "locked" is
 * the one that is the furthest from the center of the sphere
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
	
	// Is the region microscopic?
	bool bMicro;
	
	// Does this region use different diffusion coefficients?
	bool bLocalDiffusion;
	
	// Local diffusion coefficients (if applicable)
	// Length is NUM_MOL_TYPES
	double * diffusion;
	
	// Can each type of molecule flow?
	// Length is NUM_MOL_TYPES
	bool * bFlow;
	
	// Is the definition for this molecule's flow different from the global setting?
	// Length is NUM_MOL_TYPES
	bool * bFlowLocal;
	
	// What type of flow does each type of molecule experience?
	// Types are defined in global_param.h
	// Length is NUM_MOL_TYPES
	unsigned short * flowType;
	
	// Flow vector/parameters for each type of molecule
	// Size is NUM_MOL_TYPES x VECTOR_LENGTH (size of VECTOR_LENGTH depends on the type of flow)
	double ** flowVector;
	
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
	int type;
	
	// Surface type. Default is NO_SURFACE, which corresponds to REGION_NORMAL
	int surfaceType;
	
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
	
	// Does each molecule type diffuse?
	// Length is NUM_MOL_TYPES
	bool * bDiffuse;
	
	// Is at least one type of molecule carried by flow?
	// If true and region is mesoscopic, then additional information
	// must be tracked for each subvolume.
	bool bFlow;
	
	// Flow parameters for each type of molecule.
	// Intended to apply for microscopic time steps of a "typical" length
	// Size is NUM_MOL_TYPES x Y, where Y depends on the type of flow
	double ** flowConstant;
	
	// What plane is the region in?
	// Applies particularly to 2D regions.
	// Values are defined in global_param.h
	short plane;
	
	// Number of faces. 0 if region is not a surface
	short numFace;
	
	// What shape do the region's subvolumes have?
	// May not match region shape if region is a surface
	short subShape;
	
	// Dimension and effective dimension of shape
	// Dimension depends only on the shape
	// The effective dimension depends on the shape and the type
	short dimension;
	short effectiveDim;
	
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
	
	// Does a microscopic region have at least 1 mesoscopic neighbor?
	bool bHasMesoNeigh;
	
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
	// ID is from mesoscopic subvolume list and NOT GLOBAL list
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
	double (** boundSubCenterCoor)[3];
	
	// 2D array of subvolume boundary coordinates for subvolumes in current region that
	// border each neighbour. Only needed if region is mesoscopic and neighbor
	// is microscopic. Size is NUM_REGIONS x numSubRegionNeigh x 6
	double (** boundSubCoor)[6];
	
	// 3D array of directions of virtual subvolumes that neighbor each
	// subvolume in current region
	// Size is NUM_REGIONS x numSubRegionNeigh x boundSubNumFace
	unsigned short *** boundVirtualNeighDir;
	
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
	
	// ID of reaction in global chemical reaction list
	// Needed for bimolecular reactions where reactants could be in different regions
	// Size is numChemRxn
	unsigned short * globalRxnID;
	
	// Can global reaction occur in this region?
	// Needed for bimolecular reactions where reactants could be in different regions
	// Size is MAX_RXNS
	bool * bGlobalRxnID;
	
	// Is reaction reversible?
	// Size is numChemRxn
	bool * bReversible;
	
	// ID of reverse reaction if applicable
	// Size is numChemRxn
	unsigned short * reverseRxnID;
	
	// Magnitude of stoichiometry changes associated with each possible chemical reaction
	// Size is numChemRxn x NUM_MOL_TYPES
	uint64_t ** numMolChange;
	
	// Sign of stoichiometry changes associated with each possible chemical reaction
	// Size is numChemRxn x NUM_MOL_TYPES
	bool ** bMolAdd;
	
	// Number of products associated with each chemical reaction
	// Size is numChemRxn
	uint32_t * numRxnProducts;
	
	// Molecule ID of each product associated with each chemical reaction
	// Size is numChemRxn x numRxnProducts
	unsigned short ** productID;
	
	// Bool indicating whether each product is released from region surface
	// Size is numChemRxn x numRxnProducts
	bool ** bReleaseProduct;
	
	// Type of product release from surface
	// Default is PROD_PLACEMENT_LEAVE, which will leave the molecule next to
	// the surface.
	// Used only when bReleaseProduct is true
	// Size is numChemRxn
	short * releaseType;
	
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
	unsigned short * numFirstRxnWithReactant; // Length NUM_MOL_TYPES
	
	// What are the IDs of the first order reactions that have the current molecule
	// as the reactant?
	unsigned short ** firstRxnWithReactantID; // Length NUM_MOL_TYPES x numChemRxn
	
	// How many second order reactions have the current molecule as a reactant?
	unsigned short * numSecondRxnWithReactant; // Length NUM_MOL_TYPES
	
	// For bimolecular reactions, what are the indices of the two reactants?
	// Size is numChemRxn x 2
	unsigned short (* biReactants)[2];
	
	// Used in micro regions only	
	double * uniSumRate; // Length NUM_MOL_TYPES
	double ** uniCumProb; // Length NUM_MOL_TYPES x numChemRxn
	double ** uniRelativeRate; // Length NUM_MOL_TYPES x numChemRxn
	double * minRxnTimeRV; // uniform RV must be generated in the range [minRxnTimeRV,1]
						 // in order to correspond to a reaction time that is less than
						 // the region's micro time step.
						 // Length NUM_MOL_TYPES
	
	double * rBind;		// Bimolecular reaction binding radius. Length numChemRxn
	double * rBindSq;	// Bimolecular reaction binding radius squared. Length numChemRxn
	double * rUnbind;	// Bimolecular reaction Unbinding radius. Length numChemRxn
	double rBindMax;	// Largest binding radius.
						// TODO: Use for determining micro subvolume neighbors.
						 
	// Type of surface reaction probability calculation
	// Default is RXN_PROB_NORMAL, which is inaccurate for surface transition reactions
	// Length numChemRxn
	short * rxnProbType;
	
	// Diffusion coefficient used for surface reaction probability calculation
	// Length numChemRxn
	double * rxnDiffCoef;
	
	// Does an absorbing or membrane reaction exist for a given molecule?
	// Membrane reaction is for passing from the "inner" direction
	// Length NUM_MOL_TYPES
	bool * bSurfRxnIn;
	
	// ID of absorbing-type reaction if bSurfRxnIn == true
	// Length NUM_MOL_TYPES
	unsigned short * rxnInID;
	
	// Probability of absorbing or membrane reaction for a given molecule
	// Membrane reaction is for passing from the "inner" direction
	// Length NUM_MOL_TYPES
	double * surfRxnInProb;
	
	// Does a desorbing or membrane reaction exist for a given molecule?
	// Membrane reaction is for passing from the "outer" direction
	// Length NUM_MOL_TYPES
	bool * bSurfRxnOut;
	
	// ID of desorbing-type reaction if bSurfRxnIn == true
	// Length NUM_MOL_TYPES
	unsigned short * rxnOutID;
	
	// Probability of desorbing or membrane reaction for a given molecule
	// Membrane reaction is for passing from the "outer" direction
	// Length NUM_MOL_TYPES
	double * surfRxnOutProb;
	
	// Do we use the independently-calculated desorption probability in surfRxnOutProb?
	// surfRxnOutProb is only needed for reversible steady-state desorption
	// Membrane reaction is for passing from the "outer" direction
	// Length NUM_MOL_TYPES
	bool * bUseRxnOutProb;
	
	// A Priori surface reaction parameters (micro only)
	
	// Number of A Priori surface reactions that apply to given molecule type
	// Length NUM_MOL_TYPES
	unsigned short * numApmcRxn;
	
	// Regions of A priori surface reactions that apply to given molecule type
	// Size NUM_MOL_TYPES x numApmcRxn
	unsigned short ** apmcRxnRegion;
	
	// IDs of A priori surface reactions that apply to given molecule type
	// Indexing is from reacting region reaction list
	// Size NUM_MOL_TYPES x numApmcRxn
	unsigned short ** apmcRxnID;
	
	// IDs of A priori surface reactions that apply to given molecule type
	// Indexing is from global reaction list and not region where corresponding
	// reaction is defined
	// Size NUM_MOL_TYPES x numApmcRxn
	unsigned short ** apmcGlobalRxnID;
	
	// Cumulative probabilities for current A priori reactions
	// Size NUM_MOL_TYPES x numApmcRxn
	double ** uniCumProbApmc;
	
	// Intermediate parameter alpha for A Priori surface reactions
	// Needed if the absorption has a finite reaction rate
	// Length NUM_MOL_TYPES
	double * apmcAlpha;
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
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const double SUBVOL_BASE_SIZE,
	const double boundAdjError,
	struct subvolume3D * subvolArray,
	uint32_t subCoorInd[][3]);

// Free memory of region parameters
void delete_boundary_region_(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	struct region regionArray[]);

// Find volume of specified region (exluding volume of children)
double findRegionVolume(const struct region regionArray[],
	const short curRegion,
	bool bOuter);
	
// Count the cumulative number of subvolumes defined by subvol_spec
uint32_t countAllSubvolumes(const struct region regionArray[],
	const short NUM_REGIONS);

// Determine which regions of subvolumes are adjacent
void findRegionTouch(const short NUM_REGIONS,
	const struct spec_region3D subvol_spec[],
	struct region regionArray[],
	const double SUBVOL_BASE_SIZE);

// Find index of desired subvolume in list defining region's boundary with another region
uint32_t findSubInBoundaryList(const short curRegion,
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
uint32_t findNearestSub(const short curRegion,
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
	const double point[3],
	bool bIgnoreSurfaceChildren);

// Is point in a region or any of its nested children?
// If true, actualRegion will be the ID of the region with the point
bool bPointInRegionOrChild(const short curRegion,
	const struct region regionArray[],
	const double point[3],
	short * actualRegion,
	bool bSurfaceOnly);

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
	const short neighRegion,
	const double curSubBound[6],
	double boundAdjError,
	unsigned short * numFace,
	unsigned short dirArray[6]);

// Initialize the region nesting (i.e., determine each region's parent and
// children regions, if applicable)
void initializeRegionNesting(const short NUM_REGIONS,
	struct region regionArray[],
	bool * bFail);

// Determine number of subvolumes in each region
void findNumRegionSubvolumes(const short NUM_REGIONS,
	struct region regionArray[]);

// Allocate and initialize the flow parameters in each region
void initializeRegionFlow(const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	struct region regionArray[],	
	const struct spec_region3D subvol_spec[]);

// Allocate memory for each region's neighbors
void allocateRegionNeighbors(const short NUM_REGIONS,
	struct region regionArray[]);

// Final check of region overlap and correct adjacency
void validateRegions(const short NUM_REGIONS,
	const struct region regionArray[],
	bool * bFail);

// Check whether the children of a surface region cover all of the
// region's outer boundary
void checkSurfaceRegionChildren(const short curRegion,
	const short NUM_REGIONS,
	const struct region regionArray[],
	bool * bFail);

// Does a one-sided surface exist between 2 rectangular boundaries?
// Such a surface would prevent the boundaries from being neighbors
bool bSurfaceBetweenBoundaries(const struct region regionArray[],
	const short NUM_REGIONS,
	const short region1,
	const short region2,
	const double boundary1[6],
	const double boundary2[6],
	unsigned short * surfaceRegion);

// Does a boundary intersect a specified region?
// Determined by measuring volume of intersection.
// Intersection should be rectangular. Region can be spherical if it
// fully contains or is contained by the other boundary
// Round children must be entirely inside or outside intersection
bool bIntersectRegion(const short curRegion,
	const struct region regionArray[],
	const int boundary2Type,
	const double boundary2[]);

#endif // REGION_H
