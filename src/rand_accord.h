/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * rand_accord.h - interface to random number generation functions
 *
 * Last revised for AcCoRD v0.6 (public beta, 2016-05-30)
 *
 * Revision history:
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - file created to centralize calls to random number generators
 * - replacing RNG via Mersenne Twister with the PCG algorithm
 * (http://www.pcg-random.org/)
*/

#ifndef RAND_ACCORD_H
#define RAND_ACCORD_H

#include <math.h>
#include <stdbool.h> // for C++ bool naming, requires C99
#include <stdlib.h> // llabs definition (otherwise compiler warning for inttypes)
#include <inttypes.h> // for extended integer type macros
#include "pcg_basic.h" // For PRNG based on PCG

//
// Data Type Declarations
//

//
// Function Declarations
//

// Initialize the random number generator
void rngInitialize(const uint32_t SEED);

// Return a uniform random number between 0 and 1 INCLUSIVE
double generateUniform();

// Return a triangular random number
double generateTriangular();

// Return a normal random number
double generateNormal(const double mean,
	const double variance);
	
// Return a Poisson random number
uint64_t generatePoisson(double mean);
	
// Return an Exponential random number
double generateExponential(double mean);

#endif // RAND_ACCORD_H
