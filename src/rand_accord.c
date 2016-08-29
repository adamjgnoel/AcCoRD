/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * rand_accord.c - interface to random number generation functions
 *
 * Last revised for AcCoRD v0.7.0.1 (public beta, 2016-08-30)
 *
 * Revision history:
 *
 * Revision v0.7.0.1 (public beta, 2016-08-30)
 * - modified range of uniform RV generation to be [0, 1) instead of [0,1].
 * 	Allowing 1 was causing problems when RV was being used for indexing
 *
 * Revision v0.6 (public beta, 2016-05-30)
 * - file created to centralize calls to random number generators
 * - replacing RNG via Mersenne Twister with the PCG algorithm
 * (http://www.pcg-random.org/)
*/

#include "rand_accord.h" // for "Public" declarations

// Initialize the random number generator
void rngInitialize(const uint32_t SEED)
{	
	// 2nd argument is sequence selector. Sequence value doesn't really matter.
	// Sequence value chosen is same as that used in PCG demo code
	pcg32_srandom(SEED, 54u);
}

// Return a uniform random number between 0 and 1 INCLUSIVE
double generateUniform()
{	
	// Using conversion which gives us range of [0, 1)
	return (double) pcg32_random()/(UINT32_MAX + 1.);
}

// Return a triangular random number
// Value will be in range [0,2]
double generateTriangular()
{
	double uniRV = (double) pcg32_random()/UINT32_MAX;
	
	if(uniRV < 0.5)
		return sqrt(2*uniRV);
	else
		return (2-sqrt(2*(1-uniRV)));
}

// Return a normal random number via Box-Muller transform
// Transform will generate 2 normal RVs; 1 is kept to return in the following call
double generateNormal(const double mean,
	const double std)
{
	double magnitude;
	static double offset;
	double xVal;
	static double yVal; // static so that we can return it in a future call
	static bool bVal2 = false;
	
	if(bVal2)
	{ // We previously generated a second value
		bVal2 = false;
		return mean + std * yVal * offset;
	}
	bVal2 = true;
	
	do
	{
		// Generate 2 uniform RVs within (-1,1)
		xVal = 2.*generateUniform() - 1.0;
		yVal = 2.*generateUniform() - 1.0;
		magnitude = xVal*xVal + yVal*yVal;
	} while (magnitude > 1.0 || magnitude == 0.);
	
	// Convert uniform RV to normal
	offset = sqrt(-2.0*log(magnitude)/magnitude);
	
	// Return xVal (yVal will be returned in next call)
	return mean + std*xVal*offset;
}
	
// Return a Poisson random number
uint64_t generatePoisson(double mean)
{
	uint64_t poissonVal = 0;
	double tSum;
	
	while(1)
	{
		tSum += generateExponential(mean);
		if(tSum >= 1.)
			break;
		poissonVal++;
	}
	
	return poissonVal;
}

// Return an Exponential random number
double generateExponential(double mean)
{
	double randVal;
		
	do
	{
		randVal = generateUniform();
	} while(randVal == 0.0); // Avoid 0 case
	return -mean * log(randVal);
}