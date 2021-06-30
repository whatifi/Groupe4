#include "RandomNumberGenerator.h"
#include <cmath>

double RandomNumberGenerator::uniformRandomNumber()
{
	return ((double)(rand()) + 1.) / ((double)(RAND_MAX)+1.);
}

// Box-Muller method
double RandomNumberGenerator::normalRandomNumber()
{
	double u1 = uniformRandomNumber();
	double u2 = uniformRandomNumber();
	return cos(8.*atan(1.)*u2)*sqrt(-2.*log(u1));
}