#ifndef RANDOMNUMBERGENERATOR_H
#define RANDOMNUMBERGENERATOR_H

class RandomNumberGenerator
{
public:
	// Those functions generate uniforms U[0,1] and normals N(0,1)
	static double uniformRandomNumber();
	static double normalRandomNumber();
};

#endif