#include "Particle.h"

void MPSWaterParticle::SurfaceAdjudge(float a, float g, float l0,float thoNow)
{
	if (pressure < (a * thoNow * g * l0))
		isSurface = true;
	else
		isSurface = false;
}

void MPSWaterParticle::OldSurfaceAdjudge(float b, float nNow)
{
	if (nNow < (b*n0))
		isSurface = true;
	else
		isSurface = false;
}
