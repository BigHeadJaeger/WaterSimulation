#include "Particle.h"

void MPSWaterParticle::SurfaceAdjudge(float a, float g, float l0)
{
	if (pressure < (a * tho * g * l0))
		isSurface = true;
	else
		isSurface = false;
}
