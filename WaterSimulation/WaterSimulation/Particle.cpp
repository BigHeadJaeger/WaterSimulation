#include "Particle.h"

void Particle::SetInitialN0(vector<vec3> r, int currentIndex)
{
	MPSToolFun* tool = MPSToolFun::GetMPSTool();
	n0 = tool->DensityN(r, currentIndex);
}

void Particle::UpdateAdjoin(vector<Particle>& particles, float range)
{
	adjoinParticleIndex.clear();
	for (int i = 0; i < particles.size(); i++)
	{
		if (i!=index)
		{
			vec3 pos1 = particles[i].position;
			if (length(pos1 - position) <= range)
				adjoinParticleIndex.push_back(i);
		}
	}
}

void MPSWaterParticle::SurfaceAdjudge(float a, float tho, float g, float l0)
{
	if (pressure < (a * tho * g * l0))
		isSurface = true;
	else
		isSurface = false;
}
