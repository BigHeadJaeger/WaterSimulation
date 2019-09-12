#include "Particle.h"

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
