#include "ParticleGroup.h"

void MPSWaterParticleGroup::InitParticles()
{

}
void MPSWaterParticleGroup::SetInitialN0()
{

	MPSToolFun* tool = MPSToolFun::GetMPSTool();
	for (int i = 0; i < particles.size(); i++)
	{
		vector<vec3> neighborPos;
		for (int j = 0; j < particles[i].adjoinParticleIndex.size(); j++)
		{
			neighborPos.push_back(particles[particles[i].adjoinParticleIndex[j]].position);
		}
		particles[i].n0 = tool->DensityN(neighborPos, i);
	}
	
}

void MPSWaterParticleGroup::UpdateAdjoin(float range)
{
	for (int i = 0; i < particles.size(); i++)
	{
		particles[i].adjoinParticleIndex.clear();
		for (int j = 0; j < particles.size(); j++)
		{
			if (j != i)
			{
				vec3 pos1 = particles[j].position;
				if (distance(pos1, particles[i].position) <= range)
					particles[i].adjoinParticleIndex.push_back(j);
			}
		}
	}

}