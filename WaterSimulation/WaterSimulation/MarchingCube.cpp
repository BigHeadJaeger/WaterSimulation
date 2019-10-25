#include "MarchingCube.h"

void MarchingCube::GetMeshData(vector<vec3>& sourcePoints, vector<float>& verticesInfo, bool& provideNormal, bool& provideTex)
{
	sourceData = &sourcePoints;
}


MarchingCube* MarchingCube::instance = NULL;