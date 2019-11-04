#include "MarchingCube.h"

GLfloat MarchingCube::Sample(GLfloat fX, GLfloat fY, GLfloat fZ, float r)
{
	double result = 0.0;
	for (int i = 0; i < sourceData->size(); i++)
	{
		float fDx, fDy, fDz;
		fDx = fX - (*sourceData)[i].x;
		fDy = fY - (*sourceData)[i].y;
		fDz = fZ - (*sourceData)[i].z;
		result += r*r / (fDx * fDx + fDy * fDy + fDz * fDz);			//1为球的半径的平方，当前点在球外时，这个式子返回值大于1
	}

	return result;
}

void MarchingCube::MarchCube(GLfloat fX, GLfloat fY, GLfloat fZ, GLfloat fScale,vector<float>& verticesInfo)
{
	GLint iCorner, iVertex, iVertexTest, iEdge, iTriangle, iFlagIndex, iEdgeFlags;
	GLfloat fOffset;
	vec3 sColor;
	GLfloat afCubeValue[8];
	vec3 asEdgeVertex[12];
	vec3 asEdgeNorm[12];

	//Make a local copy of the values at the cube's corners
	for (iVertex = 0; iVertex < 8; iVertex++)
	{
		afCubeValue[iVertex] = Sample(fX + a2fVertexOffset[iVertex][0] * fScale,			//遍历传进来的点(fx,fy,fy)是立方体的角坐标，
			fY + a2fVertexOffset[iVertex][1] * fScale,										//此处通过换算得到当前小立方体的8个顶点
			fZ + a2fVertexOffset[iVertex][2] * fScale, radius);
	}

	//Find which vertices are inside of the surface and which are outside
	iFlagIndex = 0;
	for (iVertexTest = 0; iVertexTest < 8; iVertexTest++)
	{
		if (afCubeValue[iVertexTest] <= fTargetValue)
			iFlagIndex |= 1 << iVertexTest;
	}

	//Find which edges are intersected by the surface
	iEdgeFlags = aiCubeEdgeFlags[iFlagIndex];

	//If the cube is entirely inside or outside of the surface, then there will be no intersections
	if (iEdgeFlags == 0)
	{
		return;
	}

	//Find the point of intersection of the surface with each edge
	//Then find the normal to the surface at those points
	for (iEdge = 0; iEdge < 12; iEdge++)
	{
		//if there is an intersection on this edge
		if (iEdgeFlags & (1 << iEdge))
		{
			fOffset = GetOffset(afCubeValue[a2iEdgeConnection[iEdge][0]],
				afCubeValue[a2iEdgeConnection[iEdge][1]], fTargetValue);

			asEdgeVertex[iEdge].x = fX + (a2fVertexOffset[a2iEdgeConnection[iEdge][0]][0] + fOffset * a2fEdgeDirection[iEdge][0]) * fScale;
			asEdgeVertex[iEdge].y = fY + (a2fVertexOffset[a2iEdgeConnection[iEdge][0]][1] + fOffset * a2fEdgeDirection[iEdge][1]) * fScale;
			asEdgeVertex[iEdge].z = fZ + (a2fVertexOffset[a2iEdgeConnection[iEdge][0]][2] + fOffset * a2fEdgeDirection[iEdge][2]) * fScale;

			GetNormal(asEdgeNorm[iEdge], asEdgeVertex[iEdge].x, asEdgeVertex[iEdge].y, asEdgeVertex[iEdge].z);
		}
	}

	//Draw the triangles that were found.  There can be up to five per cube
	for (iTriangle = 0; iTriangle < 5; iTriangle++)
	{
		if (a2iTriangleConnectionTable[iFlagIndex][3 * iTriangle] < 0)
			break;

		for (iCorner = 0; iCorner < 3; iCorner++)
		{
			iVertex = a2iTriangleConnectionTable[iFlagIndex][3 * iTriangle + iCorner];

			//vGetColor(sColor, asEdgeVertex[iVertex], asEdgeNorm[iVertex]);
			//glColor3f(sColor.fX, sColor.fY, sColor.fZ);
			verticesInfo.push_back(asEdgeVertex[iVertex].x);
			verticesInfo.push_back(asEdgeVertex[iVertex].y);
			verticesInfo.push_back(asEdgeVertex[iVertex].z);
			verticesInfo.push_back(asEdgeNorm[iVertex].x);
			verticesInfo.push_back(asEdgeNorm[iVertex].y);
			verticesInfo.push_back(asEdgeNorm[iVertex].z);
			verticesInfo.push_back(0);
			verticesInfo.push_back(0);
			//glNormal3f(asEdgeNorm[iVertex].fX, asEdgeNorm[iVertex].fY, asEdgeNorm[iVertex].fZ);
			//glVertex3f(asEdgeVertex[iVertex].fX, asEdgeVertex[iVertex].fY, asEdgeVertex[iVertex].fZ);
		}
	}
}

GLfloat MarchingCube::GetOffset(GLfloat fValue1, GLfloat fValue2, GLfloat fValueDesired)
{
	double fDelta = fValue2 - fValue1;

	if (fDelta == 0.0)
	{
		return 0.5;
	}
	return (fValueDesired - fValue1) / fDelta;
}

void MarchingCube::GetNormal(vec3& rfNormal, GLfloat fX, GLfloat fY, GLfloat fZ)
{
	rfNormal.x = Sample(fX - 0.01, fY, fZ, radius) - Sample(fX + 0.01, fY, fZ, radius);
	rfNormal.y = Sample(fX, fY - 0.01, fZ, radius) - Sample(fX, fY + 0.01, fZ, radius);
	rfNormal.z = Sample(fX, fY, fZ - 0.01, radius) - Sample(fX, fY, fZ + 0.01, radius);
	NormalizeVector(rfNormal, rfNormal);
}

void MarchingCube::NormalizeVector(vec3& rfVectorResult, vec3& rfVectorSource)
{
	GLfloat fOldLength;
	GLfloat fScale;

	fOldLength = sqrtf((rfVectorSource.x * rfVectorSource.x) +
		(rfVectorSource.y * rfVectorSource.y) +
		(rfVectorSource.z * rfVectorSource.z));

	if (fOldLength == 0.0)
	{
		rfVectorResult.x = rfVectorSource.x;
		rfVectorResult.y = rfVectorSource.y;
		rfVectorResult.z = rfVectorSource.z;
	}
	else
	{
		fScale = 1.0 / fOldLength;
		rfVectorResult.x = rfVectorSource.x * fScale;
		rfVectorResult.y = rfVectorSource.y * fScale;
		rfVectorResult.z = rfVectorSource.z * fScale;
	}
}

void MarchingCube::GetMeshData(vector<vec3>& sourcePoints, vector<float>& verticesInfo, bool& provideNormal, bool& provideTex)
{
	sourceData = &sourcePoints;
	//此处先固定为只提供normal而没有tex
	provideNormal = true;
	provideTex = false;
	GLint iX, iY, iZ;
	for (iX = 0; iX < iDataSetSize; iX++)
		for (iY = 0; iY < iDataSetSize; iY++)
			for (iZ = 0; iZ < iDataSetSize; iZ++)
			{
				MarchCube(iX * fStepSize, iY * fStepSize, iZ * fStepSize, fStepSize, verticesInfo);
			}
}