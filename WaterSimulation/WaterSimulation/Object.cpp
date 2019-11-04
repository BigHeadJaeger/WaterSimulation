#include "Object.h"

void MeshObject::GetVertexDataArray(vector<float>& data)
{
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		for (Mesh::FaceVertexIter fv_ccwit = mesh.fv_iter(*f_it); fv_ccwit.is_valid(); fv_ccwit++)
		{
			data.push_back(mesh.point(*fv_ccwit).data()[0]);
			data.push_back(mesh.point(*fv_ccwit).data()[1]);
			data.push_back(mesh.point(*fv_ccwit).data()[2]);
			data.push_back(mesh.normal(*fv_ccwit).data()[0]);
			data.push_back(mesh.normal(*fv_ccwit).data()[1]);
			data.push_back(mesh.normal(*fv_ccwit).data()[2]);
			data.push_back(mesh.texcoord2D(*fv_ccwit).data()[0]);
			data.push_back(mesh.texcoord2D(*fv_ccwit).data()[1]);
		}
	}
}

void MeshObject::readObjFile(string fileName)
{
	mesh.request_vertex_normals();


	OpenMesh::IO::Options opt;
	if (!OpenMesh::IO::read_mesh(mesh, fileName), opt)
	{
		cout << "init mesh failed" << endl;
	}

	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		mesh.request_face_normals();
		mesh.update_normals();
		mesh.release_face_normals();

	}

	mesh.request_vertex_texcoords2D();

	//UE4ShaderData* temp = dynamic_cast<UE4ShaderData*>(shaderData);
	//if (temp == NULL)
	//{
	//	cout << "shaderData ptr convert to UE4ShaderData fail" << endl;
	//	return;
	//}
	//if (!opt.check(OpenMesh::IO::Options::VertexTexCoord))
	//	temp->bUseTexture = false;
	//else
	//{
	//	temp->bUseTexture = true;
	//	mesh.request_vertex_texcoords2D();
	//}

}

void MeshObject::InitBufferData()
{
	shaderData->drawUnitNumber = mesh.n_faces() * 3;
	vector<float> data;
	GetVertexDataArray(data);
	shaderData->InitVertexBuffer(data, true, false);
}

void MeshObject::Update(float dt)
{
	shaderData->UpdateMatrix(transformation);
}

void MeshObject::Draw()
{
	renderer->Render(shaderData);
}

void Object::SetRenderer(RENDERERTYPE type)
{
	switch (type)
	{
	case UE4RENDERER:
		renderer = UE4Renderer::GetRenderer();
		delete shaderData;
		shaderData = new UE4ShaderData();
		break;
	case SIMPLERENDER:
		renderer = SimpleRenderer::GetRenderer();
		delete shaderData;
		shaderData = new SimpleShaderData();
		break;
	case MPSRENDERER:
		break;
	default:
		break;
	}
}

void Metaball::SetSourcePoints(vec3 firstPos, int w, int h, int d)
{
	initPos = firstPos;
	CubeDistribute(sourcePoints, firstPos, vec3(0.2), w, h, d);
}

void Metaball::InitBufferData()
{
	int pointCount = 0;
	vector<float> verticesInfo;
	bool provideNormal;
	bool provideTex;
	marchingCube.GetMeshData(sourcePoints, verticesInfo, provideNormal, provideTex);
	//当ball运动都超出最开始的定的边界时就不会有顶点信息产生
	if (verticesInfo.size() != 0)
	{
		pointCount = verticesInfo.size() / 8;
		shaderData->drawUnitNumber = pointCount;
		shaderData->InitVertexBuffer(verticesInfo, provideNormal, provideTex);
	}
	//pointCount = verticesInfo.size() / 8;
	//shaderData->drawUnitNumber = pointCount;
	//shaderData->InitVertexBuffer(verticesInfo, provideNormal, provideTex);
}

void Metaball::Update(float dt)
{
	cout << 1 / dt << endl;
	temp += dt;
	//对sourcePoint进行坐标变化
	float fOffset = 1.0 + sinf(temp);
	//cout << fOffset << endl;

	for (int i = 0; i < sourcePoints.size(); i++)
	{
		int a = i % 3;
		switch (a)
		{
		case 0:
			sourcePoints[i].x = initPos.x * fOffset;
			//cout << sourcePoints[i].x << endl;
			break;
		case 1:
			sourcePoints[i].y = initPos.y * fOffset;
			break;
		case 2:
			sourcePoints[i].z = initPos.z * fOffset;
			break;
		default:
			break;
		}
	}

	//每一帧都需要重新初始化顶点buffer
	//if (test % 2 == 0)
	//{
	//	
	//}
	InitBufferData();
	shaderData->UpdateMatrix(transformation);
}

void Metaball::Draw()
{
	renderer->Render(shaderData);
}
