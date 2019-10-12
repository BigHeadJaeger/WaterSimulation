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


	UE4ShaderData* temp = dynamic_cast<UE4ShaderData*>(shaderData);
	if (temp == NULL)
	{
		cout << "shaderData ptr convert to UE4ShaderData fail" << endl;
		return;
	}
	if (!opt.check(OpenMesh::IO::Options::VertexTexCoord))
		temp->bUseTexture = false;
	else
	{
		temp->bUseTexture = true;
		mesh.request_vertex_texcoords2D();
	}

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
