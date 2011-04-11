/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sderstrm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "QuadricDecimationMesh.h"

const QuadricDecimationMesh::VisualizationMode QuadricDecimationMesh::QuadricIsoSurfaces = NewVisualizationMode("Quadric Iso Surfaces");

void QuadricDecimationMesh::Initialize()
{
	DecimationMesh::Initialize();
	/*
	HalfEdgeMesh::Initialize();

  // Allocate memory for the quadric array
	unsigned int numVerts = (unsigned int) mMeshData.GetNumVertexIndices();
  std::streamsize width = std::cerr.precision(); // store stream precision
  for (unsigned int i = 0; i < numVerts; i++) {

    // Compute quadric for vertex i here
	  Vertex vert = mMeshData.GetVertex(i);

	  if(vert == mMeshData.nullVertex)
		  continue;

	  vert->Quadric = createQuadricForVert(vert);


    // Calculate initial error, should be numerically close to 0

	Vector3 v0 = mMeshData.GetVertex(i)->Position;
    Vector4 v(v0[0],v0[1],v0[2],1);
    Matrix4 m = vert->Quadric;

    float error = v*(m*v);
    std::cerr << std::scientific << std::setprecision(2) << error << " ";
  }
  std::cerr << std::setprecision(width) << std::fixed; // reset stream precision

  // Run the initialize for the parent class to initialize the edge collapses
  DecimationMesh::Initialize();*/
}

/*! \lab2 Implement the computeCollapse here */
/*!
 * \param[in,out] collapse The edge collapse object to (re-)compute, DecimationMesh::EdgeCollapse
 */
void QuadricDecimationMesh::computeCollapse(Edge& e)
{
  // Compute collapse->position and collapse->cost here
  // based on the quadrics at the edge endpoints
	/*
	assert(e != HageMesh::nullEdge);

	auto vp = mMeshData.GetEdgeVertices(e);

	Matrix4 QEdge = vp[0]->Quadric + vp[1]->Quadric;
	Matrix4 QEdgeInv = Matrix4(QEdge.Row(0),QEdge.Row(1),QEdge.Row(2),Vector4(0.0f,0.0f,0.0f,1.0f)).Invert();

	if(QEdgeInv.IsNaN())
	{
		Vector4 vM = Vector4((vp[0]->Position + vp[1]->Position)/2.0f,1.0f);
		Vector4 v1 = Vector4(vp[0]->Position,1.0f);
		Vector4 v2 = Vector4(vp[1]->Position,1.0f);

		float eVM = vM*(QEdge*vM);
		float eV1 = v1*(QEdge*v1);
		float eV2 = v2*(QEdge*v2);
		if( eVM < eV1 && eVM < eV2)
		{
			e->DecimationPosition = vM.xyz();
			e->Cost = eVM;
		}
		else if( eV1 < eV2)
		{
			e->DecimationPosition = v1.xyz();
			e->Cost = eV1;
		}
		else
		{
			e->DecimationPosition = v2.xyz();
			e->Cost = eV2;
		}
	}
	else
	{
		Vector4 vIdeal = QEdgeInv*Vector4(0.0f,0.0f,0.0f,1.0f);
		e->DecimationPosition = vIdeal.xyz();
		e->Cost = vIdeal*(QEdge*vIdeal);
	}*/

}

/*! After each edge collapse the vertex properties need to be updated */
void QuadricDecimationMesh::updateVertexProperties(Vertex v)
{/*
	if(v == mMeshData.nullVertex)
		return;

	DecimationMesh::updateVertexProperties(v);
	v->Quadric = createQuadricForVert(v);*/
}

/*!
 * \param[in] indx vertex index, points into HalfEdgeMesh::mVerts
 */
QuadricDecimationMesh::Matrix4 QuadricDecimationMesh::createQuadricForVert(Vertex v) const{
  Matrix4 Q = Matrix4::Zero();
  /*
  Face f = HageMesh::nullFace;
  
  while( (f = mMeshData.GetNextVertexFace(v,f)) != HageMesh::nullFace )
  {
	  Q += createQuadricForFace(f);
  }*/

  // The quadric for a vertex is the sum of all the quadrics for the adjacent faces
  // Tip: Matrix4x4 has an operator +=
  return Q;
}

/*!
 * \param[in] indx face index, points into HalfEdgeMesh::mFaces
 */
QuadricDecimationMesh::Matrix4 QuadricDecimationMesh::createQuadricForFace(Face f) const{

  // Calculate the quadric (outer product of plane parameters) for a face
  // here using the formula from Garland and Heckbert
	Vector4 planeEq = getPlaneEquation(f);
	return Matrix4::OuterProduct(planeEq,planeEq);
}


void QuadricDecimationMesh::Render()
{
  DecimationMesh::Render();

  glEnable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);

  static GLUquadric* quad = gluNewQuadric();

  if (mVisualizationMode == QuadricIsoSurfaces)
    {
      // Apply transform
      glPushMatrix(); // Push modelview matrix onto stack

	  for(HageMesh::IndexType i = 0; i < mMeshData.GetNumVertexIndices(); ++i)
	  {
		  Vertex v = mMeshData.GetVertex(i);
		  if(v == HageMesh::nullVertex)
			  continue;

		  Matrix4 R = v->Quadric.Cholesky();
		  if(R.IsNaN())
			  continue;
		  Matrix4 InvR = R.Invert();
		  if(InvR.IsNaN())
			  continue;

		  glPushMatrix(); 
		  const Vector3 color = Vector3(0.0f,1.0f,0.0f);
		  glColor3f(color.x,color.y,color.z);
		  glMultMatrixf( InvR.Transpose().v );

		  gluSphere(quad,1.0f,10,10);
		  glPopMatrix();
	  }

      // Restore modelview matrix
      glPopMatrix();
    }
}

