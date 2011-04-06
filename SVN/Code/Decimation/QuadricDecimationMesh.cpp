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
  // Allocate memory for the quadric array
	unsigned int numVerts = mMeshData.GetNumVertexIndices();
  mQuadrics.reserve(numVerts);
  std::streamsize width = std::cerr.precision(); // store stream precision
  for (unsigned int i = 0; i < numVerts; i++) {

    // Compute quadric for vertex i here
    mQuadrics.push_back(createQuadricForVert(i));


    // Calculate initial error, should be numerically close to 0

	Vector3 v0 = mMeshData.GetVertex(i)->Position;
    Vector4 v(v0[0],v0[1],v0[2],1);
    Matrix4 m = mQuadrics.back();

    float error = v*(m*v);
    std::cerr << std::scientific << std::setprecision(2) << error << " ";
  }
  std::cerr << std::setprecision(width) << std::fixed; // reset stream precision

  // Run the initialize for the parent class to initialize the edge collapses
  DecimationMesh::Initialize();
}

/*! \lab2 Implement the computeCollapse here */
/*!
 * \param[in,out] collapse The edge collapse object to (re-)compute, DecimationMesh::EdgeCollapse
 */
void QuadricDecimationMesh::computeCollapse(EdgeCollapse * collapse)
{
  // Compute collapse->position and collapse->cost here
  // based on the quadrics at the edge endpoints

  std::cerr << "computeCollapse in QuadricDecimationMesh not implemented.\n";
}

/*! After each edge collapse the vertex properties need to be updated */
void QuadricDecimationMesh::updateVertexProperties(unsigned int ind)
{
  DecimationMesh::updateVertexProperties(ind);
  mQuadrics[ind] = createQuadricForVert(ind);
}

/*!
 * \param[in] indx vertex index, points into HalfEdgeMesh::mVerts
 */
QuadricDecimationMesh::Matrix4 QuadricDecimationMesh::createQuadricForVert(unsigned int indx) const{
  Matrix4 Q(
	  Vector4(0,0,0,0),
	  Vector4(0,0,0,0),
	  Vector4(0,0,0,0),
	  Vector4(0,0,0,0));

  // The quadric for a vertex is the sum of all the quadrics for the adjacent faces
  // Tip: Matrix4x4 has an operator +=
  return Q;
}

/*!
 * \param[in] indx face index, points into HalfEdgeMesh::mFaces
 */
QuadricDecimationMesh::Matrix4 QuadricDecimationMesh::createQuadricForFace(unsigned int indx) const{

  // Calculate the quadric (outer product of plane parameters) for a face
  // here using the formula from Garland and Heckbert
  return Matrix4();
}


void QuadricDecimationMesh::Render()
{
  DecimationMesh::Render();

  glEnable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);

  if (mVisualizationMode == QuadricIsoSurfaces)
    {
      // Apply transform
      glPushMatrix(); // Push modelview matrix onto stack

      // Implement the quadric visualization here
      std::cout << "Quadric visualization not implemented" << std::endl;

      // Restore modelview matrix
      glPopMatrix();
    }
}

