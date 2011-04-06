/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Söderström (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "SimpleDecimationMesh.h"

void SimpleDecimationMesh::computeCollapse(EdgeCollapse * collapse)
{
  // The new vertex position is implicitly stored as the
  // position halfway along the edge. The cost is computed as
  // the vertex-to-vertex distance between the new vertex
  // and the old vertices at the edge's endpoints
	auto pair = mMeshData.GetEdgeVertices(mMeshData.GetEdge(collapse->halfEdge));
	Vector3& v0 = pair[0]->Position;
	Vector3& v1 = pair[1]->Position;
  collapse->position = (v0 + v1)*0.5;
  collapse->cost = (collapse->position - v0).length();
}
