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

void SimpleDecimationMesh::Initialize()
{
	mMeshData.SetCustomWeightingFunction(
		[&] (Edge e) -> bool 
		{ 
			computeCollapse(e);
			return true;
		}
		);
	DecimationMesh::Initialize();
}

void SimpleDecimationMesh::computeCollapse(Edge& e)
{
  // The new vertex position is implicitly stored as the
  // position halfway along the edge. The cost is computed as
  // the vertex-to-vertex distance between the new vertex
  // and the old vertices at the edge's endpoints
	auto pair = mMeshData.GetEdgeVertices(e);
	Vector3& v0 = pair[0]->Position;
	Vector3& v1 = pair[1]->Position;
	e->DecimationPosition = (v0 + v1)*0.5;
	e->Cost = (e->DecimationPosition - v0).length();
}
