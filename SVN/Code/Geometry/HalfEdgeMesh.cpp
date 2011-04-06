/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sˆderstrˆm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "Geometry/HalfEdgeMesh.h"

const unsigned int HalfEdgeMesh::BORDER = (std::numeric_limits<unsigned int>::max)();
const unsigned int HalfEdgeMesh::UNINITIALIZED = (std::numeric_limits<unsigned int>::max)()-1;

const static bool bMergeLoad = true;

HalfEdgeMesh::HalfEdgeMesh()
{
}

HalfEdgeMesh::~HalfEdgeMesh()
{
}

/*! \lab1 Implement the addFace */
/*!
 * \param[in] v1 vertex 1, Vector3<float>
 * \param[in] v2 vertex 2, Vector3<float>
 * \param[in] v3 vertex 3, Vector3<float>
 */
bool HalfEdgeMesh::AddFace(const std::vector<::Vector3<float> > &verts){
  // Add your code here

	Vertex v1 = HageMesh::nullVertex,v2 = HageMesh::nullVertex,v3 = HageMesh::nullVertex;

	if(bMergeLoad)
	{
		static std::map<::Vector3<float>,Vertex> vectorMap;
		auto f1 = vectorMap.find(verts[0]);
		auto f2 = vectorMap.find(verts[1]);
		auto f3 = vectorMap.find(verts[2]);

		if(f1!=vectorMap.end())
			v1 = f1->second;
		else
		{
			v1 = mMeshData.InsertVertex();
			vectorMap.insert(std::pair<::Vector3<float>,Vertex>(verts[0],v1));
		}
		if(f2!=vectorMap.end())
			v2 = f2->second;
		else
		{
			v2 = mMeshData.InsertVertex();
			vectorMap.insert(std::pair<::Vector3<float>,Vertex>(verts[1],v2));
		}
		if(f3!=vectorMap.end())
			v3 = f3->second;
		else
		{
			v3 = mMeshData.InsertVertex();
			vectorMap.insert(std::pair<::Vector3<float>,Vertex>(verts[2],v3));
		}
	}
	else
	{
		v1 = mMeshData.InsertVertex();
		v2 = mMeshData.InsertVertex();
		v3 = mMeshData.InsertVertex();
	}
	
	
	mMeshData.GetVertexData(v1).Position = ToLocal(verts[0]);
	mMeshData.GetVertexData(v2).Position = ToLocal(verts[1]);
	mMeshData.GetVertexData(v3).Position = ToLocal(verts[2]);

	assert(mMeshData.InsertFace(mMeshData.MakeTriple(v1,v2,v3)) != HageMesh::nullFace);

  // Add the vertices of the face/triangle

  // Add all half-edge pairs

  // Connect inner ring

  // Finally, create the face, don't forget to set the normal (which should be normalized)

  // All half-edges share the same left face (previously added)

  // Optionally, track the (outer) boundary half-edges
  // to represent non-closed surfaces
  return true;
}



bool HalfEdgeMesh::AddVertex(const ::Vector3<float>& v, unsigned int &indx)
{
	return false;
}


/*! Proceeds to check if the mesh is valid. All indices are inspected and
 * checked to see that they are initialized. The method checks: mEdges, mFaces and mVerts.
 * Also checks to see if all verts have a neighborhood using the findNeighbourFaces method.
 */
void HalfEdgeMesh::Validate()
{
	mMeshData.DebugValidateMesh();
	
	/*
  std::vector<HalfEdge>::iterator iterEdge = mEdges.begin();
  std::vector<HalfEdge>::iterator iterEdgeEnd = mEdges.end();
  while (iterEdge != iterEdgeEnd) {
    if ((*iterEdge).face == UNINITIALIZED ||
        (*iterEdge).next == UNINITIALIZED ||
        (*iterEdge).pair == UNINITIALIZED ||
        (*iterEdge).prev == UNINITIALIZED ||
        (*iterEdge).vert == UNINITIALIZED)
        std::cerr << "HalfEdge " << iterEdge - mEdges.begin() << " not properly initialized" << std::endl;

    iterEdge++;
  }
  std::cerr << "Done with edge check (checked " << GetNumEdges() << " edges)" << std::endl;

  std::vector<Face>::iterator iterTri = mFaces.begin();
  std::vector<Face>::iterator iterTriEnd = mFaces.end();
  while (iterTri != iterTriEnd) {
    if ((*iterTri).edge == UNINITIALIZED)
        std::cerr << "Tri " << iterTri - mFaces.begin() << " not properly initialized" << std::endl;

    iterTri++;
  }
  std::cerr << "Done with face check (checked " << GetNumFaces() << " faces)" << std::endl;

  std::vector<Vertex>::iterator iterVertex = mVerts.begin();
  std::vector<Vertex>::iterator iterVertexEnd = mVerts.end();
  while (iterVertex != iterVertexEnd) {
    if ((*iterVertex).edge == UNINITIALIZED)
        std::cerr << "Vertex " << iterVertex - mVerts.begin() << " not properly initialized" << std::endl;

    iterVertex++;
  }
  std::cerr << "Done with vertex check (checked " << GetNumVerts() << " vertices)" << std::endl;

  std::cerr << "Looping through triangle neighborhood of each vertex... ";
  iterVertex = mVerts.begin();
  iterVertexEnd = mVerts.end();
  int emptyCount = 0;
  std::vector<unsigned int> problemVerts;
  while (iterVertex != iterVertexEnd) {
    std::vector<unsigned int> foundFaces = FindNeighborFaces(iterVertex - mVerts.begin());
    std::vector<unsigned int> foundVerts = FindNeighborVertices(iterVertex - mVerts.begin());
    if (foundFaces.empty() || foundVerts.empty())
      emptyCount++;
    std::set<unsigned int> uniqueFaces(foundFaces.begin(), foundFaces.end());
    std::set<unsigned int> uniqueVerts(foundVerts.begin(), foundVerts.end());
        if ( foundFaces.size() != uniqueFaces.size() ||
         foundVerts.size() != uniqueVerts.size() )
      problemVerts.push_back(iterVertex - mVerts.begin());
    iterVertex++;
  }
  std::cerr << std::endl << "Done: " << emptyCount << " isolated vertices found" << std::endl;
  if(problemVerts.size()){
    std::cerr << std::endl << "Found " << problemVerts.size() << " duplicate faces in vertices: ";
    std::copy(problemVerts.begin(), problemVerts.end(), std::ostream_iterator<unsigned int>(std::cerr, ", "));
    std::cerr << "\n";
  }
  std::cerr << std::endl << "The mesh has genus " << Genus() << ", and consists of " << Shells() << " shells.\n";*/
}

/*! \lab1 Implement the FindNeighborVertices */
/*! Loops over the neighborhood of a vertex and collects all the vertices sorted counter clockwise.
 * \param [in] vertexIndex  the index to vertex, unsigned int
 * \return a vector containing the indices to all the found vertices.
 */
std::vector<unsigned int> HalfEdgeMesh::FindNeighborVertices(unsigned int vertexIndex) const
{
  std::vector<unsigned int> oneRing;
  // Add your code here

  Vertex v = mMeshData.GetVertex(vertexIndex);

  Edge e = mMeshData.GetFirstVertexEdge(v);

  while( e != HageMesh::nullEdge )
  {
	  Vertex neighbor = mMeshData.GetVertex(v,e);

	  oneRing.push_back( (unsigned int)mMeshData.GetIndex(neighbor) );

	  assert( mMeshData.GetEdge(mMeshData.MakePair(v,neighbor)) != HageMesh::nullEdge);

	  e = mMeshData.GetNextVertexEdge(v,e);
  }

  return oneRing;
}

/*! \lab1 Implement the FindNeighborFaces */
/*! Loops over the neighborhood of a vertex and collects all the faces sorted counter clockwise.
 * \param [in] vertexIndex  the index to vertex, unsigned int
 * \return a vector containing the indices to all the found faces.
 */
std::vector<unsigned int> HalfEdgeMesh::FindNeighborFaces(unsigned int vertexIndex) const
{
  std::vector<unsigned int> foundFaces;

  // Add your code here
  return foundFaces;
}


/*! \lab1 Implement the curvature */
float HalfEdgeMesh::VertexCurvature(unsigned int vertexIndex) const
{
  // Copy code from SimpleMesh or compute more accurate estimate
    std::vector<unsigned int> oneRing = FindNeighborVertices(vertexIndex);
  assert(oneRing.size() != 0);

  unsigned int curr, next;
  const Vector3 &vi = mMeshData.GetVertex(vertexIndex)->Position;
  float angleSum = 0;
  float area = 0;
  for(unsigned int i=0; i<oneRing.size(); i++){
    // connections
    curr = oneRing.at(i);
    if(i < oneRing.size() - 1 )
      next = oneRing.at(i+1);
    else
      next = oneRing.front();

    // find vertices in 1-ring according to figure 5 in lab text
    // next - beta
    const Vector3 &nextPos = mMeshData.GetVertex(next)->Position;
    const Vector3 &vj = mMeshData.GetVertex(curr)->Position;

    angleSum +=  acos( (vj-vi)*(nextPos-vi) / ( (vj-vi).length()*(nextPos-vi).length() ) );

    // compute areas
    area += ((vj-vi) % (nextPos-vi)).length()*.5;
  }

  return ( 2*M_PI - angleSum ) / area;
}

float HalfEdgeMesh::FaceCurvature(unsigned int faceIndex) const
{
  // NB Assumes vertex curvature already computed

	HageMesh::VertexTriple vertices = mMeshData.GetFaceVertices(mMeshData.GetFace(faceIndex));


  return (vertices[0]->Curvature + vertices[1]->Curvature + vertices[2]->Curvature) / 3.f;
}

Vector3<float> HalfEdgeMesh::FaceNormal(unsigned int faceIndex) const
{

	HageMesh::VertexTriple vertices = mMeshData.GetFaceVertices(mMeshData.GetFace(faceIndex));

  const Vector3 &p1 = vertices[0]->Position;
  const Vector3 &p2 = vertices[1]->Position;
  const Vector3 &p3 = vertices[2]->Position;

  const Vector3 e1 = p2-p1;
  const Vector3 e2 = p3-p1;

  return ToGlobal((e1%e2).normalize());
}

Vector3<float> HalfEdgeMesh::VertexNormal(unsigned int vertexIndex) const
{

	Face current = mMeshData.GetFirstVertexFace(mMeshData.GetVertex(vertexIndex));

	float fNumNormals = 0.0f;
	Vector3 NormalSum(0.0f,0.0f,0.0f);

	while(current != HageMesh::nullFace)
	{
		NormalSum += current->Normal;
		fNumNormals+=1.0f;
		current = mMeshData.GetNextVertexFace(mMeshData.GetVertex(vertexIndex),current);
	}

	assert(fNumNormals != 0.0f);
  // Add your code here
  return ToGlobal(NormalSum * (1.0f/fNumNormals));
}



void HalfEdgeMesh::mergeMeshVertices()
{
	Vector3 first = mMeshData.GetVertexData(mMeshData.GetVertex(0)).Position;
	Vector3 min=first,max=first;
	for(int i = 1; i < mMeshData.GetNumVertexIndices(); ++i)
	{
		Vector3 v = mMeshData.GetVertexData(mMeshData.GetVertex(i)).Position;
		if(v.x < min.x)
			min.x = v.x;
		if(v.x > max.x)
			max.x = v.x;
		if(v.y < min.y)
			min.y = v.y;
		if(v.y > max.y)
			max.y = v.y;
		if(v.z < min.z)
			min.z = v.z;
		if(v.z > max.z)
			max.z = v.z;
	}

	typedef HAGE::SpatialTree<HageMesh::Vertex> myTree;

	myTree tree(myTree::SplitterType(min,max));

	for(int i = 0; i < mMeshData.GetNumVertexIndices(); ++i)
	{
		HageMesh::Vertex v = mMeshData.GetVertex(i);
		if(v!=HageMesh::nullVertex)
			tree.Insert(v);
	}

	const myTree::TreeNode* cur = tree.TopNode();
	while(cur)
	{
		if(cur->IsLeaf())
		{
			auto elements = cur->LeafElements();
			for(size_t i = 0; i < elements.size(); ++i)
				for(size_t i2 = i+1; i2 < elements.size(); ++i2)
				{
					HageMesh::Vertex v1 = elements[i];
					HageMesh::Vertex v2 = elements[i2];
					if(v1 != HageMesh::nullVertex && v2 != HageMesh::nullVertex)
					{
						Vector3 pos1 = mMeshData.GetVertexData(v1).Position;
						Vector3 pos2 = mMeshData.GetVertexData(v2).Position;
						float distance_sq = !( pos1 - pos2);
						if(pos1 == pos2)
						{
							//mMeshData.DebugValidateMesh();
							//std::cerr << "Merging " << mMeshData.GetIndex(v1) << " and " << mMeshData.GetIndex(v2) << "\n";
							if(!mMeshData.MergeVertex(mMeshData.MakePair(v1,v2)))
								continue;
							else
								break;
						}
					}
				}
		}
		cur = cur->GetNextNode();
	}
	
	mMeshData.DebugValidateMesh();

	std::cerr << "Num Edges pre Merge:" << mMeshData.GetNumEdgeIndices() << "\n";
	std::cerr << "Num Vertices pre Merge:" << mMeshData.GetNumVertexIndices() << "\n";

	mMeshData.Compact();

	std::cerr << "Num Edges post Merge:" << mMeshData.GetNumEdgeIndices() << "\n";
	std::cerr << "Num Vertices post Merge:" << mMeshData.GetNumVertexIndices() << "\n";
}


void HalfEdgeMesh::Initialize() {

	//Merge Vertices
	mergeMeshVertices();
  Validate();
  Update();
}


void HalfEdgeMesh::Update() {
	
  // Calculate and store all differentials and area

  // First update all face normals and triangle areas
	for(unsigned int i = 0; i < mMeshData.GetNumFaceIndices(); i++){
	  mMeshData.GetFace(i)->Normal = ToLocal(FaceNormal(i));
  }
  // Then update all vertex normals and curvature
  for(unsigned int i = 0; i < mMeshData.GetNumVertexIndices(); i++){
    // Vertex normals are just weighted averages
    mMeshData.GetVertex(i)->Normal = ToLocal(VertexNormal(i));
  }

  // Then update vertex curvature
  for(unsigned int i = 0; i < mMeshData.GetNumVertexIndices(); i++){
    mMeshData.GetVertex(i)->Curvature = VertexCurvature(i);
    //    std::cerr <<   mVerts.at(i).curvature << "\n";
  }

  // Finally update face curvature
  for(unsigned int i = 0; i < mMeshData.GetNumFaceIndices(); i++){
	  mMeshData.GetFace(i)->Curvature = FaceCurvature(i);
  }

  std::cerr << "Area: " << Area() << ".\n";
  std::cerr << "Volume: " << Volume() << ".\n";

  // Update vertex and face colors
  if (mVisualizationMode == CurvatureVertex) {

	  
    float minCurvature = (std::numeric_limits<float>::max)();
	float maxCurvature = -(std::numeric_limits<float>::max)();

	for(unsigned int i = 0; i < mMeshData.GetNumVertexIndices(); i++){
		if (minCurvature > mMeshData.GetVertex(i)->Curvature)  minCurvature = mMeshData.GetVertex(i)->Curvature;
		if (maxCurvature < mMeshData.GetVertex(i)->Curvature)  maxCurvature = mMeshData.GetVertex(i)->Curvature;
	}


    std::cerr << "Mapping color based on vertex curvature with range [" << minCurvature << "," << maxCurvature << "]" << std::endl;
    
	for(unsigned int i = 0; i < mMeshData.GetNumVertexIndices(); i++){
		mMeshData.GetVertex(i)->Color = ToLocal(mColorMap->Map(mMeshData.GetVertex(i)->Curvature, minCurvature, maxCurvature));
	}
  }
  else if (mVisualizationMode == CurvatureFace) {

	  
    float minCurvature = (std::numeric_limits<float>::max)();
	float maxCurvature = -(std::numeric_limits<float>::max)();

	for(unsigned int i = 0; i < mMeshData.GetNumFaceIndices(); i++){
		if (minCurvature > mMeshData.GetFace(i)->Curvature)  minCurvature = mMeshData.GetFace(i)->Curvature;
		if (maxCurvature < mMeshData.GetFace(i)->Curvature)  maxCurvature = mMeshData.GetFace(i)->Curvature;
	}


    std::cerr << "Mapping color based on vertex curvature with range [" << minCurvature << "," << maxCurvature << "]" << std::endl;
    
	for(unsigned int i = 0; i < mMeshData.GetNumFaceIndices(); i++){
		mMeshData.GetFace(i)->Color = ToLocal(mColorMap->Map(mMeshData.GetFace(i)->Curvature, minCurvature, maxCurvature));
	}
  }
}


/*! \lab1 Implement the area */
float HalfEdgeMesh::Area() const
{
  float area = 0;
  // Add code here
  std::cerr << "Area calculation not implemented for half-edge mesh!\n";
  return area;
}

/*! \lab1 Implement the volume */
float HalfEdgeMesh::Volume() const
{
  float volume = 0;
  // Add code here
  std::cerr << "Volume calculation not implemented for half-edge mesh!\n";
  return volume;
}

/*! \lab1 Calculate the number of shells  */
int HalfEdgeMesh::Shells() const
{
  return 1;
}

/*! \lab1 Implement the genus */
int HalfEdgeMesh::Genus() const
{
  // Add code here
  std::cerr << "Genus calculation not implemented for half-edge mesh!\n";
  return 0;
}

void HalfEdgeMesh::Render()
{
  glEnable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);

  // Apply transform
  glPushMatrix(); // Push modelview matrix onto stack

  // Convert transform-matrix to format matching GL matrix format
  // Load transform into modelview matrix
  glMultMatrixf(  mTransform.ToGLMatrix().GetArrayPtr() );

  if (mWireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  // Draw geometry
  glBegin(GL_TRIANGLES);
  const int numTriangles = (int)mMeshData.GetNumFaces();
  for (int i = 0; i < numTriangles; i++){

	Face face = mMeshData.GetFace(i);

	HageMesh::VertexTriple vt = mMeshData.GetFaceVertices(face);

    if (mVisualizationMode == CurvatureVertex) {
		for(int v=0;v<3;v++)
		{
		  glColor3fv(vt[v]->Color.c);
		  glNormal3fv(vt[v]->Normal.c);
		  glVertex3fv(vt[v]->Position.c);
		}
    }
    else {
      glColor3fv(face->Color.c);
      glNormal3fv(face->Normal.c);

		for(int v=0;v<3;v++)
		{
		  glVertex3fv(vt[v]->Position.c);
		}
    }

  }
  glEnd();

  // Mesh normals by courtesy of Richard Khoury
  if (mShowNormals)
  {
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    const int numTriangles = (int)mMeshData.GetNumFaces();
    for (int i = 0; i < numTriangles; i++){

     
	Face face = mMeshData.GetFace(i);
	HageMesh::VertexTriple vt = mMeshData.GetFaceVertices(face);


	Vector3 faceStart = (vt[0]->Position + vt[1]->Position + vt[2]->Position) / 3.0f;
      Vector3 faceEnd = faceStart + face->Normal*0.1;

      glColor3f(1,0,0); // Red for face normal
      glVertex3fv(faceStart.c);
      glVertex3fv(faceEnd.c);

      glColor3f(0,1,0); // Vertex normals in Green
	  
		for(int v=0;v<3;v++)
		{
		  glVertex3fv(vt[v]->Position.c);
		  glVertex3fv((vt[v]->Position + vt[v]->Normal*0.1f).c);
		}
    }
    glEnd();
    glEnable(GL_LIGHTING);

  }


  if (mWireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // Restore modelview matrix
  glPopMatrix();

  GLObject::Render();

}

