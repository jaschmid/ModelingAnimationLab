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
#include "Decimation/DecimationMesh.h"
#include "Util/ColorMap.h"
#include "GUI/GUI.h"
#include <cassert>


const DecimationMesh::VisualizationMode DecimationMesh::CollapseCost = NewVisualizationMode("Collapse cost");


void DecimationMesh::Initialize()
{

  // Allocate memory for the references from half-edge
  // to edge collapses
  mHalfEdge2EdgeCollapse.reserve(mMeshData.GetNumEdgeIndices());
  mHalfEdge2EdgeCollapse.assign(mMeshData.GetNumEdgeIndices(), NULL);

  // Loop through the half-edges (we know they are stored
  // sequentially) and create an edge collapse operation
  // for each pair
  unsigned int numCollapses = mMeshData.GetNumEdgeIndices();
  for (unsigned int i = 0; i < numCollapses; i++) {
    EdgeCollapse * collapse = new EdgeCollapse();

    // Connect the edge collapse with the half-edge pair
    collapse->halfEdge = i;

    mHalfEdge2EdgeCollapse[i] = collapse;

    // Compute the cost and push it to the heap
    computeCollapse(collapse);
    mHeap.push(collapse);
  }
  //mHeap.print(std::cout);

  HalfEdgeMesh::Initialize();
}


void DecimationMesh::Update()
{
  // Calculate and store all differentials and area

  // First update all face normals and triangle areas
  for(unsigned int i = 0; i < mMeshData.GetNumFaceIndices(); i++){
    if (mMeshData.GetFace(i) != HageMesh::nullFace)
		mMeshData.GetFace(i)->Normal = ToLocal(FaceNormal(i));
  }
  // Then update all vertex normals and curvature
  for(unsigned int i = 0; i < mMeshData.GetNumVertexIndices(); i++){
    // Vertex normals are just weighted averages
    if (mMeshData.GetVertex(i) != HageMesh::nullVertex)
      mMeshData.GetVertex(i)->Normal = ToLocal(VertexNormal(i));
  }

  // Then update vertex curvature
  for(unsigned int i = 0; i < mMeshData.GetNumVertexIndices(); i++){
    if (mMeshData.GetVertex(i) != HageMesh::nullVertex)
      mMeshData.GetVertex(i)->Curvature = VertexCurvature(i);
    //    std::cerr <<   mVerts.at(i).curvature << "\n";
  }

  // Finally update face curvature
  for(unsigned int i = 0; i < mMeshData.GetNumFaceIndices(); i++){
    if (mMeshData.GetFace(i) != HageMesh::nullFace)
      mMeshData.GetFace(i)->Curvature = FaceCurvature(i);
  }

//  std::cerr << "Area: " << Area() << ".\n";
//  std::cerr << "Volume: " << Volume() << ".\n";

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

bool DecimationMesh::decimate(unsigned int targetFaces)
{
  // We can't collapse down to less than two faces
  if (targetFaces < 2) targetFaces = 2;

  // Keep collapsing one edge at a time until the target is reached
  // or the heap is empty (when we have no possible collapses left)
  while (mMeshData.GetNumFaces()  > targetFaces && !mHeap.isEmpty())
    decimate();

  // Return true if target is reached
  std::cout << "Collapsed mesh to " << mMeshData.GetNumFaces() << " faces" << std::endl;
  return mMeshData.GetNumFaces()  == targetFaces;
}


bool DecimationMesh::decimate()
{
  EdgeCollapse * collapse = static_cast<EdgeCollapse *>(mHeap.pop());
  if (collapse == NULL) return false;

  // Stop the collapse when we only have two triangles left
  // (the smallest entity representable)
  if (mMeshData.GetNumFaces()  == 2) return false;

  
  Vertex v = HageMesh::nullVertex;

  while(mHeap.size() > 0)
  {
	  Edge e = mMeshData.GetEdge(collapse->halfEdge);
	  auto vp = mMeshData.GetEdgeVertices(e);

	  if(e != HageMesh::nullEdge && mMeshData.MergeVertex(vp))
	  {
		  v=vp[1];
		  v->Position = collapse->position;
		  //delete collapse;
		  return true;
	  }
	  else
	  {
		  auto last = collapse;
		collapse = static_cast<EdgeCollapse *>(mHeap.pop());
		//delete last;
		if (collapse == NULL) return false;
	  }
  }

  //mHeap.print(std::cout);

  return true;
}


void DecimationMesh::updateVertexProperties(unsigned int ind)
{
  std::vector<unsigned int> neighborFaces = FindNeighborFaces(ind);

  // Approximate vertex normal
  Vector3 n(0,0,0);

  const unsigned int numCandidates = neighborFaces.size();
  for (unsigned int i = 0; i < numCandidates; i++){

    n += ToLocal(FaceNormal(neighborFaces[i]));
  }

  mMeshData.GetVertex(ind)->Normal = n.normalize();
}


void DecimationMesh::updateFaceProperties(unsigned int ind)
{

  mMeshData.GetFace(ind)->Normal = ToLocal(FaceNormal(ind));
}


void DecimationMesh::Cleanup()
{
//  HalfEdgeMesh mesh;
//  *this = mesh;
}


void DecimationMesh::Render()
{
	/*
  glEnable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);

  // Apply transform
  glPushMatrix(); // Push modelview matrix onto stack

  // Convert transform-matrix to format matching GL matrix format
  // Load transform into modelview matrix
  glMultMatrixf(mTransform.ToGLMatrix().GetArrayPtr());

  if (mWireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  // Draw geometry
  const unsigned int numTriangles = mMeshData.GetNumFaceIndices();
  for (unsigned int i = 0; i < numTriangles; i++){

    if (isFaceCollapsed(i)) continue;

    // Render without notations
    Face & f = mFaces[i];

    HalfEdge* edge = &mEdges[f.edge];

    Vertex & v1 = mVerts[edge->vert];
    edge = &mEdges[edge->next];

    Vertex & v2 = mVerts[edge->vert];
    edge = &mEdges[edge->next];

    Vertex & v3 = mVerts[edge->vert];

    // Render with notations
    //  Uncomment this block, and comment the block above
    //    to render with notations
    /*
        Face & f = mFaces[i];

        char buffer[10];
        glColor3f(1.0, 0.0, 0.0);
        HalfEdge* edge = &mEdges[mFaces[i].edge];

        // draw face
        sprintf(buffer, "f%i\n", i);
        Vector3<float> vec = (mVerts[edge->vert].pos + mVerts[mEdges[edge->pair].vert].pos)*0.5;
        vec += 0.5 * (mVerts[mEdges[edge->prev].vert].pos - vec);
        drawText(vec, buffer);

        // draw e1
        sprintf(buffer, "e%i\n", mFaces[i].edge);
        vec = (mVerts[edge->vert].pos + mVerts[mEdges[edge->pair].vert].pos)*0.5;
        vec += 0.1 * (mVerts[mEdges[edge->prev].vert].pos - vec);
        drawText(vec, buffer);

        // draw v1
        Vertex& v1 = mVerts[edge->vert];
        sprintf(buffer, "v%i\n", edge->vert);
        drawText(vec, buffer);

        sprintf(buffer, "e%i\n", edge->next);
        edge = &mEdges[edge->next];

        // draw e2
        vec = (mVerts[edge->vert].pos + mVerts[mEdges[edge->pair].vert].pos)*0.5;
        vec += 0.1 * (mVerts[mEdges[edge->prev].vert].pos - vec);
        drawText(vec, buffer);

        // draw v2
        Vertex& v2 = mVerts[edge->vert];
        sprintf(buffer, "v%i\n", edge->vert);
        drawText(vec, buffer);

        sprintf(buffer, "e%i\n", edge->next);
        edge = &mEdges[edge->next];

        // draw e3
        vec = (mVerts[edge->vert].pos + mVerts[mEdges[edge->pair].vert].pos)*0.5;
        vec += 0.1 * (mVerts[mEdges[edge->prev].vert].pos - vec);
        drawText(vec, buffer);

        // draw v3
        Vertex& v3 = mVerts[edge->vert];
        sprintf(buffer, "v%i\n", edge->vert);
        drawText(vec, buffer);
     */
	
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
  const int numTriangles = (int)mMeshData.GetNumFaceIndices();
  for (int i = 0; i < numTriangles; i++){
	  
	Face face = mMeshData.GetFace(i);
	if(face == HageMesh::nullFace)
		continue;

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
    const int numTriangles = (int)mMeshData.GetNumFaceIndices();
    for (int i = 0; i < numTriangles; i++){

     
	Face face = mMeshData.GetFace(i);
	if(face == HageMesh::nullFace)
		continue;
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

void DecimationMesh::drawText(const Vector3 & pos, const char * str)
{
  glRasterPos3f(pos[0], pos[1], pos[2]);
  for (unsigned int i = 0; str[i] != '\n'; i++)
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, str[i]);
}
