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
#ifndef __HALF_EDGE_H__
#define __HALF_EDGE_H__

#ifdef S_OK
#undef S_OK
#endif

#ifdef E_FAIL
#undef E_FAIL
#endif

#ifdef Message
#undef Message
#endif

#ifdef MessageHelper
#undef MessageHelper
#endif

#define HAGE_NO_MEMORY_ALLOC_REPLACEMENT
#include <HAGE.h>

#include "Geometry/Mesh.h"
#include "Util/ObjIO.h"
#include <map>
#include <limits>
#include <cassert>
#include <set>

/*! \brief A half edge triangle mesh class.
 * A mesh data structure with fast access to neighborhood entities.
 * The neighborhood links are stored through indices.
 */

class HalfEdgeMesh : public Mesh {
private:

public:
  HalfEdgeMesh();
  ~HalfEdgeMesh();
  virtual void Update();
  virtual void Initialize();
  
  //! Adds a triangle to the mesh.
  virtual bool AddFace(const std::vector<Vector3<float> > &verts);

  //! Calculates the area of the mesh
  virtual float Area() const;

  //! Calculates the volume of the mesh
  virtual float Volume() const;

  //! Calculates the genus of a mesh
  virtual int Genus() const;

  //! Calculates the number of shells
  virtual int Shells() const;

  //! Calculates the curvature at a vertex
  virtual float VertexCurvature(unsigned int vertexIndex) const;

  //! Calculates the curvature at a vertex
  virtual float FaceCurvature(unsigned int faceIndex) const;

  //! Calculates the normal at a face
  virtual Vector3<float> FaceNormal(unsigned int faceIndex) const;

  //! Calculates the normal at a vertex
  virtual Vector3<float> VertexNormal(unsigned int vertexIndex) const;

  //! Checks to see if the mesh is valid
  void Validate();

  virtual void Render();

protected:
	
  //! Adds a vertex to the mesh
  virtual bool AddVertex(const Vector3<float>& v, unsigned int &indx) ;

  void mergeMeshVertices();

  //! Denotes a reference to a border, only for face pointers
  const static unsigned int BORDER;
  //! Denotes a reference to a non-existing object
  const static unsigned int UNINITIALIZED;

  /*! \brief The core half-edge struct
   *  Implements the linked data structure edge type
   */

  typedef HAGE::Vector3<> Vector3;
  

  struct VertexData
  {
	  Vector3 Normal;
	  Vector3 Position;
	  Vector3 Color;
	  float Curvature;

	  
	operator Vector3&()
	{
		return Position;
	}
	operator const Vector3&() const
	{
		return Position;
	}
  };

  struct FaceData
  {
	  Vector3 Normal;
	  Vector3 Color;
	  float Curvature;
  };

  ::Vector3<float> ToGlobal(const Vector3& v) const
  {
	  return ::Vector3<float>(v.x,v.y,v.z);
  }

  Vector3 ToLocal(const ::Vector3<float>& v) const
  {
	  return Vector3(v[0],v[1],v[2]);
  }

  typedef HAGE::CEditableMesh<VertexData,FaceData,void> HageMesh;
  typedef HageMesh::Vertex Vertex;
  typedef HageMesh::Face Face;
  typedef HageMesh::Edge Edge;

  HageMesh mMeshData;
  

  float FaceArea(HageMesh::IndexType i) const;
  unsigned int GetNumEdges() const { return (unsigned int)mMeshData.GetNumEdgeIndices(); }
  
  
  //! Finds all triangles that includes a given vertex.
  virtual std::vector<unsigned int> FindNeighborFaces(unsigned int vertexIndex) const;

  //! Finds all vertices that includes a given vertex.
  virtual std::vector<unsigned int> FindNeighborVertices(unsigned int vertexIndex) const;

  /*
  struct HalfEdge {
    HalfEdge() : vert(UNINITIALIZED), face(UNINITIALIZED), next(UNINITIALIZED),
                 prev(UNINITIALIZED), pair(UNINITIALIZED) { }
    unsigned int vert;  //!< index into mVerts (the origin vertex)
    unsigned int face;  //!< index into mFaces
    unsigned int next;  //!< index into mEdges
    unsigned int prev;  //!< index into mEdges
    unsigned int pair;  //!< index into mEdges
  };*/

  /*! \brief A vertex is a point and an edge index
   *//*
  struct Vertex : public Mesh::Vertex {
    Vertex() : edge(UNINITIALIZED) { }
    unsigned int edge;     //!< index into mEdges
  };
  */
  /*! \brief A face has a normal and an index to an edge
   *//*
  struct Face : public Mesh::Face {
    Face() : edge(UNINITIALIZED) { }
    unsigned int edge; //!< index into mEdges
  };

   struct OrderedPair{
    unsigned int p1, p2;
    OrderedPair(unsigned int i1, unsigned int i2) {
      p1 = std::min(i1, i2);
      p2 = std::max(i1, i2);
    }
    bool operator < (const OrderedPair & rhs) const{
      if (this->p1 < rhs.p1){
        return true;
      } else if(this->p1 == rhs.p1){
        if(this->p2 < rhs.p2){
          return true;
        }
      }
      return false;
    }
  };

  //! Return the edge at index i
  HalfEdge& e(unsigned int i) { return mEdges.at(i); }
  const HalfEdge& e(unsigned int i) const { return mEdges.at(i); }
  //! Return the face at index i
  Face& f(unsigned int i) { return mFaces.at(i); }
  const Face& f(unsigned int i) const { return mFaces.at(i); }
  //! Return the Vertex at index i
  Vertex& v(unsigned int i) { return mVerts.at(i); }
  const Vertex v(unsigned int i) const { return mVerts.at(i); }
  //! Return number of vertices
  unsigned int GetNumVerts() const { return mVerts.size(); }
  //! Return number of faces
  unsigned int GetNumFaces() const { return mFaces.size(); }
  //! Return number of edges
  unsigned int GetNumEdges() const { return mEdges.size(); }
  */
  virtual void Dilate(float amount);
  virtual void Erode(float amount);
  virtual void Smooth(float amount);

  virtual bool save(std::ostream &os){
    os << "# HalfEdgeMesh obj streamer\n# M&A 2008\n\n";
    os << "# Vertices\n";
	for(unsigned int i=0; i<mMeshData.GetNumVertexIndices(); i++){
		os << "v " << ToGlobal(mMeshData.GetVertex(i)->Position) << " " <<  ToGlobal(mMeshData.GetVertex(i)->Position) << " " <<  ToGlobal(mMeshData.GetVertex(i)->Position) << "\n";
    }
    os << "\n# Faces\n";
    for(unsigned int i=0; i<mMeshData.GetNumFaceIndices(); i++){
		HageMesh::VertexTriple vt = mMeshData.GetFaceVertices(mMeshData.GetFace(i));
      os << "f " << mMeshData.GetIndex(vt[0])+1 << " "
         << mMeshData.GetIndex(vt[1])+1
         << " " <<  mMeshData.GetIndex(vt[2])+1 << "\n";
    }
    return os.good();
  }
};

#endif
