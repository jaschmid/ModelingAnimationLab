#ifndef _strange_dubdivmesh_
#define _strange_dubdivmesh_

#include "AdaptiveLoopSubdivisionMesh.h"

class StrangeSubdivisionMesh : public AdaptiveLoopSubdivisionMesh
{
protected:
  bool Subdividable(const Face& f){
    // Every 4th face is not subdividable - kinda strange!
    // Do something more interesting...
    return false;
  }

};

#endif
