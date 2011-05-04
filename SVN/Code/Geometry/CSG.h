
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
#ifndef __CSG_H__
#define __CSG_H__

#include "Geometry/Implicit.h"


/*! \brief CSG Operator base class */
class CSG_Operator : public Implicit
{

protected:
  //! Constructor
  CSG_Operator(Implicit * l, Implicit * r) : left(l), right(r) { }

  //! Pointers to left and right child nodes
  Implicit * left, * right;
};


/*! \brief Union boolean operation */
class Union : public CSG_Operator
{
public:
  Union(Implicit * l, Implicit * r) : CSG_Operator(l, r) {
    // Compute the resulting (axis aligned) bounding box from
    // the left and right children
    mBox = BoxUnion(l->GetBoundingBox(), r->GetBoundingBox());
  }

  virtual float GetValue(float x, float y, float z) const {
    // The coordinates (x,y,z) are passed in from world space,
    // remember to transform them into object space
    // (Hint: Implicit::TransformW2O()). This
    // is needed because the CSG operators are also implicit geometry
    // and can be transformed like all implicit surfaces.
    // Then, get values from left and right children and perform the
    // boolean operation.
	Implicit::TransformW2O(x,y,z);
	float l = left->GetValue(x,y,z);
	float r = right->GetValue(x,y,z);
    return std::min<float>(l,r);
  }
};


/*! \brief Intersection boolean operation */
class Intersection : public CSG_Operator
{
public:
  Intersection(Implicit * l, Implicit * r) : CSG_Operator(l, r) {
    mBox = BoxIntersection(l->GetBoundingBox(), r->GetBoundingBox());
  }

  virtual float GetValue(float x, float y, float z) const {
	Implicit::TransformW2O(x,y,z);
	float l = left->GetValue(x,y,z);
	float r = right->GetValue(x,y,z);
    return std::max<float>(l,r);
  }
};


/*! \brief Difference boolean operation */
class Difference : public CSG_Operator
{
public:
  Difference(Implicit * l, Implicit * r) : CSG_Operator(l, r) {
    mBox = l->GetBoundingBox();
  }

  virtual float GetValue(float x, float y, float z) const {
	Implicit::TransformW2O(x,y,z);
	float l = left->GetValue(x,y,z);
	float r = right->GetValue(x,y,z);
    return std::max<float>(l,-r);
  }
};


/*! \brief BlendedUnion boolean operation */
class BlendedUnion : public CSG_Operator
{
public:
  BlendedUnion(Implicit * l, Implicit * r, int blend) : CSG_Operator(l, r), mBlend(blend) {
    mBox = BoxUnion(l->GetBoundingBox(), r->GetBoundingBox());
  }

  virtual float GetValue(float x, float y, float z) const {
	Implicit::TransformW2O(x,y,z);
	float l = exp(-left->GetValue(x,y,z));
	float r = exp(-right->GetValue(x,y,z));
	float res = pow(pow(l,mBlend) + pow(r,mBlend),1.0f/(float)mBlend);
	return -log(res);
  }

protected :
  int mBlend;
};


/*! \brief BlendedIntersection boolean operation */
class BlendedIntersection : public CSG_Operator
{
public:
  BlendedIntersection(Implicit * l, Implicit * r, int blend) : CSG_Operator(l, r), mBlend(blend) {
    mBox = BoxUnion(l->GetBoundingBox(), r->GetBoundingBox());
  }

  virtual float GetValue(float x, float y, float z) const {
	Implicit::TransformW2O(x,y,z);
	float l = exp(-left->GetValue(x,y,z));
	float r = exp(-right->GetValue(x,y,z));
	float res= pow(pow(l,-mBlend) + pow(r,-mBlend),-1.0f/(float)mBlend);
	return -log(res);
  }

protected :
  int mBlend;
};


/*! \brief BlendedDifference boolean operation */
class BlendedDifference : public CSG_Operator
{
public:
  BlendedDifference(Implicit * l, Implicit * r, int blend) : CSG_Operator(l, r), mBlend(blend) {
    mBox = l->GetBoundingBox();
  }

  virtual float GetValue(float x, float y, float z) const {
	Implicit::TransformW2O(x,y,z);
	float l = exp(-left->GetValue(x,y,z));
	float r = exp(right->GetValue(x,y,z));
	float res = pow(pow(l,-mBlend) + pow(r,-mBlend),-1.0f/(float)mBlend);
	return -log(res);
  }

protected :
  int mBlend;
};


#endif
