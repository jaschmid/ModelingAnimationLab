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
#define HAGE_NO_MEMORY_ALLOC_REPLACEMENT
#include <HAGE.h>

#include "Quadric.h"


Quadric::Quadric(const Matrix4x4<float> & q){
  this->mQuadric = q;
}

Quadric::~Quadric(){}

/*!
 * Evaluation of world coordinates are done through either transformation
 * of the world-coordinates by mWorld2Obj, or transformation of the quadric
 * coefficient matrix by GetTransform() ONCE (see Section 2.2 in lab text).
 */
float Quadric::GetValue(float x, float y, float z) const
{
	Implicit::TransformW2O(x,y,z);

	static_assert(sizeof(HAGE::Matrix4<HAGE::f32>) == sizeof(Matrix4x4<float>),"Matrix sizes must be identical");
	static_assert(sizeof(HAGE::Vector3<HAGE::f32>) == sizeof(Vector3<float>),"Vector sizes must be identical");
	using namespace HAGE;
	Matrix4<f32>& h_matrix = *(Matrix4<f32>*)&mQuadric;

	HAGE::Vector4<f32> v(x,y,z,1.0f);

	f32 result = v*(h_matrix*v);

	return result;
}

/*!
 * Use the quadric matrix to evaluate the gradient.
 */

Vector3<float> Quadric::GetGradient(float x, float y, float z) const
{
	
	static_assert(sizeof(HAGE::Matrix4<HAGE::f32>) == sizeof(Matrix4x4<float>),"Matrix sizes must be identical");
	static_assert(sizeof(HAGE::Vector3<HAGE::f32>) == sizeof(Vector3<float>),"Vector sizes must be identical");
	using namespace HAGE;
	Matrix4<f32>& h_matrix = *(Matrix4<f32>*)&mQuadric;

	HAGE::Vector4<f32> v(x,y,z,0.0f);

	return *(::Vector3<float>*)&((h_matrix*v).xyz());
}

