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
	/*
	float A = mQuadric(0,0);
	float B = mQuadric(1,0);
	assert(mQuadric(0,1) == B);
	float C = mQuadric(2,0);
	assert(mQuadric(0,2) == B);
	float D = mQuadric(3,0);
	assert(mQuadric(0,3) == D);
	float E = mQuadric(1,1);
	float F = mQuadric(2,1);
	assert(mQuadric(1,2) == F);
	float G = mQuadric(3,1);
	assert(mQuadric(1,3) == G);
	float H = mQuadric(2,2);
	float I = mQuadric(3,2);
	assert(mQuadric(2,3) == I);
	float J = mQuadric(3,3);

	f32 alternative = x*x*A + 2.0f * x*y*B + 2.0f * x*z*C 
				+ 2.0f*x*D + y*y*E + 2.0f * y*z*F
				+ 2.0f*y*G + z*z*H + 2.0f * z*I
				+ J;

	assert(alternative == result);*/

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

	HAGE::Vector4<f32> v(x,y,z,1.0f);

	return *(::Vector3<float>*)&((h_matrix*v).xyz());
}

