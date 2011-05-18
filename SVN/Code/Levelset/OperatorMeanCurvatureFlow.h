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
#ifndef __operatormeancurvatureflow_h__
#define __operatormeancurvatureflow_h__

#include "Levelset/LevelSetOperator.h"

/*! \brief A level set operator that does mean curvature flow.
 *
 * This class implements level set propagation in the normal direction
 * as defined by the mean curvature flow \f$\kappa\f$ in the following PDE
 *
 *  \f[
 *  \dfrac{\partial \phi}{\partial t} + \alpha \kappa|\nabla \phi| = 0
 *  \f]
 */
//! \lab4 Implement mean curvature flow
class OperatorMeanCurvatureFlow : public LevelSetOperator
{
protected:
  //! Scaling parameter, affects time step constraint
  float mAlpha;
public :

  OperatorMeanCurvatureFlow(LevelSet * LS, float alpha=.9f)
    : LevelSetOperator(LS), mAlpha(alpha) { }

  virtual float ComputeTimestep()
  {
    // Compute and return a stable timestep
	  return mLS->GetDx()*mLS->GetDx() / (6.0f*mAlpha);
  }

  virtual void Propagate(float time)
  {
    // Determine timestep for stability
    float dt = ComputeTimestep();

    // Propagate level set with stable timestep dt
    // until requested time is reached
    for (float elapsed = 0; elapsed < time;) {

      if (dt > time-elapsed)
        dt = time-elapsed;
      elapsed += dt;

      //IntegrateEuler(dt);
      IntegrateRungeKutta(dt);
    }
  }


  virtual float Evaluate(unsigned int i, unsigned int j, unsigned int k)
  {
    // Compute the rate of change (dphi/dt)

	float dx = mLS->DiffXpm(i,j,k);
	float dy = mLS->DiffYpm(i,j,k);
	float dz = mLS->DiffZpm(i,j,k);
	float dxy = mLS->Diff2XYpm(i,j,k);
	float dxz = mLS->Diff2ZXpm(i,j,k);
	float dyz = mLS->Diff2YZpm(i,j,k);
	float dxx = mLS->Diff2Xpm(i,j,k);
	float dyy = mLS->Diff2Ypm(i,j,k);
	float dzz = mLS->Diff2Zpm(i,j,k);

	Vector3<float> gradient(dx,dy,dz);
	Vector4<float> ex_gradient(gradient[0],gradient[1],gradient[2],1.0f);

	const float hessian[4][4] = {
		{	dxx,	dxy,	dxz,	0.0f},
		{	dxy,	dyy,	dyz,	0.0f},
		{	dxz,	dyz,	dzz ,	0.0f},
		{	0.0f,	0.0f,	0.0f,	1.0f}
	};

	float hTrace = hessian[0][0] + hessian[1][1] + hessian[2][2];
	//the - 1 is just because we have a 4d vector so the dot product has a 1 at the end we need to get rid of
	//3x3 matrix is all we need here really
	float hMult = ex_gradient*(Matrix4x4<float>(hessian)*ex_gradient) - 1.0f;
	float gradLength = gradient.Length();

	float gradLengthSq = gradient*gradient;
    float curv = -(hMult - gradLengthSq*hTrace)/(2.0f* gradLengthSq*gradLength);

	//float curv2 = ( (dx*dx*(dyy+dzz) - 2.0f* dy*dz*dyz) + (dy*dy*(dxx+dzz) - 2.0f* dx*dz*dxz) + (dz*dz*(dxx+dyy) - 2.0f* dx*dy*dxy) ) / (2.0f*pow((dx*dx+dy*dy+dz*dz),3.0f/2.0f));

	return mAlpha*curv*gradLength;
  }


};

#endif
