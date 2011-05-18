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
#ifndef __operatoradvect_h__
#define __operatoradvect_h__

#include "Util/Stopwatch.h"
#include "Levelset/LevelSetOperator.h"
#include "Math/Function3D.h"
#include "Math/Matrix4x4.h"

/*! \brief A level set operator that does external advection
 *
 * This class implements level set advectionr in an external vector field by the PDE
 *
 *  \f$
 *  \dfrac{\partial \phi}{\partial t} + \mathbf{V}(\mathbf{x})\cdot \nabla \phi = 0
 *  \f$
 */
//! \lab4 Implement advection in external vector field
class OperatorAdvect : public LevelSetOperator
{
protected :
  Function3D<Vector3<float> > * mVectorField;

public :

  OperatorAdvect(LevelSet * LS, Function3D<Vector3<float> > * vf) : LevelSetOperator(LS), mVectorField(vf) { }

  virtual float ComputeTimestep()
  {
    // Compute and return a stable timestep
    // (Hint: Function3D::GetMaxValue())
	  Vector3<float> maxVector = mVectorField->GetMaxValue();

	  return mLS->GetDx() / maxVector.Length();
  }

  virtual void Propagate(float time)
  {
    
	 Stopwatch watch;
	 watch.start();
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

	
	watch.stop();
	std::cerr << "Time for dilate/erode = " << watch.read() << std::endl;
  }

  virtual float Evaluate(unsigned int i, unsigned int j, unsigned int k)
  {
    // Compute the rate of change (dphi/dt)

    // Remember that the point (i,j,k) is given in grid coordinates, while
    // the velocity field used for advection needs to be sampled in
    // world coordinates (x,y,z). You can use LevelSet::TransformGridToWorld()
    // for this task.
	  Vector3<float> world_location(i,j,k);
	  mLS->TransformGridToWorld(world_location[0],world_location[1],world_location[2]);
	  Vector3<float> velocity = mVectorField->GetValue(world_location[0],world_location[1],world_location[2]);
	  
	  float dx,dy,dz;

	  if(velocity[0] < 0.0f)
		  dx = mLS->DiffXp(i,j,k);
	  else
		  dx = mLS->DiffXm(i,j,k);

	  if(velocity[1] < 0.0f)
		  dy = mLS->DiffYp(i,j,k);
	  else
		  dy = mLS->DiffYm(i,j,k);

	  if(velocity[2] < 0.0f)
		  dz = mLS->DiffZp(i,j,k);
	  else
		  dz = mLS->DiffZm(i,j,k);

	  return -(dx*velocity[0]+dy*velocity[1]+dz*velocity[2]);
  }

};

#endif
