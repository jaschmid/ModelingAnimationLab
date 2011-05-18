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
#ifndef __operatordilateerode_h__
#define __operatordilateerode_h__

#include "Util/Stopwatch.h"
#include "Levelset/LevelSetOperator.h"

/*! \brief A level set operator that does erosion or dilation.
 *
 * This class implements level set propagation in the normal direction
 * as defined by
 *  \f$
 *  \dfrac{\partial \phi}{\partial t}+ F(\mathbf{x})|\nabla \phi| = 0
 *  \f$
 * where the sign of F dictates erosion (c<0), or dilation (c>0).
 */
//! \lab4 Implement erosion and dilation
class OperatorDilateErode : public LevelSetOperator
{
protected :
  //! The constant speed function
  float mF;

public :

  OperatorDilateErode(LevelSet * LS, float f) : LevelSetOperator(LS), mF(f) { }

  virtual float ComputeTimestep()
  {
    // Compute and return a stable timestep
	  
	  return mLS->GetDx() / fabs (mF);
  }

  virtual void Propagate(float time)
  {
	 
	 Stopwatch watch;
	 watch.start();
	std::cerr << "Volume pre dilate/erode = " << mLS->ComputeVolume() << std::endl;

    // Determine timestep for stability
    float dt = ComputeTimestep();

    // Propagate level set with stable timestep dt
    // until requested time is reached
    for (float elapsed = 0; elapsed < time;) {

      if (dt > time-elapsed)
        dt = time-elapsed;
      elapsed += dt;

      // Integrate level set function in time using Euler integration
      //IntegrateEuler(dt);
      IntegrateRungeKutta(dt);
    }

	std::cerr << "Volume post dilate/erode = " << mLS->ComputeVolume() << std::endl;

	watch.stop();
	std::cerr << "Time for dilate/erode = " << watch.read() << std::endl;
  }


  virtual float Evaluate(unsigned int i, unsigned int j, unsigned int k)
  {
    // Compute the rate of change (dphi/dt)
	  float x2,y2,z2;
	  LevelSetOperator::Godunov(i,j,k,mF,x2,y2,z2);
	  return -mF*(x2+y2+z2);
  }


};

#endif
