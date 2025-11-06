#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"

#include <cmath>

size_t TOTAL_STEPS = 200000;

/// @brief Crossed 4-bar linkage inverse kinematics
/// @param gc 9D generalized coordinate vector (fills in gc[4], gc[5] from gc[3])
void cfbIK(Eigen::VectorXd &gc){

  // calculate gc[5]

  // Quadratic expression: 160000*cosa1^2 - 602400*cosa1 + 160000*sina1^2 + 567009
  double cosa1 = std::cos(gc[3] - 0.0224990394);
  double sina1 = std::sin(gc[3] - 0.0224990394);

  double quad_expr = 
      160000.0 * cosa1 * cosa1 
    - 602400.0 * cosa1 
    + 160000.0 * sina1 * sina1 
    + 567009.0;
    
  // Square root expression
  double sqrt_inner = 3.2620735e11 * cosa1 
                    + 3.4269867e11 * cosa1 * sina1 * sina1 
                    - 7.317723e11 * cosa1 * cosa1 
                    + 3.4269867e11 * cosa1 * cosa1 * cosa1 
                    - 4.5511111e10 * cosa1 * cosa1 * cosa1 * cosa1 
                    - 8.6642058e10 * sina1 * sina1 
                    - 4.5511111e10 * sina1 * sina1 * sina1 * sina1 
                    - 9.1022222e10 * cosa1 * cosa1 * sina1 * sina1 
                    + 1.6657711e11;
  
  double sqrt_term = std::sqrt(sqrt_inner);
  
  // Calculate Y component (first argument of atan2)
  double term1 = 8.2463372e15 * cosa1 * cosa1 
                - 3.104746e16 * cosa1 
                + 8.2463372e15 * sina1 * sina1 
                + 1.7126637e16;
  
  double sinb1_p1 = -(8.0843973e-17 * term1) / sina1;
  
  double term2_inner = 1.5461882e13 * sina1 * sqrt_term 
                      - 3.0229392e19 * cosa1 
                      - 3.2985349e18 * cosa1 * sina1 * sina1 
                      + 1.8628476e19 * cosa1 * cosa1 
                      - 3.2985349e18 * cosa1 * cosa1 * cosa1 
                      + 6.2094919e18 * sina1 * sina1 
                      + 1.2896358e19;
  
  double sinb1_p2 = -(8.0843973e-17 * (400.0 * cosa1 - 753.0) * term2_inner) 
                  / (sina1 * quad_expr);
  
  double sinb1 = sinb1_p1 + sinb1_p2;
  
  // Calculate X component (second argument of atan2)
  double x_numerator = 156472.34 * cosa1 
                      + 0.5 * sina1 * sqrt_term 
                      - 294559.18;
  
  double cosb1 = 0.66666667 * cosa1 
            - (1.0 * x_numerator) / quad_expr 
            - 1.255;

  gc[5] = std::atan2(sinb1, cosb1) - 0.9998834205;
  if(gc[5] < -M_PI){
      gc[5] += 2.0 * M_PI;
  }

  // calculate gc[4]
  double l1 = 60.0, l2 = 45.0, l3 = 85.4748;
  double a = M_PI/2 - (gc[3] - 0.0224990394);
  double b = std::acos((l1*sina1 + l2*sinb1)/l3);

  gc[4] = 2.72271 - (a + b);

}

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);
  
  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  // auto raipal_R = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_R.urdf");
  // auto raipal_L = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_L.urdf");
  auto raipal_R = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_R.urdf");
  auto raipal_L = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_L.urdf");

  // unpowered joint indices: 4/5

  raipal_R->setName("raipal_R");
  raipal_L->setName("raipal_L");

  // auto raipal_R = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal_R/urdf/raipal_R_disabled.urdf");

  // auto raipal_R = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/head/urdf/head.urdf");

  // raipal_R -> setComputeInverseDynamics(true);
  std::cout << "robot was loaded!" << std::endl;

  // world.addGround();
  world.setTimeStep(0.001);

  // Declare variables (should be in private section)
  int gcDim_, gvDim_, nJoints_ = 7;   // unpowered joint indices: 4/5
  Eigen::VectorXd gc_init_, gv_init_, gc_, gv_, pTarget_R_, pTarget_L_, dTarget_;
  //   int obDim_ = 0, actionDim_ = 0;

  gcDim_ = raipal_R->getGeneralizedCoordinateDim();
  gvDim_ = raipal_R->getDOF();
  std::cout << "gcDim: " << gcDim_ << ", gvDim: " << gvDim_ << std::endl;
  // nJoints_ = gvDim_ - 6; //for floating base
  // nJoints_ = 2;

  Eigen::VectorXd gc(gcDim_), gv(gvDim_);

  // nominal positions & velocity
  gc.setZero();
  gv.setZero();

  // gc << 
  //   0.0, 0.0;

  // gv << 
  //   0.0, 0.0;

  /// initialize containers
  gc_.setZero(gcDim_);
  gc_init_.setZero(gcDim_);
  gv_.setZero(gvDim_);
  gv_init_.setZero(gvDim_);

  pTarget_R_.setZero(gcDim_);
  pTarget_L_.setZero(gcDim_);
  dTarget_.setZero(gvDim_);

  /// set pd gains
  Eigen::VectorXd jointPgain(gvDim_), jointDgain(gvDim_);
  jointPgain.setOnes();
  jointPgain *= 200;
  // jointPgain[1] = 1000.0; // adduction needs stronger gravity compensation
  jointDgain.setOnes();
  jointDgain *= 20.0;

  // unpowered joint indices: 4/5
  jointPgain[4] = 0.0;
  jointPgain[5] = 0.0;
  jointDgain[4] = 0.0;
  jointDgain[5] = 0.0;

  // jointDgain.tail(nJoints_).setConstant(1);

  raipal_R->setPdGains(jointPgain, jointDgain);
  raipal_R->setGeneralizedForce(Eigen::VectorXd::Zero(gvDim_));

  raipal_L->setPdGains(jointPgain, jointDgain);
  raipal_L->setGeneralizedForce(Eigen::VectorXd::Zero(gvDim_));

  // utils::gcRandomize(gc);
  // gc[2] = gc[2] + 3;
  // utils::gvRandomize(gv,15);

  raipal_R->setState(gc, gv);
  raipal_L->setState(gc, gv);

  server.launchServer();
  server.focusOn(raipal_R);

  /// if you are using an old version of Raisim, you need this line
  world.integrate1();

  auto jointLimits_R = raipal_R->getJointLimits();
  auto jointLimits_L = raipal_L->getJointLimits();

  std::cout << "Joint Limits:" << std::endl;
  for(int i=0; i<jointLimits_R.size(); i++){
    std::cout << "Joint " << i << ": R[" << jointLimits_R[i][0] << ", " << jointLimits_R[i][1] << "]" << " L[" << jointLimits_L[i][0] << ", " << jointLimits_L[i][1] << "]" << std::endl;
  }

  // COUNTDOUWN
  for (int sec=3; sec>0; sec--){
    std::cout << "Starting in [" << sec << "]..." << std::endl;
    raisim::USLEEP(1000000);
  }

  // SIM LOOP
  size_t current_step = 0;
  std::cout << "START!" << std::endl;
  server.startRecordingVideo("raipal_urdf_demo.mp4");

  double incr = 0.001;

  std::cout << "Nominal gc[3]: " << gc_[3] << ", gc[5]: " << gc_[5] << std::endl;

  for (size_t t = 0; t<2000; t++){
    RS_TIMED_LOOP(world.getTimeStep()*2e6)
    server.integrateWorldThreadSafe();
    raipal_R->getState(gc, gv);
    
    // analyze step here
    // std::cout<<"STEP " << t << "/" << TOTAL_STEPS << std::endl;

    gc_[3] += incr;
    if (gc_[3] >= jointLimits_R[3][1] || gc_[3] <= jointLimits_R[3][0]) {
      incr = -incr;
      std::cout << "Reversing direction at step " << t << ", gc[3]: " << gc_[3] << ", gc[5]: " << gc_[5] << std::endl;
    }

    cfbIK(gc_);
    // std::cout << "step " << t << "/" << TOTAL_STEPS << " :" << gc_[3] << " " << gc_[5] << std::endl;
    raipal_R->setState(gc_,gv_);
    raipal_L->setState(gc_,gv_);
  }

  gc_.setZero();
  gv_.setZero();
  incr = 0.003;
  raipal_R->setState(gc_,gv_);
  raipal_L->setState(gc_,gv_);

  Eigen::VectorXd gc_IK(gcDim_);

  for (size_t t = 0; t<10000; t++){
    RS_TIMED_LOOP(world.getTimeStep()*2e6)
    server.integrateWorldThreadSafe();

    gc_[3] += incr;
    if (gc_[3] >= jointLimits_R[3][1] + 0.1 || gc_[3] <= jointLimits_R[3][0] - 0.1) {
      incr = -incr;
    }

    // std::cout << "step " << t << "/" << TOTAL_STEPS << " :" << gc_[3] << " " << gc_[5] << std::endl;
    raipal_R->setPdTarget(gc_,gv_);
    raipal_L->setPdTarget(gc_,gv_);

    if(t%10 == 0){
      raipal_R->getState(gc, gv);
      gc_IK = gc;
      cfbIK(gc_IK);
      double err_4 = (gc_IK[4] - gc[4]) * 180 / M_PI;
      double err_5 = (gc_IK[5] - gc[5]) * 180 / M_PI;
      if(err_4 > 0.01 || err_5 > 0.01){
        std::cout << "Warning: Large IK error at step " << t << std::endl;
        std::cout << "[" << t << "] gc[3]: " << gc[3] << ": gc[4]: " << gc_IK[4] << " err: " << err_4 << "deg, gc[5]: " << gc_IK[5] << " err: " << err_5 << "deg" << std::endl;
      }
    }

  }

  server.stopRecordingVideo();
  server.killServer();

  std::cout<<"SIMULATION COMPLETE"<<std::endl;

  return 0;
}

