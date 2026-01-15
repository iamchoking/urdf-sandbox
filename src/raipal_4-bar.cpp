#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"

#include <cmath>
#include <iomanip>

#include <raipal_kinematics/raipal_kinematics.hpp>

namespace rk = raipal::kinematics;

size_t TOTAL_STEPS = 200000;

/// @brief Crossed 4-bar linkage inverse kinematics (closed-form full solution, ~1e-9 rad error)
/// @param gc 9D generalized coordinate vector (fills in gc[4], gc[5] from gc[3])
void cfbFKAnalytic(Eigen::VectorXd &gc){

  // calculate gc[5] from gc[3]

  // Quadratic expression: 160000*cosa1^2 - 602400*cosa1 + 160000*sina1^2 + 567009
  double cosa1 = std::cos(gc[3] - 0.02249827895);
  double sina1 = std::sin(gc[3] - 0.02249827895);

  double quad_expr = 
      160000.0 * cosa1 * cosa1 
    - 602400.0 * cosa1 
    + 160000.0 * sina1 * sina1 
    + 567009.0;
    
  // Square root expression
  double sqrt_inner = 
      3.2620769e11 * cosa1 
    + 3.4269867e11 * cosa1 * sina1 * sina1 
    - 7.3177239e11 * cosa1 * cosa1 
    + 3.4269867e11 * cosa1 * cosa1 * cosa1 
    - 4.5511111e10 * cosa1 * cosa1 * cosa1 * cosa1 
    - 8.664215e10  * sina1 * sina1 
    - 4.5511111e10 * sina1 * sina1 * sina1 * sina1 
    - 9.1022222e10 * cosa1 * cosa1 * sina1 * sina1 
    + 1.6657692e11;
  
  double sqrt_term = std::sqrt(sqrt_inner);
  
  // Calculate Y component (first argument of atan2)
  double term1 = 
      4.9478023e16 * cosa1 * cosa1 
    - 1.8628476e17 * cosa1 
    + 4.9478023e16 * sina1 * sina1 
    + 1.0275987e17;
  
  double sinb1_p1 = -(1.3473996e-17 * term1) / sina1;
  
  double term2_inner = 
      9.2771294e13 * sina1 * sqrt_term 
    - 1.8137637e20 * cosa1 
    - 1.9791209e19 * cosa1 * sina1 * sina1 
    + 1.1177085e20 * cosa1 * cosa1 
    - 1.9791209e19 * cosa1 * cosa1 * cosa1 
    + 3.7256952e19 * sina1 * sina1 
    + 7.7378183e19;
  
  double sinb1_p2 = 
    -(1.3473996e-17 * (400.0 * cosa1 - 753.0) * term2_inner) / (sina1 * quad_expr);
  
  double sinb1 = sinb1_p1 + sinb1_p2;
  
  // Calculate X component (second argument of atan2)
  double x_numerator = 
      156472.23 * cosa1 
    + 0.5 * sina1 * sqrt_term 
    - 294558.97;

  double cosb1 = 
      0.66666667 * cosa1 
    - x_numerator / quad_expr 
    - 1.255;

  gc[5] = std::atan2(sinb1, cosb1) - 0.99988266;
  if(gc[5] < -M_PI){
      gc[5] += 2.0 * M_PI;
  }

  // calculate gc[4]
  double l1 = 60.0, l2 = 45.0, l3 = 85.47477864;

  double a = M_PI/2 - (gc[3] - 0.02249827895);
  double b = std::acos((l1*sina1 + l2*sinb1)/l3);

  gc[4] = 2.722713633 - (a + b);

}

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);
  
  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  // auto raipal_R = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_R.urdf");
  // auto raipal_L = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_L.urdf");
  auto raipal_R = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_stub-10_R.urdf");
  auto raipal_L = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_stub-10_L.urdf");

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

  // get joint limits
  auto jointLimits_R = raipal_R->getJointLimits();
  auto jointLimits_L = raipal_L->getJointLimits();

  std::cout << "Joint Limits:" << std::endl;
  for(int i=0; i<jointLimits_R.size(); i++){
    std::cout << "Joint " << i << ": R[" << jointLimits_R[i][0] << ", " << jointLimits_R[i][1] << "]" << " L[" << jointLimits_L[i][0] << ", " << jointLimits_L[i][1] << "]" << std::endl;
  }

  // ground-truth test
  double gtInputDeg[] = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,89.1304};

  for (double input : gtInputDeg){
    gc[3] = input * M_PI / 180.0;
    cfbFKAnalytic(gc);
    std::cout << std::setprecision(10) << std::endl;
    std::cout << "Input: " << input << " deg, Output: " << gc[5] * 180.0 / M_PI << " deg" << std::endl;
  }

  // COUNTDOUWN
  std::cout << "====== SWEEP TEST ======" << std::endl;
  for (int sec=3; sec>0; sec--){
    std::cout << "Starting in [" << sec << "]..." << std::endl;
    raisim::USLEEP(1000000);
  }
  std::cout << "START!" << std::endl;

  // SIM LOOP
  size_t current_step = 0;
  server.startRecordingVideo("raipal_urdf_demo.mp4");

  double incr = 0.001;

  std::cout << "Nominal gc[3]: " << gc_[3] << ", gc[5]: " << gc_[5] << std::endl;
  Eigen::VectorXd gc__(gcDim_);

  for (size_t t = 0; t<4000; t++){
    RS_TIMED_LOOP(world.getTimeStep()*2e6)
    server.integrateWorldThreadSafe();
    
    // analyze step here
    // std::cout<<"STEP " << t << "/" << TOTAL_STEPS << std::endl;

    gc_[3] += incr;
    if (gc_[3] >= jointLimits_R[3][1] || gc_[3] <= jointLimits_R[3][0]) {
      incr = -incr;
      std::cout << "Reversing direction at step " << t << ", gc[3]: " << gc_[3] * 180.0 / M_PI << "deg, gc[5]: " << gc_[5] * 180.0 / M_PI << "deg" << std::endl;
    }

    gc__ = gc_;

    auto curTime = std::chrono::high_resolution_clock::now();
    cfbFKAnalytic(gc_);
    if(t%100 == 0){
      std::cout << "Closed-form IK took : " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - curTime).count() << " ns" << std::endl;
      curTime = std::chrono::high_resolution_clock::now();

      rk::cfbFK(gc__);

      std::cout << "Polynomial IK took : " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - curTime).count() << " ns" << std::endl;
      curTime = std::chrono::high_resolution_clock::now();

      std::cout << std::setprecision(10);
      std::cout << "Polynomial IK approx. Error: gc[5]: " << (gc__[5] - gc_[5]) * 180.0 / M_PI << " deg / gc[4]: " << (gc__[4] - gc_[4]) * 180.0 / M_PI << " deg"<< std::endl << std::endl;
    }

    // std::cout << "step " << t << "/" << TOTAL_STEPS << " :" << gc_[3] << " " << gc_[5] << std::endl;
    raipal_R->setState(gc_,gv_);
    raipal_L->setState(gc_,gv_);
  }

  // COUNTDOUWN
  std::cout << "====== Loop Closure Comparison ======" << std::endl;
  for (int sec=3; sec>0; sec--){
    std::cout << "Starting in [" << sec << "]..." << std::endl;
    raisim::USLEEP(1000000);
  }
  std::cout << "START!" << std::endl;

  gc_.setZero();
  gv_.setZero();
  incr = 0.003;
  raipal_R->setState(gc_,gv_);
  raipal_L->setState(gc_,gv_);

  Eigen::VectorXd gc_IK(gcDim_);
  Eigen::VectorXd gv_IK(gvDim_);
  double maxErrorPos = 0.0;
  double avgErrorPos = 0.0;
  double maxErrorGV4 = 0.0;
  double avgErrorGV4 = 0.0;
  double maxErrorVel = 0.0;
  double avgErrorVel = 0.0;

  // Variables to track max errors within 100-step windows
  double windowMaxErrGc4 = 0.0;
  double windowMaxErrGc5 = 0.0;
  double windowMaxErrGv4 = 0.0;
  double windowMaxErrGv5 = 0.0;

  size_t simSteps = 3000;

  for (size_t t = 0; t<simSteps; t++){
    RS_TIMED_LOOP(world.getTimeStep()*2e6)
    server.integrateWorldThreadSafe();

    gc_[3] += incr;

    gc_[1] -= incr*2;
    gc_[0] += incr*3;
    gc_[7] += incr*3;
    gc_[8] += incr*3;

    if (gc_[3] >= jointLimits_R[3][1] + 0.3 || gc_[3] <= jointLimits_R[3][0] - 0.3) {
      incr = -incr;
    }

    // std::cout << "step " << t << "/" << TOTAL_STEPS << " :" << gc_[3] << " " << gc_[5] << std::endl;
    raipal_R->setPdTarget(gc_,gv_);
    raipal_L->setPdTarget(gc_,gv_);

    raipal_R->getState(gc, gv);
    gc_IK = gc;
    gv_IK = gv;
    rk::cfbFK(gc_IK, gv_IK);

    // all errors calculated in degrees
    double err_gc4 = (gc_IK[4] - gc[4]) * 180 / M_PI;
    double err_gc5 = (gc_IK[5] - gc[5]) * 180 / M_PI;
    double err_gv4 = (gv_IK[4] - gv[4]) * 180 / M_PI;
    double err_gv5 = (gv_IK[5] - gv[5]) * 180 / M_PI;
    
    maxErrorPos = std::max(maxErrorPos, std::abs(err_gc5));
    avgErrorPos += std::abs(err_gc5) / simSteps;

    maxErrorVel = std::max(maxErrorVel, std::abs(err_gv5));
    avgErrorVel += std::abs(err_gv5) / simSteps;

    // Track maximum errors within the current 100-step window
    windowMaxErrGc4 = std::max(windowMaxErrGc4, std::abs(err_gc4));
    windowMaxErrGc5 = std::max(windowMaxErrGc5, std::abs(err_gc5));
    windowMaxErrGv4 = std::max(windowMaxErrGv4, std::abs(err_gv4));
    windowMaxErrGv5 = std::max(windowMaxErrGv5, std::abs(err_gv5));

    if(t%100 == 0){
      // if(windowMaxErrGc4 > 0.01 || windowMaxErrGc5 > 0.01){
      //   std::cout << "Warning: Large position (gc) error at step " << t << std::endl;
      //   std::cout << "  Max gc[4] err in window: " << windowMaxErrGc4 << " deg, Max gc[5] err in window: " << windowMaxErrGc5 << " deg" << std::endl;
      // }

      // if(windowMaxErrGv4 > 10 || windowMaxErrGv5 > 10){
      if(true){
        std::cout << "Warning: Large velocity (gv) error at step " << t << std::endl;
        std::cout << "  Max gv[4] err in window: " << windowMaxErrGv4 << " deg/s, Max gv[5] err in window: " << windowMaxErrGv5 << " deg/s" << std::endl;
      }

      // Reset window max values for next 100-step window
      windowMaxErrGc4 = 0.0;
      windowMaxErrGc5 = 0.0;
      windowMaxErrGv4 = 0.0;
      windowMaxErrGv5 = 0.0;
    }
  }

  std::cout << "Maximum gc[5] error: " << maxErrorPos << " deg" << std::endl;
  std::cout << "Average gc[5] error: " << avgErrorPos << " deg" << std::endl;

  std::cout << "Maximum gv[5] error: " << maxErrorVel << " deg/s" << std::endl;
  std::cout << "Average gv[5] error: " << avgErrorVel << " deg/s" << std::endl;

  server.stopRecordingVideo();
  server.killServer();

  std::cout<<"SIMULATION COMPLETE"<<std::endl;

  return 0;
}

