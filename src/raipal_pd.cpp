#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"
#include <chrono>

#include "frame_timer.hpp"
#include "table_printer.hpp"

double SIM_TIME = 30.0;
// raisim::ControlMode::Type MODE = raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE;
raisim::ControlMode::Type MODE = raisim::ControlMode::FORCE_AND_TORQUE;

/// @brief Crossed 4-bar linkage forward kinematics calculation. Calulates index 4/5 of gc based on index 3.
/// @param gc 9D generalized coordinate vector (with index 4/5 as 0)
void cfbFK(Eigen::VectorXd &gc){
  static double coeffGC5[17] = {7.487348861572753, -88.17268027058002, 470.7676589474978, -1507.949571146181, 3232.924175366478, -4903.587201465159, 5428.248858581333, -4464.681332539945, 2751.745273804874, -1265.589676288956, 418.4687974259759, -86.36058832217348, 6.540591054554072, -1.820824957906633, 2.514881046035373, 1.022698355939166, 0.000002629580961635315};
  static double coeffGC4[17] = {-4.599352352699746, 54.74231689702444, -295.2311950384566, 953.4731854791298, -2052.475077842654, 3099.83081925634, -3363.608943070999, 2636.542385991094, -1477.728449862623, 579.0562233124353, -157.9024440428828, 35.05094543077983, -7.773109367008759, -0.656408504334005, 0.798469056970581, 2.098468226080274, -0.000001063969759856979};
  static double l0 = 112.95, l1 = 60.0, l2 = 45.0, l3 = 85.47477864;
  static double a0 = -0.02249827895, b0 = 0.99988266, c0 = 2.722713633;

  // gc[5] = 0;
  // gc[4] = 0;
  double gc3_pow = 1.0;
  for(int i=16; i>=0; i--){
    gc[5] += coeffGC5[i] * gc3_pow;
    gc[4] += coeffGC4[i] * gc3_pow;
    gc3_pow *= gc[3];
  }

}

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);
  
  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  // auto raipal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal.urdf");
  auto raipal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_saber_R.urdf");
  auto raipalTarget = server.addVisualArticulatedSystem("Target", std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_saber_R.urdf", 1.0,0.0,0.0,0.5);
  // unpowered joint indices: 4/5

  raipal->setName("raipal");
  raipal->setControlMode(MODE);

  // auto raipal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_disabled.urdf");

  // auto raipal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/head/urdf/head.urdf");

  // raipal -> setComputeInverseDynamics(true);
  std::cout << "robot was loaded!" << std::endl;

  // world.addGround();
  world.setTimeStep(0.001);

  // Declare variables (should be in private section)
  int gcDim_, gvDim_, nJoints_ = 7;   // unpowered joint indices: 4/5
  Eigen::VectorXd gc_init_, gv_init_, gc_, gv_, pTarget_, dTarget_;
  //   int obDim_ = 0, actionDim_ = 0;

  gcDim_ = raipal->getGeneralizedCoordinateDim();
  gvDim_ = raipal->getDOF();
  std::cout << "gcDim: " << gcDim_ << ", gvDim: " << gvDim_ << std::endl;
  // nJoints_ = gvDim_ - 6; //for floating base
  // nJoints_ = 2;

  Eigen::VectorXd gc(gcDim_), gv(gvDim_);
  Eigen::VectorXd gf(gvDim_);

  // nominal positions & velocity
  gc.setZero();
  gv.setZero();
  gf.setZero();

  // gc << 
  //   0.0, 0.0;

  // gv << 
  //   0.0, 0.0;

  /// initialize containers
  gc_.setZero(gcDim_);
  gc_init_.setZero(gcDim_);
  gv_.setZero(gvDim_);
  gv_init_.setZero(gvDim_);

  pTarget_.setZero(gcDim_);
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

  if(MODE == raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE){
    raipal->setPdGains(jointPgain, jointDgain);
  } else{
    // when you call setPdGains, you also change the mode to PD_PLUS_FEEDFORWARD_TORQUE!
    // raipal->setPdGains(Eigen::VectorXd::Zero(gvDim_), Eigen::VectorXd::Zero(gvDim_));
  }
  raipal->setGeneralizedForce(Eigen::VectorXd::Zero(gvDim_));

  // utils::gcRandomize(gc);
  // gc[2] = gc[2] + 3;
  // utils::gvRandomize(gv,15);

  raipal->setState(gc, gv);

  server.launchServer();
  server.focusOn(raipal);

  FrameTimer ft(1.0);
  for (int sec=3; sec>0; sec--){
    ft.tick();
    std::cout << "Starting in [" << sec << "]..." << std::endl;
  }
  ft.end();

  /// Generate a simple trajectory
  // size_t TRAJECTORY_STEPS = 750;
  // std::vector<Eigen::VectorXd> waypoints;
  // auto jointLimits_R = raipal->getJointLimits();
  // std::cout << "Joint Limits:" << std::endl;
  // for(int i=0; i<jointLimits_R.size(); i++){
  //   if(i == 4 || i == 5){continue;} // skip unpowered joints
  //   Eigen::VectorXd negative_full = Eigen::VectorXd::Zero(gcDim_);
  //   Eigen::VectorXd positive_full = Eigen::VectorXd::Zero(gcDim_);
  //   negative_full[i] = jointLimits_R[i][0];
  //   positive_full[i] = jointLimits_R[i][1];

  //   waypoints.push_back(positive_full);
  //   waypoints.push_back(negative_full);
  //   // waypoints.push_back(gc_init_);

  //   std::cout << "  Joint " << i << ": [" << jointLimits_R[i][0] << ", " << jointLimits_R[i][1] << "]" << std::endl;
  // }
  // waypoints.push_back(gc_init_);
  // waypoints[waypoints.size()-1][0] =  0.1;

  // Eigen::VectorXd finalPose_R(gcDim_);
  // finalPose_R << 0.37,-0.70,0.69,0.8,0,0,-0.1,0.25,0.79;

  // // SIM LOOP
  // size_t current_step = 0;

  std::cout << "START!" << std::endl;

  pTarget_ << 1, 0, 0, 1, 0, 0, 0, 0, 0.1;

  if(MODE == raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE){
    dTarget_.setZero();
    raipal->setPdTarget(pTarget_,dTarget_);
  }

  cfbFK(pTarget_);
  raipalTarget->setGeneralizedCoordinate(pTarget_);

  ft.enable(world.getTimeStep());

  size_t TOTAL_STEPS = static_cast<size_t>(SIM_TIME / world.getTimeStep());
  size_t t = 0;

  bprinter::TablePrinter tp(&std::cout);
  tp.AddColumn("#",5);
  tp.AddColumn("P Target",10);
  tp.AddColumn("P Actual",10);
  tp.AddColumn("D Target",10);
  tp.AddColumn("D Actual",10);
  tp.AddColumn("Force",10);

  for (t = 0; t<TOTAL_STEPS; t++){
    ft.tick();
    // ft.step(world.getTimeStep());
    
    raipal->getState(gc, gv);
    
    // analyze step here
    // std::cout<<"STEP " << t << "/" << TOTAL_STEPS << std::endl;

    // set pd targets here

    // if(t%TRAJECTORY_STEPS == 0){
    //   if(current_step < waypoints.size()){
    //     pTarget_ = waypoints[current_step];
    //     std::cout << "Moving to next trajectory point: " << current_step + 1 << "/" << waypoints.size() << ": " << pTarget_.transpose() << std::endl;
    //   }
    //   else if(current_step == waypoints.size()){
    //     pTarget_ = finalPose_R;
    //     std::cout << "Moving to final pose           " << std::endl;
    //   }
    //   else if(current_step > waypoints.size() + 2){
    //     // std::cout << "Final timestep: " << t << std::endl;
    //     break;
    //   }
    //   raipal->setPdTarget(pTarget_,dTarget_);
    //   current_step++;
    // }
    if(MODE == raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE){
      raipal->setPdTarget(pTarget_,dTarget_);
    } else {
      // FORCE_AND_TORQUE mode
      Eigen::VectorXd tau(gvDim_);
      tau.setZero();
      // simple pd control to compute feedforward torque
      for(int i=0; i<gcDim_; i++){
        double p_err = pTarget_[i] - gc[i];
        double d_err = dTarget_[i] - gv[i];
        tau[i] = jointPgain[i] * p_err + jointDgain[i] * d_err;
      }
      raipal->setGeneralizedForce(tau);
    }

    server.integrateWorldThreadSafe();
    gf = raipal->getGeneralizedForce().e();

    if((t >= 0 && t < 10) || (t >= 3000 && t < 3010) || (t >= 6000 && t < 6010)){
      std::cout << "MODE: " << (raipal->getControlMode() == raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE ? "PD+FF" : "F & T") << "/ Time: " << world.getWorldTime() << std::endl;

      tp.PrintHeader();
      for (int i=0; i<gcDim_; i++){
        tp << i
            << pTarget_[i]
            << gc[i]
            << dTarget_[i]
            << gv[i]
            << gf[i];
      }
      tp.PrintFooter();
    }
  }
  
  std::cout << "Nominal Time: " << t * world.getTimeStep() << " seconds" << std::endl;
  std::cout << "Elapsed Time: " << ft.end() << " seconds" << std::endl;

  server.killServer();

  std::cout<<"SIMULATION COMPLETE"<<std::endl;
  
  return 0;
}
