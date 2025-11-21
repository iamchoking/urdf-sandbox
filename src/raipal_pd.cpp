#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"
#include <chrono>

#include "frame_timer.hpp"
#include "table_printer.hpp"

// double SIM_DT = 0.0025;
double SIM_DT = 0.001;
double MAX_TIME = 15.0;
// double LPF_CUTOFF_FREQ = 100.0; // Hz (negative value for no cutoff)
double LPF_CUTOFF_FREQ = 20.0; // Hz (negative value for no cutoff)
// double LPF_CUTOFF_FREQ = 100.0; // Hz (negative value for no cutoff)
bool RANDOM_DT = false;

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
  auto raipal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_stub-0_R.urdf");
  auto raipalTarget = server.addVisualArticulatedSystem("Target", std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_stub-0_R.urdf", 1.0,0.0,0.0,0.5);
  // unpowered joint indices: 4/5

  raipal->setName("raipal");
  raipal->setControlMode(MODE);

  // auto raipal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_disabled.urdf");

  // auto raipal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/head/urdf/head.urdf");

  // raipal -> setComputeInverseDynamics(true);
  std::cout << "robot was loaded!" << std::endl;

  // world.addGround();
  world.setTimeStep(SIM_DT);

  // Declare variables (should be in private section)
  int gcDim, gvDim, nJoints = 7;   // unpowered joint indices: 4/5
  gcDim = raipal->getGeneralizedCoordinateDim();
  gvDim = raipal->getDOF();
  std::cout << "gcDim: " << gcDim << ", gvDim: " << gvDim << std::endl;

  Eigen::VectorXd gc(gcDim), gv(gvDim);
  Eigen::VectorXd gf(gvDim);
  gc.setZero();gv.setZero();gf.setZero();

  // low-pass filter implementation
  double lpf_alpha;
  if(LPF_CUTOFF_FREQ < 0.0){lpf_alpha = 1.0;}
  else{
    lpf_alpha = 2.0 * M_PI * LPF_CUTOFF_FREQ * SIM_DT / (2.0 * M_PI * LPF_CUTOFF_FREQ * SIM_DT + 1.0);
    std::cout << "Low-pass filter applied with cutoff frequency: " << LPF_CUTOFF_FREQ << " Hz, alpha: " << lpf_alpha << std::endl;
  }
  Eigen::VectorXd gcRaw(gcDim), gvRaw(gvDim);

  /// initial (nominal) pose
  Eigen::VectorXd gc_init, gv_init, pTarget_, dTarget_;
  gc_init.setZero(gcDim) ; gv_init.setZero(gvDim);
  pTarget_.setZero(gcDim); dTarget_.setZero(gvDim);

  /// set pd gains
  Eigen::VectorXd jointPgain(gvDim), jointDgain(gvDim);
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

  // jointDgain.tail(nJoints).setConstant(1);

  if(MODE == raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE){
    raipal->setPdGains(jointPgain, jointDgain);
  } else{
    // when you call setPdGains, you also change the mode to PD_PLUS_FEEDFORWARD_TORQUE!
    // raipal->setPdGains(Eigen::VectorXd::Zero(gvDim), Eigen::VectorXd::Zero(gvDim));
  }
  raipal->setGeneralizedForce(Eigen::VectorXd::Zero(gvDim));

  // utils::gcRandomize(gc);
  // gc[2] = gc[2] + 3;
  // utils::gvRandomize(gv,15);

  raipal->setState(gc_init, gv_init);

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
  //   Eigen::VectorXd negative_full = Eigen::VectorXd::Zero(gcDim);
  //   Eigen::VectorXd positive_full = Eigen::VectorXd::Zero(gcDim);
  //   negative_full[i] = jointLimits_R[i][0];
  //   positive_full[i] = jointLimits_R[i][1];

  //   waypoints.push_back(positive_full);
  //   waypoints.push_back(negative_full);
  //   // waypoints.push_back(gc_init);

  //   std::cout << "  Joint " << i << ": [" << jointLimits_R[i][0] << ", " << jointLimits_R[i][1] << "]" << std::endl;
  // }
  // waypoints.push_back(gc_init);
  // waypoints[waypoints.size()-1][0] =  0.1;

  // Eigen::VectorXd finalPose_R(gcDim);
  // finalPose_R << 0.37,-0.70,0.69,0.8,0,0,-0.1,0.25,0.79;

  // // SIM LOOP
  // size_t current_step = 0;

  std::cout << "START!" << std::endl;

  // pTarget_ << 1, 0, 0, 1, 0, 0, 0, 0, 0.1;

  if(MODE == raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE){
    dTarget_.setZero();
    raipal->setPdTarget(pTarget_,dTarget_);
  }

  cfbFK(pTarget_);
  raipalTarget->setGeneralizedCoordinate(pTarget_);

  ft.enable(world.getTimeStep());

  size_t t = 0;

  bprinter::TablePrinter tp(&std::cout);
  tp.AddColumn("#",5);
  tp.AddColumn("P Target",10);
  tp.AddColumn("P Actual",10);
  tp.AddColumn("D Target",10);
  tp.AddColumn("D Actual",10);
  tp.AddColumn("Force",10);

  double random_dt = SIM_DT;
  while (ft.getElapsedTime() < MAX_TIME) {
    ft.tick(random_dt);
    // randomize simulation timestep (low pass filter)!
    raipal->getState(gcRaw, gvRaw);
    
    if(MODE == raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE){
      // no low-pass filtering
      gc = gcRaw;
      gv = gvRaw;
    } else {
      // FORCE_AND_TORQUE mode
      // apply low-pass filter
      gc = lpf_alpha * gcRaw + (1.0 - lpf_alpha) * gc;
      gv = lpf_alpha * gvRaw + (1.0 - lpf_alpha) * gv;
    }
    
    if(MODE == raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE){
      raipal->setPdTarget(pTarget_,dTarget_);
    } else {
      // FORCE_AND_TORQUE mode
      Eigen::VectorXd tau(gvDim);
      tau.setZero();
      // simple pd control to compute feedforward torque
      for(int i=0; i<gcDim; i++){
        double p_err = pTarget_[i] - gc[i];
        double d_err = dTarget_[i] - gv[i];
        tau[i] = jointPgain[i] * p_err + jointDgain[i] * d_err;
      }
      raipal->setGeneralizedForce(tau);
    }

    server.integrateWorldThreadSafe();
    gf = raipal->getGeneralizedForce().e();

    // if((t >= 0 && t < 10) || (t >= 3000 && t < 3010) || (t >= 6000 && t < 6010)){
    if(true){
      std::cout << "MODE: " << (raipal->getControlMode() == raisim::ControlMode::PD_PLUS_FEEDFORWARD_TORQUE ? "PD+FF" : "F & T") << " / Time: " << world.getWorldTime() << " / dt: " << world.getTimeStep() << " / alpha: " << lpf_alpha << std::endl;

      tp.PrintHeader();
      for (int i=0; i<gcDim; i++){
        tp << i
            << pTarget_[i]
            << gcRaw[i]
            << dTarget_[i]
            << gvRaw[i]
            << gf[i];
      }
      tp.PrintFooter();
    }

    // randomize dt (0.001 ~ 0.0025)
    if(RANDOM_DT){
      random_dt = 0.001 + (0.0015 * ((double) rand() / RAND_MAX));
      lpf_alpha = 2.0 * M_PI * LPF_CUTOFF_FREQ * random_dt / (2.0 * M_PI * LPF_CUTOFF_FREQ * random_dt + 1.0);

      world.setTimeStep(random_dt);
      // ft.step(world.getTimeStep());
    }
  }
  
  std::cout << "Nominal Time: " << world.getWorldTime() << " seconds" << std::endl;
  std::cout << "Elapsed Time: " << ft.end() << " seconds" << std::endl;

  server.killServer();

  std::cout<<"SIMULATION COMPLETE"<<std::endl;
  
  return 0;
}
