#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"

size_t TOTAL_STEPS = 200000;

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);
  
  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

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
  jointDgain *= 30.0;

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

  for (int sec=3; sec>0; sec--){
    std::cout << "Starting in [" << sec << "]..." << std::endl;
    raisim::USLEEP(1000000);
  }
  size_t TRAJECTORY_STEPS = 1000;
  std::vector<Eigen::VectorXd> trajectory_L, trajectory_R;

  auto jointLimits_R = raipal_R->getJointLimits();
  std::cout << "Joint Limits:" << std::endl;
  for(int i=0; i<jointLimits_R.size(); i++){
    if(i == 4 || i == 5){continue;} // skip unpowered joints
    Eigen::VectorXd negative_full = Eigen::VectorXd::Zero(gcDim_);
    Eigen::VectorXd positive_full = Eigen::VectorXd::Zero(gcDim_);
    negative_full[i] = jointLimits_R[i][0];
    positive_full[i] = jointLimits_R[i][1];

    trajectory_R.push_back(positive_full);
    trajectory_R.push_back(negative_full);
    trajectory_R.push_back(gc_init_);
    
    trajectory_L.push_back(negative_full);
    trajectory_L.push_back(positive_full);
    trajectory_L.push_back(gc_init_);

    std::cout << "  Joint " << i << ": [" << jointLimits_R[i][0] << ", " << jointLimits_R[i][1] << "]" << std::endl;
  }
  trajectory_R[trajectory_R.size()-1][0] =  0.1;
  trajectory_L[trajectory_L.size()-1][0] = -0.1;


  Eigen::VectorXd finalPose_R(gcDim_);
  finalPose_R << 0.37,-0.70,0.69,0.8,0,0,-0.1,0.25,0.79;
  Eigen::VectorXd finalPose_L(gcDim_);
  // finalPose_L << 0.37,-0.70,0.69,0.8,0,0,-0.25,0.25,0.5;
  finalPose_L << 0.37,-0.70,0.69,0.8,0,0,-0.1,0.25,0.79;

  // SIM LOOP
  size_t current_step = 0;
  std::cout << "START!" << std::endl;
  server.startRecordingVideo("raipal_urdf_demo.mp4");
  for (size_t t = 0; t<TOTAL_STEPS; t++){
    RS_TIMED_LOOP(world.getTimeStep()*2e6)
    server.integrateWorldThreadSafe();
    raipal_R->getState(gc, gv);
    
    // analyze step here
    // std::cout<<"STEP " << t << "/" << TOTAL_STEPS << std::endl;

    // set pd targets here
    if(t%TRAJECTORY_STEPS == 0){
      if(current_step < trajectory_R.size()){
        pTarget_R_ = trajectory_R[current_step];
        pTarget_L_ = trajectory_L[current_step];
        std::cout << "Moving to next trajectory point: " << current_step << "/" << trajectory_R.size() << std::endl;
      }
      else if(current_step == trajectory_R.size()){
        pTarget_R_ = finalPose_R;
        pTarget_L_ = finalPose_L;
        std::cout << "Trajectory Finished. Moving to final pose." << std::endl;
      }
      else if(current_step > trajectory_R.size() + 2){
        break;
      }
      raipal_R->setPdTarget(pTarget_R_,dTarget_);
      raipal_L->setPdTarget(pTarget_L_,dTarget_);
      current_step++;
    }

  }

  server.stopRecordingVideo();
  server.killServer();

  std::cout<<"SIMULATION COMPLETE"<<std::endl;

  return 0;
}
