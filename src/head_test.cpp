//e
// Created by Jemin Hwangbo on 2022/04/08.
//

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"


size_t TOTAL_STEPS = 20000;
size_t SPRING_GC_IDX = 1; // this is a 1D version. Usually the spring gc idx is 7
double PRELOAD_N = 1000;
double SPRING_CONST_N_M = 15000;

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);
  
  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  auto head = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/HEAD-URDF/urdf/head.urdf");
  // head -> setComputeInverseDynamics(true);
  // std::cout << "robot was loaded!" << std::endl;

  // world.addGround();
  world.setTimeStep(0.001);

  Eigen::VectorXd gc(head->getGeneralizedCoordinateDim()), gv(head->getDOF());

  gc << 
    0.0, 0.0;

  gv << 
    0.0, 0.0;


  // Declare variables (should be in private section)
  int gcDim_, gvDim_, nJoints_ = 2;
  Eigen::VectorXd gc_init_, gv_init_, gc_, gv_, pTarget_, dTarget_;
  //   int obDim_ = 0, actionDim_ = 0;

  gcDim_ = head->getGeneralizedCoordinateDim();
  gvDim_ = head->getDOF();
  // nJoints_ = gvDim_ - 6; //for floating base
  nJoints_ = 2;

  /// initialize containers
  gc_.setZero(gcDim_);
  gc_init_.setZero(gcDim_);
  gv_.setZero(gvDim_);
  gv_init_.setZero(gvDim_);
  pTarget_.setZero(gcDim_);
  dTarget_.setZero(gvDim_);

  /// set pd gains
  Eigen::VectorXd jointPgain(gvDim_), jointDgain(gvDim_);
  jointPgain.setZero();
  jointPgain.tail(nJoints_).setConstant(50.0);
  jointDgain.setZero();
  jointDgain.tail(nJoints_).setConstant(1);

  head->setPdGains(jointPgain, jointDgain);
  head->setGeneralizedForce(Eigen::VectorXd::Zero(gvDim_));

  // utils::gcRandomize(gc);
  // gc[2] = gc[2] + 3;
  // utils::gvRandomize(gv,15);

  head->setState(gc, gv);
  server.launchServer();
  server.focusOn(head);

  /// if you are using an old version of Raisim, you need this line
  world.integrate1();

  for (int sec=3; sec>0; sec--){
    std::cout << "Starting in [" << sec << "]..." << std::endl;
    raisim::USLEEP(1000000);
  }
  std::cout << "START!" << std::endl;

  // SIM LOOP

  // diagnostic variables
  double z;
  raisim::Vec<3> pos;

  // controller variables
  int latch = 0;

  for (size_t t = 0; t<TOTAL_STEPS; t++){
    RS_TIMED_LOOP(world.getTimeStep()*2e6)
    server.integrateWorldThreadSafe();
    head->getState(gc, gv);
    
    // analyze step here
    std::cout<<"STEP " << t << "/" << TOTAL_STEPS << std::endl;
    
    z = head->getBodyCOM_W()[2][2];
    std::cout << "  CoM Height: " << z*1000 << "mm" << std::endl;
    head->getFramePosition("chin_d455",pos);
    std::cout << "  d455 position  : " << pos.e().transpose() << std::endl;
    head->getFramePosition("chin_mid360",pos);
    std::cout << "  mid360 position  : " << pos.e().transpose() << std::endl;

    // set pd targets here
    // pTarget_ << sin(M_PI*(t%4000/2000.0)), cos(M_PI*(t%4000/2000.0))/2.0; //circling routine
    pTarget_ << sin(M_PI*(t%4000/2000.0)), sin(M_PI*(t%4000/1000.0))/2.0; //figure-8 routine

    head->setPdTarget(pTarget_,dTarget_);

  }

  server.killServer();

  std::cout<<"SIMULATION COMPLETE"<<std::endl;

  return 0;
}
