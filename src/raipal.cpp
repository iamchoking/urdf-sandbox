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

  auto raipal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_R.urdf");
  // auto raipal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_disabled.urdf");

  // auto raipal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/head/urdf/head.urdf");

  // raipal -> setComputeInverseDynamics(true);
  std::cout << "robot was loaded!" << std::endl;

  // world.addGround();
  world.setTimeStep(0.001);

  // Declare variables (should be in private section)
  int gcDim_, gvDim_, nJoints_ = 2;
  Eigen::VectorXd gc_init_, gv_init_, gc_, gv_, pTarget_, dTarget_;
  //   int obDim_ = 0, actionDim_ = 0;

  gcDim_ = raipal->getGeneralizedCoordinateDim();
  gvDim_ = raipal->getDOF();
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

  pTarget_.setZero(gcDim_);
  dTarget_.setZero(gvDim_);

  /// set pd gains
  Eigen::VectorXd jointPgain(gvDim_), jointDgain(gvDim_);
  jointPgain.setZero();
  // jointPgain.tail(nJoints_).setConstant(50.0);
  jointDgain.setZero();
  // jointDgain.tail(nJoints_).setConstant(1);

  raipal->setPdGains(jointPgain, jointDgain);
  raipal->setGeneralizedForce(Eigen::VectorXd::Zero(gvDim_));

  // utils::gcRandomize(gc);
  // gc[2] = gc[2] + 3;
  // utils::gvRandomize(gv,15);

  raipal->setState(gc, gv);
  server.launchServer();
  server.focusOn(raipal);

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
    raipal->getState(gc, gv);
    
    // analyze step here
    std::cout<<"STEP " << t << "/" << TOTAL_STEPS << std::endl;

    
    // z = raipal->getBodyCOM_W()[2][2];
    // std::cout << "  CoM Height: " << z*1000 << "mm" << std::endl;
    // raipal->getFramePosition("chin_d455",pos);
    // std::cout << "  d455 position  : " << pos.e().transpose() << std::endl;
    // raipal->getFramePosition("chin_mid360",pos);
    // std::cout << "  mid360 position  : " << pos.e().transpose() << std::endl;

    // set pd targets here
    // // pTarget_ << sin(M_PI*(t%4000/2000.0)), cos(M_PI*(t%4000/2000.0))/2.0; //circling routine
    // pTarget_ << sin(M_PI*(t%4000/2000.0)), sin(M_PI*(t%4000/1000.0))/2.0; //figure-8 routine

    // raipal->setPdTarget(pTarget_,dTarget_);

  }

  server.killServer();

  std::cout<<"SIMULATION COMPLETE"<<std::endl;

  return 0;
}
