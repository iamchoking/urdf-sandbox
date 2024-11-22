//
// Created by Jemin Hwangbo on 2022/04/08.
//

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"

// #include "exercise_3_20190673.hpp"
// START HEADER FILE
/// do not change the name of the method

inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!


  return Eigen::MatrixXd::Ones(18,18);
}

// Header File End

bool analyzeStep(const Eigen::VectorXd& gc, size_t t, raisim::RaisimServer* server, raisim::ArticulatedSystem* anymal){
  std::cout << "STEP[" << t << "]";
  /// TEMPLATE (do some testing here)

  std::cout << std::endl;
  /// TEMPLATE (add return condition here)
  return true;
}

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  // kinova
  auto anymal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/anymal_c/urdf/anymal.urdf");

  world.addGround();
  world.setTimeStep(0.001);

  // kinova configuration
  Eigen::VectorXd gc(anymal->getGeneralizedCoordinateDim()), gv(anymal->getDOF());
  gc << 0, 0, 0.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8; /// Jemin: I'll randomize the gc, gv when grading
  gv << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8;

  utils::gcRandomize(gc);
  utils::gvRandomize(gv);
  anymal->setState(gc, gv);
  server.launchServer();

  /// if you are using an old version of Raisim, you need this line
  world.integrate1();
  bool correct = analyzeStep(gc,0,&server,anymal);

  for (int sec=5; sec>0; sec--){
    std::cout << "Dropping in [" << sec << "]..." << std::endl;
    raisim::USLEEP(1000000);
  }
  std::cout << "DROP!" << std::endl;

  // std::cout<<"mass matrix should be \n"<< anymal->getMassMatrix().e()<<std::endl;
  for (size_t i = 0; i<2000; i++){
    RS_TIMED_LOOP(world.getTimeStep()*2e6);
    if(i%10 == 0){
      correct = analyzeStep(gc,i,&server,anymal) && correct;
    }
    server.integrateWorldThreadSafe();
    anymal->getState(gc, gv);
  }

  server.killServer();

  if(correct) {
    std::cout<<"TEST PASSED"<<std::endl;
  } else {
    std::cout<<"TEST FAILED"<<std::endl;
  }

  return 0;
}
