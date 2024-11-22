//
// Created by Jemin Hwangbo on 2022/03/17.
//


#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"

#include "exercise_2_20190673.hpp"
// excercise1_20190673.hpp

// hpp end

bool analyzeStep(Eigen::VectorXd gc, Eigen::VectorXd gv, size_t t,const raisim::Vec<3>& vTrueVec,const raisim::Vec<3>& wTrueVec,const raisim::Vec<3>& pos, raisim::RaisimServer* server){
  std::cout << "STEP[" << t << "]";


  bool vCorrect;
  bool wCorrect;

  // auto footLink = setupLinks();
  // footLink->calculateKinematics(gc);

  Eigen::Vector3d vTrue = vTrueVec.e();
  Eigen::Vector3d wTrue = wTrueVec.e();
  Eigen::Vector3d vCalc = getFootLinearVelocity (gc,gv);
  Eigen::Vector3d wCalc = getFootAngularVelocity(gc,gv);

  // graphics

  auto debugO = server->getVisualObject("debug_O");
  auto debugVTrue = server->getVisualObject("debug_vTrue");
  auto debugVCalc = server->getVisualObject("debug_vCalc");

  debugO->setPosition(pos.e());
  debugVTrue->setPosition(pos.e()+vTrue*0.2);
  debugVCalc->setPosition(pos.e()+vCalc*0.2);

  double vDiffNorm = (vCalc-vTrue).norm();
  double wDiffNorm = (wCalc-wTrue).norm();

  if(vDiffNorm < 1e-8){std::cout << "[O";}
  else{std::cout << "[X";}
  if(wDiffNorm < 1e-8){std::cout << "O]";}
  else{std::cout << "X]";}

  std::cout << " v: " << vCalc.transpose() << "(<> " << vTrue.transpose() << ", diff: " << vDiffNorm << ")";
  std::cout << " w: " << wCalc.transpose() << "(<> " << wTrue.transpose() << ", diff: " << wDiffNorm << ")";

  std::cout << std::endl;

  if (vDiffNorm > 1e-8 || wDiffNorm > 1e-8){return false;}
  else                                     {return true;}
}

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();
  world.setTimeStep(0.001);

  auto debugO     = server.addVisualSphere("debug_O", 0.1);
  auto debugVTrue = server.addVisualSphere("debug_vTrue", 0.03);
  auto debugVCalc = server.addVisualSphere("debug_vCalc", 0.015);
  debugO->setColor(1,1,1,1);
  debugVTrue->setColor(1,0,0,1);
  debugVCalc->setColor(0,1,0,1);

  // a1
  // anymal
  auto anymal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/anymal_c/urdf/anymal.urdf");
  anymal->setName("anymal");
  server.focusOn(anymal);

  // a1 configuration
  Eigen::VectorXd gc(anymal->getGeneralizedCoordinateDim());
  Eigen::VectorXd gv(anymal->getDOF());

  gc << 0, 0, 10.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8;
  gv << 0.1, 0.2, 0.3, 0.1, 0.4, 0.3, 0.1,0.1,0.1, 0.2,0.2,0.2, 0.3,0.3,0.3, 0.4,0.4,0.4;

  utils::gcRandomize(gc);
  utils::gvRandomize(gv,5);

  anymal->setState(gc, gv);

  // visualization
  server.launchServer();
  raisim::Vec<3> footVel, footAngVel, pos, posRef;

  // auto target_name = "base_top_shell";// use "base_top_shell" for base
  // auto target_name = "LH_HAA";
  auto target_name = "LH_shank_fixed_LH_FOOT";

  auto ref_name = "base_top_shell";

  std::cout << "Setup Complete!" << std::endl;
  // TODO: do step analyis

  anymal->getFramePosition(ref_name,posRef);
  anymal->getFramePosition(target_name, pos);
  anymal->getFrameVelocity(target_name, footVel);
  anymal->getFrameAngularVelocity(target_name, footAngVel);
  bool answerCorrect = analyzeStep(gc,gv,42,footVel,footAngVel,pos,&server);

  for (int sec=5; sec>0; sec--){
    std::cout << "Dropping in [" << sec << "]..." << std::endl;
    // TODO; visualize
    raisim::USLEEP(1000000);
  }
  std::cout << "DROP!" << std::endl;

  for (size_t i=0; i<2000; i++) {
    RS_TIMED_LOOP(world.getTimeStep()*2e6);

    anymal->getFramePosition(ref_name,posRef);
    anymal->getFramePosition(target_name, pos);
    anymal->getFrameVelocity(target_name, footVel);
    anymal->getFrameAngularVelocity(target_name, footAngVel);

    if (i % 10 == 0){
      answerCorrect = (analyzeStep(gc,gv,i,footVel,footAngVel,pos,&server) && answerCorrect);
      // std::cout << "true pos: " <<  pos.e().transpose() << std::endl;
      // std::cout << "true wriee: "<< (pos.e()-posRef.e()).transpose() << std::endl;
    }
    server.integrateWorldThreadSafe();
    anymal->getState(gc, gv);
  }

  server.killServer();

  if(answerCorrect) {
    std::cout<<"TEST PASSED"<<std::endl;
  } else {
    std::cout<<"TEST FAILED"<<std::endl;
  }

  return 0;
}
