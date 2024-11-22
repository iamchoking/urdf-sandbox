//e
// Created by Jemin Hwangbo on 2022/04/08.
//

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"

#include "exercise_3_20190673.hpp"
// START HEADER FILE

// HEADER END

bool analyzeStep(const Eigen::VectorXd& gc, size_t t, raisim::RaisimServer* server, raisim::ArticulatedSystem* anymal){
  std::cout << "STEP[" << t << "]" << std::endl;
  /// TEMPLATE (do some testing here)

  auto r = initRobot();
  // std::cout << "STEP[" << t << "] ROBOT INIT" << std::endl;

  r->calculateKinematics(gc);
  // std::cout << "STEP[" << t << "] ROBOT KIN" << std::endl;

  r->calculateCompositeInertia();
  // std::cout << "STEP[" << t << "] ROBOT COMPI" << std::endl;

  // Eigen::MatrixXd MCalc = r->calculateMassMatrix();
  Eigen::MatrixXd MCalc = getMassMatrix(gc);
  // std::cout << "STEP[" << t << "] ROBOT MASSMAT" << std::endl;

  Eigen::MatrixXd MTrue = anymal->getMassMatrix().e(); // required for other calculations

  auto inertias = anymal->getInertia();
  auto masses = anymal->getMass();
  auto coms   = anymal->getBodyCOM_B();
  
  auto names  = anymal->getBodyNames();
  auto comWs  = anymal->getBodyCOM_W();

  auto compositeInertias = anymal->getCompositeInertia();
  auto compositeMasses   = anymal->getCompositeMass();
  auto compositeComs     = anymal->getCompositeCOM();

  // for (size_t i = 0; i < inertias.size() ;i++){
  //   std::cout << "RAISIM: [" << i <<"] " << names[i] << std::endl;
  //   std::cout << "m: " << masses[i] << "  com: " << coms[i].e().transpose() << std::endl;
  //   std::cout << inertias[i] <<std::endl;
  // }

  std::string target = "base";

  auto l = r->getLinkByName(target);
  size_t bodyIdx = anymal->getBodyIdx(target);
  // marking positions
  server->getVisualObject("debug_mine")->setPosition(l->compI.com.originPos);
  server->getVisualObject("debug_raisim")->setPosition(compositeComs[bodyIdx].e());

  // std::cout << "------COMPOSITE-INERTIA------"<<std::endl;
  // std::cout << "MINE: [" << l->name << "]" << std::endl;
  // std::cout << "m: " << l->compI.m << "  com: " << l->compI.com.originPos.transpose() << std::endl;
  // std::cout << l->compI.I <<std::endl;
  //
  // std::cout << "RAISIM: [" << bodyIdx <<"] " << names[bodyIdx] << std::endl;
  // std::cout << "m: " << compositeMasses[bodyIdx] << "  com: " << compositeComs[bodyIdx].e().transpose() << std::endl;
  // std::cout << compositeInertias[bodyIdx] <<std::endl;
  //

  std::cout << "------MASS-(SUB)MATRIX ------" << std::endl;
  std::cout << "MINE" << std::endl;
  // std::cout << MCalc.block<12,12>(6,6) << std::endl;
  std::cout << MCalc << std::endl;

  std::cout << "RAISIM" << std::endl;
  // std::cout << MTrue.block<12,12>(6,6) << std::endl;
  std::cout << MTrue << std::endl;

  // std::cout << "Error : " << (MCalc - MTrue).block<12,12>(6,6).norm() << std::endl;
  auto err = (MCalc - MTrue).norm();
  std::cout << "Error : " << err << std::endl;

 std::cout << "------[SANITY-CHECK]------" << std::endl;
  raisim::Vec<3> posBase, posLF, posRF, posLH, posRH;

  posBase = anymal->getBasePosition();
  anymal->getFramePosition("LF_shank_fixed_LF_FOOT", posLF);
  anymal->getFramePosition("RF_shank_fixed_RF_FOOT", posRF);
  anymal->getFramePosition("LH_shank_fixed_LH_FOOT", posLH);
  anymal->getFramePosition("RH_shank_fixed_RH_FOOT", posRH);

  std::cout << "BASE POS (MINE)  : " << r->getPos(gc,"base").transpose() << std::endl;
  std::cout << "BASE POS (RAISIM): " << posBase.e().transpose() << std::endl;
  std::cout<<std::endl;

  std::cout << "BASE COM (MINE)    : " << r->getLinkByName("base")->compI.com.originPos.transpose() << std::endl;
  std::cout << "BASE COM (RAISIM)  : " << compositeComs[0].e().transpose() << std::endl;
  std::cout<<std::endl;

  std::cout << "FOOT POS (MINE)  : " << std::endl;
  std::cout << "   LF: " << r->getPos(gc,"LF_FOOT").transpose();
  std::cout << "   RF: " << r->getPos(gc,"RF_FOOT").transpose();
  std::cout << std::endl;
  std::cout << "   LH: " << r->getPos(gc,"LH_FOOT").transpose();
  std::cout << "   RH: " << r->getPos(gc,"RH_FOOT").transpose();
  std::cout << std::endl;

  std::cout << "FOOT POS (RAISIM): " << std::endl;
  std::cout << "   LF: " << posLF.e().transpose();
  std::cout << "   RF: " << posRF.e().transpose();
  std::cout << std::endl;
  std::cout << "   LH: " << posLH.e().transpose();
  std::cout << "   RH: " << posRH.e().transpose();
  std::cout << std::endl ;

  std::cout << std::endl;
  /// TEMPLATE (add return condition here)
  return (err < 1e-12);
}

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  auto anymal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/anymal_c/urdf/anymal.urdf");

  world.addGround();
  world.setTimeStep(0.001);

  // debug spheres
  auto debugRAISIM  = server.addVisualSphere("debug_raisim", 0.05);
  auto debugMINE    = server.addVisualSphere("debug_mine", 0.06); //when coincident, the sphere looks green
  debugRAISIM ->setColor(1,0,0,1);
  debugMINE   ->setColor(0,1,0,1);

  Eigen::VectorXd gc(anymal->getGeneralizedCoordinateDim()), gv(anymal->getDOF());
  gc << 0, 0, 0.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8; /// Jemin: I'll randomize the gc, gv when grading
  gv << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8;

  utils::gcRandomize(gc);
  gc[2] = gc[2] + 3;
  utils::gvRandomize(gv,15);
  anymal->setState(gc, gv);
  server.launchServer();
  server.focusOn(anymal);

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
    RS_TIMED_LOOP(world.getTimeStep()*2e6)
    if(i%10 == 0){
      correct = correct && analyzeStep(gc,i,&server,anymal);
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
