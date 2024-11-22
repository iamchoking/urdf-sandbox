#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"

#include "fin-prep_20190673.hpp"
// HEADER START

// HEADER END

/// CHECKING SCRIPT START

bool END_IF_FAIL = true;
bool INIT_RANDOMIZE = false;
bool SUPPLY_GF      = true;
bool MAX_ERROR = 1e-9;

size_t SIM_STEPS = 4000;

bool analyzeStep(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& gf, size_t t, raisim::RaisimServer* server, raisim::ArticulatedSystem* rsrobot){
  std::cout << "STEP[" << t << "]" << std::endl;
  /// TEMPLATE (do some testing here)

  auto r = initCartpole();

  r->setState(gc,gv); // consolidated gc,gv into state var.s (ex4)
  r->setForce(gf);
  // r->calculateAll();
  r->calculateKinematics();
  r->calculateDiffKinematics();
  
  r->calculateCompositeInertia();
  r->calculateMassMatrixCRBA();

  r->calculateNonlinearTermRNE();

  r->calculateUdotABA();

  Eigen::MatrixXd MTrue = rsrobot->getMassMatrix().e(); // required for other mass-related calculations
  Eigen::MatrixXd bTrue = rsrobot->getNonlinearities({0,0,-9.81}).e(); // required for other calculations
  Eigen::MatrixXd aTrue = MTrue.inverse() * (gf-bTrue);

  auto inertias = rsrobot->getInertia();
  auto masses   = rsrobot->getMass();
  auto coms     = rsrobot->getBodyCOM_B();
  auto wComs    = rsrobot->getBodyCOM_W();
  
  auto names    = rsrobot->getBodyNames();
  auto comWs    = rsrobot->getBodyCOM_W();

  auto compositeInertias = rsrobot->getCompositeInertia();
  auto compositeMasses   = rsrobot->getCompositeMass();
  auto compositeComs     = rsrobot->getCompositeCOM();

  // for (size_t i = 0; i < inertias.size() ;i++){
  //   std::cout << "RAISIM: [" << i <<"] " << names[i] << std::endl;
  //   std::cout << "m: " << masses[i] << "  com: " << coms[i].e().transpose() << std::endl;
  //   std::cout << inertias[i] <<std::endl;
  // }

  auto err = 0;
  // std::cout << "Error : " << err << std::endl;


  std::cout << "------[SANITY-CHECK]------" << std::endl;

    raisim::Vec<3> rslink1, rslink2, rslink3, rslink4, rsDebug;
    rsrobot->getFramePosition("slider", rslink1);
    rsrobot->getFramePosition("bar_revolute", rslink2);
    rsrobot->getFramePosition("bar_revolute2", rslink3);

    raisim::Vec<3> rsvlink1, rsvlink2, rsvlink3, rsvlink4, rsvDebug;
    rsrobot->getFrameVelocity("slider", rsvlink1);
    rsrobot->getFrameVelocity("bar_revolute", rsvlink2);
    rsrobot->getFrameVelocity("bar_revolute2", rsvlink3);

    Eigen::Vector3d mylink1, mylink2, mylink3, mylink4, myee, myDebug;
    mylink1 = r->getPos("slider");
    mylink2 = r->getPos("rod");
    mylink3 = r->getPos("rod2");
    myee    = r->getPos("ee");

    Eigen::Vector3d myvlink1, myvlink2, myvlink3, myvlink4, myvee, myvDebug;
    myvlink1 = r->getVel("slider");
    myvlink2 = r->getVel("rod");
    myvlink3 = r->getVel("rod2");
    myvee    = r->getVel("ee");

    // Kinematics related
      // std::cout << "LINK POSITIONS : " << std::endl;
      // std::cout << "   1: " << mylink1.transpose() << " (err: " << (mylink1 - rslink1.e()).norm() << ")" << std::endl;
      // std::cout << "   2: " << mylink2.transpose() << " (err: " << (mylink2 - rslink2.e()).norm() << ")" << std::endl;
      // std::cout << "   3: " << mylink3.transpose() << " (err: " << (mylink3 - rslink3.e()).norm() << ")" << std::endl;
      // std::cout << "   4: " << mylink4.transpose() << " (err: " << (mylink4 - rslink4.e()).norm() << ")" << std::endl;


    // Diff. kinematics related
      // std::cout << "link1 VELOCITY : " << myvlink1.transpose() << " (err: " << (myvlink1 - rsvlink1.e()).norm() << ")" << std::endl;
      // std::cout << "LINK VELOCITIES : " << std::endl;
      // std::cout << "   1: " << myvlink1.transpose() << " (err: " << (myvlink1 - rsvlink1.e()).norm() << ")" << std::endl;
      // std::cout << "   2: " << myvlink2.transpose() << " (err: " << (myvlink2 - rsvlink2.e()).norm() << ")" << std::endl;
      // std::cout << "   3: " << myvlink3.transpose() << " (err: " << (myvlink3 - rsvlink3.e()).norm() << ")" << std::endl;
      // std::cout << "   4: " << myvlink4.transpose() << " (err: " << (myvlink4 - rsvlink4.e()).norm() << ")" << std::endl;

    // Composite Mass Related
      // std::cout << "ROBOT COM (MINE)    : " << r->getLinkByName("link1")->compI.com.originPos.transpose() << std::endl;
      // std::cout << "ROBOT COM (RAISIM)  : " << compositeComs[0].e().transpose() << std::endl;
      // std::cout<<std::endl;

      // std::cout << "MASS-MATRIX (RAISIM) :" << std::endl;
      // std::cout << MTrue << std::endl;

      // std::cout << "MASS-MATRIX (MINE) :" << std::endl;
      // std::cout << r->M << std::endl;

      // std::cout << "INVERTED MASS MATRIX : " << std::endl;
      // std::cout << r->M.inverse() << std::endl;

    // Nonlinear term related
      // std::cout << "Nonlinear Term (MINE)   : " << r->b.transpose() << std::endl;
      // std::cout << "Nonlinear Term (RAISIM) : " << bTrue.transpose() << std::endl;
      // std::cout << "Difference : " << (r->b - bTrue).transpose() << std::endl;
      // std::cout << std::endl;

      // raisim::Vec<3> rsaBase, rsaLF, rsaRF, rsaLH, rsaRH, rsaDebug;
      // // rsrobot->getFrameAcceleration("base_to_base_inertia"  , rsaBase);
      // rsrobot->getFrameAcceleration("base_to_base_inertia",rsaBase);
      // rsrobot->getFrameAcceleration("LF_shank_fixed_LF_FOOT", rsaLF  );
      // rsrobot->getFrameAcceleration("RF_shank_fixed_RF_FOOT", rsaRF  );
      // rsrobot->getFrameAcceleration("LH_shank_fixed_LH_FOOT", rsaLH  );
      // rsrobot->getFrameAcceleration("RH_shank_fixed_RH_FOOT", rsaRH  );

      // this is only useful in the first step, and gf is set to zero!
      // std::cout << "[[[Accelerations]]]" << std::endl;
      // std::cout << "input force: " << rsrobot->getGeneralizedForce().e().transpose() << std::endl;
      // std::cout << "raisim says base acc is: " << rsaBase.e().transpose() << std::endl;

      // rsrobot->getFrameAcceleration("base_LF_HAA",rsaDebug);
      // std::cout << "raisim says base_LF_HAA acc is  : " << rsaDebug.e().transpose()   << std::endl;
      // rsrobot->getFrameAcceleration("LF_HAA",rsaDebug);
      // std::cout << "raisim says LF_HAA acc is  : " << rsaDebug.e().transpose()   << std::endl;
      // std::cout << "raisim says LF_FOOT acc is  : " << rsaLF.e().transpose()   << std::endl;

      // rsrobot->getFrameAcceleration("base_RH_HAA",rsaDebug);
      // std::cout << "raisim says base_RH_HAA acc is  : " << rsaDebug.e().transpose()   << std::endl;
      // rsrobot->getFrameAcceleration("RH_HAA",rsaDebug);
      // std::cout << "raisim says RH_HAA acc is  : " << rsaDebug.e().transpose()   << std::endl;

      // std::cout << "raisim says LF_FOOT acc is  : " << rsaLF.e().transpose()   << std::endl;
      // std::cout << "raisim says RF_FOOT acc is  : " << rsaRF.e().transpose()   << std::endl;
      // std::cout << "raisim says LH_FOOT acc is  : " << rsaLH.e().transpose()   << std::endl;
      // std::cout << "raisim says RH_FOOT acc is  : " << rsaRH.e().transpose()   << std::endl;

      // std::cout << std::endl;

      // rsrobot->getFrameAcceleration("RH_KFE",rsaDebug);
      // std::cout << "raisim says RH shank acc is  : " << rsaDebug.e().transpose()   << std::endl;

    // acceleration (ABA) related
      // auto aCalc = r->udot;
      std::cout << "------ COMPUTED ACCELERATION ------" << std::endl;
      std::cout << "Input gf: " << gf.transpose() << std::endl << std::endl;
      std::cout << "MINE" << std::endl;
      // std::cout << "(easy)" << aCalc.transpose() << std::endl;
      std::cout << r->udot.transpose() << std::endl;
      std::cout << computeGeneralizedAcceleration(gc,gv,gf).transpose() << std::endl;

      std::cout << "RAISIM" << std::endl;
      std::cout << aTrue.transpose() << std::endl;

      // // std::cout << "Error : " << (bCalc - bTrue).block<12,12>(6,6).norm() << std::endl;

    // setting relevant check-visuals
      server->getVisualObject("p1_mine")->setPosition(mylink1);
      server->getVisualObject("p2_mine")->setPosition(mylink2);
      server->getVisualObject("p3_mine")->setPosition(mylink3);
      // server->getVisualObject("p4_mine")->setPosition(mylink4);

      server->getVisualObject("p0_mine")->setPosition(myee);

      // double vScale = 0.1;

      // server->getVisualObject("v1_mine")->setPosition(mylink1);
      // server->getVisualObject("v1_mine")->setOrientation(quatPointing(myvlink1));
      // server->getVisualObject("v1_mine")->setCylinderSize(0.1,vScale*((myvlink1).norm()));
      // server->getVisualObject("v2_mine")->setPosition(mylink2);
      // server->getVisualObject("v2_mine")->setOrientation(quatPointing(myvlink2));
      // server->getVisualObject("v2_mine")->setCylinderSize(0.1,vScale*((myvlink2).norm()));
      // server->getVisualObject("v3_mine")->setPosition(mylink3);
      // server->getVisualObject("v3_mine")->setOrientation(quatPointing(myvlink3));
      // server->getVisualObject("v3_mine")->setCylinderSize(0.1,vScale*((myvlink3).norm()));
      // server->getVisualObject("v4_mine")->setPosition(mylink4);
      // server->getVisualObject("v4_mine")->setOrientation(quatPointing(myvlink4));
      // server->getVisualObject("v4_mine")->setCylinderSize(0.1,vScale*((myvlink4).norm()));
    

    // // Sanity Check essentials
    double ePos    = (mylink1 - rslink1.e()).norm()+(mylink2 - rslink2.e()).norm()+(mylink3 - rslink3.e()).norm()+(mylink4 - rslink4.e()).norm();
    double eVel    = (myvlink1 - rsvlink1.e()).norm()+(myvlink2 - rsvlink2.e()).norm()+(myvlink3 - rsvlink3.e()).norm()+(myvlink4 - rsvlink4.e()).norm();
    double eMass   = (MTrue - r->M).norm();
    double eNonlin = (bTrue - r->b).norm();
    double eAcc    = (aTrue - r->udot).norm();
    std::cout << "FOOT POSITIONS ERROR : "<< ePos    << std::endl;
    std::cout << "FOOT VELOCITIES ERROR: "<< eVel    << std::endl;
    std::cout << "MASS MATRIX ERROR    : "<< eMass   << std::endl;
    std::cout << "NONLINEARITY ERROR   : "<< eNonlin << std::endl;
    std::cout << "ACCELERATION ERROR   : "<< eAcc    << std::endl;

    // if(ePos + eVel + eMass + eNonlin + eAcc < MAX_ERROR){std::cout << std::endl;}
    // else{std::cout << "!!SANITY CHECK FAILED!!" << std::endl;return false;}

  if(err >= MAX_ERROR){std::cout << "STEP FAILED" << std::endl;return false;}
  return true;
}

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  raisim::ArticulatedSystem* rsrobot;

  // rsrobot = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/2DrobotArm/robot_3D.urdf"); //vanilla urdf
  rsrobot = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/cartPole/doubleCartPole.urdf");

  rsrobot -> setComputeInverseDynamics(true); // required to get ground truth in accleleration

  // world.addGround();
  world.setTimeStep(0.001);

  // some position spheres
  // auto p1RAISIM = server.addVisualSphere("p1_raisim");
  auto p0MINE   = server.addVisualSphere("p0_mine", 0.08);
  auto p1MINE   = server.addVisualSphere("p1_mine", 0.03);
  auto p2MINE   = server.addVisualSphere("p2_mine", 0.03);
  auto p3MINE   = server.addVisualSphere("p3_mine", 0.03);
  auto p4MINE   = server.addVisualSphere("p4_mine", 0.03);
  p0MINE -> setColor(.1,.1,.1,1);
  p1MINE -> setColor(.5,.5,.5,1);
  p2MINE -> setColor(.5,.5,.5,1);
  p3MINE -> setColor(.5,.5,.5,1);
  p4MINE -> setColor(.5,.5,.5,1);

  // some velocity arrows
  // auto v0MINE   = server.addVisualArrow("v0_mine", 0.1, 0.5, 1,1,1,1);
  // auto v1MINE   = server.addVisualArrow("v1_mine", 0.1, 0.5, 1,1,1,1);
  // auto v2MINE   = server.addVisualArrow("v2_mine", 0.1, 0.5, 1,1,1,1);
  // auto v3MINE   = server.addVisualArrow("v3_mine", 0.1, 0.5, 1,1,1,1);
  // auto v4MINE   = server.addVisualArrow("v4_mine", 0.1, 0.5, 1,1,1,1);

  // debug sphere
  auto debugRAISIM  = server.addVisualSphere("debug_raisim", 0.05);
  auto debugMINE    = server.addVisualSphere("debug_mine", 0.051); //when coincident, the sphere looks green
  debugRAISIM ->setColor(1,0,0,1);
  debugMINE   ->setColor(0,1,0,1);

  Eigen::VectorXd gc(rsrobot->getGeneralizedCoordinateDim()), gv(rsrobot->getDOF()), gf(rsrobot->getDOF());
  gc << 0.0, 0.1, 0.2; /// Jemin: I'll randomize the gc, gv when grading
  gv << 0.1, 0.2, 0.3;
  gf << 0.2, 0.3, 0.4;
  if(INIT_RANDOMIZE){
    utils::gvRandomize(gc,0.3);
    utils::gvRandomize(gv,2);
    utils::gvRandomize(gf,3);
  }
  if(!SUPPLY_GF){gf.setZero();}

  rsrobot->setState(gc, gv);
  rsrobot->setGeneralizedForce(gf);

  server.launchServer();
  server.focusOn(rsrobot);

  /// if you are using an old version of Raisim, you need this line
  world.integrate1();
  bool correct = analyzeStep(gc,gv,gf,0,&server,rsrobot);

  for (int sec=5; sec>0; sec--){
    std::cout << "Dropping in [" << sec << "]..." << std::endl;
    raisim::USLEEP(1000000);
  }
  std::cout << "DROP!" << std::endl;

  // std::cout<<"mass matrix should be \n"<< rsrobot->getMassMatrix().e()<<std::endl;
  for (size_t i = 0; i<SIM_STEPS; i++){
    RS_TIMED_LOOP(world.getTimeStep()*2e6)
    if(i%10 == 0){
      correct = analyzeStep(gc,gv,gf,i,&server,rsrobot) && correct;
      if(!correct && END_IF_FAIL){break;}

      if(SUPPLY_GF && INIT_RANDOMIZE){
        utils::gvRandomize(gf,0.5);
      }
    }

    server.integrateWorldThreadSafe();

    rsrobot->getState(gc, gv);
  }

  server.killServer();

  if(correct) {
    std::cout<<"TEST PASSED"<<std::endl;
  } else {
    std::cout<<"TEST FAILED"<<std::endl;
  }

  return 0;
}

