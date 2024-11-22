//e
// Created by Jemin Hwangbo on 2022/04/08.
//

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"

#include "exercise_4_20190673.hpp"
// HEADER START
// HEADER END

/// CHECKING SCRIPT START

bool analyzeStep(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, size_t t, raisim::RaisimServer* server, raisim::ArticulatedSystem* anymal){
  std::cout << "STEP[" << t << "]" << std::endl;
  /// TEMPLATE (do some testing here)

  auto r = initRobot();
  r->setState(gc,gv); // consolidated gc,gv into state var.s (ex4)
  r->calculateKinematics(); // ex1 part
  // std::cout << "kin" << std::endl;
  r->calculateDiffKinematics(); // ex2 part
  //
  r->calculateCompositeInertia(); // ex3 part
  // std::cout << "comp" << std::endl;
  r->calculateMassMatrixCRBA();
  // std::cout << "mass" << std::endl;
  // DONT'DO [r->calculateAccelerations();] !! (condition-dependent)
  // std::cout << "acc" << std::endl;
  r->calculateNonlinearTermRNE();
  // std::cout << "robot's size is" << r->links.size() << std::endl;

  // Eigen::MatrixXd bCalc = r->calculateMassMatrixCRBA();
  Eigen::VectorXd bCalc = getNonlinearities(gc,gv);
  // bCalc = r->b;
  // std::cout << "STEP[" << t << "] ROBOT MASSMAT" << std::endl;

  auto MTrue = anymal->getMassMatrix(); // required for other mass-related calculations
  Eigen::MatrixXd bTrue = anymal->getNonlinearities({0,0,-9.81}).e(); // required for other calculations

  auto inertias = anymal->getInertia();
  auto masses = anymal->getMass();
  auto coms   = anymal->getBodyCOM_B();
  auto wComs  = anymal->getBodyCOM_W();
  
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

  // std::cout << "------COMPOSITE-INERTIA------"<<std::endl;
  // std::cout << "MINE: [" << l->name << "]" << std::endl;
  // std::cout << "m: " << l->compI.m << "  com: " << l->compI.com.originPos.transpose() << std::endl;
  // std::cout << l->compI.I <<std::endl;
  //
  // std::cout << "RAISIM: [" << bodyIdx <<"] " << names[bodyIdx] << std::endl;
  // std::cout << "m: " << compositeMasses[bodyIdx] << "  com: " << compositeComs[bodyIdx].e().transpose() << std::endl;
  // std::cout << compositeInertias[bodyIdx] <<std::endl;
  //

  // std::cout << "actual acc: " << anymal->getGeneralizedAcceleration() << std::endl;


  std::cout << "------ NONLINEAR TERM ------" << std::endl;
  std::cout << "MINE" << std::endl;
  // std::cout << bCalc.block<12,12>(6,6) << std::endl;
  std::cout << bCalc.transpose() << std::endl;

  std::cout << "RAISIM" << std::endl;
  // std::cout << bTrue.block<12,12>(6,6) << std::endl;
  std::cout << bTrue.transpose() << std::endl;

  // std::cout << "Error : " << (bCalc - bTrue).block<12,12>(6,6).norm() << std::endl;
  auto err = (bCalc - bTrue).norm();
  std::cout << "Error : " << err << std::endl;



  // anymal->getFramePosition("LH_HAA", rsDebug);
  // std::cout << "[DEBUG] Position: " << rsDebug.e().transpose() << std::endl;

  // std::cout << "BASE POS (MINE)  : " << r->getPos("base").transpose() << std::endl;
  // std::cout << "BASE POS (RAISIM): " << rsBase.e().transpose() << std::endl;
  // std::cout<<std::endl;

  // std::cout << "ROBOT COM (MINE)    : " << r->getLinkByName("base")->compI.com.originPos.transpose() << std::endl;
  // std::cout << "ROBOT COM (RAISIM)  : " << compositeComs[0].e().transpose() << std::endl;
  // std::cout<<std::endl;


  std::cout << "------[SANITY-CHECK]------" << std::endl;

    raisim::Vec<3> rsBase, rsLF, rsRF, rsLH, rsRH, rsDebug;
    rsBase = anymal->getBasePosition();
    anymal->getFramePosition("LF_shank_fixed_LF_FOOT", rsLF);
    anymal->getFramePosition("RF_shank_fixed_RF_FOOT", rsRF);
    anymal->getFramePosition("LH_shank_fixed_LH_FOOT", rsLH);
    anymal->getFramePosition("RH_shank_fixed_RH_FOOT", rsRH);

    raisim::Vec<3> rsvBase, rsvLF, rsvRF, rsvLH, rsvRH, rsvDebug;
    rsvBase = gv.segment(0,3);
    anymal->getFrameVelocity("LF_shank_fixed_LF_FOOT", rsvLF);
    anymal->getFrameVelocity("RF_shank_fixed_RF_FOOT", rsvRF);
    anymal->getFrameVelocity("LH_shank_fixed_LH_FOOT", rsvLH);
    anymal->getFrameVelocity("RH_shank_fixed_RH_FOOT", rsvRH);

    raisim::Vec<3> rsaBase, rsaLF, rsaRF, rsaLH, rsaRH, rsaDebug;
    // anymal->getFrameAcceleration("base_to_base_inertia"  , rsaBase);
    anymal->getFrameAcceleration("base_to_base_inertia",rsaBase);
    anymal->getFrameAcceleration("LF_shank_fixed_LF_FOOT", rsaLF  );
    anymal->getFrameAcceleration("RF_shank_fixed_RF_FOOT", rsaRF  );
    anymal->getFrameAcceleration("LH_shank_fixed_LH_FOOT", rsaLH  );
    anymal->getFrameAcceleration("RH_shank_fixed_RH_FOOT", rsaRH  );


    // std::cout << "raisim says base acc is: " << rsaBase.e().transpose() << std::endl;
    // anymal->getFrameAcceleration("LF_KFE",rsaDebug);
    // std::cout << "raisim says LF shank acc is  : " << rsaDebug.e().transpose()   << std::endl;
    // anymal->getFrameAcceleration("RH_KFE",rsaDebug);
    // std::cout << "raisim says RH shank acc is  : " << rsaDebug.e().transpose()   << std::endl;
    // anymal->getFrameAcceleration("LF_HAA",rsaDebug);
    // std::cout << "raisim says LF hip acc is  : " << rsaDebug.e().transpose()   << std::endl;
    // anymal->getFrameAcceleration("RH_HAA",rsaDebug);
    // std::cout << "raisim says RH hip acc is  : " << rsaDebug.e().transpose()   << std::endl;


    //CURSED
      // int id = 8;
      // anymal->setComputeInverseDynamics(true);
      // std::cout <<"raisim says force at LF shank is " << anymal->getForceAtJointInWorldFrame(id).e().transpose() << " " << anymal->getTorqueAtJointInWorldFrame(id).e().transpose() <<std::endl;

      // id = 7;
      // std::cout <<"raisim says force at LF thigh is " << anymal->getForceAtJointInWorldFrame(id).e().transpose() << " " << anymal->getTorqueAtJointInWorldFrame(id).e().transpose() <<std::endl;
    //


    Eigen::Vector3d myBase, myLF, myRF, myLH, myRH, myDebug;
    myBase = r->getPos("base");
    myLF   = r->getPos("LF_FOOT");
    myRF   = r->getPos("RF_FOOT");
    myLH   = r->getPos("LH_FOOT");
    myRH   = r->getPos("RH_FOOT");

    Eigen::Vector3d myvBase, myvLF, myvRF, myvLH, myvRH, myvDebug;
    myvBase = r->getVel("base");
    myvLF   = r->getVel("LF_FOOT");
    myvRF   = r->getVel("RF_FOOT");
    myvLH   = r->getVel("LH_FOOT");
    myvRH   = r->getVel("RH_FOOT");

    // CURSED
    // Eigen::Vector3d myaBase, myaLF, myaRF, myaLH, myaRH, myaDebug;
    // Eigen::VectorXd gvdot = (anymal -> getGeneralizedAcceleration()).e();
    // Eigen::VectorXd worldAcc;
    // worldAcc.resize(6);
    // worldAcc << 0,0,9.81,0,0,0;
    // r->calculateAccelerations(worldAcc,gvdot);
    // myaBase = r->getAcc("base");
    // myaLF   = r->getAcc("LF_FOOT");
    // myaRF   = r->getAcc("RF_FOOT");
    // myaLH   = r->getAcc("LH_FOOT");
    // myaRH   = r->getAcc("RH_FOOT");
    // r->resetAcc();

    // anymal->getFrameVelocity("base_LF_HAA", rsvDebug);
    // myvDebug = r->getLinkByName("LF_HIP")->fullV;

    // anymal->getFrameAngularVelocity("base_LF_HAA",rsvDebug);
    // myvDebug = r->getLinkByName("LF_HIP")->worldOm;

    // anymal->getFrameAngularVelocity("LF_HAA",rsvDebug);
    // myvDebug = r->getLinkByName("LF_HIP")->fullOm;
    // std::cout << "Debug VELOCITY : " << myvDebug.transpose() << " (err: " << (myvDebug - rsvDebug.e()).norm() << ")" << std::endl;


    // CURSED
    // anymal->getFrameAcceleration("LF_HAA",rsaDebug);
    // myaDebug = r->getLinkByName("LF_HIP")->fullA;
    // std::cout << "Debug ACCELERATION : " << myaDebug.transpose() << " (err: " << (myaDebug - rsaDebug.e()).norm() << ")" << std::endl;

    // std::cout << "FOOT POSITIONS : " << std::endl;
    // std::cout << "   LF: " << myLF.transpose() << " (err: " << (myLF - rsLF.e()).norm() << ")" << std::endl;
    // std::cout << "   RF: " << myRF.transpose() << " (err: " << (myRF - rsRF.e()).norm() << ")" << std::endl;
    // std::cout << "   LH: " << myLH.transpose() << " (err: " << (myLH - rsLH.e()).norm() << ")" << std::endl;
    // std::cout << "   RH: " << myRH.transpose() << " (err: " << (myRH - rsRH.e()).norm() << ")" << std::endl;

    // std::cout << "Base VELOCITY : " << myvBase.transpose() << " (err: " << (myvBase - rsvBase.e()).norm() << ")" << std::endl;
    // std::cout << "FOOT VELOCITIES : " << std::endl;
    // std::cout << "   LF: " << myvLF.transpose() << " (err: " << (myvLF - rsvLF.e()).norm() << ")" << std::endl;
    // std::cout << "   RF: " << myvRF.transpose() << " (err: " << (myvRF - rsvRF.e()).norm() << ")" << std::endl;
    // std::cout << "   LH: " << myvLH.transpose() << " (err: " << (myvLH - rsvLH.e()).norm() << ")" << std::endl;
    // std::cout << "   RH: " << myvRH.transpose() << " (err: " << (myvRH - rsvRH.e()).norm() << ")" << std::endl;

    // CURSED
    // std::cout << "Base ACCELERATION : " << myaBase.transpose() << " (err: " << (myaBase - rsaBase.e()).norm() << ")" << std::endl;
    // std::cout << "Base ACCELERATION*: " << rsaBase.e().transpose() << std::endl;

    // std::cout << "FOOT Acclerations : " << std::endl;
    // std::cout << "   LF: " << myaLF.transpose() << " (err: " << (myaLF - rsaLF.e()).norm() << ")" << std::endl;
    // std::cout << "   RF: " << myaRF.transpose() << " (err: " << (myaRF - rsaRF.e()).norm() << ")" << std::endl;
    // std::cout << "   LH: " << myaLH.transpose() << " (err: " << (myaLH - rsaLH.e()).norm() << ")" << std::endl;
    // std::cout << "   RH: " << myaRH.transpose() << " (err: " << (myaRH - rsaRH.e()).norm() << ")" << std::endl;

    // setting relevant check-visuals
      server->getVisualObject("p0_mine")->setPosition(myBase);
      server->getVisualObject("p1_mine")->setPosition(myLF);
      server->getVisualObject("p2_mine")->setPosition(myRF);
      server->getVisualObject("p3_mine")->setPosition(myLH);
      server->getVisualObject("p4_mine")->setPosition(myRH);

      double vScale = 0.1;
      server->getVisualObject("v0_mine")->setPosition(myBase);
      server->getVisualObject("v0_mine")->setOrientation(quatPointing(myvBase));
      server->getVisualObject("v0_mine")->setCylinderSize(0.1,vScale*(myvBase.norm()));

      server->getVisualObject("v1_mine")->setPosition(myLF);
      server->getVisualObject("v1_mine")->setOrientation(quatPointing(myvLF - myvBase));
      server->getVisualObject("v1_mine")->setCylinderSize(0.1,vScale*((myvLF - myvBase).norm()));
      server->getVisualObject("v2_mine")->setPosition(myRF);
      server->getVisualObject("v2_mine")->setOrientation(quatPointing(myvRF - myvBase));
      server->getVisualObject("v2_mine")->setCylinderSize(0.1,vScale*((myvRF - myvBase).norm()));
      server->getVisualObject("v3_mine")->setPosition(myLH);
      server->getVisualObject("v3_mine")->setOrientation(quatPointing(myvLH - myvBase));
      server->getVisualObject("v3_mine")->setCylinderSize(0.1,vScale*((myvLH - myvBase).norm()));
      server->getVisualObject("v4_mine")->setPosition(myRH);
      server->getVisualObject("v4_mine")->setOrientation(quatPointing(myvRH - myvBase));
      server->getVisualObject("v4_mine")->setCylinderSize(0.1,vScale*((myvRH - myvBase).norm()));
    //

    double ePos = (myLF - rsLF.e()).norm()+(myRF - rsRF.e()).norm()+(myLH - rsLH.e()).norm()+(myRH - rsRH.e()).norm();
    double eVel = (myvLF - rsvLF.e()).norm()+(myvRF - rsvRF.e()).norm()+(myvLH - rsvLH.e()).norm()+(myvRH - rsvRH.e()).norm();
    double eMass= (MTrue.e() - r->M).norm();
    std::cout << "FOOT POSITIONS ERROR : "<< ePos  << std::endl;
    std::cout << "FOOT VELOCITIES ERROR: "<< eVel  << std::endl;
    std::cout << "MASS MATRIX ERROR    : "<< eMass << std::endl;

    if(std::max(std::max(ePos,eVel),eMass) < 1e-12){std::cout << std::endl;}
    else{std::cout << "!!SANITY CHECK FAILED!!" << std::endl;}


  /// TEMPLATE (add return condition here)
  if(err >= 1e-10){std::cout << "STEP FAILED" << std::endl;}
  return (err < 1e-10);
}

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  auto anymal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/anymal_c/urdf/anymal.urdf");
  anymal -> setComputeInverseDynamics(true); // required to get ground truth in accleleration

  world.addGround();
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
  auto v0MINE   = server.addVisualArrow("v0_mine", 0.1, 0.5, 1,1,1,1);
  auto v1MINE   = server.addVisualArrow("v1_mine", 0.1, 0.5, 1,1,1,1);
  auto v2MINE   = server.addVisualArrow("v2_mine", 0.1, 0.5, 1,1,1,1);
  auto v3MINE   = server.addVisualArrow("v3_mine", 0.1, 0.5, 1,1,1,1);
  auto v4MINE   = server.addVisualArrow("v4_mine", 0.1, 0.5, 1,1,1,1);

  // debug sphere
  auto debugRAISIM  = server.addVisualSphere("debug_raisim", 0.05);
  auto debugMINE    = server.addVisualSphere("debug_mine", 0.06); //when coincident, the sphere looks green
  debugRAISIM ->setColor(1,0,0,1);
  debugMINE   ->setColor(0,1,0,1);

  Eigen::VectorXd gc(anymal->getGeneralizedCoordinateDim()), gv(anymal->getDOF());
  gc << 0, 0, 0.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8;
  gv << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8;

  utils::gcRandomize(gc);
  gc[2] = gc[2] + 3;
  utils::gvRandomize(gv,15);
  anymal->setState(gc, gv);
  server.launchServer();
  server.focusOn(anymal);

  /// if you are using an old version of Raisim, you need this line
  world.integrate1();
  bool correct = analyzeStep(gc,gv,0,&server,anymal);

  for (int sec=5; sec>0; sec--){
    std::cout << "Dropping in [" << sec << "]..." << std::endl;
    raisim::USLEEP(1000000);
  }
  std::cout << "DROP!" << std::endl;

  // std::cout<<"mass matrix should be \n"<< anymal->getMassMatrix().e()<<std::endl;
  for (size_t i = 0; i<2000; i++){
    RS_TIMED_LOOP(world.getTimeStep()*2e6)
    if(i%10 == 0){
      correct = analyzeStep(gc,gv,i,&server,anymal) && correct;
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
