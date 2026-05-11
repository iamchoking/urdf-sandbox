#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <raipal_kinematics/raipal_cfb.hpp>
#include "raisimRaipal/Raipal.hpp"
#include "frame_timer.hpp"
#include "random_coordinates.hpp"

double PLAYBACK_SPEED = 1.0;
double SIM_TIMESTEP = 0.0001;
bool RANDOM_SEED = true;

size_t TEST1_NUM_POSES = 0;  // random pose test
double TEST2_DURATION  = 0.0;  // pendulum test
double TEST3_DURATION  = 0.0;  // elbow drop test (~1.0s)
double TEST4_DURATION  = 0.0;  // sine-wave joint-side test (~5.0s)
double TEST5_DURATION  = 0.0;  // sine-wave actuator-side test
double TEST6_DURATION  = 10.0;  // random actuator-side target test

namespace rk9 = raipal::kinematics;

int main(int argc, char* argv[]) {
  if(RANDOM_SEED){utils::setSeed(static_cast<unsigned>(std::time(nullptr)));}
  else {utils::setSeed(0);}

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  auto raipal9 = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal9/urdf/raipal_stub-0_R.urdf");
  // auto raipal7 = new ArticulatedRaipal(world.addArticulatedSystem(
  //   std::string(_MAKE_STR(RESOURCE_DIR)) +  "/raipal/urdf/raipal_stub-0_L.urdf")
  //   ,{3}, {-1.0}
  // );

  auto raipal7 = Raipal(world.addArticulatedSystem(
    std::string(_MAKE_STR(RESOURCE_DIR)) +  "/raipal/urdf/raipal_stub-0_L.urdf")
    ,{3}, {-1}
  );

  // raipal -> setComputeInverseDynamics(true);
  std::cout << "both models loaded!" << std::endl;  

  // world.addGround();
  world.setTimeStep(SIM_TIMESTEP);
  FrameTimer testTimer(world.getTimeStep() / PLAYBACK_SPEED, false);

  // Declare variables (should be in private section)
  int gcDim_, gvDim_;

  raipal9->setPdGains(Eigen::VectorXd::Zero(9), Eigen::VectorXd::Zero(9));
  raipal7->setPdGains(Eigen::VectorXd::Zero(7), Eigen::VectorXd::Zero(7));
  
  server.launchServer();
  // If you created `ArticulatedRaipal` as a pointer, use ->get()
  // for downstream compatibility with raisim::ArticulatedSystem*.
  // server.focusOn(raipal7->get()); 
  
  server.focusOn(raipal7);
  
  /// if you are using an old version of Raisim, you need this line
  // world.integrate1();
  
  // placeholder variables
  Eigen::VectorXd gc9(9), gv9(9), pTarget9(9), dTarget9(9), pGain9(9), dGain9(9);
  Eigen::VectorXd gc7(7), gv7(7), pTarget7(7), dTarget7(7), pGain7(7), dGain7(7);
  Eigen::VectorXd gc7Actuator(7), gv7Actuator(7);
  auto mirroredPositionDiff = [](const Eigen::VectorXd& gc9, const Eigen::VectorXd& gc7, const Eigen::VectorXd& gc7Actuator) {
    Eigen::VectorXd diff(8);
    diff.head(3) = gc9.head(3) + gc7.head(3);
    diff(3) = gc9(3) + gc7Actuator(3);
    diff(4) = gc9(5) + gc7(3);
    diff.tail(3) = gc9.tail(3) + gc7.tail(3);
    return diff.cwiseAbs();
  };
  auto mirroredEndEffectorDiff = [&]() {
    raisim::Vec<3> rightEndEffector;
    raisim::Vec<3> leftEndEffector;
    raipal9->getFramePosition("RE_end_effector_fixed", rightEndEffector);
    raipal7->getFramePosition("LE_end_effector_fixed", leftEndEffector);
    Eigen::Vector3d mirroredLeftEndEffector = leftEndEffector.e();
    mirroredLeftEndEffector.y() *= -1.0;
    return (rightEndEffector.e() - mirroredLeftEndEffector).norm();
  };
  auto runPositionTrackingTest = [&](double duration, auto&& updatePdTargets) {
    const size_t testSteps = (size_t)(duration / world.getTimeStep());

    for (int sec=3; sec>0; sec--){
      std::cout << "Starting in [" << sec << "]..." << std::endl;
      raisim::USLEEP(1000000);
    }
    std::cout << "START!" << std::endl;

    double maxPositionDiff = 0.0;
    double avgPositionDiff = 0.0;
    double maxEndEffectorDiff = 0.0;
    double avgEndEffectorDiff = 0.0;
    const size_t printEverySteps = std::max<size_t>(1, (size_t)(0.1 / world.getTimeStep()));

    testTimer.reset();
    for (size_t t = 0; t<testSteps; t++){
      testTimer.tick();

      updatePdTargets(t, testSteps);

      raipal7->updateRaipal();
      server.integrateWorldThreadSafe();
      raipal7->resetUpdateFlags();

      raipal9->getState(gc9, gv9);
      raipal7->getState(gc7, gv7);
      raipal7->getActuatorState(gc7Actuator, gv7Actuator);

      const Eigen::VectorXd positionDiff = mirroredPositionDiff(gc9, gc7, gc7Actuator);
      maxPositionDiff = std::max(maxPositionDiff, positionDiff.maxCoeff());
      avgPositionDiff += positionDiff.mean() / (double)testSteps;

      const double endEffectorDiff = mirroredEndEffectorDiff();
      maxEndEffectorDiff = std::max(maxEndEffectorDiff, endEffectorDiff);
      avgEndEffectorDiff += endEffectorDiff / (double)testSteps;

      if (t % printEverySteps == 0 || t + 1 == testSteps) {
        std::cout
          << "[STEP " << t << "]"
          // << "\n  target9: " << pTarget9.transpose()
          // << "\n  target7: " << pTarget7.transpose()
          << " positionDiff max: " << positionDiff.maxCoeff() * 180.0 / M_PI << " deg"
          << ", positionDiff avg: " << positionDiff.mean() * 180.0 / M_PI << " deg"
          << ", eeDiff: " << endEffectorDiff * 1e3 << " mm"
          << std::endl;
      }
    }
    testTimer.end();

    std::cout << "Position diff: " << 
      " max: " << maxPositionDiff * 180.0 / M_PI << " deg" << 
      " avg: " << avgPositionDiff * 180.0 / M_PI << " deg" << std::endl;

    std::cout << "End-effector pos. diff: " <<
      " max: " << maxEndEffectorDiff * 1e3 << " mm" <<
      " avg: " << avgEndEffectorDiff * 1e3 << " mm" << std::endl;
  };

  raipal9->setState(Eigen::VectorXd::Zero(9), Eigen::VectorXd::Zero(9));
  raipal7->setState(Eigen::VectorXd::Zero(7), Eigen::VectorXd::Zero(7));

  auto jointLimits9 = raipal9->getJointLimits();
  auto jointLimits7 = raipal7->getJointLimits();
  Eigen::VectorXd jointLimits9Lower, jointLimits9Upper, jointLimits9Range;

  utils::convertJointLimits(jointLimits9Lower, jointLimits9Upper, jointLimits9Range, jointLimits9);

  ///////////////// TEST0: SIMPLE DIAGNOSTICS /////////////////
  std::cout << "right gcDim: " << raipal9->getGeneralizedCoordinateDim() << std::endl;
  std::cout << "left  gcDim: " << raipal7->getGeneralizedCoordinateDim() << std::endl;
  std::cout << "Mass Matrix Diagonal" << std::endl;
  std::cout << raipal9->getMassMatrix().e().diagonal().transpose() << std::endl;
  std::cout << raipal7->getMassMatrix().e().diagonal().transpose() << std::endl;

  ///////////////// TEST1: RANDOM POSE TEST /////////////////
  if(TEST1_NUM_POSES > 0){
    std::cout << "=== Random Pose Test ===" << std::endl;
    for (int sec=3; sec>0; sec--){
      std::cout << "Starting in [" << sec << "]..." << std::endl;
      raisim::USLEEP(1000000);
    }
    std::cout << "START!" << std::endl;
  }
  else {
    std::cout << "No random pose test, skipping..." << std::endl;
  }

  for (size_t pose_idx = 0; pose_idx < TEST1_NUM_POSES; pose_idx++){
    utils::sampleJointPose(gc9, jointLimits9, 0.1);
    rk9::cfbForward(gc9);
    
    gc7 << -gc9.head(3), -gc9.tail(4);

    raipal9->setState(gc9, Eigen::VectorXd::Zero(9));
    raipal7->setState(gc7, Eigen::VectorXd::Zero(7));

    std::cout << "Pose " << pose_idx << std::endl;
    std::cout << "  gc9: " << gc9.transpose() << std::endl;
    std::cout << "  gc7: " << gc7.transpose() << std::endl;
    raisim::USLEEP(1000000);
  }

  ///////////////// TEST2: PENDULUM TEST /////////////////
  size_t test2Steps = (size_t)(TEST2_DURATION/world.getTimeStep());

  if(TEST2_DURATION > 0.0){
    std::cout << "=== Pendulum Test ===" << std::endl;

    utils::sampleJointPose(gc9, jointLimits9, 0.1);
    gc9(0) = M_PI/2;
    gc9(2) = 0.0;
    // Eigen::VectorXd gc9 = Eigen::VectorXd::Zero(9);
    gv9 = Eigen::VectorXd::Random(9) * 10.0; // random initial velocity

    rk9::cfbForward(gc9,gv9);

    gc7 << -gc9.head(3), -gc9.tail(4);
    gv7 << -gv9.head(3), -gv9.tail(4);

    raipal9->setState(gc9, gv9);
    raipal7->setState(gc7, gv7);

    for (int sec=3; sec>0; sec--){
      std::cout << "Starting in [" << sec << "]..." << std::endl;
      raisim::USLEEP(1000000);
    }
    std::cout << "START!" << std::endl;
  }
  else {
    std::cout << "No pendulum test, skipping..." << std::endl;
  }

  testTimer.reset();
  for (size_t t = 0; t<test2Steps; t++){
    testTimer.tick();

    raipal7->updateRaipal();
    server.integrateWorldThreadSafe();
    raipal7->resetUpdateFlags();
  }
  testTimer.end();

  size_t test3Steps = (size_t)(TEST3_DURATION/world.getTimeStep());

  if(TEST3_DURATION > 0.0){
    std::cout << "=== Elbow Drop Test ===" << std::endl;

    gc9 = Eigen::VectorXd::Zero(9);
    gv9 = Eigen::VectorXd::Zero(9);

    gc9(3) = jointLimits9[3][1] - 0.05; // fully flexed elbow actuator
    rk9::cfbForward(gc9, gv9);

    gc7 << -gc9.head(3), -gc9.tail(4);
    gv7 << -gv9.head(3), -gv9.tail(4);

    raipal9->setPdTarget(gc9,gv9);
    raipal7->setPdTarget(gc7,gv7);
    raipal9->setState(gc9, gv9);
    raipal7->setState(gc7, gv7);


    pGain9 << 1000, 1000, 1000, 0, 0, 0, 1000, 1000, 1000;
    dGain9 << 100 , 100 , 100 , 0, 0, 0, 100 , 100 , 100 ;
    pGain7 << pGain9.head(3), pGain9.tail(4);
    dGain7 << dGain9.head(3), dGain9.tail(4);

    raipal9->setPdGains(pGain9, dGain9);
    raipal7->setPdGains(pGain7, dGain7);

    std::cout << "Initial gc9: " << gc9.transpose() << std::endl;
    std::cout << "Initial gc7: " << gc7.transpose() << std::endl;

    for (int sec=3; sec>0; sec--){
      std::cout << "Starting in [" << sec << "]..." << std::endl;
      raisim::USLEEP(1000000);
    }
    std::cout << "START!" << std::endl;
  }
  else {
    std::cout << "No elbow drop test, skipping..." << std::endl;
  }

  Eigen::VectorXd gc9Drop(9), gv9Drop(9), gc7Drop(7), gv7Drop(7);
  double maxAbsElbowDiff = 0.0;
  double rightDropTime = std::numeric_limits<double>::quiet_NaN();
  double leftDropTime = std::numeric_limits<double>::quiet_NaN();
  const double test3StartTime = world.getWorldTime();

  raipal9->getState(gc9Drop, gv9Drop);
  raipal7->getState(gc7Drop, gv7Drop);

  double previousTime = 0.0;
  double previousRightElbow = gc9Drop(5);
  double previousLeftElbowMirrored = -gc7Drop(3);

  if (previousRightElbow <= 0.0) {
    rightDropTime = 0.0;
  }
  if (previousLeftElbowMirrored <= 0.0) {
    leftDropTime = 0.0;
  }

  testTimer.reset();
  for (size_t t = 0; t<test3Steps; t++){
    testTimer.tick();

    raipal7->updateRaipal();
    server.integrateWorldThreadSafe();
    raipal7->resetUpdateFlags();

    raipal9->getState(gc9Drop, gv9Drop);
    raipal7->getState(gc7Drop, gv7Drop);

    const double currentTime = world.getWorldTime() - test3StartTime;
    const double rightElbow = gc9Drop(5);
    const double leftElbowMirrored = -gc7Drop(3);
    const double elbowDiff = rightElbow - leftElbowMirrored;
    maxAbsElbowDiff = std::max(maxAbsElbowDiff, std::abs(elbowDiff));

    if (std::isnan(rightDropTime) && previousRightElbow > 0.0 && rightElbow <= 0.0) {
      const double alpha = previousRightElbow / (previousRightElbow - rightElbow);
      rightDropTime = previousTime + alpha * (currentTime - previousTime);
      std::cout << "Right elbow dropped at t: " << rightDropTime << std::endl;
    }
    if (std::isnan(leftDropTime) && previousLeftElbowMirrored > 0.0 && leftElbowMirrored <= 0.0) {
      const double alpha = previousLeftElbowMirrored / (previousLeftElbowMirrored - leftElbowMirrored);
      leftDropTime = previousTime + alpha * (currentTime - previousTime);
      std::cout << "Left elbow dropped at t: " << leftDropTime << std::endl;
    }

    if (t % 100 == 0 || t + 1 == test3Steps) {
      std::cout
        << "t: " << currentTime
        << ", right: " << rightElbow
        << ", left: " << leftElbowMirrored
        << ", diff: " << elbowDiff * 180.0 / M_PI << " deg"
        << std::endl;
    }

    previousTime = currentTime;
    previousRightElbow = rightElbow;
    previousLeftElbowMirrored = leftElbowMirrored;
  }
  testTimer.end();

  if(TEST3_DURATION > 0.0){
    std::cout << "Max abs elbow diff: " << maxAbsElbowDiff * 180.0 / M_PI << " deg" << std::endl;
    if (std::isnan(rightDropTime) || std::isnan(leftDropTime)) {
      std::cout
        << "Drop time difference unavailable"
        << " (right: " << (std::isnan(rightDropTime) ? -1.0 : rightDropTime)
        << ", left: " << (std::isnan(leftDropTime) ? -1.0 : leftDropTime)
        << ")" << std::endl;
    }
    else {
      std::cout
        << "Drop time difference (right - left): "
        << rightDropTime - leftDropTime
        << " s" << std::endl;
    }
  }

  
  if(TEST4_DURATION > 0.0){
    std::cout << "=== Sine-Wave Joint-Side Test ===" << std::endl;

    Eigen::VectorXd sweepCenter9(9), sweepAmplitude9(9), sweepLimits9(9), padRatio9(9), minAmplitudeRatio9(9), maxAmplitudeRatio9(9);

    padRatio9           << 0.1, 0.1, 0.1, 0.0, 0.0, 0.1, 0.1, 0.1, 0.1; // all joints have the same padding ratio
    minAmplitudeRatio9  << 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0; // only elbow joint has non-zero minimum amplitude ratio
    maxAmplitudeRatio9  << 0.1, 0.1, 0.0, 0.0, 0.0, 0.4, 0.1, 0.1, 0.1; // all joints have the same maximum amplitude ratio

    utils::sampleJointSweep(
      sweepCenter9,
      sweepAmplitude9,
      jointLimits9,
      padRatio9,
      minAmplitudeRatio9,
      maxAmplitudeRatio9
    );

    rk9::cfbBackward(sweepCenter9);

    gc9 = sweepCenter9;
    gv9.setZero();

    gc7 << -gc9.head(3), -gc9.tail(4);
    gv7 << -gv9.head(3), -gv9.tail(4);

    pGain9 << 100, 100, 100, 0, 0, 100, 100, 100, 100;
    dGain9 << 10 , 10 , 10 , 0, 0, 10 , 10 , 10 , 10 ;
    pGain7 << pGain9.head(3), pGain9.tail(4);
    dGain7 << dGain9.head(3), dGain9.tail(4);

    raipal9->setPdGains(pGain9, dGain9);
    raipal7->setPdGains(pGain7, dGain7);

    raipal9->setState(gc9, gv9);
    raipal7->setState(gc7, gv7);

    pTarget9 = gc9;
    dTarget9 = gv9;
    pTarget7 = gc7;
    dTarget7 = gv7;

    raipal9->setPdTarget(pTarget9, dTarget9);
    raipal7->setPdTarget(pTarget7, dTarget7);

    std::cout << "Initial gc9: " << gc9.transpose() << std::endl;
    std::cout << "Initial gc7: " << gc7.transpose() << std::endl;

    std::array<double,2> freq = {1.0, 5.0};
    double theta = 0.0;
    runPositionTrackingTest(TEST4_DURATION, [&](size_t t, size_t testSteps) {
      const double currentFreq = freq[0] + (freq[1] - freq[0]) * ((double)t / (double)testSteps);
      theta += 2.0 * M_PI * currentFreq * world.getTimeStep();

      pTarget9 = sweepCenter9 + std::sin(theta) * sweepAmplitude9;
      rk9::cfbBackward(pTarget9);
      pTarget7 << -pTarget9.head(3), -pTarget9.tail(4);

      raipal9->setPdTarget(pTarget9, dTarget9);
      raipal7->setPdTarget(pTarget7, dTarget7);
    });
  }
  else {
    std::cout << "No sine-wave joint-side test, skipping..." << std::endl;
  }

  
  if(TEST5_DURATION > 0.0){
    std::cout << "=== Sine-Wave Actuator-Side Test ===" << std::endl;

    Eigen::VectorXd sweepCenter9(9), sweepAmplitude9(9), sweepLimits9(9), padRatio9(9), minAmplitudeRatio9(9), maxAmplitudeRatio9(9);

    padRatio9           << 0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.1, 0.1, 0.1; // all joints have the same padding ratio
    minAmplitudeRatio9  << 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0; // only elbow joint has non-zero minimum amplitude ratio
    maxAmplitudeRatio9  << 0.2, 0.2, 0.2, 0.4, 0.0, 0.0, 0.2, 0.2, 0.2; // all joints have the same maximum amplitude ratio

    utils::sampleJointSweep(
      sweepCenter9,
      sweepAmplitude9,
      jointLimits9,
      padRatio9,
      minAmplitudeRatio9,
      maxAmplitudeRatio9
    );

    rk9::cfbForward(sweepCenter9);

    gc9 = sweepCenter9;
    gv9.setZero();

    gc7 << -gc9.head(3), -gc9.tail(4);
    gv7 << -gv9.head(3), -gv9.tail(4);

    raipal9->setState(gc9, gv9);
    raipal7->setState(gc7, gv7);

    pTarget9 = gc9;
    dTarget9 = gv9;
    pTarget7 = gc7;
    dTarget7 = gv7;

    raipal9->setPdTarget(pTarget9, dTarget9);
    raipal7->setPdTarget(pTarget7, dTarget7);

    raipal7->setCfbTargetFromActuator();

    pGain9 << 100, 100, 100, 100, 0, 0, 100, 100, 100;
    dGain9 << 10 , 10 , 10 , 10 , 0, 0, 10 , 10 , 10 ;
    pGain7 << pGain9.head(4), pGain9.tail(3);
    dGain7 << dGain9.head(4), dGain9.tail(3);

    raipal9->setPdGains(pGain9, dGain9);
    raipal7->setActuatorPdGains(pGain7, dGain7);


    std::cout << "Initial gc9: " << gc9.transpose() << std::endl;
    std::cout << "Initial gc7: " << gc7.transpose() << std::endl;

    std::array<double,2> freq = {1.0, 5.0};
    double theta = 0.0;
    runPositionTrackingTest(TEST5_DURATION, [&](size_t t, size_t testSteps) {
      const double currentFreq = freq[0] + (freq[1] - freq[0]) * ((double)t / (double)testSteps);
      theta += 2.0 * M_PI * currentFreq * world.getTimeStep();

      pTarget9 = sweepCenter9 + std::sin(theta) * sweepAmplitude9;
      rk9::cfbForward(pTarget9);
      pTarget7 << -pTarget9.head(4), -pTarget9.tail(3);

      raipal9->setPdTarget(pTarget9, dTarget9);
      raipal7->setActuatorPdTarget(pTarget7, dTarget7);
    });
  }
  else {
    std::cout << "No sine-wave actuator-side test, skipping..." << std::endl;
  }

  if(TEST6_DURATION > 0.0){
    std::cout << "=== Random Actuator-Side Target Test ===" << std::endl;

    Eigen::VectorXd padRatio9(9);
    padRatio9 << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    utils::sampleJointPose(gc9, jointLimits9, padRatio9);
    rk9::cfbForward(gc9);
    gv9.setZero();

    gc7 << -gc9.head(3), -gc9.tail(4);
    gv7 << -gv9.head(3), -gv9.tail(4);

    raipal9->setState(gc9, gv9);
    raipal7->setState(gc7, gv7);

    pTarget9 = gc9;
    dTarget9 = gv9;
    pTarget7 = gc7;
    dTarget7 = gv7;

    raipal9->setPdTarget(pTarget9, dTarget9);
    raipal7->setPdTarget(pTarget7, dTarget7);

    raipal7->setCfbTargetFromActuator();

    pGain9 << 100, 100, 100, 100, 0, 0, 100, 100, 100;
    dGain9 << 10 , 10 , 10 , 10 , 0, 0, 10 , 10 , 10 ;
    pGain7 << pGain9.head(4), pGain9.tail(3);
    dGain7 << dGain9.head(4), dGain9.tail(3);

    raipal9->setPdGains(pGain9, dGain9);
    raipal7->setActuatorPdGains(pGain7, dGain7);

    std::cout << "Initial gc9: " << gc9.transpose() << std::endl;
    std::cout << "Initial gc7: " << gc7.transpose() << std::endl;

    const double randomTargetFrequency = 5.0;
    const size_t randomTargetSteps = std::max<size_t>(
      1,
      (size_t)std::round((1.0 / randomTargetFrequency) / world.getTimeStep())
    );
    Eigen::VectorXd randomTargetStart9 = gc9;
    Eigen::VectorXd randomTargetEnd9(9);
    utils::sampleJointPose(randomTargetEnd9, jointLimits9, padRatio9);

    std::cout << "Random target frequency: " << randomTargetFrequency << " Hz" << std::endl;

    runPositionTrackingTest(TEST6_DURATION, [&](size_t t, size_t) {
      if (t > 0 && t % randomTargetSteps == 0) {
        randomTargetStart9 = randomTargetEnd9;
        utils::sampleJointPose(randomTargetEnd9, jointLimits9, padRatio9);
      }

      const double interpolationAlpha =
        (double)(t % randomTargetSteps) / (double)randomTargetSteps;
      pTarget9 = randomTargetStart9 + interpolationAlpha * (randomTargetEnd9 - randomTargetStart9);
      rk9::cfbForward(pTarget9);
      pTarget7 << -pTarget9.head(4), -pTarget9.tail(3);

      raipal9->setPdTarget(pTarget9, dTarget9);
      raipal7->setActuatorPdTarget(pTarget7, dTarget7);
    });
  }
  else {
    std::cout << "No random actuator-side target test, skipping..." << std::endl;
  }

  server.killServer();
  std::cout<<"TEST COMPLETE"<<std::endl;
  return 0;
}
