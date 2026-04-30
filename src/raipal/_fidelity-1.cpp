#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include <raipal_kinematics/raipal_cfb.hpp>
#include "ArticulatedRaipal.hpp"
#include "frame_timer.hpp"

double PLAYBACK_SPEED = 1.0;
double SIM_TIMESTEP = 0.0001;
bool RANDOM_SEED = true;

size_t TEST1_NUM_POSES = 0;  // random pose test
double TEST2_DURATION  = 0.0;  // corrected pendulum test
double TEST3_DURATION  = 0.0;  // elbow drop test
double TEST4_DURATION  = 10.0;  // sine-wave joint-side test

namespace rk9 = raipal::kinematics;

double getPlaybackTimestep(double simulationTimestep) {
  if (PLAYBACK_SPEED <= 0.0) {
    std::cout << "PLAYBACK_SPEED must be positive. Falling back to real-time playback." << std::endl;
    return simulationTimestep;
  }
  return simulationTimestep / PLAYBACK_SPEED;
}

int main(int argc, char* argv[]) {
  if(RANDOM_SEED){
    std::srand(static_cast<unsigned>(std::time(nullptr)));
  }
  else {
    std::srand(0);
  }

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  auto raipal9 = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal9/urdf/raipal_stub-0_R.urdf");
  // auto raipal7 = new ArticulatedRaipal(world.addArticulatedSystem(
  //   std::string(_MAKE_STR(RESOURCE_DIR)) +  "/raipal/urdf/raipal_stub-0_L.urdf")
  //   ,{3}, {-1.0}
  // );

  auto raipal7 = ArticulatedRaipal(world.addArticulatedSystem(
    std::string(_MAKE_STR(RESOURCE_DIR)) +  "/raipal/urdf/raipal_stub-0_L.urdf")
    ,{3}, {-1.0}
  );

  // raipal -> setComputeInverseDynamics(true);
  std::cout << "both models loaded!" << std::endl;  

  // world.addGround();
  world.setTimeStep(SIM_TIMESTEP);
  const double playbackTimestep = getPlaybackTimestep(world.getTimeStep());
  std::cout << "Playback speed: " << PLAYBACK_SPEED << "x" << std::endl;

  // Declare variables (should be in private section)
  int gcDim_, gvDim_, nJoints_;
  Eigen::VectorXd gc_init_, gv_init_, gc_, gv_, pTarget_, dTarget_, pGain_, dGain_;
  //   int obDim_ = 0, actionDim_ = 0;

  std::cout << "right gcDim: " << raipal9->getGeneralizedCoordinateDim() << std::endl;
  std::cout << "left  gcDim: " << raipal7->getGeneralizedCoordinateDim() << std::endl;

  raipal9->setPdGains(Eigen::VectorXd::Zero(9), Eigen::VectorXd::Zero(9));
  raipal7->setPdGains(Eigen::VectorXd::Zero(7), Eigen::VectorXd::Zero(7));
  
  server.launchServer();
  // If you created `ArticulatedRaipal` as a pointer, use ->get()
  // for downstream compatibility with raisim::ArticulatedSystem*.
  // server.focusOn(raipal7->get()); 
  
  server.focusOn(raipal9);
  
  /// if you are using an old version of Raisim, you need this line
  raipal7->updateRaipal();
  world.integrate1();
  raipal7->resetUpdateFlag();
  
  raipal9->setState(Eigen::VectorXd::Zero(9), Eigen::VectorXd::Zero(9));
  raipal7->setState(Eigen::VectorXd::Zero(7), Eigen::VectorXd::Zero(7));

  auto jointLimits9 = raipal9->getJointLimits();
  auto jointLimits7 = raipal7->getJointLimits();
  Eigen::VectorXd jointLimitsLower9 = Eigen::VectorXd::Zero(9);
  Eigen::VectorXd jointLimitsRange9 = Eigen::VectorXd::Zero(9);

  for (size_t i=0; i<9; i++){
    jointLimitsLower9(i) = jointLimits9[i][0];
    jointLimitsRange9(i) = jointLimits9[i][1] - jointLimits9[i][0];
  }

  std::cout << "Mass Matrix Diagonal" << std::endl;
  std::cout << raipal9->getMassMatrix().e().diagonal().transpose() << std::endl;
  std::cout << raipal7->getMassMatrix().e().diagonal().transpose() << std::endl;

  
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
    Eigen::VectorXd gc9 = jointLimitsLower9 + Eigen::VectorXd::Random(9).cwiseAbs().cwiseProduct(jointLimitsRange9);
    Eigen::VectorXd gc7 = Eigen::VectorXd::Zero(7);

    rk9::cfbForward(gc9);
    
    gc7 << -gc9.head(3), -gc9.tail(4);

    raipal9->setState(gc9, Eigen::VectorXd::Zero(9));
    raipal7->setState(gc7, Eigen::VectorXd::Zero(7));

    std::cout << "Pose " << pose_idx << std::endl;
    std::cout << "  gc9: " << gc9.transpose() << std::endl;
    std::cout << "  gc7: " << gc7.transpose() << std::endl;
    raisim::USLEEP(1000000);
  }

  // SIM LOOP


  size_t test2Steps = (size_t)(TEST2_DURATION/world.getTimeStep());

  if(TEST2_DURATION > 0.0){
    std::cout << "=== Corrected Pendulum Test ===" << std::endl;

    Eigen::VectorXd gc9 = jointLimitsLower9 + Eigen::VectorXd::Random(9).cwiseAbs().cwiseProduct(jointLimitsRange9);
    gc9(0) = M_PI/2;
    gc9(1) = -M_PI/2;
    gc9(2) = 0.0;
    // Eigen::VectorXd gc9 = Eigen::VectorXd::Zero(9);
    Eigen::VectorXd gv9 = Eigen::VectorXd::Random(9) * 10.0; // random initial velocity

    Eigen::VectorXd gc7 = Eigen::VectorXd::Zero(7);
    Eigen::VectorXd gv7 = Eigen::VectorXd::Zero(7);

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
    std::cout << "No corrected pendulum test, skipping..." << std::endl;
  }

  FrameTimer test2Timer(playbackTimestep, false);
  for (size_t t = 0; t<test2Steps; t++){
    test2Timer.tick();

    raipal7->updateRaipal();
    server.integrateWorldThreadSafe();
    raipal7->resetUpdateFlag();
  }
  test2Timer.end();


  size_t test3Steps = (size_t)(TEST3_DURATION/world.getTimeStep());

  if(TEST3_DURATION > 0.0){
    std::cout << "=== Elbow Drop Test ===" << std::endl;

    Eigen::VectorXd gc9 = Eigen::VectorXd::Zero(9);
    Eigen::VectorXd gv9 = Eigen::VectorXd::Zero(9);
    Eigen::VectorXd gc7 = Eigen::VectorXd::Zero(7);
    Eigen::VectorXd gv7 = Eigen::VectorXd::Zero(7);

    gc9(3) = jointLimits9[3][1] - 0.05; // fully flexed elbow actuator
    rk9::cfbForward(gc9, gv9);

    gc7 << -gc9.head(3), -gc9.tail(4);
    gv7 << -gv9.head(3), -gv9.tail(4);

    Eigen::VectorXd pGain9 = Eigen::VectorXd::Constant(9, 1000.0);
    Eigen::VectorXd dGain9 = Eigen::VectorXd::Constant(9, 100.0);
    Eigen::VectorXd pGain7 = Eigen::VectorXd::Constant(7, 1000.0);
    Eigen::VectorXd dGain7 = Eigen::VectorXd::Constant(7, 100.0);

    // Leave the elbow linkage free. In the 9-DOF model indices 4 and 5 are passive.
    for (int idx : {3, 4, 5}) {
      pGain9(idx) = 0.0;
      dGain9(idx) = 0.0;
    }
    pGain7(3) = 0.0;
    dGain7(3) = 0.0;

    raipal9->setPdGains(pGain9, dGain9);
    raipal7->setPdGains(pGain7, dGain7);
    raipal9->setPdTarget(Eigen::VectorXd::Zero(9), Eigen::VectorXd::Zero(9));
    raipal7->setPdTarget(Eigen::VectorXd::Zero(7), Eigen::VectorXd::Zero(7));
    raipal9->setGeneralizedForce(Eigen::VectorXd::Zero(9));

    raipal9->setState(gc9, gv9);
    raipal7->setState(gc7, gv7);

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

  FrameTimer test3Timer(playbackTimestep, false);
  for (size_t t = 0; t<test3Steps; t++){
    test3Timer.tick();

    raipal7->updateRaipal();
    server.integrateWorldThreadSafe();
    raipal7->resetUpdateFlag();

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
  test3Timer.end();

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

  size_t test4Steps = (size_t)(TEST4_DURATION/world.getTimeStep());

  if(TEST4_DURATION > 0.0){
    std::cout << "=== Sine-Wave Joint-Side Test ===" << std::endl;

    Eigen::VectorXd gc9 = jointLimitsLower9 + Eigen::VectorXd::Random(9).cwiseAbs().cwiseProduct(jointLimitsRange9);
    gc9(0) = M_PI/2;
    gc9(1) = -M_PI/2;
    gc9(2) = 0.0;
    Eigen::VectorXd gv9 = Eigen::VectorXd::Zero(9);

    const double elbowLimitMargin = M_PI / 6.0;
    const double rightElbowLower = jointLimits9[5][0];
    const double rightElbowUpper = jointLimits9[5][1];
    const double mirroredLeftElbowLower = -jointLimits7[3][1];
    const double mirroredLeftElbowUpper = -jointLimits7[3][0];
    const double sineLower = std::max(rightElbowLower, mirroredLeftElbowLower) + elbowLimitMargin;
    const double sineUpper = std::min(rightElbowUpper, mirroredLeftElbowUpper) - elbowLimitMargin;
    const double sineAmplitudeMin = elbowLimitMargin;
    const double sineAmplitudeMax = 0.5 * (sineUpper - sineLower);

    if (sineAmplitudeMax <= sineAmplitudeMin) {
      std::cout << "Sine-wave test requires more elbow joint range for a >30 deg amplitude with 30 deg limit margin." << std::endl;
      return 1;
    }

    const double sineAmplitude = sineAmplitudeMin +
      std::abs(Eigen::VectorXd::Random(1)(0)) * (sineAmplitudeMax - sineAmplitudeMin);
    const double sineCenterLower = sineLower + sineAmplitude;
    const double sineCenterUpper = sineUpper - sineAmplitude;
    const double sineCenter = sineCenterLower +
      0.5 * (Eigen::VectorXd::Random(1)(0) + 1.0) * (sineCenterUpper - sineCenterLower);
    const double sineFrequencyStart = 0.5;
    const double sineFrequencyEnd   = 2.0;
    const double sineFrequencySlope = (sineFrequencyEnd - sineFrequencyStart) / TEST4_DURATION;

    gc9(5) = sineCenter;
    gv9(5) = 0.0;

    Eigen::VectorXd gc7 = Eigen::VectorXd::Zero(7);
    Eigen::VectorXd gv7 = Eigen::VectorXd::Zero(7);
    gc7 << -gc9.head(3), -gc9.tail(4);
    gv7 << -gv9.head(3), -gv9.tail(4);

    Eigen::VectorXd pGain9 = Eigen::VectorXd::Constant(9, 1000.0);
    Eigen::VectorXd dGain9 = Eigen::VectorXd::Constant(9, 100.0);
    Eigen::VectorXd pGain7 = Eigen::VectorXd::Constant(7, 1000.0);
    Eigen::VectorXd dGain7 = Eigen::VectorXd::Constant(7, 100.0);

    pGain9(3) = 0.0;
    dGain9(3) = 0.0;
    pGain9(4) = 0.0;
    dGain9(4) = 0.0;
    pGain9(5) = 100.0;
    dGain9(5) = 10.0;
    pGain7(3) = 100.0;
    dGain7(3) = 10.0;

    
    raipal9->setPdGains(pGain9, dGain9);
    raipal7->setPdGains(pGain7, dGain7);

    rk9::cfbBackward(gc9);
    raipal9->setState(gc9, gv9);
    raipal7->setState(gc7, gv7);

    std::cout << "Sine center: " << sineCenter * 180.0 / M_PI << " deg" << std::endl;
    std::cout << "Sine amplitude: " << sineAmplitude * 180.0 / M_PI << " deg" << std::endl;
    std::cout << "Sine frequency sweep: " << sineFrequencyStart << " Hz -> " << sineFrequencyEnd << " Hz" << std::endl;
    std::cout << "Initial gc9: " << gc9.transpose() << std::endl;
    std::cout << "Initial gc7: " << gc7.transpose() << std::endl;

    for (int sec=3; sec>0; sec--){
      std::cout << "Starting in [" << sec << "]..." << std::endl;
      raisim::USLEEP(1000000);
    }
    std::cout << "START!" << std::endl;

    Eigen::VectorXd gc9Sine(9), gv9Sine(9), gc7Sine(7), gv7Sine(7);
    Eigen::VectorXd gc7ActuatorSine(7), gv7ActuatorSine(7);
    Eigen::VectorXd pTarget9 = gc9;
    Eigen::VectorXd dTarget9 = Eigen::VectorXd::Zero(9);
    Eigen::VectorXd pTarget7 = gc7;
    Eigen::VectorXd dTarget7 = Eigen::VectorXd::Zero(7);

    double maxSineElbowDiff = 0.0;
    double sumAbsSineElbowDiff = 0.0;
    double maxSineActuatorDiff = 0.0;
    double sumAbsSineActuatorDiff = 0.0;
    const double test4StartTime = world.getWorldTime();
    const size_t printEverySteps = std::max<size_t>(1, (size_t)(0.1 / world.getTimeStep()));

    FrameTimer test4Timer(playbackTimestep, false);
    for (size_t t = 0; t<test4Steps; t++){
      test4Timer.tick();

      const double currentTime = world.getWorldTime() - test4StartTime;
      const double sineFrequency = sineFrequencyStart + sineFrequencySlope * currentTime;
      const double sinePhase = 2.0 * M_PI *
        (sineFrequencyStart * currentTime + 0.5 * sineFrequencySlope * currentTime * currentTime);
      const double target = sineCenter + sineAmplitude * std::sin(sinePhase);
      const double targetVelocity = 2.0 * M_PI * sineFrequency * sineAmplitude * std::cos(sinePhase);

      pTarget9(5) = target;
      dTarget9(5) = targetVelocity;
      pTarget7(3) = -target;
      dTarget7(3) = -targetVelocity;

      raipal9->setPdTarget(pTarget9, dTarget9);
      raipal7->setPdTarget(pTarget7, dTarget7);

      raipal7->updateRaipal();
      server.integrateWorldThreadSafe();
      raipal7->resetUpdateFlag();

      raipal9->getState(gc9Sine, gv9Sine);
      raipal7->getState(gc7Sine, gv7Sine);
      raipal7->getActuatorState(gc7ActuatorSine, gv7ActuatorSine);

      const double rightElbow = gc9Sine(5);
      const double leftElbowMirrored = -gc7Sine(3);
      const double elbowDiff = rightElbow - leftElbowMirrored;
      maxSineElbowDiff = std::max(maxSineElbowDiff, std::abs(elbowDiff));
      sumAbsSineElbowDiff += std::abs(elbowDiff);

      const double rightActuator = gc9Sine(3);
      const double leftActuatorMirrored = -gc7ActuatorSine(3);
      const double actuatorDiff = rightActuator - leftActuatorMirrored;
      maxSineActuatorDiff = std::max(maxSineActuatorDiff, std::abs(actuatorDiff));
      sumAbsSineActuatorDiff += std::abs(actuatorDiff);

      if (t % printEverySteps == 0 || t + 1 == test4Steps) {
        std::cout
          << "t: " << currentTime
          << ", freq: " << sineFrequency
          << ", target: " << target
          << ", R: " << rightElbow
          << ", L: " << leftElbowMirrored
          << ", diff: " << elbowDiff * 180.0 / M_PI << " deg"
          << ", R_ACT: " << rightActuator
          << ", L_ACT: " << leftActuatorMirrored
          << ", diff_ACT: " << actuatorDiff * 180.0 / M_PI << " deg"
          << std::endl;
      }
    }
    test4Timer.end();

    std::cout << "Elbow    pos. diff: " << 
      " max: " << maxSineElbowDiff * 180.0 / M_PI << " deg" << 
      " avg: " << (sumAbsSineElbowDiff / (double)test4Steps) * 180.0 / M_PI << " deg" << std::endl;

    std::cout << "Actuator pos. diff: " << 
      " max: " << maxSineActuatorDiff * 180.0 / M_PI << " deg" << 
      " avg: " << (sumAbsSineActuatorDiff / (double)test4Steps) * 180.0 / M_PI << " deg" << std::endl;
  }
  else {
    std::cout << "No sine-wave joint-side test, skipping..." << std::endl;
  }

  server.killServer();
  std::cout<<"TEST COMPLETE"<<std::endl;
  return 0;
}
