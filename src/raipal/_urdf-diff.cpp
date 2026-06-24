#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "raisim/RaisimServer.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include <raipal_kinematics/raipal_cfb.hpp>

#include "frame_timer.hpp"
#include "raisimRaipal/Raipal.hpp"

namespace rk9 = raipal::kinematics;

static constexpr double kPlaybackSpeed = 1.0;
static constexpr double kSimTimestep = 0.001;
static constexpr double kTestDuration = 2000.0;
static constexpr double kElbowSweepHz = 0.1;
static constexpr double kElbowLower = 0.0;
static constexpr double kElbowUpper = 2.2689;

static Eigen::VectorXd makeRightPose(double elbow) {
  Eigen::VectorXd gc9(9);

  // Fixed hard-coded right-arm pose. Only the elbow output angle changes.
  gc9 << 0.35, -0.65, -1.587, 0.0, 0.0, elbow, 0.15, -0.35, 0.20;

  rk9::cfbBackward(gc9);
  return gc9;
}

static Eigen::VectorXd mirrorRightJointPoseToLeft(const Eigen::VectorXd& gc9) {
  Eigen::VectorXd gc7(7);
  gc7 << -gc9.head(3), -gc9.tail(4);
  return gc7;
}

static Eigen::VectorXd mirroredPositionDiff(
    const Eigen::VectorXd& gc9,
    const Eigen::VectorXd& gc7,
    const Eigen::VectorXd& gc7Actuator) {
  Eigen::VectorXd diff(8);
  diff.head(3) = gc9.head(3) + gc7.head(3);
  diff(3) = gc9(3) + gc7Actuator(3);
  diff(4) = gc9(5) + gc7(3);
  diff.tail(3) = gc9.tail(3) + gc7.tail(3);
  return diff.cwiseAbs();
}

static double mirroredEndEffectorDiff(
    raisim::ArticulatedSystem* raipal9,
    raisim::Raipal& raipal7) {
  raisim::Vec<3> rightEndEffector;
  raisim::Vec<3> leftEndEffector;
  raipal9->getFramePosition("RE_end_effector_fixed", rightEndEffector);
  raipal7->getFramePosition("LE_end_effector_fixed", leftEndEffector);

  Eigen::Vector3d mirroredLeftEndEffector = leftEndEffector.e();
  mirroredLeftEndEffector.y() *= -1.0;
  return (rightEndEffector.e() - mirroredLeftEndEffector).norm();
}

int main(int argc, char* argv[]) {
  (void) argc;
  (void) argv;

  raisim::World world;
  world.setTimeStep(kSimTimestep);

  raisim::Vec<3> zeroGravity;
  zeroGravity.e().setZero();
  world.setGravity(zeroGravity);

  raisim::RaisimServer server(&world);

  auto raipal9 = world.addArticulatedSystem(
      std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal9/urdf/raipal_stub-0_R.urdf");
  raisim::Raipal raipal7(
      world.addArticulatedSystem(
          std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_stub-0_L.urdf"),
      {3},
      {-1});

  raipal9->setName("raipal9_right");
  raipal7->setName("raipal7_left");
  raipal9->setPdGains(Eigen::VectorXd::Zero(9), Eigen::VectorXd::Zero(9));
  raipal7->setPdGains(Eigen::VectorXd::Zero(7), Eigen::VectorXd::Zero(7));
  raipal9->setPdTarget(Eigen::VectorXd::Zero(9), Eigen::VectorXd::Zero(9));
  raipal7->setPdTarget(Eigen::VectorXd::Zero(7), Eigen::VectorXd::Zero(7));
  raipal9->setGeneralizedForce(Eigen::VectorXd::Zero(9));
  raipal7->setGeneralizedForce(Eigen::VectorXd::Zero(7));

  const Eigen::VectorXd zeroVel9 = Eigen::VectorXd::Zero(9);
  const Eigen::VectorXd zeroVel7 = Eigen::VectorXd::Zero(7);
  Eigen::VectorXd gc9 = makeRightPose(kElbowLower);
  Eigen::VectorXd gc7 = mirrorRightJointPoseToLeft(gc9);

  raipal9->setState(gc9, zeroVel9);
  raipal7->setState(gc7, zeroVel7);

  std::cout << "=== URDF elbow diff sweep ===" << std::endl;
  std::cout << "Fixed right pose: " << gc9.transpose() << std::endl;
  std::cout << "Elbow sweep: " << kElbowLower << " to " << kElbowUpper
            << " rad at " << kElbowSweepHz << " Hz" << std::endl;
  std::cout << "Duration: " << kTestDuration << " s, dt: " << kSimTimestep << " s"
            << std::endl;

  server.launchServer();
  server.focusOn(raipal7);

  const size_t testSteps = static_cast<size_t>(kTestDuration / world.getTimeStep());
  const size_t printEverySteps =
      std::max<size_t>(1, static_cast<size_t>(0.25 / world.getTimeStep()));

  Eigen::VectorXd actualGc9(9), actualGv9(9);
  Eigen::VectorXd actualGc7(7), actualGv7(7);
  Eigen::VectorXd actuatorGc7(7), actuatorGv7(7);

  double maxPositionDiff = 0.0;
  double avgPositionDiff = 0.0;
  double maxEndEffectorDiff = 0.0;
  double avgEndEffectorDiff = 0.0;

  FrameTimer timer(world.getTimeStep() / kPlaybackSpeed, false);
  timer.reset();

  for (size_t step = 0; step <= testSteps; ++step) {
    timer.tick();

    const double time = static_cast<double>(step) * world.getTimeStep();
    const double sweep =
        0.5 * (1.0 - std::cos(2.0 * M_PI * kElbowSweepHz * time));
    const double elbow = kElbowLower + sweep * (kElbowUpper - kElbowLower);

    gc9 = makeRightPose(elbow);
    gc7 = mirrorRightJointPoseToLeft(gc9);

    raipal9->setState(gc9, zeroVel9);
    raipal7->setState(gc7, zeroVel7);

    raipal9->getState(actualGc9, actualGv9);
    raipal7->getState(actualGc7, actualGv7);
    raipal7->getActuatorState(actuatorGc7, actuatorGv7);

    const Eigen::VectorXd positionDiff =
        mirroredPositionDiff(actualGc9, actualGc7, actuatorGc7);
    const double endEffectorDiff = mirroredEndEffectorDiff(raipal9, raipal7);

    maxPositionDiff = std::max(maxPositionDiff, positionDiff.maxCoeff());
    avgPositionDiff += positionDiff.mean() / static_cast<double>(testSteps + 1);
    maxEndEffectorDiff = std::max(maxEndEffectorDiff, endEffectorDiff);
    avgEndEffectorDiff += endEffectorDiff / static_cast<double>(testSteps + 1);

    if (step % printEverySteps == 0 || step == testSteps) {
      std::cout
          << "t: " << time
          << " s, elbow target: " << elbow * 180.0 / M_PI
          << " deg, right elbow: " << actualGc9(5) * 180.0 / M_PI
          << " deg, mirrored left elbow: " << -actualGc7(3) * 180.0 / M_PI
          << " deg, pos diff max: " << positionDiff.maxCoeff() * 180.0 / M_PI
          << " deg, ee diff: " << endEffectorDiff * 1e3
          << " mm" << std::endl;
    }

    server.integrateWorldThreadSafe();
    raipal7->resetUpdateFlags();
  }

  timer.end();

  std::cout << "Position diff:"
            << " max: " << maxPositionDiff * 180.0 / M_PI << " deg"
            << " avg: " << avgPositionDiff * 180.0 / M_PI << " deg"
            << std::endl;
  std::cout << "End-effector diff:"
            << " max: " << maxEndEffectorDiff * 1e3 << " mm"
            << " avg: " << avgEndEffectorDiff * 1e3 << " mm"
            << std::endl;

  server.killServer();
  return 0;
}
