#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "raisim/RaisimServer.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "frame_timer.hpp"
#include "raisimRaipal/Raipal.hpp"

namespace {

double PLAYBACK_SPEED = 1.0;
double SIM_TIMESTEP = 0.0001;
double JOINT_LIMIT_SPRING_MARGIN = 0.1;
double JOINT_LIMIT_SPRING_CONSTANT = 10.0;

constexpr size_t NUM_TRIALS = 20;
constexpr double PD_TARGET_DURATION = 0.2;
constexpr double ZERO_GAIN_DURATION = 10.0;
constexpr double POSITION_ERROR_RANGE = 1.0;

constexpr double POSITION_GAIN = 10.0;
constexpr double VELOCITY_GAIN = 1.0;

std::string URDF_PATH = "/raipal/urdf/raipal_stub-0_L.urdf";

double getPlaybackTimestep(double simulationTimestep) {
  if (PLAYBACK_SPEED <= 0.0) {
    std::cout << "PLAYBACK_SPEED must be positive. Falling back to real-time playback." << std::endl;
    return simulationTimestep;
  }
  return simulationTimestep / PLAYBACK_SPEED;
}

Eigen::VectorXd sampleJointPositionErrorTarget(
    const Eigen::VectorXd& currentPosition,
    const std::vector<raisim::Vec<2>>& jointLimits,
    std::mt19937& rng) {
  Eigen::VectorXd pose(jointLimits.size());
  std::uniform_real_distribution<double> errorDist(
      -POSITION_ERROR_RANGE,
      POSITION_ERROR_RANGE);

  for (size_t idx = 0; idx < jointLimits.size(); ++idx) {
    const double lower = jointLimits[idx][0];
    const double upper = jointLimits[idx][1];
    pose[idx] = std::clamp(currentPosition[idx] + errorDist(rng), lower, upper);
  }

  return pose;
}

Eigen::VectorXd jointLimitCenter(const std::vector<raisim::Vec<2>>& jointLimits) {
  Eigen::VectorXd pose(jointLimits.size());
  for (size_t idx = 0; idx < jointLimits.size(); ++idx) {
    pose[idx] = 0.5 * (jointLimits[idx][0] + jointLimits[idx][1]);
  }
  return pose;
}

void addJointLimitSpringTorque(
    Eigen::VectorXd& torque,
    const Eigen::VectorXd& position,
    const std::vector<raisim::Vec<2>>& jointLimits) {
  if (JOINT_LIMIT_SPRING_MARGIN <= 0.0 || JOINT_LIMIT_SPRING_CONSTANT == 0.0) {
    return;
  }

  const size_t count = std::min(
      jointLimits.size(),
      std::min(static_cast<size_t>(position.size()), static_cast<size_t>(torque.size())));

  for (size_t idx = 0; idx < count; ++idx) {
    const double lower = jointLimits[idx][0];
    const double upper = jointLimits[idx][1];
    if (!std::isfinite(lower) || !std::isfinite(upper) || upper <= lower) {
      continue;
    }

    const double margin = std::min(JOINT_LIMIT_SPRING_MARGIN, 0.5 * (upper - lower));
    const double lowerSpringStart = lower + margin;
    const double upperSpringStart = upper - margin;
    const double q = position[idx];

    if (q < lowerSpringStart) {
      torque[idx] += JOINT_LIMIT_SPRING_CONSTANT * (lowerSpringStart - q);
    } else if (q > upperSpringStart) {
      torque[idx] -= JOINT_LIMIT_SPRING_CONSTANT * (q - upperSpringStart);
    }
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);
  (void)binaryPath;

  raisim::World world;
  world.setTimeStep(SIM_TIMESTEP);
  const double playbackTimestep = getPlaybackTimestep(world.getTimeStep());

  raisim::RaisimServer server(&world);

  const std::string raipalUrdf = std::string(_MAKE_STR(RESOURCE_DIR)) + URDF_PATH;
  raisim::Raipal raipal(
      world.addArticulatedSystem(raipalUrdf),
      {3},
      {-1});

  raipal->setName("raipal_gravity_compensation");
  // raipal->setCfbTargetFromActuator();
  raipal->setCfbTargetFromJoint();

  const size_t dof = raipal->getDOF();
  Eigen::VectorXd currentPosition(dof);
  Eigen::VectorXd currentVelocity(dof);
  Eigen::VectorXd zeroVelocity = Eigen::VectorXd::Zero(dof);
  Eigen::VectorXd pdPositionGain = Eigen::VectorXd::Constant(dof, POSITION_GAIN);
  Eigen::VectorXd pdVelocityGain = Eigen::VectorXd::Constant(dof, VELOCITY_GAIN);
  Eigen::VectorXd zeroGain = Eigen::VectorXd::Zero(dof);

  auto jointLimits = raipal->getCurrentJointLimits();
  std::mt19937 rng(std::random_device{}());

  Eigen::VectorXd nominalPosition = jointLimitCenter(jointLimits);
  raipal->setCurrentState(nominalPosition, zeroVelocity);
  raipal->setCurrentPdTarget(nominalPosition, zeroVelocity);
  raipal->setCurrentPdGains(zeroGain, zeroGain);

  server.launchServer();
  server.focusOn(raipal.get());
  world.integrate1();

  const size_t pdSteps = std::max<size_t>(
      1,
      static_cast<size_t>(std::round(PD_TARGET_DURATION / world.getTimeStep())));
  const size_t zeroGainSteps = std::max<size_t>(
      1,
      static_cast<size_t>(std::round(ZERO_GAIN_DURATION / world.getTimeStep())));

  std::cout << "=== 7-DOF Gravity Compensation Test ===" << std::endl;
  std::cout << "PD target duration: " << PD_TARGET_DURATION << " s" << std::endl;
  std::cout << "Zero-gain duration: " << ZERO_GAIN_DURATION << " s" << std::endl;
  std::cout << "PD gains: kp=" << POSITION_GAIN << ", kd=" << VELOCITY_GAIN << std::endl;
  std::cout << "Joint-limit spring: margin=" << JOINT_LIMIT_SPRING_MARGIN
            << " rad, k=" << JOINT_LIMIT_SPRING_CONSTANT << " Nm/rad" << std::endl;
  std::cout << "Playback speed: " << PLAYBACK_SPEED << "x" << std::endl;

  FrameTimer timer(playbackTimestep, false);
  for (size_t trial = 0; trial < NUM_TRIALS; ++trial) {
    raipal->getCurrentState(currentPosition, currentVelocity);
    const Eigen::VectorXd target =
        sampleJointPositionErrorTarget(currentPosition, jointLimits, rng);

    std::cout << "Trial " << trial
              << " target: " << target.transpose() << std::endl;

    raipal->setCurrentPdTarget(target, zeroVelocity);

    for (size_t step = 0; step < pdSteps + zeroGainSteps; ++step) {
      timer.tick();

      if (step == 0) {
        raipal->setCurrentPdGains(pdPositionGain, pdVelocityGain);
      } else if (step == pdSteps) {
        raipal->setCurrentPdGains(zeroGain, zeroGain);
      }

      const Eigen::VectorXd gravityCompensation =
          raipal->getCurrentNonlinearities(world.getGravity()).e();
      Eigen::VectorXd commandedTorque = gravityCompensation;
      raipal->getCurrentState(currentPosition, currentVelocity);
      addJointLimitSpringTorque(commandedTorque, currentPosition, jointLimits);
      raipal->setCurrentGeneralizedForce(commandedTorque);

      raipal->updateRaipal();
      server.integrateWorldThreadSafe();
      raipal->resetUpdateFlags();
    }

    raipal->getCurrentState(currentPosition, currentVelocity);
    std::cout << "  final position: " << currentPosition.transpose() << std::endl;
    std::cout << "  final velocity: " << currentVelocity.transpose() << std::endl;

    raipal->setCurrentState(nominalPosition, zeroVelocity);
    raipal->setCurrentPdTarget(nominalPosition, zeroVelocity);
    raipal->setCurrentPdGains(zeroGain, zeroGain);
    raipal->resetUpdateFlags();
  }
  timer.end();

  server.killServer();
  return 0;
}
