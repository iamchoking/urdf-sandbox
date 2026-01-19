// Copyright (c) 2026

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"

#include <Eigen/Geometry>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <raipal_kinematics/raipal_cfb.hpp>
#include <raipal_kinematics/raipal_kin.hpp>

namespace rk = raipal::kinematics;

namespace {

using RaipalKin = rk::RaipalKin;

constexpr size_t kNumSteps = 8000;
std::string kTipFrameName = "RE_tip_fixed";
const std::array<int, 7> kActuatedIndices{0, 1, 2, 3, 6, 7, 8};

Eigen::VectorXd startingConfig(
  const std::vector<raisim::Vec<2>> &jointLimits,
  std::mt19937 &gen)
{
  Eigen::VectorXd gc = Eigen::VectorXd::Zero(jointLimits.size());

  // gc << 2.124693, -0.712080, -1.646676,  0.672676,  1.293469,  1.122118, -0.298123,  0.006298, -0.762178;
  // gc << 1.16089619, -0.80705221, -1.77980876,  0.20547090,  0.45189218,  0.29890528, -0.60543294, -1.00913896, -0.53606410;
  // rk::cfbForward(gc);
  // return gc;

  std::uniform_real_distribution<double> unitDist(0.0, 1.0);
  for (int idx : kActuatedIndices) {
    const double lower = jointLimits[idx][0];
    const double upper = jointLimits[idx][1];
    gc[idx] = lower + (upper - lower) * unitDist(gen);
  }
  rk::cfbForward(gc);
  return gc;
}

double orientationAngleError(const Eigen::Matrix3d &target, const Eigen::Matrix3d &estimate) {
  const Eigen::Vector3d targetZ = target.col(2);
  const Eigen::Vector3d estimateZ = estimate.col(2);
  const double cosine = std::clamp(targetZ.dot(estimateZ), -1.0, 1.0);
  return std::acos(cosine);
}

} // namespace

int main(int argc, char *argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);
  (void) binaryPath;
  std::cout << "[main] Starting IK tester" << std::endl; //_debug

  raisim::World world;
  world.setTimeStep(0.001);

  raisim::RaisimServer server(&world);

  const std::string raipalUrdfDir = std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/";
  const std::string raipalRightUrdf = raipalUrdfDir + "raipal_stub-10_R.urdf";

  auto raipalRobot = world.addArticulatedSystem(raipalRightUrdf);
  raipalRobot->setName("raipal_ik_target");

  RaipalKin solver(raipalRightUrdf);

  auto visualRobot = server.addVisualArticulatedSystem(
    "raipal_ik_visual", raipalRightUrdf, 0.1, 0.6, 0.9, 0.9);

  Eigen::VectorXd zeroVel = Eigen::VectorXd::Zero(raipalRobot->getDOF());
  auto jointLimits = raipalRobot->getJointLimits();

  std::mt19937 rng(std::random_device{}());

  Eigen::VectorXd currentGc = startingConfig(jointLimits, rng);
  Eigen::VectorXd prevGc = currentGc;

  std::vector<double> jointCenter(jointLimits.size(), 0.0);
  std::vector<double> jointAmplitude(jointLimits.size(), 0.0);
  std::vector<double> jointPhase(jointLimits.size(), 0.0);

  const double sinePeriodSteps = 1000.0;
  const double phaseIncrement = 2.0 * M_PI / sinePeriodSteps;

  for (size_t idx = 0; idx < jointLimits.size(); ++idx) {

    const double lower = jointLimits[idx][0];
    const double upper = jointLimits[idx][1];
    const double center = 0.5 * (lower + upper);
    const double amplitude = 0.5 * (upper - lower);
    jointCenter[idx] = center;
    jointAmplitude[idx] = amplitude;
    if (amplitude > 1e-6) {
      const double ratio = std::clamp((currentGc[idx] - center) / amplitude, -1.0, 1.0);
      jointPhase[idx] = std::asin(ratio);
    } else {
      jointPhase[idx] = 0.0;
    }
  }

  raipalRobot->setState(currentGc, zeroVel);
  visualRobot->setGeneralizedCoordinate(currentGc);

  server.launchServer();
  server.focusOn(raipalRobot);
  world.integrate1();
  std::cout << "[main] Server launched" << std::endl; //_debug

  size_t successCount = 0;
  size_t failureCount = 0;
  double accumulatedTimeNs = 0.0;
  double accumulatedPosError = 0.0;
  double accumulatedOriError = 0.0;
  double maxPosError = 0.0;
  double maxOriError = 0.0;
  double maxSolveNs = 0.0;
  size_t stepsExecuted = 0;
  bool aborted = false;
  std::string abortReason;

  // COUNTDOUWN
  std::cout << "====== IK TEST ======" << std::endl;
  for (int sec=3; sec>0; sec--){
    std::cout << "Starting in [" << sec << "]..." << std::endl;
    raisim::USLEEP(1000000);
  }
  std::cout << "START!" << std::endl;

  // ACTUAL LOOP
  for (size_t step = 0; step < kNumSteps; ++step) {
    std::cout << "[main] Step " << step << std::endl; //_debug
    RS_TIMED_LOOP(world.getTimeStep() * 2e6)
    ++stepsExecuted;

    // increment joints
    const double phaseAdvance = phaseIncrement * static_cast<double>(step + 1);
    for (int idx : kActuatedIndices) {
      if (jointAmplitude[idx] > 1e-6) {
        currentGc[idx] = jointCenter[idx] + jointAmplitude[idx] * std::sin(jointPhase[idx] + phaseAdvance * (1 + 0.173*idx));
      } else {
        currentGc[idx] = jointCenter[idx];
      }
      rk::cfbForward(currentGc);

    }

    raipalRobot->setState(currentGc, zeroVel);

    raisim::Vec<3> tipPosVec;
    raisim::Mat<3, 3> tipRotMat;
    raipalRobot->getFramePosition(kTipFrameName, tipPosVec);
    raipalRobot->getFrameOrientation(kTipFrameName, tipRotMat);

    const Eigen::Vector3d targetPos = tipPosVec.e();
    const Eigen::Matrix3d targetRot = tipRotMat.e();

    static Eigen::VectorXd ikSol = prevGc; // only set the first time (and let ik do it own thing afterwards)
    // Eigen::VectorXd ikSol = prevGc; // set every time

    const auto ikStart = std::chrono::high_resolution_clock::now();
    rk::RaipalKin::IkResult ikResult;
    try {
      ikResult = solver.ik(kTipFrameName,targetPos,targetRot,ikSol);
    } catch (const std::exception &e) {
      abortReason = e.what();
      aborted = true;
      ++failureCount;
      std::cout << "[main] solver.ik threw: " << abortReason << std::endl;
      break;
    }
    std::cout << "[main] solver.ik done, success=" << ikResult.finalConverged 
              << ", iterations=" << ikResult.totalIterations << std::endl; //_debug

    const auto ikEnd = std::chrono::high_resolution_clock::now();

    const double solveNs = std::chrono::duration_cast<std::chrono::nanoseconds>(ikEnd - ikStart).count();
    maxSolveNs = std::max(maxSolveNs, solveNs);

    visualRobot->setGeneralizedCoordinate(ikSol);

    const auto fkResult = solver.fk(kTipFrameName, ikSol, false);
    const double posError = (fkResult.first - targetPos).norm();
    const double oriError = orientationAngleError(targetRot, fkResult.second);

    accumulatedTimeNs += solveNs;
    accumulatedPosError += posError;
    accumulatedOriError += oriError;
    maxPosError = std::max(maxPosError, posError);
    maxOriError = std::max(maxOriError, oriError);

    if (ikResult.finalConverged && !ikResult.posOnly) {
      ++successCount;
    } else {
      ++failureCount;
      std::cout << std::fixed << std::setprecision(6)
                << "[IK-FAILED] step " << step
                << " | converged: " << (ikResult.finalConverged ? "true" : "false")
                << " | pos err (mm): " << posError * 1e3
                << " | ori err (deg): " << oriError * 180.0 / M_PI
                << " | solve (us): " << solveNs / 1e3
                << " | posOnly: " << (ikResult.posOnly ? "true" : "false")
                << std::endl;
      // std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }

    prevGc = currentGc;
    // server.integrateWorldThreadSafe();
  }

  if (aborted) {
    std::cout << "[main] IK loop aborted after " << stepsExecuted
              << " steps due to: " << abortReason << std::endl;
  }

  const double avgSolveUs = stepsExecuted > 0
    ? accumulatedTimeNs / static_cast<double>(stepsExecuted) / 1e3
    : 0.0;
  const double avgPosErrMm = stepsExecuted > 0
    ? accumulatedPosError / static_cast<double>(stepsExecuted) * 1e3
    : 0.0;
  const double avgOriErrDeg = stepsExecuted > 0
    ? accumulatedOriError / static_cast<double>(stepsExecuted) * 180.0 / M_PI
    : 0.0;

  const double maxSolveUs = maxSolveNs / 1e3;

  std::cout << "\n=== IK Test Summary ===" << std::endl;
  std::cout << "Steps            : " << stepsExecuted << " / " << kNumSteps << std::endl;
  std::cout << "Successes        : " << successCount << std::endl;
  std::cout << "Failures         : " << failureCount << std::endl;
  std::cout << "Avg solve (us)   : " << avgSolveUs << std::endl;
  std::cout << "Max solve (us)   : " << maxSolveUs << std::endl;
  std::cout << "Max pos err (mm) : " << maxPosError * 1e3 << std::endl;
  std::cout << "Avg pos err (mm) : " << avgPosErrMm << std::endl;
  std::cout << "Max ori err (deg): " << maxOriError * 180 / M_PI << std::endl;
  std::cout << "Avg ori err (deg): " << avgOriErrDeg << std::endl;

  server.killServer();
  return 0;
}

