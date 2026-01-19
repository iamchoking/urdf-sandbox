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

namespace rk = raipal::kinematics;

namespace {

class RaipalKin {
public:
  explicit RaipalKin(const std::string &urdfPath)
    : world_(std::make_unique<raisim::World>())
  {
    robot_ = world_->addArticulatedSystem(urdfPath);
    if (!robot_) {
      throw std::runtime_error("Failed to create ArticulatedSystem from URDF: " + urdfPath);
    }
    gcDim_ = robot_->getGeneralizedCoordinateDim();
    jointLimits_ = robot_->getJointLimits();

        set(Eigen::VectorXd::Zero(9));
  }

  void set(const Eigen::VectorXd gc, const bool doForward = true){

    if (gc.size() == 7) {gc_ << gc.head(3), 0.0, 0.0, gc.tail(4);} 
    else {gc_ = gc;}

    if(doForward){rk::cfbForward(gc_);}

    robot_ -> setGeneralizedCoordinate(gc_);
  }

  std::pair<Eigen::Vector3d, Eigen::Matrix3d> fk(
    const std::string &frameName,
    const Eigen::VectorXd &gc,
    const bool doForward = true
  )
  {
    set(gc,doForward);
    return fk(frameName);
  }

  std::pair<Eigen::Vector3d, Eigen::Matrix3d> fk(const std::string &frameName)
  {
    raisim::Vec<3> position;
    raisim::Mat<3, 3> orientation;

    robot_->getFramePosition(frameName, position);
    robot_->getFrameOrientation(frameName, orientation);
    return {position.e(), orientation.e()};
  }

  struct IkSettings {
    int maxMultisteps = 10;
    int maxIterations = 10;
    int maxRetries = 10; // maximum number of null-space retries before failing

    double posTolerance = 1e-4;
    double rotTolerance = 1e-3;
    double damping = 1e-3;

    double maxAlpha = 0.5; // what proportion of the step should be incremented
    double minAlpha = 1e-5; // at this point, give up on the step

    double maxPosDistance = 1e-2; // maximum position difference to solve for (m)
    double maxRotDistance = 0.17; // maximum rotational difference to solve for (rad)
    // double maxPosDistance = 5e-4; // maximum position difference to solve for (m)
    // double maxRotDistance = 5e-3; // maximum rotational difference to solve for (rad)

    double maxJointDistance = 5.0 * M_PI / 180.0; // maximum joint angle difference (rad)

    Eigen::VectorXd gcIkCenter = (Eigen::VectorXd(7) << 0.0, 0.523, 0.0, 1.0472, 0.0, 0.0, 0.0).finished();
    double nullSpaceGain = 0.000;
  };

  void setIkSettings(const IkSettings &settings) { ikSettings_ = settings; }

  struct IkResult {
    bool finalConverged = false;
    bool posOnly = false;
    int totalIterations = 0;
    double posError = 0.0;
    double rotError = 0.0;
  };

  IkResult ik(
    const std::string &frameName,
    const Eigen::Vector3d &targetPos,
    const Eigen::Matrix3d &targetRot,
    Eigen::VectorXd &sol,
    const bool posOnly = false
  ){
    const IkSettings &cfg = ikSettings_;
    IkResult result;
    result.posOnly = posOnly;
    
    // Initialize solution
    const size_t frameIdx = robot_->getFrameIdxByName(frameName);
    if(sol != Eigen::VectorXd::Zero(9) && sol != Eigen::VectorXd::Zero(7)) {
      set(sol,false);
    }
    sol = gc_;
    
    std::cout << "\n========== IK SOLVE START ==========" << std::endl;
    std::cout << "Initial config: " << std::setprecision(6) << gc_.transpose() << std::endl;
    std::cout << "Target position: " << targetPos.transpose() << std::endl;
    
    // Multistep outer loop
    bool finalConverged = false;
    
    for (int multistepIter = 0; multistepIter < cfg.maxMultisteps; ++multistepIter) {
      std::cout << "\n--- Multistep iteration " << multistepIter << " ---" << std::endl;
      
      // Calculate current error to final target
      auto [currentPos, currentRot] = fk(frameName);
      Eigen::Vector3d totalPosError = targetPos - currentPos;
      Eigen::Vector3d totalRotError = rotationErrorQuat(currentRot, targetRot);
      
      std::cout << "Distance to target: pos=" << totalPosError.norm() 
                << "m, rot=" << totalRotError.norm() << "rad" << std::endl;
      
      // Check if we've reached the final target
      if(totalPosError.norm() < cfg.posTolerance && totalRotError.norm() < cfg.rotTolerance) {
        finalConverged = true;
        result.finalConverged = true;
        std::cout << ">>> CONVERGED TO FINAL TARGET <<<" << std::endl;
        break;
      }
      
      // Determine step target (may be intermediate)
      Eigen::Vector3d stepTargetPos;
      Eigen::Matrix3d stepTargetRot;
      double scaleFactor = 1.0;
      bool isIntermediateStep = false;
      
      if(totalPosError.norm() > cfg.maxPosDistance) {
        scaleFactor = cfg.maxPosDistance / totalPosError.norm();
        isIntermediateStep = true;
      }
      if(totalRotError.norm() > cfg.maxRotDistance) {
        scaleFactor = std::min(scaleFactor, cfg.maxRotDistance / totalRotError.norm());
        isIntermediateStep = true;
      }
      
      if(isIntermediateStep) {
        scaleFactor = std::min(0.5, scaleFactor);
        stepTargetPos = currentPos + totalPosError * scaleFactor;
        stepTargetRot = applyRotationError(currentRot, totalRotError * scaleFactor);
        std::cout << "Taking intermediate step (scale=" << scaleFactor << ")" << std::endl;
      } else {
        stepTargetPos = targetPos;
        stepTargetRot = targetRot;
        std::cout << "Taking final step to target" << std::endl;
      }
      
      // Newton iteration loop for current step
      Eigen::Vector3d dPos = stepTargetPos - currentPos;
      Eigen::Vector3d dRot = rotationErrorQuat(currentRot, stepTargetRot);
      bool stepConverged = false;
      int retries = 0;
      
      for (int iter = 0; iter < cfg.maxIterations; ++iter) {
        ++result.totalIterations;
        std::cout << "  Iter " << iter << ": error pos=" << dPos.norm() 
                  << "m, rot=" << dRot.norm() << "rad" << std::endl;
        
        // Check convergence for this step
        if(dPos.norm() < cfg.posTolerance && dRot.norm() < cfg.rotTolerance) {
          stepConverged = true;
          std::cout << "  Step converged at iteration " << iter << std::endl;
          break;
        }
        
        // Compute Joint Step
        Eigen::MatrixXd J_pos9, J_rot9, J_7, Jplus, Nproj;
        Eigen::VectorXd jointStep;

        if(!posOnly){
          J_pos9.setZero(3,9);
          J_rot9.setZero(3,9);
          J_7.setZero(6,7);
          
          robot_->getDenseFrameJacobian(frameIdx, J_pos9);
          robot_->getDenseFrameRotationalJacobian(frameIdx, J_rot9);
          
          J_7.topRows(3) << J_pos9.leftCols(3), J_pos9.rightCols(4);
          J_7.bottomRows(3) << J_rot9.leftCols(3), J_rot9.rightCols(4);
          
          // Compute joint step with null-space projection
          Eigen::VectorXd taskError(6);
          taskError << dPos, dRot;
          
          const Eigen::MatrixXd JT = J_7.transpose();
          Eigen::MatrixXd JJt = J_7 * JT;
          JJt.diagonal().array() += cfg.damping * cfg.damping;
          
          Jplus = JT * JJt.inverse();
          Nproj = Eigen::MatrixXd::Identity(7,7) - Jplus * J_7;
          
          Eigen::VectorXd gc7__(7);
          gc7__ << gc_.head(3), gc_.tail(4);
          
          jointStep = Jplus * taskError + cfg.nullSpaceGain * Nproj * (cfg.gcIkCenter - gc7__);

        } else {

          J_pos9.setZero(3,9);
          J_7.setZero(3,7);
          
          robot_->getDenseFrameJacobian(frameIdx, J_pos9);
          
          J_7 << J_pos9.leftCols(3), J_pos9.rightCols(4);
          
          // Compute joint step with null-space projection
          Eigen::VectorXd taskError(3);
          taskError << dPos;
          
          const Eigen::MatrixXd JT = J_7.transpose();
          Eigen::MatrixXd JJt = J_7 * JT;
          JJt.diagonal().array() += cfg.damping * cfg.damping;
          
          Jplus = JT * JJt.inverse();
          Nproj = Eigen::MatrixXd::Identity(7,7) - Jplus * J_7;
          
          Eigen::VectorXd gc7__(7);
          gc7__ << gc_.head(3), gc_.tail(4);
          
          jointStep = Jplus * taskError + cfg.nullSpaceGain * Nproj * (cfg.gcIkCenter - gc7__);
        }

        // Line search with null-space retry
        double alpha = cfg.maxAlpha;
        bool stepAccepted = false;
        
        while (!stepAccepted) {
          Eigen::VectorXd candidateGc = applyJointStep(jointStep, alpha);
          
          // Check feasibility: joint limits and improvement
          if(
            (candidateGc - sol).cwiseAbs().maxCoeff() <= cfg.maxJointDistance &&
            isImprovement(candidateGc, stepTargetPos, stepTargetRot, dPos, dRot, frameName, cfg, posOnly)
          ) {
            stepAccepted = true;
            sol = gc_;
            std::cout << "    Accepted step: alpha=" << alpha << std::endl;
          } else {
            alpha *= 0.5;
            if(alpha < cfg.minAlpha) {
              // Null-space perturbation retry
              if(++retries > cfg.maxRetries) {
                std::cout << "  FAILED: Max retries exceeded" << std::endl;
                stepConverged = false;
                goto finish_newton_loop;
              }
              std::cout << "    Retry " << retries << ": null-space perturbation" << std::endl;
              Eigen::VectorXd randomFactor(7);
              randomFactor.setRandom();
              Eigen::VectorXd deviation = Nproj * randomFactor;
              jointStep += deviation * (jointStep.norm() / randomFactor.norm() * 0.1);
              alpha = cfg.maxAlpha;
            }
          }
        }
      }
      
      finish_newton_loop:
      if(!stepConverged) {
        std::cout << "Step FAILED to converge" << std::endl;
        break;
      }
      
      ++multistepIter;
    }
    
    // rk::cfbBackward(sol);
    
    // Compute final errors
    auto [finalPos, finalRot] = fk(frameName);
    result.posError = (targetPos - finalPos).norm();
    result.rotError = rotationErrorQuat(finalRot, targetRot).norm();
    
    std::cout << "\n========== IK SOLVE END ==========" << std::endl;
    std::cout << "Result: " << (finalConverged ? "SUCCESS" : "FAILED") << std::endl;
    std::cout << "Final config: " << sol.transpose() << std::endl;
    std::cout << "Total iterations: " << result.totalIterations << std::endl;
    std::cout << "Final pos error: " << result.posError << "m, rot error: " << result.rotError << "rad" << std::endl;
    
    if(!finalConverged && !posOnly){
      std::cout << "\nRETRYING IK WITH POSITION ONLY" << std::endl;
      return ik(frameName, targetPos, targetRot, sol, true);
    }
    return result;
  }

private:
  Eigen::Vector3d rotationErrorQuat(
    const Eigen::Matrix3d &current,
    const Eigen::Matrix3d &target)
  {
    Eigen::Quaterniond qCurrent(current);
    Eigen::Quaterniond qTarget(target);
    // delta rotates "current" into "target"
    Eigen::Quaterniond delta = qTarget * qCurrent.conjugate();
    if (delta.w() < 0.0) {
      delta.coeffs() *= -1.0; // keep the short arc
    }
    delta.normalize();

    const double sinHalf = delta.vec().norm();
    if (sinHalf < 1e-12) {
      return Eigen::Vector3d::Zero();
    }

    const double angle = 2.0 * std::atan2(sinHalf, delta.w());
    const Eigen::Vector3d axis = delta.vec() / sinHalf;
    return angle * axis;
  }

  Eigen::Matrix3d applyRotationError(
    const Eigen::Matrix3d &current,
    const Eigen::Vector3d &error)
  {
    const double angle = error.norm();
    if (angle < 1e-12) {
      return current;
    }
    const Eigen::Vector3d axis = error / angle;
    const Eigen::Matrix3d delta = Eigen::AngleAxisd(angle, axis).toRotationMatrix();
    return delta * current;
  }

  Eigen::VectorXd applyJointStep(const Eigen::VectorXd step, const double alpha){
    std::cout << "[RaipalKin] applyJointStep step =" << step.transpose() << std::endl; //_debug
    std::cout << "[RaipalKin] applyJointStep alpha=" << alpha << std::endl; //_debug
    Eigen::VectorXd newGc(9);
    newGc <<
      gc_.head(3) + step.head(3) * alpha,
      gc_.segment(3,2), // this is to make joint-level distance calculation easier
      gc_.tail(4) + step.tail(4) * alpha
    ;

    // rk::cfbBackward(newGc); // to be safe

    return newGc;
  }

  bool isImprovement(const Eigen::VectorXd newGc, const Eigen::VectorXd &targetPos, const Eigen::MatrixXd &targetRot, Eigen::Vector3d &dPos, Eigen::Vector3d &dRot, const std::string &frameName, const IkSettings cfg, const bool posOnly = false){
    std::cout << "[RaipalKin] Improvement check" << std::endl; //_debug

    // get required dRot and dPos
    auto gcPrev = gc_;

    auto [cPos, cRot] = fk(frameName,newGc,false);
    Eigen::Vector3d cdPos = targetPos - cPos;
    Eigen::Vector3d cdRot = rotationErrorQuat(cRot,targetRot);
    std::cout << "[RaipalKin] candidate newGc=" << newGc.transpose() << std::endl; //_debug
    std::cout << "[RaipalKin] improvement candidate cdPos=" << cdPos.norm() << " cdRot=" << cdRot.norm() << std::endl; //_debug
    std::cout << "[RaipalKin] initial                dPos=" << dPos.norm() << "  dRot=" << dRot.norm() << std::endl; //_debug

    if(
      cdPos.norm()/cfg.posTolerance + (posOnly ? 0.0 : cdRot.norm()/cfg.rotTolerance) <= 
       dPos.norm()/cfg.posTolerance +  (posOnly ? 0.0 : dRot.norm()/cfg.rotTolerance)
    ){
      dPos << cdPos;
      dRot << cdRot;
      std::cout << "[RaipalKin] improvement accepted" << std::endl; //_debug
      return true;
    } else {
      set(gcPrev,false);
      std::cout << "[RaipalKin] improvement rejected" << std::endl; //_debug
      return false;
    }
  }

  void printJointLimits() const {
    if (jointLimits_.empty()) {
      std::cout << " (none)" << std::endl;
      return;
    }
    for (size_t i = 0; i < jointLimits_.size(); ++i) {
      std::cout << " [" << jointLimits_[i][0] << ", " << jointLimits_[i][1] << "]";
    }
    std::cout << std::endl;
  }

  Eigen::VectorXd gc_;
  IkSettings ikSettings_{};
  std::unique_ptr<raisim::World> world_;
  raisim::ArticulatedSystem *robot_{nullptr};
  size_t gcDim_{0};
  std::vector<raisim::Vec<2>> jointLimits_{};
  size_t totalIterationsAccum_{0};
  size_t totalSolveCount_{0};
  double minAcceptedAlphaGlobal_{std::numeric_limits<double>::infinity()};
  int maxRetriesObserved_{0};
};

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
    RaipalKin::IkResult ikResult;
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

