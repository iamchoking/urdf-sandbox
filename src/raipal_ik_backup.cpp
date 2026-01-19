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

class RaipalIkSolver {
public:
  explicit RaipalIkSolver(const std::string &urdfPath)
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
    int maxSteps = 100;
    int maxIters = 10;
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
  };

  struct IkMetrics{
    bool converged = false;
    std::string warning = "";
    size_t step   = 0;
    size_t iter   = 0;
    size_t retry  = 0;

    size_t totalIters = 0;
    double minAcceptedAlpha = std::numeric_limits<double>::infinity();
    int maxRetries = 0;
  };

  std::string ikStatus(const IkMetrics &metrics) const {
    if(metrics.step == 0 && metrics.iter == 0){
      return "[IK]";
    }
    
    std::ostringstream oss;
    oss << "[IK s" << metrics.step << " i" << metrics.iter << " r" << metrics.retry << "]";
    return oss.str();
  }

  void resetIkSettings() { 
    IkSettings settings{};
    ikSettings_ = settings; 
  }

  IkMetrics ik(
    const std::string &frameName,
    const Eigen::Vector3d &targetPos,
    const Eigen::Matrix3d &targetRot,
    Eigen::VectorXd &sol
  ){
    IkSettings cfg{};
    IkMetrics  met{};
    std::cout << ikStatus(met) << " ik() start" << std::endl; //_debug
    std::cout << "gc_: " << gc_.transpose() << std::endl; //_debug

    // initialize solution scheme
    const size_t frameIdx = robot_->getFrameIdxByName(frameName);
    std::cout << ikStatus(met) << " frameIdx: " << frameIdx << std::endl; //_debug
    if(sol != Eigen::VectorXd::Zero(9) && sol != Eigen::VectorXd::Zero(7)){set(sol);}
    sol = gc_; // this ensures that sol is 9d

    for(met.step = 0; met.step < cfg.maxSteps; ++met.step){ // [step] loop
      // Calculate distance from current to desired position
      auto [tPos, tRot] = fk(frameName);
      std::cout << ikStatus(met) << " initial FK computed" << std::endl; //_debug

      Eigen::Vector3d dPos = targetPos - tPos;
      Eigen::Vector3d dRot = rotationErrorQuat(tRot,targetRot);

      // [0] (preprocess) scale down the distance and rotation if it's too big

      Eigen::Vector3d nTargetPos;
      Eigen::Matrix3d nTargetRot;
      bool multistep = false;
      double factor = 1.0;

      if(dPos.norm() > cfg.maxPosDistance){
        std::cout << "Position distance over limit: " << dPos.norm() << " -> " << cfg.maxPosDistance << std::endl;
        factor = cfg.maxPosDistance / dPos.norm();
        multistep = true;
      }
      if(dRot.norm() > cfg.maxRotDistance){
        std::cout << "Angular distance over limit: " << dRot.norm() << " -> " << cfg.maxRotDistance << std::endl;
        factor = std::min(factor, cfg.maxRotDistance / dRot.norm());
        multistep = true;
      }

      if(multistep){
        std::cout << ikStatus(met) << " multistep triggered with factor: " << factor << std::endl; //_debug
        factor = std::min(0.5, factor); // if we are doing multistep, limit to half-step at most
        dPos = dPos * factor;
        dRot = dRot * factor;
        nTargetPos = tPos + dPos;
        nTargetRot = applyRotationError(tRot, dRot);
      } else{
        nTargetPos << targetPos;
        nTargetRot << targetRot;
      }

      // [1] iterative solver
      // objective: obtain nTargetPos and nTargetRot
      Eigen::MatrixXd J_pos9(3,9);
      Eigen::MatrixXd J_rot9(3,9);
      Eigen::MatrixXd J_7(6,7);

      for (met.iter = 0; met.iter < cfg.maxIters; ++met.iter) { // [iter] loop
        // [1-a] get jacobian
        J_pos9.setZero();
        J_rot9.setZero();
        J_7.setZero();

        robot_->getDenseFrameJacobian(frameIdx, J_pos9);
        robot_->getDenseFrameRotationalJacobian(frameIdx, J_rot9);

        // work with the zero columns removed to prevent ill-conditioned jacobians
        J_7.topRows(3)    << J_pos9.leftCols(3), J_pos9.rightCols(4);
        J_7.bottomRows(3) << J_rot9.leftCols(3), J_rot9.rightCols(4);

        std::cout << ikStatus(met) << " Jacobians retrieved" << std::endl; //_debug

        // [1-b] solve pseudoinverse
        Eigen::VectorXd taskError(6);
        taskError << dPos, dRot;

        const Eigen::MatrixXd JT = J_7.transpose();
        Eigen::MatrixXd JJt = J_7 * JT;
        JJt.diagonal().array() += cfg.damping * cfg.damping;

        // direct solution
        // Eigen::VectorXd jointStep = JT * JJt.ldlt().solve(taskError);        

        // full solution (includes null space projector)
        const Eigen::MatrixXd Jplus = (JT * JJt.inverse());
        const Eigen::MatrixXd Nproj = Eigen::MatrixXd::Identity(7,7) - Jplus * J_7; // null-space projector
        std::cout << ikStatus(met) << " Pseudoinverse components ready" << std::endl; //_debug

        Eigen::VectorXd jointStep = Jplus * taskError;
        std::cout << ikStatus(met) << " Joint Step Calculated (norm: " << jointStep.norm() << ")" << std::endl; //_debug

        // [1-c] apply joint step with varying alpha values
        for(met.retry = 0; met.retry <= cfg.maxRetries; met.retry++){
          double alpha = cfg.maxAlpha;
          Eigen::VectorXd cGc;
          while(alpha > cfg.minAlpha){ // [alpha] loop
            cGc = applyJointStep(jointStep, alpha);

            // feasibility check
            bool inLimits, isImprovement;
            inLimits = (cGc - sol).cwiseAbs().maxCoeff() > cfg.maxJointDistance;
            isImprovement = checkImprovement(cGc,nTargetPos,nTargetRot,dPos,dRot,frameName,cfg);

            if(!inLimits){
              std::cout << ikStatus(met) << " candidate exceeds joint limits" << std::endl; //_debug
            }
            else if(!isImprovement){
              std::cout << ikStatus(met) << " candidate does not improve error" << std::endl; //_debug
            }
            else {
              std::cout << ikStatus(met) << " candidate accepted with alpha: " << alpha << std::endl; //_debug
              if(alpha < met.minAcceptedAlpha){met.minAcceptedAlpha = alpha;}
              break; // exit alpha loop
            }

            alpha = alpha * 0.5; // keep decreaseing the alpha value until feasibility is met.
          } // [alpha] loop

          if(alpha < cfg.minAlpha){
            std::cout << ikStatus(met) << " Alpha value reached minimum!" << std::endl; // _debug

            Eigen::VectorXd factor(7);
            factor.setRandom();
            Eigen::VectorXd deviation = Nproj * factor;
            jointStep += deviation * (jointStep.norm() / factor.norm() * 0.1);
            alpha = cfg.maxAlpha * 2; // reset alpha
            std::cout << "[RETRYING] Changed joint step to " << jointStep.transpose() << std::endl; //_debug
          } else {
            std::cout << ikStatus(met) << " step applied with final alpha: " << alpha << std::endl; //_debug
            std::cout << ikStatus(met) << "iterated gc: " << std::setprecision(8) << gc_.transpose() << std::endl;
            break; // iteration successful. no retries needed
          }
        } // [retry] loop

        if(met.retry >= cfg.maxRetries){
          std::cout << ikStatus(met) << " retries exceeded, relaxing tolerance." << std::endl; //_debug
          cfg.posTolerance *= 2.0;
          cfg.rotTolerance *= 2.0;
          // break; // restart iteration
        }

        if(dPos.norm() < cfg.posTolerance && dRot.norm() < cfg.rotTolerance ){
          met.converged = true;
          std::cout << ikStatus(met) << " converged within loop [" << met.iter << "]" << std::endl; //_debug
          break;
        }
      } // [iter] loop

      sol = gc_;
      if(!multistep){break;}
    } // [step] loop

    if(met.step >= cfg.maxSteps){
      std::cout << ikStatus(met) << "[FAIL] max-steps reached" << std::endl; //_debug
    } else {
      std::cout << ikStatus(met) << "[SUCC] solve finished, (iterations:" << met.iter << ", retries:" << met.retry << ", multistep=" << met.step << ")" << std::endl; //_debug
    }

    rk::cfbBackward(sol);
    return met; 
  } // [ik] function

private:
  Eigen::Vector3d rotationError(const Eigen::Matrix3d &current, const Eigen::Matrix3d &target)
  {
    // Express the orientation error in world coordinates so it matches the Jacobian frame.
    const Eigen::Matrix3d delta = target * current.transpose();
    Eigen::Vector3d error;
    error <<
      delta(2, 1) - delta(1, 2),
      delta(0, 2) - delta(2, 0),
      delta(1, 0) - delta(0, 1);
    return 0.5 * error;
  }

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
    std::cout << "[applyJointStep]" << " step  =" << step.transpose() << std::endl; //_debug
    std::cout << "[applyJointStep]" << " alpha =" << alpha << std::endl; //_debug
    Eigen::VectorXd newGc(9);
    newGc <<
      gc_.head(3) + step.head(3) * alpha,
      gc_.segment(3,2), // this is to make joint-level distance calculation easier
      gc_.tail(4) + step.tail(4) * alpha
    ;

    // rk::cfbBackward(newGc); // to be safe

    return newGc;
  }

  bool checkImprovement(const Eigen::VectorXd newGc, const Eigen::VectorXd &targetPos, const Eigen::MatrixXd &targetRot, Eigen::Vector3d &dPos, Eigen::Vector3d &dRot, const std::string &frameName, const IkSettings cfg){
    std::cout << "Improvement check" << std::endl; //_debug

    // get required dRot and dPos
    auto gcPrev = gc_;

    auto [cPos, cRot] = fk(frameName,newGc,false);
    Eigen::Vector3d cdPos = targetPos - cPos;
    Eigen::Vector3d cdRot = rotationErrorQuat(cRot,targetRot);
    std::cout << "[checkImprovement]" << " candidate newGc=" << newGc.transpose() << std::endl; //_debug
    std::cout << "[checkImprovement]" << " improvement candidate cdPos=" << cdPos.norm() << " cdRot=" << cdRot.norm() << std::endl; //_debug
    std::cout << "[checkImprovement]" << " initial                dPos=" << dPos.norm() << "  dRot=" << dRot.norm() << std::endl; //_debug

    if(
      cdPos.norm()/cfg.posTolerance + cdRot.norm()/cfg.rotTolerance <= 
       dPos.norm()/cfg.posTolerance +  dRot.norm()/cfg.rotTolerance
    ){
      dPos << cdPos;
      dRot << cdRot;
      std::cout << "[checkImprovement]" << " improvement accepted" << std::endl; //_debug
      return true;
    } else {
      set(gcPrev);
      std::cout << "[checkImprovement]" << " improvement rejected" << std::endl; //_debug
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
  std::unique_ptr<raisim::World> world_;
  raisim::ArticulatedSystem *robot_{nullptr};
  size_t gcDim_{0};
  std::vector<raisim::Vec<2>> jointLimits_{};

  IkSettings ikSettings_{};
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

  RaipalIkSolver solver(raipalRightUrdf);

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

    const auto ikStart = std::chrono::high_resolution_clock::now();
    RaipalIkSolver::IkMetrics metrics;
    try {
      metrics = solver.ik(kTipFrameName,targetPos,targetRot,ikSol);
    } catch (const std::exception &e) {
      abortReason = e.what();
      aborted = true;
      ++failureCount;
      std::cout << "[main] solver.ik threw: " << abortReason << std::endl;
      break;
    }
    std::cout << "[main] solver.ik done, success=" << metrics.converged << std::endl; //_debug

    const auto ikEnd = std::chrono::high_resolution_clock::now();

    const double solveNs = std::chrono::duration_cast<std::chrono::nanoseconds>(ikEnd - ikStart).count();
    maxSolveNs = std::max(maxSolveNs, solveNs);

    visualRobot->setGeneralizedCoordinate(ikSol);

    const auto fkResult = solver.fk(kTipFrameName, ikSol);
    const double posError = (fkResult.first - targetPos).norm();
    const double oriError = orientationAngleError(targetRot, fkResult.second);

    accumulatedTimeNs += solveNs;
    accumulatedPosError += posError;
    accumulatedOriError += oriError;
    maxPosError = std::max(maxPosError, posError);
    maxOriError = std::max(maxOriError, oriError);

    if (metrics.converged) {
      ++successCount;

      if (step % 100 == 0) {
        std::cout << std::fixed << std::setprecision(3)
                  << "[IK] step " << step
                  << " | pos err (mm): " << posError * 1e3
                  << " | ori err (deg): " << oriError * 180.0 / M_PI
                  << " | solve (us): " << solveNs / 1e3
                  << std::endl;
      }
    } else {
      ++failureCount;
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

