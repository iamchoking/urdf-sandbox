#pragma once

#include "raisim/World.hpp"
#include "raipal_solution/cfbSolution.hpp"

#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

#define RAIPAL_FORWARD(method)                             \
  template <typename... Args>                              \
  decltype(auto) method(Args&&... args) {                  \
    return robot_->method(std::forward<Args>(args)...);    \
  }                                                        \
                                                           \
  template <typename... Args>                              \
  decltype(auto) method(Args&&... args) const {            \
    return robot_->method(std::forward<Args>(args)...);    \
  }

#define RAIPAL_SET_UPDATE(method)                          \
  template <typename... Args>                              \
  void method(Args&&... args) {                            \
    robot_->method(std::forward<Args>(args)...);           \
    updateRaipal(true);                                    \
  }

class articulatedRaipal {

  public:
  // constructor from existing articulated system
  explicit articulatedRaipal(
    raisim::ArticulatedSystem* robot,
    std::vector<size_t> cfbIndices    = {3  },
    std::vector<double> cfbDirections = {1.0}
  ): 
    robot_(robot),
    rotorInertia_(robot_->getRotorInertia()),
    maxTorque_(robot_->getActuationUpperLimits().e()),
    minTorque_(robot_->getActuationLowerLimits().e()),
    maxVelocity_(robot_->getJointVelocityLimits())
  {
    dof_ = robot_->getDOF();
    for(size_t i = 0; i < cfbIndices.size(); i++){
      actuatorStates_.push_back({cfbIndices[i], cfbDirections[i]});
    }
    tauFF_ = robot_->getFeedForwardGeneralizedForce().e();
    tauCFB_.setZero(tauFF_.size());
  }

  // casting / conversion for upstream compatiblity
  raisim::ArticulatedSystem* get() { return robot_; }
  const raisim::ArticulatedSystem* get() const { return robot_; }

  operator raisim::ArticulatedSystem*() { return robot_; }
  operator const raisim::ArticulatedSystem*() const { return robot_; }

  // Raipal specific methods
  void updateRaipal(bool forceUpdate = false){
    if(updated_ && !forceUpdate) return;

    Eigen::VectorXd gc,gv,ga;
    robot_->getState(gc, gv);
    ga = robot_->getGeneralizedAcceleration().e();

    Eigen::VectorXd pTarget(robot_->getGeneralizedCoordinateDim());
    Eigen::VectorXd dTarget(robot_->getDOF());
    Eigen::VectorXd pGain(robot_->getDOF());
    Eigen::VectorXd dGain(robot_->getDOF());
    robot_->getPdTarget(pTarget, dTarget);
    robot_->getPdGains(pGain, dGain);

    tauFF_ = robot_->getFeedForwardGeneralizedForce().e() - tauCFB_;
    tauCFB_.setZero();

    for(auto& as : actuatorStates_){
      // CFB solution
      double a, da, dda;
      cfb::evalCfb(cfb::fromJoint, as.sign * gc[as.idx], a, da, dda);

      // ------- Actuator-side update -------
      // kinematic state
      as.pos = as.sign * a;
      as.vel = as.sign * da * gv[as.idx];
      as.acc = as.sign * (dda * gv[as.idx] * gv[as.idx] + da * ga[as.idx]);

      cfb::evalCfb(cfb::fromJoint, as.sign * pTarget[as.idx], as.pTarget);
      as.dTarget = as.sign * da * dTarget[as.idx];

      // TODO: this may change if the gain is set in terms of actuator gain
      as.pGain   = pGain[as.idx] / (da * da);
      as.dGain   = dGain[as.idx] / (da * da);

      // ------- Joint-side update -------
      maxVelocity_[as.idx]  = maxVelocity_[as.idx-1] / da;
      rotorInertia_[as.idx] = (rotorInertia_[as.idx-1] + cfb::actuatorInertia) * da * da;

      tauCFB_[as.idx] = -as.sign * (rotorInertia_[as.idx-1] + cfb::actuatorInertia) * da * dda * gv[as.idx] * gv[as.idx];
      maxTorque_[as.idx]    = maxTorque_[as.idx-1] * da;
      minTorque_[as.idx]    = minTorque_[as.idx-1] * da;
    }

    robot_->setJointVelocityLimits(maxVelocity_);
    robot_->setActuationLimits(maxTorque_ + tauCFB_, minTorque_ + tauCFB_);
    robot_->setRotorInertia(rotorInertia_);
    robot_->setGeneralizedForce(tauFF_ + tauCFB_);

    updated_ = true;
  }

  void resetUpdateFlag(bool updated = false) { updated_ = updated; }

  // actuator-side methods
  void getActuatorState(Eigen::VectorXd &genco, Eigen::VectorXd &genvel){
    robot_->getState(genco, genvel);
    for(auto& as : actuatorStates_){
      genco[as.idx]  = as.pos;
      genvel[as.idx] = as.vel;
    }
  }

  void getActuatorPdTarget(Eigen::VectorXd &pTarget, Eigen::VectorXd &dTarget){
    robot_->getPdTarget(pTarget, dTarget);
    for(auto& as : actuatorStates_){
      pTarget[as.idx] = as.pTarget;
      dTarget[as.idx] = as.dTarget;
    }
  }

  void getActuatorPdGains(Eigen::VectorXd &pGain, Eigen::VectorXd &dGain){
    robot_->getPdGains(pGain, dGain);
    for(auto& as : actuatorStates_){
      pGain[as.idx] = as.pGain;
      dGain[as.idx] = as.dGain;
    }
  }

  // getters that need special care
  /**
   * @return the upper joint torque/force limit*/
  [[nodiscard]] const raisim::VecDyn getActuationUpperLimits() const {
    raisim::VecDyn tauUpper;
    tauUpper = maxTorque_;
    return tauUpper;
  }

  /**
   * @return the lower joint torque/force limit*/
  [[nodiscard]] const raisim::VecDyn getActuationLowerLimits() const {
    raisim::VecDyn tauLower;
    tauLower = minTorque_;
    return tauLower;
  }

  [[nodiscard]] const raisim::VecDyn getFeedForwardGeneralizedForce() const {
    raisim::VecDyn tauFF;
    tauFF = tauFF_;
    return tauFF;
  }

  [[nodiscard]] raisim::VecDyn getGeneralizedForce() const {
    Eigen::VectorXd gf = robot_->getGeneralizedForce().e() - tauCFB_;
    raisim::VecDyn genForce;
    genForce = gf;
    return genForce;
  }

  // BOOKMARK: start coding setters
  // setters that need special care

  // RAIPAL_SET_UPDATE(setRotorInertia)
  void setRotorInertia(const raisim::VecDyn &rotorInertia) {
    rotorInertia_ = rotorInertia;
    updateRaipal(true);
  }

  // RAIPAL_SET_UPDATE(setActuationLimits)
  void setActuationLimits(const Eigen::VectorXd &upper, const Eigen::VectorXd &lower) {
    maxTorque_ = upper;
    minTorque_ = lower;
    updateRaipal(true);
  }

  // RAIPAL_SET_UPDATE(setJointVelocityLimits)
  void setJointVelocityLimits(const Eigen::VectorXd &velLimits) {
    maxVelocity_ = velLimits;
    updateRaipal(true);
  }

  // RAIPAL_SET_UPDATE(setGeneralizedForce)      // tauFF_
  void setGeneralizedForce(const raisim::VecDyn &tau) { 
    tauFF_ = tau.e();
    tauCFB_.setZero();
    robot_->setGeneralizedForce(tauFF_);
    updateRaipal(true);
  }

  void setGeneralizedForce(const Eigen::VectorXd &tau) { 
    tauFF_ = tau;
    tauCFB_.setZero();
    robot_->setGeneralizedForce(tauFF_);
    updateRaipal(true);
  }

  // setters that require updates
  RAIPAL_SET_UPDATE(setState)                 // gc, gv
  RAIPAL_SET_UPDATE(setGeneralizedCoordinate) // gc 
  RAIPAL_SET_UPDATE(setGeneralizedVelocity)   // gv
  RAIPAL_SET_UPDATE(setPdTarget)              // pTarget, dTarget
  RAIPAL_SET_UPDATE(setPTarget)               // pTarget
  RAIPAL_SET_UPDATE(setDTarget)               // dTarget
  RAIPAL_SET_UPDATE(setPdGains)               // pGain, dGain
  RAIPAL_SET_UPDATE(setPGains)                // pGain
  RAIPAL_SET_UPDATE(setDGains)                // dGain

  // forwarded to ArticluatedSystem (no need for wrapper intervention)
  RAIPAL_FORWARD(isSecondOrderOrHigher)
  RAIPAL_FORWARD(getGeneralizedCoordinate)
  RAIPAL_FORWARD(getGeneralizedVelocity)
  RAIPAL_FORWARD(getGeneralizedAcceleration)
  RAIPAL_FORWARD(getBaseOrientation)
  RAIPAL_FORWARD(getBasePosition)
  RAIPAL_FORWARD(updateKinematics)
  RAIPAL_FORWARD(getState)
  RAIPAL_FORWARD(getMassMatrix)
  RAIPAL_FORWARD(getNonlinearities)
  RAIPAL_FORWARD(getInverseMassMatrix)
  RAIPAL_FORWARD(getCompositeCOM)
  RAIPAL_FORWARD(getCOM)
  RAIPAL_FORWARD(getCompositeInertia)
  RAIPAL_FORWARD(getCompositeMass)
  RAIPAL_FORWARD(getLinearMomentum)
  RAIPAL_FORWARD(getGeneralizedMomentum)
  RAIPAL_FORWARD(getEnergy)
  RAIPAL_FORWARD(getKineticEnergy)
  RAIPAL_FORWARD(getPotentialEnergy)
  RAIPAL_FORWARD(getAngularMomentum)
  RAIPAL_FORWARD(printOutBodyNamesInOrder)
  RAIPAL_FORWARD(printOutMovableJointNamesInOrder)
  RAIPAL_FORWARD(printOutFrameNamesInOrder)
  RAIPAL_FORWARD(getMovableJointNames)
  RAIPAL_FORWARD(getGeneralizedVelocityIndex)
  RAIPAL_FORWARD(getPosition)
  RAIPAL_FORWARD(getFrameByName)
  RAIPAL_FORWARD(getFrameByLinkName)
  RAIPAL_FORWARD(getFrameIdxByLinkName)
  RAIPAL_FORWARD(getFrameByIdx)
  RAIPAL_FORWARD(getFrameIdxByName)
  RAIPAL_FORWARD(getFrames)
  RAIPAL_FORWARD(getFramePosition)
  RAIPAL_FORWARD(getPositionInFrame)
  RAIPAL_FORWARD(getFrameOrientation)
  RAIPAL_FORWARD(getFrameVelocity)
  RAIPAL_FORWARD(getFrameAngularVelocity)
  RAIPAL_FORWARD(getFrameAcceleration)
  RAIPAL_FORWARD(getPositionInBodyCoordinate)
  RAIPAL_FORWARD(getOrientation)
  RAIPAL_FORWARD(getVelocity)
  RAIPAL_FORWARD(getAngularVelocity)
  RAIPAL_FORWARD(getSparseJacobian)
  RAIPAL_FORWARD(getSparseRotationalJacobian)
  RAIPAL_FORWARD(getTimeDerivativeOfSparseJacobian)
  RAIPAL_FORWARD(getTimeDerivativeOfSparseRotationalJacobian)
  RAIPAL_FORWARD(getDenseJacobian)
  RAIPAL_FORWARD(getDenseRotationalJacobian)
  RAIPAL_FORWARD(getDenseFrameJacobian)
  RAIPAL_FORWARD(getDenseFrameRotationalJacobian)
  RAIPAL_FORWARD(getBodyIdx)
  RAIPAL_FORWARD(getDOF)
  RAIPAL_FORWARD(getGeneralizedVelocityDim)
  RAIPAL_FORWARD(getGeneralizedCoordinateDim)
  RAIPAL_FORWARD(getBodyPose)
  RAIPAL_FORWARD(getBodyPosition)
  RAIPAL_FORWARD(getBodyOrientation)
  RAIPAL_FORWARD(getJointPos_P)
  RAIPAL_FORWARD(getJointOrientation_P)
  RAIPAL_FORWARD(getJointAxis_P)
  RAIPAL_FORWARD(getJointAxis)
  RAIPAL_FORWARD(getMass)
  RAIPAL_FORWARD(getInertia)
  RAIPAL_FORWARD(getBodyCOM_B)
  RAIPAL_FORWARD(getBodyCOM_W)
  RAIPAL_FORWARD(getCollisionBodies)
  RAIPAL_FORWARD(getCollisionBody)
  RAIPAL_FORWARD(updateMassInfo)
  RAIPAL_FORWARD(setMass)
  RAIPAL_FORWARD(getTotalMass)
  RAIPAL_FORWARD(setExternalForce)
  RAIPAL_FORWARD(setExternalTorque)
  RAIPAL_FORWARD(setExternalTorqueInBodyFrame)
  RAIPAL_FORWARD(getContactPointVel)
  RAIPAL_FORWARD(setControlMode)
  RAIPAL_FORWARD(getControlMode)
  RAIPAL_FORWARD(getPdTarget)
  RAIPAL_FORWARD(getPdGains)
  RAIPAL_FORWARD(setJointDamping)
  RAIPAL_FORWARD(computeSparseInverse)
  RAIPAL_FORWARD(massMatrixVecMul)
  RAIPAL_FORWARD(ignoreCollisionBetween)
  RAIPAL_FORWARD(getOptions)
  RAIPAL_FORWARD(getBodyNames)
  RAIPAL_FORWARD(getVisOb)
  RAIPAL_FORWARD(getVisColOb)
  RAIPAL_FORWARD(getVisObPose)
  RAIPAL_FORWARD(getVisColObPose)
  RAIPAL_FORWARD(getResourceDir)
  RAIPAL_FORWARD(getRobotDescriptionfFileName)
  RAIPAL_FORWARD(getRobotDescriptionfTopDirName)
  RAIPAL_FORWARD(getRobotDescriptionFullPath)
  RAIPAL_FORWARD(getRobotDescription)
  RAIPAL_FORWARD(exportRobotDescriptionToURDF)
  RAIPAL_FORWARD(setBasePos_e)
  RAIPAL_FORWARD(setBaseOrientation_e)
  RAIPAL_FORWARD(setBasePos)
  RAIPAL_FORWARD(setBaseOrientation)
  RAIPAL_FORWARD(setBaseVelocity)
  RAIPAL_FORWARD(setBaseAngularVelocity)
  RAIPAL_FORWARD(setCollisionObjectShapeParameters)
  RAIPAL_FORWARD(setCollisionObjectPositionOffset)
  RAIPAL_FORWARD(setCollisionObjectOrientationOffset)
  RAIPAL_FORWARD(getRotorInertia)
  RAIPAL_FORWARD(getJointType)
  RAIPAL_FORWARD(getNumberOfJoints)
  RAIPAL_FORWARD(getJoint)
  RAIPAL_FORWARD(getLink)
  RAIPAL_FORWARD(getMappingFromBodyIndexToGeneralizedVelocityIndex)
  RAIPAL_FORWARD(getMappingFromBodyIndexToGeneralizedCoordinateIndex)
  RAIPAL_FORWARD(getObjectType)
  RAIPAL_FORWARD(getBodyType)
  RAIPAL_FORWARD(setIntegrationScheme)
  RAIPAL_FORWARD(getJointLimitViolations)
  RAIPAL_FORWARD(setJointLimits)
  RAIPAL_FORWARD(getJointLimits)
  RAIPAL_FORWARD(getJointVelocityLimits)
  RAIPAL_FORWARD(clearExternalForcesAndTorques)
  RAIPAL_FORWARD(addSpring)
  RAIPAL_FORWARD(getSprings)
  RAIPAL_FORWARD(getParentVector)
  RAIPAL_FORWARD(getSensorSet)
  RAIPAL_FORWARD(getSensorSets)
  RAIPAL_FORWARD(addConstraints)
  RAIPAL_FORWARD(initializeConstraints)
  RAIPAL_FORWARD(setComputeInverseDynamics)
  RAIPAL_FORWARD(getComputeInverseDynamics)
  RAIPAL_FORWARD(getForceAtJointInWorldFrame)
  RAIPAL_FORWARD(getTorqueAtJointInWorldFrame)
  RAIPAL_FORWARD(getAllowedNumberOfInternalContactsBetweenTwoBodies)
  RAIPAL_FORWARD(setAllowedNumberOfInternalContactsBetweenTwoBodies)
  RAIPAL_FORWARD(articulatedBodyAlgorithm)
  RAIPAL_FORWARD(getMinvJT)
  RAIPAL_FORWARD(getj_MinvJT_T1D)
  RAIPAL_FORWARD(getUdot)
  RAIPAL_FORWARD(getFullDelassusAndTauStar)
  RAIPAL_FORWARD(appendJointLimits)

  void setGeneralizedCoordinate(std::initializer_list<double> jointState) {
    robot_->setGeneralizedCoordinate(jointState);
  }

  void setGeneralizedVelocity(std::initializer_list<double> jointVelocity) {
    robot_->setGeneralizedVelocity(jointVelocity);
  }

  void setGeneralizedForce(std::initializer_list<double> tau) {
    robot_->setGeneralizedForce(tau);
  }

  static void convertSparseJacobianToDense(
      const raisim::SparseJacobian& sparseJaco,
      Eigen::MatrixXd& denseJaco) {
    raisim::ArticulatedSystem::convertSparseJacobianToDense(sparseJaco, denseJaco);
  }

private:

  struct actuatorState{
    size_t idx;
    double sign;

    // kinematic state
    double pos;
    double vel;
    double acc;

    // control variables
    double pTarget;
    double dTarget;
    double pGain;
    double dGain;

    // dynamic state
    double gfTerm = 0.0;
  };

  raisim::ArticulatedSystem* robot_ = nullptr;
  size_t dof_;

  std::vector<actuatorState> actuatorStates_;

  bool updated_ = false;
  
  raisim::VecDyn rotorInertia_, maxVelocity_;
  Eigen::VectorXd maxTorque_, minTorque_;

  Eigen::VectorXd tauFF_;   // feedforward term for rest of robot
  Eigen::VectorXd tauCFB_;  // feedforward term for CFB correction
};

#undef RAIPAL_FORWARD
#undef RAIPAL_SET_UPDATE
