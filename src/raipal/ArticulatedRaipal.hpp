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
    return articulated_->method(std::forward<Args>(args)...);    \
  }                                                        \
                                                           \
  template <typename... Args>                              \
  decltype(auto) method(Args&&... args) const {            \
    return articulated_->method(std::forward<Args>(args)...);    \
  }

#define RAIPAL_UPDATE_GET(method)                          \
  template <typename... Args>                              \
  decltype(auto) method(Args&&... args) {                  \
    updateRaipal();                                        \
    return articulated_->method(std::forward<Args>(args)...);    \
  }                                                        \
                                                           \
  template <typename... Args>                              \
  decltype(auto) method(Args&&... args) const {            \
    const_cast<ArticulatedRaipal*>(this)->updateRaipal();  \
    return articulated_->method(std::forward<Args>(args)...);    \
  }

#define RAIPAL_SET(method)                                 \
  template <typename... Args>                              \
  void method(Args&&... args) {                            \
    articulated_->method(std::forward<Args>(args)...);           \
    resetUpdateFlag();                                     \
  }

class ArticulatedRaipal {
  public:

  enum class CfbTargetMode {
    FROM_JOINT,    // set target and gains in terms of joint space (actuator-side pd gains change every update)
    FROM_ACTUATOR  // set target and gains in terms of actuator space (joint-side pd gains change every update)
  };

  // constructor from existing articulated system
  explicit ArticulatedRaipal(
    raisim::ArticulatedSystem* robot,
    std::vector<size_t> cfbIndices    = {3  },
    std::vector<int> cfbDirections = {1}
  ): 
    articulated_(robot),
    rotorInertia_(articulated_->getRotorInertia()),
    maxTorque_(articulated_->getActuationUpperLimits().e()),
    minTorque_(articulated_->getActuationLowerLimits().e()),
    maxVelocity_(articulated_->getJointVelocityLimits())
  {
    dof_ = articulated_->getDOF();
    for(size_t i = 0; i < cfbIndices.size(); i++){
      actuatorStates_.push_back({cfbIndices[i], cfbDirections[i]});
    }
    tauFF_ = articulated_->getFeedForwardGeneralizedForce().e();
    tauCFB_.setZero(tauFF_.size());
  }

  // casting / conversion for upstream compatiblity
  raisim::ArticulatedSystem* get() { return articulated_; }
  const raisim::ArticulatedSystem* get() const { return articulated_; }

  operator raisim::ArticulatedSystem*() { return articulated_; }
  operator const raisim::ArticulatedSystem*() const { return articulated_; }

  ArticulatedRaipal* operator->() { return this; }
  const ArticulatedRaipal* operator->() const { return this; }

  // Raipal-specific methods
  void updateRaipal(bool forceUpdate = false){
    if(updated_ && !forceUpdate) return;

    Eigen::VectorXd gc,gv,ga;
    articulated_->getState(gc, gv);
    ga = articulated_->getGeneralizedAcceleration().e();

    Eigen::VectorXd pTarget(articulated_->getGeneralizedCoordinateDim());
    Eigen::VectorXd dTarget(articulated_->getDOF());
    Eigen::VectorXd pGain(articulated_->getDOF());
    Eigen::VectorXd dGain(articulated_->getDOF());
    articulated_->getPdTarget(pTarget, dTarget);
    articulated_->getPdGains(pGain, dGain);

    tauFF_ = articulated_->getFeedForwardGeneralizedForce().e() - tauCFB_;
    tauCFB_.setZero();

    for(auto& as : actuatorStates_){
      // CFB solution
      double a, da, dda;
      cfb::evalCfb(cfb::fromJoint, gc[as.idx], as.sign, a, da, dda);

      // ------- Actuator-side update -------
      // kinematic state
      as.pos = a;
      as.vel = da * gv[as.idx];
      as.acc = (dda * gv[as.idx] * gv[as.idx]) + da * ga[as.idx];
      
      if(cfbMode_ == CfbTargetMode::FROM_JOINT){
        // control variables (joint -> actuator)
        cfb::evalCfb(cfb::fromJoint, pTarget[as.idx], as.sign, as.pTarget);
        as.dTarget = da * dTarget[as.idx];

      }

      else if(cfbMode_ == CfbTargetMode::FROM_ACTUATOR){
        // control variables (joint -> actuator)

        cfb::evalCfb(cfb::fromActuator, as.pTarget, as.sign, pTarget(as.idx));
        dTarget(as.idx) = as.dTarget / da;
 
      }

    
      // gain setting for torque-equivalence
      double jointPerror = pTarget[as.idx] - gc[as.idx];
      double actPerror = as.pTarget - as.pos;
      if(std::abs(jointPerror) > 1e-6 && std::abs(actPerror) > 1e-6){
        // joint-error * joint-gain = actuator-error * actuator-gain * da

        // => joint-gain = actuator-error * actuator-gain * da / joint-error
        if(cfbMode_ == CfbTargetMode::FROM_ACTUATOR){pGain[as.idx] = as.pGain * da * actPerror / jointPerror;}

        // => actuator-gain = joint-gain * joint-error *  / (actuator-error * da)
        else{as.pGain = pGain[as.idx] * jointPerror / (actPerror * da);}
      }
      else{
        // if the error is too small, use naive gain mapping
        if(cfbMode_ == CfbTargetMode::FROM_ACTUATOR){pGain[as.idx] = as.pGain * (da * da);}
        else{as.pGain   = pGain[as.idx] / (da * da);}
      }

      double jointDerror = dTarget[as.idx] - gv[as.idx];
      double actDerror = as.dTarget - as.vel;
      if(std::abs(jointDerror) > 1e-6 && std::abs(actDerror) > 1e-6){
        // joint-error * joint-gain = actuator-error * actuator-gain * da
        // => joint-gain = actuator-error * actuator-gain * da / joint-error
        if(cfbMode_ == CfbTargetMode::FROM_ACTUATOR){dGain[as.idx] = as.dGain * da * actDerror / jointDerror;}

        // => actuator-gain = joint-gain * joint-error *  / (actuator-error * da)
        else{as.dGain = dGain[as.idx] * jointDerror / (actDerror * da);}        
      }
      else{
        // if the error is too small, use naive gain mapping
        if(cfbMode_ == CfbTargetMode::FROM_ACTUATOR){dGain[as.idx] = as.dGain * (da * da);}
        else{as.dGain   = dGain[as.idx] / (da * da);}
      }

      // [IMPORTANT] CFB-induced residual torque
      tauCFB_[as.idx] = - (rotorInertia_[as.idx-1] + cfb::actuatorInertia) * da * dda * gv[as.idx] * gv[as.idx];

      // ------- Joint-side update -------
      maxVelocity_[as.idx]  = maxVelocity_[as.idx-1] / da;
      rotorInertia_[as.idx] = (rotorInertia_[as.idx-1] + cfb::actuatorInertia) * da * da;

      maxTorque_[as.idx]    = maxTorque_[as.idx-1] * da;
      minTorque_[as.idx]    = minTorque_[as.idx-1] * da;
    }

    if(cfbMode_ == CfbTargetMode::FROM_ACTUATOR){
      articulated_->setPdTarget(pTarget, dTarget);
      articulated_->setPdGains(pGain,dGain);
    }

    articulated_->setJointVelocityLimits(maxVelocity_);
    articulated_->setActuationLimits(maxTorque_ + tauCFB_, minTorque_ + tauCFB_);
    articulated_->setRotorInertia(rotorInertia_);
    articulated_->setGeneralizedForce(tauFF_ + tauCFB_);

    updated_ = true;
  }

  void resetUpdateFlag(bool updated = false) { updated_ = updated; }

  // mode control
  void setCfbTargetMode(CfbTargetMode mode) { cfbMode_ = mode; }
  void setCfbTargetFromJoint() { cfbMode_ = CfbTargetMode::FROM_JOINT; }
  void setCfbTargetFromActuator() { cfbMode_ = CfbTargetMode::FROM_ACTUATOR; }

  // actuator-side methods
  void getActuatorState(Eigen::VectorXd &genco, Eigen::VectorXd &genvel){
    updateRaipal();
    articulated_->getState(genco, genvel);
    for(auto& as : actuatorStates_){
      genco[as.idx]  = as.pos;
      genvel[as.idx] = as.vel;
    }
  }

  void getActuatorPdTarget(Eigen::VectorXd &pTarget, Eigen::VectorXd &dTarget){
    updateRaipal();
    articulated_->getPdTarget(pTarget, dTarget);
    for(auto& as : actuatorStates_){
      pTarget[as.idx] = as.pTarget;
      dTarget[as.idx] = as.dTarget;
    }
  }

  void getActuatorPdGains(Eigen::VectorXd &pGain, Eigen::VectorXd &dGain){
    updateRaipal();
    articulated_->getPdGains(pGain, dGain);
    for(auto& as : actuatorStates_){
      pGain[as.idx] = as.pGain;
      dGain[as.idx] = as.dGain;
    }
  }

  void setActuatorPdTarget(const Eigen::VectorXd &pTarget, const Eigen::VectorXd &dTarget){\
    RSFATAL_IF(cfbMode_ != CfbTargetMode::FROM_ACTUATOR, "setActuatorPdTarget can only be used when CFB target mode is set to FROM_ACTUATOR");
    articulated_->setPdTarget(pTarget, dTarget);
    for(auto& as/*  */ : actuatorStates_){
      as.pTarget = pTarget[as.idx];
      as.dTarget = dTarget[as.idx];
    }
    resetUpdateFlag();
  }

  void setActuatorPdGains(const Eigen::VectorXd &pGain, const Eigen::VectorXd &dGain){
    RSFATAL_IF(cfbMode_ != CfbTargetMode::FROM_ACTUATOR, "setActuatorPdGains can only be used when CFB target mode is set to FROM_ACTUATOR");
    articulated_->setPdGains(pGain, dGain);
    for(auto& as : actuatorStates_){
      as.pGain = pGain[as.idx];
      as.dGain = dGain[as.idx];
    }
    resetUpdateFlag();
  }

  // TODOS
  // getActuatorGeneralizedForce()
  // getActuatorFeedForwardGeneralizedForce()

  // cfbMode_-aware getters and setters for control variables
  inline void setCurrentPdTarget(const Eigen::VectorXd &pTarget, const Eigen::VectorXd &dTarget){
    if(cfbMode_ == CfbTargetMode::FROM_JOINT){setPdTarget(pTarget, dTarget);}
    else if(cfbMode_ == CfbTargetMode::FROM_ACTUATOR){setActuatorPdTarget(pTarget, dTarget);}
  }

  inline void setCurrentPdGains(const Eigen::VectorXd &pGain, const Eigen::VectorXd &dGain){
    if(cfbMode_ == CfbTargetMode::FROM_JOINT){setPdGains(pGain, dGain);}
    else if(cfbMode_ == CfbTargetMode::FROM_ACTUATOR){setActuatorPdGains(pGain, dGain);}
  }

  inline void getCurrentPdTarget(Eigen::VectorXd &pTarget, Eigen::VectorXd &dTarget){
    if(cfbMode_ == CfbTargetMode::FROM_JOINT){getPdTarget(pTarget, dTarget);}
    else if(cfbMode_ == CfbTargetMode::FROM_ACTUATOR){getActuatorPdTarget(pTarget, dTarget);}
  }

  inline void getCurrentPdGains(Eigen::VectorXd &pGain, Eigen::VectorXd &dGain){
    if(cfbMode_ == CfbTargetMode::FROM_JOINT){getPdGains(pGain, dGain);}
    else if(cfbMode_ == CfbTargetMode::FROM_ACTUATOR){getActuatorPdGains(pGain, dGain);}
  }

  // cfbMode_ (FROM_JOINT, FROM_ACTUATOR) and corresponding get/setCfbMode()
  // setActuatorGeneralizedForce()

  // getters that need special care
  /**
   * @return the upper joint torque/force limit*/
  [[nodiscard]] const raisim::VecDyn getActuationUpperLimits() const {
    const_cast<ArticulatedRaipal*>(this)->updateRaipal();
    raisim::VecDyn tauUpper;
    tauUpper = maxTorque_;
    return tauUpper;
  }

  /**
   * @return the lower joint torque/force limit*/
  [[nodiscard]] const raisim::VecDyn getActuationLowerLimits() const {
    const_cast<ArticulatedRaipal*>(this)->updateRaipal();
    raisim::VecDyn tauLower;
    tauLower = minTorque_;
    return tauLower;
  }

  [[nodiscard]] const raisim::VecDyn getFeedForwardGeneralizedForce() const {
    const_cast<ArticulatedRaipal*>(this)->updateRaipal();
    raisim::VecDyn tauFF;
    tauFF = tauFF_;
    return tauFF;
  }

  [[nodiscard]] raisim::VecDyn getGeneralizedForce() const {
    const_cast<ArticulatedRaipal*>(this)->updateRaipal();
    raisim::VecDyn genForce;
    genForce = articulated_->getGeneralizedForce().e() - tauCFB_;
    return genForce;
  }

  /**
   * get the coriolis and the gravitational term
   * @param[in] gravity gravitational acceleration. You should get this value from the world.getGravity();
   * @return the coriolis and the gravitational term. Check Object/ArticulatedSystem section in the manual */
  [[nodiscard]] const raisim::VecDyn getNonlinearities(const raisim::Vec<3>& gravity) {
    raisim::VecDyn b;
    b = articulated_->getNonlinearities(gravity).e() - tauCFB_; // when all else is zero, -tauCFB_ is required to cancel out CFB-induced residual force
    return b;
  }

  // setters that need special care

  // Explicit wrappers avoid auto return deduction ordering issues.
  void getPdTarget(Eigen::VectorXd &posTarget, Eigen::VectorXd &velTarget) {
    updateRaipal();
    articulated_->getPdTarget(posTarget, velTarget);
  }

  void getPdTarget(Eigen::VectorXd &posTarget, Eigen::VectorXd &velTarget) const {
    const_cast<ArticulatedRaipal*>(this)->updateRaipal();
    articulated_->getPdTarget(posTarget, velTarget);
  }

  void getPdGains(Eigen::VectorXd &pGain, Eigen::VectorXd &dGain) {
    updateRaipal();
    articulated_->getPdGains(pGain, dGain);
  }

  void getPdGains(Eigen::VectorXd &pGain, Eigen::VectorXd &dGain) const {
    const_cast<ArticulatedRaipal*>(this)->updateRaipal();
    articulated_->getPdGains(pGain, dGain);
  }

  // RAIPAL_SET(setRotorInertia)
  void setRotorInertia(const raisim::VecDyn &rotorInertia) {
    rotorInertia_ = rotorInertia;
    resetUpdateFlag();
  }

  // RAIPAL_SET(setActuationLimits)
  void setActuationLimits(const Eigen::VectorXd &upper, const Eigen::VectorXd &lower) {
    maxTorque_ = upper;
    minTorque_ = lower;
    resetUpdateFlag();
  }

  // RAIPAL_SET(setJointVelocityLimits)
  void setJointVelocityLimits(const Eigen::VectorXd &velLimits) {
    maxVelocity_ = velLimits;
    resetUpdateFlag();
  }

  // RAIPAL_SET(setGeneralizedForce)      // tauFF_
  void setGeneralizedForce(const raisim::VecDyn &tau) { 
    tauFF_ = tau.e();
    tauCFB_.setZero();
    articulated_->setGeneralizedForce(tauFF_);
    resetUpdateFlag();
  }

  void setGeneralizedForce(const Eigen::VectorXd &tau) { 
    tauFF_ = tau;
    tauCFB_.setZero();
    articulated_->setGeneralizedForce(tauFF_);
    resetUpdateFlag();
  }

  // RAIPAL_SET(setPdTarget)              // pTarget, dTarget
  /**
   * set PD targets.
   * @param[in] posTarget position target
   * @param[in] velTarget velocity target */
  void setPdTarget(const Eigen::VectorXd &posTarget, const Eigen::VectorXd &velTarget) {
    RSFATAL_IF(cfbMode_ != CfbTargetMode::FROM_JOINT, "setPdTarget can only be used when CFB target mode is set to FROM_JOINT");
    articulated_->setPdTarget(posTarget, velTarget);
    resetUpdateFlag();
  }
  
  // RAIPAL_SET(setPdGains)               // pGain, dGain
  /**
   * set PD gains.
   * @param[in] pgain position gain (dimension == getDOF())
   * @param[in] dgain velocity gain (dimension == getDOF())*/
  void setPdGains(const Eigen::VectorXd &pgain, const Eigen::VectorXd &dgain) {
    RSFATAL_IF(cfbMode_ != CfbTargetMode::FROM_JOINT, "setPdGains can only be used when CFB target mode is set to FROM_JOINT");
    articulated_->setPdGains(pgain, dgain);
    resetUpdateFlag();
  }

  // setters that require updates
  RAIPAL_SET(setState)                 // gc, gv
  RAIPAL_SET(setGeneralizedCoordinate) // gc 
  RAIPAL_SET(setGeneralizedVelocity)   // gv

  // who uses these anyways?
  // RAIPAL_SET(setPTarget)               // pTarget
  // RAIPAL_SET(setDTarget)               // dTarget
  // RAIPAL_SET(setPGains)                // pGain
  // RAIPAL_SET(setDGains)                // dGain

  // forwarded to ArticluatedSystem (no need for wrapper intervention)
  RAIPAL_FORWARD(isSecondOrderOrHigher)
  RAIPAL_UPDATE_GET(getGeneralizedCoordinate)
  RAIPAL_UPDATE_GET(getGeneralizedVelocity)
  RAIPAL_UPDATE_GET(getGeneralizedAcceleration)
  RAIPAL_UPDATE_GET(getBaseOrientation)
  RAIPAL_UPDATE_GET(getBasePosition)
  RAIPAL_FORWARD(updateKinematics)
  RAIPAL_UPDATE_GET(getState)
  RAIPAL_UPDATE_GET(getMassMatrix)
  RAIPAL_UPDATE_GET(getNonlinearities)
  RAIPAL_UPDATE_GET(getInverseMassMatrix)
  RAIPAL_UPDATE_GET(getCompositeCOM)
  RAIPAL_UPDATE_GET(getCOM)
  RAIPAL_UPDATE_GET(getCompositeInertia)
  RAIPAL_UPDATE_GET(getCompositeMass)
  RAIPAL_UPDATE_GET(getLinearMomentum)
  RAIPAL_UPDATE_GET(getGeneralizedMomentum)
  RAIPAL_UPDATE_GET(getEnergy)
  RAIPAL_UPDATE_GET(getKineticEnergy)
  RAIPAL_UPDATE_GET(getPotentialEnergy)
  RAIPAL_UPDATE_GET(getAngularMomentum)
  RAIPAL_FORWARD(printOutBodyNamesInOrder)
  RAIPAL_FORWARD(printOutMovableJointNamesInOrder)
  RAIPAL_FORWARD(printOutFrameNamesInOrder)
  RAIPAL_UPDATE_GET(getMovableJointNames)
  RAIPAL_UPDATE_GET(getGeneralizedVelocityIndex)
  RAIPAL_UPDATE_GET(getPosition)
  RAIPAL_UPDATE_GET(getFrameByName)
  RAIPAL_UPDATE_GET(getFrameByLinkName)
  RAIPAL_UPDATE_GET(getFrameIdxByLinkName)
  RAIPAL_UPDATE_GET(getFrameByIdx)
  RAIPAL_UPDATE_GET(getFrameIdxByName)
  RAIPAL_UPDATE_GET(getFrames)
  RAIPAL_UPDATE_GET(getFramePosition)
  RAIPAL_UPDATE_GET(getPositionInFrame)
  RAIPAL_UPDATE_GET(getFrameOrientation)
  RAIPAL_UPDATE_GET(getFrameVelocity)
  RAIPAL_UPDATE_GET(getFrameAngularVelocity)
  RAIPAL_UPDATE_GET(getFrameAcceleration)
  RAIPAL_UPDATE_GET(getPositionInBodyCoordinate)
  RAIPAL_UPDATE_GET(getOrientation)
  RAIPAL_UPDATE_GET(getVelocity)
  RAIPAL_UPDATE_GET(getAngularVelocity)
  RAIPAL_UPDATE_GET(getSparseJacobian)
  RAIPAL_UPDATE_GET(getSparseRotationalJacobian)
  RAIPAL_UPDATE_GET(getTimeDerivativeOfSparseJacobian)
  RAIPAL_UPDATE_GET(getTimeDerivativeOfSparseRotationalJacobian)
  RAIPAL_UPDATE_GET(getDenseJacobian)
  RAIPAL_UPDATE_GET(getDenseRotationalJacobian)
  RAIPAL_UPDATE_GET(getDenseFrameJacobian)
  RAIPAL_UPDATE_GET(getDenseFrameRotationalJacobian)
  RAIPAL_UPDATE_GET(getBodyIdx)
  RAIPAL_UPDATE_GET(getDOF)
  RAIPAL_UPDATE_GET(getGeneralizedVelocityDim)
  RAIPAL_UPDATE_GET(getGeneralizedCoordinateDim)
  RAIPAL_UPDATE_GET(getBodyPose)
  RAIPAL_UPDATE_GET(getBodyPosition)
  RAIPAL_UPDATE_GET(getBodyOrientation)
  RAIPAL_UPDATE_GET(getJointPos_P)
  RAIPAL_UPDATE_GET(getJointOrientation_P)
  RAIPAL_UPDATE_GET(getJointAxis_P)
  RAIPAL_UPDATE_GET(getJointAxis)
  RAIPAL_UPDATE_GET(getMass)
  RAIPAL_UPDATE_GET(getInertia)
  RAIPAL_UPDATE_GET(getBodyCOM_B)
  RAIPAL_UPDATE_GET(getBodyCOM_W)
  RAIPAL_UPDATE_GET(getCollisionBodies)
  RAIPAL_UPDATE_GET(getCollisionBody)
  RAIPAL_FORWARD(updateMassInfo)
  RAIPAL_FORWARD(setMass)
  RAIPAL_UPDATE_GET(getTotalMass)
  RAIPAL_FORWARD(setExternalForce)
  RAIPAL_FORWARD(setExternalTorque)
  RAIPAL_FORWARD(setExternalTorqueInBodyFrame)
  RAIPAL_UPDATE_GET(getContactPointVel)
  RAIPAL_FORWARD(setControlMode)
  RAIPAL_UPDATE_GET(getControlMode)
  RAIPAL_FORWARD(setJointDamping)
  RAIPAL_FORWARD(computeSparseInverse)
  RAIPAL_FORWARD(massMatrixVecMul)
  RAIPAL_FORWARD(ignoreCollisionBetween)
  RAIPAL_UPDATE_GET(getOptions)
  RAIPAL_UPDATE_GET(getBodyNames)
  RAIPAL_UPDATE_GET(getVisOb)
  RAIPAL_UPDATE_GET(getVisColOb)
  RAIPAL_UPDATE_GET(getVisObPose)
  RAIPAL_UPDATE_GET(getVisColObPose)
  RAIPAL_UPDATE_GET(getResourceDir)
  RAIPAL_UPDATE_GET(getRobotDescriptionfFileName)
  RAIPAL_UPDATE_GET(getRobotDescriptionfTopDirName)
  RAIPAL_UPDATE_GET(getRobotDescriptionFullPath)
  RAIPAL_UPDATE_GET(getRobotDescription)
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
  RAIPAL_UPDATE_GET(getRotorInertia)
  RAIPAL_UPDATE_GET(getJointType)
  RAIPAL_UPDATE_GET(getNumberOfJoints)
  RAIPAL_UPDATE_GET(getJoint)
  RAIPAL_UPDATE_GET(getLink)
  RAIPAL_UPDATE_GET(getMappingFromBodyIndexToGeneralizedVelocityIndex)
  RAIPAL_UPDATE_GET(getMappingFromBodyIndexToGeneralizedCoordinateIndex)
  RAIPAL_UPDATE_GET(getObjectType)
  RAIPAL_UPDATE_GET(getBodyType)
  RAIPAL_FORWARD(setIntegrationScheme)
  RAIPAL_UPDATE_GET(getJointLimitViolations)
  RAIPAL_FORWARD(setJointLimits)
  RAIPAL_UPDATE_GET(getJointLimits)
  RAIPAL_UPDATE_GET(getJointVelocityLimits)
  RAIPAL_FORWARD(clearExternalForcesAndTorques)
  RAIPAL_FORWARD(addSpring)
  RAIPAL_UPDATE_GET(getSprings)
  RAIPAL_UPDATE_GET(getParentVector)
  RAIPAL_UPDATE_GET(getSensorSet)
  RAIPAL_UPDATE_GET(getSensorSets)
  RAIPAL_FORWARD(addConstraints)
  RAIPAL_FORWARD(initializeConstraints)
  RAIPAL_FORWARD(setComputeInverseDynamics)
  RAIPAL_UPDATE_GET(getComputeInverseDynamics)
  RAIPAL_UPDATE_GET(getForceAtJointInWorldFrame)
  RAIPAL_UPDATE_GET(getTorqueAtJointInWorldFrame)
  RAIPAL_UPDATE_GET(getAllowedNumberOfInternalContactsBetweenTwoBodies)
  RAIPAL_FORWARD(setAllowedNumberOfInternalContactsBetweenTwoBodies)
  RAIPAL_FORWARD(articulatedBodyAlgorithm)
  RAIPAL_UPDATE_GET(getMinvJT)
  RAIPAL_UPDATE_GET(getj_MinvJT_T1D)
  RAIPAL_UPDATE_GET(getUdot)
  RAIPAL_UPDATE_GET(getFullDelassusAndTauStar)
  RAIPAL_FORWARD(appendJointLimits)

  void setGeneralizedCoordinate(std::initializer_list<double> jointState) {
    articulated_->setGeneralizedCoordinate(jointState);
    resetUpdateFlag();
  }

  void setGeneralizedVelocity(std::initializer_list<double> jointVelocity) {
    articulated_->setGeneralizedVelocity(jointVelocity);
    resetUpdateFlag();
  }

  void setGeneralizedForce(std::initializer_list<double> tau) {
    articulated_->setGeneralizedForce(tau);
    resetUpdateFlag();
  }

  static void convertSparseJacobianToDense(
      const raisim::SparseJacobian& sparseJaco,
      Eigen::MatrixXd& denseJaco) {
    raisim::ArticulatedSystem::convertSparseJacobianToDense(sparseJaco, denseJaco);
  }

private:

  struct actuatorState{
    size_t idx;
    int sign;

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

  raisim::ArticulatedSystem* articulated_ = nullptr;
  size_t dof_;

  std::vector<actuatorState> actuatorStates_;
  CfbTargetMode cfbMode_ = CfbTargetMode::FROM_JOINT;

  bool updated_ = false;
  
  raisim::VecDyn rotorInertia_, maxVelocity_;
  Eigen::VectorXd maxTorque_, minTorque_;

  Eigen::VectorXd tauFF_;   // feedforward term for rest of robot
  Eigen::VectorXd tauCFB_;  // feedforward term for CFB correction
};

#undef RAIPAL_FORWARD
#undef RAIPAL_UPDATE_GET
#undef RAIPAL_SET
