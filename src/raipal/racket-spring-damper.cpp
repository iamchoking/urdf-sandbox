#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "raisim/RaisimServer.hpp"
#include "raisim/World.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace {

struct Options {
  bool launchServer = true;
  size_t steps = 10000;
  double timeStep = 0.001;
  double stiffness = 60.0;
  double damping = 2.0;
  double restLength = 0.01;
  double maxTension = 80.0;
  double sphereRadius = 0.035;
  double sphereMass = 0.05;
  bool matchTipVelocity = false;
  bool customDampedWire = false;
  std::string urdf = "/raipal/urdf/raipal_racket_L.urdf";
  std::string tipFrameName = "LE_tip_fixed";
};

bool startsWith(const std::string& value, const std::string& prefix) {
  return value.rfind(prefix, 0) == 0;
}

Options parseOptions(int argc, char* argv[]) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    if (arg == "--no-server") {
      options.launchServer = false;
    } else if (arg == "--match-tip-velocity") {
      options.matchTipVelocity = true;
    } else if (arg == "--custom-damped-wire") {
      options.customDampedWire = true;
    } else if (startsWith(arg, "--steps=")) {
      options.steps = std::stoull(arg.substr(std::string("--steps=").size()));
    } else if (startsWith(arg, "--dt=")) {
      options.timeStep = std::stod(arg.substr(std::string("--dt=").size()));
    } else if (startsWith(arg, "--k=")) {
      options.stiffness = std::stod(arg.substr(std::string("--k=").size()));
    } else if (startsWith(arg, "--c=")) {
      options.damping = std::stod(arg.substr(std::string("--c=").size()));
    } else if (startsWith(arg, "--rest=")) {
      options.restLength = std::stod(arg.substr(std::string("--rest=").size()));
    } else if (startsWith(arg, "--max-tension=")) {
      options.maxTension = std::stod(arg.substr(std::string("--max-tension=").size()));
    } else if (startsWith(arg, "--sphere-radius=")) {
      options.sphereRadius = std::stod(arg.substr(std::string("--sphere-radius=").size()));
    } else if (startsWith(arg, "--sphere-mass=")) {
      options.sphereMass = std::stod(arg.substr(std::string("--sphere-mass=").size()));
    } else if (startsWith(arg, "--urdf=")) {
      options.urdf = arg.substr(std::string("--urdf=").size());
    } else if (startsWith(arg, "--frame=")) {
      options.tipFrameName = arg.substr(std::string("--frame=").size());
    } else {
      std::cout << "Ignoring unknown option: " << arg << std::endl;
    }
  }
  return options;
}

Eigen::Vector3d eigen(const raisim::Vec<3>& value) {
  return value.e();
}

double computeCableTension(
    raisim::Object* bodyA,
    size_t bodyAIdx,
    const raisim::Vec<3>& mountAInBody,
    raisim::Object* bodyB,
    size_t bodyBIdx,
    const raisim::Vec<3>& mountBInBody,
    double restLength,
    double stiffness,
    double damping,
    double maxTension,
    double* distanceOut,
    double* stretchRateOut) {
  raisim::Vec<3> posA;
  raisim::Vec<3> posB;
  raisim::Vec<3> velA;
  raisim::Vec<3> velB;

  bodyA->getPosition(bodyAIdx, mountAInBody, posA);
  bodyB->getPosition(bodyBIdx, mountBInBody, posB);
  bodyA->getVelocity(bodyAIdx, mountAInBody, velA);
  bodyB->getVelocity(bodyBIdx, mountBInBody, velB);

  const Eigen::Vector3d delta = eigen(posB) - eigen(posA);
  const double distance = delta.norm();
  Eigen::Vector3d direction = Eigen::Vector3d::UnitX();
  if (distance > 1e-9) {
    direction = delta / distance;
  }
  const double stretchRate = (eigen(velB) - eigen(velA)).dot(direction);
  const double rawTension = stiffness * (distance - restLength) + damping * stretchRate;
  const double tension = std::clamp(rawTension, 0.0, maxTension);

  if (distanceOut) {
    *distanceOut = distance;
  }
  if (stretchRateOut) {
    *stretchRateOut = stretchRate;
  }
  return tension;
}

}  // namespace

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);
  (void)binaryPath;

  const Options options = parseOptions(argc, argv);

  raisim::World world;
  world.setTimeStep(options.timeStep);
  world.setGravity({0.0, 0.0, 0.0});

  auto* raipal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + options.urdf);
  raipal->setName("raipal_with_tip_spring");
  raipal->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);

  const auto initialCoordinateDim = raipal->getGeneralizedCoordinateDim();
  const auto initialVelocityDim = raipal->getGeneralizedVelocityDim();

  Eigen::VectorXd gc = Eigen::VectorXd::Zero(initialCoordinateDim);
  Eigen::VectorXd gv = Eigen::VectorXd::Zero(initialVelocityDim);
  if (initialVelocityDim >= 1) {
    gv[0] = 0.7;
  }
  if (initialVelocityDim >= 2) {
    gv[1] = -0.4;
  }

  raipal->setState(gc, gv);
  raipal->setGeneralizedForce(Eigen::VectorXd::Zero(initialVelocityDim));

  const size_t tipFrameIdx = raipal->getFrameIdxByName(options.tipFrameName);
  if (tipFrameIdx == size_t(-1)) {
    std::cerr << "Could not find frame: " << options.tipFrameName << std::endl;
    return 1;
  }

  const auto& tipFrame = raipal->getFrameByIdx(tipFrameIdx);
  const size_t tipBodyIdx = tipFrame.parentId;
  const raisim::Vec<3> tipMountInBody = tipFrame.position;

  raisim::Vec<3> tipWorld;
  raisim::Vec<3> tipVelocity;
  raisim::Mat<3, 3> tipOrientation;
  raipal->getPosition(tipBodyIdx, tipMountInBody, tipWorld);
  raipal->getVelocity(tipBodyIdx, tipMountInBody, tipVelocity);
  raipal->getFrameOrientation(tipFrameIdx, tipOrientation);

  auto* payload = world.addSphere(options.sphereRadius, options.sphereMass, "default", 1, 0);
  payload->setName("wire_payload_sphere");
  Eigen::Vector3d payloadOffset =
      tipOrientation.e() * Eigen::Vector3d(options.restLength, 0.0, 0.0);
  if (!options.matchTipVelocity && eigen(tipVelocity).dot(payloadOffset) > 0.0) {
    payloadOffset = -payloadOffset;
  }
  const Eigen::Vector3d initialPayloadPosition = eigen(tipWorld) + payloadOffset;
  payload->setPosition(initialPayloadPosition);
  const raisim::Vec<3> zeroVelocity = {0.0, 0.0, 0.0};
  payload->setVelocity(options.matchTipVelocity ? tipVelocity : zeroVelocity, zeroVelocity);

  const raisim::Vec<3> payloadMountInBody = {0.0, 0.0, 0.0};
  raisim::LengthConstraint* wire = nullptr;
  raisim::CompliantLengthConstraint* compliantWire = nullptr;
  raisim::CustomLengthConstraint* customWire = nullptr;
  if (options.customDampedWire) {
    customWire = world.addCustomWire(
        payload,
        0,
        payloadMountInBody,
        raipal,
        tipBodyIdx,
        tipMountInBody,
        options.restLength);
    wire = customWire;
  } else {
    compliantWire = world.addCompliantWire(
        payload,
        0,
        payloadMountInBody,
        raipal,
        tipBodyIdx,
        tipMountInBody,
        options.restLength,
        options.stiffness);
    wire = compliantWire;
  }
  wire->setName(options.customDampedWire ? "tip_custom_spring_damper" : "tip_compliant_spring");
  wire->setVisualizationWidth(0.01);
  wire->setStretchType(raisim::LengthConstraint::StretchType::STRETCH_RESISTANT_ONLY);

  std::cout << "RaiSim wire payload demo" << std::endl;
  std::cout << "  URDF: " << options.urdf << std::endl;
  std::cout << "  frame: " << options.tipFrameName
            << " on body index " << tipBodyIdx
            << " (" << tipFrame.parentName << " -> " << tipFrame.bodyName << ")" << std::endl;
  std::cout << "  qdim before/after add wire: " << initialCoordinateDim
            << " / " << raipal->getGeneralizedCoordinateDim() << std::endl;
  std::cout << "  vdim before/after add wire: " << initialVelocityDim
            << " / " << raipal->getGeneralizedVelocityDim() << std::endl;
  std::cout << "  world wires: " << world.getWires().size() << std::endl;
  std::cout << "  wire mode: "
            << (options.customDampedWire ? "custom damped, manual tension" : "compliant, Raisim-computed tension")
            << std::endl;
  std::cout << "  payload sphere: mass=" << options.sphereMass
            << " radius=" << options.sphereRadius
            << " collision mask disabled"
            << " initial velocity=" << (options.matchTipVelocity ? "tip" : "zero") << std::endl;
  std::cout << "  k=" << options.stiffness
            << " c=" << (options.customDampedWire ? options.damping : 0.0)
            << " rest=" << options.restLength
            << " max tension=" << (options.customDampedWire ? options.maxTension : 0.0) << std::endl;
  std::cout << "  Use --custom-damped-wire to compute and set tension manually." << std::endl;
  std::cout << "  Use --no-server for a headless run." << std::endl;

  raisim::RaisimServer server(&world);
  if (options.launchServer) {
    server.launchServer();
    server.focusOn(raipal);
  }

  world.integrate1();

  const size_t logEvery = std::max<size_t>(1, options.steps / 10);
  double distance = 0.0;
  double stretchRate = 0.0;
  double tension = 0.0;

  for (size_t step = 0; step < options.steps; ++step) {
    if (customWire) {
      tension = computeCableTension(
          payload,
          0,
          payloadMountInBody,
          raipal,
          tipBodyIdx,
          tipMountInBody,
          options.restLength,
          options.stiffness,
          options.damping,
          options.maxTension,
          &distance,
          &stretchRate);
      customWire->setTension(tension);
    } else {
      computeCableTension(
          payload,
          0,
          payloadMountInBody,
          raipal,
          tipBodyIdx,
          tipMountInBody,
          options.restLength,
          0.0,
          0.0,
          0.0,
          &distance,
          &stretchRate);
    }

    if (options.launchServer) {
      RS_TIMED_LOOP(1e6 * world.getTimeStep())
      server.integrateWorldThreadSafe();
    } else {
      world.integrate();
    }
    if (compliantWire) {
      tension = compliantWire->getTension().norm();
    }

    if (step % logEvery == 0 || step + 1 == options.steps) {
      std::cout << "  t=" << world.getWorldTime()
                << " distance=" << distance
                << " stretch_rate=" << stretchRate
                << " tension=" << tension << std::endl;
    }
  }

  if (options.launchServer) {
    server.killServer();
  }

  return 0;
}
