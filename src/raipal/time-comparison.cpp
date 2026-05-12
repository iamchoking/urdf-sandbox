#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "raisim/World.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>

#include <raipal_kinematics/raipal_cfb.hpp>
#include "raisimRaipal/Raipal.hpp"
#include "random_coordinates.hpp"

namespace rk9 = raipal::kinematics;

static constexpr double SIM_TIMESTEP = 0.0001;
static constexpr double TEST_DURATION = 100.0;

struct TimingStats {
  size_t count = 0;
  double mean = 0.0;
  double m2 = 0.0;
  double max = 0.0;

  void add(double value) {
    ++count;
    const double delta = value - mean;
    mean += delta / static_cast<double>(count);
    const double delta2 = value - mean;
    m2 += delta * delta2;
    if (value > max) {
      max = value;
    }
  }

  double variance() const {
    return count > 1 ? m2 / static_cast<double>(count - 1) : 0.0;
  }

  double stddev() const {
    return std::sqrt(variance());
  }

  double total() const {
    return mean * static_cast<double>(count);
  }
};

static void printStats(const std::string &label, const TimingStats &stats) {
  const double meanMs = stats.mean * 1e3;
  const double stdMs = stats.stddev() * 1e3;
  const double maxMs = stats.max * 1e3;
  const double totalMs = stats.total() * 1e3;

  std::cout
      << label << " world.integrate timing (ms)"
      << " | avg: " << meanMs
      << " | std: " << stdMs
      << " | max: " << maxMs
      << " | total: " << totalMs
      << std::endl;
}

static TimingStats runRaipal9Test() {
  raisim::World world;
  world.setTimeStep(SIM_TIMESTEP);

  auto raipal9 = world.addArticulatedSystem(
      std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal9/urdf/raipal_stub-0_R.urdf");

  Eigen::VectorXd gc9(9), gv9(9), pTarget9(9), dTarget9(9), pGain9(9), dGain9(9);

  auto jointLimits9 = raipal9->getJointLimits();

  Eigen::VectorXd sweepCenter9(9), sweepAmplitude9(9), sweepLimits9(9);
  Eigen::VectorXd padRatio9(9), minAmplitudeRatio9(9), maxAmplitudeRatio9(9);

  padRatio9           << 0.1, 0.1, 0.1, 0.1, 0.0, 0.0, 0.1, 0.1, 0.1;
  minAmplitudeRatio9  << 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0;
  maxAmplitudeRatio9  << 0.2, 0.2, 0.2, 0.4, 0.0, 0.0, 0.2, 0.2, 0.2;

  utils::sampleJointSweep(
      sweepCenter9,
      sweepAmplitude9,
      jointLimits9,
      padRatio9,
      minAmplitudeRatio9,
      maxAmplitudeRatio9);

  rk9::cfbForward(sweepCenter9);

  gc9 = sweepCenter9;
  gv9.setZero();

  raipal9->setState(gc9, gv9);

  pTarget9 = gc9;
  dTarget9 = gv9;

  raipal9->setPdTarget(pTarget9, dTarget9);

  pGain9 << 100, 100, 100, 100, 0, 0, 100, 100, 100;
  dGain9 << 10 , 10 , 10 , 10 , 0, 0, 10 , 10 , 10 ;

  raipal9->setPdGains(pGain9, dGain9);

  std::array<double, 2> freq = {1.0, 5.0};
  const size_t steps = static_cast<size_t>(TEST_DURATION / world.getTimeStep());

  TimingStats stats;
  double theta = 0.0;

  for (size_t t = 0; t < steps; ++t) {
    const double currentFreq = freq[0] + (freq[1] - freq[0]) *
        (static_cast<double>(t) / static_cast<double>(steps));
    theta += 2.0 * M_PI * currentFreq * world.getTimeStep();

    pTarget9 = sweepCenter9 + std::sin(theta) * sweepAmplitude9;
    rk9::cfbForward(pTarget9);

    raipal9->setPdTarget(pTarget9, dTarget9);

    const auto start = std::chrono::steady_clock::now();
    world.integrate();
    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed = end - start;

    stats.add(elapsed.count());
  }

  return stats;
}

static TimingStats runRaipal7Test() {
  raisim::World world;
  world.setTimeStep(SIM_TIMESTEP);

  raisim::Raipal raipal7(
      world.addArticulatedSystem(
          std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_stub-0_L.urdf"),
      {3},
      {-1});

    Eigen::VectorXd gc7(7), gv7(7), pTarget7(7), dTarget7(7), pGain7(7), dGain7(7);

    auto jointLimits7 = raipal7->getJointLimits();

    Eigen::VectorXd sweepCenter7(7), sweepAmplitude7(7);
    Eigen::VectorXd padRatio7(7), minAmplitudeRatio7(7), maxAmplitudeRatio7(7);

    padRatio7           << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    minAmplitudeRatio7  << 0.0, 0.0, 0.0, 0.3, 0.0, 0.0, 0.0;
    maxAmplitudeRatio7  << 0.2, 0.2, 0.2, 0.4, 0.2, 0.2, 0.2;

    utils::sampleJointSweep(
      sweepCenter7,
      sweepAmplitude7,
      jointLimits7,
      padRatio7,
      minAmplitudeRatio7,
      maxAmplitudeRatio7);

    gc7 = sweepCenter7;
    gv7.setZero();

  raipal7->setState(gc7, gv7);

  pTarget7 = gc7;
  dTarget7 = gv7;

  raipal7->setPdTarget(pTarget7, dTarget7);
  raipal7->setCfbTargetFromActuator();

  pGain7 << 100, 100, 100, 100, 100, 100, 100;
  dGain7 << 10 , 10 , 10 , 10 , 10 , 10 , 10 ;

  raipal7->setActuatorPdGains(pGain7, dGain7);

  std::array<double, 2> freq = {1.0, 5.0};
  const size_t steps = static_cast<size_t>(TEST_DURATION / world.getTimeStep());

  TimingStats stats;
  double theta = 0.0;

  for (size_t t = 0; t < steps; ++t) {
    const double currentFreq = freq[0] + (freq[1] - freq[0]) *
        (static_cast<double>(t) / static_cast<double>(steps));
    theta += 2.0 * M_PI * currentFreq * world.getTimeStep();

    pTarget7 = sweepCenter7 + std::sin(theta) * sweepAmplitude7;

    raipal7->setActuatorPdTarget(pTarget7, dTarget7);

    raipal7->updateRaipal();
    const auto start = std::chrono::steady_clock::now();
    world.integrate();
    const auto end = std::chrono::steady_clock::now();
    raipal7->resetUpdateFlags();

    const std::chrono::duration<double> elapsed = end - start;
    stats.add(elapsed.count());
  }

  return stats;
}

int main(int argc, char* argv[]) {
  (void)argc;
  (void)argv;

  std::cout << "Running time comparison (no visualization)..." << std::endl;

  const TimingStats stats9 = runRaipal9Test();
  const TimingStats stats7 = runRaipal7Test();

  printStats("raipal9", stats9);
  printStats("raipal7", stats7);

  return 0;
}
