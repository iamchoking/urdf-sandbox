#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"
#include <raipal_kinematics/raipal_cfb.hpp>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

size_t TOTAL_STEPS = 20000;
double logging_dt = 0.001;
bool FORWARD = false;

int main(int argc, char* argv[]) {
  
  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  Eigen::VectorXd desired_vel(3);
  Eigen::VectorXd desired_pose(7);

  desired_vel << 28.0, 0.0, 0.0; // desired velocity in base frame
  desired_pose << 1.303674,-0.871144,-0.796485,0.515208,0.973577,-0.462641,0.251565;

  Eigen::VectorXd snap_torque(9);
  snap_torque << 180.0, 180.0, 180.0, 180.0, 0.0, 0.0, 13.0, 13.0, 13.0;

  double window = 0.005; // seconds with no torque

  Eigen::VectorXd hold_p_gain(9), hold_d_gain(9);
  hold_p_gain << 400.0, 400.0, 400.0, 700.0, 0.0, 0.0, 50.0,  120.0,  75.0;
  hold_d_gain << 40.0,  40.0,  40.0,  70.0,  0.0, 0.0, 5.0,  12.0,  7.5;
  // hold_p_gain.setZero();
  // hold_d_gain.setZero();

  auto raipal_L = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal9/urdf/raipal_stub-0_L.urdf");

  raipal_L->setName("raipal");

  std::cout << "robot was loaded!" << std::endl;  

  int dof_ = raipal_L->getDOF();
  int joints_ = dof_ == 9 ? 7 : 14;

  const size_t frameIdx = raipal_L->getFrameIdxByName("LE_tip_fixed");
  Eigen::VectorXd desired_gc(9);
  Eigen::VectorXd desired_gv(9);
  ::raipal::kinematics::cfbExpand(desired_pose, desired_gc, true);

  raipal_L->setGeneralizedCoordinate(desired_gc);

  Eigen::MatrixXd J_pos9(3,9), J_pos7(3,7);  
  J_pos9.setZero();
  J_pos7.setZero();

  raipal_L->getDenseFrameJacobian(frameIdx, J_pos9);  

  ::raipal::kinematics::cfbJacobian(desired_gc,J_pos9,J_pos7);

  desired_gv = J_pos9.completeOrthogonalDecomposition().solve(desired_vel);
  std::cout << "desired gc: " << desired_gc.transpose() << std::endl;
  std::cout << "desired gv: " << desired_gv.transpose() << std::endl;

  // raipal -> setComputeInverseDynamics(true);

  // world.addGround();
  world.setTimeStep(0.0001);

  // Declare variables (should be in private section)
  Eigen::VectorXd gc_, gv_, pTarget_, dTarget_, pGain_, dGain_;

  /// initialize containers
  gc_.setZero(dof_);
  gv_.setZero(dof_);

  pTarget_.setZero(dof_);
  dTarget_.setZero(dof_);

  pGain_.setZero(dof_);
  dGain_.setZero(dof_);

  // std::cout << "Setting pd gains..." << std::endl;
  raipal_L->setPdGains(pGain_, dGain_);
  // std::cout << "Pd gains set." << std::endl;
  raipal_L->setGeneralizedForce(Eigen::VectorXd::Zero(dof_));


  server.launchServer();
  server.focusOn(raipal_L);

  /// if you are using an old version of Raisim, you need this line
  world.integrate1();

  std::cout << "buildup_simulation" << std::endl;
  const double velocity_sign = FORWARD ? 1.0 : -1.0;
  raipal_L->setState(desired_gc, velocity_sign * desired_gv);
  std::cout << "Simulation mode: " << (FORWARD ? "FORWARD" : "REVERSE")
            << " (initial gv sign: " << velocity_sign << ")" << std::endl;

  for (int sec=3; sec>0; sec--){
    std::cout << "Starting in [" << sec << "]..." << std::endl;
    raisim::USLEEP(1000000);
  }
  std::cout << "START!" << std::endl;

  auto jointLimits = raipal_L->getJointLimits();
  const double vel_stop_thresh = 5e-2;
  const double limit_eps = 1e-2;
  const double dt = world.getTimeStep();

  std::vector<Eigen::VectorXd> gc_log;
  std::vector<double> time_log;
  gc_log.reserve(TOTAL_STEPS + 1);
  time_log.reserve(TOTAL_STEPS + 1);

  double next_log_t = 0.0;
  bool trajectory_invalid = false;
  std::string invalid_reason;

  raipal_L->getState(gc_, gv_);
  gc_log.push_back(gc_);
  time_log.push_back(FORWARD ? 0.0 : -0.0);
  next_log_t += logging_dt;

  std::vector<bool> joint_stopped(dof_, false);
  std::vector<bool> joint_stopping(dof_, false);
  Eigen::VectorXd gv_prev = gv_;
  pTarget_ = gc_;
  dTarget_.setZero();

  // tau warmup
  double warmup_alpha = 0.995;
  double warmup_factor = 1.0;

  // SIM LOOP
  size_t executed_steps = 0;
  for (size_t t = 0; t<TOTAL_STEPS; t++){
    RS_TIMED_LOOP(world.getTimeStep()*2e6)

    pGain_.setZero();
    dGain_.setZero();
    dTarget_.setZero();

    for (int i = 0; i < dof_; i++) {
      const bool below_thresh = std::abs(gv_[i]) < vel_stop_thresh;
      const bool sign_changed = (gv_prev[i] * gv_[i]) < 0.0;

      const double kd_hold = (i < hold_d_gain.size()) ? hold_d_gain[i] : 0.0;
      const bool has_damping_mode = kd_hold > 1e-12;
      const double damping_entry_speed = has_damping_mode ? (snap_torque[i] / kd_hold) : 0.0;

      if (!joint_stopped[i] && !joint_stopping[i] && has_damping_mode && (std::abs(gv_[i]) <= damping_entry_speed)) {
        joint_stopping[i] = true;
        std::cout << "Joint " << i << " entering damping-stop at t=" << (t * dt)
                  << "s, dq=" << gv_[i] << ", |dq|_entry=" << damping_entry_speed << std::endl;
      }

      if (!joint_stopped[i] && (below_thresh || sign_changed)) {
        joint_stopped[i] = true;
        joint_stopping[i] = false;
        pTarget_[i] = gc_[i];
        std::cout << "Joint " << i << " stopped at t=" << (t * dt)
                  << "s, q=" << gc_[i] << ", dq=" << gv_[i]
                  << (sign_changed ? " (sign change detected)" : "") << std::endl;
      }

      if (!joint_stopped[i]) {
        // Keep non-stopped joints neutral in PD by default.
        pTarget_[i] = gc_[i];

        // In stopping mode, use damping-only PD to avoid abrupt accelerations.
        if (joint_stopping[i]) {
          pGain_[i] = 0.0;
          dGain_[i] = hold_d_gain[i];
        }
        continue;
      }

      // Hold stopped joints with prescribed gains.
      pGain_[i] = hold_p_gain[i];
      dGain_[i] = hold_d_gain[i];
    }
    raipal_L->setPdGains(pGain_, dGain_);
    raipal_L->setPdTarget(pTarget_, dTarget_);

    Eigen::VectorXd tau = Eigen::VectorXd::Zero(dof_);
    const double sim_t = t * dt;
    for (int i = 0; i < dof_ && i < snap_torque.size(); i++) {
      const double max_tau = snap_torque[i] * (1 - warmup_factor);
      if (sim_t < window || max_tau <= 0.0) {
        tau[i] = 0.0;
        continue;
      }

      if (joint_stopped[i]) {
        tau[i] = 0.0;
        continue;
      }

      if (joint_stopping[i]) {
        tau[i] = 0.0;
        continue;
      }

      if (std::abs(gv_[i]) > vel_stop_thresh) {
        tau[i] = -std::copysign(max_tau, gv_[i]);
      } else {
        tau[i] = 0.0;
      }

      tau[i] = std::clamp(tau[i], -max_tau, max_tau);
    }
    raipal_L->setGeneralizedForce(tau);

    server.integrateWorldThreadSafe();
    raipal_L->getState(gc_, gv_);
    executed_steps = t + 1;

    // analyze step here
    if (t % 100 == 0) {
      std::cout << "STEP " << t << "/" << TOTAL_STEPS << std::endl;
    }

    // joint-limit safety check
    const int nCheck = std::min<int>(jointLimits.size(), dof_);
    for (int i = 0; i < nCheck; i++) {
      const double lower = jointLimits[i][0];
      const double upper = jointLimits[i][1];
      const bool at_or_outside_limit = (gc_[i] <= lower + limit_eps) || (gc_[i] >= upper - limit_eps);
      if (at_or_outside_limit && std::abs(gv_[i]) > vel_stop_thresh) {
        trajectory_invalid = true;
        invalid_reason = "joint " + std::to_string(i) + " hit limit with nonzero velocity (" + std::to_string(gv_[i]) + ")";
        break;
      }
    }
    if (trajectory_invalid) {
      std::cout << "Trajectory invalid: " << invalid_reason << std::endl;
      break;
    }

    const double current_t = executed_steps * dt;
    if (current_t + 1e-12 >= next_log_t) {
      gc_log.push_back(gc_);
      time_log.push_back(FORWARD ? current_t : -current_t);
      next_log_t += logging_dt;
    }

    if (gv_.cwiseAbs().maxCoeff() < vel_stop_thresh) {
      std::cout << "All joint speeds below " << vel_stop_thresh
                << " at step " << executed_steps << ". Finishing early." << std::endl;
      break;
    }

    gv_prev = gv_;
    warmup_factor *= warmup_alpha;
  }

  // Save logged trajectory
  const std::string csv_path = "raipal_throw_gc_log.csv";
  std::ofstream csv(csv_path);
  if (csv.is_open()) {
    csv << "time";
    for (int i = 0; i < dof_; i++) csv << ",q" << i;
    csv << "\n";

    if (FORWARD) {
      for (size_t k = 0; k < gc_log.size(); k++) {
        csv << time_log[k];
        for (int i = 0; i < dof_; i++) csv << "," << gc_log[k][i];
        csv << "\n";
      }
    } else {
      for (size_t k = gc_log.size(); k-- > 0;) {
        csv << time_log[k];
        for (int i = 0; i < dof_; i++) csv << "," << gc_log[k][i];
        csv << "\n";
      }
    }
    csv.close();
    std::cout << "Saved trajectory snapshots to: " << csv_path
              << " (" << gc_log.size() << " rows)" << std::endl;
  } else {
    std::cout << "Failed to write CSV: " << csv_path << std::endl;
  }

  server.killServer();

  std::cout << "SIMULATION COMPLETE"
            << (trajectory_invalid ? " (INVALID TRAJECTORY)" : "")
            << std::endl;

  return 0;
}
