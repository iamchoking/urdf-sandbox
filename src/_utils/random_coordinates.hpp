# pragma once
#include <Eigen/Core>

namespace utils{

void setSeed(unsigned int seed){
  std::srand(seed);
}

void gcRandomize(Eigen::VectorXd& gc,double gamma = 1){
  int dims = gc.size();

  Eigen::VectorXd noise = Eigen::VectorXd::Random(dims);

  Eigen::Vector4d quat_raw = Eigen::Vector4d::Random(4);
  Eigen::Vector3d quat_xyz = quat_raw.segment(1,3)/quat_raw.segment(1,3).norm();
  double quat_angle = quat_raw(0)*M_PI;
  Eigen::Vector4d quat_rand;
  quat_rand << cos(quat_angle/2),sin(quat_angle/2)*quat_xyz;
  // std::cout << "random quaternion" << quat_rand.transpose() << " (norm: " << quat_rand.norm() << ")" << std::endl;
  // std::cout << "quat axis" << quat_xyz.transpose() << " (norm: " << quat_xyz.norm() << ")"  << std::endl;
  // std::cout << "quat raw" << quat_raw.transpose() << " (norm: " << quat_raw.norm() << ")"  << std::endl;

  gc = gc+gamma*noise;
  gc(3) = quat_rand(0);
  gc(4) = quat_rand(1);
  gc(5) = quat_rand(2);
  gc(6) = quat_rand(3);
}

void gvRandomize(Eigen::VectorXd& gv, double gamma = 1){
  int dims = gv.size();

  Eigen::VectorXd noise = Eigen::VectorXd::Random(dims);

  gv = gv+gamma*noise;
}

void eigenRandomize(Eigen::VectorXd& vec,const Eigen::VectorXd& center,const double gamma = 1){
  int dims = center.size();
  vec.setZero(dims);

  Eigen::VectorXd noise = Eigen::VectorXd::Random(dims);

  vec = center+gamma*noise;
}

void eigenRandomize(Eigen::VectorXd& vec, const Eigen::VectorXd& center,const Eigen::VectorXd& gamma){
  int dims = center.size();
  vec.setZero(dims);
  Eigen::VectorXd noise = Eigen::VectorXd::Random(dims).cwiseAbs();
  vec = center+gamma.cwiseProduct(noise);
}

void eigenRandomizeRange(Eigen::VectorXd& vec, const Eigen::VectorXd& lower, const Eigen::VectorXd& upper){
  const Eigen::VectorXd center = ((upper + lower) / 2.0).eval();
  const Eigen::VectorXd gamma  = ((upper - lower) / 2.0).eval();
  return eigenRandomize(vec, center, gamma);
}

void convertJointLimits(Eigen::VectorXd& lower, Eigen::VectorXd& upper, Eigen::VectorXd& range, const std::vector<raisim::Vec<2UL>>& jointLimits){
  int dims = jointLimits.size();
  lower.setZero(dims);
  upper.setZero(dims);
  range.setZero(dims);
  for (size_t i=0; i<dims; i++){
    lower[i] = jointLimits[i][0];
    upper[i] = jointLimits[i][1];
    range[i] = upper[i] - lower[i];
  }
}

void sampleJointPose(Eigen::VectorXd &sample, const std::vector<raisim::Vec<2UL>>& jointLimits, double padRatio = 0.0){
  Eigen::VectorXd lower(jointLimits.size()), upper(jointLimits.size());
  for (size_t i=0; i<9; i++){
    lower[i] = jointLimits[i][0] + padRatio * (jointLimits[i][1] - jointLimits[i][0]);
    upper[i] = jointLimits[i][1] - padRatio * (jointLimits[i][1] - jointLimits[i][0]);
  }
  eigenRandomizeRange(sample, lower, upper);
}

void sampleJointPose(Eigen::VectorXd &sample, const std::vector<raisim::Vec<2UL>>& jointLimits, Eigen::VectorXd padRatio, Eigen::VectorXd center){
  Eigen::VectorXd lower(jointLimits.size()), upper(jointLimits.size());
  for (size_t i=0; i<9; i++){
    if(padRatio[i] < 0 || padRatio[i] > 0.5){
      lower[i] = center[i];
      upper[i] = center[i];
      continue;
    }
    lower[i] = jointLimits[i][0] + padRatio[i] * (jointLimits[i][1] - jointLimits[i][0]);
    upper[i] = jointLimits[i][1] - padRatio[i] * (jointLimits[i][1] - jointLimits[i][0]);
  }
  eigenRandomizeRange(sample, lower, upper);
}

void sampleJointPose(Eigen::VectorXd &sample, const std::vector<raisim::Vec<2UL>>& jointLimits, Eigen::VectorXd padRatio){
  sampleJointPose(sample, jointLimits, padRatio, Eigen::VectorXd::Zero(padRatio.size()));
}

// center + sin(t) * amplitude yieilds a joint sweep of:
// amplitude ratio within [minAmplitudeRatio, maxAmplitudeRatio]
// amplitude is maximum and minimum point within joint limits & pad ratio
void sampleJointSweep(
  Eigen::VectorXd& center, Eigen::VectorXd& amplitude,
  std::vector<raisim::Vec<2UL>> jointLimits,
  Eigen::VectorXd padRatio,
  Eigen::VectorXd minAmplitudeRatio, Eigen::VectorXd maxAmplitudeRatio
){
  if((maxAmplitudeRatio*2 + padRatio*2).maxCoeff() > 1.0){
    std::cerr << "Error: maxAmplitudeRatio*2 + padRatio*2 should be less than or equal to 1.0" << std::endl;
    return;
  }

  Eigen::VectorXd amplitudeRatio;
  eigenRandomizeRange(amplitudeRatio, minAmplitudeRatio, maxAmplitudeRatio);

  Eigen::VectorXd possibleMin,possibleMax, possibleRange;
  convertJointLimits(possibleMin, possibleMax, possibleRange, jointLimits);
  possibleMin += possibleRange.cwiseProduct(padRatio + amplitudeRatio);
  possibleMax -= possibleRange.cwiseProduct(padRatio + amplitudeRatio);

  eigenRandomizeRange(center, possibleMin, possibleMax);
  amplitude = possibleRange.cwiseProduct(amplitudeRatio);

}

// center + sin(t) * amplitude yieilds a joint sweep of:
// amplitude ratio within [minAmplitudeRatio, maxAmplitudeRatio]
// amplitude is maximum and minimum point within joint limits & pad ratio
void sampleJointSweep(
  Eigen::VectorXd& center, Eigen::VectorXd& amplitude,
  std::vector<raisim::Vec<2UL>> jointLimits,
  Eigen::VectorXd padRatio,
  Eigen::VectorXd minAmplitudeRatio
){
  Eigen::VectorXd maxAmplitudeRatio = Eigen::VectorXd::Ones(padRatio.size()) - padRatio*2;
  sampleJointSweep(center, amplitude, jointLimits, padRatio, minAmplitudeRatio, maxAmplitudeRatio);
}

void sampleJointSweep(
  Eigen::VectorXd& center, Eigen::VectorXd& amplitude,
  std::vector<raisim::Vec<2UL>> jointLimits,
  Eigen::VectorXd padRatio
){
  Eigen::VectorXd minAmplitudeRatio = Eigen::VectorXd::Zero(padRatio.size());
  sampleJointSweep(center, amplitude, jointLimits, padRatio, minAmplitudeRatio);
}

} // namespace utils

