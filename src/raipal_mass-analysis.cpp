#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"

size_t TOTAL_STEPS = 20000;

Eigen::VectorXd condense(const Eigen::VectorXd& v9) {
  Eigen::VectorXd v7(7);
  v7 << v9.head(3), v9.tail(4);
  return v7;
}

size_t condenseIdx(size_t idx) {
  if (idx <= 2) return idx;
  else return idx - 2; // skip the spring coordinate
}

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);
  
  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  // auto raipal_R = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_stub-0_R.urdf");
  // auto raipal_L = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_stub-0_L.urdf");

  // raipal_R->setName("raipal_R");
  // raipal_L->setName("raipal_L");

  auto raipal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/urdf/raipal_stub-0_L.urdf");
  raipal->setName("raipal");

  auto raipal_alt = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/raipal/alternative/raipal-ALT_stub-0_R.urdf");
  raipal_alt->setName("raipal_alt");

  std::cout << "raipal has " << raipal->getGeneralizedCoordinateDim() << " generalized coordinates" << std::endl;
  std::cout << "raipal_alt has " << raipal_alt->getGeneralizedCoordinateDim() << " generalized coordinates" << std::endl;
  // raipal -> setComputeInverseDynamics(true);
  std::cout << "robots were loaded!" << std::endl;

  std::vector<Eigen::VectorXd> poses;
  std::vector<std::vector<size_t>> indices;

  Eigen::VectorXd pose(9);
  pose << 
    0.0,-0.52359,0.0,
    0.0,0.0,
    0.0,
    0.0,0.0,0.0
  ;
  poses.push_back(pose);
  indices.push_back({0,1});  

  pose <<
    0.0,0.0,0.0,
    0.0,0.0,
    1.57,
    0.0,0.0,0.0
  ;
  poses.push_back(pose);
  indices.push_back({2,5});

  pose <<
    0.0,0.0,0.0,
    0.0,0.0,
    0.0,
    0.0,1.57,0.0
  ;
  poses.push_back(pose);
  indices.push_back({6});  

  pose <<
    0.0,0.0,0.0,
    0.0,0.0,
    0.0,
    0.0,0.0,0.0
  ;
  poses.push_back(pose);
  indices.push_back({7,8});

  std::cout << "Inertia Poses Created" << std::endl;

  server.launchServer();
  server.focusOn(raipal);
  for (int sec=5; sec>0; sec--){
    std::cout << "Starting in [" << sec << "]..." << std::endl;
    raisim::USLEEP(1000000);
  }
  std::cout << "START!" << std::endl;

  for(size_t i=0; i<poses.size(); i++){
    raipal    ->setState(poses[i], Eigen::VectorXd::Zero(9));
    raipal_alt->setState(condense(poses[i]), Eigen::VectorXd::Zero(7));
    std::cout << "[Pose " << i << "] " << poses[i].transpose() << std::endl;
    auto mass_diag     = raipal->getMassMatrix().e().diagonal();
    auto mass_diag_alt = raipal_alt->getMassMatrix().e().diagonal();
    std::cout << "diag(M)       : " << mass_diag.transpose() << std::endl;
    std::cout << "diag(M) (alt) : " << mass_diag_alt.transpose() << std::endl;
    for(auto joint_idx : indices[i]){
      double inertia = mass_diag[joint_idx];
      double inertia_alt = mass_diag_alt[condenseIdx(joint_idx)];

      if(joint_idx == 5){ // worst case elbow reflected inertia
        inertia = inertia / 0.53 + mass_diag[3];
      }

      std::cout << "  Joint " << joint_idx << "\n" <<
        "    inertia: " << inertia << " / alt: " << inertia_alt << " (" << 100.0*inertia_alt/inertia << "%)" << "\n" <<
        "    max acceleration: " << (joint_idx <= 5 ? 183 : 12.6) / inertia << " rad/s^2" << "\n";
    }
    raisim::USLEEP(1000000);
  }

  server.killServer();

  std::cout<<"SIMULATION COMPLETE"<<std::endl;

  return 0;
}
