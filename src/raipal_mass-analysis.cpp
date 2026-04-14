#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x
#include "raisim/RaisimServer.hpp"
#include "random_coordinates.hpp"

size_t TOTAL_STEPS = 20000;

bool DO_ACC = true;
bool DO_END = false;

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
  server.launchServer();
  server.focusOn(raipal);

  Eigen::VectorXd pose(9);
  if(DO_ACC){
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
  }

  // endpoint force calculation
  if(DO_END){

    poses.clear();
    indices.clear();
    pose.setZero(9);

    pose << 
      1.5708,-0.52359,0.0,
      0.0,0.0,
      0.0,
      0.0,0.0,0.0
    ;
    poses.push_back(pose);
    indices.push_back({0});  

    pose << 
      0.0,-1.5709,0.0,
      0.0,0.0,
      0.0,
      0.0,0.0,0.0
    ;
    poses.push_back(pose);
    indices.push_back({1});  

    pose <<
      0.0,-1.5709,0.0,
      0.0,0.0,
      1.5708,
      0.0,0.0,0.0
    ;
    poses.push_back(pose);
    indices.push_back({2});

    pose <<
      0.0,0.0,0.0,
      0.0,0.0,
      1.5708,
      0.0,0.0,0.0
    ;
    poses.push_back(pose);
    indices.push_back({5});

    pose <<
      0.0,0.0,0.0,
      0.0,0.0,
      1.57,
      0.0,1.57,0.0
    ;
    poses.push_back(pose);
    indices.push_back({6});

    pose <<
      0.0,0.0,0.0,
      0.0,0.0,
      0.0,
      0.0,1.57,0.0
    ;
    poses.push_back(pose);
    indices.push_back({7});

    pose <<
      0.0,0.0,0.0,
      0.0,0.0,
      0.0,
      0.0,0.0,1.57
    ;
    poses.push_back(pose);
    indices.push_back({8});

    std::cout << "Endpoint Poses Created" << std::endl;
    std::cout << "Gravity: " << world.getGravity().e().transpose() << std::endl;

    auto bodyNames = raipal->getBodyNames();
    for(auto name : bodyNames){
      std::cout << "Body: " << name << std::endl;
    }

    for (int sec=5; sec>0; sec--){
      std::cout << "Starting in [" << sec << "]..." << std::endl;
      raisim::USLEEP(1000000);
    }
    std::cout << "START!" << std::endl;

    for(size_t i=0; i<poses.size(); i++){
      raipal    ->setState(poses[i], Eigen::VectorXd::Zero(9));
      raipal_alt->setState(condense(poses[i]), Eigen::VectorXd::Zero(7));

      std::cout << "[Pose " << i << "] " << poses[i].transpose() << std::endl;

      auto b = raipal->getNonlinearities(world.getGravity()).e();
      std::cout << "Nonlinearities: " << b.transpose() << std::endl;

      raisim::Vec<3> pos;
      raipal->getFramePosition("LE_tip_fixed", pos);

      std::cout << "Tip Position: " << pos.e().transpose() << std::endl;

      for(auto joint_idx : indices[i]){
        double tau = joint_idx <= 5 ? 183 : 12.6;
        if(joint_idx == 5){tau = tau*0.5313;}
        tau -= b[joint_idx];
        
        raisim::Vec<3> pos_B;
        raipal->getPositionInBodyCoordinate(joint_idx+1, pos, pos_B);
        Eigen::Vector3d pos_e;
        pos_e << pos_B[0], pos_B[1], 0.0;
        std::cout << "  pos_B: " << pos_B.e().transpose() << std::endl;
        double r = joint_idx <= 5 ? pos_e.norm() : pos_e.norm() - 0.08;

        double f_end = tau/r;

        std::cout << "  Joint " << joint_idx << "\n" <<
          "    tau: " << tau << "\n" <<
          "    r: " << r << "\n" <<
          "    f_end: " << f_end << "\n";
      }
      raisim::USLEEP(1000000);
    }

  }
  server.killServer();

  std::cout<<"SIMULATION COMPLETE"<<std::endl;

  return 0;
}
