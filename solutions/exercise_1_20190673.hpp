//
// Created by Jemin Hwangbo on 2022/03/17.
//

#ifndef ME553_2022_SOLUTIONS_EXERCISE1_20190673_HPP_
#define ME553_2022_SOLUTIONS_EXERCISE1_20190673_HPP_

#include <Eigen/Core>
#include <utility>
#include <vector>

using namespace std;
using namespace Eigen;

// Kinematic chain (of interest):
// (world)
// (implicit floating joint) (gc[0]~gc[6])

// "BASE"
// base
// <base_LH_HAA> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="-0.2999 0.104 0.0"/>

// "LH_BASE"
// LH_HAA
// <<LH_HAA>> (revolute) <axis xyz="-1 0 0"/> (gc[13])

// "HIP"
// LH_HIP
// <LH_HIP_LH_hip_fixed> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="0 0 0"/>
// LH_hip_fixed
// <LH_hip_fixed_LH_HFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0599 0.08381 0.0"/>
// LH_HFE
// <<LH_HFE>> (revolute) <axis xyz="1 0 0"/> (gc[14])

// "THIGH"
// LH_THIGH
// <LH_THIGH_LH_thigh_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
// LH_thigh_fixed
// <LH_thigh_fixed_LH_KFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0 0.1003 -0.285"/>
// LH_KFE
// <<LH_KFE>> (revolute) <axis xyz="1 0 0"/> (gc[15])

// "SHANK"
// LH_SHANK
// <LH_SHANK_LH_shank_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
// LH_shank_fixed
// <LH_shank_fixed_LH_FOOT> (fixed) <origin rpy="0 0 0" xyz="-0.08795 0.01305 -0.33797"/>
// LH_FOOT <-- (objective joint origin)

Eigen::Matrix3d quatToRot(Eigen::Vector4d q){
  // ***we assume q is normalized and valid (q.norm() > 0)
  // from formula
  Eigen::Matrix3d R;
  R << 1 - 2 * (q(2)*q(2) + q(3)*q(3)), 2 * (q(1)*q(2) - q(3)*q(0))    , 2 * (q(1)*q(3) + q(2)*q(0))    ,
    2 * (q(1)*q(2) + q(3)*q(0))    , 1 - 2 * (q(1)*q(1) + q(3)*q(3)), 2 * (q(2)*q(3) - q(1)*q(0))    ,
    2 * (q(1)*q(3) - q(2)*q(0))    , 2 * (q(2)*q(3) + q(1)*q(0))    , 1 - 2 * (q(1)*q(1) + q(2)*q(2));
  return R;
}

Eigen::Matrix3d rpyToRot(double r,double p,double y){
  Eigen::Matrix3d Rx, Ry, Rz;

  // Individual rotation matrices
  //roll: x axis
  Rx << 1, 0, 0,
    0, cos(r), -sin(r),
    0, sin(r), cos(r);
  // pitch: y axis
  Ry << cos(p), 0, sin(p),
    0, 1, 0,
    -sin(p), 0, cos(p);
  // yaw: z axis
  Rz << cos(y), -sin(y), 0,
    sin(y), cos(y), 0,
    0, 0, 1;

  // Combine rotations (roll, then pitch, then yaw)
  return Rz * Ry * Rx; // Passive Rotation + urdf convention
}

Eigen::Matrix3d rpyToRot(Eigen::Vector3d rpy){
  return rpyToRot(rpy(0),rpy(1),rpy(2));
}

class Trans {
public:
  char typ; //'f' fixed, 'r' revolute (to be added)
  Eigen::Vector3d originPos;
  Eigen::Matrix<double, 3, 3> originRot;
  Eigen::Vector3d axis;
  int gcIdx; // index of gc to use in kinematics

  Trans(const char t, Eigen::Vector3d xyz, Eigen::Vector3d rpy, Eigen::Vector3d ax, int gcIndex = -1) {
    typ = t;
    originPos = std::move(xyz);
    originRot = rpyToRot(std::move(rpy));
    axis = std::move(ax);
    gcIdx = gcIndex;

    if(typ == 'f'){
      axis << 0,0,0;
      gcIdx = -1;
    }
  }

  // evaluate the full transform for the given generalized coord.
  void evalTrans(const Eigen::VectorXd &gc,Eigen::Vector3d &r,Eigen::Matrix<double,3,3> &R){
    // urdf convention does the "origin move" first, then the rotation.
    if (typ == 'f'){r = originPos;R = originRot;return;}
    else if (typ == 'r'){
      r = originPos;

      Eigen::Vector4d quat; // rotation about an axis can be expressed as quat.
      quat << cos(gc[gcIdx]/2),sin(gc[gcIdx]/2)*(axis/axis.norm());
      R = originRot * quatToRot(quat);
      return;
    }
  }
};

class Link{
public:
  std::string name; //name
  char typ; //'b': base link, 'a': articulated link
  Link* parent; //pointer to parent link
  vector<Trans*> transforms;
  bool calc; // 'true': Kinematics Calculated 'false': Kinematics not calculated yet
  Eigen::Vector3d r; // location of the final transformation (expressed in world frame)
  Eigen::Matrix3d R; // orientation of the final transformation (expressed in world frame)

  Link(const std::string& n,const char t,Link* p){
    name = n;
    typ = t;
    parent = p;
    calc = false;
    r = Eigen::Vector3d::Zero();
    R = Eigen::Matrix3d::Identity();
  }

  Link(const std::string& n,const char t): Link(n,t,nullptr){}

  ~Link() {
    // Delete dynamically allocated memory within the Link object
    for (auto & transform : transforms) {
      delete transform;
    }
    transforms.clear();
  }

  void addTrans(Trans* t){
    transforms.push_back(t);
  }

  void calculateKinematics(const Eigen::VectorXd& gc){
    if(typ == 'b'){ //base case (no parent link, get r / R from base-pos)
      // std::cout<< "base floating body: " << gc.transpose() << std::endl;
      r = gc.segment(0,3);
      R = quatToRot(gc.segment(3,4));
    }
    else if(!parent->calc){
      parent->calculateKinematics(gc);
      r = parent->r;
      R = parent->R;
    } //recursive calls for articulated links

    // traverse and calculate through transforms
    Eigen::Vector3d  tempr;
    Eigen::Matrix<double,3,3> tempR;

    for (auto trans:transforms){
      trans ->evalTrans(gc,tempr,tempR);
      r = r+R*tempr;
      R = R*tempR;
    }
  }

};

Link* solveLinks (const Eigen::VectorXd& gc) {

  // all the "hard-coding" is done in this function.

  Eigen::Vector3d xyz_;
  Eigen::Vector3d rpy_;
  Eigen::Vector3d ax_;

  // note: though this can be further simplified, this is left as-is to mimic the workflow in [anymal.urdf].
  // "BASE"
  // base
  auto base    = new Link("BASE",'b');
  // <base_LH_HAA> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="-0.2999 0.104 0.0"/>
  xyz_ << -0.2999,0.104,0.0;
  rpy_ << -2.61799387799,0,-3.14159265359;
  // rpy_ << -2.61799387799,1.4,-0.69;
  ax_  << 0,0,0;
  base -> addTrans(new Trans('f',xyz_,rpy_,ax_,-1));

  // "LH_BASE"
  auto lhBase  = new Link("LH_BASE",'a',base);
  // LH_HAA
  // <<LH_HAA>> (revolute) <axis xyz="-1 0 0"/> (gc[13])
  xyz_ << 0,0.0;
  rpy_ << 0,0,0;
  // ax_  << -1,1,3;
  ax_  << -1,0,0;
  lhBase -> addTrans(new Trans('r',xyz_,rpy_,ax_,13));

  // "HIP"
  auto lhHip   = new Link("LH_HIP",'a',lhBase);
  // LH_HIP
  // <LH_HIP_LH_hip_fixed> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="0 0 0"/>
  xyz_ << 0,0,0;
  rpy_ << -2.61799387799,0,-3.14159265359;
  ax_  << 0,0,0;
  lhHip -> addTrans(new Trans('f',xyz_,rpy_,ax_,-1));
  // LH_hip_fixed
  // <LH_hip_fixed_LH_HFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0599 0.08381 0.0"/>
  xyz_ << -0.0599,0.08381,0.0;
  rpy_ << 0,0,1.57079632679;
  ax_  << 0,0,0;
  lhHip -> addTrans(new Trans('f',xyz_,rpy_,ax_,-1));
  // LH_HFE
  // <<LH_HFE>> (revolute) <axis xyz="1 0 0"/> (gc[14])
  xyz_ << 0,0,0;
  rpy_ << 0,0,0;
  ax_  << 1,0,0;
  lhHip -> addTrans(new Trans('r',xyz_,rpy_,ax_,14));

  // "THIGH"
  auto lhThigh = new Link("LH_THIGH",'a',lhHip);
  // LH_THIGH
  // <LH_THIGH_LH_thigh_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  xyz_ << 0,0,0;
  rpy_ << 0,0,-1.57079632679;
  ax_  << 0,0,0;
  lhThigh -> addTrans(new Trans('f',xyz_,rpy_,ax_,-1));
  // LH_thigh_fixed
  // <LH_thigh_fixed_LH_KFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0 0.1003 -0.285"/>
  xyz_ << 0,0.1003,-0.285;
  rpy_ << 0,0,1.57079632679;
  ax_  << 0,0,0;
  lhThigh -> addTrans(new Trans('f',xyz_,rpy_,ax_,-1));
  // LH_KFE
  // <<LH_KFE>> (revolute) <axis xyz="1 0 0"/> (gc[15])
  xyz_ << 0,0,0;
  rpy_ << 0,0,0;
  ax_  << 1,0,0;
  lhThigh -> addTrans(new Trans('r',xyz_,rpy_,ax_,15));

  // "SHANK"
  auto lhShank = new Link("LH_SHANK",'a',lhThigh);
  // LH_SHANK
  // <LH_SHANK_LH_shank_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  xyz_ << 0,0,0;
  rpy_ << 0,0,-1.57079632679;
  ax_  << 0,0,0;
  lhShank -> addTrans(new Trans('f',xyz_,rpy_,ax_,-1));
  // LH_shank_fixed
  // <LH_shank_fixed_LH_FOOT> (fixed) <origin rpy="0 0 0" xyz="-0.08795 0.01305 -0.33797"/>
  xyz_ << -0.08795,0.01305,-0.33797;
  rpy_ << 0,0,0;
  ax_  << 0,0,0;
  lhShank -> addTrans(new Trans('f',xyz_,rpy_,ax_,-1));
  // LH_FOOT <-- (objective joint origin)

  lhShank -> calculateKinematics(gc);

  return lhShank;
}

Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {
  Link* ee = solveLinks(gc);

  // auto final_pos = ee -> r;
  // auto final_ori = ee -> R;
  // std::cout << "CALCULATION RESULTS" << std::endl;
  // std::cout << "CALC_POS:" << std::endl << final_pos.transpose() << std::endl;
  // std::cout << "CALC_ORI:" << std::endl << final_ori << std::endl;

  return ee -> r;
}

#endif
