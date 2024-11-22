//
// Created by Jemin Hwangbo on 2022/03/17.
//

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

// #include "exercise1_20190673.hpp"
#include "raisim/RaisimServer.hpp"
#include "raisim/World.hpp"
#include "random_coordinates.hpp"
#include <iostream>

// #include "exercise2_20190673.hpp"
// excercise1_20190673.hpp
#include <Eigen/Core>
#include <utility>
#include <vector>
#include <iostream>

// SETUP
  // Kinematic chain (of interest):
  // (world)
  // (implicit floating joint) (gc[0]~gc[6])
  //
  // "BASE"
  // base
  // <base_LH_HAA> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="-0.2999 0.104 0.0"/>
  // 
  // "LH_BASE"
  // LH_HAA
  // <<LH_HAA>> (revolute) <axis xyz="-1 0 0"/> (gc[13])
  //
  // "HIP"
  // LH_HIP
  // <LH_HIP_LH_hip_fixed> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="0 0 0"/>
  // LH_hip_fixed
  // <LH_hip_fixed_LH_HFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0599 0.08381 0.0"/>
  // LH_HFE
  // <<LH_HFE>> (revolute) <axis xyz="1 0 0"/> (gc[14])
  //
  // "THIGH"
  // LH_THIGH
  // <LH_THIGH_LH_thigh_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  // LH_thigh_fixed
  // <LH_thigh_fixed_LH_KFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0 0.1003 -0.285"/>
  // LH_KFE
  // <<LH_KFE>> (revolute) <axis xyz="1 0 0"/> (gc[15])
  //
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

Eigen::Matrix3d skew3d(const Eigen::Vector3d& w){
  Eigen::Matrix3d R;
  R <<   0 ,-w[2], w[1],
       w[2],   0 , w[0],
      -w[1], w[0],   0 ;
  return R;
}

class Trans {
public:
  char typ; //'f' fixed, 'r' revolute (to be added) ('p' prismatic)
  Eigen::Vector3d originPos;
  Eigen::Matrix3d originRot;
  Eigen::Vector3d axis;
  int gcIdx; // index of gc to use in kinematics
  int gvIdx; // index of gv to use in diff. kinematics

  Trans(const char t, Eigen::Vector3d xyz, Eigen::Vector3d rpy, Eigen::Vector3d ax, int gcIndex = -1, int gvIndex = -1) {
    typ = t;
    originPos = std::move(xyz);
    originRot = rpyToRot(std::move(rpy));
    axis = std::move(ax);
    gcIdx = gcIndex;
    gvIdx = gvIndex;

    if(typ == 'f'){
      initFixed();
    }
  }

  Trans(Eigen::Vector3d r, Eigen::Matrix3d R){ //simple creation (result of evalTrans)
    initFixed();
    originPos = r;
    originRot = R;
  }

  Trans(){ // trivial trans
    initFixed();
    originPos << 0,0,0;
    originRot = Eigen::Matrix3d::Identity();
  }

  void setXyz (double x,double y,double z){originPos << x,y,z;}
  void setRpy (double r,double p,double y){originRot = rpyToRot(r,p,y);}
  void setAxis(double ax,double ay,double az){axis << ax,ay,az;}
  void setTyp (char newTyp){typ=newTyp;}
  void setIdx (int gcIndex = -1,int gvIndex = -1){gcIdx = gcIndex;gvIdx = gvIndex;}

  void initFixed(){ // pattern for fixed transformations
    typ = 'f';
    gcIdx = -1;
    gvIdx = -1;
    axis << 0,0,0;
  }

  void attachTrans(const Trans& newT){ //attach transformation to the end of this one (modifier)
    if(typ != 'f'){ //error: middle actuation
      throw std::invalid_argument("middle-actuation detected!");
    }
    originPos = originPos + originRot*newT.originPos;
    originRot = originRot*newT.originRot;
    if(newT.typ != 'f'){
      if(typ != 'f'){ //error: double actuation! (this is actually impossible)
        throw std::invalid_argument("double-actuation detected!");
      }
      typ = newT.typ;
      axis = newT.axis;
      gcIdx = newT.gcIdx;
      gvIdx = newT.gvIdx;
    }
  }

  // evaluate the full transform for the given generalized coord.
  // returns a fixed transformation
  Trans* evalTrans(const Eigen::VectorXd &gc){
    // urdf convention does the "origin move" first, then the actuation wrt axis.
    Eigen::Vector3d newPos = originPos;
    Eigen::Matrix3d newRot = originRot;


    if(typ == 'r'){
      Eigen::Vector4d quat; // rotation about an axis can be expressed as quat.
      quat << cos(gc[gcIdx]/2),sin(gc[gcIdx]/2)*(axis/axis.norm());
      newRot = originRot * quatToRot(quat);

    }
    else if(typ == 'p'){ //UNTESTED
      newPos = originPos + gc[gcIdx] * (originRot*axis);
    }
    return new Trans(newPos,newRot);
  }
};

//TODO: throw warning/error if an articulated link does not have a gcindex

class Link{
public:
  std::string name; //name
  char typ; //'b': base link, 'a': articulated link, 'e': end-effector (no actuation)
  Link* parent; //pointer to parent link
  // ex2: !! changed this from "final" to "root" references!!

  Trans bodyT; // full transform from (parent^2<>parent) to (parent<>this) (r,R,p,gcIdx,gvIdx) ("assume parent is fixed to world")
  // transform order (urdf convention) r -> R -> p
    // translate by r (expressed in (p^2<>p) frame)
    // rotate by R (expressed in (p^2<>p) frame)
    // actuate wrt p by (gc[gcIdx]) (expressed in new frame(r -> R))

  // state variable (another transform)
  Trans worldT;
  bool calcKin; // 'true': Kinematics Calculated. 'false': Kinematics not calculated yet
  //*: non-constant


  Link(const std::string& n,const char linkTyp,Link* p){
    name = n;
    typ = linkTyp;
    parent = p;
    calcKin = false;

    bodyT  = *(new Trans());
    worldT = *(new Trans());
  }

  Link(const std::string& n,const char t): Link(n,t,nullptr){}

  void addTrans(const Trans& newT){
    // add an additional transformation at the end of the link
    // only modifies constant properties
    bodyT.attachTrans(newT);
  }

  Trans* calculateKinematics(const Eigen::VectorXd& gc){
    if(typ == 'b'){ //base case (no parent link, get wr / wR from base-pos)
      // std::cout<< "base floating body: " << gc.transpose() << std::endl;
      worldT.originPos = gc.segment(0,3);
      worldT.originRot = quatToRot(gc.segment(3,4));
      return worldT.evalTrans(gc);
    }
    else if(!parent->calcKin){
      worldT = *(parent->calculateKinematics(gc));
    } //recursive calls for articulated links
    // std::cout << "parent pos" << worldT.originPos << std::endl;
    // std::cout << "parent ori" << worldT.originRot << std::endl;
    // std::cout << "parent typ" << worldT.typ << std::endl;

    worldT.attachTrans(bodyT); // just attach my own transform!
    // TODO warning for un-actuated worldT's (in the case of 'a' links)

    // although each links have to keep axis information (for J), the output for kinematics must be evaluated!
    return worldT.evalTrans(gc);
  }

  void calculateJacobian(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv){

  }

};

Link* setupLinks() {

  // all the "hard-coding" is done in this function.
  Trans tempT = *(new Trans());

  // note: though this can be further simplified, this is left as-is to mimic the workflow in [anymal.urdf].
  // "BASE"
  auto base    = new Link("BASE",'b');
  // base

  // "HIP"
  auto lhHip = new Link("LH_HIP",'a',base);
  // base
  // <base_LH_HAA> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="-0.2999 0.104 0.0"/>
  tempT.setXyz (-0.2999,0.104,0.0);
  tempT.setRpy (-2.61799387799,0,-3.14159265359);
  tempT.setAxis(0,0,0);
  tempT.setTyp ('f');
  lhHip -> addTrans(tempT);
  // LH_HAA

  // <<LH_HAA>> (revolute) <axis xyz="-1 0 0"/> (gc[13])
  tempT.setXyz (0,0,0);
  tempT.setRpy (0,0,0);
  tempT.setAxis(-1,0,0);
  tempT.setTyp ('r');
  tempT.setIdx (13,12);
  lhHip -> addTrans(tempT);
  // LH_HIP

  // "THIGH"
  auto lhThi   = new Link("LH_THIGH",'a',lhHip);
  // LH_HIP
  // <LH_HIP_LH_hip_fixed> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="0 0 0"/>
  tempT.setXyz (0,0,0);
  tempT.setRpy (-2.61799387799,0,-3.14159265359);
  tempT.setAxis(0,0,0);
  tempT.setTyp ('f');
  lhThi -> addTrans(tempT);
  // LH_hip_fixed
  // <LH_hip_fixed_LH_HFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0599 0.08381 0.0"/>
  tempT.setXyz (-0.0599,0.08381,0.0);
  tempT.setRpy (0,0,1.57079632679);
  tempT.setAxis(0,0,0);
  tempT.setTyp ('f');
  lhThi -> addTrans(tempT);
  // LH_HFE
  // <<LH_HFE>> (revolute) <axis xyz="1 0 0"/> (gc[14])
  tempT.setXyz (0,0,0);
  tempT.setRpy (0,0,0);
  tempT.setAxis(1,0,0);
  tempT.setTyp ('r');
  tempT.setIdx (14,13);
  lhThi -> addTrans(tempT);

  // "SHANK"
  auto lhSha = new Link("LH_SHANK",'a',lhThi);
  // LH_THIGH
  // <LH_THIGH_LH_thigh_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  tempT.setXyz (0,0,0);
  tempT.setRpy (0,0,-1.57079632679);
  tempT.setAxis(0,0,0);
  tempT.setTyp ('f');
  lhSha -> addTrans(tempT);
  // LH_thigh_fixed
  // <LH_thigh_fixed_LH_KFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0 0.1003 -0.285"/>
  tempT.setXyz (0,0.1003,-0.285);
  tempT.setRpy (0,0,1.57079632679);
  tempT.setAxis(0,0,0);
  tempT.setTyp ('f');
  lhSha -> addTrans(tempT);
  // LH_KFE
  // <<LH_KFE>> (revolute) <axis xyz="1 0 0"/> (gc[15])
  tempT.setXyz (0,0,0);
  tempT.setRpy (0,0,0);
  tempT.setAxis(1,0,0);
  tempT.setTyp ('r');
  tempT.setIdx (15,14);
  lhSha -> addTrans(tempT);

  // "FOOT"
  auto lhFoot = new Link("LH_FOOT",'e',lhSha);
  // LH_SHANK
  // <LH_SHANK_LH_shank_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  tempT.setXyz (0,0,0);
  tempT.setRpy (0,0,-1.57079632679);
  tempT.setAxis(0,0,0);
  tempT.setTyp ('f');
  lhFoot -> addTrans(tempT);
  // LH_shank_fixed
  // <LH_shank_fixed_LH_FOOT> (fixed) <origin rpy="0 0 0" xyz="-0.08795 0.01305 -0.33797"/>
  tempT.setXyz (-0.08795,0.01305,-0.33797);
  tempT.setRpy (0,0,0);
  tempT.setAxis(0,0,0);
  tempT.setTyp ('f');
  lhFoot -> addTrans(tempT);
  // LH_FOOT <-- (objective joint origin)

  return lhFoot;
}

int main(int argc, char* argv[]) {
  // create raisim world
  raisim::World world; // physics world 
  raisim::RaisimServer server(&world); // visualization server
  // world.addGround(); //we aren't integrating. we don't really need ground

  // anymal
  auto anymal = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/anymal_c/urdf/anymal_TEST.urdf");
  anymal->setName("anymal");
  server.focusOn(anymal);
  

  // anymal configuration
  // for randomness

  Eigen::VectorXd gc(anymal->getGeneralizedCoordinateDim());

  gc <<
    0.5,  //[ 0]Body Pos z
    0,    //[ 1]Body Pos y
    0.54, //[ 2]Body Pos z
    1.0,  //[ 3]Body Ori w
    0.0,  //[ 4]Body Ori x
    0.0,  //[ 5]Body Ori y
    0.0,  //[ 6]Body Ori z
    -0.03,//[ 7][LH leg] Hip Abduction / Adduction
    0.4,  //[ 8][LH leg] Hip Flexion   / Extension
    -0.8, //[ 9][LH leg] Knee Flexion  / Extension
    -0.03,//[10][RF leg]
    0.4,
    -0.8,
    -0.03, //[13][LH leg] (Left "Hind")
    -0.4,
    0.8,
    -0.03,//[16][RH leg] (Right "Hind")
    -0.4,
    0.8
    ;

  utils::gcRandomize(gc,1);
  // std::srand((unsigned int) time(nullptr));
  // int gcDims = int(anymal->getGeneralizedCoordinateDim());
  // Eigen::VectorXd joint_noise = Eigen::VectorXd::Random(gcDims);
  // Eigen::Vector4d quat_raw = Eigen::Vector4d::Random(4);
  // Eigen::Vector3d quat_xyz = quat_raw.segment(1,3)/quat_raw.segment(1,3).norm();
  // double quat_angle = quat_raw(0)*M_PI;

  // Eigen::Vector4d quat_rand;
  // quat_rand << cos(quat_angle/2),sin(quat_angle/2)*quat_xyz;
  // // std::cout << "random quaternion" << quat_rand.transpose() << " (norm: " << quat_rand.norm() << ")" << std::endl;
  // // std::cout << "quat axis" << quat_xyz.transpose() << " (norm: " << quat_xyz.norm() << ")"  << std::endl;
  // // std::cout << "quat raw" << quat_raw.transpose() << " (norm: " << quat_raw.norm() << ")"  << std::endl;

  // Eigen::VectorXd gcInput = jointNominalConfig+joint_noise;
  // gcInput(3) = quat_rand(0);
  // gcInput(4) = quat_rand(1);
  // gcInput(5) = quat_rand(2);
  // gcInput(6) = quat_rand(3);
  Eigen::VectorXd gcInput = gc;
  // auto gcInput = jointNominalConfig;
  std::cout << "Input gc: " << gcInput.transpose() << std::endl;

  anymal->setGeneralizedCoordinate(gcInput);
  anymal->updateKinematics();

  // debug sphere

  auto debugO = server.addVisualSphere("debug_O", 0.02);
  auto debugX = server.addVisualSphere("debug_X", 0.015);
  auto debugY = server.addVisualSphere("debug_Y", 0.015);
  auto debugZ = server.addVisualSphere("debug_Z", 0.015);

  debugO->setColor(1,1,1,1);
  debugX->setColor(1,0,0,1);
  debugY->setColor(0,1,0,1);
  debugZ->setColor(0,0,1,1);

  // SOLUTION VERIFICATION
  auto L = setupLinks();
  auto result = L->calculateKinematics(gcInput);
  Eigen::Vector3d ex{0.2,0,0};
  Eigen::Vector3d ey{0,0.2,0};
  Eigen::Vector3d ez{0,0,0.2};

  Eigen::Vector3d solO = result->originPos;
  Eigen::Vector3d solX = solO + result->originRot * ex;
  Eigen::Vector3d solY = solO + result->originRot * ey;
  Eigen::Vector3d solZ = solO + result->originRot * ez;

  debugO->setPosition(solO);
  debugX->setPosition(solX);
  debugY->setPosition(solY);
  debugZ->setPosition(solZ);

  auto final_pos = result->originPos;
  auto final_ori = result->originRot;

  std::cout << "CALCULATION RESULTS" << std::endl;
  std::cout << "CALC_POS:" << std::endl << final_pos.transpose() << std::endl;
  std::cout << "CALC_ORI:" << std::endl << final_ori << std::endl;
  delete L; // prevent memory leak


  // solution sphere
  // auto target_name = "base_top_shell";// use "base_top_shell" for base
  // auto target_name = "base_LH_HAA";
  // auto target_name = "LH_HAA";
  auto target_name = "LH_shank_fixed_LH_FOOT";

  auto answerSphere = server.addVisualSphere("answer_sphere", 0.04);
  answerSphere->setColor(0,1,0,1);
  raisim::Vec<3> pos;
  raisim::Mat<3,3> ori{0,0,0,0,0,0,0,0,0};
  anymal->getFramePosition(target_name, pos);
  anymal->getFrameOrientation(target_name,ori);
  answerSphere->setPosition(pos.e());

  std::cout << std::endl;
  // status
  std::cout << "ACTUAL [target: " << target_name << "]" << std::endl;
  std::cout << "POS: " << std::endl << pos.e().transpose() << std::endl;
  std::cout << "ORI :" << std::endl << ori << std::endl;

  auto pos_diff = solO - pos.e();
  std::cout << "POS DIFF. : " << pos_diff.transpose() << "( norm: " << pos_diff.norm() << ")" <<std::endl;

  std::cout << "this was also tested on VSCode!!!!~" << std::endl;

  // base
  std::cout << std::endl;
  std::cout << "ADDITIONAL CALCS" << std::endl;
  raisim::Mat<3,3> ori2{0,0,0,0,0,0,0,0,0};
  raisim::quatToRotMat(gcInput.segment(3,4),ori2);
  raisim::Vec<4> quat;
  anymal->getBaseOrientation(quat);
  std::cout << "Converted Base Orientation :" << std::endl << ori2 << std::endl;
  std::cout << "Base Orientation in quat: " << quat.e().transpose() << std::endl;
  // visualization
  server.launchServer();
  for (int i=0; i<2000000; i++)
    std::this_thread::sleep_for(std::chrono::microseconds(1000));

  server.killServer();
}
