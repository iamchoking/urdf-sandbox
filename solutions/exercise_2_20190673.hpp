#pragma once
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
       w[2],   0 ,-w[0],
      -w[1], w[0],   0 ;
  // just one typo here set me back 3 hours... :(
  return R;
}

class Trans {
public:
  char typ; //'f' fixed, 'r' revolute (to be added: 'p' prismatic)
  Eigen::Vector3d originPos;
  Eigen::Matrix3d originRot;
  Eigen::Vector3d axis;
  int gcIdx; // index of gc to use in kinematics
  int gvIdx; // index of gv to use in diff. kinematics (usually gcIdx-1)

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

  // "set" functions for streamlining hard-coding
  void setXyz (double x,double y,double z){originPos << x,y,z;}
  void setRpy (double r,double p,double y){originRot = rpyToRot(r,p,y);}
  void setAxis(double ax,double ay,double az){axis << ax,ay,az;}
  void setTyp (char newTyp = 'f'){typ=newTyp;}
  void setIdx (int gcIndex = -1,int gvIndex = -1){gcIdx = gcIndex;gvIdx = gvIndex;}

  [[maybe_unused]] void setProfile(double* xyz,double* rpy,char newTyp = 'f',double* axis = nullptr,int gcIndex = -1,int gvIndex = -1){
    setXyz (xyz[1] ,xyz[2] ,xyz[3] );
    setRpy (rpy[1] ,rpy[2] ,rpy[3] );
    setAxis(axis[1],axis[2],axis[3]);
    setTyp(newTyp);
    setIdx(gcIndex,gvIndex);
  }

  void setProfile(double x,double y,double z, double R,double P,double Y,char newTyp = 'f',double ax=0,double ay=0,double az=0,int gcIndex = -1,int gvIndex = -1){
    setXyz (x,y,z);
    setRpy (R,P,Y);
    setAxis(ax,ay,az);
    setTyp(newTyp);
    setIdx(gcIndex,gvIndex);
    if(typ == 'f'){initFixed();}
  }

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
      calcKin = true;
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
    calcKin = true;
    return worldT.evalTrans(gc);
  }

  Eigen::MatrixXd calculateVJ(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& wree){
    // wree: the position of end effector (point of question) in world coordinates.
    if (!calcKin){calculateKinematics(gc);}

    Eigen::MatrixXd vJ;
    if(typ == 'b'){
      vJ.resize(3,gv.size());
      vJ.setZero();
      vJ.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
      vJ.block<3,3>(0,3) = -skew3d(wree-worldT.originPos);
      // std::cout << "computed arm : " << (wree-worldT.originPos).transpose() << std::endl;
      // std::cout << "BASE vJ" << std::endl << vJ << std::endl;
    }
    else if(typ == 'a' || typ == 'e'){
      vJ = parent->calculateVJ(gc,gv,wree);
      if (worldT.gvIdx >= 0){
        if (worldT.typ == 'r'){
          vJ.block<3,1>(0,worldT.gvIdx) = skew3d(worldT.originRot*worldT.axis)*(wree-worldT.originPos);
          // std::cout << "computed axis: " << (worldT.originRot*worldT.axis).transpose() << std::endl;
          // std::cout << "computed arm : " << (wree-worldT.originPos).transpose() << std::endl;
          // std::cout << "computed block: " << vJ.block<3,1>(0,worldT.gvIdx).transpose() << std::endl;
        } //TODO: prismatic jacobian
      }

      // vJ.block<3,1>(0,worldT.gvIdx) = skew3d(worldT.axis)*(wree-worldT.originPos);
    }

    return vJ;
  }

  Eigen::MatrixXd calculateWJ(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& wree){
    // wree: the position of end effector (point of question) in world coordinates.
    if (!calcKin){calculateKinematics(gc);}

    Eigen::MatrixXd wJ;
    if(typ == 'b'){
      wJ.resize(3,gv.size());
      wJ.setZero();
      // wJ.block<3,3>(0,0) = Eigen::Matrix3d::Zero()); //redundant
      wJ.block<3,3>(0,3) = Eigen::Matrix3d::Identity();
    }
    else if(typ == 'a' || typ == 'e'){
      wJ = parent->calculateWJ(gc,gv,wree);
      if (worldT.gvIdx >= 0){
        if (worldT.typ == 'r'){
          wJ.block<3,1>(0,worldT.gvIdx) = worldT.originRot*worldT.axis;
        } //TODO: prismatic jacobian
      }

      // wJ.block<3,1>(0,worldT.gvIdx) = skew3d(worldT.axis)*(wree-worldT.originPos);
    }

    return wJ;
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
  // order: x y z   R P Y  (type = f) (ax ay az) (gcIdx gvIdx)
  tempT.setProfile(-0.2999,0.104,0.0,  -2.61799387799,0,-3.14159265359);
  lhHip -> addTrans(tempT);
  // LH_HAA

  // <<LH_HAA>> (revolute) <axis xyz="-1 0 0"/> (gc[13])
  tempT.setProfile(0.0,0.0,0.0,  0.0,0.0,0.0, 'r', -1,0.0,0.0, 13,12);
  lhHip -> addTrans(tempT);
  // LH_HIP

  // "THIGH"
  auto lhThi   = new Link("LH_THIGH",'a',lhHip);
  // LH_HIP
  // <LH_HIP_LH_hip_fixed> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="0 0 0"/>
  tempT.setProfile(0.0,0.0,0.0,  -2.61799387799,0.0,-3.14159265359);
  lhThi -> addTrans(tempT);
  // LH_hip_fixed
  // <LH_hip_fixed_LH_HFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0599 0.08381 0.0"/>
  tempT.setProfile(-0.0599,0.08381,0.0,  0.0,0.0,1.57079632679);
  lhThi -> addTrans(tempT);
  // LH_HFE
  // <<LH_HFE>> (revolute) <axis xyz="1 0 0"/> (gc[14])
  tempT.setProfile(0.0,0.0,0.0,  0.0,0.0,0.0, 'r', 1,0.0,0.0, 14,13);
  lhThi -> addTrans(tempT);

  // "SHANK"
  auto lhSha = new Link("LH_SHANK",'a',lhThi);
  // LH_THIGH
  // <LH_THIGH_LH_thigh_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0.0,0.0,0.0,  0.0,0.0,-1.57079632679);
  lhSha -> addTrans(tempT);
  // LH_thigh_fixed
  // <LH_thigh_fixed_LH_KFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0 0.1003 -0.285"/>
  tempT.setProfile(0.0,0.1003,-0.285,  0.0,0.0,1.57079632679);
  lhSha -> addTrans(tempT);
  // LH_KFE
  // <<LH_KFE>> (revolute) <axis xyz="1 0 0"/> (gc[15])
  tempT.setProfile(0.0,0.0,0.0,  0.0,0.0,0.0, 'r', 1,0.0,0.0, 15,14);
  lhSha -> addTrans(tempT);

  // "FOOT"
  auto lhFoot = new Link("LH_FOOT",'e',lhSha);
  // LH_SHANK
  // <LH_SHANK_LH_shank_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0.0,0.0,0.0,  0.0,0.0,-1.57079632679);
  lhFoot -> addTrans(tempT);
  // LH_shank_fixed
  // <LH_shank_fixed_LH_FOOT> (fixed) <origin rpy="0 0 0" xyz="-0.08795 0.01305 -0.33797"/>
  tempT.setProfile(-0.08795,0.01305,-0.33797,  0.0,0.0,0.0);
  lhFoot -> addTrans(tempT);
  // LH_FOOT <-- (objective joint origin)

  // return base;
  // return lhHip;
  return lhFoot;
}

/// do not change the name of the method
inline Eigen::Vector3d getFootLinearVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Link* L = setupLinks();
    Eigen::Vector3d wree = L->calculateKinematics(gc)->originPos;
    Eigen::MatrixXd vJ = L->calculateVJ(gc,gv,wree);

    return vJ*gv; /// replace this
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Link* L = setupLinks();
    Eigen::Vector3d wree = L->calculateKinematics(gc)->originPos;
    Eigen::MatrixXd wJ = L->calculateWJ(gc,gv,wree);

    return wJ*gv; /// replace this
}

// hpp end