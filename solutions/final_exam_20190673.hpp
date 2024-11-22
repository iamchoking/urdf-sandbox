#pragma once

#include <Eigen/Core>
#include <utility>
#include <vector>
#include <iostream>

/// useful tool functions

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

Eigen::Matrix3d hcInertia(double ixx, double ixy, double ixz, double iyy, double iyz, double izz){
  Eigen::Matrix3d I;
  I << ixx,ixy,ixz,
    ixy,iyy,iyz,
    ixz,iyz,izz
    ;
  return I;
}

[[maybe_unused]] Eigen::MatrixXd hc6x6(double m,double ixx, double ixy, double ixz, double iyy, double iyz, double izz){
  Eigen::MatrixXd I;
  I.resize(6,6);
  I.setZero();
  I.block<3,3>(0,0) = m * Eigen::Matrix3d::Identity();
  I.block<3,3>(3,3) = hcInertia(ixx,ixy,ixz,iyy,iyz,izz);

  return I;
}

Eigen::Vector4d rotToQuat(Eigen::Matrix3d rot){ // this is used only for visualization and nowhere else!
  Eigen::Quaterniond q(rot);
  return q.coeffs();
}

Eigen::Vector4d quatPointing(Eigen::Vector3d dir){ // this is used only for visualization and nowhere else!
  //quaternion that points the z vector to wherever you want
  // dir.normalize();
  Eigen::Vector3d axis(1,0,0);
  Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(axis,dir);
  return q.coeffs();
  // Eigen::Vector3d axis;
  // axis << 1,0,0;
  // dir.normalize();

  // Eigen::Vector3d v = skew3d(axis) * dir;
  // double s = v.norm(); //sine of angle
  // double c = (axis.transpose() * dir); //cosine of angle

  // Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + skew3d(v) + (1.0/(1+c))*(skew3d(v)*skew3d(v));
  // // R.setIdentity();
  // // R = rpyToRot(M_PI/2,0,0);
  // // R = rpyToRot(0,M_PI/2,0);
  // // R = rpyToRot(0,0,M_PI/2);
  // return rotToQuat(R);
}

class Trans {
public:
  char typ; //'f' fixed, 'r' revolute, 'b' (floating) base, (to be added: 'p' prismatic)
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
    if(typ == 'b'){ // floating base!
      gcIdx = -2; //(actually it is 0~6)
      gvIdx = -2; //(actually it is 0~5)
    }
  }

  Trans(const Eigen::Vector3d& r, const Eigen::Matrix3d& R){ //simple fixed creation (result of evalTrans)
    typ = 'f';
    gcIdx = -1;
    gvIdx = -1;
    axis << 0,0,0;
    // initFixed();
    originPos = r;
    originRot = R;
  }

  Trans(){ // trivial trans
    typ = 'f';
    gcIdx = -1;
    gvIdx = -1;
    axis << 0,0,0;
    // initFixed();
    originPos << 0,0,0;
    originRot = Eigen::Matrix3d::Identity();
  }

  // "set" functions for streamlining hard-coding
  void setXyz (double x,double y,double z){originPos << x,y,z;}
  void setRpy (double r,double p,double y){originRot = rpyToRot(r,p,y);}
  void setAxis(double ax,double ay,double az){axis << ax,ay,az;}
  void setTyp (char newTyp = 'f'){typ=newTyp;}
  void setIdx (int gcIndex = -1,int gvIndex = -1){gcIdx = gcIndex;gvIdx = gvIndex;}

  // axis given as 3 numbers
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
    // if(newT.typ == 'b'){ //error: attaching floating base
    //   throw std::invalid_argument("cannot attach floating base!");
    // }
    // std::cout << "rotation det: " << originRot.determinant() << std::endl;
    originPos = originPos + originRot*newT.originPos;
    originRot = originRot*newT.originRot;
    if(newT.typ != 'f'){
      // if(typ != 'f'){ //error: double actuation! (this is actually impossible)
      //   throw std::invalid_argument("double-actuation detected!");
      // }
      typ = newT.typ;
      axis = newT.axis;
      gcIdx = newT.gcIdx;
      gvIdx = newT.gvIdx;
    }
  }

  // evaluate the full transform for the given generalized coord (pure function).
  // returns a fixed transformation
  Trans evalTrans(const Eigen::VectorXd &gc){
    // urdf convention does the "origin move" first, then the actuation wrt axis.
    Eigen::Vector3d newPos = originPos;
    Eigen::Matrix3d newRot = originRot;

    if(typ == 'r'){
      Eigen::Vector4d quat; // rotation about an axis can be expressed as quat.
      quat << cos(gc[gcIdx]/2),sin(gc[gcIdx]/2)*(axis/axis.norm());
      newRot = originRot * quatToRot(quat);
    }
    else if(typ == 'p'){ //!!: prismatic joint kinematics
      newPos = originPos + gc[gcIdx] * (originRot*axis);
    }
    else if (typ == 'b'){ // floating base (migrated from Link class to here at ex4)
      newPos = originPos + gc.segment(0,3);
      // newRot = originRot + quatToRot(gc.segment(3,4)); //THIS LINE TOOK 5 HOURS TO CORRECT (it should be *, not +)
      newRot = originRot * quatToRot(gc.segment(3,4)); //THIS LINE TOOK 5 HOURS TO CORRECT (it should be *, not +)
    }
    return {newPos,newRot}; // equivalent to Trans(newPos,newRot)
  }
};

class Inertia{
public:
  double m;
  Eigen::Matrix3d I;
  Trans com;

  bool isEmpty;

  Inertia(double mass,const Eigen::Matrix3d& momentInertia,const Trans& centerOfMass){ // standard init
    m = mass;
    I = momentInertia;
    com = centerOfMass;
    isEmpty = false;

    // rotate I right away if the com's rotation is not trivial
    straighten();

  }

  Inertia(const Inertia& i){
    m = i.m;
    I = i.I;
    com = i.com;
    isEmpty = i.isEmpty;
  }

  Inertia(){ // empty init
    m = 0.0;
    I = Eigen::Matrix3d::Zero();
    com = Trans();
    isEmpty = true;
    // clear();
  }

  void straighten(){
    if (!com.originRot.isIdentity()) {I = com.originRot * I * com.originRot.transpose();}
    com.originRot = Eigen::Matrix3d::Identity();
  }

  [[maybe_unused]] void clear(){
    m = 0.0;
    I = Eigen::Matrix3d::Zero();
    com = Trans();
    isEmpty = true;
  }

  void merge(Inertia* i2){
    // i2 must be in the same coordinate frame as i2
    straighten();
    i2->straighten();
    if(i2->isEmpty){return;}
    double m1 = m;
    double m2 = i2 -> m;
    Eigen::Vector3d rcom = (m1*com.originPos + m2*(i2->com.originPos))/(m1+m2);
    Eigen::Vector3d r1 = com.originPos-rcom;
    Eigen::Vector3d r2 = i2->com.originPos - rcom;

    m = m + i2->m;
    I = I + i2->I - m1* skew3d(r1) * skew3d(r1) - m2*skew3d(r2)* skew3d(r2);
    com.originPos = rcom;
  }

  void addInertia(Inertia* i){
    if(isEmpty){
      m = i->m;
      I = i->I;
      com = i->com;
      isEmpty = false;
    }
    else{
      merge(i);
    }
  }

  void setProfile(double x,double y,double z,double R,double P,double Y,
                  double mass,double ixx, double ixy, double ixz, double iyy, double iyz, double izz){
    com.setProfile(x,y,z,R,P,Y);
    m = mass;
    I = hcInertia(ixx,ixy,ixz,iyy,iyz,izz);
    isEmpty=false;
    straighten(); //important to straighten the orientation
  }

  /// (pure function) returns another inertia that is expressed on another frame
  [[nodiscard]] Inertia expressedIn(Trans expT) const{ // express inertia in this frame
    expT.attachTrans(com); //expT comes first (attach com to expT)
    return {m,I,expT};
    // return Inertia(m,I,expT); // (equivalent)
  }
};

class Link{
public:
  std::string name; //name
  char typ; //'b': base link, 'a': articulated link, 'e': end-effector (no actuation)
  Link* parent; //pointer to parent link
  std::vector<Link*> children; //pointer(s) to children link

  // ex2: !! changed this from "final" to "root" references!!

  //CONSTANTS
  Trans bodyT; // full transform from (parent^2<>parent) to (parent<>this) (r,R,p,gcIdx,gvIdx) ("assume parent^2 is fixed to world")
  Inertia bodyI; // body inertia (of the current body (expressed with body origin as origin)
  // transform order (urdf convention) r -> R -> p
  // translate by r (expressed in (p^2<>p) frame)
  // rotate by R (expressed in (p^2<>p) frame)
  // actuate wrt p by (gc[gcIdx]) (expressed in new frame(r -> R))

  // STATE VARIABLES (another transform)
  // Kinematic
  Trans worldT;    // transform from world to i  (world --> frame attached to parent (used for jacobian)) (has "axis" property!)
  Trans fullT;     // transform from world to i' (world --> frame attached to self)

  // Diff.kinematics
  Eigen::Vector3d worldV;  // translational velocity of worldT frame (i, not i')
  Eigen::Vector3d worldOm; // angular velocity of worldT frame (i, not i')

  Eigen::Vector3d fullV;  // translational velocity of fullT frame (i')
  Eigen::Vector3d fullOm; // angular velocity of fullT frame (i')

  // Accleration <-- this value is different based on calculation conditions, and reset for every calculation.
  Eigen::VectorXd worldA;
  Eigen::VectorXd fullA;
  
  Eigen::VectorXd fi; // generalized force!

  // Inertial
  Inertia fullI;  // transform body inertia to world inertia
  Inertia compI;   // composite inertia including own body and all supporting bodies (in world coordinates)

  Eigen::MatrixXd Ma;    // articulated body inertia matrix
  Eigen::VectorXd ba;    // articulated body bias force

  // flags
  bool calcKin;  // kinematics calculated (root->leaf)
  bool calcComp; // composite inertia calculated (leaf -> root)
  
  bool calcDiffKin; // (optional) differnential kinematics calulated (root -> leaf) 

  bool calcAcc; // (**) this value is often reset independent to "resetCalculations"
  bool calcNE;  // newton-euler eom solved (again, reset indep to "resetCalulations")



  Link(const std::string& n,const char linkTyp){
    name = n;
    typ = linkTyp;
    calcKin = false; // pass 1: kinematics calculated (root->leaf)
    calcComp = false; // pass 2: composite inertia calculated (leaf -> root)
    calcDiffKin = false;

    calcAcc = false;
    calcNE  = false;

    bodyT  = Trans();
    bodyI  = Inertia();

    worldT = Trans();
    fullT  = Trans();

    worldV  << 0,0,0;
    worldOm << 0,0,0;
    fullV   << 0,0,0;
    fullOm  << 0,0,0;

    worldA.resize(6);
    fullA.resize(6);
    worldA << 0,0,0,0,0,0;
    fullA  << 0,0,0,0,0,0;

    fi.resize(6);
    fi << 0,0,0,0,0,0;

    Ma.resize(6,6);
    Ma.setZero();
    ba.resize(6);
    ba.setZero();

    fullI = Inertia();
    compI  = Inertia();

    parent = nullptr; //this will be modified within the Robot Class
  }

  void addTrans(const Trans& newT){
    // add a transformation at the end of the link
    // only modifies constant properties
    bodyT.attachTrans(newT);
  }

  void addInertia(Inertia newI){
    bodyI.addInertia(&newI);
  }

  // Trans propLinkKin(const Eigen::VectorXd& gc){ //returns frame for (i') frame (instead of (i))
  //   return worldT.evalTrans(gc);
  // }

  void calcLinkKin(const Eigen::VectorXd& gc){
    if( (typ != 'b') && (!parent->calcKin) ){
      throw(std::invalid_argument("Parent Kinematics Not Calculated!"));
    }

    if(typ == 'b'){ //base case (no parent link, get wr / wR from base-pos)
      worldT = Trans();
    }
    else{
      // std::cout << "parent pos" << worldT.originPos << std::endl;
      // std::cout << "parent ori" << worldT.originRot << std::endl;
      // std::cout << "parent typ" << worldT.typ << std::endl;
      // worldT = parent->propLinkKin(gc);
      worldT = parent->fullT; //get full transform of parent (transform of parent's "parent joint" (attached to parent))
      // although each links have to keep axis information (for J), the output for kinematics must be evaluated!
    }
    worldT.attachTrans(bodyT); // just attach my own transform!
    // std::cout << "worldT resolved to " << worldT.originPos.transpose() << std::endl;
    fullT = worldT.evalTrans(gc); // fullT is attached to myself
    // std::cout << "fullT resolved to " << fullT.originPos.transpose() << std::endl;
    resolveFullI();
    calcKin = true;
    // return;
  }

  /// only call this inside calcLinkKin!
  void resolveFullI(){
    fullI = bodyI.expressedIn(fullT);
  }

  /// get veolcity of worldT frame
  void calcLinkDiffKin(const Eigen::VectorXd& gv){
    if(!calcKin){throw(std::invalid_argument("This link's kinematics was not calculated yet..."));}
    if(typ == 'b'){
      worldV = Eigen::Vector3d::Zero(); // remember, the "i" joint of floating base is fixed to ground.
      worldOm = Eigen::Vector3d::Zero();
      if(worldT.typ == 'b'){ // floating base case
        fullV = gv.segment(0,3);
        fullOm = gv.segment(3,3);
      }
      else{ //fixed base case
        fullV = worldV;
        fullOm = worldOm;
      }
      // std::cout << "base velocity calculated as " << fullV.transpose() << " / " << fullOm.transpose() << std::endl;
    }
    else{
      // get i frame velocity (fixed to parent)
      worldV = parent->fullV + skew3d(parent->fullOm)*(worldT.originPos - parent->fullT.originPos); // vb = va + wx(rab)
      worldOm = parent->fullOm; // rotation doesn't change

      // get i' frame velocity (actuated by speed gv[worldT.gvIdx] along worldT.axis)
      if(worldT.typ == 'r'){
        // std::cout << "revolute joint ang vel evaluation from " << worldOm.transpose();
        fullV = worldV; //translational velocity doesn't change
        fullOm = worldOm + gv[worldT.gvIdx] * (worldT.originRot * worldT.axis);
        // std::cout << " to " << fullOm.transpose() << std::endl;
      }
      else if(worldT.typ == 'p'){
        // ** remember: since the origins are no longer coincident, worldOm also affects fullV! (this took 3 hours to find)
        fullV = worldV + gv[worldT.gvIdx] * (worldT.originRot * worldT.axis) + skew3d(worldOm)*(fullT.originPos - worldT.originPos);
        fullOm = worldOm;
        // std::cout << "prismatic joint velocity" << std::endl;
        // std::cout << "(i) frame: v: " << worldV.transpose() << " / om: " << worldOm.transpose() << std::endl;
        // std::cout << "axis: " << (worldT.originRot * worldT.axis).transpose() << " / velocity: " << gv[worldT.gvIdx] << std::endl;
        // std::cout << "(i')frame: v: " << fullV.transpose() << " / om: " << fullOm.transpose() << std::endl;
        // std::cout << std::endl;
        /// !!: prismatic diff.kin
      }
      else if(worldT.typ == 'f'){ // fixed joints simply transfer over the velocities.
        fullV = worldV;
        fullOm = worldOm;
      }

      // std::cout << name << " velocity calculated as " << fullV.transpose() << " / " << fullOm.transpose() << std::endl;
    }
    calcDiffKin = true;
  }


  /// add influence of this link to the provided Jp
  void augmentJp(Eigen::MatrixXd& Jp, const Eigen::VectorXd& wree){
    // wree: the position of end effector (point of question) in world coordinates.
    if (!calcKin){throw(std::invalid_argument("This link's kinematics is not calculated!"));}

    if(typ == 'b'){
      Jp.block<3,3>(0,0) = Eigen::Matrix3d::Identity();
      Jp.block<3,3>(0,3) = -skew3d(wree-worldT.originPos);
      // std::cout << "computed arm : " << (wree-worldT.originPos).transpose() << std::endl;
      // std::cout << "BASE vJ" << std::endl << vJ << std::endl;
    }
    else if(typ == 'a' || typ == 'e'){
      if (worldT.gvIdx >= 0){
        if (worldT.typ == 'r'){
          Jp.block<3,1>(0,worldT.gvIdx) = skew3d(worldT.originRot*worldT.axis)*(wree-worldT.originPos);
          // std::cout << "computed axis: " << (worldT.originRot*worldT.axis).transpose() << std::endl;
          // std::cout << "computed arm : " << (wree-worldT.originPos).transpose() << std::endl;
          // std::cout << "computed block: " << vJ.block<3,1>(0,worldT.gvIdx).transpose() << std::endl;
        }
        else if (worldT.typ == 'p'){
          Jp.block<3,1>(0,worldT.gvIdx) = worldT.originRot*worldT.axis;
        }
        // UNTESTED: prismatic positional jacobian
      }

      // vJ.block<3,1>(0,worldT.gvIdx) = skew3d(worldT.axis)*(wree-worldT.originPos);
    }
  }

  void augmentJa(Eigen::MatrixXd Ja, const Eigen::VectorXd& wree){
    // wree: the position of end effector (point of question) in world coordinates.
    if (!calcKin){throw(std::invalid_argument("This link's kinematics is not calculated!"));}

    if(typ == 'b'){
      // wJ.block<3,3>(0,0) = Eigen::Matrix3d::Zero()); //redundant
      Ja.block<3,3>(0,3) = Eigen::Matrix3d::Identity();
    }
    else if(typ == 'a' || typ == 'e'){
      if (worldT.gvIdx >= 0){
        if (worldT.typ == 'r'){
          Ja.block<3,1>(0,worldT.gvIdx) = worldT.originRot*worldT.axis;
        } 
        else if (worldT.typ == 'p'){
          Ja.block<3,1>(0,worldT.gvIdx) = Eigen::Vector3d::Zero();
        }
        //UNTESTED: prismatic angular jacobian
      }

      // Ja.block<3,1>(0,worldT.gvIdx) = skew3d(worldT.axis)*(wree-worldT.originPos);
    }
  }

  void calcCompInertia(){
    for(auto child:children){
      if(!(child->calcComp) || !(child->calcKin)){throw(std::invalid_argument("The children's composite inertia is not calculated!"));}
    }
    compI = fullI;
    for(auto child:children){
      compI.merge(&(child->compI));
    }
    // std::cout << "[" << name << "] Composite Inertia Calculated" << std::endl;

    calcComp = true;
    // return;
  }

  void calcLinkAcc(const Eigen::VectorXd& gv,const Eigen::VectorXd& gvdot){ // gvdot is often not a state but a condition for calculation. so it is passed as an arg on each call.
    if(!calcKin    ){throw(std::invalid_argument("Kinematics not calculated!!"));}
    if(!calcDiffKin){throw(std::invalid_argument("Velocities not calculated!"));}

    if( (typ != 'b') && (!parent->calcAcc) ){
      throw(std::invalid_argument("Parent Acceleration Not Calculated!"));
    }
    if(calcAcc){return;}

    // propagating through link ( (i-1)' --> (i) )
    if(typ == 'b'){ //worldA for the base is often set by outside functions.
      // std::cout << "calculating base acceleration with world acceleration of: " << worldA.transpose() << std::endl;
    }
    if(typ != 'b'){
      Eigen::Vector3d r = (worldT.originPos - parent->fullT.originPos);
      // worldA: attached to parent
      worldA.segment(0,3) << (parent->fullA).segment(0,3) + skew3d( (parent->fullA).segment(3,3) )*r + skew3d(parent->fullOm)*(skew3d(parent->fullOm)*r);
      worldA.segment(3,3) << (parent->fullA).segment(3,3);
    }

    // propagating through joint ( i --> i' )
    if(typ == 'b'){ //floating base case
      if(worldT.typ == 'b'){fullA = worldA + gvdot.segment(0,6);}
      else{fullA = worldA;} //fixed base case
    }
    else if (typ == 'a'){

      if(worldT.typ == 'r'){
        fullA = worldA;
        fullA = fullA + getSdot() * gv[worldT.gvIdx] + getS() * gvdot[worldT.gvIdx];
      }
      else if(worldT.typ == 'p'){
        // again, since fullA and worldA's origins are not coincident, we need to compensate its effect
        // first "move over" worldA's origin
        Eigen::Vector3d rp = fullT.originPos - worldT.originPos;
        fullA.segment(0,3) << worldA.segment(0,3) + skew3d( worldA.segment(3,3) )* rp + skew3d(worldOm)*(skew3d(worldOm)*rp);
        fullA.segment(3,3) << worldA.segment(3,3);
        // std::cout << "[" << name << "] prismatic acceleration propagation: intermediate Acc: " << fullA.transpose() << std::endl;
        // std::cout << "Sdot: " << getSdot().transpose() << " S:" << getS().transpose() << std::endl;
        // std::cout << "gv  : " << gv[worldT.gvIdx] << "  gvdot: " << gvdot[worldT.gvIdx] << std::endl;
        fullA = fullA + getSdot() * gv[worldT.gvIdx] + getS() * gvdot[worldT.gvIdx];
        // std::cout << "fullA: " << fullA.transpose() << std::endl;
      }
      else if(worldT.typ == 'f'){
        fullA = worldA;
      }
    }
    else{ // fixed joints have no axis / actuation (same with worldA)
      fullA = worldA;
    }

    calcAcc = true;
    if(fullI.m <= 1e-12 || worldT.gvIdx < 0){ //massless or no-axis case.
      fi << 0,0,0,0,0,0;
      calcNE = true;
    }
  }

  void calcLinkNE(bool gravity = false,Eigen::VectorXd extF = Eigen::VectorXd::Zero(6)){ //reset with acc (often uses same condition)
    // gravity: add gravity term when calculating bias force (b)
    // extF is assumed to be correctly transformed
    if(!calcKin    ){throw(std::invalid_argument("Kinematics not calculated!!"));}
    if(!calcDiffKin){throw(std::invalid_argument("Velocities not calculated!"));}

    if( (typ != 'b') && (!parent->calcAcc) ){throw(std::invalid_argument("Parent Acceleration Not Calculated!"));}
    for(auto child:children){
      if(! (child->calcNE)){
        throw(std::invalid_argument("Child NE is not yet solved!"));
      }
    }

    Eigen::Vector3d r;
    Eigen::Vector3d tempF;
    Eigen::Vector3d tempTau;
    Eigen::VectorXd tempSpatialF(6);

    // std::cout << "[" << name << "] calculating NE with acceleration: " << fullA.transpose() << std::endl;

    for(auto child:children){
      tempSpatialF.setZero();
      // r = child->worldT.originPos - fullT.originPos;
      r = child->fullT.originPos - fullT.originPos;      // this line took 4 hours to find...
      tempF = (child->fi).segment(0,3);
      tempTau = (child->fi).segment(3,3) + skew3d(r)*tempF;
      // tempTau = (child->fi).segment(3,3) - skew3d(r)*tempF;
      tempSpatialF << tempF,tempTau;
      extF += tempSpatialF;
      // std::cout << "  [" << name << "] added child force (" << child->name << ") r = <" << r.transpose() << "> dist = " << r.norm() << ": " << tempSpatialF.transpose() << std::endl;
    }
    // std::cout << "[" << name << "] total external force: " << extF.transpose() << std::endl;
    
    fi = getSpatialInertia() * fullA + getSpatialBiasVector(gravity) + extF; // the famous NE!

    // std::cout << "[" << name << "] calculated force at joint: " << fi.transpose() << std::endl;
    // std::cout << std::endl;
    calcNE = true;
  }

  void calcLinkABI(const Eigen::VectorXd& gv, const Eigen::VectorXd& gf){
    Ma = getSpatialInertia();
    // ba = getSpatialBiasVector(true); // [opt] include gravity in bias vector
    ba = getSpatialBiasVector();

    Eigen::MatrixXd Xbp(6,6);
    Eigen::MatrixXd Xbpdot(6,6); //...huh?
    Eigen::VectorXd wp(6);
    Eigen::MatrixXd Sb;
    Eigen::MatrixXd Sbdot;
    Eigen::MatrixXd Mb(6,6);
    // Eigen::VectorXd taub(1);
    int bIdx;

    for(auto c:children){ // base is never included as children!
      if(c->worldT.gvIdx < 0){continue;}
      Xbp.setIdentity();
      Xbp.block<3,3>(3,0) = skew3d(c->fullT.originPos - fullT.originPos);
      // Xbpdot.setIdentity();
      Xbpdot.setZero();
      Xbpdot.block<3,3>(3,0) = skew3d( skew3d(fullOm) * (c->fullT.originPos - fullT.originPos) ); //sketchy line
      wp << fullV,fullOm;
      Mb = c->Ma;
      bIdx = c->worldT.gvIdx;

      Sb = c->getS();
      Sbdot = c->getSdot();

      Ma += Xbp*(Mb)*( - Sb * (Sb.transpose() * Mb * Sb).inverse() * (Sb.transpose() * Mb * Xbp.transpose() ) + Xbp.transpose());

      ba += Xbp*(Mb *( (Sb * (Sb.transpose() * Mb * Sb).inverse()) * (gf.segment(bIdx,1) - Sb.transpose() * Mb * (Sbdot* gv[bIdx] + Xbpdot.transpose()*wp) - Sb.transpose() * c->ba) + Sbdot*gv[bIdx] + Xbpdot.transpose()*wp ) + c->ba);
    }

    // std::cout << "[" << name << "] Calculated Ma:" << std::endl << Ma << std::endl << "ba: " << ba.transpose() << std::endl;

  }

  Eigen::VectorXd calcAccFromABI(const Eigen::VectorXd& gv, const Eigen::VectorXd& gf){
    // we need to get ubdot (returned), fullA (6x1 vector)
    Eigen::VectorXd ubdot;
    if(typ == 'b'){ubdot.resize(6);}
    else{ubdot.resize(1);}

    Eigen::VectorXd wp(6), wpdot(6);
    Eigen::MatrixXd Xbp(6,6),Xbpdot(6,6);

    if(typ == 'b'){
      wp.setZero();
      wpdot.setZero();
      Xbp.setZero();
      Xbpdot.setZero();

      ubdot = (getS().transpose()*Ma*getS()).inverse() * (gf.segment(0,6) - getS().transpose()*Ma*(getSdot() * gv.segment(0,6)) - getS().transpose()*ba);
      // std::cout << "ubdot calculated size: " << ubdot.rows() << "x" << ubdot.cols() << std::endl;
      fullA = getS()*ubdot + getSdot()*gv.segment(0,6);

    }
    else{
      wp << parent->fullV,parent->fullOm;
      wpdot << parent->fullA;
      Xbp.setIdentity();
      Xbp.block<3,3>(3,0) = skew3d(fullT.originPos - parent->fullT.originPos);
      // Xbpdot.setIdentity();
      Xbpdot.setZero();
      Xbpdot.block<3,3>(3,0) = skew3d( skew3d(parent->fullOm) * (fullT.originPos - parent->fullT.originPos) ); //sketchy line

      ubdot = (getS().transpose()*Ma*getS()).inverse() * (gf.segment(worldT.gvIdx,1) - getS().transpose()*Ma*(getSdot() * gv[worldT.gvIdx] + Xbpdot.transpose()*wp + Xbp.transpose()*wpdot) - getS().transpose()*ba);
      // std::cout << "ubdot calculated size: " << ubdot.rows() << "x" << ubdot.cols() << std::endl;
      fullA = getS()*ubdot + getSdot()*gv[worldT.gvIdx] + Xbpdot.transpose()*wp + Xbp.transpose()*wpdot;
    }
    // std::cout << "ubdot for [" << name << "]: " << ubdot.transpose() << std::endl;
    calcAcc = true;
    return ubdot;
  }

  // "get" functions
  [[nodiscard]] Eigen::MatrixXd getS() const{ //retrieve motion subspace matrix
    Eigen::MatrixXd S;
    if(typ == 'b'){ //base: full space
      S.resize(6,6);
      S.setIdentity();
      return S;
    }

    if(typ == 'a'){
      S.resize(6,1);
      if(worldT.typ == 'r'){
        // sketchy line
        S << 0,0,0,worldT.originRot*worldT.axis; //worldT.axis is expressed in terms of worldT!
        return S;
      }
      else if(worldT.typ == 'p'){
        S << worldT.originRot*worldT.axis,0,0,0;
        return S;
      }
      // !!: prismatic motion subspace
    }

    if(typ == 'e'){ //no dof (this should not be possible)
      S.resize(1,1);
      S.setZero();
      return S;
    }
    throw(std::invalid_argument("Cannot determine motion subspace matrix"));
  }

  [[nodiscard]] Eigen::MatrixXd getSdot() const{ //retrieve motion subspace matrix //needs to be pre-declared (used by getLinkAcc)
    Eigen::MatrixXd Sdot;
    auto S = getS();
    if(typ == 'b'){ //base: full space
      Sdot.resize(6,6);
      // sketchy line
      // **axis of base floating joint is not moving
      Sdot.setZero();
      return Sdot;
    }

    if(typ == 'a'){
      Sdot.resize(6,1);
      if(worldT.typ == 'r'){
        // sketchy line
        Sdot << 0,0,0, skew3d(worldOm) * (worldT.originRot*worldT.axis); //worldT.axis is expressed in terms of worldT frame!
        return Sdot;
      }
      else if(worldT.typ == 'p'){
        // sketchy line
        Sdot << 2* skew3d(worldOm) * (worldT.originRot*worldT.axis),0,0,0; //worldT.axis is expressed in terms of worldT frame!
        // Sdot.setZero();
        // I have absolutely no idea why I have to multiply 2 here...
        return Sdot;
      }
      // !!SKETCHY!!: prismatic motion subspace derivative
    }

    if(typ == 'e'){ //no dof (this should not be possible)
      Sdot.resize(1,1);
      Sdot.setZero();
      return Sdot;
    }
    throw(std::invalid_argument("Cannot determine motion subspace derivative matrix"));
  }

  Eigen::MatrixXd getSpatialInertia(){
    if (!calcKin){throw(std::invalid_argument("This link's kinematics is not calculated!"));}

    Eigen::MatrixXd Mi;
    Mi.resize(6,6);
    Mi.setZero();

    Eigen::Matrix3d rx = skew3d(fullI.com.originPos - fullT.originPos);
    Mi.block<3,3>(0,0) = fullI.m * Eigen::Matrix3d::Identity();
    Mi.block<3,3>(3,0) = fullI.m*rx;
    Mi.block<3,3>(0,3) = -fullI.m*rx;
    Mi.block<3,3>(3,3) = fullI.I - fullI.m*rx*rx;
    // std::cout << "Mi for " << name << std::endl;
    // std::cout << Mi << std::endl;

    return Mi;    
  }

  Eigen::MatrixXd getSpatialCompositeInertia(){
    if (!calcKin){throw(std::invalid_argument("This link's kinematics is not calculated!"));}
    if (!calcComp){throw(std::invalid_argument("This link's composite inertia is not calculated!"));}

    Eigen::MatrixXd Mi;
    Mi.resize(6,6);
    Mi.setZero();
    Eigen::Matrix3d rx = skew3d(compI.com.originPos-fullT.originPos);
    Mi.block<3,3>(0,0) = compI.m * Eigen::Matrix3d::Identity();
    Mi.block<3,3>(3,0) = compI.m*rx;
    Mi.block<3,3>(0,3) = -compI.m*rx;
    Mi.block<3,3>(3,3) = compI.I - compI.m*rx*rx;
    // std::cout << "Mi for " << name << std::endl;
    // std::cout << Mi << std::endl;

    return Mi;
  }

  Eigen::VectorXd getSpatialBiasVector(bool gravity = false){
    Eigen::VectorXd bi;
    bi.resize(6);
    Eigen::Vector3d r = fullI.com.originPos-fullT.originPos;

    bi.segment(0,3) << fullI.m * (skew3d(fullOm) * (skew3d(fullOm) * r));
    bi.segment(3,3) << skew3d(fullOm) * (fullI.I - fullI.m * (skew3d(r) * skew3d(r) ) ) * fullOm;

    if(gravity){ // add an upward force of mg at com ()
      // std::cout << "[" << name << "] gravity term changed b from " << std::endl << bi.transpose() << std::endl;
      Eigen::VectorXd bg(6);
      Eigen::Vector3d mg;
      mg << 0,0,-9.81*fullI.m;
      bg << -mg,(skew3d(-r)*mg);
      bi += bg;
      // std::cout << "to " << std::endl << bi.transpose() << std::endl;
    }

    return bi;
  }

};

class Robot{
public:
  std::vector<Link*> links; //follows "parent-first" convention
  size_t gvDim; // (==dof)
  size_t dof;
  size_t gcDim;

  // state variables
  Eigen::VectorXd gc; //gen. coord.
  Eigen::VectorXd gv; //gen. vel.

  Eigen::MatrixXd M; // mass matrix
  Eigen::VectorXd b; // bias force

  Eigen::VectorXd gf; // generalized force

  Eigen::VectorXd udot; // acceleration calculated from forward dynamics (gf = M*fdUdot + b)

  // flags
  bool calcKin;  // kinematics calculated
  bool calcDiffKin;
  bool calcComp; // composite inertias calculated
  bool calcAcc; // acceleration (<- reset often!)

  Link* root;

  Robot(size_t gcDimensions,size_t gvDimensions){
    gcDim = gcDimensions;
    gvDim = gvDimensions;
    gc.resize(gcDim);
    gv.resize(gvDim);
    dof = int(gvDimensions);
    root = nullptr;
    calcKin = false;
    calcComp = false;
    calcDiffKin = false;
    calcAcc = false;

    M.resize(int(gvDim),int(gvDim));
    M.setZero();

    b.resize(int(gvDim));
    b.setZero();

    gf.resize(gvDim);
    gf.setZero();

    udot.resize(gvDim);
    udot.setZero();
  }

  int findLinkIdx(Link* l){
    auto it = std::find(links.begin(),links.end(),l);
    if (it == links.end()){
      return -1;
    }
    else{
      return int(std::distance(links.begin(),it));
    }
  }

  int findLinkIdx(const std::string& n){
    for(size_t i = 0;i<links.size();i++){
      if(links[i]->name == n){return int(i);}
    }
    return -1;
  }

  int findLinkIdx(const int& gvIndex){ //find index with gv index
    bool fBase = (root->worldT.typ == 'b');
    if(fBase && gvIndex < 6){return 0;} //floating base link case
    for(size_t i = 0;i<links.size();i++){
      if(links[i]->bodyT.gvIdx == gvIndex){return int(i);}
    }
    return -1;
  }

  Link* getLinkByName(const std::string& n){
    return links[findLinkIdx(n)];
  }

  Link* getLinkByGvIdx(const int& gvIndex){
    if(gvIndex >= gvDim || gvIndex < 0){
      throw(std::invalid_argument("This index is out of range! >> " + std::to_string(gvIndex)));
    }
    return links[findLinkIdx(gvIndex)];
  }

  void addLink(Link* l,Link* p = nullptr){
    // validity check
    if ( ( p == nullptr ) && ( !links.empty() ) ){throw std::invalid_argument("double-base detected!");}
    else if( (p) && findLinkIdx(p) == -1){throw std::invalid_argument("parent not found!");}
    links.push_back(l);
    l->parent = p;
    if( p != nullptr){p->children.push_back(l);}
    else{root=l;}

    // std::cout << "initilized link [" << l->name << "] with rotation matrix: " << l->bodyT.originRot << std::endl;
    // std::cout << "added link " << l->name << std::endl;
  }

  [[maybe_unused]] void addLink(Link* l,const std::string& pName){
    if(pName.empty()){return addLink(l,nullptr);}
    if(findLinkIdx(pName) < 0){throw std::invalid_argument("no such parent with name" + pName);}
    return addLink(l,links[findLinkIdx(pName)]);
  }

  void setState(const Eigen::VectorXd& gcIn, const Eigen::VectorXd& gvIn){
    if(gcIn.size() != gcDim){throw std::invalid_argument("gc dimensions do not match.");}
    if(gvIn.size() != gvDim){throw std::invalid_argument("gv dimensions do not match.");}
    gc = gcIn;
    gv = gvIn;
    resetCalculations();
  }

  void setState(const Eigen::VectorXd& gcIn){ // sometimes, only gc is given.
    if(gcIn.size() != gcDim){throw std::invalid_argument("gc dimensions do not match.");}
    gc = gcIn;
    gv.setZero();
    resetCalculations();
  }  
  
  void setForce(const Eigen::VectorXd& gfIn){
    if(gfIn.size() != gvDim){throw std::invalid_argument("gf dimensions do not match.");}
    gf = gfIn;
    // you don't have to reset calculations here!
  }

  // important stuff
  // all the following functions assume a "parent-first" indexing of [links]
  void resetCalculations(){
    calcKin     = false;
    calcComp    = false;
    calcDiffKin = false;
    calcAcc     = false;
    for(auto l:links){
      l -> calcKin     = false;
      l -> calcComp    = false;
      l -> calcDiffKin = false;
      l -> calcAcc     = false;
      l -> calcNE      = false;
    }
    M.setZero();
  }

  void resetAcc(){
    calcAcc = false;
    for(auto l:links){l->calcAcc = false;l->calcNE = false;} // the calcNE flag only exists for links (not the robot)
  }

  void calculateKinematics(){
    if(calcKin){return;}
    for(auto l:links){
      l->calcLinkKin(gc);
      // std::cout << "kinematics for " << l->name << " calculated" << std::endl;
    }
    calcKin = true;
  }

  void calculateDiffKinematics(){
    if(!calcKin){throw std::invalid_argument("Kinematics not yet calculated");}
    if(calcDiffKin){return;}
    for(auto l:links){
      l->calcLinkDiffKin(gv);
      // std::cout << "kinematics for " << l->name << " calculated" << std::endl;
    }
    calcDiffKin = true;
  }

  /// Positional Jacobian (modifier) (for a point attached to linkName, with coordinate wree)
  void calculateJp(Eigen::MatrixXd& Jp,const std::string& linkName,const Eigen::Vector3d &wree){
    // initialize
    if(!calcKin){calculateKinematics();}
    Jp.resize(3,long(gvDim));
    Jp.setZero();

    Link* l = links[findLinkIdx(linkName)];
    while(l->parent != nullptr){
      l->augmentJp(Jp,wree);
      l = l->parent;
    }
    //final time for base
    l->augmentJp(Jp,wree);
  }

  /// Angular Jacobian
  void calculateJa(Eigen::MatrixXd& Ja,const std::string& linkName,const Eigen::Vector3d &wree){
    // initialize
    if(!calcKin){calculateKinematics();}  
    Ja.resize(3,long(gvDim));
    Ja.setZero();

    Link* l = links[findLinkIdx(linkName)];
    while(l->parent != nullptr){
      l->augmentJa(Ja,wree);
      l = l->parent;
    }
    //final time for base
    l->augmentJa(Ja,wree);
  }

  [[maybe_unused]] void calculateJ(Eigen::MatrixXd& J,const std::string& linkName,const Eigen::Vector3d &wree){
    J.resize(6,long(gvDim));
    J.setZero();

    Eigen::MatrixXd Jp;
    Eigen::MatrixXd Ja;
    calculateJp(Jp,linkName,wree);
    calculateJa(Ja,linkName,wree);

    J << Jp,Ja;
  }

  void calculateCompositeInertia(){
    if(!calcKin){throw(std::invalid_argument("This robot's kinematics was not yet calculated!"));}
    // reverse calculation!
    for(int idx = int(links.size())-1;idx >= 0;--idx){
      links[idx]->calcCompInertia();
      // std::cout << "Comp Inertia for " << links[idx] -> name << " finished" << std::endl;
    }
    calcComp = true;
  }

  [[maybe_unused]] bool isConnected(int i, int j){ //checks if gvIndex i and gvIndex j share a path to root.
    if((findLinkIdx(int(i)) == -1) || (findLinkIdx(int(j)) == -1)){return false;}
    if(i > j){ std::tie(i, j) = std::make_tuple(j, i);} //ensures j >= i;
    if(i == j){return true;}
    if(i <= 5){return true;} // base is connected to everything
    Link* tempL = getLinkByGvIdx(int(j)); //start with j and move up
    while(tempL->parent != nullptr){
      if(tempL->parent->worldT.gvIdx == i){return true;}
      tempL = tempL->parent;
    }
    return false;
  }

  void calculateMassMatrixCRBA(){
    if(!calcKin ){throw(std::invalid_argument("[calc M] This robot's kinematics was not yet calculated!"));}
    if(!calcComp){throw(std::invalid_argument("[calc M] This robot's composite inertia isn't populated!"));}

    // row i and column j
    Eigen::MatrixXd tempM; //temp variable for M_ij

    Eigen::MatrixXd tempX; // the [I;-rx;0;I] matrix
    tempX.resize(6,6);

    Link* li;
    Link* lj;

    int i;

    bool fBase = (root->worldT.typ == 'b');

    for(int j = int(gvDim);j >= 0;--j){ // from leaf to root
      if(fBase && (j>0 && j<6)){continue;} //floating base: (1~5) is redundant
      if(findLinkIdx(j) < 0){continue;} //not actuated
      lj = getLinkByGvIdx(j);
      li = lj;
      while(li != nullptr){
        i = li->worldT.gvIdx;
        if(fBase && li->typ == 'b'){i = 0;}
        if(i<0){break;} // fixed base case ("a parent is not actuated")

        tempX.setIdentity();
        if(i != j){tempX.block<3,3>(0,3) = -skew3d(lj->fullT.originPos - li->fullT.originPos);}
        tempM = lj->getS().transpose() * lj->getSpatialCompositeInertia() * tempX * li->getS();

        // std::cout << "M_" << i << j << std::endl;
        // std::cout << tempM << std::endl;

        if(fBase && (i == 0)){ // only applicable for floating base relations
          if(j == 0){ // 6x6 case
            M.block<6,6>(i,j) = tempM;
          }
          else{ // 6x1 case
            // std::cout << tempM << std::endl;
            M.block<6,1>(i,j) = tempM.transpose(); // !! careful !! (this line (without transpose) caused errors in debug mode)
            M.block<1,6>(j,i) = tempM;
          }
          // std::cout << "written as floating base case" << std::endl;
        }
        else{
          // std::cout << "putting in M_" << i << j << " = " << tempM << std::endl;
          M.block<1,1>(i,j) = tempM;
          if(i != j){M.block<1,1>(j,i) = tempM.transpose();}
        }
        li = li->parent;

      }
    }
  }

  void calculateAccelerations(Eigen::VectorXd worldAcc = Eigen::Vector3d::Zero(), Eigen::VectorXd gvdot = Eigen::Vector3d::Zero()){

    if(worldAcc.isZero()){
      worldAcc.resize(6);
      worldAcc.setZero();
    }

    if(gvdot.isZero()){
      gvdot.resize(gvDim);
      gvdot.setZero();
      // std::cout << "Calculating acceleration for gvdot = 0 case" << std::endl;
    }
    else{
      // std::cout << "supplied gvdot: " << gvdot.transpose() << std::endl;
    }



    root -> worldA = worldAcc;
    for (auto l:links){l -> calcLinkAcc(gv,gvdot);}
    calcAcc = true;
  }

  void calculateNonlinearTermRNE(bool gravity = true){
    //setting fictitous accleration values
    Eigen::VectorXd worldAcc(6);
    Eigen::VectorXd gvdot;
    worldAcc << 0,0,9.81*gravity,0,0,0; // [opt1] gravity as acceleration
    // worldAcc.setZero(); // [opt2] gravity as force
    bool fBase = root->worldT.typ == 'b';

    gvdot.resize(gvDim);
    gvdot.setZero();

    calculateAccelerations(worldAcc,gvdot);

    // std::cout << "calculated base acc is " << root->fullA.transpose() << std::endl;
    // std::cout << "calculated LF_HIP worldA is " << getLinkByName("LF_HIP")->worldA.transpose() << std::endl;
    // std::cout << "calculated LF_HIP fullA is " << getLinkByName("LF_HIP")->fullA.transpose() << std::endl;
    // std::cout << "calculated RH_HIP worldA is " << getLinkByName("RH_HIP")->worldA.transpose() << std::endl;
    // std::cout << "calculated RH_HIP fullA is " << getLinkByName("RH_HIP")->fullA.transpose() << std::endl;

    // std::cout << "calculated LF_FOOT fullA(=worldA) is " << getLinkByName("LF_FOOT")->worldA.transpose() << std::endl;
    // std::cout << "calculated RF_FOOT fullA(=worldA) is " << getLinkByName("RF_FOOT")->worldA.transpose() << std::endl;
    // std::cout << "calculated LH_FOOT fullA(=worldA) is " << getLinkByName("LH_FOOT")->worldA.transpose() << std::endl;
    // std::cout << "calculated RH_FOOT fullA(=worldA) is " << getLinkByName("RH_FOOT")->worldA.transpose() << std::endl;


    for(int i = int(gvDim) - 1; i >= 0;--i){
      if(fBase && i <= 5){ // base case
        auto link = root;
        // std::cout << "calculating NE for " << link->name << std::endl;
        link -> calcLinkNE(); // [opt1] gravity as acceleration
        // link -> calcLinkNE(gravity); // [opt2] gravity as force
        b.segment(0,6) << (link->getS()).transpose() * link->fi;
        break;
      }
      auto link = getLinkByGvIdx(i);
      // std::cout << "calculating NE for " << link->name << std::endl;
      link -> calcLinkNE(); // [opt1] gravity as acceleration
      // link -> calcLinkNE(gravity); // [opt2] gravity as force
      b.segment(i,1) << (link->getS()).transpose() * link->fi;
    }

    // std::cout << "calculated b: " << b.transpose() << std::endl;
    resetAcc();
  }

  void calculateUdotABA(){
    // pass 01: (root -> leaf) get joint positions and velocities
    // calculateKinematics();
    // calculateDiffKinematics();
    bool fBase = root->worldT.typ == 'b';
    
    if(!calcKin){throw(std::invalid_argument("[ABA] Kinematics not yet calculated"));}
    if(!calcDiffKin){throw(std::invalid_argument("[ABA] Velocites not yet calculated"));}

    Eigen::VectorXd a;
    a.resize(gvDim);
    // pass 02: (leaf -> root) get ABI for each link (with mass)
    for(int i = int(gvDim) - 1; i >= 0;--i){
      if(fBase && i < 5){break;} // floating base case
      auto link = getLinkByGvIdx(i);
      link->calcLinkABI(gv,gf);
    }
    // pass 03: calculate body acceleration
    a.setZero();
    for(int i = 0;i<gvDim;i++){
      if(fBase && i < 5){continue;} // floating base case
      // std::cout << "calc for index " << i <<std::endl;
      if(fBase && i == 5){
        a.segment(0,6) = root->calcAccFromABI(gv,gf);
        continue;
      }
      a[i] = (getLinkByGvIdx(i)->calcAccFromABI(gv,gf)) [0];

    }
    // std::cout << "ABA finished" << std::endl;

    // finally, subtract gravitational acceleration
    if(fBase){a[2] -= 9.81;}
    // std::cout << std::endl;
    // std::cout << "Calculated acceleration is " << std::endl << a.transpose() << std::endl;
    udot = a;
  }

  // the full service
  void calculateAll(){ // populate all calculated state variables (pos, vel, M, b)
    calculateKinematics(); //ex1 part
    calculateDiffKinematics(); // ex2 part
    
    calculateCompositeInertia(); //ex3 part
    calculateMassMatrixCRBA();

    calculateNonlinearTermRNE(); //ex4 part

    calculateUdotABA();
  }

  // end-level "get" methods
  Inertia* getInertia(const std::string &linkName){
    auto l = links[findLinkIdx(linkName)];
    if(!(l->calcKin)){throw(std::invalid_argument("Link Kinematics not yet calculated!"));}
    return &(l->fullI);
  }

  Eigen::Vector3d getPos(const std::string &linkName){
    if(!calcKin){calculateKinematics();}
    return getLinkByName(linkName)->fullT.originPos;
  }

  Eigen::Vector3d getVel(const std::string &linkName){
    if(!calcKin){calculateKinematics();}
    if(!calcDiffKin){calculateDiffKinematics();}
    return getLinkByName(linkName) -> fullV;
  }

  Eigen::Vector3d getAcc(const std::string &linkName){
    if(!calcKin){calculateKinematics();}
    if(!calcDiffKin){calculateDiffKinematics();}
    if(!calcAcc){throw(std::invalid_argument("Acceleration not calculated!"));}

    return getLinkByName(linkName) -> fullA.segment(0,3);
  }

};


// (2024-SPRING-FIN) Cartpole hard-coding
void initCartpoleTrans(Robot& robot){
  // HARD-CODING: Link transformations.

  Trans tempT = Trans();

  //////////////////////      BASE     ////////////////////////
  // "base"
  auto sliderBar    = new Link("sliderBar",'b'); // <--base
  tempT.setProfile(0,0,5.0,  0,0,0, 'f');
  sliderBar -> addTrans(tempT);
  // std::cout << "base transformation set as " << base -> bodyT.typ << std::endl;
  robot.addLink(sliderBar);

  auto slider    = new Link("slider",'a');
  tempT.setProfile(0,0,0, 0,0,0, 'p', 1,0,0, 0,0);
  slider -> addTrans(tempT);
  robot.addLink(slider,sliderBar);

  auto rod    = new Link("rod",'a');
  tempT.setProfile(0,0,0, 0,0,0, 'r', 0,1,0, 1,1);
  rod -> addTrans(tempT);
  robot.addLink(rod,slider);

  auto rod2    = new Link("rod2",'a');
  tempT.setProfile(0,0,0.8, 0,0,0, 'r', 0,1,0, 2,2);
  rod2 -> addTrans(tempT);
  robot.addLink(rod2,rod);

  auto ee      = new Link("ee",'e');
  tempT.setProfile(0,0,0.9, 0,0,0);
  ee -> addTrans(tempT);
  robot.addLink(ee,rod2);
}

void initCartpoleInertia(Robot& robot){
  auto tempT  = Trans(); // the "cumulating" trans
  auto tempTT = Trans(); // the "adding" trans
  auto tempI  = Inertia();

  // mass value="6.222"/> <inertia ixx="0.017938806" ixy="0.00387963" ixz="0.001500772" iyy="0.370887745" iyz="0.370887745" izz="0.372497653"/>
  tempI.setProfile(0,0,0, 0,0,0,  0,  1.0,0,0,1,0,1);
  robot.getLinkByName("sliderBar")->addInertia(tempI);

  tempI.setProfile(0,0,0, 0,0,0,  2,  2,0,0,1,0,2);
  robot.getLinkByName("slider")->addInertia(tempI);

  tempI.setProfile(0,0,0.5, 0,0,0,  5,  1.0,0,0,1,0,1);
  robot.getLinkByName("rod")->addInertia(tempI);

  tempI.setProfile(0,0,0.5, 0,0,0,  5,  1.0,0,0,1,0,1);
  robot.getLinkByName("rod2")->addInertia(tempI);

}

Robot* initCartpole(){
  static Robot robot = Robot(3,3);
  if(!robot.links.empty()){
    robot.resetCalculations();
  }
  else{
    initCartpoleTrans(robot);
    initCartpoleInertia(robot);
  }
  return &robot;
}

/// do not change the name of the method
inline Eigen::VectorXd computeGeneralizedAcceleration (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& gf) {
  auto r = initCartpole();
  r->setState(gc,gv);
  r->setForce(gf);

  // r->calculateAll();

  r->calculateKinematics();
  r->calculateDiffKinematics();

  r->calculateNonlinearTermRNE();

  r->calculateUdotABA();

  return r->udot;
}