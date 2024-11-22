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
  dir.normalize();
  Eigen::Vector3d z(0,0,1);
  Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(z,dir);
  return q.coeffs();
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
    else if(typ == 'p'){ //TODO: UNTESTED
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

  // declared methods
  // [[nodiscard]] Eigen::MatrixXd getSdot();

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
    // std::cout << "fullT position is: " << fullT.originPos.transpose() << std::endl;
    resolveFullI();
    calcKin = true;
    // return;
  }

  /// only call this inside calcLinkKin!
  void resolveFullI(){
    fullI = bodyI.expressedIn(fullT);
  }

  /// get veolcity of worldT frame
  void calcLinkDiffkin(const Eigen::VectorXd& gv){
    if(!calcKin){throw(std::invalid_argument("This link's kinematics was not calculated yet..."));}
    if(typ == 'b'){
      worldV = Eigen::Vector3d::Zero(); // remember, the "i" joint of floating base is fixed to ground.
      worldOm = Eigen::Vector3d::Zero();
      fullV = gv.segment(0,3);
      fullOm = gv.segment(3,3);
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
        /// TODO prismatic diff.kin
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
        } //TODO: prismatic jacobian
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
        } //TODO: prismatic jacobian
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

    if(typ == 'b'){ //worldA for the base is often set by outside functions.
      // std::cout << "calculating base acceleration with world acceleration of: " << worldA.transpose() << std::endl;
    }
    if(typ != 'b'){
      Eigen::Vector3d r = (worldT.originPos - parent->fullT.originPos);
      // worldA: attached to parent
      worldA.segment(0,3) << (parent->fullA).segment(0,3) + skew3d( (parent->fullA).segment(3,3) )*r + skew3d(worldOm)*(skew3d(worldOm)*r);
      worldA.segment(3,3) << (parent->fullA).segment(3,3);
    }

    if(typ == 'b'){
      fullA = worldA + gvdot.segment(0,6);
    }
    else if (typ == 'a'){
      fullA = worldA + getSdot() * gv[worldT.gvIdx] + getS() * gvdot[worldT.gvIdx];
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
      r = child->worldT.originPos - fullT.originPos;
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
      // TODO: prismatic motion subspace
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
        Sdot << 0,0,0, skew3d(worldOm) * (worldT.originRot*worldT.axis); //worldT.axis is expressed in terms of worldT!
        return Sdot;
      }
      // TODO: prismatic motion subspace
    }

    if(typ == 'e'){ //no dof (this should not be possible)
      Sdot.resize(1,1);
      Sdot.setZero();
      return Sdot;
    }
    throw(std::invalid_argument("Cannot determine motion subspace matrix"));
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
  Eigen::VectorXd b; // TODO bias force

  Eigen::VectorXd gf; // generalized force

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
    if(gvIndex < 6){return 0;} //base link case
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
      l->calcLinkDiffkin(gv);
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
    for(int j = int(gvDim);j >= 0;--j){ // from leaf to root
      if(j>0 && j<6){continue;} //base cases redundant
      if(findLinkIdx(j) < 0){continue;} //"not coded yet"
      lj = getLinkByGvIdx(j);
      li = lj;
      while(li != nullptr){
        i = li->worldT.gvIdx;
        if(li->typ == 'b'){i = 0;}

        tempX.setIdentity();
        if(i != j){tempX.block<3,3>(0,3) = -skew3d(lj->fullT.originPos - li->fullT.originPos);}
        tempM = lj->getS().transpose() * lj->getSpatialCompositeInertia() * tempX * li->getS();
        if(i == 0){
          if(j == 0){ // 6x6 case
            M.block<6,6>(i,j) = tempM;
          }
          else{ // 6x1 case
            // std::cout << tempM << std::endl;
            M.block<6,1>(i,j) = tempM.transpose(); // !! careful !! (this line (without transpose) caused errors in debug mode)
            M.block<1,6>(j,i) = tempM;
          }
        }
        else{
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
    // worldAcc << 0,0,9.81*gravity,0,0,0; // [opt1] gravity as acceleration
    worldAcc.setZero(); // [opt2] gravity as force

    gvdot.resize(gvDim);
    gvdot.setZero();

    calculateAccelerations(worldAcc,gvdot);

    // std::cout << "calculated base acc is" << root->fullA.transpose() << std::endl;
    // std::cout << "calculated LF   acc is" << getLinkByName("LF_FOOT")->fullA.transpose() << std::endl;
    // std::cout << "calculated LF_HIP acc is" << getLinkByName("LF_HIP")->fullA.transpose() << std::endl;
    // std::cout << "calculated LF_HIP acc is" << getLinkByName("LF_HIP")->worldA.transpose() << std::endl;


    for(int i = int(gvDim) - 1; i >= 0;--i){
      if(i <= 5){ // base case
        auto link = root;
        // std::cout << "calculating NE for " << link->name << std::endl;
        // link -> calcLinkNE(); // [opt1] gravity as acceleration
        link -> calcLinkNE(gravity); // [opt2] gravity as force
        b.segment(0,6) << (link->getS()).transpose() * link->fi;
        break;
      }
      auto link = getLinkByGvIdx(i);
      // std::cout << "calculating NE for " << link->name << std::endl;
      // link -> calcLinkNE(); // [opt1] gravity as acceleration
      link -> calcLinkNE(gravity); // [opt2] gravity as force
      b.segment(i,1) << (link->getS()).transpose() * link->fi;
    }
    // std::cout << "calculated b: " << b.transpose() << std::endl;
    resetAcc();
  }



  Eigen::VectorXd computeAccelerationABA(){
    // pass 01: (root -> leaf) get joint positions and velocities
    // calculateKinematics();
    // calculateDiffKinematics();
    if(!calcKin){throw(std::invalid_argument("[ABA] Kinematics not yet calculated"));}
    if(!calcDiffKin){throw(std::invalid_argument("[ABA] Velocites not yet calculated"));}

    Eigen::VectorXd a;
    a.resize(gvDim);
    // pass 02: (leaf -> root) get ABI for each link (with mass)
    for(int i = int(gvDim) - 1; i >= 5;--i){ // we end at 5 since index 0-5 is the base link
      auto link = getLinkByGvIdx(i);
      link->calcLinkABI(gv,gf);
    }
    // pass 03: calculate body acceleration
    a.setZero();
    for(int i = 5;i<gvDim;i++){
      // std::cout << "calc for index " << i <<std::endl;
      if(i == 5){
        a.segment(0,6) = root->calcAccFromABI(gv,gf);
        continue;
      }
      a[i] = (getLinkByGvIdx(i)->calcAccFromABI(gv,gf)) [0];

    }
    // std::cout << "ABA finished" << std::endl;

    // finally, subtract gravitational acceleration
    a[2] -= 9.81;
    std::cout << std::endl;

    // std::cout << "Calculated acceleration is " << std::endl << a.transpose() << std::endl;
    return a;
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

//HARD-CODE #1: adding transformation relations
void initRobotTrans(Robot& robot) {
  // HARD-CODING: Link transformations.

  Trans tempT = Trans();
  //////////////////////      BASE     ////////////////////////
  // "base"
  auto base    = new Link("base",'b');
  tempT.setProfile(0,0,0,  0,0,0, 'b');
  base -> addTrans(tempT);
  // std::cout << "base transformation set as " << base -> bodyT.typ << std::endl;
  robot.addLink(base);
  // std::cout << "base transformation set as " << base -> bodyT.typ << std::endl;
  tempT = Trans();
  /////////////////////////////////////////////////////////////

  ////////////////////// LF LEG START  ////////////////////////
  auto LFHip = new Link("LF_HIP",'a');
  // <base_LF_HAA> (fixed) <origin rpy="2.61799387799 0 0.0" xyz="0.2999 0.104 0.0"/>
  tempT.setProfile(0.2999, 0.104, 0.0,   2.61799387799, 0, 0.0);
  LFHip -> addTrans(tempT);

  // <<LF_HAA>> (revolute) <origin rpy="0 0 0" xyz="0 0 0"/> <axis xyz="1 0 0/>" (gc[7])
  tempT.setProfile(0, 0, 0,   0, 0, 0, 'r',  1, 0, 0, 7, 6);
  LFHip -> addTrans(tempT);
  robot.addLink(LFHip,base);

  auto LFThigh = new Link("LF_THIGH",'a');
  // <LF_HIP_LF_hip_fixed> (fixed) <origin rpy="-2.61799387799 0 0.0" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   -2.61799387799, 0, 0.0);
  LFThigh -> addTrans(tempT);

  // <LF_hip_fixed_LF_HFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="0.0599 0.08381 0.0"/>
  tempT.setProfile(0.0599, 0.08381, 0.0,   0, 0, 1.57079632679);
  LFThigh -> addTrans(tempT);

  // <<LF_HFE>> (revolute) <origin rpy="0 0 0" xyz="0 0 0"/> <axis xyz="1 0 0/>" (gc[8])
  tempT.setProfile(0, 0, 0,   0, 0, 0, 'r',  1, 0, 0, 8, 7);
  LFThigh -> addTrans(tempT);
  robot.addLink(LFThigh,LFHip);

  auto LFShank = new Link("LF_SHANK",'a');
  // <LF_THIGH_LF_thigh_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, -1.57079632679);
  LFShank -> addTrans(tempT);

  // <LF_thigh_fixed_LF_KFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="0.0 0.1003 -0.285"/>
  tempT.setProfile(0.0, 0.1003, -0.285,   0, 0, 1.57079632679);
  LFShank -> addTrans(tempT);

  // <<LF_KFE>> (revolute) <origin rpy="0 0 0" xyz="0 0 0"/> <axis xyz="1 0 0/>" (gc[9])
  tempT.setProfile(0, 0, 0,   0, 0, 0, 'r',  1, 0, 0, 9, 8);
  LFShank -> addTrans(tempT);
  robot.addLink(LFShank,LFThigh);

  auto LFFoot = new Link("LF_FOOT",'e');
  // <LF_shank_LF_shank_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, -1.57079632679);
  LFFoot -> addTrans(tempT);

  // <LF_shank_fixed_LF_FOOT> (fixed) <origin rpy="0 0 0" xyz="0.08795 0.01305 -0.33797"/>
  tempT.setProfile(0.08795, 0.01305, -0.33797,   0, 0, 0);
  LFFoot -> addTrans(tempT);
  robot.addLink(LFFoot,LFShank);

  ////////////////////// LF LEG FINISH ////////////////////////

  ////////////////////// RF LEG START  ////////////////////////

  auto RFHip = new Link("RF_HIP",'a');
  // <base_RF_HAA> (fixed) <origin rpy="-2.61799387799 0 0.0" xyz="0.2999 -0.104 0.0"/>
  tempT.setProfile(0.2999, -0.104, 0.0,   -2.61799387799, 0, 0.0);
  RFHip -> addTrans(tempT);

  // <<RF_HAA>> (revolute) <origin rpy="0 0 0" xyz="0 0 0"/> <axis xyz="1 0 0/>" (gc[10])
  tempT.setProfile(0, 0, 0,   0, 0, 0, 'r',  1, 0, 0, 10, 9);
  RFHip -> addTrans(tempT);
  robot.addLink(RFHip,base);

  auto RFThigh = new Link("RF_THIGH",'a');
  // <RF_HIP_RF_hip_fixed> (fixed) <origin rpy="2.61799387799 0 0.0" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   2.61799387799, 0, 0.0);
  RFThigh -> addTrans(tempT);

  // <RF_hip_fixed_RF_HFE> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0.0599 -0.08381 0.0"/>
  tempT.setProfile(0.0599, -0.08381, 0.0,   0, 0, -1.57079632679);
  RFThigh -> addTrans(tempT);

  // <<RF_HFE>> (revolute) <origin rpy="0 0 0" xyz="0 0 0"/> <axis xyz="-1 0 0/>" (gc[11])
  tempT.setProfile(0, 0, 0,   0, 0, 0, 'r',  -1, 0, 0, 11, 10);
  RFThigh -> addTrans(tempT);
  robot.addLink(RFThigh,RFHip);

  auto RFShank = new Link("RF_SHANK",'a');
  // <RF_THIGH_RF_thigh_fixed> (fixed) <origin rpy="0 0 1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, 1.57079632679);
  RFShank -> addTrans(tempT);

  // <RF_thigh_fixed_RF_KFE> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0.0 -0.1003 -0.285"/>
  tempT.setProfile(0.0, -0.1003, -0.285,   0, 0, -1.57079632679);
  RFShank -> addTrans(tempT);

  // <<RF_KFE>> (revolute) <origin rpy="0 0 0" xyz="0 0 0"/> <axis xyz="-1 0 0/>" (gc[12])
  tempT.setProfile(0, 0, 0,   0, 0, 0, 'r',  -1, 0, 0, 12, 11);
  RFShank -> addTrans(tempT);
  robot.addLink(RFShank,RFThigh);

  auto RFFoot = new Link("RF_FOOT",'e');
  // <RF_shank_RF_shank_fixed> (fixed) <origin rpy="0 0 1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, 1.57079632679);
  RFFoot -> addTrans(tempT);

  // <RF_shank_fixed_RF_FOOT> (fixed) <origin rpy="0 0 0" xyz="0.08795 -0.01305 -0.33797"/>
  tempT.setProfile(0.08795, -0.01305, -0.33797,   0, 0, 0);
  RFFoot -> addTrans(tempT);
  robot.addLink(RFFoot,RFShank);

  ////////////////////// RF LEG FINISH ////////////////////////

  ////////////////////// LH LEG START  ////////////////////////
  // (base)

  // "HIP"
  auto LHHip = new Link("LH_HIP",'a');
  robot.addLink(LHHip,base);
  // base
  // <base_LH_HAA> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="-0.2999 0.104 0.0"/>
  // order: x y z   R P Y  (type = f) (ax ay az) (gcIdx gvIdx)
  tempT.setProfile(-0.2999,0.104,0.0,  -2.61799387799,0,-3.14159265359);
  LHHip -> addTrans(tempT);
  // LH_HAA

  // <<LH_HAA>> (revolute) <axis xyz="-1 0 0"/> (gc[13])
  tempT.setProfile(0.0,0.0,0.0,  0.0,0.0,0.0, 'r', -1,0.0,0.0, 13,12);
  LHHip -> addTrans(tempT);
  // LH_HIP

  // "THIGH"
  auto LHThi   = new Link("LH_THIGH",'a');
  robot.addLink(LHThi,LHHip);
  // LH_HIP
  // <LH_HIP_LH_hip_fixed> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="0 0 0"/>
  tempT.setProfile(0.0,0.0,0.0,  -2.61799387799,0.0,-3.14159265359);
  LHThi -> addTrans(tempT);
  // LH_hip_fixed
  // <LH_hip_fixed_LH_HFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0599 0.08381 0.0"/>
  tempT.setProfile(-0.0599,0.08381,0.0,  0.0,0.0,1.57079632679);
  LHThi -> addTrans(tempT);
  // LH_HFE
  // <<LH_HFE>> (revolute) <axis xyz="1 0 0"/> (gc[14])
  tempT.setProfile(0.0,0.0,0.0,  0.0,0.0,0.0, 'r', 1,0.0,0.0, 14,13);
  LHThi -> addTrans(tempT);

  // "SHANK"
  auto LHShank = new Link("LH_SHANK",'a');
  robot.addLink(LHShank,LHThi);
  // LH_THIGH
  // <LH_THIGH_LH_thigh_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0.0,0.0,0.0,  0.0,0.0,-1.57079632679);
  LHShank -> addTrans(tempT);
  // LH_thigh_fixed
  // <LH_thigh_fixed_LH_KFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0 0.1003 -0.285"/>
  tempT.setProfile(0.0,0.1003,-0.285,  0.0,0.0,1.57079632679);
  LHShank -> addTrans(tempT);
  // LH_KFE
  // <<LH_KFE>> (revolute) <axis xyz="1 0 0"/> (gc[15])
  tempT.setProfile(0.0,0.0,0.0,  0.0,0.0,0.0, 'r', 1,0.0,0.0, 15,14);
  LHShank -> addTrans(tempT);

  // "FOOT"
  auto LHFoot = new Link("LH_FOOT",'e');
  robot.addLink(LHFoot,LHShank);
  // LH_SHANK
  // <LH_SHANK_LH_shank_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0.0,0.0,0.0,  0.0,0.0,-1.57079632679);
  LHFoot -> addTrans(tempT);
  // LH_shank_fixed
  // <LH_shank_fixed_LH_FOOT> (fixed) <origin rpy="0 0 0" xyz="-0.08795 0.01305 -0.33797"/>
  tempT.setProfile(-0.08795,0.01305,-0.33797,  0.0,0.0,0.0);
  LHFoot -> addTrans(tempT);
  // LH_FOOT

  ////////////////////// LH LEG FINISH ////////////////////////

  ////////////////////// RH LEG START  ////////////////////////

  auto RHHip = new Link("RH_HIP",'a');
  // <base_RH_HAA> (fixed) <origin rpy="2.61799387799 0 -3.14159265359" xyz="-0.2999 -0.104 0.0"/>
  tempT.setProfile(-0.2999, -0.104, 0.0,   2.61799387799, 0, -3.14159265359);
  RHHip -> addTrans(tempT);

  // <<RH_HAA>> (revolute) <origin rpy="0 0 0" xyz="0 0 0"/> <axis xyz="-1 0 0/>" (gc[16])
  tempT.setProfile(0, 0, 0,   0, 0, 0, 'r',  -1, 0, 0, 16, 15);
  RHHip -> addTrans(tempT);
  robot.addLink(RHHip,base);

  auto RHThigh = new Link("RH_THIGH",'a');
  // <RH_HIP_RH_hip_fixed> (fixed) <origin rpy="2.61799387799 0 -3.14159265359" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   2.61799387799, 0, -3.14159265359);
  RHThigh -> addTrans(tempT);

  // <RH_hip_fixed_RH_HFE> (fixed) <origin rpy="0 0 -1.57079632679" xyz="-0.0599 -0.08381 0.0"/>
  tempT.setProfile(-0.0599, -0.08381, 0.0,   0, 0, -1.57079632679);
  RHThigh -> addTrans(tempT);

  // <<RH_HFE>> (revolute) <origin rpy="0 0 0" xyz="0 0 0"/> <axis xyz="-1 0 0/>" (gc[17])
  tempT.setProfile(0, 0, 0,   0, 0, 0, 'r',  -1, 0, 0, 17, 16);
  RHThigh -> addTrans(tempT);
  robot.addLink(RHThigh,RHHip);

  auto RHShank = new Link("RH_SHANK",'a');
  // <RH_THIGH_RH_thigh_fixed> (fixed) <origin rpy="0 0 1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, 1.57079632679);
  RHShank -> addTrans(tempT);

  // <RH_thigh_fixed_RH_KFE> (fixed) <origin rpy="0 0 -1.57079632679" xyz="-0.0 -0.1003 -0.285"/>
  tempT.setProfile(-0.0, -0.1003, -0.285,   0, 0, -1.57079632679);
  RHShank -> addTrans(tempT);

  // <<RH_KFE>> (revolute) <origin rpy="0 0 0" xyz="0 0 0"/> <axis xyz="-1 0 0/>" (gc[18])
  tempT.setProfile(0, 0, 0,   0, 0, 0, 'r',  -1, 0, 0, 18, 17);
  RHShank -> addTrans(tempT);
  robot.addLink(RHShank,RHThigh);

  auto RHFoot = new Link("RH_FOOT",'e');
  // <RH_shank_RH_shank_fixed> (fixed) <origin rpy="0 0 1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, 1.57079632679);
  RHFoot -> addTrans(tempT);

  // <RH_shank_fixed_RH_FOOT> (fixed) <origin rpy="0 0 0" xyz="-0.08795 -0.01305 -0.33797"/>
  tempT.setProfile(-0.08795, -0.01305, -0.33797,   0, 0, 0);
  RHFoot -> addTrans(tempT);
  robot.addLink(RHFoot,RHShank);

  ////////////////////// RH LEG FINISH ////////////////////////

}

//HARD-CODE #2: adding mass and inertia
void initRobotInertia(Robot& robot){
  auto tempT  = Trans(); // the "cumulating" trans
  auto tempTT = Trans(); // the "adding" trans
  auto tempI  = Inertia();

  ////////////////////// BASE   START   ////////////////////////
  // (base == base_inertia)
  // [base_inertia] <origin rpy="0 0 0" xyz="-0.018 -0.002 0.024"/>
  // mass value="6.222"/> <inertia ixx="0.017938806" ixy="0.00387963" ixz="0.001500772" iyy="0.370887745" iyz="0.370887745" izz="0.372497653"/>
  tempI.setProfile(-0.018,-0.002,0.024,  0,0,0,  6.222,  0.017938806, 0.00387963, 0.001500772, 0.370887745, 6.8963e-05, 0.372497653);
  robot.getLinkByName("base")->addInertia(tempI);

  // (base --"base_face_front"--> face_front)
  // <base_face_front> (fixed) <origin rpy="0 0 -0.0" xyz="0.4145 0 0"/>
  tempT.setProfile(0.4145, 0, 0,   0, 0, -0.0);
  // [face_front] <origin rpy="0 0 0" xyz="0.042 -0.001 0.004"/>
  // mass value="0.73"/> <inertia ixx="0.005238611" ixy="1.7609e-05" ixz="7.2167e-05" iyy="0.002643098" iyz="0.002643098" izz="0.004325938"/>
  tempI.setProfile(0.042,-0.001,0.004,  0,0,0,  0.73,  0.005238611, 1.7609e-05, 7.2167e-05, 0.002643098, 1.9548e-05, 0.004325938);
  robot.getLinkByName("base")->addInertia(tempI.expressedIn(tempT));

  // (base --"base_face_rear"--> face_rear)
  // <base_face_rear> (fixed) <origin rpy="0 0 3.14159265359" xyz="-0.4145 0 0"/>
  tempT.setProfile(-0.4145, 0, 0,   0, 0, 3.14159265359);
  // [face_rear] <origin rpy="0 0 0" xyz="0.042 -0.001 0.004"/>
  // mass value="0.73"/> <inertia ixx="0.005238611" ixy="1.7609e-05" ixz="7.2167e-05" iyy="0.002643098" iyz="0.002643098" izz="0.004325938"/>
  tempI.setProfile(0.042,-0.001,0.004,  0,0,0,  0.73,  0.005238611, 1.7609e-05, 7.2167e-05, 0.002643098, 1.9548e-05, 0.004325938);
  robot.getLinkByName("base")->addInertia(tempI.expressedIn(tempT));

  // (base --"base_battery"--> battery)
  // <base_battery> (fixed) <origin rpy="0 0 0" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, 0);
  // [battery] <origin rpy="0 0 0" xyz="-0.00067 -0.00023 -0.03362"/>
  // mass value="5.53425"/> <inertia ixx="0.00749474794" ixy="0.00016686282" ixz="7.82763e-05" iyy="0.0722338913" iyz="0.0722338913" izz="0.07482717535"/>
  tempI.setProfile(-0.00067,-0.00023,-0.03362,  0,0,0,  5.53425,  0.00749474794, 0.00016686282, 7.82763e-05, 0.0722338913, 1.42902e-06, 0.07482717535);
  robot.getLinkByName("base")->addInertia(tempI.expressedIn(tempT));

  // (base --"base_to_docking_hatch_cover"--> docking_hatch_cover)
  // <base_to_docking_hatch_cover> (fixed) <origin rpy="0 0 0" xyz="0.343 0.0 -0.07"/>
  tempT.setProfile(0.343, 0.0, -0.07,   0, 0, 0);
  // [docking_hatch_cover] <origin rpy="0 0 0" xyz="-0.003 0.0 0.005"/>
  // mass value="0.065"/> <inertia ixx="0.00063283" ixy="0.0" ixz="3.45e-07" iyy="0.00110971" iyz="0.00110971" izz="0.00171883"/>
  tempI.setProfile(-0.003,0.0,0.005,  0,0,0,  0.065,  0.00063283, 0.0, 3.45e-07, 0.00110971, 0.0, 0.00171883);
  robot.getLinkByName("base")->addInertia(tempI.expressedIn(tempT));

  // (base --"base_to_lidar_cage"--"lidar_cage_to_lidar"--> lidar)
  // <base_to_lidar_cage> (fixed) <origin rpy="0 0 0" xyz="-0.364 0.0 0.0735"/>
  tempT.setProfile(-0.364, 0.0, 0.0735,   0, 0, 0);
  // <lidar_cage_to_lidar> (fixed) <origin rpy="0.0 0.0 -1.57079632679" xyz="0.0 0.0 0.0687"/>
  tempTT.setProfile(0.0, 0.0, 0.0687,   0.0, 0.0, -1.57079632679);
  tempT.attachTrans(tempTT);
  // [lidar] <origin rpy="0.0 0.0 0.0" xyz="-0.012 0.001 -0.008"/>
  // mass value="0.695"/> <inertia ixx="0.000846765" ixy="6.9565e-05" ixz="0.00027111" iyy="0.001367583" iyz="0.001367583" izz="0.001363673"/>
  tempI.setProfile(-0.012,0.001,-0.008,  0.0,0.0,0.0,  0.695,  0.000846765, 6.9565e-05, 0.00027111, 0.001367583, 5.8984e-05, 0.001363673);
  robot.getLinkByName("base")->addInertia(tempI.expressedIn(tempT));

  // (base --"base_hatch"--> hatch)
  // <base_hatch> (fixed) <origin rpy="0 0 0" xyz="0.0 0.0 0.0"/>
  tempT.setProfile(0.0, 0.0, 0.0,   0, 0, 0);
  // [hatch] <origin rpy="0 0 0" xyz="0.116 0.0 0.0758"/>
  // mass value="0.142"/> <inertia ixx="0.001" ixy="0.001" ixz="0.001" iyy="0.001" iyz="0.001" izz="0.001"/>
  tempI.setProfile(0.116,0.0,0.0758,  0,0,0,  0.142,  0.001, 0.001, 0.001, 0.001, 0.001, 0.001);
  robot.getLinkByName("base")->addInertia(tempI.expressedIn(tempT));

  ////////////////////// BASE   FINISH  ////////////////////////

  ////////////////////// LF LEG START   ////////////////////////
  //---attached to BASE ---
  // (base --"base_LF_HAA" --> LF_HAA) (base added later)
  // <base_LF_HAA> (fixed) <origin rpy="2.61799387799 0 0.0" xyz="0.2999 0.104 0.0"/>
  tempT.setProfile(0.2999, 0.104, 0.0,   2.61799387799, 0, 0.0);
  // [LF_HAA] <origin rpy="0 0 0" xyz="-0.063 7e-05 0.00046"/>
  // mass value="2.04"/> <inertia ixx="0.001053013" ixy="4.527e-05" ixz="8.855e-05" iyy="0.001805509" iyz="0.001805509" izz="0.001765827"/>
  tempI.setProfile(-0.063,7e-05,0.00046,  0,0,0,  2.04,  0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827);
  robot.getLinkByName("base")->addInertia(tempI.expressedIn(tempT));

  //---attached to LF_HIP ---
  // (LF_HIP --"LF_HIP_LF_hip_fixed"--> LF_hip_fixed --> "LF_hip_fixed_LF_HFE" --> LF_HFE)
  // [LF_HIP] <origin rpy="0 0 0" xyz="0 0 0"/>
  // mass value="0.001"/> <inertia ixx="0.000001" ixy="0.0" ixz="0.0" iyy="0.000001" iyz="0.000001" izz="0.000001"/>
  tempI.setProfile(0,0,0,  0,0,0,  0.001,  0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001);
  robot.getLinkByName("LF_HIP")->addInertia(tempI);

  // <LF_HIP_LF_hip_fixed> (fixed) <origin rpy="-2.61799387799 0 0.0" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   -2.61799387799, 0, 0.0);
  // [LF_hip_fixed] <origin rpy="0 0 0" xyz="0.048 0.008 -0.003"/>
  // mass value="0.74"/> <inertia ixx="0.001393106" ixy="8.4012e-05" ixz="2.3378e-05" iyy="0.003798579" iyz="0.003798579" izz="0.003897509"/>
  tempI.setProfile(0.048,0.008,-0.003,  0,0,0,  0.74,  0.001393106, 8.4012e-05, 2.3378e-05, 0.003798579, 7.1319e-05, 0.003897509);
  robot.getLinkByName("LF_HIP")->addInertia(tempI.expressedIn(tempT));

  // <LF_hip_fixed_LF_HFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="0.0599 0.08381 0.0"/>
  tempTT.setProfile(0.0599, 0.08381, 0.0,   0, 0, 1.57079632679);
  tempT.attachTrans(tempTT);
  // [LF_HFE] <origin rpy="0 0 0" xyz="-0.063 7e-05 0.00046"/>
  // mass value="2.04"/> <inertia ixx="0.001053013" ixy="4.527e-05" ixz="8.855e-05" iyy="0.001805509" iyz="0.001805509" izz="0.001765827"/>
  tempI.setProfile(-0.063,7e-05,0.00046,  0,0,0,  2.04,  0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827);
  robot.getLinkByName("LF_HIP")->addInertia(tempI.expressedIn(tempT));

  //---attached to LF_THIGH ---
  // (LF_THIGH --"LF_THIGH_LF_thigh_fixed"--> LF_thigh_fixed --> "LF_thigh_fixed_LF_KFE" --> LF_KFE)
  // [LF_THIGH] <origin rpy="0 0 0" xyz="0 0 0"/>
  // mass value="0.001"/> <inertia ixx="0.000001" ixy="0.0" ixz="0.0" iyy="0.000001" iyz="0.000001" izz="0.000001"/>
  tempI.setProfile(0,0,0,  0,0,0,  0.001,  0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001);
  robot.getLinkByName("LF_THIGH")->addInertia(tempI);

  // <LF_THIGH_LF_thigh_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, -1.57079632679);
  // [LF_thigh_fixed] <origin rpy="0 0 0" xyz="0.0 0.018 -0.169"/>
  // mass value="1.03"/> <inertia ixx="0.018644469" ixy="5.2e-08" ixz="1.0157e-05" iyy="0.019312599" iyz="0.019312599" izz="0.002838361"/>
  tempI.setProfile(0.0,0.018,-0.169,  0,0,0,  1.03,  0.018644469, 5.2e-08, 1.0157e-05, 0.019312599, 0.002520077, 0.002838361);
  robot.getLinkByName("LF_THIGH")->addInertia(tempI.expressedIn(tempT));

  // <LF_thigh_fixed_LF_KFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="0.0 0.1003 -0.285"/>
  tempTT.setProfile(0.0, 0.1003, -0.285,   0, 0, 1.57079632679);
  tempT.attachTrans(tempTT);
  // [LF_KFE] <origin rpy="0 0 0" xyz="-0.063 7e-05 0.00046"/>
  // mass value="2.04"/> <inertia ixx="0.001053013" ixy="4.527e-05" ixz="8.855e-05" iyy="0.001805509" iyz="0.001805509" izz="0.001765827"/>
  tempI.setProfile(-0.063,7e-05,0.00046,  0,0,0,  2.04,  0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827);
  robot.getLinkByName("LF_THIGH")->addInertia(tempI.expressedIn(tempT));

  //---attached to LF_SHANK ---
  // (LF_SHANK --"LF_SHANK_LF_shank_fixed"--> LF_shank_fixed --> "LF_shank_fixed_LF_FOOT" --> LF_FOOT)

  // [LF_SHANK] <origin rpy="0 0 0" xyz="0 0 0"/>
  // mass value="0.001"/> <inertia ixx="0.000001" ixy="0.0" ixz="0.0" iyy="0.000001" iyz="0.000001" izz="0.000001"/>
  tempI.setProfile(0,0,0,  0,0,0,  0.001,  0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001);
  robot.getLinkByName("LF_SHANK")->addInertia(tempI);

  // <LF_shank_LF_shank_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, -1.57079632679);
  // [LF_shank_fixed] <origin rpy="0 0 0" xyz="0.03463 0.00688 0.00098"/>
  // mass value="0.33742"/> <inertia ixx="0.00032748005" ixy="2.142561e-05" ixz="1.33942e-05" iyy="0.00110974122" iyz="0.00110974122" izz="0.00089388521"/>
  tempI.setProfile(0.03463,0.00688,0.00098,  0,0,0,  0.33742,  0.00032748005, 2.142561e-05, 1.33942e-05, 0.00110974122, 7.601e-08, 0.00089388521);
  robot.getLinkByName("LF_SHANK")->addInertia(tempI.expressedIn(tempT));

  // <LF_shank_fixed_LF_FOOT> (fixed) <origin rpy="0 0 0" xyz="0.08795 0.01305 -0.33797"/>
  tempTT.setProfile(0.08795, 0.01305, -0.33797,   0, 0, 0);
  tempT.attachTrans(tempTT);
  // [LF_FOOT] <origin rpy="0 0 0" xyz="0.00948 -0.00948 0.1468"/>
  // mass value="0.25"/> <inertia ixx="0.00317174097" ixy="2.63048e-06" ixz="6.815581e-05" iyy="0.00317174092" iyz="0.00317174092" izz="8.319196e-05"/>
  tempI.setProfile(0.00948,-0.00948,0.1468,  0,0,0,  0.25,  0.00317174097, 2.63048e-06, 6.815581e-05, 0.00317174092, 6.815583e-05, 8.319196e-05);
  robot.getLinkByName("LF_SHANK")->addInertia(tempI.expressedIn(tempT));
  ////////////////////// LF LEG FINISH  ////////////////////////

  ////////////////////// RF LEG START   ////////////////////////

  //---attached to BASE ---
  // (base --"base_RF_HAA" --> RF_HAA) (base added later)
  // <base_RF_HAA> (fixed) <origin rpy="-2.61799387799 0 0.0" xyz="0.2999 -0.104 0.0"/>
  tempT.setProfile(0.2999, -0.104, 0.0,   -2.61799387799, 0, 0.0);
  // [RF_HAA] <origin rpy="0 0 0" xyz="-0.063 7e-05 0.00046"/>
  // mass value="2.04"/> <inertia ixx="0.001053013" ixy="4.527e-05" ixz="8.855e-05" iyy="0.001805509" iyz="0.001805509" izz="0.001765827"/>
  tempI.setProfile(-0.063,7e-05,0.00046,  0,0,0,  2.04,  0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827);
  robot.getLinkByName("base")->addInertia(tempI.expressedIn(tempT));

  //---attached to RF_HIP ---
  // (RF_HIP --"RF_HIP_RF_hip_fixed"--> RF_hip_fixed --> "RF_hip_fixed_RF_HFE" --> RF_HFE)
  // [RF_HIP] <origin rpy="0 0 0" xyz="0 0 0"/>
  // mass value="0.001"/> <inertia ixx="0.000001" ixy="0.0" ixz="0.0" iyy="0.000001" iyz="0.000001" izz="0.000001"/>
  tempI.setProfile(0,0,0,  0,0,0,  0.001,  0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001);
  robot.getLinkByName("RF_HIP")->addInertia(tempI);

  // <RF_HIP_RF_hip_fixed> (fixed) <origin rpy="2.61799387799 0 0.0" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   2.61799387799, 0, 0.0);
  // [RF_hip_fixed] <origin rpy="0 0 0" xyz="0.048 -0.008 -0.003"/>
  // mass value="0.74"/> <inertia ixx="0.001393106" ixy="-8.4012e-05" ixz="2.3378e-05" iyy="0.003798579" iyz="0.003798579" izz="0.003897509"/>
  tempI.setProfile(0.048,-0.008,-0.003,  0,0,0,  0.74,  0.001393106, -8.4012e-05, 2.3378e-05, 0.003798579, -7.1319e-05, 0.003897509);
  robot.getLinkByName("RF_HIP")->addInertia(tempI.expressedIn(tempT));

  // <RF_hip_fixed_RF_HFE> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0.0599 -0.08381 0.0"/>
  tempTT.setProfile(0.0599, -0.08381, 0.0,   0, 0, -1.57079632679);
  tempT.attachTrans(tempTT);
  // [RF_HFE] <origin rpy="0 0 0" xyz="-0.063 7e-05 0.00046"/>
  // mass value="2.04"/> <inertia ixx="0.001053013" ixy="4.527e-05" ixz="8.855e-05" iyy="0.001805509" iyz="0.001805509" izz="0.001765827"/>
  tempI.setProfile(-0.063,7e-05,0.00046,  0,0,0,  2.04,  0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827);
  robot.getLinkByName("RF_HIP")->addInertia(tempI.expressedIn(tempT));

  //---attached to RF_THIGH ---
  // (RF_THIGH --"RF_THIGH_RF_thigh_fixed"--> RF_thigh_fixed --> "RF_thigh_fixed_RF_KFE" --> RF_KFE)
  // [RF_THIGH] <origin rpy="0 0 0" xyz="0 0 0"/>
  // mass value="0.001"/> <inertia ixx="0.000001" ixy="0.0" ixz="0.0" iyy="0.000001" iyz="0.000001" izz="0.000001"/>
  tempI.setProfile(0,0,0,  0,0,0,  0.001,  0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001);
  robot.getLinkByName("RF_THIGH")->addInertia(tempI);

  // <RF_THIGH_RF_thigh_fixed> (fixed) <origin rpy="0 0 1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, 1.57079632679);
  // [RF_thigh_fixed] <origin rpy="0 0 0" xyz="0.0 -0.018 -0.169"/>
  // mass value="1.03"/> <inertia ixx="0.018644469" ixy="-5.2e-08" ixz="1.0157e-05" iyy="0.019312599" iyz="0.019312599" izz="0.002838361"/>
  tempI.setProfile(0.0,-0.018,-0.169,  0,0,0,  1.03,  0.018644469, -5.2e-08, 1.0157e-05, 0.019312599, -0.002520077, 0.002838361);
  robot.getLinkByName("RF_THIGH")->addInertia(tempI.expressedIn(tempT));

  // <RF_thigh_fixed_RF_KFE> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0.0 -0.1003 -0.285"/>
  tempTT.setProfile(0.0, -0.1003, -0.285,   0, 0, -1.57079632679);
  tempT.attachTrans(tempTT);
  // [RF_KFE] <origin rpy="0 0 0" xyz="-0.063 7e-05 0.00046"/>
  // mass value="2.04"/> <inertia ixx="0.001053013" ixy="4.527e-05" ixz="8.855e-05" iyy="0.001805509" iyz="0.001805509" izz="0.001765827"/>
  tempI.setProfile(-0.063,7e-05,0.00046,  0,0,0,  2.04,  0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827);
  robot.getLinkByName("RF_THIGH")->addInertia(tempI.expressedIn(tempT));

  //---attached to RF_SHANK ---
  // (RF_SHANK --"RF_SHANK_RF_shank_fixed"--> RF_shank_fixed --> "RF_shank_fixed_RF_FOOT" --> RF_FOOT)

  // [RF_SHANK] <origin rpy="0 0 0" xyz="0 0 0"/>
  // mass value="0.001"/> <inertia ixx="0.000001" ixy="0.0" ixz="0.0" iyy="0.000001" iyz="0.000001" izz="0.000001"/>
  tempI.setProfile(0,0,0,  0,0,0,  0.001,  0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001);
  robot.getLinkByName("RF_SHANK")->addInertia(tempI);

  // <RF_shank_RF_shank_fixed> (fixed) <origin rpy="0 0 1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, 1.57079632679);
  // [RF_shank_fixed] <origin rpy="0 0 0" xyz="0.03463 -0.00688 0.00098"/>
  // mass value="0.33742"/> <inertia ixx="0.00032748005" ixy="-2.142561e-05" ixz="1.33942e-05" iyy="0.00110974122" iyz="0.00110974122" izz="0.00089388521"/>
  tempI.setProfile(0.03463,-0.00688,0.00098,  0,0,0,  0.33742,  0.00032748005, -2.142561e-05, 1.33942e-05, 0.00110974122, -7.601e-08, 0.00089388521);
  robot.getLinkByName("RF_SHANK")->addInertia(tempI.expressedIn(tempT));

  // <RF_shank_fixed_RF_FOOT> (fixed) <origin rpy="0 0 0" xyz="0.08795 -0.01305 -0.33797"/>
  tempTT.setProfile(0.08795, -0.01305, -0.33797,   0, 0, 0);
  tempT.attachTrans(tempTT);
  // [RF_FOOT] <origin rpy="0 0 0" xyz="0.00948 0.00948 0.1468"/>
  // mass value="0.25"/> <inertia ixx="0.00317174097" ixy="-2.63048e-06" ixz="6.815581e-05" iyy="0.00317174092" iyz="0.00317174092" izz="8.319196e-05"/>
  tempI.setProfile(0.00948,0.00948,0.1468,  0,0,0,  0.25,  0.00317174097, -2.63048e-06, 6.815581e-05, 0.00317174092, -6.815583e-05, 8.319196e-05);
  robot.getLinkByName("RF_SHANK")->addInertia(tempI.expressedIn(tempT));
  ////////////////////// RF LEG FINISH  ////////////////////////

  ////////////////////// LH LEG START  ////////////////////////
  //---attached to BASE ---
  // (base --"base_LH_HAA" --> LH_HAA) (base added later)
  // <base_LH_HAA> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="-0.2999 0.104 0.0"/>
  tempT.setProfile(-0.2999,0.104,0.0,  -2.61799387799,0,-3.14159265359);
  // [LH_HAA]
  // <origin rpy="0 0 0" xyz="-0.063 7e-05 0.00046"/>
  // <mass value="2.04"/>
  // <inertia ixx="0.001053013" ixy="4.527e-05" ixz="8.855e-05" iyy="0.001805509" iyz="9.909e-05" izz="0.001765827"/>
  tempI.setProfile(-0.063,7e-05,0.00046,  0,0,0,  2.04,  0.001053013,4.527e-05,8.855e-05,0.001805509,9.909e-05,0.001765827);
  robot.getLinkByName("base")->addInertia(tempI.expressedIn(tempT));

  //---attached to LH_HIP ---
  // (LH_HIP --"LH_HIP_LH_hip_fixed"--> LH_hip_fixed --> "LH_hip_fixed_LH_HFE" --> LH_HFE)
  // [LH_HIP]
  // <origin rpy="0 0 0" xyz="0 0 0"/>
  // <mass value="0.001"/>
  // <inertia ixx="0.000001" ixy="0.0" ixz="0.0" iyy="0.000001" iyz="0.0" izz="0.000001"/>
  tempI.setProfile(0,0,0,  0,0,0,  0.001,  0.000001,0.0,0.0,0.000001,0.0,0.000001);
  robot.getLinkByName("LH_HIP")->addInertia(tempI);

  // <LH_HIP_LH_hip_fixed> (fixed) <origin rpy="-2.61799387799 0 -3.14159265359" xyz="0 0 0"/>
  tempT.setProfile(0.0,0.0,0.0,  -2.61799387799,0.0,-3.14159265359);
  // [LH_hip_fixed]
  // <origin rpy="0 0 0" xyz="-0.048 0.008 -0.003"/>
  // <mass value="0.74"/>
  // <inertia ixx="0.001393106" ixy="-8.4012e-05" ixz="-2.3378e-05" iyy="0.003798579" iyz="7.1319e-05" izz="0.003897509"/>
  tempI.setProfile(-0.048,0.008,-0.003,  0,0,0,  0.74,0.001393106,-8.4012e-05,-2.3378e-05,0.003798579,7.1319e-05,0.003897509);
  robot.getLinkByName("LH_HIP")->addInertia(tempI.expressedIn(tempT));

  // <LH_hip_fixed_LH_HFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0599 0.08381 0.0"/>
  tempTT.setProfile(-0.0599,0.08381,0.0,  0.0,0.0,1.57079632679);
  tempT.attachTrans(tempTT);
  // [LH_HFE]
  // <origin rpy="0 0 0" xyz="-0.063 7e-05 0.00046"/>
  // <mass value="2.04"/>
  // <inertia ixx="0.001053013" ixy="4.527e-05" ixz="8.855e-05" iyy="0.001805509" iyz="9.909e-05" izz="0.001765827"/>
  tempI.setProfile(-0.063,7e-05,0.00046,  0,0,0,  2.04,  0.001053013,4.527e-05,8.855e-05,0.001805509,9.909e-05,0.001765827);
  robot.getLinkByName("LH_HIP")->addInertia(tempI.expressedIn(tempT));

  //---attached to LH_THIGH ---
  // (LH_THIGH --"LH_THIGH_LH_thigh_fixed"--> LH_thigh_fixed --> "LH_thigh_fixed_LH_KFE" --> LH_KFE)
  //[LH_THIGH]
  // <origin rpy="0 0 0" xyz="0 0 0"/>
  // <mass value="0.001"/>
  // <inertia ixx="0.000001" ixy="0.0" ixz="0.0" iyy="0.000001" iyz="0.0" izz="0.000001"/>
  tempI.setProfile(0,0,0,0,0,0,0.001,0.000001,0.0,0.0,0.000001,0.0,0.000001);
  robot.getLinkByName("LH_THIGH")->addInertia(tempI);

  // <LH_THIGH_LH_thigh_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0.0,0.0,0.0,  0.0,0.0,-1.57079632679);
  // [LH_thigh_fixed]
  // <origin rpy="0 0 0" xyz="-0.0 0.018 -0.169"/>
  // <mass value="1.03"/>
  // <inertia ixx="0.018644469" ixy="-5.2e-08" ixz="-1.0157e-05" iyy="0.019312599" iyz="0.002520077" izz="0.002838361"/>
  tempI.setProfile(0,0.018,-0.169, 0,0,0,  1.03,  0.018644469,-5.2e-08,-1.0157e-05,0.019312599,0.002520077,0.002838361);
  robot.getLinkByName("LH_THIGH")->addInertia(tempI.expressedIn(tempT));

  // <LH_thigh_fixed_LH_KFE> (fixed) <origin rpy="0 0 1.57079632679" xyz="-0.0 0.1003 -0.285"/>
  tempTT.setProfile(0.0,0.1003,-0.285,  0.0,0.0,1.57079632679);
  tempT.attachTrans(tempTT);
  // [LH_KFE]
  // <origin rpy="0 0 0" xyz="-0.063 7e-05 0.00046"/>
  // <mass value="2.04"/>
  // <inertia ixx="0.001053013" ixy="4.527e-05" ixz="8.855e-05" iyy="0.001805509" iyz="9.909e-05" izz="0.001765827"/>
  tempI.setProfile(-0.063,7e-05,0.00046,0,0,0,2.04 ,0.001053013,4.527e-05,8.855e-05,0.001805509,9.909e-05,0.001765827);
  robot.getLinkByName("LH_THIGH")->addInertia(tempI.expressedIn(tempT));

  //---attached to LH_SHANK ---
  // (LH_SHANK --"LH_SHANK_LH_shank_fixed"--> LH_shank_fixed --> "LH_shank_fixed_LH_FOOT" --> LH_FOOT)

  // [LH_SHANK] <origin rpy="0 0 0" xyz="0 0 0"/>
  // mass value="0.001"/> <inertia ixx="0.000001" ixy="0.0" ixz="0.0" iyy="0.000001" iyz="0.000001" izz="0.000001"/>
  tempI.setProfile(0,0,0,  0,0,0,  0.001,  0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001);
  robot.getLinkByName("LH_SHANK")->addInertia(tempI);

  // <LH_shank_LH_shank_fixed> (fixed) <origin rpy="0 0 -1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0.0, 0, 0,   0, 0, -1.57079632679);
  // [LH_shank_fixed] <origin rpy="0 0 0" xyz="-0.03463 0.00688 0.00098"/>
  // mass value="0.33742"/> <inertia ixx="0.00032748005" ixy="-2.142561e-05" ixz="-1.33942e-05" iyy="0.00110974122" iyz="0.00110974122" izz="0.00089388521"/>
  tempI.setProfile(-0.03463,0.00688,0.00098,  0,0,0,  0.33742,  0.00032748005, -2.142561e-05, -1.33942e-05, 0.00110974122, 7.601e-08, 0.00089388521);
  robot.getLinkByName("LH_SHANK")->addInertia(tempI.expressedIn(tempT));

  // <LH_shank_fixed_LH_FOOT> (fixed) <origin rpy="0 0 0" xyz="-0.08795 0.01305 -0.33797"/>
  tempTT.setProfile(-0.08795, 0.01305, -0.33797,   0, 0, 0);
  tempT.attachTrans(tempTT);
  // [LH_FOOT] <origin rpy="0 0 0" xyz="-0.00948 -0.00948 0.1468"/>
  // mass value="0.25"/> <inertia ixx="0.00317174097" ixy="-2.63048e-06" ixz="-6.815581e-05" iyy="0.00317174092" iyz="0.00317174092" izz="8.319196e-05"/>
  tempI.setProfile(-0.00948,-0.00948,0.1468,  0,0,0,  0.25,  0.00317174097, -2.63048e-06, -6.815581e-05, 0.00317174092, 6.815583e-05, 8.319196e-05);
  robot.getLinkByName("LH_SHANK")->addInertia(tempI.expressedIn(tempT));

  //---attached to LH_FOOT --- (<- "ghost" link (nothing attached)) ---
  ////////////////////// LH LEG FINISH ////////////////////////

  ////////////////////// RH LEG START   ////////////////////////
  //---attached to BASE ---
  // (base --"base_RH_HAA" --> RH_HAA) (base added later)
  // <base_RH_HAA> (fixed) <origin rpy="2.61799387799 0 -3.14159265359" xyz="-0.2999 -0.104 0.0"/>
  tempT.setProfile(-0.2999, -0.104, 0.0,   2.61799387799, 0, -3.14159265359);
  // [RH_HAA] <origin rpy="0 0 0" xyz="-0.063 7e-05 0.00046"/>
  // mass value="2.04"/> <inertia ixx="0.001053013" ixy="4.527e-05" ixz="8.855e-05" iyy="0.001805509" iyz="0.001805509" izz="0.001765827"/>
  tempI.setProfile(-0.063,7e-05,0.00046,  0,0,0,  2.04,  0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827);
  robot.getLinkByName("base")->addInertia(tempI.expressedIn(tempT));

  //---attached to RH_HIP ---
  // (RH_HIP --"RH_HIP_RH_hip_fixed"--> RH_hip_fixed --> "RH_hip_fixed_RH_HFE" --> RH_HFE)
  // [RH_HIP] <origin rpy="0 0 0" xyz="0 0 0"/>
  // mass value="0.001"/> <inertia ixx="0.000001" ixy="0.0" ixz="0.0" iyy="0.000001" iyz="0.000001" izz="0.000001"/>
  tempI.setProfile(0,0,0,  0,0,0,  0.001,  0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001);
  robot.getLinkByName("RH_HIP")->addInertia(tempI);

  // <RH_HIP_RH_hip_fixed> (fixed) <origin rpy="2.61799387799 0 -3.14159265359" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   2.61799387799, 0, -3.14159265359);
  // [RH_hip_fixed] <origin rpy="0 0 0" xyz="-0.048 -0.008 -0.003"/>
  // mass value="0.74"/> <inertia ixx="0.001393106" ixy="8.4012e-05" ixz="-2.3378e-05" iyy="0.003798579" iyz="0.003798579" izz="0.003897509"/>
  tempI.setProfile(-0.048,-0.008,-0.003,  0,0,0,  0.74,  0.001393106, 8.4012e-05, -2.3378e-05, 0.003798579, -7.1319e-05, 0.003897509);
  robot.getLinkByName("RH_HIP")->addInertia(tempI.expressedIn(tempT));

  // <RH_hip_fixed_RH_HFE> (fixed) <origin rpy="0 0 -1.57079632679" xyz="-0.0599 -0.08381 0.0"/>
  tempTT.setProfile(-0.0599, -0.08381, 0.0,   0, 0, -1.57079632679);
  tempT.attachTrans(tempTT);
  // [RH_HFE] <origin rpy="0 0 0" xyz="-0.063 7e-05 0.00046"/>
  // mass value="2.04"/> <inertia ixx="0.001053013" ixy="4.527e-05" ixz="8.855e-05" iyy="0.001805509" iyz="0.001805509" izz="0.001765827"/>
  tempI.setProfile(-0.063,7e-05,0.00046,  0,0,0,  2.04,  0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827);
  robot.getLinkByName("RH_HIP")->addInertia(tempI.expressedIn(tempT));

  //---attached to RH_THIGH ---
  // (RH_THIGH --"RH_THIGH_RH_thigh_fixed"--> RH_thigh_fixed --> "RH_thigh_fixed_RH_KFE" --> RH_KFE)
  // [RH_THIGH] <origin rpy="0 0 0" xyz="0 0 0"/>
  // mass value="0.001"/> <inertia ixx="0.000001" ixy="0.0" ixz="0.0" iyy="0.000001" iyz="0.000001" izz="0.000001"/>
  tempI.setProfile(0,0,0,  0,0,0,  0.001,  0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001);
  robot.getLinkByName("RH_THIGH")->addInertia(tempI);

  // <RH_THIGH_RH_thigh_fixed> (fixed) <origin rpy="0 0 1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, 1.57079632679);
  // [RH_thigh_fixed] <origin rpy="0 0 0" xyz="-0.0 -0.018 -0.169"/>
  // mass value="1.03"/> <inertia ixx="0.018644469" ixy="5.2e-08" ixz="-1.0157e-05" iyy="0.019312599" iyz="0.019312599" izz="0.002838361"/>
  tempI.setProfile(-0.0,-0.018,-0.169,  0,0,0,  1.03,  0.018644469, 5.2e-08, -1.0157e-05, 0.019312599, -0.002520077, 0.002838361);
  robot.getLinkByName("RH_THIGH")->addInertia(tempI.expressedIn(tempT));

  // <RH_thigh_fixed_RH_KFE> (fixed) <origin rpy="0 0 -1.57079632679" xyz="-0.0 -0.1003 -0.285"/>
  tempTT.setProfile(-0.0, -0.1003, -0.285,   0, 0, -1.57079632679);
  tempT.attachTrans(tempTT);
  // [RH_KFE] <origin rpy="0 0 0" xyz="-0.063 7e-05 0.00046"/>
  // mass value="2.04"/> <inertia ixx="0.001053013" ixy="4.527e-05" ixz="8.855e-05" iyy="0.001805509" iyz="0.001805509" izz="0.001765827"/>
  tempI.setProfile(-0.063,7e-05,0.00046,  0,0,0,  2.04,  0.001053013, 4.527e-05, 8.855e-05, 0.001805509, 9.909e-05, 0.001765827);
  robot.getLinkByName("RH_THIGH")->addInertia(tempI.expressedIn(tempT));

  //---attached to RH_SHANK ---
  // (RH_SHANK --"RH_SHANK_RH_shank_fixed"--> RH_shank_fixed --> "RH_shank_fixed_RH_FOOT" --> RH_FOOT)

  // [RH_SHANK] <origin rpy="0 0 0" xyz="0 0 0"/>
  // mass value="0.001"/> <inertia ixx="0.000001" ixy="0.0" ixz="0.0" iyy="0.000001" iyz="0.000001" izz="0.000001"/>
  tempI.setProfile(0,0,0,  0,0,0,  0.001,  0.000001, 0.0, 0.0, 0.000001, 0.0, 0.000001);
  robot.getLinkByName("RH_SHANK")->addInertia(tempI);

  // <RH_shank_RH_shank_fixed> (fixed) <origin rpy="0 0 1.57079632679" xyz="0 0 0"/>
  tempT.setProfile(0, 0, 0,   0, 0, 1.57079632679);
  // [RH_shank_fixed] <origin rpy="0 0 0" xyz="-0.03463 -0.00688 0.00098"/>
  // mass value="0.33742"/> <inertia ixx="0.00032748005" ixy="2.142561e-05" ixz="-1.33942e-05" iyy="0.00110974122" iyz="0.00110974122" izz="0.00089388521"/>
  tempI.setProfile(-0.03463,-0.00688,0.00098,  0,0,0,  0.33742,  0.00032748005, 2.142561e-05, -1.33942e-05, 0.00110974122, -7.601e-08, 0.00089388521);
  robot.getLinkByName("RH_SHANK")->addInertia(tempI.expressedIn(tempT));

  // <RH_shank_fixed_RH_FOOT> (fixed) <origin rpy="0 0 0" xyz="-0.08795 -0.01305 -0.33797"/>
  tempTT.setProfile(-0.08795, -0.01305, -0.33797,   0, 0, 0);
  tempT.attachTrans(tempTT);
  // [RH_FOOT] <origin rpy="0 0 0" xyz="-0.00948 0.00948 0.1468"/>
  // mass value="0.25"/> <inertia ixx="0.00317174097" ixy="2.63048e-06" ixz="-6.815581e-05" iyy="0.00317174092" iyz="0.00317174092" izz="8.319196e-05"/>
  tempI.setProfile(-0.00948,0.00948,0.1468,  0,0,0,  0.25,  0.00317174097, 2.63048e-06, -6.815581e-05, 0.00317174092, -6.815583e-05, 8.319196e-05);
  robot.getLinkByName("RH_SHANK")->addInertia(tempI.expressedIn(tempT));
  ////////////////////// RH LEG FINISH  ////////////////////////
}

//separate function keeps initialization static throughout runtime
Robot* initRobot(){
  static Robot robot = Robot(19,18);
  if(!robot.links.empty()){
    // std::cout << "using static robot object (no re-init)" << std::endl;
    robot.resetCalculations();
  }
  else {
    initRobotTrans(robot);
    initRobotInertia(robot);
    // std::cout << "robot initialized with " << robot.links.size() << " bodies" << std::endl;
  }
  return &robot;
}

// the homework function

/// do not change the name of the method
inline Eigen::VectorXd computeGeneralizedAcceleration (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& gf) {
  auto r = initRobot();
  r->setState(gc,gv);
  r->setForce(gf);

  r->calculateKinematics();
  r->calculateDiffKinematics();

  return   r->computeAccelerationABA();
}