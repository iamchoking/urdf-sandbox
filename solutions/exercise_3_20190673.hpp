#pragma once

// PASTE HERE
#include <Eigen/Core>
#include <utility>
#include <vector>
#include <iostream>

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

  Trans(const Eigen::Vector3d& r, const Eigen::Matrix3d& R){ //simple creation (result of evalTrans)
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
    else if(typ == 'p'){ //UNTESTED
      newPos = originPos + gc[gcIdx] * (originRot*axis);
    }
    return {newPos,newRot}; // equivalent to evalTrans(newPos,newRot)
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

  Trans bodyT; // full transform from (parent^2<>parent) to (parent<>this) (r,R,p,gcIdx,gvIdx) ("assume parent^2 is fixed to world")
  Inertia bodyI; // body inertia (of the current body (expressed with body origin as origin)
  // transform order (urdf convention) r -> R -> p
  // translate by r (expressed in (p^2<>p) frame)
  // rotate by R (expressed in (p^2<>p) frame)
  // actuate wrt p by (gc[gcIdx]) (expressed in new frame(r -> R))

  // state variable (another transform)
  Trans worldT;    // transform from world to i  (world --> frame attached to parent (used for jacobian))
  Trans fullT;     // transform from world to i' (world --> frame attached to self)
  Inertia worldI;  // transform body inertia to world inertia
  Inertia compI;   // composite inertia including own body and all supporting bodies (in world coordinates)

  // flags
  bool calcKin;  // pass 1: kinematics calculated (root->leaf)
  bool calcComp; // pass 2: composite inertia calculated (leaf -> root)
  //*: non-constant

  Link(const std::string& n,const char linkTyp){
    name = n;
    typ = linkTyp;
    calcKin = false; // pass 1: kinematics calculated (root->leaf)
    calcComp = false; // pass 2: composite inertia calculated (leaf -> root)

    bodyT  = Trans();
    bodyI  = Inertia();

    worldT = Trans();
    fullT  = Trans();

    worldI = Inertia();
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
    if(typ == 'b'){ //base case (no parent link, get wr / wR from base-pos)
      // std::cout<< "base floating body: " << gc.transpose() << std::endl;
      worldT.originPos = gc.segment(0,3);
      worldT.originRot = quatToRot(gc.segment(3,4));
      fullT = worldT;
      resolveWorldI();
      calcKin = true;
      return;
    }
    else if(!parent->calcKin){
      throw(std::invalid_argument("Parent Kinematics Not Calculated!"));
    }
    // std::cout << "parent pos" << worldT.originPos << std::endl;
    // std::cout << "parent ori" << worldT.originRot << std::endl;
    // std::cout << "parent typ" << worldT.typ << std::endl;

    // worldT = parent->propLinkKin(gc);
    worldT = parent->fullT; //get full transform of parent (transform of parent's "parent joint" (attached to parent))
    worldT.attachTrans(bodyT); // just attach my own transform!
    fullT = worldT.evalTrans(gc);
    // although each links have to keep axis information (for J), the output for kinematics must be evaluated!
    resolveWorldI();
    calcKin = true;
    // return;
  }

  /// only call this inside calcLinkKin!
  void resolveWorldI(){
    worldI = bodyI.expressedIn(fullT);
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
    // std::cout << "[" << name << "] input check complete" << std::endl;
    compI = worldI;
    for(auto child:children){
      compI.merge(&(child->compI));
    }

    calcComp = true;
    // return;
  }

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
        S << 0,0,0,worldT.originRot*worldT.axis; //this is expressed in terms of worldT!
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

};

class Robot{
public:
  std::vector<Link*> links; //follows "parent-first" convention
  size_t gvDim; // (==dof)
  size_t dof;
  size_t gcDim;

  // state variables
  Eigen::MatrixXd M; // mass matrix
  Eigen::VectorXd b; // TODO bias force

  // flags
  bool calcKin;
  bool calcComp;

  Link* root;

  Robot(size_t gcDimensions,size_t gvDimensions){
    gcDim = gcDimensions;
    gvDim = gvDimensions;
    dof = int(gvDimensions);
    root = nullptr;
    calcKin = false;
    calcComp = false;

    M.resize(int(gvDim),int(gvDim));
    M.setZero();
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
    return links[findLinkIdx(gvIndex)];
  }

  void addLink(Link* l,Link* p = nullptr){
    // validity check
    if ( ( p == nullptr ) && ( !links.empty() ) ){throw std::invalid_argument("double-root detected!");}
    else if( (p) && findLinkIdx(p) == -1){throw std::invalid_argument("parent not found!");}
    links.push_back(l);
    l->parent = p;
    if( p != nullptr){p->children.push_back(l);}
    else{root=l;}

    // std::cout << "added link " << l->name << std::endl;
  }

  [[maybe_unused]] void addLink(Link* l,const std::string& pName){
    if(pName.empty()){return addLink(l,nullptr);}
    if(findLinkIdx(pName) < 0){throw std::invalid_argument("no such parent with name" + pName);}
    return addLink(l,links[findLinkIdx(pName)]);
  }

  // important stuff
  // all the following functions assume a "parent-first" indexing of [links]

  void calculateKinematics(const Eigen::VectorXd& gc){
    if(calcKin){return;}
    for(auto l:links){
      l->calcLinkKin(gc);
    }
    calcKin = true;
  }

  [[maybe_unused]] void resetKinematics(){
    calcKin = false;
    calcComp = false;
    for(auto l:links){
      l -> calcKin = false;
      l -> calcComp = false;
    }
    M.setZero();
  }

  /// Positional Jacobian (modifier)
  void calculateJp(Eigen::MatrixXd& Jp,const Eigen::VectorXd& gc,const std::string& linkName,const Eigen::Vector3d &wree){
    // initialize
    if(!calcKin){calculateKinematics(gc);}
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
  void calculateJa(Eigen::MatrixXd& Ja,const Eigen::VectorXd& gc,const std::string& linkName,const Eigen::Vector3d &wree){
    // initialize
    if(!calcKin){calculateKinematics(gc);}
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

  [[maybe_unused]] void calculateJ(Eigen::MatrixXd& J,const Eigen::VectorXd& gc,const std::string& linkName,const Eigen::Vector3d &wree){
    J.resize(6,long(gvDim));
    J.setZero();

    Eigen::MatrixXd Jp;
    Eigen::MatrixXd Ja;
    calculateJp(Jp,gc,linkName,wree);
    calculateJa(Ja,gc,linkName,wree);

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

  void calculateMassMatrix(){
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
            std::cout << tempM << std::endl;
            M.block<6,1>(i,j) = tempM.transpose(); // !! careful !! (this line caused errors in debug mode)
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

  // end-level "get" methodsF
  Trans* getTrans(const std::string &linkName){
    auto l = links[findLinkIdx(linkName)];
    if(!(l->calcKin)){throw(std::invalid_argument("Link Kinematics not yet calculated!"));}
    return &(l->fullT);
  }

  Inertia* getInertia(const std::string &linkName){
    auto l = links[findLinkIdx(linkName)];
    if(!(l->calcKin)){throw(std::invalid_argument("Link Kinematics not yet calculated!"));}
    return &(l->worldI);
  }

  Eigen::Vector3d getPos(const Eigen::VectorXd &gc,const std::string &linkName){
    calculateKinematics(gc);
    return getTrans(linkName)->originPos;
  }

  [[maybe_unused]] Eigen::Vector3d getCom(const Eigen::VectorXd &gc,const std::string &linkName){
    calculateKinematics(gc);
    return getInertia(linkName)->com.originPos;
  }

};

void initRobotTrans(Robot& robot) {
  // HARD-CODING: Link transformations.

  Trans tempT = Trans();
  //////////////////////      BASE     ////////////////////////
  // "base"
  auto base    = new Link("base",'b');
  robot.addLink(base);
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

Robot* initRobot(){
  static Robot robot = Robot(19,18);
  if(!robot.links.empty()){
    // std::cout << "using static robot object (no re-init)" << std::endl;
    robot.resetKinematics();
  }
  else {
    initRobotTrans(robot);
    initRobotInertia(robot);
    // std::cout << "robot initialized with " << robot.links.size() << " bodies" << std::endl;
  }
  return &robot;
}

inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {
  auto r = initRobot();
  r->resetKinematics();
  r->calculateKinematics(gc);
  r->calculateCompositeInertia();
  r->calculateMassMatrix();

  return r->M;
}

// TO HERE