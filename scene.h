#ifndef SCENE_HEADER_FILE
#define SCENE_HEADER_FILE

#include <vector>
#include <fstream>
#include <igl/bounding_box.h>
#include <igl/readOFF.h>
#include "ccd.h"
#include "volInt.h"
#include "auxfunctions.h"

using namespace Eigen;
using namespace std;


void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p);
void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir);
void center(const void *_obj, ccd_vec3_t *dir);

//double EPSILON = (double) std::numeric_limits<double>::epsilon;
bool RANDOM_VELOCITIES = false;
float GRAVITY = 9.807;
// Standard is none
float FRICTION_KINETIC = 0.2;
float FRICTION_STATIC = 0.3;
double frictionThreshold = 0.001;

float DRAG_COEFF = 1;

// Explicit gravity vector
RowVector3d gravityVec = RowVector3d(0.0, -GRAVITY, 0.0);

//Impulse is defined as a pair <position, direction>
typedef std::pair<RowVector3d,RowVector3d> Impulse;

double dRand(double dMin, double dMax)
{
	double f = (double)rand() / RAND_MAX;
	return dMin + f * (dMax - dMin);
}

//the class the contains each individual rigid objects and their functionality
class RigidObject{
public:
  MatrixXd origV;   //original vertex positions, where COM=(0.0,0.0,0.0) - never change this!
  MatrixXd currV;   //current vertex position
  MatrixXi T;
  
  //position of object in space. We must always have that currV = QRot(origV, orientation)+ COM
  RowVector4d orientation; //current orientation
  RowVector3d COM;  //current center of mass
  Matrix3d invIT;  //Original *inverse* inertia tensor around the COM, defined in the rest state to the object (so to the canonical world system)
  
  //kinematics
  bool isFixed;  //is the object immobile
  double mass;
  RowVector3d comVelocity;  //the linear velocity of the center of mass
  RowVector3d angVelocity;  //the angular velocity of the object.
  
  //dynamics
  std::vector<Impulse> impulses;  //Gets updated by the collision process
  
  
  //checking collision between bounding boxes, and consequently the rigid objects if succeeds.
  //you do not need to update these functions (isBoxCollide and isCollide) unless you are doing a different collision
  
  bool isBoxCollide(const RigidObject& ro){
    RowVector3d VMin1=currV.colwise().minCoeff();
    RowVector3d VMax1=currV.colwise().maxCoeff();
    RowVector3d VMin2=ro.currV.colwise().minCoeff();
    RowVector3d VMax2=ro.currV.colwise().maxCoeff();
    
    //checking all axes for non-intersection of the dimensional interval
    for (int i=0;i<3;i++)
      if ((VMax1(i)<VMin2(i))||(VMax2(i)<VMin1(i)))
        return false;
    
    return true;  //all dimensional intervals are overlapping = intersection
    
  }
  
  bool isCollide(const RigidObject& ro, double& depth, RowVector3d& intNormal, RowVector3d& intPosition){

    //collision between bounding boxes
    if (!isBoxCollide(ro))
      return false;
    
    ccd_t ccd;
    ccd_vec3_t sep;
    
    CCD_INIT(&ccd);
    //sophisticated collision between convex triangle meshes
    ccd.support1       = support; // support function for first object
    ccd.support2       = support; // support function for second object
    ccd.center1         =center;
    ccd.center2         =center;
    
    ccd.first_dir       = stub_dir;
    ccd.max_iterations = 100;     // maximal number of iterations
    
    void* obj1=(void*)this;
    void* obj2=(void*)&ro;
    
    ccd_real_t _depth;
    ccd_vec3_t dir, pos;
    
    // int nonintersect = ccdGJKPenetration(obj1, obj2,&ccd, &_depth,&dir,&pos);
    int nonintersect = ccdMPRPenetration(obj1, obj2, &ccd, &_depth, &dir, &pos);
    if (nonintersect)
      return false;
    
    
    for (int i=0;i<3;i++){
      intNormal(i)=dir.v[i];
      intPosition(i)=pos.v[i];
    }
    
    depth=_depth;
    intPosition-=depth*intNormal/2.0;  //to bring it to (current) obj2
    return !nonintersect;
  }
  
  
  //return the current inverted inertia tensor around the current COM. Update it by applying the orientation
  Matrix3d getCurrInvInertiaTensor(){
	  Matrix3d R;
	  R = Q2RotMatrix(orientation).transpose() * invIT * Q2RotMatrix(orientation);
	  return R;
  }
  
  
  //Update the current position and orientation by integrating the linear and angular velocities, and update currV accordingly
  //You need to modify this according to its purpose
  void updatePosition(double timeStep){
	  //Euler step
	  //COM += comVelocity * timeStep;

	  //Fourth order runge-kutta integration to update position
	  RowVector3d p2 = COM + comVelocity * timeStep * 0.5f;
	  RowVector3d p3 = p2 + comVelocity * timeStep * 0.5f;
	  RowVector3d p4 = p3 + comVelocity * timeStep;
	  COM = (COM + 2.0f * p2 + 2.0f * p3 + p4) / 6.0f;

	  //Integrate angular velocity to update orientation 
	  RowVector4d angularVelocityContainer(0, angVelocity(0, 0), angVelocity(0, 1), angVelocity(0, 2));
	  orientation += timeStep * 0.5 * QMult(angularVelocityContainer,orientation);
	  orientation.normalize();

	  //Update the position of every vertex with new orientation and COM.
	  for (int i = 0; i < currV.rows(); i++) {
		  currV.row(i) = QRot(origV.row(i), orientation) + COM; 
	  }  
  }
  
  //Updating velocity *instantaneously*. i.e., not integration from acceleration, but as a result of a collision impulse from the "impulses" list
  //You need to modify this for that purpose.
  void updateImpulseVelocities(){
    
    if (isFixed){
      comVelocity.setZero();
      impulses.clear();
      angVelocity.setZero();
      return;
    }
    
    //update linear and angular velocity according to all impulses
    for (int i=0;i<impulses.size();i++){

		comVelocity += impulses[i].second * (1 / mass);

		RowVector3d arm = impulses[i].first - COM;
		angVelocity += getCurrInvInertiaTensor() * arm.cross(impulses[i].second).transpose();
    }
    impulses.clear();
  }
  
  
  //Updating the linear and angular velocities of the object
  //You need to modify this to integrate from acceleration in the field (basically gravity)
  void updateVelocity(double timeStep){
    
    if (isFixed)
      return;
	
	RowVector3d netAcceleration = gravityVec - (DRAG_COEFF * comVelocity * (1 / mass));

	//Fourth order runge-kutta integration
	RowVector3d v2 = comVelocity + netAcceleration * timeStep * 0.5f;
	RowVector3d v3 = v2 + netAcceleration * timeStep * 0.5f;
	RowVector3d v4 = v3 + netAcceleration * timeStep;

	comVelocity = (comVelocity + 2.0f * v2 + 2.0f * v3 + v4) / 6.0f;
  }
  
  
  //the full integration for the time step (velocity + position)
  //You need to modify this if you are changing the integration
  void integrate(double timeStep){
    updateVelocity(timeStep);
    updatePosition(timeStep);
  }
  
  
  RigidObject(const MatrixXd& _V, const MatrixXi& _T, const double density, const bool _isFixed, const RowVector3d& _COM, const RowVector4d& _orientation){
    origV=_V;
    T=_T;
    isFixed=_isFixed;
    COM=_COM;
    orientation=_orientation;
    comVelocity.setZero();
    angVelocity.setZero();
    
    RowVector3d naturalCOM;  //by the geometry of the object
    
    //initializes the original geometry (COM + IT) of the object
    getCOMandInvIT(origV, T, density, mass, naturalCOM, invIT);
    
    origV.rowwise()-=naturalCOM;  //removing the natural COM of the OFF file (natural COM is never used again)
    
    currV.resize(origV.rows(), origV.cols());
    for (int i=0;i<currV.rows();i++)
      currV.row(i)<<QRot(origV.row(i), orientation)+COM;
  }
  
  ~RigidObject(){}
};

//This class contains the entire scene operations, and the engine time loop.
class Scene{
public:
  double currTime;
  int numFullV, numFullT;
  std::vector<RigidObject> rigidObjects;
  
  //adding an objects. You do not need to update this generally
  void addRigidObject(const MatrixXd& V, const MatrixXi& T, const double density, const bool isFixed, const RowVector3d& COM, const RowVector4d orientation){
    
    RigidObject ro(V,T, density, isFixed, COM, orientation);
    rigidObjects.push_back(ro);
    numFullV+=V.rows();
    numFullT+=T.rows();
  }
  
  /*********************************************************************
   This function handles a collision between objects ro1 and ro2 when found, by assigning impulses to both objects.
   Input: RigidObjects ro1, ro2
   depth: the depth of penetration
   contactNormal: the normal of the conact measured ro1->ro2
   penPosition: a point on ro2 such that if ro2 <= ro2 + depth*contactNormal, then penPosition+depth*contactNormal is the common contact point
   CRCoeff: the coefficient of restitution
   *********************************************************************/
  void handleCollision(RigidObject& ro1, RigidObject& ro2, const double depth, const RowVector3d& contactNormal, const RowVector3d& penPosition, const double CRCoeff){

    // Interpretation resolution: move each object by inverse mass weighting, unless either is fixed, and then move the other. Remember to respect the direction of contactNormal and update penPosition accordingly.
    RowVector3d contactPosition = penPosition + depth*contactNormal;
	
	double invSummedMass;

	double ro2Depth;
	double ro1Depth;

	// If both objects are free
	if (!(ro1.isFixed || ro2.isFixed)) {
		invSummedMass = 1 / (ro1.mass + ro2.mass);

		// Heavier object moves back less
		ro1Depth = depth * ro2.mass * invSummedMass;
		ro2Depth = depth * ro1.mass * invSummedMass;

		ro1.COM -= contactNormal * ro1Depth;
		for (int i = 0; i < ro1.currV.rows(); i++) {
			ro1.currV.row(i) = QRot(ro1.origV.row(i), ro1.orientation) + ro1.COM;
		}

		ro2.COM += contactNormal * ro2Depth;
		for (int i = 0; i < ro2.currV.rows(); i++) {
			ro2.currV.row(i) = QRot(ro2.origV.row(i), ro2.orientation) + ro2.COM;
		}
	}

	// If one object is fixed 
	else if (ro1.isFixed) {
		ro2.COM += contactNormal * depth;
		for (int i = 0; i < ro2.currV.rows(); i++) {
			ro2.currV.row(i) = QRot(ro2.origV.row(i), ro2.orientation) + ro2.COM;
		}
	}

	else if (ro2.isFixed) {
		ro1.COM -= contactNormal * depth;
		for (int i = 0; i < ro1.currV.rows(); i++) {
			ro1.currV.row(i) = QRot(ro1.origV.row(i), ro1.orientation) + ro1.COM;
		}
	}

	RowVector3d fullVelocity1 = ro1.comVelocity + ro1.angVelocity.cross(contactPosition - ro1.COM);
	RowVector3d fullVelocity2 = ro2.comVelocity + ro2.angVelocity.cross(contactPosition - ro2.COM);

	RowVector3d relativeVelocity = fullVelocity1 - fullVelocity2;

    //Create impulses and push them into ro1.impulses and ro2.impulses.
	RowVector3d arm1 = contactPosition - ro1.COM;
	RowVector3d arm2 = contactPosition - ro2.COM;

	RowVector3d armCrossNormal1 = arm1.cross(contactNormal);
	RowVector3d armCrossNormal2 = arm2.cross(contactNormal);

	double augTermA = armCrossNormal1 * ro1.getCurrInvInertiaTensor() * armCrossNormal1.transpose();
	double augTermB = armCrossNormal2 * ro2.getCurrInvInertiaTensor() * armCrossNormal2.transpose();

	RowVector3d impulse = -(1.0f + CRCoeff) * ((relativeVelocity).dot(contactNormal) * (contactNormal))
		/ ((1.0f / ro1.mass) + (1.0f / ro2.mass) + augTermA + augTermB);

	// Collision impulses
 	ro1.impulses.push_back(Impulse(contactPosition, impulse));
	ro2.impulses.push_back(Impulse(contactPosition, -impulse));

	//Friction Calculations - Disabled
	RowVector3d tangent = (contactNormal.cross(relativeVelocity)).cross(contactNormal);
	tangent.normalize();

	RowVector3d armCrossTangent1 = arm1.cross(tangent);
	RowVector3d armCrossTangent2 = arm2.cross(tangent);

	double augTermFricA = armCrossTangent1 * ro1.getCurrInvInertiaTensor() * armCrossTangent1.transpose();
	double augTermFricB = armCrossTangent2 * ro2.getCurrInvInertiaTensor() * armCrossTangent2.transpose();

	RowVector3d fricImpulse;
	if (relativeVelocity.norm() <= frictionThreshold) {
		fricImpulse = -(1.0f + CRCoeff) * ((relativeVelocity).dot(tangent) * FRICTION_STATIC *(tangent))
			/ ((1.0f / ro1.mass) + (1.0f / ro2.mass) + augTermFricA + augTermFricB);
	}
	else {
		fricImpulse = -(1.0f + CRCoeff) * ((relativeVelocity).dot(tangent) * FRICTION_KINETIC *(tangent))
			/ ((1.0f / ro1.mass) + (1.0f / ro2.mass) + augTermFricA + augTermFricB);
	}
	
	//FRICTION IMPULSES -- Works as intended on decreasing total impulse, BUT collision detection breaks when enabled. small fast objects pass through (?)
	//ro1.impulses.push_back(Impulse(contactPosition, fricImpulse));
	//ro2.impulses.push_back(Impulse(contactPosition, -fricImpulse));
	  
	//updating velocities according to impulses
	ro1.updateImpulseVelocities();
	ro2.updateImpulseVelocities();
  }
  
  /*********************************************************************
   This function handles a single time step by:
   1. Integrating velocities, positions, and orientations by the timeStep
   2. detecting and handling collisions with the coefficient of restitution CRCoeff
   3. updating the visual scene in fullV and fullT
   *********************************************************************/
  void updateScene(double timeStep, double CRCoeff, MatrixXd& fullV, MatrixXi& fullT){

	if (RANDOM_VELOCITIES) {
		  initalizeRandomVelocities(20.0f);
		  RANDOM_VELOCITIES = false;
	 }
    fullV.conservativeResize(numFullV,3);
    fullT.conservativeResize(numFullT,3);
    int currVIndex=0, currFIndex=0;
    
    //integrating velocity, position and orientation from forces and previous states
    for (int i=0;i<rigidObjects.size();i++)
      rigidObjects[i].integrate(timeStep);
    
    //detecting and handling collisions when found
    //This is done exhaustively: checking every two objects in the scene.
    double depth;
    RowVector3d contactNormal, penPosition;
    for (int i=0;i<rigidObjects.size();i++)
      for (int j=i+1;j<rigidObjects.size();j++)
        if (rigidObjects[i].isCollide(rigidObjects[j],depth, contactNormal, penPosition))
          handleCollision(rigidObjects[i], rigidObjects[j],depth, contactNormal.normalized(), penPosition, CRCoeff);
    
    
    
    //Code for updating visualization meshes
    for (int i=0;i<rigidObjects.size();i++){
      fullT.block(currFIndex, 0, rigidObjects[i].T.rows(),3)=rigidObjects[i].T.array()+currVIndex;   //need to advance the indices, because every object is indexed independently
      fullV.block(currVIndex, 0, rigidObjects[i].currV.rows(),3)=rigidObjects[i].currV;
      currFIndex+=rigidObjects[i].T.rows();
      currVIndex+=rigidObjects[i].currV.rows();
    }
    currTime+=timeStep;
  }
  
  //loading a scene from the scene .txt files
  //you do not need to update this function
  bool loadScene(const std::string path, const std::string sceneFileName){
    
    ifstream sceneFileHandle;
    sceneFileHandle.open(path+std::string("/")+sceneFileName);
    if (!sceneFileHandle.is_open())
      return false;
    int numofObjects;
    
    currTime=0;
    numFullT=0;
    numFullV=0;
    sceneFileHandle>>numofObjects;
    for (int i=0;i<numofObjects;i++){
      MatrixXi objT;
      MatrixXd objV;
      std::string OFFFileName;
      bool isFixed;
      double density;
      RowVector3d COM;
      RowVector4d orientation;
      sceneFileHandle>>OFFFileName>>density>>isFixed>>COM(0)>>COM(1)>>COM(2)>>orientation(0)>>orientation(1)>>orientation(2)>>orientation(3);
      orientation.normalize();
      igl::readOFF(path+std::string("/")+OFFFileName,objV,objT);
      addRigidObject(objV,objT,density, isFixed, COM, orientation);
      cout << "COM: " << COM <<endl;
      cout << "orientation: " << orientation <<endl;
    }
    return true;
  }
  
  void initalizeRandomVelocities(float bounds) {

	  int i = 0;
	  for (auto &ro : rigidObjects) {
		  if (i != rigidObjects.size() - 1) {
			  ro.comVelocity(0, 0) = dRand(-bounds, bounds);
			  ro.comVelocity(0, 1) = dRand(-bounds, bounds * 0.5f);
			  ro.comVelocity(0, 2) = dRand(-bounds, bounds);
		  }
		  i++;
	  }
  }
  
  Scene(){}
  ~Scene(){}
};


/*****************************Auxiliary functions for collision detection. Do not need updating********************************/

/** Support function for libccd*/
void support(const void *_obj, const ccd_vec3_t *_d, ccd_vec3_t *_p)
{
  // assume that obj_t is user-defined structure that holds info about
  // object (in this case box: x, y, z, pos, quat - dimensions of box,
  // position and rotation)
  //std::cout<<"calling support"<<std::endl;
  RigidObject *obj = (RigidObject *)_obj;
  RowVector3d p;
  RowVector3d d;
  for (int i=0;i<3;i++)
    d(i)=_d->v[i]; //p(i)=_p->v[i];
  
  
  d.normalize();
  //std::cout<<"d: "<<d<<std::endl;
  
  int maxVertex=-1;
  int maxDotProd=-32767.0;
  for (int i=0;i<obj->currV.rows();i++){
    double currDotProd=d.dot(obj->currV.row(i)-obj->COM);
    if (maxDotProd < currDotProd){
      maxDotProd=currDotProd;
      //std::cout<<"maxDotProd: "<<maxDotProd<<std::endl;
      maxVertex=i;
    }
    
  }
  //std::cout<<"maxVertex: "<<maxVertex<<std::endl;
  
  for (int i=0;i<3;i++)
    _p->v[i]=obj->currV(maxVertex,i);
  
  //std::cout<<"end support"<<std::endl;
}

void stub_dir(const void *obj1, const void *obj2, ccd_vec3_t *dir)
{
  dir->v[0]=1.0;
  dir->v[1]=0.0;
  dir->v[2]=0.0;
}

void center(const void *_obj,ccd_vec3_t *center)
{
  RigidObject *obj = (RigidObject *)_obj;
  for (int i=0;i<3;i++)
    center->v[i]=obj->COM(i);
}






#endif
