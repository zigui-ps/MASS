#include "Character.h"
#include "BVH.h"
#include "DARTHelper.h"
#include "Muscle.h"
#include <tinyxml.h>
using namespace dart;
using namespace dart::dynamics;
using namespace MASS;
Character::
Character()
	:mSkeleton(nullptr),mBVH(nullptr),mTc(Eigen::Isometry3d::Identity())
{

}

void
Character::
LoadSkeleton(const std::string& path,bool create_obj)
{
	mSkeleton = BuildFromFile(path,create_obj);
	std::map<std::string,std::string> bvh_map;
	TiXmlDocument doc;
	doc.LoadFile(path);
	TiXmlElement *skel_elem = doc.FirstChildElement("Skeleton");

	for(TiXmlElement* node = skel_elem->FirstChildElement("Node");node != nullptr;node = node->NextSiblingElement("Node"))
	{
		if(node->Attribute("endeffector")!=nullptr)
		{
			std::string ee =node->Attribute("endeffector");
			if(ee == "True")
			{
				mEndEffectors.push_back(mSkeleton->getBodyNode(std::string(node->Attribute("name"))));
			}
		}
		TiXmlElement* joint_elem = node->FirstChildElement("Joint");
		if(joint_elem->Attribute("bvh")!=nullptr)
		{
			bvh_map.insert(std::make_pair(node->Attribute("name"),joint_elem->Attribute("bvh")));
		}
	}
	
	mBVH = new BVH(mSkeleton,bvh_map);
}
void
Character::
LoadMuscles(const std::string& path)
{
	TiXmlDocument doc;
	if(!doc.LoadFile(path)){
		std::cout << "Can't open file : " << path << std::endl;
		return;
	}

	TiXmlElement *muscledoc = doc.FirstChildElement("Muscle");
	for(TiXmlElement* unit = muscledoc->FirstChildElement("Unit");unit!=nullptr;unit = unit->NextSiblingElement("Unit"))
	{
		std::string name = unit->Attribute("name");
		double f0 = std::stod(unit->Attribute("f0"));
		double lm = std::stod(unit->Attribute("lm"));
		double lt = std::stod(unit->Attribute("lt"));
		double pa = std::stod(unit->Attribute("pen_angle"));
		double lmax = std::stod(unit->Attribute("lmax"));
		mMuscles.push_back(new Muscle(name,f0,lm,lt,pa,lmax));
		int num_waypoints = 0;
		for(TiXmlElement* waypoint = unit->FirstChildElement("Waypoint");waypoint!=nullptr;waypoint = waypoint->NextSiblingElement("Waypoint"))	
			num_waypoints++;
		int i = 0;
		for(TiXmlElement* waypoint = unit->FirstChildElement("Waypoint");waypoint!=nullptr;waypoint = waypoint->NextSiblingElement("Waypoint"))	
		{
			std::string body = waypoint->Attribute("body");
			Eigen::Vector3d glob_pos = string_to_vector3d(waypoint->Attribute("p"));
			if(i==0||i==num_waypoints-1)
			// if(true)
				mMuscles.back()->AddAnchor(mSkeleton->getBodyNode(body),glob_pos);
			else
				mMuscles.back()->AddAnchor(mSkeleton,mSkeleton->getBodyNode(body),glob_pos,2);

			i++;
		}
		mMuscles.back()->SetMuscle();
	}
	

}
void
Character::
LoadBVH(const std::string& path,bool cyclic)
{
	if(mBVH ==nullptr){
		std::cout<<"Initialize BVH class first"<<std::endl;
		return;
	}
	mBVH->Parse(path,cyclic);
}
void
Character::
Reset()
{
	mTc = mBVH->GetT0();
	mTc.translation()[1] = 0.0;
}
void
Character::
SetPDParameters(double kp, double kv)
{
	int dof = mSkeleton->getNumDofs();
	mKp = Eigen::VectorXd::Constant(dof,kp);	
	mKv = Eigen::VectorXd::Constant(dof,kv);	
}
Eigen::VectorXd
Character::
GetSPDForces(const Eigen::VectorXd& p_desired, const Eigen::VectorXd& v_desired)
{
	auto &skel = mSkeleton;
	Eigen::VectorXd q = skel->getPositions();
	Eigen::VectorXd dq = skel->getVelocities();
	double dt = skel->getTimeStep();
	Eigen::MatrixXf M = (skel->getMassMatrix() + Eigen::MatrixXd(dt * mKv.asDiagonal())).cast<float>();
	Eigen::MatrixXd M_inv = M.inverse().cast<double>();

	Eigen::VectorXd p_d = q + dq*dt - p_desired;
	// clamping radians to [-pi, pi], only for ball joints
	// TODO : make it for all type joints
	/*
	p_d.segment<6>(0) = Eigen::VectorXd::Zero(6);
	for (int i = 6; i < (int)skel->getNumDofs(); i += 3)
	{
		Eigen::Quaterniond q_s = Basic::DARTPositionToQuaternion(q.segment<3>(i));
		Eigen::Quaterniond dq_s = Basic::DARTPositionToQuaternion(dt * (dq.segment<3>(i)));
		Eigen::Quaterniond q_d_s = Basic::DARTPositionToQuaternion(p_desired.segment<3>(i));

		Eigen::Quaterniond p_d_s = q_d_s.inverse() * q_s * dq_s;

		Eigen::Vector3d v = Basic::QuaternionToDARTPosition(p_d_s);
		p_d.segment<3>(i) = v;
	} // */
	Eigen::VectorXd p_diff = -mKp.cwiseProduct(p_d);
	Eigen::VectorXd v_diff = -mKv.cwiseProduct(dq - v_desired);
	Eigen::VectorXd qddot = M_inv * (-skel->getCoriolisAndGravityForces() +
		p_diff + v_diff + skel->getConstraintForces());

	Eigen::VectorXd tau = p_diff + v_diff - dt * mKv.cwiseProduct(qddot);
	tau.segment<6>(0) = Eigen::VectorXd::Zero(6);

	return tau;
}
Eigen::VectorXd
Character::
GetTargetPositions(double t,double dt)
{
	// std::cout<<"GetTargetPositions"<<std::endl;
	Eigen::VectorXd p = mBVH->GetMotion(t);	
	Eigen::Isometry3d T_current = dart::dynamics::FreeJoint::convertToTransform(p.head<6>());
	T_current = mBVH->GetT0().inverse()*T_current;
	Eigen::Isometry3d T_head = mTc*T_current;
	Eigen::Vector6d p_head = dart::dynamics::FreeJoint::convertToPositions(T_head);
	p.head<6>() = p_head;
	
	
	if(mBVH->IsCyclic())
	{
		double t_mod = std::fmod(t, mBVH->GetMaxTime());
		t_mod = t_mod/mBVH->GetMaxTime();

		double r = 0.95;
		if(t_mod>r)
		{
			double ratio = 1.0/(r-1.0)*t_mod - 1.0/(r-1.0);
			Eigen::Isometry3d T01 = mBVH->GetT1()*(mBVH->GetT0().inverse());
			double delta = T01.translation()[1];
			delta *= ratio;
			p[5] += delta;
		}



		double tdt_mod = std::fmod(t+dt, mBVH->GetMaxTime());
		if(tdt_mod-dt<0.0){
			Eigen::Isometry3d T01 = mBVH->GetT1()*(mBVH->GetT0().inverse());
			Eigen::Vector3d p01 = dart::math::logMap(T01.linear());
			p01[0] =0.0;
			p01[2] =0.0;
			T01.linear() = dart::math::expMapRot(p01);

			mTc = T01*mTc;
			mTc.translation()[1] = 0.0;
		}
	}
	

	return p;
}

std::pair<Eigen::VectorXd,Eigen::VectorXd>
Character::
GetTargetPosAndVel(double t,double dt)
{
	Eigen::VectorXd p = this->GetTargetPositions(t,dt);
	Eigen::Isometry3d Tc = mTc;
	Eigen::VectorXd p1 = this->GetTargetPositions(t+dt,dt);
	mTc = Tc;

	return std::make_pair(p,(p1-p)/dt);
}