/*
 * MIT License
 * 
 * Copyright (c) 2020 Jenna Reher (jreher@caltech.edu)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/

#include <cassie_estimation/kinematics_hip_velocity_estimator.hpp>
#include <unsupported/Eigen/EulerAngles>
#include <cassie_common_toolbox/geometry.hpp>


KinematicsHipVelocityEstimator::KinematicsHipVelocityEstimator(ros::NodeHandle &nh, cassie_model::Cassie &robot, bool do_update_robot) {
    this->robot = &robot;
    this->do_update_robot = do_update_robot;
    this->config.paramChecker.init(nh.getNamespace() + "/kinematics_velocity");
    this->reconfigureService = nh.advertiseService("reconfigure_kinematics_estimation", &KinematicsHipVelocityEstimator::reconfigure, this);
    this->reconfigure();
    this->initialize();
}

void KinematicsHipVelocityEstimator::initialize() {
    this->reset();
}

void KinematicsHipVelocityEstimator::Config::reconfigure() {
    paramChecker.checkAndUpdate("/cassie/dt", this->control_rate);
    paramChecker.checkAndUpdate("lowpass_vx_dt_cutoff", this->lowpass_vx_dt_cutoff);
    paramChecker.checkAndUpdate("lowpass_vy_dt_cutoff", this->lowpass_vy_dt_cutoff);
    paramChecker.checkAndUpdate("lowpass_vz_dt_cutoff", this->lowpass_vz_dt_cutoff);
}

void KinematicsHipVelocityEstimator::reconfigure() {
    std::cout << "Polling rosparams under: " << this->config.paramChecker.node.getNamespace() << std::endl;
    this->config.reconfigure();
    this->lpVX.reconfigure(this->config.control_rate, this->config.lowpass_vx_dt_cutoff);
    this->lpVY.reconfigure(this->config.control_rate, this->config.lowpass_vy_dt_cutoff);
    this->lpVZ.reconfigure(this->config.control_rate, this->config.lowpass_vz_dt_cutoff);
}

bool KinematicsHipVelocityEstimator::reconfigure(std_srvs::Empty::Request &req, std_srvs::Empty::Response &res) {
    this->reconfigure();
    return true;
}

void KinematicsHipVelocityEstimator::reset() {
    this->velocity.setZero();
    this->lpVX.reset(0.);
    this->lpVY.reset(0.);
    this->lpVZ.reset(0.);
}

Vector3d KinematicsHipVelocityEstimator::getValue() {
    return this->velocity;
}

void KinematicsHipVelocityEstimator::update() {
    // Get system state
    VectorXd denc_left(8), denc_right(8);
    denc_left << this->robot->dq(LeftHipYaw),
                 this->robot->dq(LeftHipRoll),
                 this->robot->dq(LeftHipPitch),
                 this->robot->dq(LeftKneePitch),
                 this->robot->dq(LeftShinPitch),
                 this->robot->dq(LeftTarsusPitch),
                 this->robot->dq(LeftHeelSpring),
                 this->robot->dq(LeftFootPitch);
    denc_right << this->robot->dq(RightHipYaw),
                  this->robot->dq(RightHipRoll),
                  this->robot->dq(RightHipPitch),
                  this->robot->dq(RightKneePitch),
                  this->robot->dq(RightShinPitch),
                  this->robot->dq(RightTarsusPitch),
                  this->robot->dq(RightHeelSpring),
                  this->robot->dq(RightFootPitch);

    MatrixXd fk_left(5,1), fk_right(5,1);
    MatrixXd J_fk_left(3,22), J_fk_right(3,22), temp_J_fk_left(5,22), temp_J_fk_right(5,22);
    SymFunction::p_leftSole_constraint(fk_left, this->robot->q); //
    SymFunction::p_rightSole_constraint(fk_right, this->robot->q);
    SymFunction::J_leftSole_constraint(temp_J_fk_left, this->robot->q);
    SymFunction::J_rightSole_constraint(temp_J_fk_right, this->robot->q);
    Vector3d velocity_left = -(temp_J_fk_left.block(0,BasePosX, 3,3)).inverse()*temp_J_fk_left.block(0,BaseRotX, 3,3)*this->robot->dq.segment(BaseRotX,3) - temp_J_fk_left.block(0,LeftHipRoll, 3,8) * denc_left;
    Vector3d velocity_right = -(temp_J_fk_right.block(0,BasePosX, 3,3)).inverse()*temp_J_fk_right.block(0,BaseRotX, 3,3)*this->robot->dq.segment(BaseRotX,3) - temp_J_fk_right.block(0,RightHipRoll, 3,8) * denc_right;
    
    // Filter
    Vector3d velocity;
    velocity.setZero();
    if (this->robot->leftContact >= 0.99 && this->robot->rightContact >= 0.99)
        velocity = (this->robot->leftContact * velocity_left + this->robot->rightContact * velocity_right) / (this->robot->leftContact + this->robot->rightContact);
    else if (this->robot->leftContact >= 0.99)
        velocity = velocity_left;
    else if (this->robot->rightContact >= 0.99)
        velocity = velocity_right;
    else
        velocity << 0.9*this->lpVX.getValue(), 0.9*this->lpVY.getValue(), 0.9*this->lpVZ.getValue(); // Decay
    this->lpVX.update(velocity(0));
    this->lpVY.update(velocity(1));
    this->lpVZ.update(velocity(2));

    // Return and store
    this->velocity << this->lpVX.getValue(), this->lpVY.getValue(), this->lpVZ.getValue();
    if (this->do_update_robot){
        this->robot->dq(BasePosX) = this->velocity(0);
        this->robot->dq(BasePosY) = this->velocity(1);
        this->robot->dq(BasePosZ) = this->velocity(2);
    }
}
