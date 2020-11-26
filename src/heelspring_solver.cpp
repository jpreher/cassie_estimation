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

#include <cassie_estimation/heelspring_solver.hpp>
#include <ros/ros.h>

using namespace cassie_model;

void HeelspringSolver::Config::reconfigure() {
    paramChecker.checkAndUpdate("zero_springs_on_startup", this->zero_springs_on_startup);
    paramChecker.checkAndUpdate("num_zeroing_samples", this->num_zeroing_samples);
}

void HeelspringSolver::reconfigure() {
    this->config.reconfigure();
    this->memory.reset();
    this->cache.reset();
}

bool HeelspringSolver::reconfigure(std_srvs::Empty::Request &req, std_srvs::Empty::Response &res) {
    this->reconfigure();
    return true;
}

HeelspringSolver::HeelspringSolver(ros::NodeHandle &nh, cassie_model::Cassie &robot) {
    this->robot = &robot;
    this->config.paramChecker.init(nh.getNamespace() + "/heelspring_solver");

    this->cache.init();
    this->memory.init();

    // Service for calling reconfigures
    this->reconfigureService = nh.advertiseService("reconfigure_heelspring_solver", &HeelspringSolver::reconfigure, this);
    this->config.reconfigure();
}

void HeelspringSolver::update() {
    // Run the inverse kinematics
    VectorXd x(2);

    if (this->config.zero_springs_on_startup && (this->memory.nZeroSamplesAccrued < this->config.num_zeroing_samples)) {
        // Zero springs
        Vector4d temp;
        temp(0) = 0.0;
        temp(1) = 0.0;
        if ( this->InverseKinematics(x) >= 0 ) {
            // Increment the number of samples
            ++this->memory.nZeroSamplesAccrued;

            // Assign heel spring values
            temp(2) = x(0);
            temp(3) = x(1);

            // Assign the offset
            this->memory.spring_offset += ( temp - this->memory.spring_offset ) / this->memory.nZeroSamplesAccrued;
        }
    } else {
        // Update IK
        if ( this->InverseKinematics(x) >= 0 ) {
            // Valid Solution
            if ((fabs(x(0)) > 0.15) || (fabs(x(1)) > 0.15)) {
                // Run again
                x << this->memory.prev_x;
                this->memory.prev_x.setZero();
                if (!this->InverseKinematics(x)) {
                    // Failed a second time, just use last value
                    x << this->memory.prev_x;
                }
            }

            // Assign position to robot
            this->robot->q(LeftHeelSpring)  = x(0) - this->memory.spring_offset(2);
            this->robot->q(RightHeelSpring) = x(1) - this->memory.spring_offset(3);

            // Compute the velocity
            VectorXd dq_meas(6), dqAch(2);
            MatrixXd J(2,22), J_ach(2,2), J_meas(2,6);
            SymFunction::J_achilles_constraint(J, this->robot->q);
            J_ach << J.col(LeftHeelSpring), J.col(RightHeelSpring);
            J_meas << J.col(LeftKneePitch),
                      J.col(LeftShinPitch),
                      J.col(LeftTarsusPitch),
                      J.col(RightKneePitch),
                      J.col(RightShinPitch),
                      J.col(RightTarsusPitch);
            dq_meas << this->robot->dq(LeftKneePitch),
                       this->robot->dq(LeftShinPitch),
                       this->robot->dq(LeftTarsusPitch),
                       this->robot->dq(RightKneePitch),
                       this->robot->dq(RightShinPitch),
                       this->robot->dq(RightTarsusPitch);
            dqAch = - J_ach.inverse() * J_meas * dq_meas;

            // Assign velocity to robot
            this->robot->dq(LeftHeelSpring)  = dqAch(0);
            this->robot->dq(RightHeelSpring) = dqAch(1);
        }
    }
}

void HeelspringSolver::IKfunction(Eigen::VectorXd &x, Eigen::VectorXd &fvec) {
    this->cache.q(LeftHeelSpring) = x(0);
    this->cache.q(RightHeelSpring) = x(1);
    {
        VectorWrap fvec_(fvec);
        SymFunction::p_achilles_constraint(fvec_, this->cache.q);
    }
}

void HeelspringSolver::J_IKfunction(Eigen::VectorXd &x, Eigen::MatrixXd &J) {
    this->cache.q(LeftHeelSpring) = x(0);
    this->cache.q(RightHeelSpring) = x(1);

    MatrixXd J_temp(2,22);
    SymFunction::J_achilles_constraint(J_temp,this->cache.q);

    J.setZero();
    J(0,0) = J_temp(0,LeftHeelSpring);
    J(1,1) = J_temp(1,RightHeelSpring);
}

int HeelspringSolver::InverseKinematics(Eigen::VectorXd &x) {
    int iteration_limit = 10;
    double xtol = 1e-16;

    this->cache.q << this->robot->q;
    for (int i=0; i<6; ++i)
        this->cache.q(i) = 0.0;

    Eigen::VectorXd F(2);
    Eigen::MatrixXd J(2,2);
    double f_d = 0.0;
    x << this->memory.prev_x;

    for (int i=0; i<iteration_limit; ++i){
        // Reinitialize the state at the current guess
        this->IKfunction(x, F);

        // Check for completion
        f_d = 0.5 * F.transpose() * F;
        if (f_d < xtol) {
            this->memory.prev_x << x;
            return i;
        }

        // Compute the Jacobian
        this->J_IKfunction(x, J);

        // Perform the update
        x = x - (J.transpose() * (J * J.transpose()).inverse()) * F;
    }
    ROS_WARN("HEELSPRING IK DID NOT CONVERGE");
    return -1;
}
