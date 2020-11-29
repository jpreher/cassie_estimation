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

#include <cassie_estimation/rigidtarsus_solver.hpp>
#include <ros/ros.h>

using namespace cassie_model;

void RigidTarsusSolver::reconfigure() {
    this->memory.reset();
    this->cache.reset();
}

RigidTarsusSolver::RigidTarsusSolver(cassie_model::Cassie &robot) {
    this->robot = &robot;
    this->cache.init();
    this->memory.init();
}

void RigidTarsusSolver::reset() {
    this->cache.reset();
    this->memory.reset();
}

void RigidTarsusSolver::update() {
    // Run the inverse kinematics
    VectorXd x(2);
    if (!this->memory.is_initialized)
        this->memory.prev_x << this->robot->q(LeftTarsusPitch), this->robot->q(RightTarsusPitch);

    // Update IK
    x << this->robot->q(LeftTarsusPitch), this->robot->q(RightTarsusPitch); // start from measured
    this->InverseKinematics(x);

    // Assign position to robot
    this->memory.tar_sol << x;

    // Compute the velocity
    VectorXd dq_meas(6), dqtar(2);
    MatrixXd J(2,22), J_tar(2,2), J_meas(2,6);
    SymFunction::J_achilles_constraint(J, this->cache.q);
    J_tar << J.col(LeftTarsusPitch), J.col(RightTarsusPitch);
    J_meas << J.col(LeftKneePitch),
              J.col(LeftShinPitch),
              J.col(LeftHeelSpring),
              J.col(RightKneePitch),
              J.col(RightShinPitch),
              J.col(RightHeelSpring);
    dq_meas << this->robot->dq(LeftKneePitch),
               0., // Left Shin
               0., // Left Heel
               this->robot->dq(RightKneePitch),
               0., // Right shin
               0.; // Right heel
    this->memory.dtar_sol = - J_tar.inverse() * J_meas * dq_meas;
}

void RigidTarsusSolver::IKfunction(Eigen::VectorXd &x, Eigen::VectorXd &fvec) {
    this->cache.q(LeftTarsusPitch) = x(0);
    this->cache.q(RightTarsusPitch) = x(1);
    {
        VectorWrap fvec_(fvec);
        SymFunction::p_achilles_constraint(fvec_, this->cache.q);
    }
}

void RigidTarsusSolver::J_IKfunction(Eigen::VectorXd &x, Eigen::MatrixXd &J) {
    this->cache.q(LeftTarsusPitch) = x(0);
    this->cache.q(RightTarsusPitch) = x(1);

    MatrixXd J_temp(2,22);
    SymFunction::J_achilles_constraint(J_temp,this->cache.q);

    J.setZero();
    J(0,0) = J_temp(0,LeftTarsusPitch);
    J(1,1) = J_temp(1,RightTarsusPitch);
}

int RigidTarsusSolver::InverseKinematics(Eigen::VectorXd &x) {
    int iteration_limit = 10;
    double xtol = 1e-16;

    this->cache.q << this->robot->q;
    for (int i=0; i<6; ++i)
        this->cache.q(i) = 0.0;
    this->cache.q(LeftShinPitch) = 0.;
    this->cache.q(LeftHeelSpring) = 0.;
    this->cache.q(RightShinPitch) = 0.;
    this->cache.q(RightHeelSpring) = 0.;

    Eigen::VectorXd F(2);
    Eigen::MatrixXd J(2,2);
    double f_d = 0.0;

    for (int i=0; i<iteration_limit; ++i) {
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
    ROS_WARN("TARSUS IK DID NOT CONVERGE");
    return -1;
}
