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

#include <cassie_estimation/contact_classifier.hpp>
#include <control_utilities/limits.hpp>

using namespace cassie_model;

ContactClassifier::ContactClassifier(ros::NodeHandle &nh, cassie_model::Cassie &robot, double dt)
    : LowPassLeft(dt, 1/300), LowPassRight(dt, 1/300)
{
    this->config.paramChecker.init(nh.getNamespace() + "/contact_classifier");

    // Pointer to the robot
    this->robot = &robot;

    this->config.init();

    // Service for calling reconfigures
    this->reconfigureService = nh.advertiseService("reconfigure_contact_classifier", &ContactClassifier::reconfigure, this);
    this->reconfigure();

    this->grf.resize(6);
}

double sigmoid(double s, double A, double B, double power) {
    return 1.0 / ( 1.0 + exp(-A * pow((s - B),power)) );
}

void ContactClassifier::Config::init() {
    this->sigmoid_A = 0.35;
    this->sigmoid_B = 30;
    this->sigmoid_power = 1;
    this->dt = 0.0005;
    this->lowpass_dt_cutoff = 1/300.;
}

void ContactClassifier::Config::reconfigure() {
    std::cout << "Polling rosparams under: " << this->paramChecker.node.getNamespace() << std::endl;

    paramChecker.checkAndUpdate("/cassie/dt", this->dt);
    paramChecker.checkAndUpdate("use_sigmoid", this->use_sigmoid);
    paramChecker.checkAndUpdate("sigmoid_A", this->sigmoid_A);
    paramChecker.checkAndUpdate("sigmoid_B", this->sigmoid_B);
    paramChecker.checkAndUpdate("sigmoid_power", this->sigmoid_power);
    paramChecker.checkAndUpdate("linear_lb", this->linear_lb);
    paramChecker.checkAndUpdate("linear_ub", this->linear_ub);
    paramChecker.checkAndUpdate("lowpass_dt_cutoff", this->lowpass_dt_cutoff);
}

void ContactClassifier::reconfigure() {
    this->config.reconfigure();
    this->LowPassLeft.reconfigure(this->config.dt, this->config.lowpass_dt_cutoff);
    this->LowPassRight.reconfigure(this->config.dt, this->config.lowpass_dt_cutoff);
}

bool ContactClassifier::reconfigure(std_srvs::Empty::Request &req, std_srvs::Empty::Response &res) {
    this->reconfigure();
    return true;
}

void ContactClassifier::update() {
    MatrixXd Jl(3,7), Jr(3,7), J(6,14), JT(14,6);
    VectorXd tau(14);

    this->robot->kinematics.computeConstrainedToeJacobian(this->robot->q, Jl, Jr);

    J.setZero();
    J.block(0,0,3,7) << Jl;
    J.block(3,7,3,7) << Jr;

    tau << this->robot->torque(0), // lhr
           this->robot->torque(1), // lhy
           this->robot->torque(2), // lhp
           this->robot->torque(3), // lkp
           - this->robot->q(CassieStateEnum::LeftShinPitch) * 2300., // lsp
           - this->robot->q(CassieStateEnum::LeftHeelSpring) * 2000., // lhsp
           this->robot->torque(4), // ltp
           this->robot->torque(5), // rhr
           this->robot->torque(6), // rhy
           this->robot->torque(7), // rhp
           this->robot->torque(8), // rkp
           - this->robot->q(CassieStateEnum::RightShinPitch) * 2300., // rsp
           - this->robot->q(CassieStateEnum::RightHeelSpring) * 2000., // rhsp
           this->robot->torque(9);

    // Compute the quasi-static grf estimate in robot body frame
    this->grf = - (J.transpose()).completeOrthogonalDecomposition().solve(tau);

    // Rotate into the world
    Eigen::Matrix3d Rx, Ry, R;
    Rx << 1.,0.,0.,
          0., cos(this->robot->q(BaseRotX)), -sin(this->robot->q(BaseRotX)),
          0., sin(this->robot->q(BaseRotX)), cos(this->robot->q(BaseRotX));
    Ry << cos(this->robot->q(BaseRotY)), 0, sin(this->robot->q(BaseRotY)),
          0.,1.,0.,
          -sin(this->robot->q(BaseRotY)), 0, cos(this->robot->q(BaseRotY));
    R << Ry * Rx;

    grf.segment(0,3) = R * grf.segment(0,3);
    grf.segment(3,3) = R * grf.segment(3,3);

    // Update vertical grf lowpass
    this->LowPassLeft.update(this->grf(2));
    this->LowPassRight.update(this->grf(5));

    if (this->config.use_sigmoid) {
        // Sigmoid classifier
        this->robot->leftContact  = sigmoid(this->LowPassLeft.getValue(), this->config.sigmoid_A, this->config.sigmoid_B, this->config.sigmoid_power);
        this->robot->rightContact = sigmoid(this->LowPassRight.getValue(), this->config.sigmoid_A, this->config.sigmoid_B, this->config.sigmoid_power);
    } else {
        // Linear classifier
        this->robot->leftContact = control_utilities::clamp((this->LowPassLeft.getValue() - this->config.linear_lb) / (this->config.linear_ub - this->config.linear_lb), 0, 1);
        this->robot->rightContact = control_utilities::clamp((this->LowPassRight.getValue() - this->config.linear_lb) / (this->config.linear_ub - this->config.linear_lb), 0, 1);
    }



}
