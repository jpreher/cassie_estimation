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

#include <cassie_estimation/contact_ekf.hpp>
#include <unsupported/Eigen/EulerAngles>

void contact_ekf::Cache::init() {
    this->Fc.resize(23,23);
    this->Fk.resize(23,23);
    this->Gc.resize(23,20);
    this->Hk.resize(8,23);
    this->Rk.resize(8,8);
    this->Sk.resize(8,8);
    this->K.resize(23,8);
    this->Qscaled.resize(23,23);
    
    this->reset();
}

void contact_ekf::Cache::reset() {
    this->Fc.setZero();
    this->Fk.setZero();
    this->Gc.setZero();
    this->Hk.setZero();
    this->Rk.setZero();
    this->Sk.setZero();
    this->K.setZero();
    this->Qscaled.setZero();
}

void contact_ekf::Memory::init() {
    this->X.resize(3,10);
    this->P.resize(23,23);
    this->encoders_prev.resize(14);
    this->dencoders_prev.resize(14);

    this->reset();
}

void contact_ekf::Memory::reset() {
    this->timer.reset();
    this->t_last = 0.;
    this->filter_enabled = false;
    this->bias_initialized = false;
    this->contact_mode = 0;

    this->X.setZero();
    this->P.setIdentity();

    this->ba0.setZero();
    this->bg0.setZero();

    this->w_prev.setZero();
    this->a_prev.setZero();
    this->encoders_prev.setZero();
    this->dencoders_prev.setZero();
    this->contact_prev.setZero();

     this->v_init.setZero();
}

void contact_ekf::Config::init() {
    this->P_prior.resize(23,23);
    this->Qcontinuous.resize(20,20);
}

void contact_ekf::Config::reconfigure() {
    this->g << 0,0,-9.81;
    this->paramChecker.checkAndUpdate("/cassie/dt", this->dt);

    this->paramChecker.checkAndUpdate("gyro_noise_std",        this->gyro_noise_std);
    this->paramChecker.checkAndUpdate("gyro_bias_noise_std",   this->gyro_bias_noise_std);
    this->paramChecker.checkAndUpdate("accel_noise_std",       this->accel_noise_std);
    this->paramChecker.checkAndUpdate("accel_bias_noise_std",  this->accel_bias_noise_std);
    this->paramChecker.checkAndUpdate("contact_noise_std",     this->contact_noise_std);
    this->paramChecker.checkAndUpdate("contact_yaw_noise_std", this->contact_yaw_noise_std);
    this->paramChecker.checkAndUpdate("dcontact_noise_std",    this->dcontact_noise_std);
    this->paramChecker.checkAndUpdate("encoder_noise_std",     this->encoder_noise_std);

    this->paramChecker.checkAndUpdate("position_std",    this->position_std);
    this->paramChecker.checkAndUpdate("velocity_std",    this->velocity_std);
    this->paramChecker.checkAndUpdate("contact_std",     this->contact_std);
    this->paramChecker.checkAndUpdate("contact_yaw_std", this->contact_yaw_std);

    this->paramChecker.checkAndUpdate("prior_base_rotation_std",      this->prior_base_rotation_std);
    this->paramChecker.checkAndUpdate("prior_base_position_std",      this->prior_base_position_std);
    this->paramChecker.checkAndUpdate("prior_base_velocity_std",      this->prior_base_velocity_std);
    this->paramChecker.checkAndUpdate("prior_contact_position_std",   this->prior_contact_position_std);
    this->paramChecker.checkAndUpdate("prior_gyro_bias_std",          this->prior_gyro_bias_std);
    this->paramChecker.checkAndUpdate("prior_accel_bias_std",         this->prior_accel_bias_std);
    this->paramChecker.checkAndUpdate("prior_forward_kinematics_std", this->prior_forward_kinematics_std);
    this->paramChecker.checkAndUpdate("prior_contact_yaw_std",        this->prior_contact_yaw_std);

    this->paramChecker.checkAndUpdate("apply_post_filter",       this->apply_post_filter);
    this->paramChecker.checkAndUpdate("post_filter_dt_cutoff_x", this->post_filter_dt_cutoff_x);
    this->paramChecker.checkAndUpdate("post_filter_dt_cutoff_y", this->post_filter_dt_cutoff_y);
    this->paramChecker.checkAndUpdate("post_filter_dt_cutoff_z", this->post_filter_dt_cutoff_z);

    // Assign covariance values
    // Create the continuous process covariance matrix
    this->Qcontinuous.resize(20,20); this->Qcontinuous.setZero();
    this->Qcontinuous.diagonal() << pow(this->accel_noise_std,2)      * VectorXd::Ones(3),
                                    pow(this->gyro_noise_std,2)       * VectorXd::Ones(3),
                                    pow(this->accel_bias_noise_std,2) * VectorXd::Ones(3),
                                    pow(this->gyro_bias_noise_std,2)  * VectorXd::Ones(3),
                                    pow(this->contact_std,2)          * VectorXd::Ones(3),
                                    pow(this->contact_std,2)          * VectorXd::Ones(3),
                                    pow(this->contact_yaw_std,2)      * VectorXd::Ones(2);

    // Build prior
    this->P_prior.setZero();
    this->P_prior.diagonal() <<
            pow(this->prior_base_rotation_std,2)    * VectorXd::Ones(3), // rotation
            pow(this->prior_base_position_std,2)    * VectorXd::Ones(3), // position
            pow(this->prior_base_velocity_std,2)    * VectorXd::Ones(3), // velocity
            pow(this->prior_gyro_bias_std,2)        * VectorXd::Ones(3), // bias
            pow(this->prior_accel_bias_std,2)       * VectorXd::Ones(3),
            pow(this->prior_contact_position_std,2) * VectorXd::Ones(3), // lf_pos
            pow(this->prior_contact_position_std,2) * VectorXd::Ones(3), // rf_pos
            pow(this->prior_contact_yaw_std,2)      * VectorXd::Ones(2); // lf/rf yaw
}

contact_ekf::contact_ekf(ros::NodeHandle &nh, cassie_model::Cassie &robot, bool do_update_robot) : nh(&nh) {
    // Main constructor
    this->robot = &robot;
    this->cache.init();
    this->memory.init();
    this->config.init();
    this->config.do_update_robot = do_update_robot;

    // Parameter checker
    this->config.paramChecker.init(nh.getNamespace() + "/contact_ekf");

    // Service for calling reconfigures
    this->reconfigureService = nh.advertiseService("reconfigure_contact_ekf", &contact_ekf::reconfigure, this);
    this->reconfigure();
}

bool contact_ekf::reconfigure() {
    std::cout << "Polling rosparams under: " << this->config.paramChecker.node.getNamespace() << std::endl;
    this->config.reconfigure();
    this->lpVX.reconfigure(this->config.dt, this->config.post_filter_dt_cutoff_x);
    this->lpVY.reconfigure(this->config.dt, this->config.post_filter_dt_cutoff_y);
    this->lpVZ.reconfigure(this->config.dt, this->config.post_filter_dt_cutoff_z);

    return true;
}

bool contact_ekf::reconfigure(std_srvs::Empty::Request &req, std_srvs::Empty::Response &res) {
    return this->reconfigure();
}

void contact_ekf::initialize() {

}

void contact_ekf::reset() {
    this->cache.reset();
    this->memory.reset();
    this->lpVX.reset(0.);
    this->lpVY.reset(0.);
    this->lpVZ.reset(0.);
}

void contact_ekf::reset(Vector3d &initial_velocity) {
    this->reset();
    this->memory.v_init << initial_velocity;
}

bool contact_ekf::isEnabled() {
    return this->memory.filter_enabled;
}

void contact_ekf::getValues(Matrix3d &R, Vector3d &p, Vector3d &v, Vector3d &ba, Vector3d &bg, Vector3d &plf, Vector3d &prf, Vector2d &footYaws) {
    unpackState(this->memory.X, R, p, v, ba, bg, plf, prf, footYaws);

    if (this->config.apply_post_filter)
        v << this->lpVX.getValue(), this->lpVY.getValue(), this->lpVZ.getValue();
}

Vector3d contact_ekf::getRawVelocity() {
    Vector3d v;
    v << this->memory.X.col(4);
    return v;
}

void contact_ekf::update(double dt, VectorXd &w, VectorXd &a, VectorXd &encoders, VectorXd &dencoders, VectorXd &contact) {
    // Initialize bias
    // (does nothing if bias is already initialized)
    this->initializeBias(w, a);

    // Check dt and clamp for safety
    if (dt < 0.0001)
        dt = 0.0001;
    if (dt > 0.01)
        dt = 0.01;

    // Initialize filter
    // (does nothing if filter is already initialized)
    if ( (contact(0) >= 0.99) || (contact(1) >= 0.99) ) {
        this->initializeFilter(encoders, contact);
    }

    // Update
    if ( this->memory.filter_enabled  ) {
        // Predict state using IMU and contact measurements
        this->predict_state(dt);

        if ( (contact(0) >= 0.1) || (contact(1) >= 0.1) ) {
            // Update state using forward kinematic measurements
            this->update_forward_kinematics(w, encoders, dencoders, contact);
            this->memory.t_last = ros::Time::now().toSec();
        }
    }

    // Reset the filter if there is no contact for longer than 0.5 seconds
    if ( (ros::Time::now().toSec() - this->memory.t_last) > 0.5 )
        this->reset();

    if (this->config.apply_post_filter) {
        this->lpVX.update(this->memory.X(0,4));
        this->lpVY.update(this->memory.X(1,4));
        this->lpVZ.update(this->memory.X(2,4));
    }

    // Store
    this->memory.w_prev << w;
    this->memory.a_prev << a;
    this->memory.encoders_prev << encoders;
    this->memory.dencoders_prev << dencoders;
    this->memory.contact_prev << contact;

    // Update robot
    if (this->config.do_update_robot) {
        //this->robot->q.block(BasePosX,0,3,1)  << this->memory.X.block(0,4, 3,1);
        //this->robot->dq.block(BasePosX,0,3,1) << this->memory.X.block(0,3, 3,1);
        //this->robot->q.block(BaseRotX,0,3,1) <<
    }
}

void contact_ekf::initializeBias(VectorXd &w, VectorXd &a) {
    // TODO static bias initialization...
    // for now do no bias initialization
    this->memory.ba0.setZero();
    this->memory.ba0 << 0.0, 0.0, -0.20; // Hard coded
    this->memory.bg0.setZero();
    this->memory.bias_initialized = true;
}

void contact_ekf::initializeFilter(VectorXd &encoders, VectorXd &contact) {
    // Attempt to enable filter (successful if enable is true, and
    // at least one foot is on the ground)
    if ( !this->memory.filter_enabled ) {
        // Build a starting guess for the estimator
        Matrix3d Rx, Ry, R;
        Rx << 1.,0.,0.,
              0., cos(this->robot->q(BaseRotX)), -sin(this->robot->q(BaseRotX)),
              0., sin(this->robot->q(BaseRotX)), cos(this->robot->q(BaseRotX));
        Ry << cos(this->robot->q(BaseRotY)), 0., sin(this->robot->q(BaseRotY)),
              0.,1.,0.,
              -sin(this->robot->q(BaseRotY)), 0., cos(this->robot->q(BaseRotY));
        R << Ry * Rx;
      
        Matrix3d Rlf, Rrf;
        Vector2d footYaws = Vector2d::Zero();
        Vector3d pLF, pRF, p_init, v_init;
        p_init.setZero();
        v_init << this->memory.v_init;
        relative_foot_positions(encoders, pLF, pRF, Rlf, Rrf);
        pLF = R * pLF;
        pRF = R * pRF;
        if ( contact(0) >= 0.99 )
            p_init(2) = -pLF(2);
        if ( contact(1) >= 0.99 )
            p_init(2) = -pRF(2); // Z position from ground
        pLF(2) += p_init(2);
        pRF(2) += p_init(2);

        MatrixXd X_init = MatrixXd::Identity(3,10);
        packState(X_init, R, p_init, v_init, this->memory.ba0, this->memory.bg0, pLF, pRF, footYaws);

        this->memory.timer.restart();
        this->memory.t_last = ros::Time::now().toSec() - 0.0005; // Assume 2kHz (better not to assume...)
        this->memory.X << X_init;
        this->memory.P << this->config.P_prior;
        this->memory.filter_enabled = true;
        ROS_INFO("contact_ekf Initialized!");
        std::cout << "Initial ekf state: " << std::endl << X_init << std::endl << std::endl;
    }

    // If filter is disabled, zero everything
    if ( !this->memory.filter_enabled ) {
        this->cache.reset();
        this->memory.reset();
    }
}

void contact_ekf::predict_state(double dt) {
    // Timing
    // double dt = this->config.dt; //ros::Time::now().toSec() - this->memory.t_last;
    // this->memory.t_last = ros::Time::now().toSec();

    // Separate state vector into components
    Matrix3d R;
    Vector2d footYaws;
    Vector3d p, v, ba, bg, pRF, pLF;
    unpackState(this->memory.X, R, p, v, ba, bg, pLF, pRF, footYaws);
    Matrix3d Rt = R.transpose();

    // Bias corrected IMU information
    Vector3d w_k(this->memory.w_prev - bg);
    Vector3d a_k(this->memory.a_prev - ba);

    Vector3d phi = w_k * dt;
    Matrix3d G0 = Gamma(phi,0);
    Matrix3d G1 = Gamma(phi,1);
    Matrix3d G2 = Gamma(phi,2);

    // Strapdown IMU motion model
    Matrix3d R_pred = R * this->Exp_SO3(phi);
    Vector3d p_pred = p + v*dt + (R*G2*a_k + 0.5*this->config.g)*dt*dt;
    Vector3d v_pred = v + (R*G1*a_k + this->config.g) * dt;

    // Foot positions -- keep static during full contact, kinematics during swing
    //Foot Position Dynamics
    Matrix3d RLF_enc, RRF_enc;
    Vector3d pLF_enc(3,1), pRF_enc(3,1);
    relative_foot_positions(this->memory.encoders_prev, pLF_enc, pRF_enc, RLF_enc, RRF_enc);
    Vector3d pRF_off = p_pred + R_pred * pRF_enc;
    Vector3d pLF_off = p_pred + R_pred * pLF_enc;

    // Foot Rotations
    Matrix3d RLF_pred, RRF_pred;
    RLF_pred = R_pred * RLF_enc;
    RRF_pred = R_pred * RRF_enc;

    Vector2d footYaws_off;
    Eigen::EulerAnglesXYZd euler;
    eulerXYZ(RLF_pred, euler);
    footYaws_off(0) = euler.gamma();
    eulerXYZ(RRF_pred, euler);
    footYaws_off(1) = euler.gamma();

    // Get blended prediction of feet
    MatrixXd Qcont(20,20);
    Qcont << this->config.Qcontinuous;
    Vector2d footYaws_pred;
    Vector3d pRF_pred, pLF_pred;
    if (this->memory.contact_prev(1) >= 0.25) {
        pRF_pred = pRF;
        footYaws_pred(1) = footYaws(1);
    } else {
        pRF_pred = pRF_off;
        Qcont.block(15,15,3,3) << 10000.0 * Matrix3d::Identity();

        footYaws_pred(1) = footYaws_off(1);
        Qcont(19,19) = 10000.0;

    }
    Qcont.block(15,15,3,3) = RRF_enc * Qcont.block(15,15,3,3) * RRF_enc.transpose();
    if (this->memory.contact_prev(0) >= 0.25) {
        pLF_pred = pLF;
        footYaws_pred(0) = footYaws(0);
    } else {
        pLF_pred = pLF_off;
        Qcont.block(12,12,3,3) << 10000.0 * Matrix3d::Identity();

        footYaws_pred(0) = footYaws_off(0);
        Qcont(18,18) = 10000.0;
    }
    Qcont.block(12,12,3,3) = RLF_enc * Qcont.block(12,12,3,3) * RLF_enc.transpose();

    // Analytical State transition matrix
    double theta = sqrt(pow(phi(0),2) + pow(phi(1),2) + pow(phi(2),2));
    Matrix3d Psi1, Psi2;
    if (theta < 0.000000001) {
        Psi1 = skew(a_k) * Gamma(-phi, 2)
            + ((sin(theta)-theta*cos(theta))/(pow(theta,3)))*(skew(phi)*skew(a_k))
            - ((cos(2*theta)-4*cos(theta)+3)/(4*pow(theta,4)))*(skew(phi)*skew(a_k)*skew(phi))
            + ((4*sin(theta)+sin(2*theta)-4*theta*cos(theta)-2*theta)/(4*pow(theta,5)))*(skew(phi)*skew(a_k)*skew(phi)*skew(phi))
            + (((theta*theta)-2*theta*sin(theta)-2*cos(theta)+2)/(2*pow(theta,4)))*(skew(phi)*skew(phi)*skew(a_k))
            - ((6*theta-8*sin(theta)+sin(2*theta))/(4*pow(theta,5)))*(skew(phi)*skew(phi)*skew(a_k)*skew(phi))
            + ((2*(theta*theta)-4*theta*sin(theta)-cos(2*theta)+1)/(4*pow(theta,6)))*(skew(phi)*skew(phi)*skew(a_k)*skew(phi)*skew(phi));

        Psi2 = skew(a_k)*Gamma(-phi,3)
            - ((theta*sin(theta)+2*cos(theta)-2)/(pow(theta,4)))*(skew(phi)*skew(a_k))
            - ((6*theta-8*sin(theta)+sin(2*theta))/(8*pow(theta,5)))*(skew(phi)*skew(a_k)*skew(phi))
            - ((2*(theta*theta)+8*theta*sin(theta)+16*cos(theta)+cos(2*theta)-17)/(8*pow(theta,6)))*(skew(phi)*skew(a_k)*skew(phi)*skew(phi))
            + (((pow(theta,3))+6*theta-12*sin(theta)+6*theta*cos(theta))/(6*pow(theta,5)))*(skew(phi)*skew(phi)*skew(a_k))
            - ((6*(pow(theta,2))+16*cos(theta)-cos(2*theta)-15)/(8*pow(theta,6)))*(skew(phi)*skew(phi)*skew(a_k)*skew(phi))
            + ((4*(pow(theta,3))+6*theta-24*sin(theta)-3*sin(2*theta)+24*theta*cos(theta))/(24*pow(theta,7)))*(skew(phi)*skew(phi)*skew(a_k)*skew(phi)*skew(phi));
    } else {
        Psi1 = 0.5 * skew(a_k);
        Psi2 = (1.0/6.0) * skew(a_k);
    }

    MatrixXd eye  = MatrixXd::Identity(3,3);
    MatrixXd zero = MatrixXd::Zero(3,3);
    MatrixXd zero_col = MatrixXd::Zero(3,1);

    this->cache.Fk.setIdentity();
    //                R                          p     v       ba               bg                  plf   prf    ylf       yrf
    this->cache.Fk << G0.transpose(),            zero, zero,   zero,            -G1.transpose()*dt, zero, zero,  zero_col, zero_col,
                      -R*skew(G2*a_k)*pow(dt,2), eye,  eye*dt, -R*G2*pow(dt,2), R*Psi2*pow(dt,3),   zero, zero,  zero_col, zero_col,
                      -R*skew(G1*a_k)*dt,        zero, eye,    -R*G1*dt,        R*Psi1*pow(dt,2),   zero, zero,  zero_col, zero_col,
                      zero,                      zero, zero,   eye,             zero,               zero, zero,  zero_col, zero_col,
                      zero,                      zero, zero,   zero,            eye,                zero, zero,  zero_col, zero_col,
                      zero,                      zero, zero,   zero,            zero,               eye,  zero,  zero_col, zero_col,
                      zero,                      zero, zero,   zero,            zero,               zero, eye,   zero_col, zero_col,
                      0,0,0,                     0,0,0,0,0,0,  0,0,0,           0,0,0,              0,0,0,0,0,0, 1.0,      0,
                      0,0,0,                     0,0,0,0,0,0,  0,0,0,           0,0,0,              0,0,0,0,0,0, 0,        1.0;

    //                na    nw    nba   nbw   nclf  ncrf  nylf      nyrf
    this->cache.Gc << zero, eye,  zero, zero, zero, zero, zero_col, zero_col,
                      zero, zero, zero, zero, zero, zero, zero_col, zero_col,
                      R,    zero, zero, zero, zero, zero, zero_col, zero_col,
                      zero, zero, eye,  zero, zero, zero, zero_col, zero_col,
                      zero, zero, zero, eye,  zero, zero, zero_col, zero_col,
                      zero, zero, zero, zero, R,   zero,  zero_col, zero_col,
                      zero, zero, zero, zero, zero, R,    zero_col, zero_col,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  1,        0,
                      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  0,        1;

    // Update
    packState(this->memory.X, R_pred, p_pred, v_pred, ba, bg, pLF_pred, pRF_pred, footYaws_pred);
    this->memory.Qc = this->cache.Fk * this->cache.Gc * Qcont * this->cache.Gc.transpose() * this->cache.Fk.transpose() * dt;
    this->memory.P = this->cache.Fk * this->memory.P * this->cache.Fk.transpose() + this->memory.Qc;
}

void contact_ekf::update_forward_kinematics(VectorXd &w, VectorXd &encoders, VectorXd &dencoders, VectorXd &contact) {
    // Function to perform Right-Invariant EKF update from forward kinematic measurements
    // Compute all forward kinematics
    Matrix3d RLF_enc, RRF_enc;
    Vector3d pRF_enc, pLF_enc;
    MatrixXd JRF_enc(3,14), JLF_enc(3,14);
    MatrixXd JrotLF_enc(3,14), JrotRF_enc(3,14);
    relative_foot_positions(encoders, pLF_enc, pRF_enc, RLF_enc, RRF_enc);
    relative_foot_jacobians(encoders, JLF_enc, JRF_enc, JrotLF_enc, JrotRF_enc);

    // Predicted state
    MatrixXd zero_col = MatrixXd::Zero(3,1);
    Matrix3d zero = Matrix3d::Zero();
    Matrix3d R;
    Vector2d footYaws;
    Vector3d p, v, wa, wg, ba, bw, plf, prf;
    unpackState(this->memory.X, R, p, v, ba, bw, plf, prf, footYaws);
    Matrix3d Rt = R.transpose();
    Vector3d wmb =  w - bw;

    // Floating base
    Vector3d pLF_pred = Rt * (plf - p);
    Vector3d pRF_pred = Rt * (prf - p);
    Vector3d vLF = skew(wmb)*pLF_enc + JLF_enc * dencoders;
    Vector3d vRF = skew(wmb)*pRF_enc + JRF_enc * dencoders;
    Vector3d v_pred = -R.transpose() * v;

    // Foot Rotations
    Matrix3d Rz_LF, Rz_RF;
    Rz_LF << cos(footYaws(0)), -sin(footYaws(0)), 0,
             sin(footYaws(0)), cos(footYaws(0)),  0,
             0,                  0,               1;
    Rz_RF << cos(footYaws(1)), -sin(footYaws(1)), 0,
             sin(footYaws(1)), cos(footYaws(1)),  0,
             0,                  0,               1;
    Matrix3d RLF_pred, RRF_pred;
    RLF_pred = Rt * Rz_LF;
    RRF_pred = Rt * Rz_RF;

    Vector2d footYaws_pred;
    Eigen::EulerAnglesXYZd euler;
    eulerXYZ(RLF_pred, euler);
    footYaws_pred(0) = euler.gamma();
    eulerXYZ(RRF_pred, euler);
    footYaws_pred(1) = euler.gamma();

    // Encoder noise
    MatrixXd Renc(14,14);
    Renc.setZero();
    Renc.diagonal() << pow(this->config.encoder_noise_std,2) * VectorXd::Ones(14);
    MatrixXd Rcon(3,3);
    Rcon.setZero();
    Rcon.diagonal() << pow(this->config.contact_noise_std,2) * VectorXd::Ones(3);
    MatrixXd Rdcon(3,3);
    Rdcon.setZero();
    Rdcon.diagonal() << pow(this->config.dcontact_noise_std,2) * VectorXd::Ones(3);
    double Rzcon = pow(this->config.contact_yaw_noise_std,2);

    // Update measurement based on contact condition
    VectorXd y, h;
    if ( (contact(0) >= 0.25) && (contact(1) >= 0.25) ) {
        // Double support
        y.resize(6);
        y << pLF_enc, pRF_enc;
        h.resize(6);
        h << pLF_pred, pRF_pred;
        this->cache.Hk.resize(6,23); this->cache.Hk.setZero();
        this->cache.Rk.resize(6,6); this->cache.Rk.setZero();
        this->cache.Hk << skew(pLF_pred), -Rt,  zero, zero, zero, Rt, zero, zero_col, zero_col,
                          skew(pRF_pred), -Rt,  zero, zero, zero, zero, Rt, zero_col, zero_col;
        this->cache.Rk.block(0,0,3,3) << JLF_enc * Renc * JLF_enc.transpose() + Rcon;
        this->cache.Rk.block(3,3,3,3) << JRF_enc * Renc * JRF_enc.transpose() + Rcon;

    } else if (contact(1) >= 0.25) {
        // Right support
        y.resize(3);
        y << pRF_enc;
        h.resize(3);
        h << pRF_pred;
        this->cache.Hk.resize(3,23); this->cache.Hk.setZero();
        this->cache.Rk.resize(3,3); this->cache.Rk.setZero();
        this->cache.Hk << skew(pRF_pred), -Rt, zero, zero, zero, zero, Rt, zero_col, zero_col;
        this->cache.Rk.block(0,0,3,3) << JRF_enc * Renc * JRF_enc.transpose() + Rcon;

    } else if (contact(0) >= 0.25) {
        // Left support
        y.resize(3);
        y << pLF_enc;
        h.resize(3);
        h << pLF_pred;
        this->cache.Hk.resize(3,23); this->cache.Hk.setZero();
        this->cache.Rk.resize(3,3); this->cache.Rk.setZero();
        this->cache.Hk << skew(pLF_pred), -Rt, zero, zero, zero, Rt, zero, zero_col, zero_col;
        this->cache.Rk.block(0,0,3,3) << JLF_enc * Renc * JLF_enc.transpose() + Rcon;
    } else {
        // Don't update
        return;
    }

    // Kalman gain step
    // Standard EKF update
    VectorXd dx(23);
    MatrixXd Sk(this->cache.Rk.rows(), this->cache.Rk.cols());
    MatrixXd K(this->memory.P.rows(), this->cache.Rk.rows());
    Sk = this->cache.Hk * this->memory.P * this->cache.Hk.transpose() + this->cache.Rk;
    K = ( this->memory.P * this->cache.Hk.transpose() ) * Sk.inverse();
    dx = K * ( y - h );
    this->memory.P = (MatrixXd::Identity(23,23) - K * this->cache.Hk) * this->memory.P * (MatrixXd::Identity(23,23) - K * this->cache.Hk).transpose() + K * this->cache.Rk * K.transpose();

    // Update values in X
    Vector3d alpha = dx.block(0,0,3,1);
    R = R * Exp_SO3(alpha);
    p += dx.block(3,0,3,1);
    v += dx.block(6,0,3,1);
    ba += dx.block(9,0,3,1);
    bw += dx.block(12,0,3,1);
    plf += dx.block(15,0,3,1);
    prf += dx.block(18,0,3,1);
    footYaws += dx.block(21,0,2,1);
    packState(this->memory.X, R, p, v, ba, bw, plf, prf, footYaws);
}

void contact_ekf::unpackState(MatrixXd &X, Matrix3d &R, Vector3d &p, Vector3d &v, Vector3d &ba, Vector3d &bw, Vector3d &plf, Vector3d &prf, Vector2d &footYaws) {
    R << X.block(0,0,3,3);
    p << X.col(3);
    v << X.col(4);
    ba << X.col(5);
    bw << X.col(6);
    plf << X.col(7);
    prf << X.col(8);
    footYaws(0) = X(0,9);
    footYaws(1) = X(1,9);
}

void contact_ekf::packState(MatrixXd &X, Matrix3d &R, Vector3d &p, Vector3d &v, Vector3d &ba, Vector3d &bw, Vector3d &plf, Vector3d &prf, Vector2d &footYaws) {
    X.block(0,0,3,3) << R;
    X.col(3) << p;
    X.col(4) << v;
    X.col(5) << ba;
    X.col(6) << bw;
    X.col(7) << plf;
    X.col(8) << prf;
    X.col(9) << footYaws, 0.0;
}

void contact_ekf::relative_foot_positions(VectorXd &enc, Vector3d &plf, Vector3d &prf, Matrix3d &Rlf, Matrix3d &Rrf) {
    VectorXd q(22); q.setZero();
    for (int i=0; i<this->robot->iEncoderMap.size(); i++)
        q(this->robot->iEncoderMap(i)) = enc(i);
    Matrix3d Rx, Ry, Rz;
    MatrixXd temp(6,1);
    SymFunction::pose_leftFoot(temp,q);
    plf << temp.block(0,0,3,1);
    Rx << 1.0, 0, 0,
          0, cos(temp(3)), -sin(temp(3)),
          0, sin(temp(3)), cos(temp(3));
    Ry << cos(temp(4)), 0., sin(temp(4)),
          0.,1.,0.,
          -sin(temp(4)), 0., cos(temp(4));
    Rz << cos(temp(5)), -sin(temp(5)), 0,
          sin(temp(5)), cos(temp(5)), 0,
          0, 0, 1.0;
    Rlf = Rz * Ry * Rx;

    SymFunction::pose_rightFoot(temp,q);
    prf << temp.block(0,0,3,1);

    Rx << 1.0, 0, 0,
          0, cos(temp(3)), -sin(temp(3)),
          0, sin(temp(3)), cos(temp(3));
    Ry << cos(temp(4)), 0., sin(temp(4)),
          0.,1.,0.,
          -sin(temp(4)), 0., cos(temp(4));
    Rz << cos(temp(5)), -sin(temp(5)), 0,
          sin(temp(5)), cos(temp(5)), 0,
          0, 0, 1.0;
    Rrf = Rz * Ry * Rx;
}

void contact_ekf::relative_foot_jacobians(VectorXd &enc, MatrixXd &Jlf, MatrixXd &Jrf, MatrixXd &JrotLF_enc, MatrixXd &JrotRF_enc) {
    VectorXd q(22); q.setZero();
    for (int i=0; i<this->robot->iEncoderMap.size(); i++)
        q(this->robot->iEncoderMap(i)) = enc(i);

    MatrixXd temp(6,22);
    SymFunction::J_leftFoot(temp, q);
    for (int i=0; i<this->robot->iEncoderMap.size(); i++) {
        Jlf.block(0,i,3,1) << temp.block(0,this->robot->iEncoderMap(i),3,1);
        JrotLF_enc.block(0,i,3,1) << temp.block(3,this->robot->iEncoderMap(i),3,1);
    }

    SymFunction::J_rightFoot(temp, q);
    for (int i=0; i<this->robot->iEncoderMap.size(); i++) {
        Jrf.block(0,i,3,1) << temp.block(0,this->robot->iEncoderMap(i),3,1);
        JrotRF_enc.block(0,i,3,1) << temp.block(3,this->robot->iEncoderMap(i),3,1);
    }
}

int contact_ekf::factorial(int n)
{
    // single line to find factorial
    return (n==1 || n==0) ? 1: n * factorial(n - 1);
}

Matrix3d contact_ekf::skew(Vector3d v) {
    Matrix3d A;
    A <<  0.,   -v(2),   v(1),
          v(2),  0.   , -v(0),
         -v(1),  v(0),   0.;
    return A;
}

Matrix3d contact_ekf::Exp_SO3(Vector3d w) {
    // Computes the vectorized exponential map for SO(3)
    Eigen::Matrix3d A = skew(w);
    double theta = w.norm();
    if (theta < 1e-10) {
        return Eigen::Matrix3d::Identity();
    }
    Eigen::Matrix3d R =  Eigen::Matrix3d::Identity() + (sin(theta)/theta)*A + ((1-cos(theta))/(theta*theta))*A*A;
    return R;
}

Matrix3d contact_ekf::Gamma(Vector3d w, int n) {
    Matrix3d output = Matrix3d::Identity();
    Matrix3d R = Exp_SO3(w);
    Matrix3d A = skew(w);
    double theta = sqrt(pow(w(0),2) + pow(w(1),2) + pow(w(2),2));

    if (theta <= 0.00000000001) {
        output = (1.0 / factorial(n)) * Matrix3d::Identity();
        return output;
    }

    Matrix3d S = Matrix3d::Identity();
    for (int k = 1; k<n+1; k++) {
        // Calculate matrix power
        Matrix3d Apow = A;
        for (int j=1; j<k; j++) {
            Apow = Apow * Apow;
        }
        // Sum
        S += Apow / factorial(k);
    }

    if (n==0)
        output = R;
    else if (n % 2) {
        double powr = pow(-1,(n+1)/2) / pow(theta, n+1.0);
        output = (1.0/factorial(n)) * Matrix3d::Identity() + powr * A * (R - S);
    } else {
        double powr = pow(-1,(n)/2) / pow(theta, n);
        output = (1.0/factorial(n)) * Matrix3d::Identity() + powr * (R - S);
    }
    return output;
}
