/*
 * @author Jenna Reher (jreher@caltech.edu)
 *
 * This filter is adapted from several sources, including
 *  - Bloesch, Michael, et al. "State estimation for legged robots-consistent fusion of leg kinematics and IMU." Robotics 17 (2013): 17-24.
 *  - Reher, J, Wen-Loong Ma, and Aaron D. Ames. "Dynamic walking with compliance on a cassie bipedal robot." 2019 18th European Control Conference (ECC). IEEE, 2019.
 *  - Several functions released in an example implementation of Bloesch's method from Ross Hartley: https://github.com/UMich-BipedLab/Cassie_StateEstimation/blob/master/Estimators/Quaternion_EKF/QuaternionEKF.m
 */

#ifndef CONTACT_EKF_HPP
#define CONTACT_EKF_HPP

#include <cassie_description/cassie_model.hpp>
#include <control_utilities/filters.hpp>
#include <ros_utilities/ros_utilities.hpp>
#include <ros_utilities/timing.hpp>
#include <std_srvs/Empty.h>

using namespace Eigen;

class contact_ekf {

public:
    contact_ekf(ros::NodeHandle &nh, cassie_model::Cassie &robot, bool do_update_robot);

    void initialize();
    void reset();
    void reset(Vector3d &initial_velocity);
    void update(double dt, VectorXd &w, VectorXd &a, VectorXd &encoders, VectorXd &dencoders, VectorXd &contact);
    void getValues(Matrix3d &R, Vector3d &p, Vector3d &v, Vector3d &ba, Vector3d &bg, Vector3d &plf, Vector3d &prf, Vector2d &footYaws);
    Vector3d getRawVelocity();
    bool isEnabled();
private:
    /**
    * @brief Computations that can be reproduced given a Config and Memory object
    * Meant to only be valid during one computation frame. Anything meant to persist across time should be stored in Memory.
    */
    struct Cache {
        MatrixXd Fc;
        MatrixXd Fk;
        MatrixXd Gc;
        MatrixXd Hk;
        MatrixXd Rk;
        MatrixXd Sk;
        MatrixXd K;
        MatrixXd Qscaled;

        void init();
        void reset();
    } cache;

    /**
    * @brief Data meant to persist across data frames
    */
    struct Memory {
        ros_utilities::Timer timer = ros_utilities::Timer(true);
        double t_last;
        bool filter_enabled;
        bool bias_initialized;
        int contact_mode;

        MatrixXd X;
        MatrixXd P;
        MatrixXd Qc;

        Vector3d ba0;
        Vector3d bg0;
        Vector3d w_prev;
        Vector3d a_prev;
        VectorXd encoders_prev;
        VectorXd dencoders_prev;
        Vector2d contact_prev;
        Vector3d v_init;

        void init();
        void reset();
    } memory;

    /**
    * @brief Persistent configurations
    * Meant to be configurable parameters that are pulled from YAML or rosparam
    */
    struct Config {
        double dt;
        Vector3d g;
        bool do_update_robot;

        // Filter parameters
        double gyro_noise_std;
        double gyro_bias_noise_std;
        double accel_noise_std;
        double accel_bias_noise_std;
        double contact_noise_std;
        double contact_yaw_noise_std;
        double dcontact_noise_std;
        double encoder_noise_std;
        double contact_yaw_std;

        double position_std;
        double velocity_std;
        double contact_std;
        
        double prior_base_rotation_std;
        double prior_base_position_std;
        double prior_base_velocity_std;
        double prior_contact_position_std;
        double prior_gyro_bias_std;
        double prior_accel_bias_std;
        double prior_forward_kinematics_std;
        double prior_contact_yaw_std;

        bool apply_post_filter;
        double post_filter_dt_cutoff_x;
        double post_filter_dt_cutoff_y;
        double post_filter_dt_cutoff_z;

        // Sensor Covariances
        MatrixXd Qcontinuous;
        MatrixXd P_prior;

        // Parameter checker
        ros_utilities::ParamChecker paramChecker;

        // Methods
        void init();
        void reconfigure();

    } config;

    // Pointer to the controlling nodehandle and related ROS things
    ros::NodeHandle *nh;
    ros::ServiceServer reconfigureService;
    bool reconfigure();
    bool reconfigure(std_srvs::Empty::Request &req, std_srvs::Empty::Response &res) ;

    // Lowpass post-filter for velocity
    control_utilities::LowPassFilter lpVX = control_utilities::LowPassFilter(0.0005, 0.12);
    control_utilities::LowPassFilter lpVY = control_utilities::LowPassFilter(0.0005, 0.12);
    control_utilities::LowPassFilter lpVZ = control_utilities::LowPassFilter(0.0005, 0.12);

    // Methods
    void initializeBias(VectorXd &w, VectorXd &a);
    void initializeFilter(VectorXd &encoders, VectorXd &contact);
    void predict_state(double dt);
    void update_forward_kinematics(VectorXd &w, VectorXd &encoders, VectorXd &dencoders, VectorXd &contact);
    
    void packState(MatrixXd &X, Matrix3d &R, Vector3d &p, Vector3d &v, Vector3d &ba, Vector3d &bw, Vector3d &plf, Vector3d &prf, Vector2d &footYaws);
    void unpackState(MatrixXd &X, Matrix3d &R, Vector3d &p, Vector3d &v, Vector3d &ba, Vector3d &bw, Vector3d &plf, Vector3d &prf, Vector2d &footYaws);

    // Robot model
    bool do_update_robot;
    cassie_model::Cassie *robot;
    void relative_foot_positions(VectorXd &enc, Vector3d &plf, Vector3d &prf, Matrix3d &Rlf, Matrix3d &Rrf);
    void relative_foot_jacobians(VectorXd &enc, MatrixXd &Jlf, MatrixXd &Jrf, MatrixXd &JrotLF_enc, MatrixXd &JrotRF_enc);

    // Helper functions
    int factorial(int n);
    Matrix3d skew(Vector3d w);
    Matrix3d Exp_SO3(Vector3d w);
    Matrix3d Gamma(Vector3d w, int n);
};

#endif // CONTACT_EKF_HPP
