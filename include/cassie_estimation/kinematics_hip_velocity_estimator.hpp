/*
 * @author Jenna Reher (jreher@caltech.edu)
 */

#ifndef KINEMATICS_HIP_VELOCITY_ESTIMATOR_HPP
#define KINEMATICS_HIP_VELOCITY_ESTIMATOR_HPP

#include <cassie_description/cassie_model.hpp>
#include <control_utilities/filters.hpp>
#include <cassie_common_toolbox/smoothing.hpp>
#include <ros_utilities/ros_utilities.hpp>
#include <std_srvs/Empty.h>

class KinematicsHipVelocityEstimator {

public:
    KinematicsHipVelocityEstimator(ros::NodeHandle &nh, cassie_model::Cassie &robot, bool do_update_robot);

    void initialize();
    void reset();
    void update();
    Vector3d getValue();
    void reconfigure();

private:
    // Robot model
    bool do_update_robot;
    cassie_model::Cassie *robot;
    Vector3d velocity;

    struct Config {
        double control_rate;
        double lowpass_vx_dt_cutoff;
        double lowpass_vy_dt_cutoff;
        double lowpass_vz_dt_cutoff;

        void init() {this->reconfigure();}
        void reconfigure();
        ros_utilities::ParamChecker paramChecker;
    } config;

    ros::NodeHandle *nh;
    ros::ServiceServer reconfigureService;
    bool reconfigure(std_srvs::Empty::Request &req, std_srvs::Empty::Response &res);

    control_utilities::LowPassFilter lpVX = control_utilities::LowPassFilter(0.0005, 0.12);
    control_utilities::LowPassFilter lpVY = control_utilities::LowPassFilter(0.0005, 0.12);
    control_utilities::LowPassFilter lpVZ = control_utilities::LowPassFilter(0.0005, 0.12);
};



#endif // KINEMATICS_HIP_VELOCITY_ESTIMATOR_HPP
