/*
 * @author Jenna Reher (jreher@caltech.edu)
 */

#ifndef CONTACT_CLASSIFIER_HPP
#define CONTACT_CLASSIFIER_HPP

#include <cassie_description/cassie_model.hpp>
#include <control_utilities/filters.hpp>
#include <ros/ros.h>
#include <ros_utilities/ros_utilities.hpp>
#include <std_srvs/Empty.h>

class ContactClassifier {

public:

    VectorXd grf;

    ContactClassifier(ros::NodeHandle &nh, cassie_model::Cassie &robot, double dt);
    void update();
    void reconfigure();

private:
    control_utilities::LowPassFilter LowPassLeft;
    control_utilities::LowPassFilter LowPassRight;
    cassie_model::Cassie *robot;

    struct Config {
        double dt;
        bool use_sigmoid;
        double sigmoid_A;
        double sigmoid_B;
        double sigmoid_power;
        double lowpass_dt_cutoff;
        double linear_lb;
        double linear_ub;

        ros_utilities::ParamChecker paramChecker;

        void init();
        void reconfigure();
    } config;

    ros::ServiceServer reconfigureService;
    bool reconfigure(std_srvs::Empty::Request &req, std_srvs::Empty::Response &res);
    ros::NodeHandle *nh;
};

#endif // CONTACT_CLASSIFIER_HPP
