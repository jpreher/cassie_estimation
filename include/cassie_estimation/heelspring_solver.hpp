/*
 * @author Jenna Reher (jreher@caltech.edu)
 */

#ifndef HEELSPRING_SOLVER_HPP
#define HEELSPRING_SOLVER_HPP

#include <cassie_description/cassie_model.hpp>
#include <ros_utilities/ros_utilities.hpp>
#include <std_srvs/Empty.h>

class HeelspringSolver {

public:
    cassie_model::Cassie * robot;

    HeelspringSolver(ros::NodeHandle &nh, cassie_model::Cassie &robot);
    void update();
    void reconfigure();

private:
    struct Cache {
        VectorXd q;
        void init() {this->q.resize(22); this->reset();}
        void reset() {this->q.setZero();}
    } cache;

    struct Memory {
        Vector2d prev_x;
        Vector4d spring_offset;
        int nZeroSamplesAccrued;

        void init() {this->reset();}
        void reset() {this->prev_x.setZero(); this->nZeroSamplesAccrued=0; this->spring_offset.setZero();}
    } memory;

    struct Config {
        bool zero_springs_on_startup;
        int num_zeroing_samples;

        void init() {this->zero_springs_on_startup = false; this->num_zeroing_samples=100;}
        void reconfigure();

        ros_utilities::ParamChecker paramChecker;
    } config;

    ros::ServiceServer reconfigureService;
    bool reconfigure(std_srvs::Empty::Request &req, std_srvs::Empty::Response &res);
    ros::NodeHandle *nh;

    void IKfunction(Eigen::VectorXd &x, Eigen::VectorXd &fvec);
    void J_IKfunction(Eigen::VectorXd &x, Eigen::MatrixXd &J);
    int InverseKinematics(Eigen::VectorXd &x);

};


#endif // HEELSPRING_SOLVER_HPP
