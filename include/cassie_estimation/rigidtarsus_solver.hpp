/*
 * @author Jenna Reher (jreher@caltech.edu)
 */

#ifndef RIGIDTARSUS_SOLVER_HPP
#define RIGIDTARSUS_SOLVER_HPP

#include <cassie_description/cassie_model.hpp>
#include <ros_utilities/ros_utilities.hpp>
#include <std_srvs/Empty.h>

class RigidTarsusSolver {

public:
    cassie_model::Cassie * robot;

    RigidTarsusSolver(cassie_model::Cassie &robot);
    void update();
    void reconfigure();
    void reset();
    double getLeftRigidTarsusPosition() {return this->memory.tar_sol[0];}
    double getLeftRigidTarsusVelocity() {return this->memory.dtar_sol[0];}
    double getRightRigidTarsusPosition() {return this->memory.tar_sol[1];}
    double getRightRigidTarsusVelocity() {return this->memory.dtar_sol[1];}

private:
    struct Cache {
        VectorXd q;
        void init() {this->q.resize(22); this->reset();}
        void reset() {this->q.setZero();}
    } cache;

    struct Memory {
        Vector2d prev_x;
        bool is_initialized;
        Vector2d tar_sol;
        Vector2d dtar_sol;

        void init() {this->reset();}
        void reset() {this->prev_x.setZero(); tar_sol.setZero(); dtar_sol.setZero(); this->is_initialized=false;}
    } memory;

    void IKfunction(Eigen::VectorXd &x, Eigen::VectorXd &fvec);
    void J_IKfunction(Eigen::VectorXd &x, Eigen::MatrixXd &J);
    int InverseKinematics(Eigen::VectorXd &x);

};


#endif // RIGIDTARSUS_SOLVER_HPP
