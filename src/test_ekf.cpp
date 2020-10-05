// ROS
#include <ros/ros.h>

// Main ekf
#include <cassie_description/cassie_model.hpp>
#include <cassie_estimation/contact_ekf.hpp>
#include <cassie_estimation/contact_classifier.hpp>


// Logging
#include <string>
#include <fstream>
#include <iostream>


// Main node
int main(int argc, char *argv[])
{
    // Establish the current ROS node and associated timing
    ros::init(argc, argv, "test_node");
    ros::NodeHandle nh("/cassie/interface");

    // Output file for plotting
    VectorXd log = VectorXd::Zero(39);
    std::fstream logfile;
    std::string home=getenv("HOME");
    std::string path= home+"/datalog/ekf/ekf_test.bin";
    logfile.open(path, std::ios::out | std::ios::binary);

    // Input file from hardware
    std::string path_in= home+"/datalog/ekf/estimation_log.bin";
    std::ifstream myFile (path_in, std::ios::in | std::ios::binary);
    if (!myFile) {
        std::cout << "Couldn't load input file!!" << std::endl;
        logfile.close();
        return 0;
    }

    myFile.seekg(0, std::ios::end);
    std::streampos size = myFile.tellg();

    double *fileBuffer;
    unsigned long long sizeOfBuffer = size / sizeof(double);
    fileBuffer = new double[sizeOfBuffer];

    std::cout << "Open file is of size: " << size << std::endl;
    std::cout << "With buffer size:     " << sizeOfBuffer << std::endl;

    myFile.seekg(0, std::ios::beg);
    myFile.read(reinterpret_cast<char*>(fileBuffer), size);

    unsigned long long ndbge = 46;
    unsigned long long nlogse = sizeOfBuffer / ndbge;
    std::cout << "Total log entries:    " << nlogse << std::endl;
    VectorXd te      = VectorXd::Zero(nlogse);
    MatrixXd quat    = MatrixXd::Zero(4,nlogse);
    MatrixXd gyro    = MatrixXd::Zero(3,nlogse);
    MatrixXd accel   = MatrixXd::Zero(3,nlogse);
    MatrixXd v       = MatrixXd::Zero(3,nlogse);
    MatrixXd enc     = MatrixXd::Zero(14,nlogse);
    MatrixXd denc    = MatrixXd::Zero(14,nlogse);
    MatrixXd ach     = MatrixXd::Zero(2,nlogse);
    MatrixXd contact = MatrixXd::Zero(2,nlogse);
    for (unsigned long long i=0; i<nlogse; i++) {
        unsigned long long j = ndbge * i;
        te(i) = fileBuffer[j];
        quat.col(i)    << fileBuffer[j+1], fileBuffer[j+2], fileBuffer[j+3], fileBuffer[j+4];
        gyro.col(i)    << fileBuffer[j+5], fileBuffer[j+6], fileBuffer[j+7];
        accel.col(i)   << fileBuffer[j+8], fileBuffer[j+9], fileBuffer[j+10];
        v.col(i)       << fileBuffer[j+11], fileBuffer[j+12], fileBuffer[j+13];
        enc.col(i)     << fileBuffer[j+14], fileBuffer[j+15], fileBuffer[j+16], fileBuffer[j+17], fileBuffer[j+18], fileBuffer[j+19], fileBuffer[j+20],
                          fileBuffer[j+21], fileBuffer[j+22], fileBuffer[j+23], fileBuffer[j+24], fileBuffer[j+25], fileBuffer[j+26], fileBuffer[j+27];
        denc.col(i)    << fileBuffer[j+28], fileBuffer[j+29], fileBuffer[j+30], fileBuffer[j+31], fileBuffer[j+32], fileBuffer[j+33], fileBuffer[j+34],
                          fileBuffer[j+35], fileBuffer[j+36], fileBuffer[j+37], fileBuffer[j+38], fileBuffer[j+39], fileBuffer[j+40], fileBuffer[j+41];
        ach.col(i)     << fileBuffer[j+42], fileBuffer[j+43];
        contact.col(i) << fileBuffer[j+44], fileBuffer[j+45];
    }
    te = te - te(0) * VectorXd::Ones(te.size());

    // Build estimator
    cassie_model::Cassie robot;
    contact_ekf ekf(nh, robot, true);
    ContactClassifier contact_classifier(nh, robot, 0.0005);

    // Run and log
    for (unsigned long long i=0; i<nlogse; i++) {
        // Get values at timestep
        VectorXd w, a, encoder, dencoder, con, quatern, achilles;
        w        = gyro.col(i);
        a        = accel.col(i);
        encoder  = enc.col(i);
        dencoder = denc.col(i);
        con      = contact.col(i);
        quatern  = quat.col(i);
        achilles = ach.col(i);

        double dt = 0;
        if (i>0)
            dt = te(i) - te(i-1);

        // Update the robot model
        robot.q.setZero();
        robot.dq.setZero();
        for (int i=0; i<robot.iEncoderMap.size(); i++) {
            robot.q(robot.iEncoderMap(i)) = encoder(i);
            robot.dq(robot.iEncoderMap(i)) = dencoder(i);
        }
        robot.q(LeftHeelSpring) = achilles(0);
        robot.q(RightHeelSpring) = achilles(1);

        // Do euler angles - SUBTRACTING OUT THE YAW!!!
        Eigen::Quaterniond quat(quatern(0), quatern(1), quatern(2), quatern(3));
        Eigen::Matrix3d Ro = quat.toRotationMatrix();
        Eigen::EulerAnglesXYZd euler = Eigen::EulerAnglesXYZd::FromRotation<false, false, false>(quat);
        eulerXYZ(quat, euler);
        Eigen::Matrix3d Rz;
        Rz << cos(euler.gamma()), -sin(euler.gamma()), 0,
              sin(euler.gamma()), cos(euler.gamma()),  0,
              0,                  0,                   1;
        Ro = Rz.transpose() * Ro;
        Eigen::Quaterniond tempquat(Ro);
        eulerXYZ(tempquat, euler);
        robot.q(BaseRotX) = euler.alpha(); // roll
        robot.q(BaseRotY) = euler.beta();  // pitch
        robot.q(BaseRotZ) = euler.gamma(); // yaw

        // Do contact estimation
        contact_classifier.update();
        if (isnan((robot.leftContact))) {
            robot.leftContact = 0.;
            ROS_WARN("Left contact nan!");
        }
        if (isnan((robot.rightContact))) {
            robot.rightContact = 0.;
            ROS_WARN("Right contact nan!");
        }
        con << robot.leftContact, robot.rightContact; // override stored measurement

        // Update
        dt = 0.0005;
        ekf.update(dt, w, a, encoder, dencoder, con);

        // Extract result
        Matrix3d R;
        Vector2d footYaws;
        Vector3d pos, vel, ba, bg, plf, prf;
        ekf.getValues(R,pos,vel,ba,bg,plf,prf,footYaws);

        // Rotate yaw
        eulerXYZ(R, euler);
        Rz << cos(euler.gamma()), -sin(euler.gamma()), 0,
              sin(euler.gamma()), cos(euler.gamma()), 0,
              0, 0, 1.0;
        vel = Rz.transpose() * vel;

        // Rotate into foot frames
        VectorXd q(22); q.setZero();
        for (int i=0; i<robot.iEncoderMap.size(); i++)
            q(robot.iEncoderMap(i)) = encoder(i);
        Matrix3d Rzl, Rzr;
        MatrixXd temp(6,1);
        SymFunction::pose_leftFoot(temp,q);
        Rzl << cos(temp(5)), -sin(temp(5)), 0,
               sin(temp(5)), cos(temp(5)), 0,
               0, 0, 1.0;

        SymFunction::pose_rightFoot(temp,q);
        Rzr << cos(temp(5)), -sin(temp(5)), 0,
               sin(temp(5)), cos(temp(5)), 0,
               0, 0, 1.0;

        // Rotate into frame relative to stance foot
        if ( (con(0) >= 0.25) && (con(1) >= 0.25) ) {
            // Double support
            VectorXd vl(3), vr(3);
            vl = Rzl * vel;
            vr = Rzr * vel;
            vel = (vl + vr) / 2.0;
        } else if (con(1) >= 0.25) {
            // Right support
            vel = Rzr * vel;
        } else if (con(0) >= 0.25) {
            // Left support
            vel = Rzl * vel;
        }

        log << te(i),    // 1
               R(0,0), R(0,1), R(0,2), // 3
               R(1,0), R(1,1), R(1,2), // 3
               R(2,0), R(2,1), R(2,2), // 3
               pos,      // 3
               vel,      // 3
               ekf.getRawVelocity(), //v.col(i), // 3
               ba,       // 3
               bg,       // 3
               plf,      // 3
               prf,      // 3
               con,      // 2
               contact_classifier.grf;
        logfile.write(reinterpret_cast<char *>(log.data()), (log.size())*sizeof(double));
    }

    // Close down
    free(fileBuffer);
    logfile.close();
    myFile.close();
    return 0;
}
