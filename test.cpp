#include <stdio.h>
#include <math.h>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/SVD>

#include"CRkinematics.h"

using namespace Eigen;


int main(){
    MatrixXd out;
    CRkinematics test;
    test.dataset(9);
    test.theta << 90, -90 , 90, 0, 0, 0, 0, 0, 0;
    test.phi << 0, 0, 0, 0, 0, 0, 0, 0, 0;

    std::cout << test.theta << std::endl;
 


    std::cout << "--------------" << std::endl;
    test.inv_LSK(test.theta,test.phi);
    std::cout << test.wire_length << std::endl;
 


    std::cout << "--------------" << std::endl;
    test.fw_LSK(test.wire_length);  
    std::cout << test.beta<< std::endl;

    std::cout << "--------------" << std::endl;

    out = test.pose_CR(test.theta,test.phi); 

    std::cout << out.block(0,0,3,4) << std::endl;


    std::cout << "--------------" << std::endl;

    out = test.pose_CR2d(test.theta); 

    std::cout << out.block(0,0,3,4) << std::endl;
}