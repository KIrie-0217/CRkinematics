#include <stdio.h>
#include <math.h>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/SVD>

#include"CRkinematics.h"

using Eigen::MatrixXd;
using namespace Eigen;


void CRkinematics::dataset(int section_num){

    section = section_num;
    beta.conservativeResize(2, 9);
    A.conservativeResize(3,2);
    wire_length.conservativeResize(3,9);
    theta.conservativeResize(1,9);
    phi.conservativeResize(1,9);
    rc.conservativeResize(1,9);
    L.conservativeResize(1,9);
    offset.conservativeResize(1,9);
  
    // the dataset of robot
    offset << 20, 60, 100, 0, 40, 80, 90, 130, 170;
    offset = offset * M_PI / 180;
    rc << 4.6, 4.6, 4.6, 4.8, 4.8, 4.8, 4.6, 4.6, 4.6;
    N = 5;
    L << 120,120,120,100,100,100,75,75,75;

}

MatrixXd CRkinematics::pseudo_inverse(MatrixXd A){
    // SVDによる疑似逆行列計算
    JacobiSVD<MatrixXd> SVD(A, ComputeThinU | ComputeThinV);
    MatrixXd Ainv(3, 2);
    Ainv = SVD.matrixV() * SVD.singularValues().asDiagonal().inverse() * SVD.matrixU().transpose();

    return Ainv;
}

MatrixXd CRkinematics::SVD_solve(MatrixXd A, MatrixXd w){
    //SVDを用いた最小二乗法によるAx＝ｗのxについての解
    JacobiSVD<MatrixXd> SVD(A, ComputeThinU | ComputeThinV);
    MatrixXd b(2, 1);
    b = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(w);

    return b;
}

MatrixXd CRkinematics::inv_LSK(MatrixXd theta,MatrixXd phi){
    // inverse LSK 

    for (int i = 0; i < section;i++){
        beta(0, i) = theta(0,i)/N * M_PI /180;
        beta(1, i) = phi(0,i) / N * M_PI / 180;
    }
    for (int i = 0; i < section; i++)
    {
        
        A << std::cos(offset(0,i)), std::sin(offset(0,i)),
            std::cos(offset(0,i) + M_PI * 2 / 3), std::sin(offset(0,i) + M_PI * 2 / 3),
            std::cos(offset(0,i) + M_PI * 4 / 3), std::sin(offset(0,i) + M_PI * 4 / 3);

        
        
        wire_length.col(i) = -N * rc(0,i) * A * beta.block(0, 0, 2, i + 1).rowwise().sum();

    }

    return wire_length;
}

MatrixXd CRkinematics::fw_LSK(MatrixXd wire_length){
    // forward LSK
    for (int i = 0; i < section; i++)
    {
        A << std::cos(offset(0,i)), std::sin(offset(i)),
            std::cos(offset(0,i) + M_PI * 2 / 3), std::sin(offset(0,i) + M_PI * 2 / 3),
            std::cos(offset(0,i) + M_PI * 4 / 3), std::sin(offset(0,i) + M_PI * 4 / 3);


        A = -N * rc(0,i) * A;

        if (i == 0){
            beta.col(i) = SVD_solve(A, wire_length.col(i));

        }
        else{
            beta.col(i) = SVD_solve(A, wire_length.col(i)) - beta.block(0, 0, 2, i ).rowwise().sum();
        }
    }
    beta = beta * 180 / M_PI*N;
    theta = beta.col(0);
    phi = beta.col(1);
    return beta,theta,phi;
}

MatrixXd CRkinematics::pose_CR(MatrixXd theta,MatrixXd phi){
    //pose
    MatrixXd xyz(4,1);
    MatrixXd L_;
    L_ = L/10;
    MatrixXd simat;
    MatrixXd simat_;


    for(i=0;i<section;i++){
        simat = Rx(phi(0,i))*Tz(L_(0,i))*Ry(theta(0,1))*Tz(L_(0,1));
        if(i==0){
            simat_ = simat^5;
            pose.conservativeResize(4,3*(i+1))
            pose.block(0,0,4,1) = xyz;
            pose.block(0,3*(i+1)-2,4,1) =  simat * xyz;           
            pose.block(0,3*(i+1)-1,4,1) =  simat^3 * xyz; 
            pose.block(0,3*(i+1),4,1) =  simat^5 * xyz; 
        }
      else{
        pose.conservativeResize(4,3*(i+1))
        pose.block(0,3*(i+1)-2,4,1) =  simat_ *simat xyz;           
        pose.block(0,3*(i+1)-1,4,1) =  simat_ *simat^3 xyz; 
        pose.block(0,3*(i+1),4,1) =  simat_ *simat^5 xyz;
        simat_ = simat_ * simat^5; 

        }
    }

    std::cout << pose << std::endl;
    return pose;

} 


MatrixXd CRkinematics::pose_CR2d(MatrixXd theta){
    MatrixXd a;
    return a;
} 


Matrix4d CRkinematics::Rx(double a){
    Matrix4d R;

    R << 1,0,0,0,
        0,std::cos(a),-1*std::sin(a),0,
        0,std::sin(a),std::cos(a),0,
        0,0,0,1;

    return R;

}


Matrix4d CRkinematics::Ry(double a){
    Matrix4d R;

    R << std::cos(a),0,std::cos(a),0,
        0,1,0,0,
        -1*std::sin(a),0,std::cos(a),0,
        0,0,0,1;

    return R;
}

Matrix4d CRkinematics::Rz(double a){
    Matrix4d R;

    R << std::cos(a),-1*std::sin(a),0,0,
        std::sin(a),std::cos(a),0,0,
        0,0,1,0,
        0,0,0,1;

    return R;
}


Matrix4d CRkinematics::Tx(double h){
    Matrix4d T;

    T << 1,0,0,h,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1;

    return T;
}

Matrix4d CRkinematics::Ty(double h){
    Matrix4d T;

    T << 1,0,0,0,
        0,1,0,h,
        0,0,1,0,
        0,0,0,1;

    return T;
}

Matrix4d CRkinematics::Tz(double h){
    Matrix4d T;

    T << 1,0,0,0,
        0,1,0,0,
        0,0,1,h,
        0,0,0,1;

    return T;    
}