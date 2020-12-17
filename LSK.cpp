#include <stdio.h>
#include <math.h>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/SVD>

using Eigen::MatrixXd;
using namespace Eigen;


MatrixXd beta(2, 9);
MatrixXd offset(1, 9); //ワイヤ配置オフセット
MatrixXd A(3, 2); // 変換行列
MatrixXd rc(1, 9); //ワイヤピッチ円半径
MatrixXd L(1, 9); // sectionあたり長さ
MatrixXd wire_length(3, 9); // ワイヤ引き量
double N; // 1セクション内の1自由度あたりのセグメント数
MatrixXd theta(1, 9);　//曲げ角θ
MatrixXd phi(1, 9);　//曲げ角φ


void dataset(){
    // the dataset of robot
    ::offset << 20, 60, 100, 0, 40, 80, 90, 130, 170;
    ::offset = offset * M_PI / 180;
    ::rc << 4.6, 4.6, 4.6, 4.8, 4.8, 4.8, 4.6, 4.6, 4.6;
    ::N = 5;
    ::L = 120,120,120,100,100,100,75,75,75;

}

MatrixXd pseudo_inverse(MatrixXd A){
    // SVDによる疑似逆行列計算
    JacobiSVD<MatrixXd> SVD(A, ComputeThinU | ComputeThinV);
    MatrixXd Ainv(3, 2);
    Ainv = SVD.matrixV() * SVD.singularValues().asDiagonal().inverse() * SVD.matrixU().transpose();

    return Ainv;
}

MatrixXd SVD_solve(MatrixXd A, MatrixXd w){
    //SVDを用いた最小二乗法によるAx＝ｗのxについての解
    JacobiSVD<MatrixXd> SVD(A, ComputeThinU | ComputeThinV);
    MatrixXd b(2, 1);
    b = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(w);

    return b;
}

MatrixXd inv_LSK(MatrixXd theta,MatrixXd phi){
    // inverse LSK 

    for (int i = 0; i < 9;i++){
        beta(0, i) = theta(i)/N * M_PI /180;
        beta(1, i) = phi(i) / N * M_PI / 180;
    }
    
    for (int i = 0; i < 9; i++)
    {

        A << std::cos(offset(i)), std::sin(offset(i)),
            std::cos(offset(i) + M_PI * 2 / 3), std::sin(offset(i) + M_PI * 2 / 3),
            std::cos(offset(i) + M_PI * 4 / 3), std::sin(offset(i) + M_PI * 4 / 3);

        wire_length.col(i) = -N * rc(i) * A * beta.block(0, 0, 2, i + 1).rowwise().sum();

    }
    return wire_length;
}

MatrixXd fw_LSK(MatrixXd wire_length){
    // forward LSK
    for (int i = 0; i < 9; i++)
    {
        A << std::cos(offset(i)), std::sin(offset(i)),
            std::cos(offset(i) + M_PI * 2 / 3), std::sin(offset(i) + M_PI * 2 / 3),
            std::cos(offset(i) + M_PI * 4 / 3), std::sin(offset(i) + M_PI * 4 / 3);


        A = -N * rc(i) * A;

        if (i == 0){
            beta.col(i) = SVD_solve(A, wire_length.col(i));

        }
        else{
            beta.col(i) = SVD_solve(A, wire_length.col(i)) - beta.block(0, 0, 2, i ).rowwise().sum();
        }
    }
    beta = beta * 180 / M_PI*N;
    return beta;
}

MatrixXd pose_CR(MatrixXd beta){
    // the pose in the cartesian space of the robot from the configuration angles
}

int main(){
    dataset();
    ::theta << 90, -90 , 90, 0, 0, 0, 0, 0, 90;
    ::phi << 0, 0, 0, 0, 0, 0, 0, 0, 0;

    wire_length = inv_LSK(theta, phi);
    std::cout << wire_length << std::endl;
    beta =  fw_LSK(wire_length);
    std::cout << beta << std::endl;
}