#ifndef CRkinematics_h
#define CRkinematics_h

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/SVD>


//class
class CRkinematics{
    public:
        void dataset(int section_num);
        Eigen::MatrixXd pseudo_inverse(Eigen::MatrixXd A);
        Eigen::MatrixXd SVD_solve(Eigen::MatrixXd A, Eigen::MatrixXd w);
        Eigen::MatrixXd inv_LSK(Eigen::MatrixXd theta,Eigen::MatrixXd phi);
        Eigen::MatrixXd fw_LSK(Eigen::MatrixXd wire_length);
        Eigen::MatrixXd pose_CR(Eigen::MatrixXd theta,Eigen::MatrixXd phi);
        Eigen::MatrixXd pose_CR2d(Eigen::MatrixXd theta);


        Eigen::MatrixXd beta; //2*9

        Eigen::MatrixXd A; // 変換行列 2*3

        Eigen::MatrixXd wire_length; // ワイヤ引き量 3*9

        Eigen::MatrixXd theta; //曲げ角θ 1*9
        Eigen::MatrixXd phi; //曲げ角φ 1*9
        int section;

    private:
        Eigen::MatrixXd rc; //ワイヤピッチ円半径 1*9
        Eigen::MatrixXd L; // sectionあたり長さ 1*9
        double N; // 1セクション内の1自由度あたりのセグメント数 =5
        Eigen::MatrixXd offset; //ワイヤ配置オフセット 1*9

        Eigen::Matrix4d Rx(double a); //rx
        Eigen::Matrix4d Ry(double a); //ry
        Eigen::Matrix4d Rz(double a); //rz
        Eigen::Matrix4d Tx(double h); //tx
        Eigen::Matrix4d Ty(double h); //tx
        Eigen::Matrix4d Tz(double h); //tx
};
#endif


