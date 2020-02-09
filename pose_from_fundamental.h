#include "math\matrix.h"
#include "math\matrix_svd.h"
#include <vector>

typedef math::Matrix<double,3,3> FundamentalMatrix;
typedef math::Matrix<double,3,3> EssentialMatrix;

typedef std::pair<math::Matrix3d,math::Vec3d> pose;

/*
* description 对匹配点进行三角化得到空间三维点
* @
*/
math::Vec3d triangulation(const math::Vec2d& p1,
                            const math::Vec2d& p2,
                            const math::Matrix3d& K1,
                            const math::Matrix3d& K2,
                            const math::Vec3d& t1,
                            const math::Vec3d& t2,
                            const math::Matrix3d& R1,
                            const math::Matrix3d& R2)
{
    //构造投影矩阵
    math::Matrix<double,3,4>P1,P2;

    math::Matrix<double,3,3> KR1 = K1 * R1;
    math::Matrix<double,3,1> Kt1(*(K1 * t1));
    P1 = KR1.hstack(Kt1);

    math::Matrix<double,3,3> KR2 = K2 * R2;
    math::Matrix<double,3,1> Kt2(*(K2 * t2));
    P2 = KR2.hstack(Kt2);

    std::cout<<"P1: "<<P1<<std::endl;
    std::cout<<"P1 for fist pose should be\n"
             <<"0.972222 0 0 0\n"
             <<"0 0.972222 0 0\n"
             <<"0 0 1 0\n";

    std::cout<<"P2: "<<P2<<std::endl;
    std::cout<<"P2 for fist pose should be\n"
             <<" -0.957966 0.165734 -0.00707496 0.0774496\n"
             <<"0.164089 0.952816 0.102143 0.967341\n"
             <<"0.0250416 0.102292 -0.994439 0.0605768\n";
    

    //下面这一块整个看不懂啊。。。。贼几把难受！
    /* 构造A矩阵 */
    math::Matrix<double, 4, 4> A;
    // 对A的每一列分别进行赋值
    for(int i=0; i<4; i++){
        // 第1个点
        A(0, i) = p1[0]*P1(2, i) - P1(0, i);
        A(1, i) = p1[1]*P1(2, i) - P1(1, i);

        // 第2个点
        A(2, i) = p2[0]*P2(2, i) - P2(0, i);
        A(3, i) = p2[1]*P2(2, i) - P2(1, i);
    }

    std::cout<<"A: "<<std::endl;
    std::cout<<"A for first pose should be:\n"
             <<"-0.972222 0 0.180123 0\n"
             <<"-0 -0.972222 -0.156584 -0\n"
             <<"0.963181 -0.14443 -0.200031 -0.0648336\n"
               <<"-0.164975 -0.956437 -0.0669352 -0.969486\n";

    math::Matrix<double, 4, 4> V;
    math::matrix_svd<double, 4, 4> (A, nullptr, nullptr, &V);
    math::Vec3d X;
    X[0] = V(0, 3)/V(3, 3);
    X[1] = V(1, 3)/V(3, 3);
    X[2] = V(2, 3)/V(3, 3);
    std::cout<<"--------------------------------------"<<std::endl;
    std::cout<<"X is:\n"<<X[0]<<" "<<X[1] <<" "<<X[2]<<std::endl;
    std::cout<<"X for first pose should be:\n"
             <<"3.2043116948585566 -2.7710180887818652 17.195578538234088\n";
    return X;


}


/*
* description 判断相机姿态是否正确，方法是计算三维点在两个相机中的坐标，要求其
* z坐标为正，即三位点同时位于两个相机前方
* @param p3d
* @param R1
* @param R2
* @param t1
* @param t2
*/
bool is_correct_pose(const math::Vec3d p3d,const math::Matrix3d& R1,
                    const math::Matrix3d& R2,const math::Vec3d& t1,
                    const math::Vec3d& t2,const math::Matrix3d& K1,
                    const math::Matrix3d& K2)
{
    math::Vector<double,3> x1 = R1 * p3d + t1;
    math::Vector<double,3> x2 = R2 * p3d + t2;

    return x1[2] > 0 && x2[2] > 0;

}

pose pose_from_fundamental(const FundamentalMatrix& F,const math::Matrix3d& K1,
                            const math::Matrix3d& K2,const math::Vec2d& p1,
                            const math::Vec2d& p2)
{
    EssentialMatrix E = K2.transposed() * F * K1;

    std::cout<<"EssentialMatrix result is "<<E<<std::endl;
    std::cout<<"EssentialMatrix should be: \n"
             <<"-0.00490744 -0.0146139 0.34281\n"
             <<"0.0212215 -0.000748851 -0.0271105\n"
             <<"-0.342111 0.0315182 -0.00552454\n";
    
     /* 本质矩阵求解的是相机之间的相对姿态，第一个相机姿态可以设置为[I|0], 第二个相机的姿态[R|t]
     * 可以通过对本质矩阵进行分解来求得, E=U*S*V',其中S是进行尺度归一化之后是diag(1,1,0)
     */

    math::Matrix3d W(0.0);
    W(0,1) = -1.0; W(1,0) = 1.0; W(2,2) = 1.0;
    math::Matrix3d Wt(0.0);
    Wt(0,1) = 1.0; Wt(1,0) = -1.0; Wt(2,2) = 1.0;

    math::Matrix3d U,S,V;
    math::matrix_svd(E,&U,&S,&V);

    //保证旋转矩阵det(R) = 1 (instead of -1).
    if(math::matrix_determinant(U)<0.0)
    {
        for(int i = 0; i<3; i++)
        {
            U(i,2) = -U(i,2);
        }
    }
    if(math::matrix_determinant(V)<0.0)
        for(int i = 0; i<3; i++)
            V(i,2) = -V(i,2);

    //相机姿态一共有4种情况
    V = V.transposed();
    std::vector<std::pair<math::Matrix3d,math::Vec3d>> poses(4);
    poses[0].first = U * W * V;
    poses[1].first = U * W * V;
    poses[2].first = U * Wt * V;
    poses[3].first = U * Wt * V;

    poses[0].second = U.col(2);
    poses[1].second = -U.col(2);
    poses[2].second = U.col(2);
    poses[3].second = -U.col(2);

        std::cout<<"Result of 4 candidate camera poses shoule be \n"
    <<"R0:\n"
      <<"-0.985336 0.170469 -0.0072771\n"
     <<"0.168777 0.980039 0.105061\n"
     <<"0.0250416 0.102292 -0.994439\n"
     <<"t0:\n"
     <<" 0.0796625 0.99498 0.0605768\n"
     <<"R1: \n"
     <<"-0.985336 0.170469 -0.0072771\n"
     <<"0.168777 0.980039 0.105061\n"
     <<"0.0250416 0.102292 -0.994439\n"
     <<"t1:\n"
     <<"-0.0796625 -0.99498 -0.0605768\n"
     <<"R2: \n"
     <<"0.999827 -0.0119578 0.0142419\n"
     <<"0.0122145 0.999762 -0.0180719\n"
     <<"-0.0140224 0.0182427 0.999735\n"
     <<"t2:\n"
     <<"0.0796625 0.99498 0.0605768\n"
     <<"R3: \n"
     <<"0.999827 -0.0119578 0.0142419\n"
     <<"0.0122145 0.999762 -0.0180719\n"
     <<"-0.0140224 0.0182427 0.999735\n"
     <<"t3: \n"
     <<"-0.0796625 -0.99498 -0.0605768\n";

     //第一个相机的旋转矩阵R1设置为单位矩阵，平移向量t1设置为0
     math::Matrix3d R1;
     math::matrix_set_identity(&R1);
     math::Vec3d t1;
     t1.fill(0.0);

     //判断姿态是否合理
     bool flag[4];
     for(int i = 0; i<4; i++)
     {
         math::Vec3d p3d = triangulation(p1,p2,K1,K2,t1,poses[i].second,R1,poses[i].first);
         flag[i] = is_correct_pose(p3d,R1,poses[i].first,t1,poses[i].second,K1,K2);

     }
    //找到正确姿态
    math::Matrix3d R;
    math::Vec3d t;
    if(flag[0] || flag[1] || flag[2] || flag[3])
    {
        for(int i = 0; i<4; i++)
        {
            if(flag[i])
            {
                R = poses[i].first;
                t = poses[i].second;
                std::cout<<"-----------------\n"
                        <<"correct pose found!"<<std::endl;
            }
                    
        }
    }

    pose corr_pose;
    corr_pose.first = R;
    corr_pose.second = t;

    return corr_pose;

}
