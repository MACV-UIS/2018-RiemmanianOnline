#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include <opencv2/opencv.hpp>

#include "../Resources/GeneralFunctions.h"
#include "../Resources/Matrixmanage.h"


//#include "MatricialTrajectoryManage.h"

using namespace cv;
using namespace std;


//
//int num_traj_valides(float** sequenceKin_k, int num_tot_traj_k)
//{
//
//    int cont=0;
//    for(int j=0; j<num_tot_traj_k; j++)
//        if(sequenceKin_k[j][2]>2) cont++;
//
//    return  cont;
//}



void ToComputeFrameApp( float** MatTraj, float** MatAppTraj, int num_tot_traj_k, int heightFrame_k,
                        int widthFrame_k, string pathFileFrames, int frame_index, bool WRTT) // WRTT = with respect to trajectories
{
    float** f_label = ToCreateMatrix2D(widthFrame_k, heightFrame_k);
    ToInitMatrix2D(f_label, widthFrame_k, heightFrame_k);


        /// calcular el frame appareance
        string str_zeros= "0000000";
        if(frame_index>9 && frame_index<100)
        {
            str_zeros = "000000";
        }
        else if(frame_index>99 && frame_index<1000)
        {
            str_zeros = "00000";
        }

    stringstream ss_numframe;
    ss_numframe << frame_index;
    string frame_name = pathFileFrames + str_zeros + ss_numframe.str() + ".png";
    //cout<< "frame_name: " <<frame_name << endl;
   // string frame_name = pathFileFrames  + ss_numframe.str() + ".png";


    Mat src = imread(frame_name.c_str());
    Mat src_gray, grad;
    Mat grad_x, grad_y, grad_xx, grad_yy;
    Mat abs_grad_x, abs_grad_y, abs_grad_xx, abs_grad_yy;
    Mat abs_dst, dst; //Laplacian

    int kernel_size = 3;
    int scale = 1;
    int delta = 0;
    int ddepth = CV_16S;
    ///Gaussian Blur
    GaussianBlur( src, src, Size(kernel_size,kernel_size), 0, 0, BORDER_DEFAULT );
    /// Convert it to gray
    cvtColor( src, src_gray, CV_BGR2GRAY );
    /// Gradient X
    Sobel( src_gray, grad_x, ddepth, 1, 0, 3, scale, delta, BORDER_DEFAULT );
    /// Gradient Y
    Sobel( src_gray, grad_y, ddepth, 0, 1, 3, scale, delta, BORDER_DEFAULT );
    /// Gradient XX
    Sobel( src_gray, grad_xx, ddepth, 2, 0, 3, scale, delta, BORDER_DEFAULT );
    /// Gradient YY
    Sobel( src_gray, grad_yy, ddepth, 0, 2, 3, scale, delta, BORDER_DEFAULT );
    /// Apply Laplace function
    Laplacian( src_gray, dst, ddepth, kernel_size, scale, delta, BORDER_DEFAULT );
    ///Compute Absoulete values and scales
    convertScaleAbs( grad_x, abs_grad_x );
    convertScaleAbs( grad_y, abs_grad_y );
    convertScaleAbs( grad_xx, abs_grad_xx );
    convertScaleAbs( grad_yy, abs_grad_yy );
    convertScaleAbs( dst, abs_dst );

    //cout<< " num_tot_traj_k: "<< num_tot_traj_k<< " widthFrame_k: "<< widthFrame_k << " heightFrame_k: "<< heightFrame_k <<endl;
    for(int j=0; j< num_tot_traj_k; j++)
    {
        if(MatTraj[j][2] > 2)
        {
            //        cout << " x_val-> "<< (int)sequenceKin_k[j][0] << " y_val-> "<< (int)sequenceKin_k[j][1] << " j-> "<< j << endl;
            int x_val=0, y_val=0;
            if((int)MatTraj[j][0]>=widthFrame_k)
            {
                x_val = (int)widthFrame_k-1; // because dense  improved trajectories are based on interpolation and can exceeed the maximum
            }
            else
            {
                x_val = (int)MatTraj[j][0];
            }

            if((int)MatTraj[j][1]>=heightFrame_k)
            {
                y_val = (int)heightFrame_k-1;
            }
            else
            {
                y_val = (int)MatTraj[j][1];
            }

            f_label[x_val][y_val] =1;

            ///load features ------------------
            MatAppTraj[j][0] = x_val;  /// x
            MatAppTraj[j][1] = y_val;  /// y
            MatAppTraj[j][2] = (float)abs(grad_x.at<short>(y_val,x_val));  /// |Ix|
            MatAppTraj[j][3] = (float)abs(grad_y.at<short>(y_val,x_val));  /// |Iy|
            MatAppTraj[j][4] = (float)sqrt(pow(MatAppTraj[j][2],2) + pow(MatAppTraj[j][3],2)); /// sqrt(grad_x**2 + grad_y**2)
            MatAppTraj[j][5] = (float)abs(grad_xx.at<short>(y_val,x_val));  /// |Ixx|
            MatAppTraj[j][6] = (float)abs(grad_yy.at<short>(y_val,x_val));  /// |Iyy|
            /// Derivada cruzada /// |Ixy| ?? to add
            MatAppTraj[j][7] = (float)fastAtan2(MatAppTraj[j][2], MatAppTraj[j][3]);  /// atan |Ix|/|Iy|
//            cout<< endl;
//            cout<< " cargo! "<< endl;
        }
    }
////
//    addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad );
//    imshow( "sobel", grad );
//    waitKey(0);


//    imshow( "laplace", abs_dst );

    src.release();
    src_gray.release();
    grad.release();
    grad_x.release();
    grad_y.release();
    grad_xx.release();
    grad_yy.release();
    abs_grad_x.release();
    abs_grad_y.release();
    abs_grad_xx.release();
    abs_grad_yy.release();
    abs_dst.release();
    dst.release();

    ToEliminateMatrix2D(f_label, widthFrame_k, heightFrame_k);


}



void ToLoadMatTrajAppFeatures(int star_seq, int end_seq, string pathFileTrajectories, string pathFileFrames, int height_frame,
                              int width_frame, int num_tot_traj, int num_scales_traj, float*** MatTraj, float*** MatAppTraj, bool WRTT)
{


// leer traj per frame
// mirar si hay mas que N traj
// calcular appareance features and only left the correspondinf to the trajectories
// store in MatAppTraj

    for(int frame_index=1; frame_index< end_seq; frame_index++)
    {
//        int traj_used = num_traj_valides(MatTraj[frame_index],  num_tot_traj);
//        if(traj_used>100)
//        {
            ToComputeFrameApp( MatTraj[frame_index], MatAppTraj[frame_index], num_tot_traj,
                               height_frame, width_frame, pathFileFrames,  frame_index,  WRTT);
        //} //end if(traj_used>100)
    } // end for(int frame_index=3; frame_index< numTotalFrames_k;


}
