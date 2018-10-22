#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <vector>

#include <cv.h>
#include <highgui.h>
#include <opencv2/opencv.hpp>

#include"KinematicFeatures.h"
#include "../Resources/Matrixmanage.h"
#include "GoMKByFeat.h"

using namespace cv;
using namespace std;

#define PI  3.14159265358979323846;

//--------------gloab variables ----------------------------------------------------------------------------------------------
const int components=3; // namely the mean, variance and weight

//------------- function declaration ------------------------------------------------------------------------------------------
void compute_Velxy_MoG(vector<float> &v_x, vector<float> &v_y,  vector<float> &vec_vel_x, vector<float> &vec_vel_y);
vector< vector<float> > updateMoGFeat_Original(int numberOfGaussians, float alphaMoG, float std_init, float w_init,
                                     float mean_init, vector<float> vec_feat);
 vector< vector<float> > initCurrentMoG(int numberOfGaussians, float std_init, float w_init , float mean_init);
 vector< vector<float> > sortMoGVector(vector< vector<float> > CurrentMoG, int numberOfGaussians);
 vector< vector<float> > normalizedMoGVector(vector< vector<float> > CurrentMoG, int numberOfGaussians);
//-----------------------------------------------------------------------------------------------------------------------------


vector< vector< float > > computeMoGFeat(vector<float> vect_traj_x, vector<float> vect_traj_y, vector<float> vect_traj_z,
                                         string feature_name, int numberOfGaussians){

    vector<float> feat_vect = computeKinFeatByVector(vect_traj_x, vect_traj_y, vect_traj_z, feature_name);
    vector< vector< float > > CurrentMoG(numberOfGaussians, vector<float>(3));
    float std_init = 30.0; float alphaMoG = 0.005; float w_init = 0.05; float mean_init = feat_vect[0];
    CurrentMoG = updateMoGFeat_Original(numberOfGaussians, alphaMoG, std_init, w_init, mean_init, feat_vect);
    return CurrentMoG;


}


void computeMoGPerTraj(vector<float> vect_traj_x, vector<float> vect_traj_y, vector<float> vect_traj_z,
                                         string feature_name, int numberOfGaussians, float*** matrixMoGTraj,
                                         float alphaMoG, float std_init, float w_init, float rho_init){

    vector<float> feat_vect = computeKinFeatByVector(vect_traj_x, vect_traj_y, vect_traj_z, feature_name);
    float mean_init = feat_vect[0];

//**************** 0. Init matrix *********************
    init_MatrixMoGTraj(matrixMoGTraj,  numberOfGaussians,  std_init,  w_init,  mean_init);


//    cout<< "t: "<< 0<< " ******************** Feat val: "<< feat_vect[0] <<endl;
//    for(int nGauss=0; nGauss<numberOfGaussians; nGauss++){
//    cout<<" mean: "<< matrixMoGTraj[0][nGauss][0] << " std: "<< matrixMoGTraj[1][nGauss][0] << " w: "<< matrixMoGTraj[2][nGauss][0] <<endl ;
//    }
//    cout<<endl;



    for(int t=1; t<feat_vect.size(); t++){


        int bad_wG=1; // bandera para que solo le sume 1 a la primera gaussiana que cumpla con la condicion.

//************************1. Asignar ***********************************
        for(int nGauss=0; nGauss<numberOfGaussians; nGauss++){

            float rho=rho_init;
//            if(matrixMoGTraj[2][nGauss][t]>0){
//             rho =alphaMoG/matrixMoGTraj[2][nGauss][t]; //
//            }else{
//             rho =rho_init;//alphaMoG/matrixMoGTraj[2][nGauss][t]; //
//            }


             if(val_abs(feat_vect[t] - matrixMoGTraj[0][nGauss][t-1]) <= 2.5*matrixMoGTraj[1][nGauss][t-1] && bad_wG==1 ){ // 2.5*matrixMoGTraj[1][nGauss][t-1] - 0.01
                 matrixMoGTraj[0][nGauss][t] = val_meanMoG(matrixMoGTraj[0][nGauss][t-1], feat_vect[t], rho ); // mean
                 matrixMoGTraj[1][nGauss][t] = val_stdMoG (matrixMoGTraj[1][nGauss][t-1], matrixMoGTraj[0][nGauss][t],  feat_vect[t], rho ); // std
                 matrixMoGTraj[2][nGauss][t] = val_meanMoG(matrixMoGTraj[2][nGauss][t-1] , 1.0f,alphaMoG); // w
                 bad_wG=0;
             }
             else{
                 matrixMoGTraj[0][nGauss][t] = matrixMoGTraj[0][nGauss][t-1]; // mean
                 matrixMoGTraj[1][nGauss][t] = matrixMoGTraj[1][nGauss][t-1]; // std
                 matrixMoGTraj[2][nGauss][t] = matrixMoGTraj[2][nGauss][t-1]; // w
                matrixMoGTraj[2][nGauss][t] = val_meanMoG(matrixMoGTraj[2][nGauss][t] , 0.0, alphaMoG); // w
             }
        }



//************************ 2. normalizar  ***********************************
         norm_MatrixMoGTraj(matrixMoGTraj,  numberOfGaussians,  t);
//************************ 3. ordenar     ***********************************
         sort_MatrixMoGTraj(matrixMoGTraj,  numberOfGaussians,  t);
//************************4. incluir si es el caso  ***********************************

        if(bad_wG==1){
            matrixMoGTraj[0][numberOfGaussians-1][t]= feat_vect[t]; //feature value
            matrixMoGTraj[1][numberOfGaussians-1][t]= std_init; // large value
            matrixMoGTraj[2][numberOfGaussians-1][t]= alphaMoG; // very small value
           // cout<<"*************no entro a ninguna Gaussiana *********" << endl;
            // si no acutalizo ninguna gaussiana entonces remplaza la [ultima
        }


///************** print the MoG ***********************

//    cout<< "t: "<< t<< " ******************** Feat val: "<< feat_vect[t] <<endl;
//    for(int nGauss=0; nGauss<numberOfGaussians; nGauss++){
//    cout<<" mean: "<< matrixMoGTraj[0][nGauss][t] << " std: "<< matrixMoGTraj[1][nGauss][t] << " w: "<< matrixMoGTraj[2][nGauss][t] <<endl ;
//    }
//    cout<<endl;
//    std::cin.ignore();
    }


//*************************************


}


void norm_MatrixMoGTraj(float*** matrixMoGTraj, int numberOfGaussians, int t){

        float pesoNor =0;
        for(int i=0; i<numberOfGaussians; i++){ pesoNor= pesoNor+  matrixMoGTraj[2][i][t];}
        for(int i=0; i<numberOfGaussians; i++){ if(pesoNor>0){matrixMoGTraj[2][i][t]=matrixMoGTraj[2][i][t]/pesoNor;}else{matrixMoGTraj[2][i][t]=0;}}

}



void init_MatrixMoGTraj(float*** matrixMoGTraj, int numberOfGaussians, float std_init, float w_init, float mean_init){

    for(int y=0; y< numberOfGaussians; y++){
            matrixMoGTraj[0][y][0]=mean_init;
            matrixMoGTraj[1][y][0]=std_init;
            matrixMoGTraj[2][y][0]=w_init;
    }

}


void sort_MatrixMoGTraj(float*** matrixMoGTraj, int numberOfGaussians, int t){

    int j,i=0;
    int pos_max=0;
    for(j=0; j<numberOfGaussians-1; j++){
         pos_max =j;
        for(i=j+1; i<numberOfGaussians; i++){
            if(matrixMoGTraj[2][j][t]/matrixMoGTraj[1][j][t] < matrixMoGTraj[2][i][t]/matrixMoGTraj[1][i][t]){
                pos_max = i;
            }
        }
        if(pos_max !=j){
            float temp_mean = matrixMoGTraj[0][j][t];
            float temp_std  = matrixMoGTraj[1][j][t];
            float temp_w    = matrixMoGTraj[2][j][t];

            matrixMoGTraj[0][j][t] = matrixMoGTraj[0][pos_max][t];
            matrixMoGTraj[1][j][t] = matrixMoGTraj[1][pos_max][t];
            matrixMoGTraj[2][j][t] = matrixMoGTraj[2][pos_max][t];

            matrixMoGTraj[0][pos_max][t] = temp_mean;
            matrixMoGTraj[1][pos_max][t] = temp_std;
            matrixMoGTraj[2][pos_max][t] = temp_w;
        }
    }


}


float val_meanMoG( float val_m, float val, float alpha){
    return alpha*val + (1-alpha)*(val_m);
}

float val_stdMoG( float sigma_m, float val_m, float val, float alpha){
    return (float)sqrt(alpha*pow(val_m - val, 2) + (1-alpha)*(sigma_m));
}


float val_euclnom_MOG(float val_ax, float val_ay){
    return sqrt(pow((val_ax),2) + pow((val_ay),2) );
}

float val_Distance1D_MOG(float val_a , float val_aprev){
    return val_a - val_aprev;

}

float val_abs(float val){

    if(val>=0){
        return val;
    }else{return (-1.0*val);}

}



vector< vector<float> > updateMoGFeat_Original(int numberOfGaussians, float alphaMoG, float std_init, float w_init,
                                     float mean_init, vector<float> vec_feat){

    vector< vector<float> > CurrentMoG(numberOfGaussians, vector<float>(3));
    CurrentMoG = initCurrentMoG( numberOfGaussians, std_init, w_init , mean_init);

    /*for(int nGauss=0; nGauss<numberOfGaussians; nGauss++){
             CurrentMoG[nGauss][0] = nGauss*360/numberOfGaussians;
    }*/


    for(int t=1; t<vec_feat.size(); t++){


        int bad_wG=1; // bandera para que solo le sume 1 a la primera gaussiana que cumpla con la condicion.

        for(int nGauss=0; nGauss<numberOfGaussians; nGauss++){
            float rho =0.05;//alphaMoG/CurrentMoG[nGauss][2];


             if(val_abs(vec_feat[t]-CurrentMoG[nGauss][0]) <= 0.01 && bad_wG==1 ){ //2.5*(CurrentMoG[nGauss][1])  ,
                 CurrentMoG[nGauss][0] = val_meanMoG(CurrentMoG[nGauss][0], vec_feat[t], rho ); // mean
                 CurrentMoG[nGauss][1] = val_stdMoG(CurrentMoG[nGauss][1], CurrentMoG[nGauss][0],  vec_feat[t], rho ); // std
                 CurrentMoG[nGauss][2] = val_meanMoG(CurrentMoG[nGauss][2] , 1.0,alphaMoG); // w
                 bad_wG=0;
             }
             else{
                CurrentMoG[nGauss][2] = val_meanMoG(CurrentMoG[nGauss][2] , 0.0,alphaMoG); // w
             }
        }

        CurrentMoG = normalizedMoGVector(CurrentMoG, numberOfGaussians);

        CurrentMoG = sortMoGVector(CurrentMoG, numberOfGaussians);//ordenar por el factor w/sigma
        if(bad_wG==1){
            CurrentMoG[numberOfGaussians-1][0]= vec_feat[t]; //feature value
            CurrentMoG[numberOfGaussians-1][1]= 30.0; // large value
            CurrentMoG[numberOfGaussians-1][2]= 0.005; // very small value
           // cout<<"*************no entro a ninguna Gaussiana *********" << endl;
            // si no acutalizo ninguna gaussiana entonces remplaza la [ultima
        }

    }
    return CurrentMoG;
}

vector< vector<float> > normalizedMoGVector(vector< vector<float> > CurrentMoG, int numberOfGaussians){
    float pesoNor =0;
    for(int i=0; i<numberOfGaussians; i++){ pesoNor= pesoNor+  CurrentMoG[i][2];}
    for(int i=0; i<numberOfGaussians; i++){ CurrentMoG[i][2]=CurrentMoG[i][2]/pesoNor;}
 return CurrentMoG;
}

vector< vector<float> > sortMoGVector(vector< vector<float> > CurrentMoG, int numberOfGaussians){
    int j,i=0;
    int pos_max=0;
    for(j=0; j<numberOfGaussians-1; j++){
         pos_max =j;
        for(i=j+1; i<numberOfGaussians; i++){
            if(CurrentMoG[j][2]/CurrentMoG[j][1] < CurrentMoG[i][2]/CurrentMoG[i][1]){
                pos_max = i;
            }
        }
        if(pos_max !=j){
            float temp_mean = CurrentMoG[j][0];
            float temp_std  = CurrentMoG[j][1];
            float temp_w    = CurrentMoG[j][2];

            CurrentMoG[j][0] = CurrentMoG[pos_max][0];
            CurrentMoG[j][1] = CurrentMoG[pos_max][1];
            CurrentMoG[j][2] = CurrentMoG[pos_max][2];

            CurrentMoG[pos_max][0] = temp_mean;
            CurrentMoG[pos_max][1] = temp_std;
            CurrentMoG[pos_max][2] = temp_w;
        }
    }

    return CurrentMoG;
}

 vector< vector<float> > initCurrentMoG(int numberOfGaussians, float std_init, float w_init , float mean_init){

        vector< vector<float> > CurrentMoG(numberOfGaussians, vector<float>(3));
        for(int nGauss=0; nGauss<numberOfGaussians; nGauss++){
             CurrentMoG[nGauss][0] = mean_init; // mean
             CurrentMoG[nGauss][1] = std_init; // std
             CurrentMoG[nGauss][2] = w_init; // w
        }
        return CurrentMoG;
}

