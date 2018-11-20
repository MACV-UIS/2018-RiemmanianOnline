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

#include "RecursiveCov.h"



#include "../Resources/GeneralFunctions.h"
#include "../Resources/Matrixmanage.h"
#include "../TrajManage/MatricialTrajectoryManage.h"


#include "../Gradien_descent_Mean/Interface_NR3_Float.h"
#include "../Fast_Covariance/Fast_Covariance.h"




int num_traj_valides(float** sequenceKin_k, int num_tot_traj_k)
{

    int cont=0;
    for(int j=0; j<num_tot_traj_k; j++)
        if(sequenceKin_k[j][2]>2) cont++;

    return  cont;
}





//float val_mean( float val_m, float val, float alpha){
//    return alpha*val + (1-alpha)*(val_m);
//}
void compute_Rec_Cov_Mean(float** cov_mean_prev, float ** cov_current, float ** cov_mean_actual, float  alpha_val, float distAB, int NFEATURES){

float** LoG_Cov = ToCreateMatrix2D(NFEATURES, NFEATURES);
ToInitMatrix2D(LoG_Cov, NFEATURES, NFEATURES);
float_Log_A_resp_I(cov_current, LoG_Cov,  NFEATURES, distAB);

for(int i= 0; i<NFEATURES; i++)
    for(int j=0; j<NFEATURES; j++)
        cov_mean_actual[i][j] = alpha_val*LoG_Cov[i][j] + (1-alpha_val)*cov_mean_prev[i][j];


ToEliminateMatrix2D(LoG_Cov, NFEATURES, NFEATURES);

}

//
//float val_std( float std_m, float val_m, float val, float alpha){
//    return (float)sqrt(pow(val_m - alpha*val, 2) + (1-alpha)*(std_m));
//}

void compute_Rec_Cov_Var(float** cov_var_prev,  float ** cov_mean_actual,  float ** cov_current, float ** cov_var_actual, float  alpha_val, float distAB, int NFEATURES){

float** LoG_Cov = ToCreateMatrix2D(NFEATURES, NFEATURES);
ToInitMatrix2D(LoG_Cov, NFEATURES, NFEATURES);
float_Log_A_resp_I(cov_current, LoG_Cov,  NFEATURES, distAB);

for(int i= 0; i<NFEATURES; i++)
    for(int j=0; j<NFEATURES; j++)
        cov_var_actual[i][j] =   pow(cov_mean_actual[i][j] - alpha_val*LoG_Cov[i][j], 2)   + (1-alpha_val)*cov_var_prev[i][j];


ToEliminateMatrix2D(LoG_Cov, NFEATURES, NFEATURES);

}

float max_v( float  val_1 , float val_2 ){
 if(val_1> val_2)return val_1;else return val_2;
}
float min_v( float  val_1 , float val_2 ){
 if(val_1< val_2)return val_1;else return val_2;
}

//float val_minimun(float val_min, float angle, float alpha){
// return alpha*angle + (1-alpha)*min_value(angle, val_min);
//}

void compute_Rec_Cov_Min(float** cov_min_prev, float ** cov_current, float ** cov_min_actual, float  alpha_val, float distAB, int NFEATURES){

float** LoG_Cov = ToCreateMatrix2D(NFEATURES, NFEATURES);
ToInitMatrix2D(LoG_Cov, NFEATURES, NFEATURES);
float_Log_A_resp_I(cov_current, LoG_Cov,  NFEATURES, distAB);

for(int i= 0; i<NFEATURES; i++)
    for(int j=0; j<NFEATURES; j++)
        cov_min_actual[i][j] = alpha_val*LoG_Cov[i][j] + (1-alpha_val)* min_v(LoG_Cov[i][j], cov_min_prev[i][j]);


ToEliminateMatrix2D(LoG_Cov, NFEATURES, NFEATURES);

}

//float val_maximun(float val_max, float angle, float alpha){
// return alpha*angle + (1-alpha)*max_value(angle, val_max);
//}

void compute_Rec_Cov_Max(float** cov_max_prev, float ** cov_current, float ** cov_max_actual, float  alpha_val, float distAB, int NFEATURES){

float** LoG_Cov = ToCreateMatrix2D(NFEATURES, NFEATURES);
ToInitMatrix2D(LoG_Cov, NFEATURES, NFEATURES);
float_Log_A_resp_I(cov_current, LoG_Cov,  NFEATURES, distAB);

for(int i= 0; i<NFEATURES; i++)
    for(int j=0; j<NFEATURES; j++)
        cov_max_actual[i][j] = alpha_val*LoG_Cov[i][j] + (1-alpha_val)* max_v(LoG_Cov[i][j], cov_max_prev[i][j]);


ToEliminateMatrix2D(LoG_Cov, NFEATURES, NFEATURES);

}



void compute_Rec_Cov_Dif(float** cov_max, float ** cov_min, float ** cov_diff_actual, float  alpha_val, float distAB, int NFEATURES){



for(int i= 0; i<NFEATURES; i++)
    for(int j=0; j<NFEATURES; j++)
        cov_diff_actual[i][j] = cov_max[i][j]  - cov_min[i][j];


}








//

void Intercambiar_num(float** K_near, int i, int j)
{
    float tmp_dist;
    float tmp_numFr;
    float tmp_index;

    tmp_dist   = K_near[i][1];
    tmp_numFr  = K_near[i][0];
    tmp_index  = K_near[i][2];

    K_near[i][1] = K_near[j][1];
    K_near[i][0] = K_near[j][0];
    K_near[i][2] = K_near[j][2];

    K_near[j][1] = tmp_dist;
    K_near[j][0] = tmp_numFr;
    K_near[j][2] = tmp_index;

}


void quickSort_inv(float**  K_near, int left, int right)
{

    int i = left, j = right;
    float pivot = K_near[(left + right) / 2][1];
    //cout<<" pivot "<< pivot << endl;

    /** Partición */
    while (i <= j)
    {
        while (K_near[i][1] > pivot)i++;
        while (K_near[j][1] < pivot)j--;
        if (i <= j)
        {
            Intercambiar_num(K_near, i, j);
            i++;
            j--;
        }
    }

    /** Recursividad */
    if (left < j)
    {
        quickSort_inv(K_near, left, j);
    }
    if (i < right)
    {
        quickSort_inv(K_near, i, right);
    }



}


void Copy_Cov(float** cov_source, float** cov_des, int NFEATURES)
{
    for(int i=0; i<NFEATURES; i++)
        for(int j=0; j<NFEATURES; j++)
            cov_des[i][j] = cov_source[i][j];

}


void Recursive_mean(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                    vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha)
{

    // The experiments were carried out with LED
    ofstream outdata_mean((name_file +"mean.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_mean.precision(10);
    outdata_mean<<  fixed;

    ofstream outdata_var((name_file + "var.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_var.precision(10);
    outdata_var<<  fixed;

    ofstream outdata_min((name_file + "min.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_min.precision(10);
    outdata_min<<  fixed;

    ofstream outdata_max((name_file + "max.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_max.precision(10);
    outdata_max<<  fixed;

    ofstream outdata_der((name_file + "der.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_der.precision(10);
    outdata_der<<  fixed;


    float** submat = ToCreateMatrix2D(NFEATURES, NFEATURES);
    float* EigenValSym = ToCreateMatrix1D(NFEATURES);

    for(int k=0; k< vec_videos.size(); k++)
    {



        cout<< k<< " )  processing video: "<< vec_videos[k]  << endl;
        string path_textTrajectories_k = str_path_trajectories + vec_videos[k] + ".txt";



        int numTotalFrames_k, heightFrame_k, widthFrame_k, num_tot_traj_k;
        infoVideoFromTrajecAllScales(path_textTrajectories_k, numTotalFrames_k,
                                     heightFrame_k, widthFrame_k, num_tot_traj_k, 1); //3 pqrq ENSTA TRAJ

        int tot_seq_matrix_k = numTotalFrames_k+1; // mas dos porque en las dos ultimas posiciones se escribe el inicion y final de la trajectoria. Cabezera info.
        int tot_feat_matrix_k = 12;//((int)featLabels.size()*(int)alphaScales.size())+3; 9 features + x,y,z

        cout<< "<----- creating mat: ["<< tot_seq_matrix_k<<"]["<< num_tot_traj_k<<"]["<< tot_feat_matrix_k<<"]    ------>"<<endl;
        float*** sequenceKin_k = ToCreateMatrix3D(tot_seq_matrix_k, num_tot_traj_k, tot_feat_matrix_k);
        ToInitMatrix3D(sequenceKin_k, tot_seq_matrix_k, num_tot_traj_k, tot_feat_matrix_k);

        ToLoadMatTrajOnlyKin(0, numTotalFrames_k,path_textTrajectories_k,
                             heightFrame_k, widthFrame_k, num_tot_traj_k, 1, sequenceKin_k);

        cout<< " ---------------- numero de Frames -> "<< numTotalFrames_k << endl;

        float **covmat_prev = ToCreateMatrix2D(NFEATURES, NFEATURES); //Covariance Matrix buffer
        ToInitMatrix2D(covmat_prev  , NFEATURES, NFEATURES);
        int cont_frameCov =0;

        int cont_tot_DifCov =0;
        float *** arr_cov  = ToCreateMatrix3D(numTotalFrames_k, NFEATURES, NFEATURES);
        float**   arr_info_first = ToCreateMatrix2D(numTotalFrames_k, 3);

        ToInitMatrix3D(arr_cov,numTotalFrames_k, NFEATURES, NFEATURES );
        ToInitMatrix2D(arr_info_first, numTotalFrames_k, 3);


        for(int frame_index=3; frame_index< numTotalFrames_k; frame_index++)
        {
            int traj_used = num_traj_valides(sequenceKin_k[frame_index],  num_tot_traj_k);
            //cout << " numero de frame: "<< frame_index << " numTotalFrames_k: "<< numTotalFrames_k << " trajectories used "<< traj_used << endl;
            if(traj_used>100)
            {
                float **covmat_K = ToCreateMatrix2D(NFEATURES, NFEATURES); //Covariance Matrix buffer
                ToInitMatrix2D(covmat_K, NFEATURES, NFEATURES);
                //cout<<" entro cov"<< endl;
                ToComputeCovmat(covmat_K, sequenceKin_k[frame_index],  NFEATURES, num_tot_traj_k, heightFrame_k, widthFrame_k);

                //cin.ignore();
                //ToPrintMatrix2D(covmat_K, NFEATURES, NFEATURES);

                //cout<<" salio cov "<< endl;
                ToInitMatrix2D(submat, NFEATURES, NFEATURES);
                ToInitMatrix1D(EigenValSym, NFEATURES);
                //cout<< " entro eig"<<endl;
                float_eigValues(covmat_K, EigenValSym, NFEATURES);
                //cout<< " salio eig "<< endl;
                if(symPos_label(EigenValSym, NFEATURES)== true)
                {

                    if (cont_frameCov>1)
                    {
                        float val_dist =0;
                        val_dist =  float_Frobenius_Dist(covmat_K, covmat_prev, NFEATURES);

                        arr_info_first[cont_tot_DifCov][0] = frame_index;
                        arr_info_first[cont_tot_DifCov][1] = val_dist;
                        arr_info_first[cont_tot_DifCov][2] = cont_tot_DifCov;




                        for( int m=0; m<NFEATURES; m++)
                        {
                            for(int n=0; n<NFEATURES; n++)
                            {
                                covmat_prev[m][n] = covmat_K[m][n];
                                arr_cov[cont_tot_DifCov][m][n] = covmat_K[m][n];

                            }
                        }

                        cont_tot_DifCov++;


                    }
                    else
                    {

                        for( int m=0; m<NFEATURES; m++)
                            for(int n=0; n<NFEATURES; n++)
                                covmat_prev[m][n] = covmat_K[m][n];

                    }

                    cont_frameCov++;


                }


                ToEliminateMatrix2D(covmat_K, NFEATURES, NFEATURES); // no la pasa porque es del vector
            }

        }

        ToEliminateMatrix2D(covmat_prev, NFEATURES, NFEATURES);
        ToEliminateMatrix3D(sequenceKin_k, tot_seq_matrix_k, num_tot_traj_k, tot_feat_matrix_k);


        //cont_tot_DifCov++;
        cout<< " cont_tot_DifCov: "<< cont_tot_DifCov << endl;


///-----------------------------------------------
        if(cont_tot_DifCov >num_mean)
        {

            //cout<< " entro "<< endl;
            int label_k = atoi(returnLabelAction_KTH(vec_videos[k], vec_activities ).c_str()); //-1 para que la primera sea cero

            outdata_mean<< label_k<< " ";
            outdata_var<< label_k<< " ";
            outdata_min<< label_k<< " ";
            outdata_max<< label_k<< " ";
            outdata_der<< label_k<< " ";



//            //cout<< " hizo lo de  "<< endl;
            float*** array_mean = ToCreateMatrix3D(cont_tot_DifCov, NFEATURES, NFEATURES);
            ToInitMatrix3D(array_mean, cont_tot_DifCov, NFEATURES, NFEATURES);

            float*** array_var = ToCreateMatrix3D(cont_tot_DifCov, NFEATURES, NFEATURES);
            ToInitMatrix3D(array_var, cont_tot_DifCov, NFEATURES, NFEATURES);

            float*** array_min = ToCreateMatrix3D(cont_tot_DifCov, NFEATURES, NFEATURES);
            ToInitMatrix3D(array_min, cont_tot_DifCov, NFEATURES, NFEATURES);

            float*** array_max = ToCreateMatrix3D(cont_tot_DifCov, NFEATURES, NFEATURES);
            ToInitMatrix3D(array_max, cont_tot_DifCov, NFEATURES, NFEATURES);

            float*** array_dif = ToCreateMatrix3D(cont_tot_DifCov, NFEATURES, NFEATURES);
            ToInitMatrix3D(array_dif, cont_tot_DifCov, NFEATURES, NFEATURES);


            float**   arr_info_mean = ToCreateMatrix2D(numTotalFrames_k, 2);
            ToInitMatrix2D(arr_info_mean, numTotalFrames_k, 2);


            float alpha_val = pow(2, -1*init_alpha);
            cout<< "alpha val: "<<alpha_val << endl;
            //Copy_Cov(arr_cov[0], array_mean[0],  NFEATURES);
            float distAB=0;
            float_Log_A_resp_I(arr_cov[0], array_mean[0],  NFEATURES, distAB);
            float_Log_A_resp_I(arr_cov[0], array_max[0],  NFEATURES, distAB);
            float_Log_A_resp_I(arr_cov[0], array_dif[0],  NFEATURES, distAB);


            for(int i_mean = 1; i_mean < cont_tot_DifCov; i_mean++)
            {

                compute_Rec_Cov_Mean(array_mean[i_mean-1], arr_cov[i_mean], array_mean[i_mean],   alpha_val,  distAB,   NFEATURES);
                compute_Rec_Cov_Var(array_var[i_mean-1], array_mean[i_mean], arr_cov[i_mean], array_var[i_mean],   alpha_val,  distAB,  NFEATURES);
                compute_Rec_Cov_Min(array_min[i_mean-1], arr_cov[i_mean], array_min[i_mean], alpha_val, distAB, NFEATURES);
                compute_Rec_Cov_Max(array_max[i_mean-1], arr_cov[i_mean], array_max[i_mean], alpha_val, distAB, NFEATURES);
                compute_Rec_Cov_Dif(array_max[i_mean], array_min[i_mean], array_dif[i_mean], alpha_val,  distAB,  NFEATURES);

                arr_info_mean[i_mean][0] = i_mean;
            }


            int cont_P =1;
            for(int i=0; i< NFEATURES; i++){
                for(int j=i; j< NFEATURES; j++){
                        outdata_mean<< cont_P << ":"<<array_mean[cont_tot_DifCov-1][i][j] << " ";
                        outdata_var<< cont_P++ << ":"<<array_var[cont_tot_DifCov-1][i][j] << " ";
                        outdata_min<< cont_P++ << ":"<<array_min[cont_tot_DifCov-1][i][j] << " ";
                        outdata_max<< cont_P++ << ":"<<array_max[cont_tot_DifCov-1][i][j] << " ";
                        outdata_der<< cont_P++ << ":"<<array_dif[cont_tot_DifCov-1][i][j] << " ";
                        cont_P++;

                        }}
            outdata_mean<<endl;
            outdata_var<<endl;
            outdata_min<<endl;
            outdata_max<<endl;
            outdata_der<<endl;


///---------------------------------------------------
            ToEliminateMatrix3D(array_mean, cont_tot_DifCov, NFEATURES, NFEATURES);
            ToEliminateMatrix3D(array_var, cont_tot_DifCov, NFEATURES, NFEATURES);
            ToEliminateMatrix3D(array_min, cont_tot_DifCov, NFEATURES, NFEATURES);
            ToEliminateMatrix3D(array_max, cont_tot_DifCov, NFEATURES, NFEATURES);
            ToEliminateMatrix3D(array_dif, cont_tot_DifCov, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(arr_info_mean, numTotalFrames_k, 2);
        }
        ToEliminateMatrix3D(arr_cov, numTotalFrames_k, NFEATURES, NFEATURES);
        ToEliminateMatrix2D(arr_info_first, numTotalFrames_k, 3);
    }


    ToEliminateMatrix2D(submat, NFEATURES, NFEATURES);
    ToEliminateMatrix1D(EigenValSym, NFEATURES);


            outdata_mean.close();
            outdata_var.close();
            outdata_min.close();
            outdata_max.close();
            outdata_der.close();



}

bool symPos_label(float* EigenValSym, int NFEATURES)
{

    bool label = true;
    for(int z=0; z<NFEATURES; z++)
    {
        if( EigenValSym[z] <= 0.001)
        {
            label= false;
        }
    }

    return label;

}


void ToComputeCovmat(float** covmat_K, float** sequenceKin_k, int NFEATURES,int num_tot_traj_k, int heightFrame_k, int widthFrame_k)
{

    float ***f_mat = ToCreateMatrix3D(NFEATURES, widthFrame_k, heightFrame_k);
    ToInitMatrix3D(f_mat, NFEATURES, widthFrame_k, heightFrame_k);

    float** f_label = ToCreateMatrix2D(widthFrame_k, heightFrame_k);
    ToInitMatrix2D(f_label, widthFrame_k, heightFrame_k);

     //cout<< " num_tot_traj_k: "<< num_tot_traj_k<< " widthFrame_k: "<< widthFrame_k << " heightFrame_k: "<< heightFrame_k <<endl;
    for(int j=0; j< num_tot_traj_k; j++)
    {
        if(sequenceKin_k[j][2] > 2)
        {

    //        cout << " x_val-> "<< (int)sequenceKin_k[j][0] << " y_val-> "<< (int)sequenceKin_k[j][1] << " j-> "<< j << endl;
            int x_val=0, y_val=0;
            if((int)sequenceKin_k[j][0]>=widthFrame_k){
                x_val = (int)widthFrame_k-1; // because dense  improved trajectories are based on interpolation and can exceeed the maximum
                }else{
                x_val = (int)sequenceKin_k[j][0];
                }

            if((int)sequenceKin_k[j][1]>=heightFrame_k){
                y_val = (int)heightFrame_k-1;
                }else{
                y_val = (int)sequenceKin_k[j][1];
                }

           f_label[x_val][y_val] =1;
            for(int feat_index=0; feat_index<NFEATURES; feat_index++)
            {
  //              cout<< sequenceKin_k[j][feat_index] <<" ";
                f_mat[feat_index][x_val][y_val] = sequenceKin_k[j][feat_index];


            }
//            cout<< endl;
//            cout<< " cargo! "<< endl;
        }
    }

    //cout<< " lleno fmat "<< endl;

    //cout<< " cont_feat_test: "<< cont_feat_test << endl;
    int nmatrices = ((NFEATURES * (NFEATURES + 1)) / 2 );

    float*** product_matrix_Q = ToCreateMatrix3D(nmatrices, widthFrame_k, heightFrame_k); //Feature Matrix
    ToInitMatrix3D(product_matrix_Q, nmatrices,  widthFrame_k, heightFrame_k); //10 x 10 region to measure

    computeIntegral_matrices(f_mat,  NFEATURES,  heightFrame_k,  widthFrame_k, product_matrix_Q );

    // cout<<" calculo bien las integral para ir a cov "<< endl;

    covariance(0, 0, widthFrame_k-1, heightFrame_k-1, f_mat, product_matrix_Q, covmat_K, f_label,  NFEATURES );
    ToEliminateMatrix3D(f_mat, NFEATURES,  widthFrame_k, heightFrame_k);
    ToEliminateMatrix3D(product_matrix_Q, nmatrices, widthFrame_k, heightFrame_k);
    ToEliminateMatrix2D(f_label, widthFrame_k, heightFrame_k);

//    ToPrintMatrix2D(covmat_K, NFEATURES, NFEATURES);
//    cin.ignore();
    // cout<<  " \n ******* cov matrix computed !************** "<< endl;


    //ToPrintMatrix2D(covmat_K, NFEATURES, NFEATURES);
    //cin.ignore();
}
