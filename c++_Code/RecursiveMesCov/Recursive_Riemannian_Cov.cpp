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

#include "Recursive_Riemannian_Cov.h"
#include "RecursiveCov.h"



#include "../Resources/GeneralFunctions.h"
#include "../Resources/Matrixmanage.h"
#include "../TrajManage/MatricialTrajectoryManage.h"


#include "../Gradien_descent_Mean/Interface_NR3_Float.h"
#include "../Fast_Covariance/Fast_Covariance.h"







//float val_mean( float val_m, float val, float alpha){
//    return alpha*val + (1-alpha)*(val_m);
//}
void compute_Rec_Riemmannian_Cov_Mean(float** cov_mean_prev, float ** cov_current, float ** cov_mean_actual, float  alpha_val, float distAB,  int NFEATURES){

float** LoG_Cov = ToCreateMatrix2D(NFEATURES, NFEATURES);
ToInitMatrix2D(LoG_Cov, NFEATURES, NFEATURES);

float_Log_A_resp_B_Alpha(cov_current, cov_mean_prev, LoG_Cov, NFEATURES, distAB, alpha_val);
//ToPrintMatrix2D(LoG_Cov, NFEATURES, NFEATURES);

float_Exp_A_resp_B_Fletcher(LoG_Cov, cov_mean_prev, cov_mean_actual,NFEATURES);
//ToPrintMatrix2D(cov_mean_actual, NFEATURES, NFEATURES);


ToEliminateMatrix2D(LoG_Cov, NFEATURES, NFEATURES);

}

//
//float val_std( float std_m, float val_m, float val, float alpha){
//    return (float)sqrt(pow(val_m - alpha*val, 2) + (1-alpha)*(std_m));
//}

void compute_Rec_Riemmannian_Cov_Var(float** cov_var_prev,  float ** cov_mean_actual,  float ** cov_current, float ** cov_var_actual, float  alpha_val, float distAB, int NFEATURES){

float** LoG_Cov = ToCreateMatrix2D(NFEATURES, NFEATURES);
ToInitMatrix2D(LoG_Cov, NFEATURES, NFEATURES);

float** LoG_var = ToCreateMatrix2D(NFEATURES, NFEATURES);
ToInitMatrix2D(LoG_var, NFEATURES, NFEATURES);

float_Log_A_resp_B_Variance(cov_current, cov_mean_actual, LoG_var, NFEATURES, distAB );
//ToPrintMatrix2D(LoG_var, NFEATURES, NFEATURES);

float_Log_A_resp_B_Alpha(LoG_var,  cov_var_prev, LoG_Cov, NFEATURES, distAB, alpha_val);

float_Exp_A_resp_B_Fletcher(LoG_Cov, cov_var_prev, cov_var_actual, NFEATURES);


ToEliminateMatrix2D(LoG_Cov, NFEATURES, NFEATURES);
ToEliminateMatrix2D(LoG_var, NFEATURES, NFEATURES);

}




void Recursive_Riemannian_mean(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                    vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha)
{

    // The experiments were carried out with LED
    ofstream outdata_mean((name_file +"mean.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_mean.precision(10);
    outdata_mean<<  fixed;

    ofstream outdata_var((name_file + "var.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_var.precision(10);
    outdata_var<<  fixed;


    ofstream outdata_var2((name_file + "var2.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_var2.precision(10);
    outdata_var2<<  fixed;

    // The experiments were carried out with LED
    ofstream outdata_meanSV((name_file +"meanSV.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_meanSV.precision(10);
    outdata_meanSV<<  fixed;

    ofstream outdata_varSV((name_file + "varSV.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_varSV.precision(10);
    outdata_varSV<<  fixed;

    ofstream outdata_var2SV((name_file + "var2SV.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_var2SV.precision(10);
    outdata_var2SV<<  fixed;





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
            outdata_var2<< label_k<< " ";


            outdata_meanSV<< label_k<< " ";
            outdata_varSV<< label_k<< " ";
            outdata_var2SV<< label_k<< " ";


//            //cout<< " hizo lo de  "<< endl;
            float*** array_mean = ToCreateMatrix3D(cont_tot_DifCov, NFEATURES, NFEATURES);
            ToInitMatrix3D(array_mean, cont_tot_DifCov, NFEATURES, NFEATURES);

            float*** array_var = ToCreateMatrix3D(cont_tot_DifCov, NFEATURES, NFEATURES);
            ToInitMatrix3D(array_var, cont_tot_DifCov, NFEATURES, NFEATURES);




            float**   arr_info_mean = ToCreateMatrix2D(numTotalFrames_k, 2);
            ToInitMatrix2D(arr_info_mean, numTotalFrames_k, 2);


            float alpha_val = pow(2, -1*init_alpha);
            cout<< "alpha val: "<<alpha_val << endl;
            Copy_Cov(arr_cov[0], array_mean[0],  NFEATURES);

            float distAB=0;
            Copy_Cov(arr_cov[0], array_mean[0],  NFEATURES);
            Copy_Cov(arr_cov[0], array_var[0],  NFEATURES);
            float_Log_A_resp_I(arr_cov[0], array_mean[0],  NFEATURES, distAB);



            for(int i_mean = 1; i_mean < cont_tot_DifCov; i_mean++)
            {

                compute_Rec_Riemmannian_Cov_Mean(array_mean[i_mean-1], arr_cov[i_mean], array_mean[i_mean],   alpha_val,  distAB,   NFEATURES);
                compute_Rec_Riemmannian_Cov_Var(array_var[i_mean-1], array_mean[i_mean], arr_cov[i_mean], array_var[i_mean],   alpha_val,  distAB,  NFEATURES);

                arr_info_mean[i_mean][0] = i_mean;
            }



            float** mean_Log = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(mean_Log, NFEATURES, NFEATURES);
            float_Log_A_resp_I(array_mean[cont_tot_DifCov-1], mean_Log,  NFEATURES, distAB);

            int cont_P =1;
            for(int i=0; i< NFEATURES; i++){
                for(int j=i; j< NFEATURES; j++){
                        outdata_mean<< cont_P << ":"<<mean_Log[i][j] << " ";

                        outdata_var<< cont_P << ":"<<array_var[cont_tot_DifCov-1][i][j] << " ";
            //outdata_var<<endl;
                        outdata_var2<< cont_P << ":"<<array_var[cont_tot_DifCov-1][i][j]/cont_tot_DifCov << " ";
            //outdata_var2<<endl;

            cont_P++;

            }}
            outdata_mean<<endl;
            outdata_var<<endl;
            outdata_var2<<endl;


            cont_P =1;
            for(int i=0; i< NFEATURES; i++){
                for(int j=i+1; j< NFEATURES; j++){
                        outdata_meanSV<< cont_P << ":"<<mean_Log[i][j] << " ";
                        outdata_varSV<< cont_P << ":"<<array_var[cont_tot_DifCov-1][i][j] << " ";
            //outdata_varSV<<endl;
                        outdata_var2SV<< cont_P << ":"<<array_var[cont_tot_DifCov-1][i][j]/cont_tot_DifCov << " ";
            //outdata_var2SV<<endl;

            cont_P++;

            }}
            outdata_meanSV<<endl;
            outdata_varSV<<endl;
            outdata_var2SV<<endl;

             ToEliminateMatrix2D(mean_Log, NFEATURES, NFEATURES);

///---------------------------------------------------
            ToEliminateMatrix3D(array_mean, cont_tot_DifCov, NFEATURES, NFEATURES);
            ToEliminateMatrix3D(array_var, cont_tot_DifCov, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(arr_info_mean, numTotalFrames_k, 2);
        }
        ToEliminateMatrix3D(arr_cov, numTotalFrames_k, NFEATURES, NFEATURES);
        ToEliminateMatrix2D(arr_info_first, numTotalFrames_k, 3);
    }


    ToEliminateMatrix2D(submat, NFEATURES, NFEATURES);
    ToEliminateMatrix1D(EigenValSym, NFEATURES);


            outdata_mean.close();
            outdata_var.close();
            outdata_var2.close();

            outdata_meanSV.close();
            outdata_varSV.close();
            outdata_var2SV.close();


}

