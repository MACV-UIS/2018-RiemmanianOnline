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
#include "ConcatRecLogCov.h"


#include "../Resources/GeneralFunctions.h"
#include "../Resources/Matrixmanage.h"
#include "../TrajManage/MatricialTrajectoryManage.h"


#include "../Gradien_descent_Mean/Interface_NR3_Float.h"
#include "../Fast_Covariance/Fast_Covariance.h"

#include "../AppareanceTrajectoryFeatures/apperFeatTraj.h"



void FromAppareance_Feat(vector<string>  vec_videos, string name_file,  string str_path_trajectories, string path_frames,
                        vector<string> vec_activities, int NFEATURES, int num_mean, bool WRTT)
{

    // The experiments were carried out with LED
    // The experiments were carried out with LED
    ofstream outdata_meanW((name_file +"app_meanW.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_meanW.precision(10);
    outdata_meanW<<  fixed;

    ofstream outdata_LogMeanW((name_file +"app_LogMeanW.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogMeanW.precision(10);
    outdata_LogMeanW<<  fixed;

    ofstream outdata_LogVar((name_file +"app_LogVar.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogVar.precision(10);
    outdata_LogVar<<  fixed;

    ofstream outdata_LogVarR((name_file +"app_LogVarR.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogVarR.precision(10);
    outdata_LogVarR<<  fixed;


    ofstream outdata_LogMeanVar((name_file +"app_LogMeanVar.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogMeanVar.precision(10);
    outdata_LogMeanVar<<  fixed;

    ofstream outdata_LogMeanVarR((name_file +"app_LogMeanVarR.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogMeanVarR.precision(10);
    outdata_LogMeanVarR<<  fixed;



    float** submat = ToCreateMatrix2D(NFEATURES, NFEATURES);
    float* EigenValSym = ToCreateMatrix1D(NFEATURES);

    for(int k=0; k< vec_videos.size(); k++)
    {

        cout<< k<< " )  processing video: "<< vec_videos[k]  << endl;
        string path_textTrajectories_k = str_path_trajectories + vec_videos[k] + ".scale";
        string path_framesVideo_k = path_frames +  vec_videos[k] + "/";


        int numTotalFrames_k, heightFrame_k, widthFrame_k, num_tot_traj_k;
        infoVideoFromTrajecAllScales(path_textTrajectories_k, numTotalFrames_k,
                                     heightFrame_k, widthFrame_k, num_tot_traj_k, 1); //3 pqrq ENSTA TRAJ

        int tot_seq_matrix_k = numTotalFrames_k+1; // mas dos porque en las dos ultimas posiciones se escribe el inicion y final de la trajectoria. Cabezera info.
        int tot_feat_matrix_k = 8;//((int)featLabels.size()*(int)alphaScales.size())+3; 9 features + x,y,z

        cout<< "<----- creating mat of appareance: ["<< tot_seq_matrix_k<<"]["<< num_tot_traj_k<<"]["<< tot_feat_matrix_k<<"]    ------>"<<endl;

        float*** MatTraj = ToCreateMatrix3D(tot_seq_matrix_k, num_tot_traj_k, 3);
        ToInitMatrix3D(MatTraj, tot_seq_matrix_k, num_tot_traj_k, 3);

        float*** MatAppTraj = ToCreateMatrix3D(tot_seq_matrix_k, num_tot_traj_k, tot_feat_matrix_k);
        ToInitMatrix3D(MatAppTraj, tot_seq_matrix_k, num_tot_traj_k, tot_feat_matrix_k);


        ToLoadMatTrajOnly_XYZ(0, numTotalFrames_k,path_textTrajectories_k,
                                 heightFrame_k, widthFrame_k, num_tot_traj_k, 1, MatTraj);


        ToLoadMatTrajAppFeatures(0, numTotalFrames_k, path_textTrajectories_k, path_framesVideo_k,
                                 heightFrame_k, widthFrame_k, num_tot_traj_k, 1, MatTraj, MatAppTraj, WRTT); //true is w.r.t trajectories


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
            int traj_used = num_traj_valides(MatTraj[frame_index],  num_tot_traj_k);
            //cout << " numero de frame: "<< frame_index << " numTotalFrames_k: "<< numTotalFrames_k << " trajectories used "<< traj_used << endl;
            if(traj_used>100)
            {
                float **covmat_K = ToCreateMatrix2D(NFEATURES, NFEATURES); //Covariance Matrix buffer
                ToInitMatrix2D(covmat_K, NFEATURES, NFEATURES);
                //cout<<" entro cov"<< endl;
                ToComputeCovmat(covmat_K, MatAppTraj[frame_index],  NFEATURES, num_tot_traj_k, heightFrame_k, widthFrame_k);

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

        ToEliminateMatrix3D(MatTraj, tot_seq_matrix_k, num_tot_traj_k, 3);
        ToEliminateMatrix3D(MatAppTraj, tot_seq_matrix_k, num_tot_traj_k, tot_feat_matrix_k);




        //cont_tot_DifCov++;
        cout<< " cont_tot_DifCov: "<< cont_tot_DifCov << endl;


///-----------------------------------------------
///-----------------------------------------------
        if(cont_tot_DifCov >num_mean)
        {

            //cout<< " entro "<< endl;
            int label_k = atoi(returnLabelAction_KTH(vec_videos[k], vec_activities ).c_str()); //-1 para que la primera sea cero

            outdata_LogMeanW<< label_k<< " ";
            outdata_meanW<< label_k<< " ";
            outdata_LogVar<< label_k<< " ";
            outdata_LogVarR<< label_k<< " ";
            outdata_LogMeanVar<< label_k<< " ";
            outdata_LogMeanVarR<< label_k<< " ";

//            //cout<< " hizo lo de  "<< endl;
            float** _mean = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(_mean, NFEATURES, NFEATURES);

            float** VecMap_mean = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(VecMap_mean, NFEATURES, NFEATURES);

            float** Log_mean = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(Log_mean, NFEATURES, NFEATURES);

            float** Log_var2 = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(Log_var2, NFEATURES, NFEATURES);
            float** Log_var4 = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(Log_var4, NFEATURES, NFEATURES);


            float distAB =0.0;
            float_grad_mean_Array_Fletcher(arr_cov, _mean, cont_tot_DifCov, NFEATURES);
            float_compute_grad_vectorMapping(arr_cov, VecMap_mean, cont_tot_DifCov, NFEATURES);

//            cin.ignore();
//            cout<< " covariance media original" << endl;
//            ToPrintMatrix2D(_mean, NFEATURES, NFEATURES);

//            cout<< " Vect map mean "<< endl;

            for(int i=0; i< NFEATURES; i++)
                for(int j=i+1; j< NFEATURES; j++)
                    VecMap_mean[i][j] = 1.4142*VecMap_mean[i][j];

//            ToPrintMatrix2D(VecMap_mean, NFEATURES, NFEATURES);
//            cout<< " Log mean Weighted"<< endl;

            float_Log_A_resp_I(_mean, Log_mean,  NFEATURES, distAB);

///----------------------------------------------------------------------------
                for(int i_cov = 0; i_cov < cont_tot_DifCov; i_cov++)
                {
                    float** Log_temp_v2 = ToCreateMatrix2D(NFEATURES, NFEATURES);
                    ToInitMatrix2D(Log_temp_v2, NFEATURES, NFEATURES);

                    float** Log_temp_v4 = ToCreateMatrix2D(NFEATURES, NFEATURES);
                    ToInitMatrix2D(Log_temp_v4, NFEATURES, NFEATURES);

                    float_Log_A_resp_B_Variance(arr_cov[i_cov], _mean, Log_temp_v2,  NFEATURES,  distAB);
                    float_Log_A_resp_I(arr_cov[i_cov], Log_temp_v4,  NFEATURES, distAB);

                    for(int m=0; m<NFEATURES; m++)
                    {
                        for(int n=0; n<NFEATURES; n++)
                        {
                            Log_var2[m][n] =Log_var2[m][n] + Log_temp_v2[m][n];
                            Log_var4[m][n] = Log_var4[m][n] + pow( Log_temp_v4[m][n] -Log_mean[m][n], 2 );
                        }
                    }

                    ToEliminateMatrix2D(Log_temp_v2, NFEATURES, NFEATURES);
                    ToEliminateMatrix2D(Log_temp_v4, NFEATURES, NFEATURES);

                }

///----------------------------------------------------------------------------
//            ToPrintMatrix2D(Log_mean, NFEATURES, NFEATURES);

            ///------ outdata_meanCV ---
            int cont_P =1;
            for(int i=0; i< NFEATURES; i++){
                for(int j=i; j< NFEATURES; j++){
                    outdata_LogMeanW<< cont_P << ":"<<VecMap_mean[i][j] << " ";
                    outdata_meanW<< cont_P << ":"<<Log_mean[i][j] << " ";
                    outdata_LogVar<< cont_P << ":"<<Log_var2[i][j] << " ";
                    outdata_LogVarR<< cont_P << ":"<<Log_var4[i][j]/cont_tot_DifCov << " ";
                    outdata_LogMeanVar<< cont_P << ":"<<Log_mean[i][j] << " ";
                    outdata_LogMeanVarR<< cont_P << ":"<<Log_mean[i][j] << " ";
                    cont_P++;
                    }}
            for(int i=0; i< NFEATURES; i++){
                for(int j=i; j< NFEATURES; j++){
                outdata_LogMeanVar<< cont_P << ":"<<Log_var2[i][j] << " ";
                outdata_LogMeanVarR<< cont_P << ":"<<Log_var4[i][j]/cont_tot_DifCov << " ";
                cont_P++;
                }}

            outdata_LogMeanW<<endl;
            outdata_meanW<<endl;
            outdata_LogVar<<endl;
            outdata_LogVarR<<endl;

            outdata_LogMeanVar<<endl;
            outdata_LogMeanVarR<<endl;


            ToEliminateMatrix2D(_mean, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(Log_mean, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(VecMap_mean, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(Log_var2, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(Log_var4, NFEATURES, NFEATURES);


        }
        ToEliminateMatrix3D(arr_cov, numTotalFrames_k, NFEATURES, NFEATURES);
        ToEliminateMatrix2D(arr_info_first, numTotalFrames_k, 3);
    }


    ToEliminateMatrix2D(submat, NFEATURES, NFEATURES);
    ToEliminateMatrix1D(EigenValSym, NFEATURES);

    outdata_LogMeanW.close();
    outdata_meanW.close();
    outdata_LogVar.close();
    outdata_LogVarR.close();


}
