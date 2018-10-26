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








void Concat_Rec_Log_CoV(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                    vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha)
{

    // The experiments were carried out with LED
    ofstream outdata_mean_var((name_file + "mean_var.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_mean_var.precision(10);
    outdata_mean_var<<  fixed;

    ofstream outdata_mean_min((name_file + "mean_min.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_mean_min.precision(10);
    outdata_mean_min<<  fixed;

    ofstream outdata_mean_max((name_file + "mean_max.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_mean_max.precision(10);
    outdata_mean_max<<  fixed;

    ofstream outdata_mean_der((name_file + "mean_der.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_mean_der.precision(10);
    outdata_mean_der<<  fixed;

    ofstream outdata_conc_all((name_file + "all.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_conc_all.precision(10);
    outdata_conc_all<<  fixed;



    float** submat = ToCreateMatrix2D(NFEATURES, NFEATURES);
    float* EigenValSym = ToCreateMatrix1D(NFEATURES);

    for(int k=0; k< vec_videos.size(); k++)
    {



        cout<< k<< " )  processing video: "<< vec_videos[k]  << endl;
        string path_textTrajectories_k = str_path_trajectories + vec_videos[k] + ".scale";



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

            outdata_mean_var<< label_k<< " ";
            outdata_mean_min<< label_k<< " ";
            outdata_mean_max<< label_k<< " ";
            outdata_mean_der<< label_k<< " ";
            outdata_conc_all<< label_k<< " ";





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
                        outdata_mean_var<< cont_P << ":"<<array_mean[cont_tot_DifCov-1][i][j] << " ";
                        outdata_mean_min<< cont_P << ":"<<array_mean[cont_tot_DifCov-1][i][j] << " ";
                        outdata_mean_max<< cont_P << ":"<<array_mean[cont_tot_DifCov-1][i][j] << " ";
                        outdata_mean_der<< cont_P << ":"<<array_mean[cont_tot_DifCov-1][i][j] << " ";
                        outdata_conc_all<< cont_P << ":"<<array_mean[cont_tot_DifCov-1][i][j] << " ";
                        cont_P++;
                        }}


            for(int i=0; i< NFEATURES; i++){
                for(int j=i; j< NFEATURES; j++){
                        outdata_mean_var<< cont_P << ":"<<array_var[cont_tot_DifCov-1][i][j] << " ";
                        outdata_mean_min<< cont_P << ":"<<array_min[cont_tot_DifCov-1][i][j] << " ";
                        outdata_mean_max<< cont_P << ":"<<array_max[cont_tot_DifCov-1][i][j] << " ";
                        outdata_mean_der<< cont_P << ":"<<array_dif[cont_tot_DifCov-1][i][j] << " ";
                        outdata_conc_all<< cont_P << ":"<<array_var[cont_tot_DifCov-1][i][j] << " ";
                        cont_P++;
                        }}

             for(int i=0; i< NFEATURES; i++)
                for(int j=i; j< NFEATURES; j++)
                    outdata_conc_all<< cont_P++ << ":"<<array_min[cont_tot_DifCov-1][i][j] << " ";

             for(int i=0; i< NFEATURES; i++)
                for(int j=i; j< NFEATURES; j++)
                    outdata_conc_all<< cont_P++ << ":"<<array_max[cont_tot_DifCov-1][i][j] << " ";

             for(int i=0; i< NFEATURES; i++)
                for(int j=i; j< NFEATURES; j++)
                    outdata_conc_all<< cont_P++ << ":"<<array_dif[cont_tot_DifCov-1][i][j] << " ";



            outdata_mean_var<<endl;
            outdata_mean_min<<endl;
            outdata_mean_max<<endl;
            outdata_mean_der<<endl;
            outdata_conc_all<<endl;


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



            outdata_mean_var.close();
            outdata_mean_min.close();
            outdata_mean_max.close();
            outdata_mean_der.close();
            outdata_conc_all.close();



}




void Concat_Scal_Rec_Log_CoV(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                    vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha, int num_scales)
{

    // The experiments were carried out with LED

    ofstream outdata_mean((name_file + "mean.txt").c_str(), fstream::out); //1 means numero de divisiones del video
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
        string path_textTrajectories_k = str_path_trajectories + vec_videos[k] + ".scale";



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


            int cont_P =1;
            for(int iter_alpha=0; iter_alpha< num_scales; iter_alpha++){
//            //cout<< " hizo lo de  "<< endl;
            float alpha_val = pow(2, -1*(init_alpha+iter_alpha));
            cout<< "alpha val: "<<alpha_val << endl;

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

            }






            for(int i=0; i< NFEATURES; i++){
                for(int j=i; j< NFEATURES; j++){
                        outdata_mean<< cont_P << ":"<<array_mean[cont_tot_DifCov-1][i][j]<< " ";
                        outdata_var<< cont_P << ":"<<array_var[cont_tot_DifCov-1][i][j] << " ";
                        outdata_min<< cont_P << ":"<<array_min[cont_tot_DifCov-1][i][j] << " ";
                        outdata_max<< cont_P << ":"<<array_max[cont_tot_DifCov-1][i][j] << " ";
                        outdata_der<< cont_P << ":"<<array_dif[cont_tot_DifCov-1][i][j] << " ";
                        cont_P++;
                        }}



///---------------------------------------------------
            ToEliminateMatrix3D(array_mean, cont_tot_DifCov, NFEATURES, NFEATURES);
            ToEliminateMatrix3D(array_var, cont_tot_DifCov, NFEATURES, NFEATURES);
            ToEliminateMatrix3D(array_min, cont_tot_DifCov, NFEATURES, NFEATURES);
            ToEliminateMatrix3D(array_max, cont_tot_DifCov, NFEATURES, NFEATURES);
            ToEliminateMatrix3D(array_dif, cont_tot_DifCov, NFEATURES, NFEATURES);
            }

            outdata_mean<<endl;
            outdata_var<<endl;
            outdata_min<<endl;
            outdata_max<<endl;
            outdata_der<<endl;


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

