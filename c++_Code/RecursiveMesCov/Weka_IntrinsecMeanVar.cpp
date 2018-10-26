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

#include "Weka_IntrinsecMeanVar.h"
#include "RecursiveCov.h"


#include "../Resources/GeneralFunctions.h"
#include "../Resources/Matrixmanage.h"
#include "../TrajManage/MatricialTrajectoryManage.h"


#include "../Gradien_descent_Mean/Interface_NR3_Float.h"
#include "../Fast_Covariance/Fast_Covariance.h"




void Intrinsec_mean_Weka(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                        vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha)
{

    // The experiments were carried out with LED
//    ofstream outdata_meanW((name_file +"meanW.txt").c_str(), fstream::out); //1 means numero de divisiones del video
//    outdata_meanW.precision(10);
//    outdata_meanW<<  fixed;

    ofstream outdata_LogMeanW((name_file +"LogMeanW.arff").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogMeanW.precision(10);
    outdata_LogMeanW<<  fixed;

    ofstream outdata_LogVar((name_file +"LogVar.arff").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogVar.precision(10);
    outdata_LogVar<<  fixed;

    ofstream outdata_LogVarR((name_file +"LogVarR.arff").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogVarR.precision(10);
    outdata_LogVarR<<  fixed;


    ofstream outdata_LogMeanVar((name_file +"LogMeanVar.arff").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogMeanVar.precision(10);
    outdata_LogMeanVar<<  fixed;

    ofstream outdata_LogMeanVarR((name_file +"LogMeanVarR.arff").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogMeanVarR.precision(10);
    outdata_LogMeanVarR<<  fixed;

    /// Weka File Configuration


    outdata_LogMeanW    <<"%1. Title: intrinsec mean covariance features\n "<<endl;
    outdata_LogVar      <<"%1. Title: intrinsec var covariance features\n "<<endl;
    outdata_LogVarR     <<"%1. Title: intrinsec var-log covariance features\n "<<endl;
    outdata_LogMeanVar  <<"%1. Title: intrinsec mean + var covariance features\n "<<endl;
    outdata_LogMeanVarR <<"%1. Title: intrinsec mean + log(var) covariance features\n "<<endl;

    outdata_LogMeanW    <<"%2. author: fmartinezc \n % Date: 23  Mars\n\n "<<endl;
    outdata_LogVar      <<"%2. author: fmartinezc \n % Date: 23  Mars\n\n "<<endl;
    outdata_LogVarR     <<"%2. author: fmartinezc \n % Date: 23  Mars\n\n "<<endl;
    outdata_LogMeanVar  <<"%2. author: fmartinezc \n % Date: 23  Mars\n\n "<<endl;
    outdata_LogMeanVarR <<"%2. author: fmartinezc \n % Date: 23  Mars\n\n "<<endl;

    outdata_LogMeanW    <<"@relation IntrinsecMeanVar"<<endl;
    outdata_LogVar      <<"@relation IntrinsecMeanVar"<<endl;
    outdata_LogVarR     <<"@relation IntrinsecMeanVar"<<endl;
    outdata_LogMeanVar  <<"@relation IntrinsecMeanVar"<<endl;
    outdata_LogMeanVarR <<"@relation IntrinsecMeanVar"<<endl;


    int numCovFeat = NFEATURES*(NFEATURES + 1)/2;
    for(int nameF =0; nameF<numCovFeat; nameF++)
    {
    stringstream ss_numF; ss_numF << nameF;
    outdata_LogMeanW    << "@attribute F_" + ss_numF.str()  + " numeric" <<endl;
    outdata_LogVar      << "@attribute F_" + ss_numF.str()  + " numeric" <<endl;
    outdata_LogVarR     << "@attribute F_" + ss_numF.str()  + " numeric" <<endl;
    outdata_LogMeanVar  << "@attribute F_" + ss_numF.str()  + " numeric" <<endl;
    outdata_LogMeanVarR << "@attribute F_" + ss_numF.str()  + " numeric" <<endl;
    }

    outdata_LogMeanW    << "@attribute class {boxing handclapping handwaving jogging running walking}"<<endl;
    outdata_LogVar      << "@attribute class {boxing handclapping handwaving jogging running walking}"<<endl;
    outdata_LogVarR     << "@attribute class {boxing handclapping handwaving jogging running walking}"<<endl;
    outdata_LogMeanVar  << "@attribute class {boxing handclapping handwaving jogging running walking}"<<endl;
    outdata_LogMeanVarR << "@attribute class {boxing handclapping handwaving jogging running walking}"<<endl;


    outdata_LogMeanW    << "@data"<<endl;
    outdata_LogVar      << "@data"<<endl;
    outdata_LogVarR     << "@data"<<endl;
    outdata_LogMeanVar  << "@data"<<endl;
    outdata_LogMeanVarR << "@data"<<endl;


    float** submat = ToCreateMatrix2D(NFEATURES, NFEATURES);
    float* EigenValSym = ToCreateMatrix1D(NFEATURES);

    for(int k=0; k< vec_videos.size(); k++)
    {

        cout<< k<< " )  processing video: "<< vec_videos[k]  << endl;
        string path_textTrajectories_k = str_path_trajectories + vec_videos[k] + ".scale";
        //cin.ignore();


        int numTotalFrames_k, heightFrame_k, widthFrame_k, num_tot_traj_k;
        infoVideoFromTrajecAllScales(path_textTrajectories_k, numTotalFrames_k,
                                     heightFrame_k, widthFrame_k, num_tot_traj_k, 1); //3 pqrq ENSTA TRAJ

        int tot_seq_matrix_k = numTotalFrames_k+1; // mas dos porque en las dos ultimas posiciones se escribe el inicion y final de la trajectoria. Cabezera info.
       // int tot_feat_matrix_k = 12;//((int)featLabels.size()*(int)alphaScales.size())+3; 9 features + x,y,z

        cout<< "<----- creating mat: ["<< tot_seq_matrix_k<<"]["<< num_tot_traj_k<<"]["<< NFEATURES<<"]    ------>"<<endl;
        float*** sequenceKin_k = ToCreateMatrix3D(tot_seq_matrix_k, num_tot_traj_k, NFEATURES);
        ToInitMatrix3D(sequenceKin_k, tot_seq_matrix_k, num_tot_traj_k, NFEATURES);

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
                //cout<<" entro cov en frame "<<  frame_index <<endl;
                ToComputeCovmat(covmat_K, sequenceKin_k[frame_index],  NFEATURES, num_tot_traj_k, heightFrame_k, widthFrame_k);

//                cin.ignore();
//                ToPrintMatrix2D(covmat_K, NFEATURES, NFEATURES);

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
        ToEliminateMatrix3D(sequenceKin_k, tot_seq_matrix_k, num_tot_traj_k, NFEATURES);


        //cont_tot_DifCov++;
        //cout<< " cont_tot_DifCov: "<< cont_tot_DifCov << endl;


///-----------------------------------------------
        if(cont_tot_DifCov >num_mean)
        {


//            //cout<< " hizo lo de  "<< endl;
            float** _mean = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(_mean, NFEATURES, NFEATURES);

            float** Log_mean = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(Log_mean, NFEATURES, NFEATURES);

            float** Log_var2 = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(Log_var2, NFEATURES, NFEATURES);
            float** Log_var4 = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(Log_var4, NFEATURES, NFEATURES);


            float distAB =0.0;
            float_grad_mean_Array_Fletcher(arr_cov, _mean, cont_tot_DifCov, NFEATURES);
//            float_compute_grad_vectorMapping(arr_cov, VecMap_mean, cont_tot_DifCov, NFEATURES);

//            cin.ignore();
//            cout<< " covariance media original" << endl;
//            ToPrintMatrix2D(_mean, NFEATURES, NFEATURES);

//            cout<< " Vect map mean "<< endl;


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

///------------------------ To write weka file  -----------






//            ToPrintMatrix2D(Log_mean, NFEATURES, NFEATURES);

            ///------ outdata_meanCV ---
            //int cont_P =1;
            for(int i=0; i< NFEATURES; i++){
                for(int j=i; j< NFEATURES; j++){
                    outdata_LogMeanW<< Log_mean[i][j] << ",";
                    outdata_LogVar<< Log_var2[i][j] << ",";
                    outdata_LogVarR<< Log_var4[i][j]/cont_tot_DifCov << ",";
                    outdata_LogMeanVar<<Log_mean[i][j] << ",";
                    outdata_LogMeanVarR<<Log_mean[i][j] << ",";
                    //cont_P++;
                    }}


            for(int i=0; i< NFEATURES; i++){
                for(int j=i; j< NFEATURES; j++){
                outdata_LogMeanVar<<Log_var2[i][j] << ",";
                outdata_LogMeanVarR<<Log_var4[i][j]/cont_tot_DifCov << ",";
                //cont_P++;
                }}


           // cout<< " entro "<< endl;
             string  label_k = returnNameLabelAction_KTH(vec_videos[k], vec_activities );

            outdata_LogMeanW<< label_k<<endl;
            outdata_LogVar<< label_k<<endl;
            outdata_LogVarR<< label_k<<endl;
            outdata_LogMeanVar<< label_k<<endl;
            outdata_LogMeanVarR<< label_k<<endl;



            ToEliminateMatrix2D(_mean, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(Log_mean, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(Log_var2, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(Log_var4, NFEATURES, NFEATURES);


        }
        ToEliminateMatrix3D(arr_cov, numTotalFrames_k, NFEATURES, NFEATURES);
        ToEliminateMatrix2D(arr_info_first, numTotalFrames_k, 3);
    }


    ToEliminateMatrix2D(submat, NFEATURES, NFEATURES);
    ToEliminateMatrix1D(EigenValSym, NFEATURES);

    outdata_LogMeanW.close();
    //outdata_meanW.close();
    outdata_LogVar.close();
    outdata_LogVarR.close();

}

