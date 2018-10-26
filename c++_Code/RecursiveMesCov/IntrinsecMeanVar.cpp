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

#include "IntrinsecMeanVar.h"
#include "RecursiveCov.h"


#include "../Resources/GeneralFunctions.h"
#include "../Resources/Matrixmanage.h"
#include "../TrajManage/MatricialTrajectoryManage.h"


#include "../Gradien_descent_Mean/Interface_NR3_Float.h"
#include "../Fast_Covariance/Fast_Covariance.h"



void Intrinsec_Mean_Var_tDiv_Iterative(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                                       vector<string> vec_activities, int NFEATURES, int num_mean, int t_Div, int max_iter)
{

    // The experiments were carried out with LED
    ofstream outdata_meanCV((name_file +"_meanCV.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_meanCV.precision(10);
    outdata_meanCV<<  fixed;

    ofstream outdata_meanSV((name_file +"_meanSV.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_meanSV.precision(10);
    outdata_meanSV<<  fixed;

    ofstream outdata_var2CV((name_file + "_var2CV.txt").c_str(), fstream::out); //suma de logup cuadrado; peor el cuadrado es en la diagonal
    outdata_var2CV.precision(10);
    outdata_var2CV<<  fixed;

    ofstream outdata_var2SV((name_file + "_var2SV.txt").c_str(), fstream::out); //suma de logup cuadrado; peor el cuadrado es en la diagonal
    outdata_var2SV.precision(10);
    outdata_var2SV<<  fixed;

    ofstream outdata_var3CV((name_file + "_var3CV.txt").c_str(), fstream::out); // como var2 pero normalizando entre N
    outdata_var3CV.precision(10);
    outdata_var3CV<<  fixed;

    ofstream outdata_var3SV((name_file + "_var3SV.txt").c_str(), fstream::out); // como var2 pero normalizando entre N
    outdata_var3SV.precision(10);
    outdata_var3SV<<  fixed;

    ofstream outdata_var4CV((name_file + "_var4CV.txt").c_str(), fstream::out); // pasar mean and var a log space identity y sumar diferencias cudradas
    outdata_var4CV.precision(10);
    outdata_var4CV<<  fixed;

    ofstream outdata_var4SV((name_file + "_var4SV.txt").c_str(), fstream::out); // pasar mean and var a log space identity y sumar diferencias cudradas
    outdata_var4SV.precision(10);
    outdata_var4SV<<  fixed;


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
        int cont_frameCov =0;

        int cont_tot_DifCov =0;
        float *** arr_cov  = ToCreateMatrix3D(numTotalFrames_k, NFEATURES, NFEATURES);
        ToInitMatrix3D(arr_cov,numTotalFrames_k, NFEATURES, NFEATURES );


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
                        for( int m=0; m<NFEATURES; m++)
                        {
                            for(int n=0; n<NFEATURES; n++)
                            {
                                arr_cov[cont_tot_DifCov][m][n] = covmat_K[m][n];

                            }
                        }

                        cont_tot_DifCov++;
                    cont_frameCov++;

                }


                ToEliminateMatrix2D(covmat_K, NFEATURES, NFEATURES); // no la pasa porque es del vector
            }

        }

        ToEliminateMatrix3D(sequenceKin_k, tot_seq_matrix_k, num_tot_traj_k, tot_feat_matrix_k);
        //cont_tot_DifCov++;
        cout<< " cont_tot_DifCov: "<< cont_tot_DifCov << endl;
/////-----------------------------------------------
        if(cont_tot_DifCov >num_mean)
        {
            int label_k = atoi(returnLabelAction_KTH(vec_videos[k], vec_activities ).c_str()); //-1 para que la primera sea cero

            outdata_meanCV<< label_k<< " ";
            outdata_meanSV<< label_k<< " ";
            outdata_var2CV<< label_k<< " ";
            outdata_var2SV<< label_k<< " ";
            outdata_var3CV<< label_k<< " ";
            outdata_var3SV<< label_k<< " ";
            outdata_var4CV<< label_k<< " ";
            outdata_var4SV<< label_k<< " ";

            int div_fix=(cont_tot_DifCov/t_Div);
            cout<< " divisiones que corresponden: "<< div_fix << endl;
            int cont_P =1;

            float *** vec_means = ToCreateMatrix3D(t_Div,NFEATURES, NFEATURES);
            ToInitMatrix3D(vec_means, t_Div,NFEATURES, NFEATURES);

            float *** Log_vec_means = ToCreateMatrix3D(t_Div,NFEATURES, NFEATURES);
            ToInitMatrix3D(Log_vec_means, t_Div,NFEATURES, NFEATURES);

            float *** Log_vec_var2 = ToCreateMatrix3D(t_Div,NFEATURES, NFEATURES);
            ToInitMatrix3D(Log_vec_var2, t_Div,NFEATURES, NFEATURES);

            float *** Log_vec_var4 = ToCreateMatrix3D(t_Div,NFEATURES, NFEATURES);
            ToInitMatrix3D(Log_vec_var4, t_Div,NFEATURES, NFEATURES);


            cout<< " indices de inicializacion: "<< endl;
            for( int cont_fix =0; cont_fix< t_Div; cont_fix++ )
            {
                cout<< (int)((cont_tot_DifCov/(t_Div+1))*(cont_fix+1))<< " ";
                Copy_Cov(arr_cov[(int)((cont_tot_DifCov/(t_Div+1))*(cont_fix+1))], vec_means[cont_fix], NFEATURES);

            }
            cout<< endl;

            int cont_iter =0;
            float** arra_cont = ToCreateMatrix2D(t_Div,2);
            float**** array_toComp_means = ToCreateMatrix4D(t_Div, cont_tot_DifCov, NFEATURES, NFEATURES);
            while(cont_iter < max_iter )
            {
                cout<<" k-means iteraction: "<<cont_iter << endl;
                ToInitMatrix4D(array_toComp_means, t_Div, cont_tot_DifCov, NFEATURES, NFEATURES);
                ToInitMatrix2D(arra_cont, t_Div, 2); // array con "contadores dinamicos"
                for(int l=0; l<t_Div; l++)
                    arra_cont[l][1]=l;



                for(int i_cov =0; i_cov<cont_tot_DifCov; i_cov++ )
                {
                    ///2) compare each cov with the four k and assign to the closer
                    float dist_min = 10000;
                    int k_index = 0;
                    for(int i_mean =0; i_mean < t_Div; i_mean++)
                    {
                        float dist_temp = float_AffineInv_Dist(arr_cov[i_cov], vec_means[i_mean],  NFEATURES);
                        //cout<< "i_mean: "<< i_mean << " dist_temp "<< dist_temp;
                        if(dist_temp< dist_min)
                        {
                            dist_min = dist_temp;
                            k_index=i_mean;
                        }
                    }
                    ///----prueba de distancia seleccionado
//                    cout<< " prueba indice de media seleccionada: "<< k_index << endl;
//                    cin.ignore();
                    ///3) Assign the cov to the close mean
                    Copy_Cov( arr_cov[i_cov], array_toComp_means[k_index][(int)arra_cont[k_index][0]++], NFEATURES);
                    /// PROBAR SI EL CONTADOR SI ESTA CRECIENDO
//                    cin.ignore();
//                    cout<< " test  the coun updating "<< endl;
//                    cout<< " selected k "<< k << endl;
//                    ToPrintMatrix1D(arra_cont, t_Div);
                }
                ///4) compute the mean of the cov --------------------------------------------

                float distAB =0.0;
//                 cout<< " selected k's " << endl;
//                 ToPrintMatrix1D(arra_cont, t_Div);
//                 cout<< endl;

                for( int cont_means =0; cont_means< t_Div; cont_means++ )
                {
                    if( arra_cont[cont_means][0] > 1){
                    float_grad_mean_Array_Fletcher( array_toComp_means[cont_means] , vec_means[cont_means] , arra_cont[cont_means][0], NFEATURES);
                    float_Log_A_resp_I(vec_means[cont_means], Log_vec_means[cont_means],  NFEATURES, distAB);
                    ///------- means ---
//                    cin.ignore();
//                    cout<< cont_means <<") prueba de mean  "<< endl;
//                    ToPrintMatrix2D(vec_means[cont_means], NFEATURES, NFEATURES);

                    for(int j=0; j< (int)arra_cont[cont_means][0]; j++)
                    {
                        //------------------------------------------------------------------
                        float** Log_temp_v2 = ToCreateMatrix2D(NFEATURES, NFEATURES);
                        ToInitMatrix2D(Log_temp_v2, NFEATURES, NFEATURES);

                        float** Log_temp_v4 = ToCreateMatrix2D(NFEATURES, NFEATURES);
                        ToInitMatrix2D(Log_temp_v4, NFEATURES, NFEATURES);

                        float_Log_A_resp_B_Variance(array_toComp_means[cont_means][j], vec_means[cont_means], Log_temp_v2,  NFEATURES,  distAB);
                        float_Log_A_resp_I(array_toComp_means[cont_means][j], Log_temp_v4,  NFEATURES, distAB);

                        for(int m=0; m<NFEATURES; m++)
                       {
                           for(int n=0; n<NFEATURES; n++)
                            {
                                Log_vec_var2[cont_means][m][n] = Log_vec_var2[cont_means][m][n] + Log_temp_v2[m][n];
                                Log_vec_var4[cont_means][m][n] = Log_vec_var4[cont_means][m][n] + pow( Log_temp_v4[m][n] - Log_vec_means[cont_means][m][n], 2 );
                            }
                        }
//                        /// otra prueba de construccion de cov
//                        cin.ignore();
//                        cout<< " print cov temporal "<< endl;
//                        ToPrintMatrix2D(Log_vec_var2[cont_means], NFEATURES, NFEATURES);

                        ToEliminateMatrix2D(Log_temp_v2, NFEATURES, NFEATURES);
                        ToEliminateMatrix2D(Log_temp_v4, NFEATURES, NFEATURES);
                        //------------------------------------------------------------------
                    }
                    }else{
                        Copy_Cov(arr_cov[(int)(cont_tot_DifCov/2)], vec_means[cont_means], NFEATURES);
                    }

                }

                cont_iter++;
            }

            ///-------------- Ordenando las matrices ------

             quickSort_cov(arra_cont,  0, t_Div-1);

        ///--------------------------copiado de matrices en orden ------

            float*** Log_means_sort = ToCreateMatrix3D(t_Div,NFEATURES, NFEATURES);
            float*** Log_var2_sort  = ToCreateMatrix3D(t_Div,NFEATURES, NFEATURES);
            float*** Log_var4_sort  = ToCreateMatrix3D(t_Div,NFEATURES, NFEATURES);

            ToInitMatrix3D(Log_means_sort, t_Div,NFEATURES, NFEATURES);
            ToInitMatrix3D(Log_var2_sort, t_Div,NFEATURES, NFEATURES);
            ToInitMatrix3D(Log_var4_sort, t_Div,NFEATURES, NFEATURES);


            for(int a =0; a< t_Div; a++){
                Copy_Cov(Log_vec_means[(int)arra_cont[a][1]], Log_means_sort[a], NFEATURES);
                Copy_Cov(Log_vec_var2[(int)arra_cont[a][1]], Log_var2_sort[a], NFEATURES);
                Copy_Cov(Log_vec_var4[(int)arra_cont[a][1]], Log_var4_sort[a], NFEATURES);
            }

            ToEliminateMatrix3D(Log_vec_means, t_Div,NFEATURES, NFEATURES);
            ToEliminateMatrix3D(Log_vec_var2, t_Div,NFEATURES, NFEATURES);
            ToEliminateMatrix3D(Log_vec_var4, t_Div,NFEATURES, NFEATURES);

            ///------------------------------------------------------------------------
            ///---------------------- outdata_meanCV ---------------------------

            for(int a=0; a < t_Div; a++)
            {
                for(int i=0; i< NFEATURES; i++)
                {
                    for(int j=i; j< NFEATURES; j++)
                    {
                        outdata_meanCV<< cont_P << ":"<<Log_means_sort[a][i][j] << " ";
                        outdata_var2CV<< cont_P << ":"<<Log_var2_sort[a][i][j] << " ";
                        outdata_var3CV<< cont_P << ":"<<Log_var2_sort[a][i][j]/arra_cont[a][0] << " ";
                        outdata_var4CV<< cont_P << ":"<<Log_var4_sort[a][i][j]/arra_cont[a][0] << " ";
                        cont_P++;
                    }
                }


            }

        cont_P=1;
        for(int a=0; a < t_Div; a++)
            {
                for(int i=0; i< NFEATURES; i++)
                {
                    for(int j=i+1; j< NFEATURES; j++)
                    {
                        outdata_meanSV<< cont_P << ":"<<Log_means_sort[a][i][j] << " ";
                        outdata_var2SV<< cont_P << ":"<<Log_var2_sort[a][i][j] << " ";
                        outdata_var3SV<< cont_P << ":"<<Log_var2_sort[a][i][j]/arra_cont[a][0] << " ";
                        outdata_var4SV<< cont_P << ":"<<Log_var4_sort[a][i][j]/arra_cont[a][0] << " ";
                        cont_P++;
                    }
                }

            }

        outdata_meanCV<< endl;
        outdata_meanSV<< endl;
        outdata_var2CV<< endl;
        outdata_var2SV<< endl;
        outdata_var3CV<< endl;
        outdata_var3SV<< endl;
        outdata_var4CV<< endl;
        outdata_var4SV<< endl;

        ToEliminateMatrix3D(vec_means, t_Div,NFEATURES, NFEATURES);
        ToEliminateMatrix3D(Log_means_sort, t_Div,NFEATURES, NFEATURES);
        ToEliminateMatrix3D(Log_var2_sort, t_Div,NFEATURES, NFEATURES);
        ToEliminateMatrix3D(Log_var4_sort, t_Div,NFEATURES, NFEATURES);

        ToEliminateMatrix4D(array_toComp_means, t_Div, cont_tot_DifCov, NFEATURES, NFEATURES);
        ToEliminateMatrix2D(arra_cont, t_Div,2);

}
    ToEliminateMatrix3D(arr_cov, numTotalFrames_k, NFEATURES, NFEATURES);
}


ToEliminateMatrix2D(submat, NFEATURES, NFEATURES);
ToEliminateMatrix1D(EigenValSym, NFEATURES);

outdata_meanCV.close();
outdata_meanSV.close();
outdata_var2CV.close();
outdata_var2SV.close();
outdata_var3CV.close();
outdata_var3SV.close();
outdata_var4CV.close();
outdata_var4SV.close();

}


void Intercambiar_cov(float** K_sort, int i, int j)
{
    float tmp_dist;
    float tmp_index;

    tmp_dist   = K_sort[i][0];
    tmp_index   = K_sort[i][1];

    K_sort[i][0]   = K_sort[j][0];
    K_sort[i][1]   = K_sort[j][1];


    K_sort[j][0] = tmp_dist;
    K_sort[j][1] = tmp_index;

}



void quickSort_cov(float**  K_sort,  int left, int right)
{

    int i = left, j = right;
    float pivot = K_sort[(left + right) / 2][0];
    //cout<<" pivot "<< pivot << endl;

    /** Partición */
    while (i <= j)
    {
        while (K_sort[i][0] > pivot)i++;
        while (K_sort[j][0] < pivot)j--;
        if (i <= j)
        {
            Intercambiar_cov(K_sort, i, j);
            i++;
            j--;
        }
    }

    /** Recursividad */
    if (left < j)
    {
        quickSort_cov(K_sort, left, j);
    }
    if (i < right)
    {
        quickSort_cov(K_sort, i, right);
    }



}




void Intrinsec_Mean_Var_tDiv(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                             vector<string> vec_activities, int NFEATURES, int num_mean, int t_Div)
{

    // The experiments were carried out with LED
    ofstream outdata_meanCV((name_file +"_tDiv_meanCV.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_meanCV.precision(10);
    outdata_meanCV<<  fixed;

    ofstream outdata_meanSV((name_file +"_tDiv_meanSV.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_meanSV.precision(10);
    outdata_meanSV<<  fixed;


    ofstream outdata_var2CV((name_file + "_tDiv_var2CV.txt").c_str(), fstream::out); //suma de logup cuadrado; peor el cuadrado es en la diagonal
    outdata_var2CV.precision(10);
    outdata_var2CV<<  fixed;

    ofstream outdata_var2SV((name_file + "_tDiv_var2SV.txt").c_str(), fstream::out); //suma de logup cuadrado; peor el cuadrado es en la diagonal
    outdata_var2SV.precision(10);
    outdata_var2SV<<  fixed;

    ofstream outdata_var3CV((name_file + "_tDiv_var3CV.txt").c_str(), fstream::out); // como var2 pero normalizando entre N
    outdata_var3CV.precision(10);
    outdata_var3CV<<  fixed;

    ofstream outdata_var3SV((name_file + "_tDiv_var3SV.txt").c_str(), fstream::out); // como var2 pero normalizando entre N
    outdata_var3SV.precision(10);
    outdata_var3SV<<  fixed;

    ofstream outdata_var4CV((name_file + "_tDiv_var4CV.txt").c_str(), fstream::out); // pasar mean and var a log space y sumar diferencias cudradas
    outdata_var4CV.precision(10);
    outdata_var4CV<<  fixed;

    ofstream outdata_var4SV((name_file + "_tDiv_var4SV.txt").c_str(), fstream::out); // pasar mean and var a log space y sumar diferencias cudradas
    outdata_var4SV.precision(10);
    outdata_var4SV<<  fixed;



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

            outdata_meanCV<< label_k<< " ";
            outdata_meanSV<< label_k<< " ";
            outdata_var2CV<< label_k<< " ";
            outdata_var2SV<< label_k<< " ";
            outdata_var3CV<< label_k<< " ";
            outdata_var3SV<< label_k<< " ";
            outdata_var4CV<< label_k<< " ";
            outdata_var4SV<< label_k<< " ";




            int div_fix=(cont_tot_DifCov/t_Div);
            cout<< " divisiones que corresponden: "<< div_fix << endl;
            int cont_P1 =1, cont_P2 =1, cont_P3 =1, cont_P4 =1, cont_P5 =1, cont_P6 =1, cont_P7 =1;
            int cont_P8 =1, cont_P9 =1, cont_P10 =1, cont_P11 =1, cont_P12 =1, cont_P13 =1, cont_P14 =1;


            for( int cont_fix =0; cont_fix< t_Div; cont_fix++ )
            {
                cout<< " entro a dividir"<< endl;
                float*** temp_split_vec = ToCreateMatrix3D(div_fix, NFEATURES, NFEATURES);
                ToInitMatrix3D(temp_split_vec, div_fix, NFEATURES, NFEATURES);


                int index_temp =0;
                for(int n_div= div_fix*cont_fix; n_div < div_fix*(cont_fix+1); n_div++ )
                {
                    Copy_Cov(arr_cov[n_div], temp_split_vec[index_temp],NFEATURES);
//                cout<< " cov copiada"<< endl;
//                cin.ignore();
//                ToPrintMatrix2D();
                    index_temp++;
                }

//            //cout<< " hizo lo de  "<< endl;
                float** _mean = ToCreateMatrix2D(NFEATURES, NFEATURES);
                ToInitMatrix2D(_mean, NFEATURES, NFEATURES);
                float** Log_mean = ToCreateMatrix2D(NFEATURES, NFEATURES);
                ToInitMatrix2D(Log_mean, NFEATURES, NFEATURES);

                float** _var2 = ToCreateMatrix2D(NFEATURES, NFEATURES);
                ToInitMatrix2D(_var2, NFEATURES, NFEATURES);
                float** Log_var2 = ToCreateMatrix2D(NFEATURES, NFEATURES);
                ToInitMatrix2D(Log_var2, NFEATURES, NFEATURES);

                float** _var4 = ToCreateMatrix2D(NFEATURES, NFEATURES);
                ToInitMatrix2D(_var4, NFEATURES, NFEATURES);
                float** Log_var4 = ToCreateMatrix2D(NFEATURES, NFEATURES);
                ToInitMatrix2D(Log_var4, NFEATURES, NFEATURES);


                float distAB =0.0;
                float_grad_mean_Array_Fletcher(temp_split_vec, _mean, div_fix, NFEATURES);
                float_Log_A_resp_I(_mean, Log_mean,  NFEATURES, distAB);

                for(int i_cov = 0; i_cov < div_fix; i_cov++)
                {
                    float** Log_temp_v2 = ToCreateMatrix2D(NFEATURES, NFEATURES);
                    ToInitMatrix2D(Log_temp_v2, NFEATURES, NFEATURES);

                    float** Log_temp_v4 = ToCreateMatrix2D(NFEATURES, NFEATURES);
                    ToInitMatrix2D(Log_temp_v4, NFEATURES, NFEATURES);

                    float_Log_A_resp_B_Variance(temp_split_vec[i_cov], _mean, Log_temp_v2,  NFEATURES,  distAB);
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


                ///------ outdata_meanCV ---

                for(int i=0; i< NFEATURES; i++)
                    for(int j=i; j< NFEATURES; j++)
                        outdata_meanCV<< cont_P1++ << ":"<<Log_mean[i][j] << " ";

                for(int i=0; i< NFEATURES; i++)
                    for(int j=i+1; j< NFEATURES; j++)
                        outdata_meanSV<< cont_P2++ << ":"<<Log_mean[i][j] << " ";


                ///---var2
                for(int i=0; i< NFEATURES; i++)
                    for(int j=i; j< NFEATURES; j++)
                        outdata_var2CV<< cont_P5++ << ":"<<Log_var2[i][j] << " ";


                for(int i=0; i< NFEATURES; i++)
                    for(int j=i+1; j< NFEATURES; j++)
                        outdata_var2SV<< cont_P6++ << ":"<<Log_var2[i][j] << " ";

                ///---var3

                for(int i=0; i< NFEATURES; i++)
                    for(int j=i; j< NFEATURES; j++)
                        outdata_var3CV<< cont_P7++ << ":"<<Log_var2[i][j]/cont_tot_DifCov << " ";


                for(int i=0; i< NFEATURES; i++)
                    for(int j=i+1; j< NFEATURES; j++)
                        outdata_var3SV<< cont_P8++ << ":"<<Log_var2[i][j]/cont_tot_DifCov << " ";

                ///---var4
                for(int i=0; i< NFEATURES; i++)
                    for(int j=i; j< NFEATURES; j++)
                        outdata_var4CV<< cont_P9++ << ":"<<Log_var4[i][j]/cont_tot_DifCov << " ";


                for(int i=0; i< NFEATURES; i++)
                    for(int j=i+1; j< NFEATURES; j++)
                        outdata_var4SV<< cont_P10++ << ":"<<Log_var4[i][j]/div_fix << " ";






                ToEliminateMatrix2D(_mean, NFEATURES, NFEATURES);
                ToEliminateMatrix2D(Log_mean, NFEATURES, NFEATURES);
                ToEliminateMatrix2D(_var2, NFEATURES, NFEATURES);
                ToEliminateMatrix2D(Log_var2, NFEATURES, NFEATURES);
                ToEliminateMatrix2D(_var4, NFEATURES, NFEATURES);
                ToEliminateMatrix2D(Log_var4, NFEATURES, NFEATURES);



///-----------------------------------------------------------------------------------
                ToEliminateMatrix3D(temp_split_vec, div_fix, NFEATURES, NFEATURES);
            }


            outdata_meanCV<< endl;
            outdata_meanSV<< endl;
            outdata_var2CV<< endl;
            outdata_var2SV<< endl;
            outdata_var3CV<< endl;
            outdata_var3SV<< endl;
            outdata_var4CV<< endl;
            outdata_var4SV<< endl;


        }
        ToEliminateMatrix3D(arr_cov, numTotalFrames_k, NFEATURES, NFEATURES);
        ToEliminateMatrix2D(arr_info_first, numTotalFrames_k, 3);
    }


    ToEliminateMatrix2D(submat, NFEATURES, NFEATURES);
    ToEliminateMatrix1D(EigenValSym, NFEATURES);

    outdata_meanCV.close();
    outdata_meanSV.close();
    outdata_var2CV.close();
    outdata_var2SV.close();
    outdata_var3CV.close();
    outdata_var3SV.close();
    outdata_var4CV.close();
    outdata_var4SV.close();

}




void Intrinsec_mean_VectoMapping(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                        vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha)
{

    // The experiments were carried out with LED
//    ofstream outdata_meanW((name_file +"meanW.txt").c_str(), fstream::out); //1 means numero de divisiones del video
//    outdata_meanW.precision(10);
//    outdata_meanW<<  fixed;

    ofstream outdata_LogMeanW((name_file +"LogMeanW.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogMeanW.precision(10);
    outdata_LogMeanW<<  fixed;

    ofstream outdata_LogVar((name_file +"LogVar.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogVar.precision(10);
    outdata_LogVar<<  fixed;

    ofstream outdata_LogVarR((name_file +"LogVarR.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogVarR.precision(10);
    outdata_LogVarR<<  fixed;


    ofstream outdata_LogMeanVar((name_file +"LogMeanVar.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogMeanVar.precision(10);
    outdata_LogMeanVar<<  fixed;

    ofstream outdata_LogMeanVarR((name_file +"LogMeanVarR.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_LogMeanVarR.precision(10);
    outdata_LogMeanVarR<<  fixed;


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

           // cout<< " entro "<< endl;
            int label_k = atoi(returnLabelAction_KTH(vec_videos[k], vec_activities ).c_str()); //-1 para que la primera sea cero
           // int label_k = atoi(returnLabelAction_UT(vec_videos[k], vec_activities ).c_str()); //-1 para que la primera sea cero

            outdata_LogMeanW<< label_k<< " ";
//            outdata_meanW<< label_k<< " ";
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
                    outdata_LogMeanW<< cont_P << ":"<<Log_mean[i][j] << " ";
                    //outdata_meanW<< cont_P << ":"<<VecMap_mean[i][j] << " ";
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
            //outdata_meanW<<endl;
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
    //outdata_meanW.close();
    outdata_LogVar.close();
    outdata_LogVarR.close();

}





void Intrinsec_Mean_Var(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                        vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha)
{

    // The experiments were carried out with LED
    ofstream outdata_meanCV((name_file +"meanCV.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_meanCV.precision(10);
    outdata_meanCV<<  fixed;

    ofstream outdata_meanSV((name_file +"meanSV.txt").c_str(), fstream::out); //1 means numero de divisiones del video
    outdata_meanSV.precision(10);
    outdata_meanSV<<  fixed;


    ofstream outdata_var2CV((name_file + "var2CV.txt").c_str(), fstream::out); //suma de logup cuadrado; peor el cuadrado es en la diagonal
    outdata_var2CV.precision(10);
    outdata_var2CV<<  fixed;

    ofstream outdata_var2SV((name_file + "var2SV.txt").c_str(), fstream::out); //suma de logup cuadrado; peor el cuadrado es en la diagonal
    outdata_var2SV.precision(10);
    outdata_var2SV<<  fixed;

    ofstream outdata_var3CV((name_file + "var3CV.txt").c_str(), fstream::out); // como var2 pero normalizando entre N
    outdata_var3CV.precision(10);
    outdata_var3CV<<  fixed;

    ofstream outdata_var3SV((name_file + "var3SV.txt").c_str(), fstream::out); // como var2 pero normalizando entre N
    outdata_var3SV.precision(10);
    outdata_var3SV<<  fixed;

    ofstream outdata_var4CV((name_file + "var4CV.txt").c_str(), fstream::out); // pasar mean and var a log space y sumar diferencias cudradas
    outdata_var4CV.precision(10);
    outdata_var4CV<<  fixed;

    ofstream outdata_var4SV((name_file + "var4SV.txt").c_str(), fstream::out); // pasar mean and var a log space y sumar diferencias cudradas
    outdata_var4SV.precision(10);
    outdata_var4SV<<  fixed;


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
        //int tot_feat_matrix_k = 12;//((int)featLabels.size()*(int)alphaScales.size())+3; 9 features + x,y,z

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
        ToEliminateMatrix3D(sequenceKin_k, tot_seq_matrix_k, num_tot_traj_k, NFEATURES);


        //cont_tot_DifCov++;
        cout<< " cont_tot_DifCov: "<< cont_tot_DifCov << endl;


///-----------------------------------------------
        if(cont_tot_DifCov >num_mean)
        {

            //cout<< " entro "<< endl;
            int label_k = atoi(returnLabelAction_KTH(vec_videos[k], vec_activities ).c_str()); //-1 para que la primera sea cero

            outdata_meanCV<< label_k<< " ";
            outdata_meanSV<< label_k<< " ";
            outdata_var2CV<< label_k<< " ";
            outdata_var2SV<< label_k<< " ";
            outdata_var3CV<< label_k<< " ";
            outdata_var3SV<< label_k<< " ";
            outdata_var4CV<< label_k<< " ";
            outdata_var4SV<< label_k<< " ";

//            //cout<< " hizo lo de  "<< endl;
            float** _mean = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(_mean, NFEATURES, NFEATURES);
            float** Log_mean = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(Log_mean, NFEATURES, NFEATURES);



            float** _var1 = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(_var1, NFEATURES, NFEATURES);
            float** Log_var1 = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(Log_var1, NFEATURES, NFEATURES);

            float** _var2 = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(_var2, NFEATURES, NFEATURES);
            float** Log_var2 = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(Log_var2, NFEATURES, NFEATURES);




            float** _var4 = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(_var4, NFEATURES, NFEATURES);
            float** Log_var4 = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(Log_var4, NFEATURES, NFEATURES);

            float** Log_var5 = ToCreateMatrix2D(NFEATURES, NFEATURES);
            ToInitMatrix2D(Log_var5, NFEATURES, NFEATURES);


            float distAB =0.0;
            float_grad_mean_Array_Fletcher(arr_cov, _mean, cont_tot_DifCov, NFEATURES);

            cin.ignore();
            cout<< " covariance original" << endl;
            ToPrintMatrix2D(_mean, NFEATURES, NFEATURES);

            cout<< " utilizando sqrt cov "<< endl;
            float_SQRT_A_resp_I_Pennec(_mean, Log_var4,  NFEATURES);
            ToPrintMatrix2D(Log_var4, NFEATURES, NFEATURES);

            for(int i=0; i< NFEATURES; i++)
                for(int j=i+1; j< NFEATURES; j++)
                    Log_var4[i][j] = 1.4142*Log_var4[i][j];

            cout<< " utilizando sqrt cov PONDERADO"<< endl;
            ToPrintMatrix2D(Log_var4, NFEATURES, NFEATURES);


            cout<< " utilizando log de cov "<< endl;
            float_Log_A_resp_I(_mean, Log_mean,  NFEATURES, distAB);
            ToPrintMatrix2D(Log_mean, NFEATURES, NFEATURES);
            //void
            cout<< " utilizando sqrt de la matrix log cov "<< endl;
            float_SQRT_A_resp_I_Pennec(Log_mean, Log_mean,  NFEATURES);
            ToPrintMatrix2D(Log_mean, NFEATURES, NFEATURES);





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






            ///------ outdata_meanCV ---
            int cont_P =1;
            for(int i=0; i< NFEATURES; i++)
                for(int j=i; j< NFEATURES; j++)
                    outdata_meanCV<< cont_P++ << ":"<<Log_mean[i][j] << " ";
            outdata_meanCV<<endl;

            cont_P =1;

            for(int i=0; i< NFEATURES; i++)
                for(int j=i+1; j< NFEATURES; j++)
                    outdata_meanSV<< cont_P++ << ":"<<Log_mean[i][j] << " ";
            outdata_meanSV<<endl;



            ///---var2
            cont_P =1;
            for(int i=0; i< NFEATURES; i++)
                for(int j=i; j< NFEATURES; j++)
                    outdata_var2CV<< cont_P++ << ":"<<Log_var2[i][j] << " ";
            outdata_var2CV<<endl;

            cont_P =1;

            for(int i=0; i< NFEATURES; i++)
                for(int j=i+1; j< NFEATURES; j++)
                    outdata_var2SV<< cont_P++ << ":"<<Log_var2[i][j] << " ";
            outdata_var2SV<<endl;

            ///---var3

            cont_P =1;
            for(int i=0; i< NFEATURES; i++)
                for(int j=i; j< NFEATURES; j++)
                    outdata_var3CV<< cont_P++ << ":"<<Log_var2[i][j]/cont_tot_DifCov << " ";
            outdata_var3CV<<endl;

            cont_P =1;

            for(int i=0; i< NFEATURES; i++)
                for(int j=i+1; j< NFEATURES; j++)
                    outdata_var3SV<< cont_P++ << ":"<<Log_var2[i][j]/cont_tot_DifCov << " ";
            outdata_var3SV<<endl;

            ///---var4
            cont_P =1;
            for(int i=0; i< NFEATURES; i++)
                for(int j=i; j< NFEATURES; j++)
                    outdata_var4CV<< cont_P++ << ":"<<Log_var4[i][j]/cont_tot_DifCov << " ";
            outdata_var4CV<<endl;

            cont_P =1;

            for(int i=0; i< NFEATURES; i++)
                for(int j=i+1; j< NFEATURES; j++)
                    outdata_var4SV<< cont_P++ << ":"<<Log_var4[i][j]/cont_tot_DifCov << " ";
            outdata_var4SV<<endl;




            ToEliminateMatrix2D(_mean, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(Log_mean, NFEATURES, NFEATURES);



            ToEliminateMatrix2D(_var2, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(Log_var2, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(_var4, NFEATURES, NFEATURES);
            ToEliminateMatrix2D(Log_var4, NFEATURES, NFEATURES);




        }
        ToEliminateMatrix3D(arr_cov, numTotalFrames_k, NFEATURES, NFEATURES);
        ToEliminateMatrix2D(arr_info_first, numTotalFrames_k, 3);
    }


    ToEliminateMatrix2D(submat, NFEATURES, NFEATURES);
    ToEliminateMatrix1D(EigenValSym, NFEATURES);

    outdata_meanCV.close();
    outdata_meanSV.close();
    outdata_var2CV.close();
    outdata_var2SV.close();
    outdata_var3CV.close();
    outdata_var3SV.close();
    outdata_var4CV.close();
    outdata_var4SV.close();

}
