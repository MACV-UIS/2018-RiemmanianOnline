#include <iostream>
#include <math.h>

#include "ManageMat_OPS.h"
#include "Geodesics.h"
#include "Riemannian_distances.h"
#include "Grad_des_Mean.h"

using namespace std;


void sum_mat( MatDoub &mat_iter,  MatDoub &mat_temp)
{
    int NFEATURES = mat_iter.nrows();
    for(int i=0; i<NFEATURES; i++)
        for(int j=0; j<NFEATURES; j++)
            mat_iter[i][j] = mat_iter[i][j] + mat_temp[i][j];


}

void div_mat(MatDoub &mat_iter, int div_index)
{

    int NFEATURES = mat_iter.nrows();
    for(int i=0; i<NFEATURES; i++)
        for(int j=0; j<NFEATURES; j++)
            mat_iter[i][j] = mat_iter[i][j]/div_index;
}

void compute_grad_mean_Pennec(vector<MatDoub> vect_Matdoub, MatDoub &mean_deg)
{

    float error_val =0.001;
    int NFEATURES = mean_deg.nrows();
    for(int l=0; l<NFEATURES; l++)mean_deg[l][l]=1.0; // comp_identity

    MatDoub mat_iter(NFEATURES, NFEATURES,0.0);
    MatDoub mat_temp(NFEATURES, NFEATURES,0.0);
    float distAB =0;
    float dist_temp =0;
    int cont=0;

    do
    {
        for(int i=0; i<vect_Matdoub.size(); i++)
        {
            compute_Log_A_resp_B_Pennec(vect_Matdoub[i], mean_deg, mat_temp, distAB);
            sum_mat(mat_iter, mat_temp);
        }

        div_mat(mat_iter, vect_Matdoub.size());
        compute_Exp_A_resp_B_Pennec(mat_iter, mean_deg, mean_deg);

//        cout<< " mean_deg  "<< endl;
//        displayMat(mean_deg);
//        cin.ignore();

        dist_temp =difussion_dist(mat_iter);
        cont++;
        //cout << " measure of error: "<<difussion_dist(mat_iter)<< " cont: " << cont++ << "total mat:  " << vect_Matdoub.size()<<  endl;
        restart_mat(mat_iter);

    }
    while( dist_temp> error_val  && cont <500 ); //mat_iter
    cout<< "numero de iteaciones: "<<cont<< " error "<< dist_temp<< endl;
}




void compute_grad_vectorMapping(vector<MatDoub> vect_Matdoub, MatDoub &vect_map)
{

    float error_val =0.001;
    int NFEATURES = vect_map.nrows();
    MatDoub mean_deg(NFEATURES, NFEATURES,0.0);

    MatDoub mean_sqrt(NFEATURES, NFEATURES,0.0);
    MatDoub mean_sqrt_trans(NFEATURES, NFEATURES,0.0);

    MatDoub mat_iter(NFEATURES, NFEATURES,0.0);
    MatDoub mat_temp(NFEATURES, NFEATURES,0.0);
    MatDoub mat_temp_multi(NFEATURES, NFEATURES,0.0);


    for(int l=0; l<NFEATURES; l++)mean_deg[l][l]=1.0; // comp_identity

    float distAB =0;
    float dist_temp =0;
    int cont=0;

    do
    {
        for(int i=0; i<vect_Matdoub.size(); i++)
        {
            compute_Log_A_resp_B_Pennec(vect_Matdoub[i], mean_deg, mat_temp, distAB);
            sum_mat(mat_iter, mat_temp);
        }

        div_mat(mat_iter, vect_Matdoub.size());

        comp_SQRT_A_resp_I_Pennec(mean_deg, mean_sqrt);
        transpose_of_any_matrix(mean_sqrt, mean_sqrt_trans);

        matrix_multiplication(mean_sqrt_trans, mat_iter, mat_temp_multi);
        matrix_multiplication(mat_temp_multi, mean_sqrt_trans, vect_map);

        compute_Exp_A_resp_B_Pennec(mat_iter, mean_deg, mean_deg);

        dist_temp =difussion_dist(mat_iter);
        cont++;
        //cout << " measure of error: "<<difussion_dist(mat_iter)<< " cont: " << cont++ << "total mat:  " << vect_Matdoub.size()<<  endl;
        restart_mat(mat_iter);

    }
    while( dist_temp> error_val  && cont <500 ); //mat_iter
    cout<< "numero de iteaciones: "<<cont<< " error "<< dist_temp<< endl;
}





void  compute_grad_mean_Fletcher(vector<MatDoub> vect_Matdoub, MatDoub &mean_deg)
{

    float error_val =0.001;
    int NFEATURES = mean_deg.nrows();
    for(int l=0; l<NFEATURES; l++)mean_deg[l][l]=1.0; // comp_identity

    MatDoub mat_iter(NFEATURES, NFEATURES,0.0);
    MatDoub mat_temp(NFEATURES, NFEATURES,0.0);
    float distAB =0;
    float dist_temp =0;
    int cont=0;

    do
    {
        // cout<< " numero total de matrices "<< vect_Matdoub.size() << endl;

        for(int i=0; i<vect_Matdoub.size(); i++)
        {
            //cout<< "determinant of: "<< compute_determinant(vect_Matdoub[i])<< endl;
            compute_Log_A_resp_B_Fletcher(vect_Matdoub[i], mean_deg, mat_temp, distAB);
            sum_mat(mat_iter, mat_temp);
        }

        div_mat(mat_iter, vect_Matdoub.size());
        compute_Exp_A_resp_B_Fletcher(mat_iter, mean_deg, mean_deg);

//        cout<< " mean_deg  "<< endl;
//        displayMat(mean_deg);
//        cin.ignore();

        dist_temp =difussion_dist(mat_iter);
        cont++;
        cout << " measure of error: "<<difussion_dist(mat_iter)<< " cont: " << cont << " total mat:  " << vect_Matdoub.size()<<  endl;
        restart_mat(mat_iter);

    }
    while( dist_temp> error_val && cont <50 ); //mat_iter
    cout<< "numero de iteaciones: "<<cont<< " error "<< dist_temp<< endl;

}

