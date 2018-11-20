#include <iostream>
#include <math.h>

#include "../NR3/nr3.h"
#include "../NR3/eigen_sym.h"
#include "../NR3/eigen_unsym.h"
#include "../NR3/ludcmp.h"

#include "ManageMat_OPS.h"
#include "Geodesics.h"

using namespace std;


///compute eigenValues
void comp_eigValues(MatDoub &input_C, VecDoub &Eigen_values)
{

    int NFEATURES = input_C.nrows();
    Symmeig s(input_C);

    for (int i = 0; i < NFEATURES; i++)
    {
        Eigen_values[i] = s.d[i];

    }

}

///compute eigenDescomposition of a Matrix C= UDU; where D diag matrix of eigenvalues
void comp_eigDescomposition(MatDoub &input_C, MatDoub &Eigen_values, MatDoub &Eigen_vectors)
{

    int NFEATURES = input_C.nrows();
    Symmeig s(input_C);

    for (int i = 0; i < NFEATURES; i++)
    {
        for(int j=0; j<NFEATURES; j++)
        {
            Eigen_vectors[i][j] = s.z[i][j];
            if(i==j) Eigen_values[i][j] = s.d[i];
        }
    }

}

///compute exponential respect to Identity or simply A= U(exp(Diag(i)))U
void comp_Exp_A_resp_I(MatDoub &input_A, MatDoub &output_A)
{

    int NFEATURES = input_A.nrows();

    MatDoub Meigen_values(NFEATURES, NFEATURES, 0.0);
    MatDoub Meigen_vectors(NFEATURES, NFEATURES, 0.0);
    MatDoub Meigen_vectors_trans(NFEATURES, NFEATURES, 0.0);
    MatDoub Meucl_temp(NFEATURES, NFEATURES, 0.0);

    comp_eigDescomposition(input_A, Meigen_values, Meigen_vectors);

    //cout<< " eigen values "<< endl;

    //displayMat(Meigen_values);

    for (int i = 0; i < NFEATURES; i++)
    {
        for(int j=0; j<NFEATURES; j++)
        {
            if(i==j)
            {
                if(Meigen_values[i][j]<10)  // 708 because protcolo IEEE-compatible type
                {
                    Meigen_values[i][j] = exp(Meigen_values[i][j]);
                }
                else
                {
                    Meigen_values[i][j] = exp(10);
                }
            }
        }
    }

//
//    cout<< " eigen values after exp "<< endl;
//
//    displayMat(Meigen_values);


    transpose_of_any_matrix(Meigen_vectors, Meigen_vectors_trans);
    matrix_multiplication(Meigen_vectors, Meigen_values, Meucl_temp);
    matrix_multiplication(Meucl_temp, Meigen_vectors_trans, output_A);

}

///compute logarithm respect to Identity or simply A= U(log(Diag(i)))U
void comp_Log_A_resp_I(MatDoub &input_A, MatDoub &output_A, float &distAB)
{

    int NFEATURES = input_A.nrows();
    MatDoub Meigen_values(NFEATURES, NFEATURES, 0.0);
    MatDoub Meigen_vectors(NFEATURES, NFEATURES, 0.0);
    MatDoub Meigen_vectors_trans(NFEATURES, NFEATURES, 0.0);
    MatDoub Meucl_temp(NFEATURES, NFEATURES, 0.0);

    comp_eigDescomposition(input_A, Meigen_values, Meigen_vectors);
    distAB =0;

    for(int j=0; j<NFEATURES; j++)
    {

        if(Meigen_values[j][j] > 0)
        {
            Meigen_values[j][j] = log(Meigen_values[j][j]);
            distAB = distAB + pow(Meigen_values[j][j],2) ;
        }
        else
        {
            Meigen_values[j][j] = log(0.001);
            distAB = distAB + pow(Meigen_values[j][j],2) ;
            cout<< " ****************** !!!!! ATTENTION LOG (0) ou <0 Log_A_resp_I!!!! ***************** "<< endl;
        }

    }



    transpose_of_any_matrix(Meigen_vectors, Meigen_vectors_trans);
    matrix_multiplication(Meigen_vectors, Meigen_values, Meucl_temp);
    matrix_multiplication(Meucl_temp, Meigen_vectors_trans, output_A);

}

///compute sqrt respect to Identity or simply A= U(exp(1/2(log(Diag(i)))))U
void comp_SQRT_A_resp_I_Pennec(MatDoub &input_A, MatDoub &output_A)
{

    int NFEATURES = input_A.nrows();
    MatDoub Meigen_values(NFEATURES, NFEATURES, 0.0);
    MatDoub Meigen_vectors(NFEATURES, NFEATURES, 0.0);
    MatDoub Meigen_vectors_trans(NFEATURES, NFEATURES, 0.0);
    MatDoub Meucl_temp(NFEATURES, NFEATURES, 0.0);

    comp_eigDescomposition(input_A, Meigen_values, Meigen_vectors);

    for (int i = 0; i < NFEATURES; i++)
    {
        for(int j=0; j<NFEATURES; j++)
        {

            if(i==j)
            {
                if(log(Meigen_values[i][j])/2 < 10)
                {
                    Meigen_values[i][j] = exp(log(Meigen_values[i][j])/2);
                }
                else
                {
                    Meigen_values[i][j] = exp(10);
                }

            }
        }
    }

    transpose_of_any_matrix(Meigen_vectors, Meigen_vectors_trans);
    matrix_multiplication(Meigen_vectors, Meigen_values, Meucl_temp);
    matrix_multiplication(Meucl_temp, Meigen_vectors_trans, output_A);


}


void comp_SQRT_A_resp_I_Fletcher(MatDoub &input_A, MatDoub &output_A)
{

    int NFEATURES = input_A.nrows();
    MatDoub Meigen_values(NFEATURES, NFEATURES, 0.0);
    MatDoub Meigen_vectors(NFEATURES, NFEATURES, 0.0);
    MatDoub Meigen_vectors_trans(NFEATURES, NFEATURES, 0.0);
    MatDoub Meucl_temp(NFEATURES, NFEATURES, 0.0);

    comp_eigDescomposition(input_A, Meigen_values, Meigen_vectors);

    for (int i = 0; i < NFEATURES; i++)
    {
        for(int j=0; j<NFEATURES; j++)
        {

            if(i==j)
            {
                if(Meigen_values[i][j]>0)
                    Meigen_values[i][j] = sqrt(Meigen_values[i][j]);
            }
        }
    }

    matrix_multiplication(Meigen_vectors, Meigen_values, output_A);


}

/// compute exp_b(A)
void compute_Exp_A_resp_B_Pennec(MatDoub &input_A, MatDoub &input_B, MatDoub &output_A)
{
    int NFEATURES = input_B.nrows();
    MatDoub Mat_sqrtB(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_sqrtB_inv(NFEATURES, NFEATURES, 0.0);

    comp_SQRT_A_resp_I_Pennec(input_B, Mat_sqrtB);
    inverse(Mat_sqrtB, Mat_sqrtB_inv);

    MatDoub Mat_tem(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_input_exp(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_out_exp(NFEATURES, NFEATURES, 0.0);

    matrix_multiplication(Mat_sqrtB_inv, input_A, Mat_tem);
    matrix_multiplication(Mat_tem, Mat_sqrtB_inv, Mat_input_exp);

    comp_Exp_A_resp_I(Mat_input_exp, Mat_out_exp);

    MatDoub Mat_temp_2(NFEATURES, NFEATURES, 0.0);

    matrix_multiplication(Mat_sqrtB, Mat_out_exp, Mat_temp_2);
    matrix_multiplication(Mat_temp_2, Mat_sqrtB, output_A);

}

/// compute exp_b(A) according to Fletcher
void compute_Exp_A_resp_B_Fletcher(MatDoub &input_A, MatDoub &input_B, MatDoub &output_A)
{

    int NFEATURES = input_B.nrows();
    MatDoub Mat_sqrtB(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_sqrtB_inv(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_sqrtB_trans(NFEATURES, NFEATURES, 0.0);



    comp_SQRT_A_resp_I_Fletcher(input_B, Mat_sqrtB);
    inverse(Mat_sqrtB, Mat_sqrtB_inv);
    transpose_of_any_matrix(Mat_sqrtB_inv, Mat_sqrtB_trans);

    MatDoub Mat_temp_Y(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_Y(NFEATURES, NFEATURES, 0.0);

    matrix_multiplication(Mat_sqrtB_inv, input_A, Mat_temp_Y);
    matrix_multiplication(Mat_temp_Y, Mat_sqrtB_trans, Mat_Y);

    MatDoub Meigen_values_Y(NFEATURES, NFEATURES, 0.0);
    MatDoub Meigen_vectors_Y(NFEATURES, NFEATURES, 0.0);
    comp_eigDescomposition(Mat_Y, Meigen_values_Y, Meigen_vectors_Y);

    for(int j=0; j<NFEATURES; j++)
    {
        if(Meigen_values_Y[j][j]  < 10)
        {
            Meigen_values_Y[j][j] = exp(Meigen_values_Y[j][j]);
        }
        else
        {
            Meigen_values_Y[j][j] = exp(10);
        }
    }

    MatDoub Mat_gv(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_gv_trans(NFEATURES, NFEATURES, 0.0);


    matrix_multiplication(Mat_sqrtB, Meigen_vectors_Y, Mat_gv);
    transpose_of_any_matrix(Mat_gv, Mat_gv_trans);

    MatDoub Mat_gve_temp(NFEATURES, NFEATURES, 0.0);
    matrix_multiplication(Mat_gv, Meigen_values_Y, Mat_gve_temp);
    matrix_multiplication(Mat_gve_temp, Mat_gv_trans, output_A);

}


/// compute Log_b(A) according to Fletcher
void compute_Log_A_resp_B_Variance(MatDoub &input_A, MatDoub &input_B, MatDoub &output_A, float &distAB)
{

    int NFEATURES = input_B.nrows();
    MatDoub Mat_sqrtB(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_sqrtB_inv(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_sqrtB_trans(NFEATURES, NFEATURES, 0.0);



    comp_SQRT_A_resp_I_Fletcher(input_B, Mat_sqrtB);
    inverse(Mat_sqrtB, Mat_sqrtB_inv);
    transpose_of_any_matrix(Mat_sqrtB_inv, Mat_sqrtB_trans);

    MatDoub Mat_temp_Y(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_Y(NFEATURES, NFEATURES, 0.0);

    matrix_multiplication(Mat_sqrtB_inv, input_A, Mat_temp_Y);
    matrix_multiplication(Mat_temp_Y, Mat_sqrtB_trans, Mat_Y);

    MatDoub Meigen_values_Y(NFEATURES, NFEATURES, 0.0);
    MatDoub Meigen_vectors_Y(NFEATURES, NFEATURES, 0.0);
    comp_eigDescomposition(Mat_Y, Meigen_values_Y, Meigen_vectors_Y);

    distAB = 0;

    for(int j=0; j<NFEATURES; j++)
    {
        if(Meigen_values_Y[j][j]  >0)
        {
            Meigen_values_Y[j][j] = pow(log(Meigen_values_Y[j][j]),2);// log(pow(Meigen_values_Y[j][j],2));
            distAB = distAB + pow(Meigen_values_Y[j][j],2) ;
        }
        else
        {
            Meigen_values_Y[j][j] = log(0.001);
            distAB = distAB + pow(Meigen_values_Y[j][j],2) ;
            //cout<< "error: valor para log <= 0 Log_A_resp_B_ " << endl;
            //cout<< " viene de la matrix: "<< endl;
            displayMat(input_A);
            //cout << " con valores eigen values de: "<< endl;
            displayMat(Meigen_values_Y);
            //cin.ignore();


        }
    }

    MatDoub Mat_gv(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_gv_trans(NFEATURES, NFEATURES, 0.0);


    matrix_multiplication(Mat_sqrtB, Meigen_vectors_Y, Mat_gv);
    transpose_of_any_matrix(Mat_gv, Mat_gv_trans);

    MatDoub Mat_gve_temp(NFEATURES, NFEATURES, 0.0);
    matrix_multiplication(Mat_gv, Meigen_values_Y, Mat_gve_temp);
    matrix_multiplication(Mat_gve_temp, Mat_gv_trans, output_A);

}



/// compute Log_b(A) according to Fletcher but weighted according to an alpha
void compute_Log_A_resp_B_Alpha(MatDoub &input_A, MatDoub &input_B, MatDoub &output_A, float &distAB, float alpha_val)
{

    int NFEATURES = input_B.nrows();
    MatDoub Mat_sqrtB(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_sqrtB_inv(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_sqrtB_trans(NFEATURES, NFEATURES, 0.0);



    comp_SQRT_A_resp_I_Fletcher(input_B, Mat_sqrtB);
    inverse(Mat_sqrtB, Mat_sqrtB_inv);
    transpose_of_any_matrix(Mat_sqrtB_inv, Mat_sqrtB_trans);

    MatDoub Mat_temp_Y(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_Y(NFEATURES, NFEATURES, 0.0);

    matrix_multiplication(Mat_sqrtB_inv, input_A, Mat_temp_Y);
    matrix_multiplication(Mat_temp_Y, Mat_sqrtB_trans, Mat_Y);

    MatDoub Meigen_values_Y(NFEATURES, NFEATURES, 0.0);
    MatDoub Meigen_vectors_Y(NFEATURES, NFEATURES, 0.0);
    comp_eigDescomposition(Mat_Y, Meigen_values_Y, Meigen_vectors_Y);

    distAB = 0;

    for(int j=0; j<NFEATURES; j++)
    {
        if(Meigen_values_Y[j][j]  >0)
        {
            Meigen_values_Y[j][j] = alpha_val*log(Meigen_values_Y[j][j]) ;// log(pow(Meigen_values_Y[j][j],2));
            distAB = distAB + pow(Meigen_values_Y[j][j],2) ;
        }
        else
        {
            Meigen_values_Y[j][j] = log(0.001);
            distAB = distAB + pow(Meigen_values_Y[j][j],2) ;
            //cout<< "error: valor para log <= 0 Log_A_resp_B_ " << endl;
            //cout<< " viene de la matrix: "<< endl;
            displayMat(input_A);
            //cout << " con valores eigen values de: "<< endl;
            displayMat(Meigen_values_Y);
            //cin.ignore();


        }
    }

    MatDoub Mat_gv(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_gv_trans(NFEATURES, NFEATURES, 0.0);


    matrix_multiplication(Mat_sqrtB, Meigen_vectors_Y, Mat_gv);
    transpose_of_any_matrix(Mat_gv, Mat_gv_trans);

    MatDoub Mat_gve_temp(NFEATURES, NFEATURES, 0.0);
    matrix_multiplication(Mat_gv, Meigen_values_Y, Mat_gve_temp);
    matrix_multiplication(Mat_gve_temp, Mat_gv_trans, output_A);

}


/// compute Log_b(A) according to Fletcher
void compute_Log_A_resp_B_Fletcher(MatDoub &input_A, MatDoub &input_B, MatDoub &output_A, float &distAB)
{

    int NFEATURES = input_B.nrows();
    MatDoub Mat_sqrtB(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_sqrtB_inv(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_sqrtB_trans(NFEATURES, NFEATURES, 0.0);



    comp_SQRT_A_resp_I_Fletcher(input_B, Mat_sqrtB);
    inverse(Mat_sqrtB, Mat_sqrtB_inv);
    transpose_of_any_matrix(Mat_sqrtB_inv, Mat_sqrtB_trans);

    MatDoub Mat_temp_Y(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_Y(NFEATURES, NFEATURES, 0.0);

    matrix_multiplication(Mat_sqrtB_inv, input_A, Mat_temp_Y);
    matrix_multiplication(Mat_temp_Y, Mat_sqrtB_trans, Mat_Y);

    MatDoub Meigen_values_Y(NFEATURES, NFEATURES, 0.0);
    MatDoub Meigen_vectors_Y(NFEATURES, NFEATURES, 0.0);
    comp_eigDescomposition(Mat_Y, Meigen_values_Y, Meigen_vectors_Y);

    distAB = 0;

    for(int j=0; j<NFEATURES; j++)
    {
        if(Meigen_values_Y[j][j]  >0)
        {
            Meigen_values_Y[j][j] =     (Meigen_values_Y[j][j]);
            distAB = distAB + pow(Meigen_values_Y[j][j],2) ;
        }
        else
        {
            Meigen_values_Y[j][j] = log(0.001);
            distAB = distAB + pow(Meigen_values_Y[j][j],2) ;
            //cout<< "error: valor para log <= 0 Log_A_resp_B_ " << endl;
            //cout<< " viene de la matrix: "<< endl;
            displayMat(input_A);
            //cout << " con valores eigen values de: "<< endl;
            displayMat(Meigen_values_Y);
            //cin.ignore();


        }
    }

    MatDoub Mat_gv(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_gv_trans(NFEATURES, NFEATURES, 0.0);


    matrix_multiplication(Mat_sqrtB, Meigen_vectors_Y, Mat_gv);
    transpose_of_any_matrix(Mat_gv, Mat_gv_trans);

    MatDoub Mat_gve_temp(NFEATURES, NFEATURES, 0.0);
    matrix_multiplication(Mat_gv, Meigen_values_Y, Mat_gve_temp);
    matrix_multiplication(Mat_gve_temp, Mat_gv_trans, output_A);

}



/// compute log_b(A)
void compute_Log_A_resp_B_Pennec(MatDoub &input_A, MatDoub &input_B, MatDoub &output_A, float &distAB)
{
    int NFEATURES = input_B.nrows();
    MatDoub Mat_sqrtB(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_sqrtB_inv(NFEATURES, NFEATURES, 0.0);

    comp_SQRT_A_resp_I_Pennec(input_B, Mat_sqrtB);
    inverse(Mat_sqrtB, Mat_sqrtB_inv);

    MatDoub Mat_tem(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_input_exp(NFEATURES, NFEATURES, 0.0);
    MatDoub Mat_out_exp(NFEATURES, NFEATURES, 0.0);

    matrix_multiplication(Mat_sqrtB_inv, input_A, Mat_tem);
    matrix_multiplication(Mat_tem, Mat_sqrtB_inv, Mat_input_exp);


    comp_Log_A_resp_I(Mat_input_exp, Mat_out_exp, distAB );

    MatDoub Mat_temp_2(NFEATURES, NFEATURES, 0.0);

    matrix_multiplication(Mat_sqrtB, Mat_out_exp, Mat_temp_2);
    matrix_multiplication(Mat_temp_2, Mat_sqrtB, output_A);

}
