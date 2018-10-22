#include "../NR3/nr3.h"
#include "ManageMat_OPS.h"
#include "Geodesics.h"
#include "Interface_NR3_Float.h"
#include "Grad_des_Mean.h"
#include "Riemannian_distances.h"

///compute eigenValues of a Matrix C
void float_eigValues(float** input_C, float* output_C, int NFEATURES)
{

    MatDoub MD_input_C(NFEATURES, NFEATURES, 0.0);
    VecDoub MD_output_C(NFEATURES, 0.0);

    float2Mat(input_C, MD_input_C, NFEATURES, NFEATURES);
    comp_eigValues(MD_input_C, MD_output_C);
    vec2float(MD_output_C, output_C, NFEATURES);


}

///compute eigenDescomposition of a Matrix C= UDU; where D diag matrix of eigenvalues
void float_eigDescomposition(float** input_C, float** Eigen_values, float** Eigen_vectors, int NFEATURES)
{

    MatDoub MD_input_C(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_Eigen_values(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_Eigen_vectors(NFEATURES, NFEATURES, 0.0);

    float2Mat(input_C, MD_input_C,  NFEATURES,  NFEATURES);
    comp_eigDescomposition(MD_input_C, MD_Eigen_values, MD_Eigen_vectors);

    mat2float(MD_Eigen_values, Eigen_values,  NFEATURES,  NFEATURES);
    mat2float(MD_Eigen_vectors, Eigen_vectors,  NFEATURES,  NFEATURES);


}

///compute exponential respect to Identity or simply A= U(exp(Diag(i)))U
void float_Exp_A_resp_I(float** input_A, float** output_A,  int NFEATURES)
{

    MatDoub MD_input_A(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_output_A(NFEATURES, NFEATURES, 0.0);

    float2Mat(input_A, MD_input_A,  NFEATURES,  NFEATURES);
    comp_Exp_A_resp_I(MD_input_A, MD_output_A);

//cout << " version MatDoub Input "<< " NFEATURES: "<< NFEATURES << endl;
//
// displayMat(MD_input_A);
//
//cout << " version MatDoub Output"<< endl;
//
// displayMat(MD_output_A);

    mat2float(MD_output_A, output_A,  NFEATURES,  NFEATURES);

}

///compute logarithm respect to Identity or simply A= U(log(Diag(i)))U
void float_Log_A_resp_I(float** input_A, float** output_A, int NFEATURES, float &distAB)
{

    MatDoub MD_input_A(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_output_A(NFEATURES, NFEATURES, 0.0);

    float2Mat(input_A, MD_input_A,  NFEATURES,  NFEATURES);
    comp_Log_A_resp_I(MD_input_A, MD_output_A, distAB);
    mat2float(MD_output_A, output_A,  NFEATURES,  NFEATURES);

}

///compute sqrt respect to Identity or simply A= U(exp(1/2(log(Diag(i)))))U
void float_SQRT_A_resp_I_Pennec(float** input_A, float** output_A, int NFEATURES)
{

    MatDoub MD_input_A(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_output_A(NFEATURES, NFEATURES, 0.0);

    float2Mat(input_A, MD_input_A, NFEATURES, NFEATURES);
    comp_SQRT_A_resp_I_Pennec(MD_input_A, MD_output_A);
    mat2float(MD_output_A, output_A, NFEATURES, NFEATURES);

}

///compute sqrt respect to Identity according to Fletcher
void float_SQRT_A_resp_I_Fletcher(float** input_A, float** output_A, int NFEATURES)
{

    MatDoub MD_input_A(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_output_A(NFEATURES, NFEATURES, 0.0);

    float2Mat(input_A, MD_input_A, NFEATURES, NFEATURES);
    comp_SQRT_A_resp_I_Fletcher(MD_input_A, MD_output_A);
    mat2float(MD_output_A, output_A, NFEATURES, NFEATURES);

}

/// compute exp_b(A)
void float_Exp_A_resp_B_Pennec(float** input_A, float** input_B, float** output_A, int NFEATURES)
{

    MatDoub MD_input_A(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_input_B(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_output_A(NFEATURES, NFEATURES, 0.0);

    float2Mat(input_A, MD_input_A,  NFEATURES,  NFEATURES);
    float2Mat(input_B, MD_input_B,  NFEATURES,  NFEATURES);

    compute_Exp_A_resp_B_Pennec(MD_input_A, MD_input_B, MD_output_A);
    mat2float(MD_output_A, output_A, NFEATURES,  NFEATURES);

}


/// compute exp_b(A) according to Fletcher
void float_Exp_A_resp_B_Fletcher(float** input_A, float** input_B, float** output_A, int NFEATURES)
{

    MatDoub MD_input_A(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_input_B(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_output_A(NFEATURES, NFEATURES, 0.0);

    float2Mat(input_A, MD_input_A,  NFEATURES,  NFEATURES);
    float2Mat(input_B, MD_input_B,  NFEATURES,  NFEATURES);

    compute_Exp_A_resp_B_Fletcher(MD_input_A, MD_input_B, MD_output_A);
    mat2float(MD_output_A, output_A, NFEATURES,  NFEATURES);

}


/// compute Log_b(A) according to Fletcher
void float_Log_A_resp_B_Fletcher(float** input_A, float** input_B, float** output_A, int NFEATURES, float &distAB)
{

    MatDoub MD_input_A(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_input_B(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_output_A(NFEATURES, NFEATURES, 0.0);

    float2Mat(input_A, MD_input_A,  NFEATURES,  NFEATURES);
    float2Mat(input_B, MD_input_B,  NFEATURES,  NFEATURES);

    compute_Log_A_resp_B_Fletcher(MD_input_A, MD_input_B, MD_output_A, distAB);
    mat2float(MD_output_A, output_A, NFEATURES,  NFEATURES);

}


/// compute varaince = ||Log_b(A)||2 according to Fletcher
void float_Log_A_resp_B_Variance(float** input_A, float** input_B, float** output_A, int NFEATURES, float &distAB)
{

    MatDoub MD_input_A(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_input_B(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_output_A(NFEATURES, NFEATURES, 0.0);

    float2Mat(input_A, MD_input_A,  NFEATURES,  NFEATURES);
    float2Mat(input_B, MD_input_B,  NFEATURES,  NFEATURES);

    compute_Log_A_resp_B_Variance(MD_input_A, MD_input_B, MD_output_A, distAB);
    mat2float(MD_output_A, output_A, NFEATURES,  NFEATURES);

}



/// compute varaince = Log_b(\alpha*A) according to Fletcher
void float_Log_A_resp_B_Alpha(float** input_A, float** input_B, float** output_A, int NFEATURES, float &distAB, float alpha_val)
{

    MatDoub MD_input_A(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_input_B(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_output_A(NFEATURES, NFEATURES, 0.0);

    float2Mat(input_A, MD_input_A,  NFEATURES,  NFEATURES);
    float2Mat(input_B, MD_input_B,  NFEATURES,  NFEATURES);

    compute_Log_A_resp_B_Alpha(MD_input_A, MD_input_B, MD_output_A, distAB, alpha_val);
    mat2float(MD_output_A, output_A, NFEATURES,  NFEATURES);

}



/// compute log_b(A)
void float_Log_A_resp_B_Pennec(float** input_A, float** input_B, float** output_A, int NFEATURES, float &distAB)
{

    MatDoub MD_input_A(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_input_B(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_output_A(NFEATURES, NFEATURES, 0.0);

    float2Mat(input_A, MD_input_A,  NFEATURES,  NFEATURES);
    float2Mat(input_B, MD_input_B,  NFEATURES,  NFEATURES);

    compute_Log_A_resp_B_Pennec(MD_input_A, MD_input_B, MD_output_A, distAB);
    mat2float(MD_output_A, output_A,  NFEATURES,  NFEATURES);

}

void float_grad_mean_Pennec(vector<float**> vect_covMat, float** mean_deg, int NFEATURES)
{


    vector<MatDoub> vect_Matdoub;
    MatDoub MD_mean_deg(NFEATURES, NFEATURES, 0.0);

    for(int i=0; i<vect_covMat.size(); i++)
    {
        MatDoub MD_temp(NFEATURES, NFEATURES, 0.0);
        float2Mat(vect_covMat[i], MD_temp,  NFEATURES,  NFEATURES);
        vect_Matdoub.push_back(MD_temp);
    }



    compute_grad_mean_Pennec( vect_Matdoub,  MD_mean_deg);
    mat2float(MD_mean_deg, mean_deg,  NFEATURES,  NFEATURES);

}



void float_compute_grad_vectorMapping(float*** array_covMat, float** vect_map, int sizeArray, int NFEATURES){


    vector<MatDoub> vect_Matdoub;
    MatDoub MD_vect_map(NFEATURES, NFEATURES, 0.0);

    for(int i=0; i<sizeArray; i++)
    {
        MatDoub MD_temp(NFEATURES, NFEATURES, 0.0);
        float2Mat(array_covMat[i], MD_temp,  NFEATURES,  NFEATURES);
        vect_Matdoub.push_back(MD_temp);
    }



    compute_grad_vectorMapping( vect_Matdoub,  MD_vect_map);
    mat2float(MD_vect_map, vect_map,  NFEATURES,  NFEATURES);



}


void float_grad_mean_Fletcher(vector<float**> vect_covMat, float** mean_deg, int NFEATURES)
{

    vector<MatDoub> vect_Matdoub;
    MatDoub MD_mean_deg(NFEATURES, NFEATURES, 0.0);

    for(int i=0; i<vect_covMat.size(); i++)
    {
        MatDoub MD_temp(NFEATURES, NFEATURES, 0.0);
        float2Mat(vect_covMat[i], MD_temp,  NFEATURES,  NFEATURES);
        vect_Matdoub.push_back(MD_temp);
    }



    compute_grad_mean_Fletcher( vect_Matdoub,  MD_mean_deg);
    mat2float(MD_mean_deg, mean_deg,  NFEATURES,  NFEATURES);


}


void float_grad_mean_Array_Fletcher(float*** array_covMat, float** mean_deg, int sizeArray, int NFEATURES)
{

    vector<MatDoub> vect_Matdoub;
    MatDoub MD_mean_deg(NFEATURES, NFEATURES, 0.0);

    for(int i=0; i<sizeArray; i++)
    {
        MatDoub MD_temp(NFEATURES, NFEATURES, 0.0);
        float2Mat(array_covMat[i], MD_temp,  NFEATURES,  NFEATURES);
        vect_Matdoub.push_back(MD_temp);
    }



    compute_grad_mean_Fletcher( vect_Matdoub,  MD_mean_deg);
    mat2float(MD_mean_deg, mean_deg,  NFEATURES,  NFEATURES);


}


float float_AffineInv_Dist(float** cov_A, float** cov_B, int NFEATURES)
{
    MatDoub MD_A(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_B(NFEATURES, NFEATURES, 0.0);

    float2Mat(cov_A, MD_A,  NFEATURES,  NFEATURES);
    float2Mat(cov_B, MD_B,  NFEATURES,  NFEATURES);

    float dist =  comp_AffineInv_Dist_AID(MD_A, MD_B);
    return dist;
}


float float_Frobenius_Dist(float** cov_A, float** cov_B, int NFEATURES)
{

    MatDoub MD_A(NFEATURES, NFEATURES, 0.0);
    MatDoub MD_B(NFEATURES, NFEATURES, 0.0);

    float2Mat(cov_A, MD_A,  NFEATURES,  NFEATURES);
    float2Mat(cov_B, MD_B,  NFEATURES,  NFEATURES);

    float dist =  comp_Frobenius_Dist_LED(MD_A, MD_B);

    return dist;
}

float float_difussion_dist(float** cov_A, int NFEATURES){

MatDoub MD_A(NFEATURES, NFEATURES, 0.0);
float2Mat(cov_A, MD_A,  NFEATURES,  NFEATURES);

float dist =  difussion_dist(MD_A);
return dist;
}

