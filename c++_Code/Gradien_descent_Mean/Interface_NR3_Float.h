#ifndef INTERNR3FLOAT_H
#define INTERNR3FLOAT_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;

///compute eigenValues of a Matrix C
void float_eigValues(float** input_C, float* output_C, int NFEATURES);

///compute eigenDescomposition of a Matrix C= UDU; where D diag matrix of eigenvalues
void float_eigDescomposition(float** input_C, float** Eigen_values, float** Eigen_vectors, int NFEATURES);
///compute exponential respect to Identity or simply A= U(exp(Diag(i)))U
void float_Exp_A_resp_I(float** input_A, float** output_A,  int NFEATURES);

///compute logarithm respect to Identity or simply A= U(log(Diag(i)))U
void float_Log_A_resp_I(float** input_A, float** output_A, int NFEATURES, float &distAB);

///compute sqrt respect to Identity or simply A= U(exp(1/2(log(Diag(i)))))U
void float_SQRT_A_resp_I_Pennec(float** input_A, float** output_A, int NFEATURES);
void float_SQRT_A_resp_I_Fletcher(float** input_A, float** output_A, int NFEATURES);

/// compute exp_b(A)
void float_Exp_A_resp_B_Pennec(float** input_A, float** input_B, float** output_A, int NFEATURES);
void float_Exp_A_resp_B_Fletcher(float** input_A, float** input_B, float** output_A, int NFEATURES);

/// compute log_b(A)
void float_Log_A_resp_B_Pennec(float** input_A, float** input_B, float** output_A, int NFEATURES, float &distAB);
void float_Log_A_resp_B_Fletcher(float** input_A, float** input_B, float** output_A, int NFEATURES, float &distAB);


/// compute gradient descendent mean
void float_grad_mean_Pennec(vector<float**> vect_covMat, float** mean_deg, int NFEATURES);
void float_grad_mean_Fletcher(vector<float**> vect_covMat, float** mean_deg, int NFEATURES);

void float_grad_mean_Array_Fletcher(float*** array_covMat, float** mean_deg, int sizeArray, int NFEATURES);


/// compute the vector mapping
void float_compute_grad_vectorMapping(float*** array_covMat, float** vect_map, int sizeArray, int NFEATURES);


/// compute variance
void float_Log_A_resp_B_Variance(float** input_A, float** input_B, float** output_A, int NFEATURES, float &distAB);
void float_Log_A_resp_B_Alpha(float** input_A, float** input_B, float** output_A, int NFEATURES, float &distAB, float alpha_val);


/// compute distances
float float_Frobenius_Dist(float** cov_A, float** cov_B,int NFEATURES);
float float_AffineInv_Dist(float** cov_A, float** cov_B,int  NFEATURES);
float float_difussion_dist(float** cov_A, int NFEATURES);

#endif
