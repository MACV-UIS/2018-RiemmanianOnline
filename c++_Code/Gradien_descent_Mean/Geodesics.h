#ifndef GEODESICS_H
#define GEODESICS_H

///compute eigenValues of a Matrix C
void comp_eigValues(MatDoub &input_C, VecDoub &Eigen_values);

///compute eigenDescomposition of a Matrix C= UDU; where D diag matrix of eigenvalues
void comp_eigDescomposition(MatDoub &input_C, MatDoub &Eigen_values, MatDoub &Eigen_vectors);

///compute exponential respect to Identity or simply A= U(exp(Diag(i)))U
void comp_Exp_A_resp_I(MatDoub &input_A, MatDoub &output_A);

///compute logarithm respect to Identity or simply A= U(log(Diag(i)))U
void comp_Log_A_resp_I(MatDoub &input_A, MatDoub &output_A, float &distAB);

///compute sqrt respect to Identity or simply A= U(exp(1/2(log(Diag(i)))))U
void comp_SQRT_A_resp_I_Pennec(MatDoub &input_A, MatDoub &output_A);
void comp_SQRT_A_resp_I_Fletcher(MatDoub &input_A, MatDoub &output_A);
/// compute exp_b(A)
void compute_Exp_A_resp_B_Pennec(MatDoub &input_A, MatDoub &input_B, MatDoub &output_A);
void compute_Exp_A_resp_B_Fletcher(MatDoub &input_A, MatDoub &input_B, MatDoub &output_A);
/// compute log_b(A)
void compute_Log_A_resp_B_Pennec(MatDoub &input_A, MatDoub &input_B, MatDoub &output_A, float &distAB);
void compute_Log_A_resp_B_Fletcher(MatDoub &input_A, MatDoub &input_B, MatDoub &output_A, float &distAB);

/// compute variance
void compute_Log_A_resp_B_Variance(MatDoub &input_A, MatDoub &input_B, MatDoub &output_A, float &distAB);

void compute_Log_A_resp_B_Alpha(MatDoub &input_A, MatDoub &input_B, MatDoub &output_A, float &distAB, float alpha_val);
#endif // GEODESICS_H



