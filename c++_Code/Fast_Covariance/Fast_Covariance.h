/**---------------------------------------------------------------------------------------
Author: fabio.martinez-carillo@ensta-paristech.fr
Description: implementation of fast covariance implementation.
Date: 08/12/2015
----------------------------------------------------------------------------------------*/

#ifndef FAST_COVARIANCE_H
#define FAST_COVARIANCE_H

void computeIntegral_matrices(float*** feature_matrix, int NFEATURES, int height, int width, float*** product_matrix_Q );
void covariance(int LT_x, int LT_y, int RB_x, int RB_y, float ***fsum, float*** prsum, float **covmat, float ** f_label, int NFEATURES );
float value_sumimage(int lefttop_x, int lefttop_y, int rightbot_x, int rightbot_y, float **sum_image);
#endif // FAST_COVARIANCE_H
