#include "Fast_Covariance.h"
#include <iostream>
using namespace std;



///* To compute to P = feature_matrix y Q= product_matrix_Q
void computeIntegral_matrices(float*** feature_matrix, int NFEATURES, int height, int width, float*** product_matrix_Q ){

///----------------------------------------------------------------------------------
/// **************** Matrix multiplication  ******************
///----------------------------------------------------------------------------------
	int pos = 0;
	for (int k = 0; k < NFEATURES; k++) {
		for (int l = k; l < NFEATURES; l++) {
			for (int i = 0; i < width; i++) {
				for (int j = 0; j <height ; j++) {
					product_matrix_Q[pos][i][j] = (feature_matrix[l][i][j])
							* (feature_matrix[k][i][j]);
				}
			}
			pos++;
		}
	}


///----------------------------------------------------------------------------------
/// **************** integral image of features region matrix = Q  ******************
///----------------------------------------------------------------------------------

	for (int x = 0; x < ((NFEATURES * (NFEATURES + 1)) / 2); x++) {
		for (int i = 1; i <width; i++)
			product_matrix_Q[x][i][0] += product_matrix_Q[x][i - 1][0];
		for (int j = 1; j < height; j++)
			product_matrix_Q[x][0][j] += product_matrix_Q[x][0][j - 1];

		for (int i = 1; i < width; i++)
			for (int j = 1; j < height; j++)
				product_matrix_Q[x][i][j] += (product_matrix_Q[x][i][j - 1]
						+ product_matrix_Q[x][i - 1][j])
						- product_matrix_Q[x][i - 1][j - 1];
	}

///----------------------------------------------------------------------------------
/// **************** integral image of features region matrix = Ptr ******************
///----------------------------------------------------------------------------------
	for (int x = 0; x < NFEATURES; x++) {
		for (int i = 1; i < width; i++)
			feature_matrix[x][i][0] += feature_matrix[x][i - 1][0];
		for (int j = 1; j < height; j++)
			feature_matrix[x][0][j] += feature_matrix[x][0][j - 1];

		for (int i = 1; i < width; i++)
			for (int j = 1; j < height; j++)
				feature_matrix[x][i][j] += (feature_matrix[x][i][j - 1]
						+ feature_matrix[x][i - 1][j]) - feature_matrix[x][i
						- 1][j - 1];
	}
///-----------------------------------------------------------------------------
}


///* To compute Cov

void covariance(int LT_x, int LT_y, int RB_x, int RB_y, float ***fsum, float*** prsum, float **covmat, float ** f_label, int NFEATURES ) {


	int width  = RB_x - LT_x;
	int height = RB_y - LT_y;

    ///------------------------------------------------------
    ///-------      number of availables samples ------------
    ///-------------------------------------------------------
        int n_samples = 0;
        for(int x= LT_x; x< RB_x; x++)
            for(int y= LT_y; y< RB_y; y++)
                if(f_label[x][y]==1)n_samples++;





    //cout<< "number of samples: " <<  n_samples<< " R width: " << width << " R height: "<< height <<endl<< endl;
    ///--------------------------------------------------------
    ///-------------- Covariance computation ------------------
    ///--------------------------------------------------------
	float* mean = new float[NFEATURES];
	int pos = 0;

	for (int i = 0; i < NFEATURES; i++) {
		mean[i] = value_sumimage(LT_x, LT_y, RB_x, RB_y, fsum[i]);
		//mean[i] = mean[i]/n_samples; // if I want the mean value
	}

	for (int i = 0; i < NFEATURES; i++) {
		for (int j = i; j < NFEATURES; j++) {
			covmat[i][j] = value_sumimage(LT_x, LT_y, RB_x, RB_y, prsum[pos]) - ((mean[i])
					* (mean[j])/n_samples );
			covmat[i][j]= covmat[i][j]/((n_samples)-1);
			if (i != j) {  //copy of values
				covmat[j][i] = covmat[i][j];
			}
			pos++;
		}
	}


}




float value_sumimage(int lefttop_x, int lefttop_y, int rightbot_x, int rightbot_y, float **sum_image) {
	//This code will give the value  of
	//summation of the values of that particular feature
	// NOTE you have to enter co-ordinates of the left top pixel
	// and right bottom pixel
	//both pixels are included in the summation
	float value = 0.0;

    value = sum_image[rightbot_x][rightbot_y]
				+ sum_image[lefttop_x][lefttop_y]
				- sum_image[rightbot_x][lefttop_y]
				- sum_image[lefttop_x][rightbot_y];



//	if ( lefttop_y  < 1 && lefttop_x  >= 1) {
//		value = sum_image[rightbot_x][rightbot_y]
//				- sum_image[lefttop_x - 1][rightbot_y]; //-1 porque el index empieza en zero y va hasta width and height
//	}


//	if ((lefttop_x - 1) < 0 && (lefttop_y - 1) >= 0) {
//		value = sum_image[rightbot_x][rightbot_y]
//				- sum_image[rightbot_x][lefttop_y - 1];
//
//	}
//
//
//	if ((lefttop_x - 1) < 0 && (lefttop_y - 1) < 0) {
//		value = sum_image[rightbot_x][rightbot_y];
//
//	}
//
//
//	if ((lefttop_x - 1) >= 0 && (lefttop_y - 1) >= 0) {
//		value = sum_image[rightbot_x][rightbot_y]
//				+ sum_image[lefttop_x - 1][lefttop_y - 1]
//				- sum_image[rightbot_x][lefttop_y - 1]
//				- sum_image[lefttop_x- 1][rightbot_y];
//
//
//	}
	return value;
}
