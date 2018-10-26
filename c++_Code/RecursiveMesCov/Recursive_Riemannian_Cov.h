#ifndef RECURSIVE_RIEMANNIAN_COV_H
#define RECURSIVE_RIEMANNIAN_COV_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;



void compute_Rec_Riemmannian_Cov_Mean(float** cov_mean_prev, float ** cov_current, float ** cov_mean_actual, float  alpha_val, float distAB, int NFEATURES);
void compute_Rec_Riemmannian_Cov_Var(float** cov_var_prev,  float ** cov_mean_actual,  float ** cov_current, float ** cov_var_actual, float  alpha_val, float distAB, int NFEATURES);

void Recursive_Riemannian_mean(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                    vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha);



#endif
