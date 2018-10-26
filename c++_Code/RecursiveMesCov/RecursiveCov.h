#ifndef RECURSIVECOV_H
#define RECURSIVECOV_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

int num_traj_valides(float** sequenceKin_k, int num_tot_traj_k);

void ToComputeCovmat(float** covmat_K, float** sequenceKin_k, int NFEATURES,int num_tot_traj_k, int heightFrame_k, int widthFrame_k);
bool symPos_label(float* EigenValSym, int NFEATURES);

void Copy_Cov(float** cov_source, float** cov_des, int NFEATURES);
void compute_Rec_Cov_Mean(float** cov_mean_prev, float ** cov_current, float ** cov_mean_actual, float  alpha_val, float distAB, int NFEATURES);
void Recursive_mean(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                    vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha);


void compute_Rec_Cov_Var(float** cov_var_prev,  float ** cov_mean_actual,  float ** cov_current, float ** cov_var_actual, float  alpha_val, float distAB, int NFEATURES);

void compute_Rec_Cov_Min(float** cov_min_prev, float ** cov_current, float ** cov_min_actual, float  alpha_val, float distAB, int NFEATURES);

void compute_Rec_Cov_Max(float** cov_max_prev, float ** cov_current, float ** cov_max_actual, float  alpha_val, float distAB, int NFEATURES);

void compute_Rec_Cov_Dif(float** cov_max, float ** cov_min, float ** cov_diff_actual, float  alpha_val, float distAB, int NFEATURES);

#endif
