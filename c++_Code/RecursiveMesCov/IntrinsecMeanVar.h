#ifndef INTRINSECMEANVAR_H
#define INTRINSECMEANVAR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

void Intrinsec_Mean_Var(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                    vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha);

void Intrinsec_Mean_Var_tDiv(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                    vector<string> vec_activities, int NFEATURES, int num_mean, int t_Div);

void Intrinsec_Mean_Var_tDiv_Iterative(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                                       vector<string> vec_activities, int NFEATURES, int num_mean, int t_Div, int max_iter);


void Intrinsec_mean_VectoMapping(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                        vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha);


void Intercambiar_cov(float** K_sort, int i, int j);
void quickSort_cov(float** K_sort,  int left, int right);

#endif // INTRINSECMEANVAR_H
