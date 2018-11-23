#ifndef CONCATRECLOGCOV_H
#define CONCATRECLOGCOV_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;

void Concat_Rec_Log_CoV(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                    vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha);
void Framelevel_Rec_Log_CoV(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                    vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha);
void Concat_Scal_Rec_Log_CoV(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                    vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha, int num_scales);

#endif
