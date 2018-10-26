#ifndef RECURSIVE_WEKA_INTRINSECMEANVAR_H
#define RECURSIVE_WEKA_INTRINSECMEANVAR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;


void Intrinsec_mean_Weka(vector<string>  vec_videos, string name_file,  string str_path_trajectories,
                        vector<string> vec_activities, int NFEATURES, int num_mean, float init_alpha);

#endif // RECURSIVE_WEKA_INTRINSECMEANVAR_H
