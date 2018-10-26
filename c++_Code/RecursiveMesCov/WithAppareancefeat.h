#ifndef WITHAPPAREANCEFEAT_H
#define WITHAPPAREANCEFEAT_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;


void FromAppareance_Feat(vector<string>  vec_videos, string name_file,  string str_path_trajectories, string path_frames,
                        vector<string> vec_activities, int NFEATURES, int num_mean, bool WRTT);
#endif
