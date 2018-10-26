#ifndef APPERFEATTRAJ_H
#define APPERFEATTRAJ_H
/// compute apareance features from the video sequence to after compute the

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;


void ToComputeFrameApp( float** MatTraj, float** MatAppTraj, int NFEATURES, int num_tot_traj_k, int heightFrame_k, int widthFrame_k, string pathFileFrames, int frame_index);


void ToLoadMatTrajAppFeatures(int star_seq, int end_seq, string pathFileTrajectories, string pathFileFrames, int height_frame,
                       int width_frame, int num_tot_traj, int num_scales_traj, float*** MatTraj, float*** MatAppTraj,  bool WRTT);
#endif // APPERFEATTRAJ_H
