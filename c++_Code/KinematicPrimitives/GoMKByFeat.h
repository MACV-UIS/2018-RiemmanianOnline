#ifndef GOMKINEMATICFEATURES_H
#define GOMKINEMATICFEATURES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <vector>

#include <cv.h>
#include <highgui.h>
#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

#define PI  3.14159265358979323846;


vector< vector< float > > computeMoGFeat(vector<float> vect_traj_x, vector<float> vect_traj_y, vector<float> vect_traj_z,
                                         string feature_name, int numberOfGaussians);

void computeMoGPerTraj(vector<float> vect_traj_x, vector<float> vect_traj_y, vector<float> vect_traj_z,
                                         string feature_name, int numberOfGaussians, float*** matrixMoGTraj,
                                         float alphaMoG, float std_init, float w_init, float rho_init );

void init_MatrixMoGTraj(float*** matrixMoGTraj,
                          int numberOfGaussians, float std_init, float w_init, float mean_init);

void sort_MatrixMoGTraj(float*** matrixMoGTraj, int numberOfGaussians, int t);

void norm_MatrixMoGTraj(float*** matrixMoGTraj, int numberOfGaussians, int t);

float val_abs(float val);

float val_meanMoG( float val_m, float val, float alpha);
float val_stdMoG( float sigma_m, float val_m, float val, float alpha);

#endif // GOMKINEMATICFEATURES_H
