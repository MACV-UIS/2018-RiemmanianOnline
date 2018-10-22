
#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;


///--------------------------------  KTH -----------------------------------------
//string root                       = "/home/fmartinezc/POSDOC/DATA/KTH/";
//string root                     = "/people/fmartinez/POSDOC/DATA/KTH/";
string root                     = "/home/fmartinezc/main/datasets/KTH/";
string str_path_train           = root +  "/ConfigurationFiles/KTH_trainingSeq1.txt";
string str_path_test            = root +  "/ConfigurationFiles/KTH_testSeq1.txt";
string str_path_valid           = root +  "/ConfigurationFiles/KTH_validationSeq1.txt";
string str_path_activities      = root +  "/ConfigurationFiles/activities.txt";
string str_path_frames          = root + "/FramesSeq/";


//DenseTraj
//EnstaTraj
//ImpTraj

string str_path_trajectories    = root + "//";// "/DenseTrajSubSeq/FILTERED/";//"/EnstaTraj/Th15/";//"/ImprovedDenseTrajSubSeq/FILTERED/";//"/DenseTrajSubSeq/FILTERED/";// "
string str_path_expAct          = root +  "/ConfigurationFiles/KTH_ExpActionsSeq1.txt";


#endif // CONFIG_H
