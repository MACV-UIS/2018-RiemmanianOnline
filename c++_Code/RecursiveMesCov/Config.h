
#ifndef CONFIG_H
#define CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

using namespace std;


///--------------------------------  KTH -----------------------------------------
//string root                       = "/home/fmartinezc/POSDOC/DATA/KTH/";
////string root                     = "/people/fmartinez/POSDOC/DATA/KTH/";
////string root                     = "/home/fmartinezc/main/datasets/KTH/";
//string str_path_train           = "/home/macv/DATOS/fabio/ConfigurationFiles/KTH_trainingSeq1.txt";
//string str_path_test            = "/home/macv/DATOS/fabio/ConfigurationFiles/KTH_testSeq1.txt";
//string str_path_valid           = "/home/macv/DATOS/fabio/ConfigurationFiles/KTH_validationSeq1.txt";
//string str_path_activities      = "/home/macv/DATOS/fabio/ConfigurationFiles/activities.txt";
//string str_path_frames          = root + "/FramesSeq/";
//
////                                   "/EnstaTraj/Th30/";
//                                     //"/EnstaTraj/Th15/"
////                                   // "/DenseTrajSubSeq/DENSE/ -  FILTERED";// "
//                                    // ImprovedDenseTrajSubSeq/DENSE
//
//string str_path_trajectories    = "/home/macv/DATOS/fabio/improved/";// "/DenseTrajSubSeq/FILTERED/";//"/EnstaTraj/Th15/";//"/ImprovedDenseTrajSubSeq/FILTERED/";//"/DenseTrajSubSeq/FILTERED/";// "
//string str_path_expAct          = root +  "/ConfigurationFiles/KTH_ExpActionsSeq1.txt";


///--------------------------------- UT ------------------------------------------------------

//

string root                     = "/home/macv/DATOS/gustavo/datasets/UT-interaction/";
//string str_path_train           = root +  "ConfigurationFiles/segmented_set1_rand.txt";
string str_path_train            = root +  "/ConfigurationFiles/segmented_set2_rand.txt";

string str_path_activities      = root +  "ConfigurationFiles/activities.txt";
//string str_path_frames          = root + "UT_segmented_set1/FramesSeq/";
//
////                                   "/EnstaTraj/Th30/";
//                                     //"/EnstaTraj/Th15/"
////                                   // "/DenseTraj/DENSE/ -  FILTERED";// "
//                                    // ImprovedDenseTraj/DENSE
//
//string str_path_trajectories    = root + "ImprovedDenseTraj/segmented_set1/FILTERED/";// "
string str_path_trajectories    = root + "ImprovedDenseTraj/segmented_set2/FILTERED/";// "
////string str_path_expAct;             //= root +  "UT_segmented_set1/ConfigurationFiles/KTH_ExpActionsSeq1.txt";


#endif // CONFIG_H
