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

#include "Config.h"
#include "RecursiveCovStat.h"


#include "../Resources/GeneralFunctions.h"

using namespace std;

vector<string> vec_tes_train_val, vec_activities, vec_examplesActions;
vector<string> vec_test, vect_training_val;

void ToLoadFiles();

int KIN_NFEATURES=12;///incluyendo x, y,t //((int)featLabels.size()*(int)alphaScales.size())+3; 9 features + x,y,z
int APP_NFEATURES=6;
/// La idea principal es contruir los experimentos para el paper, desarrollando unicamente los datos que seran utilizados


int main()
{
    ToLoadFiles();
    cout<<" vec test: "<< vec_test.size() << " vec train: "<< vect_training_val.size() << " vec activities: "<< vec_activities.size() << endl;

    /// To compute training complete sequences
    //Concat_Rec_Log_CoV(vec_test, "test_concMes_DF_25_",  str_path_trajectories, vec_activities,  KIN_NFEATURES, 1, 5.0);
    //Concat_Rec_Log_CoV(vect_training_val, "train_concMes_DF_25_",  str_path_trajectories, vec_activities,  KIN_NFEATURES, 1, 5.0);

    Concat_Rec_Log_CoV_splitLearning(vect_training_val, "Gtrain_concMes_IFs_25_",  str_path_trajectories, vec_activities,  KIN_NFEATURES, 1, 5.0);
    Concat_Rec_Log_CoV_splitLearning(vec_test, "Gtest_concMes_IFs_24_",  str_path_trajectories, vec_activities,  KIN_NFEATURES, 1, 4.0);


//    Concat_Scal_Rec_Log_CoV_WithApp(vec_test, "test_concScal_Th_15_22_5_",  str_path_trajectories, vec_activities,  NFEATURES, 1, 2.0, 5);
//    Concat_Scal_Rec_Log_CoV_WithApp(vect_training_val, "train_concScal_Th15_22_5_",  str_path_trajectories, vec_activities,  NFEATURES, 1, 2.0, 5);

    /// To compute training segemented sequences with ratio of alpha


    /// To build test complete sequences

    /// To compute the prediction online using complete sequences as model - max label
    /// To compute the prediction online using segemented sequences as model - max label

    /// To compute global index for

    cout << "Termino todo ok!" << endl;
    return 0;
}

void ToLoadFiles()
{
    redFileSequences(str_path_train, vec_tes_train_val);
    redFileSequences(str_path_test, vec_tes_train_val);
    redFileSequences(str_path_valid, vec_tes_train_val);

    redFileSequences(str_path_test, vec_test);

    redFileSequences(str_path_train, vect_training_val);
    redFileSequences(str_path_valid, vect_training_val);

    redFileSequences(str_path_activities, vec_activities);
    redFileSequences(str_path_expAct, vec_examplesActions);
}
