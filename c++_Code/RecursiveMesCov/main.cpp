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
#include "RecursiveCov.h"

#include "IntrinsecMeanVar.h"
#include "Weka_IntrinsecMeanVar.h"
#include "Recursive_Riemannian_Cov.h"
#include "ConcatRecLogCov.h"
#include "WithAppareancefeat.h"

#include "../Resources/GeneralFunctions.h"


using namespace std;

vector<string> vec_tes_train_val, vec_activities, vec_examplesActions;
vector<string> vec_test, vect_training_val;

void ToLoadFiles();
int NFEATURES=12; ///incluyendo x, y,t //((int)featLabels.size()*(int)alphaScales.size())+3; 9 features + x,y,z

int main()
{

    ToLoadFiles();
    cout<<" vec test: "<< vec_test.size() << " vec train: "<< vect_training_val.size() << " vec activities: "<< vec_activities.size() << endl;

//    Recursive_Riemannian_mean(vec_test, "test_Riemannian_Th_30_23",  str_path_trajectories, vec_activities,  NFEATURES, 1, 3.0);
//    Recursive_Riemannian_mean(vect_training_val, "train_Riemannian_Th30_23",  str_path_trajectories, vec_activities,  NFEATURES, 1, 3.0);

//    Recursive_mean(vec_test, "test_IF_26",  str_path_trajectories, vec_activities,  NFEATURES, 1, 10.0);
//    Recursive_mean(vect_training_val, "train_IF_26",  str_path_trajectories, vec_activities,  NFEATURES, 1, 10.0);

   // Intrinsec_Mean_Var(vec_test, "mv_test_Th_30_",  str_path_trajectories, vec_activities,  NFEATURES, 3, 5.0);
   // Intrinsec_Mean_Var(vect_training_val, "mv_train_Th_30_",  str_path_trajectories, vec_activities,  NFEATURES, 3, 5.0);

//    Intrinsec_Mean_Var_tDiv(vec_test, "tDiv_test_Th_10_",  str_path_trajectories, vec_activities,  NFEATURES, 8, 4);
//    Intrinsec_Mean_Var_tDiv(vect_training_val, "tDiv_train_Th_10_",  str_path_trajectories, vec_activities,  NFEATURES, 8, 4);

//      Intrinsec_Mean_Var_tDiv_Iterative(vect_training_val, "train_tIter_Th_30",  str_path_trajectories, vec_activities,  NFEATURES, 10, 4, 20);
//      Intrinsec_Mean_Var_tDiv_Iterative(vec_test, "test__tIter_Th_30",  str_path_trajectories, vec_activities,  NFEATURES, 10, 4, 20);


//      Concat_Rec_Log_CoV(vec_test, "test_concMes_DF_26_",  str_path_trajectories, vec_activities,  NFEATURES, 1, 10.0);
//      Concat_Rec_Log_CoV(vect_training_val, "train_concMes_DF_26_",  str_path_trajectories, vec_activities,  NFEATURES, 1, 10.0);


      Concat_Scal_Rec_Log_CoV(vec_test, "test_concScal_Th_15_22_5_",  str_path_trajectories, vec_activities,  NFEATURES, 1, 3.0, 2);
      Concat_Scal_Rec_Log_CoV(vect_training_val, "train_concScal_Th15_22_5_",  str_path_trajectories, vec_activities,  NFEATURES, 1, 3.0, 2);


     //Intrinsec_mean_VectoMapping(vec_test, "test_VectMap_IF_",  str_path_trajectories, vec_activities,  NFEATURES, 3, 5.0);
     //Intrinsec_mean_VectoMapping(vect_training_val, "train_VectMap_IF_",  str_path_trajectories, vec_activities,  NFEATURES, 3, 5.0);
///     Intrinsec_mean_VectoMapping(vec_test, "test_VectMap_Div_",  str_path_trajectories, vec_activities,  NFEATURES, 3, 5.0);
///     Intrinsec_mean_VectoMapping(vect_training_val, "train_VectMap_Div_",  str_path_trajectories, vec_activities,  NFEATURES, 3, 5.0);



///-------------------  not going to use appearance for now -------------------

//    FromAppareance_Feat(vec_test, "test_App_ID_",  str_path_trajectories, str_path_frames,
//                        vec_activities,  NFEATURES, 3, true);
//
//    FromAppareance_Feat(vect_training_val, "training_App_ID_",  str_path_trajectories, str_path_frames,
//                        vec_activities,  NFEATURES, 3, true);


///--------------------  UT --------------------------------
//     Intrinsec_mean_VectoMapping(vec_test, "test_XYt-TS-NxNy-k_DF_",  str_path_trajectories, vec_activities,  NFEATURES, 3, 5.0);
//     Intrinsec_mean_VectoMapping(vect_training_val, "train_XYt-TS-NxNy-k_DF_",  str_path_trajectories, vec_activities,  NFEATURES, 3, 5.0);

     //XYTxTy
//     Intrinsec_mean_Weka(vec_test, "weka_test_DF_",  str_path_trajectories, vec_activities,  NFEATURES, 3, 5.0);
//     Intrinsec_mean_Weka(vect_training_val, "weka_train_DF_",  str_path_trajectories, vec_activities,  NFEATURES, 3, 5.0);

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
    //redFileSequences(str_path_expAct, vec_examplesActions);
}
