#include <iostream>
#include <math.h>

#include "ManageMat_OPS.h"
#include "Geodesics.h"
#include "Riemannian_distances.h"

using namespace std;

float comp_AffineInv_Dist_AID(MatDoub &input_A, MatDoub &input_B)
{

    int NFEATURES = input_A.nrows();

    MatDoub L(NFEATURES, NFEATURES), Li(NFEATURES, NFEATURES), L_transpose(NFEATURES,NFEATURES), L_transpose_inverse(NFEATURES, NFEATURES), M2(NFEATURES, NFEATURES);


    try
    {
        choleskyd(input_A, L);
    }
    catch (const char* str)
    {
        cout << "Error caught: " << str << endl;
//		displayMat(A, "CovMat");
//		displayMat(B, "Model");
        return 100000;
    }
    catch (int str)
    {
        cout << "Error caught: " << str << endl;
//		displayMat(A, "CovMat");
//		displayMat(B, "Model");
        //    cout << "Similarity measure = " << 100000 <<endl;
        return 100000;
    }

    /*cout<<"L is :"<<endl;
    displayMat(L);*/
    inverse(L, Li);
    /*cout<<"L_inverse is "<<endl;
     displayMat(Li);*/

    transpose_of_any_matrix(L, L_transpose);
    /*cout<<"L transpose is"<<endl;
     displayMat(L_transpose);*/

    inverse(L_transpose, L_transpose_inverse);
    /*cout<<"L transpose inverse is"<<endl;
     displayMat(L_transpose_inverse);*/

    //M2 = Li * B * L_transpose_inverse;
    matrix_multiplication(Li, input_B, M2);
    matrix_multiplication(M2, L_transpose_inverse, M2);

    VecDoub Eigen_values(NFEATURES);
    comp_eigValues(M2, Eigen_values);


    float similarity = 0;

    for (int q = 0; q < NFEATURES; q++)
    {
        //cout<<p[q]<<endl;
        similarity += (log(Eigen_values[q]) * log(Eigen_values[q]));
    }
    similarity = sqrt(similarity);

    return similarity;
}


float comp_Frobenius_Dist_LED(MatDoub &input_A, MatDoub &input_B)
{

    int NFEATURES = input_A.nrows();
    MatDoub Log_A(NFEATURES, NFEATURES, 0.0);
    MatDoub Log_B(NFEATURES, NFEATURES, 0.0);

    float distAB =0;

    comp_Log_A_resp_I(input_A, Log_A, distAB);
    comp_Log_A_resp_I(input_B, Log_B, distAB);

    float similarity =0;

    for (int i = 0; i < NFEATURES; i++)
    {
        for (int j = 0; j < NFEATURES; j++)
        {
            similarity = similarity + pow((Log_A[i][j] - Log_B[i][j]),2) ;
        }
    }

    return sqrt(similarity);

}


float difussion_dist(MatDoub &input_A)
{

    int NFEATURES = input_A.nrows();
    VecDoub Eigen_values(NFEATURES);
    comp_eigValues(input_A, Eigen_values);
    float similarity = 0;

    for(int i=0; i<NFEATURES; i++) similarity = similarity + pow(Eigen_values[i],2);

    return similarity;
}
