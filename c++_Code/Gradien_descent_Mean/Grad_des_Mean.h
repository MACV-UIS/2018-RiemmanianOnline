#ifndef GRADESMEAN_H
#define GRADESMEAN_H

#include "ManageMat_OPS.h"

void compute_grad_mean_Pennec(vector<MatDoub> vect_Matdoub, MatDoub &mean_deg);
void compute_grad_mean_Fletcher(vector<MatDoub> vect_Matdoub, MatDoub &mean_deg);

void compute_grad_vectorMapping(vector<MatDoub> vect_Matdoub, MatDoub &vect_map);
#endif
