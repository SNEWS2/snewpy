
#include "EMEWS.h"

#ifndef jacobians_H
#define jacobians_H

MATRIX<double,10,8> J(std::array<double,NY> Y);
void J(std::array<double,NY> Y,MATRIX<double,10,8> &j);

MATRIX<double,8,10> JInverse(std::array<double,NY> Y);
void JInverse(std::array<double,NY> Y,MATRIX<double,8,10> &j);

#endif
