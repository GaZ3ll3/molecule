/*
 *
 *  Copyright (C) 2016, Yimin Zhong, University of Texas at Austin.
 *
 *  Quadrature rules for singular integral
 *
 */

#ifndef MOLECULE_QUADRATUE_H
#define MOLECULE_QUADRATUE_H

#include "mapping.h"

inline scalar_t singularAngle(point& a, point &b, int id, vector<point>& centers) {
    return angle(a, b, centers[id]);
}

#endif //MOLECULE_QUADRATUE_H
