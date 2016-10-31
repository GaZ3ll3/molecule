//
// Created by Yimin Zhong on 10/29/16.
//

#ifndef MOLECULE_UTILS_H
#define MOLECULE_UTILS_H

#include <iostream>
#include <fstream>
#include <sstream>

#include <cmath>
#include <cstring>
#include <cassert>

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <utility>
#include <functional>

#include <chrono>
#include <iomanip>

#include "cblas.h"

#ifdef RUN_OMP
#include "omp.h"
#endif

typedef int index_t;
typedef double scalar_t;
typedef bool bool_t;

#define EPS 1e-14

#ifndef DISP
#define RUN(s, func){\
func;\
}
#endif


#ifdef DISP
#define RUN(s, func){ \
std::chrono::steady_clock::time_point begin =std::chrono::steady_clock::now(); \
func;\
std::chrono::steady_clock::time_point end =  end = std::chrono::steady_clock::now();\
std::cout << std::setw(15)<< s << " "  << std::setprecision(5) << std::setw(8) << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000.0 << " seconds"<<std::endl;\
}
#endif

using std::unordered_set;
using std::vector;
using std::queue;
using std::max;
using std::min;

#endif //MOLECULE_UTILS_H
