#include <cmath>
#include <cstring>
#include <sstream>
#include <queue>
#include <fstream>
#include <chrono>
#include "omp.h"
#include "getopt.h"
#ifdef ED
    #include "./tree/AVtreeString.h"
#else
    #include "./tree/AVtree.h"
#endif
#include <sys/stat.h>                                                     
#include <iostream> 
#include <functional>
#include <iomanip>


using namespace std;
using namespace std::chrono;

#define LINEAR_SCAN                                                 0
#define STANDARD                                                    1
#define STANDARD_MEDIOCRE                                           2
#define MEDIOCRE_CACHE                                              3

#define KNNMediocre                                                 4
#define KNNCrackOnKth                                               5