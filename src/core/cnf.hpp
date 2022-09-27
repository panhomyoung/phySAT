/* phySAT: This is a SAT solver based on Semi-Tonser Product
 * Copyright (C) 2022 */

/**
 * @file cnf.hpp
 *
 * @brief Semi-Tensor Product based SAT solver for CNF solving
 *
 * @author Hongyang Pan
 * @since  2022/09/27
*/

#pragma once

#include "matrix.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>
#include <time.h>

using namespace std;

namespace phySAT
{
    void stp_cnf(vector<int> &expression, vector<string> &tt, vector<MatrixXi> &mtxvec);
}
