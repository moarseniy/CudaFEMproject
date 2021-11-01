#ifndef TESTS_H
#define TESTS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

#include "Tools.h"
#include "femstruct.h"

void test_EbePCG_diag(std::vector<Element> &elems, std::vector<float*> &m_e, std::vector<float*> &b_e, float soln[18]);

#endif
