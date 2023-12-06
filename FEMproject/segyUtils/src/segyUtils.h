
#ifndef SEGY_UTILS_H_
#define SEGY_UTILS_H_

#include <string>
#include <vector>
#include <set>


#include <segy.h>
#include <Tools.h>

#include <matrix_pack.h>

void convert(std::string out_path, std::string txt_filename, const ReceiversType::eType type, Matrix &nodes, size_t dof, size_t elementsCount, size_t nodesCount);

#endif  // SEGY_UTILS_H_
