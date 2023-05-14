
#include <Tools.h>

namespace fs {

void normalizeString(std::string &s) {
  if (s.length() <= 1) {
    return;
  }

  for (size_t i = s.length() - 1; i >= 1; --i) {
    if (s[i] == '/' && s[i - 1]  == '/') {
      s.erase(i, 1);
    }
  }
}

std::string joinString(const std::string &src1,
                const std::string &src2) {
  std::string str_path = src1 + "/" + src2;
  normalizeString(str_path);
  return str_path;
}

} // namespace fs

void _my_assert(const bool expr, const char *file, const int line, const char *func) {
  if (!expr) {
    fprintf(stderr, "Assertion failed in %s \nFile %s(%i)\n", func, file, line);
    exit(EXIT_FAILURE);
//    throw std::runtime_error(std::string("Assertion failed in ") + func + std::string(" \nFile ") + file + std::string("(") + std::to_string(line) + std::string(")\n"));
  }
}

const std::string ReceiversTypeToInfoName(const ReceiversType::eType type, const size_t dof) {
  switch (type) {
  case ReceiversType::DISPLACEMENT: {
    if (dof == 0)
      return "Displacement Along X";
    else if (dof == 1)
      return "Displacement Along Y";
    else if (dof == 2)
      return "Displacement Along Z";
    break;
  }
  case ReceiversType::VELOCITY: {
    if (dof == 0)
      return "Velocity Along X";
    else if (dof == 1)
      return "Velocity Along Y";
    else if (dof == 2)
      return "Velocity Along Z";
    break;
  }
  case ReceiversType::ACCELERATION: {
    if (dof == 0)
      return "Acceleration Along X";
    else if (dof == 1)
      return "Acceleration Along Y";
    else if (dof == 2)
      return "Acceleration Along Z";
    break;
  }
  case ReceiversType::PRINCIPAL_STRESS: {
    if (dof == 0)
      return "Principal Stress Along X";
    else if (dof == 1)
      return "Principal Stress Along Y";
    else if (dof == 2)
      return "Principal Stress Along Z";
    break;
  }
  case ReceiversType::PRESSURE:
    return "Pressure";
  default:
    break;
  }

  return "Undefined";
}

//ResultsType::eType ReceiversTypeToResultsType(const ReceiversType::eType type) {
//  switch (type) {
//  case ReceiversType::DISPLACEMENT:
//    return ResultsType::RECEIVER_DISPLACEMENT;
//  case ReceiversType::VELOCITY:
//    return ResultsType::RECEIVER_VELOCITY;
//  case ReceiversType::PRINCIPAL_STRESS:
//    return ResultsType::RECEIVER_PRINCIPAL_STRESS;
//  case ReceiversType::PRESSURE:
//    return ResultsType::RECEIVER_PRESSURE;
//  case ReceiversType::ACCELERATION:
//    return ResultsType::RECEIVER_ACCELERATION;
//  default:
//    break;
//  }

//  return ResultsType::COUNT;
//}
