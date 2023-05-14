#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <string>
#include <ctime>

#include <vector>
#include <array>

namespace fs {
  std::string joinString(const std::string &src1,
                         const std::string &src2);
} // namespace fs

class Timer{
public:
  Timer(std::string namefunc) {
    this->start_time = clock();
    this->namefunc = namefunc;
  }
  ~Timer() {
    std::cout << "Time(" << namefunc << "): "<< clock() - start_time<< " ms\n";
  }
private:
  int start_time;
  std::string namefunc;
};

#ifdef TOOLS_TIMER
#define CheckRunTime(nameFunc) \
  Timer timer(nameFunc);
#else
#define CheckRunTime(nameFunc)
#endif

#define CheckAssert(a) _my_assert(a, __FILE__, __LINE__, __FUNCTION__)
void _my_assert(const bool expr, const char *file, const int line, const char *func);

namespace ReceiversType {
  typedef enum {
    DISPLACEMENT = 0,
    VELOCITY = 1,
    PRINCIPAL_STRESS = 2,
    PRESSURE = 3,
    ACCELERATION = 4,
    COUNT
  } eType;
}

struct RECEIVER {
  ReceiversType::eType type;
  std::array<bool, 3> dofs = { false, false, false };
  std::vector<int> glob_nodes;
  std::vector<std::array<double,3>> coords;
  int id = -1;
  int output_every_step = 1;
};

#ifdef _MSC_VER
#define native_to_big(v) _byteswap_ulong(v); // Visual C++ 32 bit
#else
#define __builtin_bswap32(v); // GCC 32 bit
#endif

const std::string ReceiversTypeToInfoName(const ReceiversType::eType type, const size_t dof);
//ResultsType::eType ReceiversTypeToResultsType(const ReceiversType::eType type);

#endif
