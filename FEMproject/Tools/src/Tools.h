#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <string>
#include <ctime>

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

#endif
