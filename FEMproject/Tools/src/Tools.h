#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <ctime>

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


#endif
