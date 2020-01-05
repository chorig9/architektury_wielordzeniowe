#include "constants.hpp"
#include "simulate.cu"

#include <iostream>
#include <chrono>

template<typename TimeT = std::chrono::milliseconds>
struct measure
{
    template<typename F, typename ...Args>
    static typename TimeT::rep execution(F&& func, Args&&... args)
    {
        auto start = std::chrono::steady_clock::now();
        std::forward<decltype(func)>(func)(std::forward<Args>(args)...);
        auto duration = std::chrono::duration_cast< TimeT> 
                            (std::chrono::steady_clock::now() - start);
        return duration.count();
    }
};

int main(int argc, char *argv[])
{
    simulation sim(10000, time(NULL));

    size_t iterations = 3000;
    if (argc >= 2)
        iterations = std::stoull(argv[1]);

    std::cout << measure<>::execution([&]{
        sim.step(iterations);
    }) << std::endl;
}