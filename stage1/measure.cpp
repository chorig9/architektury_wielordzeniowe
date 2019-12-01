#include "constants.hpp"
#include "simulate.hpp"

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
    simulation sim(4000, time(NULL));

    size_t iterations = 1000;
    if (argc >= 2)
        iterations = std::stoull(argv[1]);

    std::cout << measure<>::execution([&]{
        for (size_t i = 0; i < iterations; i++)
            sim.step();
    }) << std::endl;
}