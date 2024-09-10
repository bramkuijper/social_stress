#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <random>
#include "parameters.hpp"

class Individual
{
    private:

    public:
        bool is_attacked{false};
        double stress_hormone[2]{0.0,0.0}; // stress hormone
        double damage{0.0};
        double v[2]{0.0,0.0}; // vigilance

        // constructor
        Individual(Parameters const &params);

        // copy constructor
        Individual(Individual const &other); 

        // birth constructor
        Individual(Individual const &mum,
                Individual const &dad,
                Parameters const &param,
                std::mt19937 &rng_r);

        // mutation function
        double mutate(double to_mutate,
                double const mutation_prob,
                double const mutation_sd,
                std::mt19937 &rng_r);

        // assignment operator
        void operator=(Individual const &other);
};

#endif
