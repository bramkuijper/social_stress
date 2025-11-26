#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <random>
#include "parameters.hpp"

class Individual
{
    private:

    public:
        bool is_alive{true};
        bool is_attacked{false};
        double baseline_influx[2]{0.0,0.0}; 
        double stress_influx[2]{0.0,0.0}; 
        double vigilance_influx[2]{0.0,0.0}; 
        double removal[2]{0.0,0.0};

// TODO EG: Evolvable baseline vigilance alleles (a_v). No b_stressvigilance trait yet as v is baseline only
        double v[2]{0.0,0.0}; // baeline vigilance alleles

        // state variables
        double damage{0.0}; // 
        double stress_hormone{0.0}; // how much stress hormone does an individual have now

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
