#ifndef _STRESS_SOCIAL_HPP
#define _STRESS_SOCIAL_HPP

#include <string>
#include <vector>
#include <fstream>
#include <random>
#include "parameters.hpp"
#include "patch.hpp"

class StressSocial
{
    private:
        // stuff for random number generation
        unsigned time_step{0};

        // parameters for this simulation
        Parameters param{};
        
        // random device which is used to generate
        // proper random seeds
        std::random_device rd;
        
        // store the random seed
        // we need to store this so that we can output the
        // random seed, so that we could 'replay' the exact
        // same sequence of random numbers for debugging purposes etc
        unsigned int seed;

        // random number generator
        std::mt19937 rng_r;

        // uniform distribution to compare against probabilities
        std::uniform_real_distribution<double> uniform;

        // uniform distribution to get random patch
        std::uniform_int_distribution<int> patch_sampler;
        
        // a metapopulation of patches containing individuals
        std::vector <Patch> metapopulation;

        // functions for the life cycle: survival, replacement, etc
        void initialize_patches();
        void predator_visit();
        void switch_predator_status();

        void survive(); 
        void reproduce(); 


        // functions for data output
        std::ofstream data_file;

        void write_parameters();
        void write_data_headers();
        void write_data();

    public:
        StressSocial(Parameters const &param);

};

#endif
