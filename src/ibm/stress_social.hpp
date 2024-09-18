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
        std::uniform_int_distribution<unsigned> patch_sampler;

        // uniform distribution to sample breeders from patch
        std::uniform_int_distribution<unsigned> take_random_breeder;
        
        // a metapopulation of patches containing individuals
        std::vector <Patch> metapopulation;

        // functions for the life cycle: survival, replacement, etc
        void initialize_patches();
        void predator_visit();
        void switch_predator_status();

        double calculate_group_vigilance(Patch const &current_patch);
        double attack_survival(double const h);
        double update_stress_hormone();

        // life history functions
        void survive_damage_vigilance();
        void reproduce(); 
        void update_damage();

        // damage and vigilance-investment 
        // related mortality
        double mu(double const damage, 
                double const vigilance);

        // functions for data output
        std::ofstream data_file;

        void write_parameters();
        void write_data_headers();
        void write_data();

    public:
        StressSocial(Parameters const &param);

};

#endif
