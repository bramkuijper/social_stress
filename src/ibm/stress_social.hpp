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
        
        // data file for output 
        std::ofstream data_file;

        // some data members to keep track of numbers of attacks, mortalities etc
        unsigned int n_attacked{0};
        unsigned int n_death_damage{0};
        unsigned int n_death_predator{0};

        // functions for the life cycle: survival, replacement, etc
        void initialize_patches();
        void predator_visit();
        void switch_predator_status();

        double calculate_group_vigilance(Patch const &current_patch);
        double attack_survival(double const h);
        void update_stress_hormone();

        // life history functions
        void survive_damage_vigilance();
        void reproduce(); 

        // damage and vigilance-investment 
        // related mortality
        double mu(double const damage, 
                double const vigilance);

        // NOTE: Not sure if this needs including twice?
        void write_parameters();
        void write_data_headers();
        void write_data();

        void write_distribution();

    public:
        StressSocial(Parameters const &param);

};

#endif
