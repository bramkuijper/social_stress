#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

enum Sex
{
    female = 0,
    male = 1
};

enum PatchType
{
    NP = 0,
    P = 1
};

// all parameters with their default values
class Parameters
{
    public:
        // first element female, second element male
        unsigned n{20};

        // dispersal
        double d{0.1};

        unsigned npatches{100}; 

        unsigned max_time{3};

        // switch rates between predator present in patch
        // vs absent
        double s[2]{0.5,0.5};

        // attack probability in a patch where a predator
        // is present

        double p_attack{0.1};

        // power of how fecundity decreases with vigilance
        double fecundity_power{1.0};

        // min max hormone
        double hmin{0.0};
        double hmax{10.0};

        // max damage level
        double dmax{10.0};

        double survival_power{1.0};

        double init_v{0.0};
        double init_stress_hormone_level{0.0};

        // base name for the file
        std::string file_name{"sim_stress_social"};

        // mutation rates
        double mu_baseline{0.01};
        double mu_stress_influx{0.01};
        double mu_vigilance_influx{0.01};
        double mu_removal{0.01};
        double mu_v{0.01};
        double sdmu{0.01};

        // mortality rates 
        double m0{0.1};
        double md{1.0}; // weighting of damage-related mortality
        double mv{1.0}; // weighting of vigilance-investment-related mortality
};

#endif
