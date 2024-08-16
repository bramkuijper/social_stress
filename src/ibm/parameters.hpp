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
        int n[2]{5,5};

        // dispersal
        double d[2]{0.1,0.9};

        int npatches{100}; 

        unsigned max_time{3};

        // switch rates between predator present in patch
        // vs absent
        double s[2]{0.5,0.5};

        // attack probability in a patch where a predator
        // is present

        double p_attack{0.1};

        // base name for the file
        std::string file_name{"sim_stress_social"};
};

#endif
