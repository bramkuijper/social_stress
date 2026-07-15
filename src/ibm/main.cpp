#include <cassert>
#include <iostream>
#include "stress_social.hpp"



int main(int argc, char **argv)
{

    if (argc != 18)
    {
        std::cerr
            << "Usage: " << argv[0]
            << " file_name npatches n max_time s_np s_p md mv"
            << " p_mig p_attack fecundity_power hmax init_v"
            << " init_stress_hormone g k vigilance\n"
            << "vigilance must be 1 (on) or 0 (off)\n";
    
        return 1;
    }

    Parameters pars; // Assigning parameter order
    pars.file_name = argv[1]; // file name
    pars.npatches = std::stoi(argv[2]); // number of patches
    pars.n = std::stoi(argv[3]); // individuals per patch
    pars.max_time = std::stoul(argv[4]); // max time
    pars.s[NP] = std::stod(argv[5]); // switch rate NP to P
    pars.s[P] = std::stod(argv[6]); // switch rate P to NP
    pars.md = std::stod(argv[7]); // weight of damage-related mortality
    pars.mv = std::stod(argv[8]); // weight of vigilance-related mortality
    pars.p_mig = std::stod(argv[9]); // migration probability
    pars.p_attack = std::stod(argv[10]); // probability of being attacked when predator present
    pars.fecundity_power = std::stod(argv[11]); // power of fecundity cost of vigilance
    pars.hmax = std::stod(argv[12]); // maximum stress hormone level
    pars.init_v = std::stod(argv[13]); // initial vigilance
    pars.init_stress_hormone_level = std::stod(argv[14]); // initial stress hormone level
    pars.g = std::stod(argv[15]); // damage removal per timestep
    pars.k = std::stod(argv[16]); // increase in damage due to hormone != optimum
    pars.vigilance = std::stoi(argv[17]) != 0; // vigilance on/off: 1 = on, 0 = off
    
    StressSocial sim_object(pars);
    return 0;

}
