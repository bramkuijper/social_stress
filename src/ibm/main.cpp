#include <cassert>
#include "stress_social.hpp"

int main(int argc, char **argv)
{
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
    pars.init_v = std::stod(argv[11]); // initial vigilance
    pars.init_stress_hormone_level = std::stod(argv[12]); // initial stress hormone level
    pars.g = std::stod(argv[13]); // damage removal per timestep
    pars.k = std::stod(argv[14]); // increase ini damage due to hormone != optimum
    
    StressSocial sim_object(pars);

}
