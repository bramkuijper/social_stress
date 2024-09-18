#include <cassert>
#include "stress_social.hpp"

int main(int argc, char **argv)
{
    Parameters pars; // Assigning parameter order
    pars.file_name = argv[1];
    pars.s[NP] = std::stod(argv[2]); // String to double
    pars.s[P] = std::stod(argv[3]);
    pars.npatches = std::stoi(argv[4]); // String to integer

    StressSocial sim_object(pars);

}
