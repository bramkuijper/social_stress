#include "stress_social.hpp"

int main(int argc, char **argv)
{
    Parameters pars;

    pars.file_name = argv[1];

    StressSocial sim_object(pars);

}
