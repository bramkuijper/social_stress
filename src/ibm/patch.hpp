#ifndef _PATCH_HPP 
#define _PATCH_HPP

#include <vector>
#include "individual.hpp"

// TODO how are these breeders actually initialized

class Patch
{
    public:
        bool predator_patch{false};

        double V{0.0}; // group-level vigilance

        std::vector <Individual> breeders;

        std::vector <Individual> breeders_surviving;

        Patch(Parameters const &params); // Declaration of default constructor patch
        Patch(Patch const &other);
        void operator=(Patch const &other);

};

#endif 
