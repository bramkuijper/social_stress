#ifndef _PATCH_HPP 
#define _PATCH_HPP

#include <vector>
#include "individual.hpp"

// TODO how are these breeders actually initialized

class Patch
{
    public:
        bool predator_patch{false};

        std::vector <Individual> breeders[2];

        std::vector <Individual> breeders_surviving[2];

        Patch();


};

#endif 
