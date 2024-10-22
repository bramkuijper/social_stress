#include "patch.hpp"

Patch::Patch(Parameters const &params) :
    breeders{params.n, Individual(params)}
{}

Patch::Patch(Patch const &other):
    predator_patch{other.predator_patch},
    V{other.V},
    breeders{other.breeders},
    breeders_surviving{other.breeders_surviving}
{}

void Patch::operator=(Patch const &other)
{
    predator_patch = other.predator_patch;
    V = other.V;
    breeders = other.breeders;
    breeders_surviving = other.breeders_surviving;
}
