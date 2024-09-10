#include <algorithm>
#include "individual.hpp"

// main constructor
Individual::Individual(Parameters const &params) :
    stress_hormone{params.init_stress_hormone_level/2, params.init_stress_hormone_level/2},
    v{params.init_v/2, params.init_v/2},
{}

// copy constructor
Individual::Individual(Individual const &other) :
    is_attacked{other.is_attacked},
    stress_hormone{other.stress_hormone[0], other.stress_hormone[1]},
    damage{other.damage},
    v{other.v[0],other.v[1]}
{}

// birth constructor
Individual::Individual(Individual const &mum,
                Individual const &dad,
                Parameters const &param,
                std::mt19937 &rng_r) :
{
    std::bernoulli_distribution segregator{0.5};

    stress_hormone[0] = mutate(mum.stress_hormone[segregator(rng_r)], param.mu_h, param.sdmu);
    stress_hormone[0] = std::clamp(stress_hormone[0], param.hmin, param.hmax);

    stress_hormone[1] = mutate(dad.stress_hormone[segregator(rng_r)], param.mu_h, param.sdmu);
    stress_hormone[1] = std::clamp(stress_hormone[1], param.hmin, param.hmax);

    v[0] = mutate(mum.v[segregator(rng_r)], param.mu_v, param.sdmu);
    v[0] = std::clamp(v[0], 0, 1.0);
    v[1] = mutate(dad.v[segregator(rng_r)], param.mu_v, param.sdmu);
    v[1] = std::clamp(v[1], 0, 1.0);

} // birth constructor
  
void Individual::operator=(Individual const &other)
{
    is_attacked = other.is_attacked;
    damage = other.damage;

    for (unsigned allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        stress_hormone[allele_idx] = other.stress_hormone[allele_idx];
        v[allele_idx] = other.v[allele_idx];
    }
} // end operator=()

// mutate according to a continuum-of-alleles model
double Individual::mutate(double to_mutate,
        double const mutation_prob,
        double const mutation_sd,
        std::mt19937 &rng_r)
{
    std::uniform_real_distribution<double> uniform{0.0,1.0};

    if (mutation_prob < uniform(rng_r))
    {
        std::normal_distribution<double> mutational_effect_size{0.0,mutation_sd};

        // our new allelic value
        to_mutate += mutational_effect_size(rng_r);
    }

    return(to_mutate);
} // end mutate 


